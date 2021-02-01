#include "casm/enumerator/io/dof_space_analysis.hh"

#include <boost/filesystem.hpp>

#include "casm/app/DirectoryStructure.hh"
#include "casm/casm_io/Log.hh"
#include "casm/clex/SimpleStructureTools.hh"
#include "casm/crystallography/SimpleStructure.hh"
#include "casm/crystallography/io/SimpleStructureIO.hh"
#include "casm/enumerator/ConfigEnumInput.hh"
#include "casm/enumerator/DoFSpace_impl.hh"
#include "casm/enumerator/io/json/DoFSpace.hh"
#include "casm/global/definitions.hh"
#include "casm/symmetry/json_io.hh"

namespace CASM {

namespace DoFSpaceIO {

OutputImpl::Error::Error(Index _state_index, std::string _identifier,
                         DoFKey _dof_key, std::string _what, jsonParser _data)
    : state_index(_state_index),
      identifier(_identifier),
      dof_key(_dof_key),
      what(_what),
      data(_data) {}

/// Write symmetry groups (lattice point group, factor_group,
/// crystal_point_group)
void OutputImpl::write_symmetry(Index state_index,
                                std::string const &identifier,
                                ConfigEnumInput const &config_enum_input,
                                std::vector<PermuteIterator> const &group) {
  Lattice config_lattice = config_enum_input.configuration().ideal_lattice();
  SymGroup lattice_point_group{SymGroup::lattice_point_group(config_lattice)};
  SymGroup factor_group =
      make_sym_group(group.begin(), group.end(), config_lattice);
  SymGroup crystal_point_group =
      make_point_group(group.begin(), group.end(), config_lattice);
  this->write_symmetry(state_index, identifier, config_enum_input,
                       lattice_point_group, factor_group, crystal_point_group);
}

void OutputImpl::write_dof_space_error(
    make_symmetry_adapted_dof_space_error const &e, Index state_index,
    DoFSpace const &dof_space, std::string const &identifier,
    ConfigEnumInput const &config_enum_input,
    std::optional<VectorSpaceSymReport> const &sym_report) {
  jsonParser json;
  to_json(dof_space, json, identifier, config_enum_input, sym_report);
  errors().emplace_back(state_index, identifier, dof_space.dof_key(), e.what(),
                        json);
}

/// Write dof space analysis error information to
/// <current_path>/dof_space_errors.json
void OutputImpl::write_errors() const {
  if (!errors().size()) {
    return;
  }

  fs::path output_path = fs::current_path() / "dof_space_errors.json";
  Index i = 1;
  while (fs::exists(output_path)) {
    output_path = fs::current_path() / (std::string("dof_space_errors.") +
                                        std::to_string(i) + ".json");
    ++i;
  }

  jsonParser json = jsonParser::array();
  for (auto const e : errors()) {
    jsonParser tjson;
    tjson["state_index"] = e.state_index;
    tjson["identifier"] = e.identifier;
    tjson["dof"] = e.dof_key;
    tjson["what"] = e.what;
    tjson["data"] = e.data;
    json.push_back(tjson);
  }
  CASM::log() << "Encountered " << errors().size() << " errors" << std::endl;

  fs::ofstream outfile{output_path};
  json.print(outfile);
  CASM::log() << "Writing: " << output_path << std::endl;
}

void DirectoryOutput::write_symmetry(Index state_index,
                                     std::string const &identifier,
                                     ConfigEnumInput const &config_enum_input,
                                     SymGroup const &lattice_point_group,
                                     SymGroup const &factor_group,
                                     SymGroup const &crystal_point_group) {
  _check_config(state_index, identifier, config_enum_input);
  fs::path output_dir = _output_dir(state_index, identifier, config_enum_input);

  // Write lattice point group
  {
    jsonParser json;
    write_symgroup(lattice_point_group, json);
    fs::path output_path = output_dir / "lattice_point_group.json";
    fs::ofstream outfile{output_path};
    json.print(outfile);
    CASM::log() << "Writing: " << output_path << std::endl;
  }

  // Write factor group
  {
    jsonParser json;
    write_symgroup(factor_group, json);
    fs::path output_path = output_dir / "factor_group.json";
    fs::ofstream outfile{output_path};
    json.print(outfile);
    CASM::log() << "Writing: " << output_path << std::endl;
  }

  // Write crystal point group
  {
    jsonParser json;
    write_symgroup(crystal_point_group, json);
    fs::path output_path = output_dir / "crystal_point_group.json";
    fs::ofstream outfile{output_path};
    json.print(outfile);
    CASM::log() << "Writing: " << output_path << std::endl;
  }
}

// void DirectoryOutput::write_config(
//   Index state_index,
//   std::string const &identifier,
//   ConfigEnumInput const &config_enum_input) {
//
//   _check_config(state_index, identifier, config_enum_input);
//   fs::path output_dir = _output_dir(state_index, identifier,
//   config_enum_input);
//
//   fs::path output_path = output_dir / "config.json";
//   fs::ofstream outfile {output_path};
//   jsonParser json;
//   to_json(config_enum_input.configuration(), json);
//   json.print(outfile);
//   CASM::log() << "Writing: " << output_path << std::endl << std::endl;
// }

void DirectoryOutput::write_structure(
    Index state_index, std::string const &identifier,
    ConfigEnumInput const &config_enum_input) {
  _check_config(state_index, identifier, config_enum_input);
  fs::path output_dir = _output_dir(state_index, identifier, config_enum_input);

  fs::path output_path = output_dir / "structure.json";
  fs::ofstream outfile{output_path};
  jsonParser json;
  to_json(make_simple_structure(config_enum_input.configuration()), json);
  json.print(outfile);
  CASM::log() << "Writing: " << output_path << std::endl << std::endl;
}

void DirectoryOutput::write_dof_space(
    Index state_index, DoFSpace const &dof_space, std::string const &identifier,
    ConfigEnumInput const &config_enum_input,
    std::optional<VectorSpaceSymReport> const &sym_report) {
  _check_config(state_index, identifier, config_enum_input);
  fs::path output_dir = _output_dir(state_index, identifier, config_enum_input);

  jsonParser json;
  to_json(dof_space, json, identifier, config_enum_input, sym_report);

  std::string filename = "dof_space_" + dof_space.dof_key() + ".json";
  fs::path output_path = output_dir / filename;
  json.write(output_path);
  CASM::log() << "Writing: " << output_path << std::endl << std::endl;
}

SymmetryDirectoryOutput::SymmetryDirectoryOutput(DirectoryStructure const &dir)
    : m_dir(dir) {}

/// For SymmetryDirectoryOutput, configurations must exist in database and this
/// will throw otherwise
void SymmetryDirectoryOutput::_check_config(
    Index state_index, std::string const &identifier,
    ConfigEnumInput const &config_enum_input) {
  Configuration const &configuration = config_enum_input.configuration();
  // These should not occur for now. They should be prevented by
  // require_database_configurations.
  // TODO: add support for other dof space analyses
  if (configuration.id() == "none") {
    throw std::runtime_error(
        "Error in output_dof_space: Cannot output to symmetry directory, "
        "configuration does not exist in database. Choose a different output "
        "method.");
  }
  if (identifier != configuration.name()) {
    throw std::runtime_error(
        "Error in output_dof_space: Cannot output to symmetry directory, "
        "unknown name error. Choose a different output method.");
  }
  if (config_enum_input.sites().size() != configuration.size()) {
    throw std::runtime_error(
        "Error in output_dof_space: Cannot output to symmetry directory, "
        "incomplete site selection. Choose a different output method.");
  }
}

/// For SymmetryDirectoryOutput return output directory path in symmetry dir
fs::path SymmetryDirectoryOutput::_output_dir(
    Index state_index, std::string const &identifier,
    ConfigEnumInput const &config_enum_input) {
  fs::path output_dir =
      m_dir.symmetry_dir(config_enum_input.configuration().name());
  fs::create_directories(output_dir);
  return output_dir;
}

SequentialDirectoryOutput::SequentialDirectoryOutput(fs::path output_dir)
    : m_output_dir(output_dir) {
  if (fs::exists(m_output_dir / "dof_space")) {
    throw std::runtime_error(
        "Error in output_dof_space: \"dof_space\" directory "
        "already exists. Will not overwrite.");
  }
}

/// For SequentialDirectoryOutput, any input state is allowed
void SequentialDirectoryOutput::_check_config(
    Index state_index, std::string const &identifier,
    ConfigEnumInput const &config_enum_input) {
  return;
}

/// For SequentialDirectoryOutput return output directory path
fs::path SequentialDirectoryOutput::_output_dir(
    Index state_index, std::string const &identifier,
    ConfigEnumInput const &config_enum_input) {
  std::string dirname = std::string("state.") + std::to_string(state_index);
  fs::path output_dir = m_output_dir / "dof_space" / dirname;
  fs::create_directories(output_dir);
  return output_dir;
}

CombinedJsonOutput::CombinedJsonOutput(fs::path output_dir)
    : m_combined_json(jsonParser::array()), m_output_dir(output_dir) {
  fs::path output_path = m_output_dir / "dof_space.json";
  if (fs::exists(output_path)) {
    throw std::runtime_error(
        "Error in output_dof_space: \"dof_space.json\" file "
        "already exists. Will not overwrite.");
  }
}

CombinedJsonOutput::~CombinedJsonOutput() {
  fs::path output_path = m_output_dir / "dof_space.json";
  fs::ofstream outfile{output_path};
  m_combined_json.print(outfile);
  CASM::log() << "Writing: " << output_path << std::endl;
}

void CombinedJsonOutput::write_symmetry(
    Index state_index, std::string const &identifier,
    ConfigEnumInput const &config_enum_input,
    SymGroup const &lattice_point_group, SymGroup const &factor_group,
    SymGroup const &crystal_point_group) {
  jsonParser &output_json = _output_json(state_index);
  write_symgroup(lattice_point_group, output_json["lattice_point_group"]);
  write_symgroup(factor_group, output_json["factor_group"]);
  write_symgroup(crystal_point_group, output_json["crystal_point_group"]);
}

// void CombinedJsonOutput::write_config(
//   Index state_index,
//   std::string const &identifier,
//   ConfigEnumInput const &config_enum_input) {
//
//   jsonParser &output_json = _output_json(state_index);
//   to_json(config_enum_input.configuration(), output_json["config"]);
// }

void CombinedJsonOutput::write_structure(
    Index state_index, std::string const &identifier,
    ConfigEnumInput const &config_enum_input) {
  jsonParser &output_json = _output_json(state_index);
  to_json(make_simple_structure(config_enum_input.configuration()),
          output_json["structure"]);
}

void CombinedJsonOutput::write_dof_space(
    Index state_index, DoFSpace const &dof_space, std::string const &identifier,
    ConfigEnumInput const &config_enum_input,
    std::optional<VectorSpaceSymReport> const &sym_report) {
  jsonParser &output_json =
      _output_json(state_index)["dof_space"][dof_space.dof_key()];
  to_json(dof_space, output_json, identifier, config_enum_input, sym_report);
}

jsonParser &CombinedJsonOutput::_output_json(Index state_index) {
  while (m_combined_json.size() < state_index) {
    m_combined_json.push_back(jsonParser::object());
  }
  return m_combined_json[state_index - 1];
}

void output_dof_space(Index state_index, std::string const &identifier,
                      ConfigEnumInput const &input_state,
                      DoFSpaceAnalysisOptions const &options,
                      OutputImpl &output) {
  Log &log = CASM::log();
  log.begin(identifier);
  log.increase_indent();

  std::vector<PermuteIterator> group = make_invariant_subgroup(input_state);

  if (options.write_symmetry) {
    Lattice config_lattice = input_state.configuration().ideal_lattice();
    SymGroup lattice_point_group{SymGroup::lattice_point_group(config_lattice)};
    SymGroup factor_group =
        make_sym_group(group.begin(), group.end(), config_lattice);
    SymGroup crystal_point_group =
        make_point_group(group.begin(), group.end(), config_lattice);
    output.write_symmetry(state_index, identifier, input_state,
                          lattice_point_group, factor_group,
                          crystal_point_group);
  }

  // if(options.write_config) {
  //   output.write_config(state_index, identifier, config_enum_input);
  // }

  if (options.write_structure) {
    output.write_structure(state_index, identifier, input_state);
  }

  for (DoFKey const &dof : options.dofs) {
    log << "Working on: " << identifier << " " << dof << std::endl;
    DoFSpace dof_space = make_dof_space(dof, input_state);
    std::optional<VectorSpaceSymReport> report;
    try {
      if (options.sym_axes) {
        dof_space = make_symmetry_adapted_dof_space(
            dof_space, input_state, group, options.calc_wedge, report);
      }
    } catch (make_symmetry_adapted_dof_space_error &e) {
      log << "Error: " << e.what() << std::endl;
      log << "skipping: " << identifier << " " << dof << std::endl << std::endl;
      output.write_dof_space_error(e, state_index, dof_space, identifier,
                                   input_state, report);
      continue;
    }
    output.write_dof_space(state_index, dof_space, identifier, input_state,
                           report);
  }

  log.decrease_indent();
  log << std::endl;
}

void dof_space_analysis(
    std::vector<std::pair<std::string, ConfigEnumInput>> const &named_inputs,
    DoFSpaceAnalysisOptions const &options, OutputImpl &output) {
  // For each state, for each DoF type specified, perform analysis and write
  // files.
  Index state_index{0};
  for (auto const &named_input : named_inputs) {
    std::string identifier = named_input.first;
    ConfigEnumInput const &input_state = named_input.second;
    output_dof_space(state_index, identifier, input_state, options, output);

    ++state_index;
  }
  output.write_errors();
}

}  // namespace DoFSpaceIO

}  // namespace CASM
