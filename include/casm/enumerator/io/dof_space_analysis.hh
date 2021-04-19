#ifndef CASM_enumerator_io_json_dof_space_analysis
#define CASM_enumerator_io_json_dof_space_analysis

#include <boost/filesystem/path.hpp>

#include "casm/casm_io/json/jsonParser.hh"
#include "casm/crystallography/DoFDecl.hh"
#include "casm/global/definitions.hh"
#include "casm/symmetry/PermuteIterator.hh"

namespace CASM {

class ConfigEnumInput;
class DirectoryStructure;
class DoFSpace;
class SymGroup;
class make_symmetry_adapted_dof_space_error;

namespace SymRepTools_v2 {
struct VectorSpaceSymReport;
}

namespace DoFSpaceIO {

class OutputImpl {
 public:
  /// Provide state_index, indentifier, and dof_key for any errors
  struct Error {
    Error(Index _state_index, std::string _identifier, DoFKey _dof_key,
          std::string _what, jsonParser _data);
    Index state_index;
    std::string identifier;
    DoFKey dof_key;
    std::string what;
    jsonParser data;
  };

  /// Write symmetry groups (lattice point group, factor_group,
  /// crystal_point_group)
  void write_symmetry(Index state_index, std::string const &identifier,
                      ConfigEnumInput const &config_enum_input,
                      std::vector<PermuteIterator> const &group);

  /// Write symmetry groups (lattice point group, factor_group,
  /// crystal_point_group)
  ///
  /// - lattice_point_group: point group of supercell lattice
  /// - factor_group: subgroup of supercell factor group that keeps
  /// configuration DoF invariant
  ///   and selected sites invariant
  /// - crystal_point_group: factor_group excluding translations
  virtual void write_symmetry(Index state_index, std::string const &identifier,
                              ConfigEnumInput const &config_enum_input,
                              SymGroup const &lattice_point_group,
                              SymGroup const &factor_group,
                              SymGroup const &crystal_point_group) = 0;

  // /// Write input state configuration
  // virtual void write_config(
  //   Index state_index,
  //   std::string const &identifier,
  //   ConfigEnumInput const &config_enum_input) = 0;

  /// Write input state structure
  virtual void write_structure(Index state_index, std::string const &identifier,
                               ConfigEnumInput const &config_enum_input) = 0;

  /// Write dof space analysis
  virtual void write_dof_space(
      Index state_index, DoFSpace const &dof_space,
      std::string const &identifier, ConfigEnumInput const &config_enum_input,
      std::optional<SymRepTools_v2::VectorSpaceSymReport> const
          &sym_report) = 0;

  /// Write dof space analysis error information
  void write_dof_space_error(
      make_symmetry_adapted_dof_space_error const &e, Index state_index,
      DoFSpace const &dof_space, std::string const &identifier,
      ConfigEnumInput const &config_enum_input,
      std::optional<SymRepTools_v2::VectorSpaceSymReport> const &sym_report);

  std::vector<Error> const &errors() const { return m_errors; }
  std::vector<Error> &errors() { return m_errors; }

  /// Write dof space analysis error information to
  /// <current_path>/dof_space_errors.json
  void write_errors() const;

  /// If true, output current status to log
  virtual bool output_status() const { return true; }

 private:
  std::vector<Error> m_errors;
};

class DirectoryOutput : public OutputImpl {
 public:
  using OutputImpl::write_symmetry;

  void write_symmetry(Index state_index, std::string const &identifier,
                      ConfigEnumInput const &config_enum_input,
                      SymGroup const &lattice_point_group,
                      SymGroup const &factor_group,
                      SymGroup const &crystal_point_group) override;

  // void write_config(
  //   Index state_index,
  //   std::string const &identifier,
  //   ConfigEnumInput const &config_enum_input) override;

  void write_structure(Index state_index, std::string const &identifier,
                       ConfigEnumInput const &config_enum_input) override;

  void write_dof_space(Index state_index, DoFSpace const &dof_space,
                       std::string const &identifier,
                       ConfigEnumInput const &config_enum_input,
                       std::optional<SymRepTools_v2::VectorSpaceSymReport> const
                           &sym_report) override;

 private:
  /// For SymmetryDirectoryOutput, configurations must exist in database and
  /// this will throw otherwise
  virtual void _check_config(Index state_index, std::string const &identifier,
                             ConfigEnumInput const &config_enum_input) = 0;

  /// For SymmetryDirectoryOutput vs SequentialDirectoryOutput return output
  /// directory path
  virtual fs::path _output_dir(Index state_index, std::string const &identifier,
                               ConfigEnumInput const &config_enum_input) = 0;
};

/// Implementation that outputs to <casm_project>/symmetry/analysis/<configname>
class SymmetryDirectoryOutput : public DirectoryOutput {
 public:
  SymmetryDirectoryOutput(DirectoryStructure const &dir);

 private:
  /// For SymmetryDirectoryOutput, configurations must exist in database and
  /// this will throw otherwise
  void _check_config(Index state_index, std::string const &identifier,
                     ConfigEnumInput const &config_enum_input) override;

  /// For SymmetryDirectoryOutput return output directory path in symmetry dir
  fs::path _output_dir(Index state_index, std::string const &identifier,
                       ConfigEnumInput const &config_enum_input) override;

  DirectoryStructure const &m_dir;
};

/// Implementation that outputs to <output_dir>/dof_space/state.<index>
class SequentialDirectoryOutput : public DirectoryOutput {
 public:
  SequentialDirectoryOutput(fs::path output_dir);

 private:
  /// For SequentialDirectoryOutput, any input state is allowed
  void _check_config(Index state_index, std::string const &identifier,
                     ConfigEnumInput const &config_enum_input) override;

  /// For SequentialDirectoryOutput return output directory path
  fs::path _output_dir(Index state_index, std::string const &identifier,
                       ConfigEnumInput const &config_enum_input) override;

  fs::path m_output_dir;
};

/// Implementation that outputs all results to one JSON array
class CombinedJsonOutput : public OutputImpl {
 public:
  /// Does not write to file on destruction
  CombinedJsonOutput();

  /// Writes to <output_dir>/dof_space.json on destruction
  CombinedJsonOutput(fs::path output_dir);

  ~CombinedJsonOutput();

  using OutputImpl::write_symmetry;

  void write_symmetry(Index state_index, std::string const &identifier,
                      ConfigEnumInput const &config_enum_input,
                      SymGroup const &lattice_point_group,
                      SymGroup const &factor_group,
                      SymGroup const &crystal_point_group) override;

  // void write_config(
  //   Index state_index,
  //   std::string const &identifier,
  //   ConfigEnumInput const &config_enum_input) override;

  void write_structure(Index state_index, std::string const &identifier,
                       ConfigEnumInput const &config_enum_input) override;

  void write_dof_space(Index state_index, DoFSpace const &dof_space,
                       std::string const &identifier,
                       ConfigEnumInput const &config_enum_input,
                       std::optional<SymRepTools_v2::VectorSpaceSymReport> const
                           &sym_report) override;

  jsonParser const &combined_json() const { return m_combined_json; }

  /// Only output status if writing to file
  bool output_status() const override { return !m_output_dir.empty(); }

 private:
  jsonParser &_output_json(Index state_index);

  jsonParser m_combined_json;

  fs::path m_output_dir;
};

struct DoFSpaceAnalysisOptions {
  std::vector<std::string> dofs;
  bool sym_axes = true;
  bool calc_wedge = false;
  bool write_symmetry = true;
  bool write_structure = true;
};

void output_dof_space(Index state_index, std::string const &identifier,
                      ConfigEnumInput const &input_state,
                      DoFSpaceAnalysisOptions const &options,
                      OutputImpl &output);

void dof_space_analysis(
    std::vector<std::pair<std::string, ConfigEnumInput>> const &named_inputs,
    DoFSpaceAnalysisOptions const &options, OutputImpl &output);
}  // namespace DoFSpaceIO

}  // namespace CASM

#endif
