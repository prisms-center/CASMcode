#include "casm/app/bset.hh"

#include <boost/filesystem.hpp>

#include "casm/app/ClexDescription.hh"
#include "casm/app/DirectoryStructure.hh"
#include "casm/app/ProjectSettings.hh"
#include "casm/basis_set/DoFTraits.hh"
#include "casm/casm_io/Log.hh"
#include "casm/clex/ClexBasisSpecs.hh"
#include "casm/clex/ClexBasis_impl.hh"
#include "casm/clex/Clexulator.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/io/stream/ClexBasis_stream_io_impl.hh"
#include "casm/clusterography/ClusterSpecs_impl.hh"
#include "casm/clusterography/io/OrbitPrinter_impl.hh"
#include "casm/completer/Handlers.hh"

namespace CASM {

namespace Completer {

BsetOption::BsetOption() : OptionHandlerBase("bset") {}

void BsetOption::initialize() {
  add_help_suboption();
  add_coordtype_suboption();

  m_desc.add_options()(
      //
      "update,u", "Update basis set")(
      //
      "orbits", "Pretty-print orbit prototypes")(
      //
      "functions", "Pretty-print prototype cluster functions for each orbit")(
      //
      "clusters", "Pretty-print all clusters")(
      //
      "clex", po::value<std::string>()->value_name(ArgHandler::clex()),
      "Select a non-default basis set to print by specifying the name of a "
      "cluster expansion using the basis set")(
      //
      "print-invariant-group",
      "Use with --orbits or --clusters to print the group that leaves clusters "
      "invariant.")(
      //
      "print-equivalence-map",
      "Use with --clusters to print the orbit equivalence map (sets of "
      "symmetry operations that map the prototype cluster onto an equivalent "
      "cluster / cosets of the invariant group).")(
      //
      "align",
      "Use with --functions to print aligned functions ready to copy and paste "
      "into a latex document")(
      //
      "no-compile",
      "Use with --update to generate basis.json, clust.json, and "
      "the basis function source code but not compile it.")(
      //
      "only-compile",
      "Use with --update to keep the existing basis.json, clust.json, and "
      "basis function source code but to re-compile it if the .so and .o files "
      "do not exist.")(
      //
      "json",
      "Use with --orbits, --clusters, or --functions to print using "
      "JSON format.")(
      //
      "force,f", "Force overwrite");
  return;
}
}  // namespace Completer

namespace bset_impl {

// might be useful for other casm commands...
template <typename CommandType>
ClexDescription get_clex_description(const CommandType &cmd) {
  const auto &vm = cmd.vm();
  const auto &primclex = cmd.primclex();

  if (!vm.count("clex")) {
    return primclex.settings().default_clex();
  } else {
    const auto &set = primclex.settings();
    auto it =
        set.cluster_expansions().find(vm["clex"].template as<std::string>());
    if (it == set.cluster_expansions().end()) {
      throw std::runtime_error{"Invalid --clex value: " +
                               vm["clex"].template as<std::string>()};
    }
    return it->second;
  }
}

/// Check for generated bset files, print messages when found to CASM::log()
bool any_existing_files(const std::string &bset, const BsetCommand &cmd) {
  const DirectoryStructure &dir = cmd.primclex().dir();
  const ProjectSettings &set = cmd.primclex().settings();

  bool _any_existing_files = false;
  for (const auto &p : dir.bset_data(set.project_name(), bset)) {
    if (fs::exists(p)) {
      if (!_any_existing_files) {
        log().custom("Found existing files");
        _any_existing_files = true;
      }
      log() << "found: " << p << "\n";
    }
  }
  return _any_existing_files;
}

void check_force(const std::string &bset, const BsetCommand &cmd) {
  const auto &vm = cmd.vm();

  /// check for existing fiels
  if (any_existing_files(bset, cmd)) {
    if (vm.count("force")) {
      log() << "Using --force. Will overwrite existing files.\n" << std::endl;
    } else {
      throw CASM::runtime_error{
          "Exiting due to existing files.  Use --force to force overwrite.",
          ERR_EXISTING_FILE};
    }
  }
}

/// Implements `casm bset --update`
void update_bset(const BsetCommand &cmd) {
  const auto &vm = cmd.vm();
  const auto &primclex = cmd.primclex();

  auto basis_set_name = get_clex_description(cmd).bset;

  // force compilation by default, opt out with --no-compile, or only compile
  // with --only-compile
  if (vm.count("only-compile")) {
    auto clexulator = primclex.clexulator(basis_set_name);
    return;
  }

  if (!primclex.has_basis_set_specs(basis_set_name)) {
    auto bspecs_path = primclex.dir().bspecs(basis_set_name);
    throw CASM::runtime_error{
        "'bspecs.json' file not found at: " + bspecs_path.string(),
        ERR_MISSING_INPUT_FILE};
  }

  check_force(basis_set_name, cmd);

  try {
    write_basis_set_data(
        primclex.shared_prim(), primclex.settings(), basis_set_name,
        primclex.basis_set_specs(basis_set_name), primclex.nlist());
  } catch (std::exception &e) {
    throw CASM::runtime_error{e.what(), ERR_INVALID_INPUT_FILE};
  }

  // force compilation by default, opt out with --no-compile
  if (!vm.count("no-compile")) {
    auto clexulator = primclex.clexulator(basis_set_name);
  }
}

/// This functor implements a template method so that clusters can be printed
/// for any orbit type, as determined at runtime from the JSON file parameters
template <typename PrinterType>
class OrbitPrinterAdapter {
 public:
  OrbitPrinterAdapter(Log &_log, OrbitPrinterOptions _opt, bool _json);

  template <typename OrbitVecType>
  void operator()(OrbitVecType const &orbits) const;

 private:
  Log &m_log;
  PrinterType m_printer;
  bool m_json;
};

/// Implements `casm bset --orbits --clusters --functions` (any combination is
/// allowed)
void print_bset(const BsetCommand &cmd) {
  auto const &vm = cmd.vm();
  Log &log = CASM::log();
  std::string basis_set_name = get_clex_description(cmd).bset;
  auto const &primclex = cmd.primclex();

  auto const &shared_prim = primclex.shared_prim();
  auto const &basis_set_specs = primclex.basis_set_specs(basis_set_name);
  auto const &cluster_specs = *basis_set_specs.cluster_specs;

  std::vector<IntegralCluster> prototypes;
  jsonParser clust_json{primclex.dir().clust(basis_set_name)};
  read_clust(std::back_inserter(prototypes), clust_json, *shared_prim);

  OrbitPrinterOptions orbit_printer_options;
  orbit_printer_options.coord_type = cmd.opt().coordtype_enum();
  orbit_printer_options.indent_space = 4;
  if (orbit_printer_options.coord_type != INTEGRAL) {
    orbit_printer_options.delim = 0;
  }
  if (vm.count("print-invariant-group")) {
    orbit_printer_options.print_invariant_group = true;
  }
  if (vm.count("print-equivalence-map")) {
    orbit_printer_options.print_equivalence_map = true;
  }

  bool output_json = vm.count("json");

  if (vm.count("orbits")) {
    OrbitPrinterAdapter<ProtoSitesPrinter> printer{log, orbit_printer_options,
                                                   output_json};
    for_all_orbits(cluster_specs, prototypes, printer);
  }
  if (vm.count("clusters")) {
    OrbitPrinterAdapter<FullSitesPrinter> printer{log, orbit_printer_options,
                                                  output_json};
    for_all_orbits(cluster_specs, prototypes, printer);
  }
  if (vm.count("functions")) {
    bool align = vm.count("align");
    if (align) {
      orbit_printer_options.itemize_orbits = true;
    }
    ClexBasisFunctionPrinter printer{
        log,   shared_prim,           basis_set_specs,
        align, orbit_printer_options, output_json};
    for_all_orbits(cluster_specs, prototypes, printer);
  }
}
}  // namespace bset_impl

const std::string BsetCommand::name = "bset";

BsetCommand::BsetCommand(const CommandArgs &_args, Completer::BsetOption &_opt)
    : APICommand<Completer::BsetOption>(_args, _opt) {}

BsetCommand::~BsetCommand() {}

int BsetCommand::vm_count_check() const {
  if (!in_project()) {
    err_log().error("No casm project found");
    err_log() << std::endl;
    return ERR_NO_PROJ;
  }

  if (vm().count("update")) {
    return 0;
  } else if (vm().count("orbits") || vm().count("clusters") ||
             vm().count("functions")) {
    return 0;
  } else {
    err_log() << "Error in 'casm bset'. Choose --update or one or more of "
                 "--orbits, --clusters, ----functions."
              << std::endl;
    return ERR_INVALID_ARG;
  }

  return 0;
}

int BsetCommand::help() const {
  log() << opt().desc() << std::endl;
  return 0;
}

int BsetCommand::desc() const {
  log() << "\n";
  log() << opt().desc() << std::endl;

  log().indent() << "DESCRIPTION\n" << std::endl;
  log().increase_indent();
  log().indent() << "Generate and inspect cluster basis functions. A "
                    "bspecs.json file should be available at\n";

  log().increase_indent();
  log().indent() << "$ROOT/basis_set/$current_bset/bspecs.json\n";
  log().decrease_indent();

  log().indent() << std::endl;
  log().indent() << "For a description of the bspecs.json file use one of:"
                 << std::endl;
  log().increase_indent();
  log() << "'casm format --bpsecs'" << std::endl;
  log() << "'casm format --local_bpsecs'" << std::endl;
  log().decrease_indent();

  log().decrease_indent();
  log() << std::endl;

  return 0;
}

int BsetCommand::run() const {
  const auto &vm = this->vm();

  try {
    if (vm.count("update")) {
      bset_impl::update_bset(*this);
    } else if (vm.count("orbits") || vm.count("clusters") ||
               vm.count("functions")) {
      bset_impl::print_bset(*this);
    } else {
      err_log() << "Error in 'casm bset'. Choose --update or one or more of "
                   "--orbits, --clusters, ----functions."
                << std::endl;
      return ERR_INVALID_ARG;
    }
  } catch (CASM::runtime_error &e) {
    err_log() << e.what() << std::endl;
    return e.code();
  } catch (std::exception &e) {
    err_log() << e.what() << std::endl;
    return ERR_UNKNOWN;
  }

  return 0;
}

namespace bset_impl {

template <typename PrinterType>
OrbitPrinterAdapter<PrinterType>::OrbitPrinterAdapter(Log &_log,
                                                      OrbitPrinterOptions _opt,
                                                      bool _json)
    : m_log(_log), m_printer(_opt), m_json(_json) {}

template <typename PrinterType>
template <typename OrbitVecType>
void OrbitPrinterAdapter<PrinterType>::operator()(
    OrbitVecType const &orbits) const {
  if (m_json) {
    jsonParser json;
    write_clust(orbits.begin(), orbits.end(), json, m_printer);
    m_log << json << std::endl;
  } else {
    print_clust(orbits.begin(), orbits.end(), m_log, m_printer);
  }
}

}  // namespace bset_impl
}  // namespace CASM
