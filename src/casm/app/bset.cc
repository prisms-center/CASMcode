#include "casm/app/bset.hh"

#include <boost/filesystem.hpp>

#include "casm/app/AppIO_impl.hh"
#include "casm/app/ClexDescription.hh"
#include "casm/app/DirectoryStructure.hh"
#include "casm/app/ProjectSettings.hh"
#include "casm/basis_set/DoFTraits.hh"
#include "casm/casm_io/Log.hh"
#include "casm/clex/ClexBasisSpecs.hh"
#include "casm/clex/ClexBasis_impl.hh"
#include "casm/clex/Clexulator.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clusterography/ClusterSpecs_impl.hh"
#include "casm/completer/Handlers.hh"

namespace CASM {

namespace Completer {

BsetOption::BsetOption() : OptionHandlerBase("bset") {}

void BsetOption::initialize() {
  add_help_suboption();

  m_desc.add_options()("update,u", "Update basis set")(
      "orbits", "Pretty-print orbit prototypes")(
      "functions", "Pretty-print prototype cluster functions for each orbit")(
      "clusters", "Pretty-print all clusters")(
      "clex", po::value<std::string>()->value_name(ArgHandler::clex()),
      "Name of the cluster expansion using the basis set")("force,f",
                                                           "Force overwrite");
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
  const auto &primclex = cmd.primclex();

  auto basis_set_name = get_clex_description(cmd).bset;

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

  // force compilation
  auto clexulator = primclex.clexulator(basis_set_name);
}

/// This functor implements a template method so that clusters can be printed
/// for any orbit type, as determined at runtime from the JSON file parameters
template <typename PrinterType>
class OrbitPrinterAdapter {
 public:
  OrbitPrinterAdapter(Log &_log);

  template <typename OrbitVecType>
  void operator()(OrbitVecType const &orbits) const;

 private:
  PrinterType m_printer;
  Log &m_log;
};

/// This functor implements a template method so that basis functions can be
/// printed for any orbit type, as determined at runtime from the JSON file
/// parameters
class BasisFunctionPrinter {
 public:
  BasisFunctionPrinter(Log &_log, std::shared_ptr<Structure const> _shared_prim,
                       ClexBasisSpecs const &_basis_set_specs);

  template <typename OrbitVecType>
  void operator()(OrbitVecType const &orbits) const;

 private:
  std::shared_ptr<Structure const> m_shared_prim;
  ClexBasisSpecs m_basis_set_specs;
  Log &m_log;
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

  if (vm.count("orbits")) {
    for_all_orbits(cluster_specs, prototypes,
                   OrbitPrinterAdapter<ProtoSitesPrinter>{log});
  }
  if (vm.count("clusters")) {
    for_all_orbits(cluster_specs, prototypes,
                   OrbitPrinterAdapter<FullSitesPrinter>{log});
  }
  if (vm.count("functions")) {
    BasisFunctionPrinter printer{log, shared_prim, basis_set_specs};
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
OrbitPrinterAdapter<PrinterType>::OrbitPrinterAdapter(Log &_log)
    : m_log(_log) {}

template <typename PrinterType>
template <typename OrbitVecType>
void OrbitPrinterAdapter<PrinterType>::operator()(
    OrbitVecType const &orbits) const {
  print_clust(orbits.begin(), orbits.end(), m_log, m_printer);
}

BasisFunctionPrinter::BasisFunctionPrinter(
    Log &_log, std::shared_ptr<Structure const> _shared_prim,
    ClexBasisSpecs const &_basis_set_specs)
    : m_shared_prim(_shared_prim),
      m_basis_set_specs(_basis_set_specs),
      m_log(_log) {}

template <typename OrbitVecType>
void BasisFunctionPrinter::operator()(OrbitVecType const &orbits) const {
  ParsingDictionary<DoFType::Traits> const *dof_dict = &DoFType::traits_dict();
  ClexBasis clex_basis{m_shared_prim, m_basis_set_specs, dof_dict};
  clex_basis.generate(orbits.begin(), orbits.end());

  print_site_basis_funcs(m_shared_prim, clex_basis, m_log);
  ProtoFuncsPrinter funcs_printer{clex_basis,
                                  m_shared_prim->shared_structure()};
  print_clust(orbits.begin(), orbits.end(), m_log, funcs_printer);
}
}  // namespace bset_impl
}  // namespace CASM
