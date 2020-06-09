#include <boost/filesystem.hpp>
#include "casm/app/AppIO_impl.hh"
#include "casm/app/ClexDescription.hh"
#include "casm/app/DirectoryStructure.hh"
#include "casm/app/ProjectSettings.hh"
#include "casm/app/bset.hh"
#include "casm/clex/Clexulator.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clusterography/ClusterOrbits_impl.hh"
#include "casm/clusterography/ClusterSymCompare_impl.hh"
#include "casm/clusterography/IntegralCluster_impl.hh"
#include "casm/completer/Handlers.hh"

namespace CASM {

  namespace Completer {

    BsetOption::BsetOption(): OptionHandlerBase("bset") {}

    void BsetOption::initialize() {
      add_help_suboption();

      m_desc.add_options()
      ("update,u", "Update basis set")
      ("orbits", "Pretty-print orbit prototypes")
      ("functions", "Pretty-print prototype cluster functions for each orbit")
      ("clusters", "Pretty-print all clusters")
      ("clex", po::value<std::string>()->value_name(ArgHandler::clex()), "Name of the cluster expansion using the basis set")
      ("force,f", "Force overwrite");
      return;
    }
  }


  namespace bset_impl {

    // might be useful for other casm commands...
    template<typename CommandType>
    ClexDescription get_clex_description(const CommandType &cmd) {
      const auto &vm = cmd.vm();
      const auto &primclex = cmd.primclex();

      if(!vm.count("clex")) {
        return primclex.settings().default_clex();
      }
      else {
        const auto &set = primclex.settings();
        auto it = set.cluster_expansions().find(vm["clex"].template as<std::string>());
        if(it == set.cluster_expansions().end()) {
          throw std::runtime_error {"Invalid --clex value: " + vm["clex"].template as<std::string>()};
        }
        return it->second;
      }
    }

    /// Check for generated bset files, print messages when found
    bool any_existing_files(const std::string &bset, const BsetCommand &cmd) {
      const DirectoryStructure &dir = cmd.primclex().dir();
      const ProjectSettings &set = cmd.primclex().settings();
      auto &log = cmd.args().log();

      bool _any_existing_files = false;
      for(const auto &p : dir.bset_data(set.project_name(), bset)) {
        if(fs::exists(p)) {
          if(!_any_existing_files) {
            log.custom("Found existing files");
            _any_existing_files = true;
          }
          log << "found: " << p << "\n";
        }
      }
      return _any_existing_files;
    }

    void check_force(const std::string &bset, const BsetCommand &cmd) {
      const auto &vm = cmd.vm();
      auto &log = cmd.args().log();

      /// check for existing fiels
      if(any_existing_files(bset, cmd)) {
        if(vm.count("force")) {
          log << "Using --force. Will overwrite existing files.\n" << std::endl;
        }
        else {
          throw CASM::runtime_error {"Exiting due to existing files.  Use --force to force overwrite.", ERR_EXISTING_FILE};
        }
      }
    }

    /// Implements `casm bset --update`
    void update_bset(const BsetCommand &cmd) {
      const auto &primclex = cmd.primclex();
      auto basis_set_name = get_clex_description(cmd).bset;

      if(!primclex.has_basis_set_specs(basis_set_name)) {
        auto bspecs_path = primclex.dir().bspecs(basis_set_name);
        throw CASM::runtime_error {"'bspecs.json' file not found at: " + bspecs_path.string(), ERR_MISSING_INPUT_FILE};
      }

      check_force(basis_set_name, cmd);

      jsonParser basis_set_specs;

      try {
        write_clexulator(
          primclex.shared_prim(),
          primclex.settings(),
          basis_set_name,
          primclex.basis_set_specs(basis_set_name),
          primclex.nlist());
      }
      catch(std::exception &e) {
        err_log() << e.what() << std::endl;
        throw CASM::runtime_error {e.what(), ERR_INVALID_INPUT_FILE};
      }

      // force compilation
      auto clexulator = primclex.clexulator(basis_set_name);
    }

    template<typename OrbitVecType>
    void print_bset(const BsetCommand &cmd, const OrbitVecType &orbits) {
      const auto &vm = cmd.vm();
      const auto &primclex = cmd.primclex();
      auto clex_desc = get_clex_description(cmd);
      auto &log = cmd.args().log();

      if(vm.count("orbits")) {
        print_clust(orbits.begin(), orbits.end(), log, ProtoSitesPrinter());
      }
      if(vm.count("clusters")) {
        print_clust(orbits.begin(), orbits.end(), log, FullSitesPrinter());
      }
      if(vm.count("functions")) {
        const auto &shared_prim = primclex.shared_prim();
        const auto &clex_basis = primclex.clex_basis(clex_desc.bset);

        print_site_basis_funcs(shared_prim, clex_basis, log);
        print_clust(
          orbits.begin(),
          orbits.end(),
          log,
          ProtoFuncsPrinter {clex_basis, shared_prim->shared_structure()});
      }
    }

    /// Implements `casm bset --orbits --clusters --functions` (any combination is allowed)
    void print_bset(const BsetCommand &cmd) {
      const auto &primclex = cmd.primclex();
      auto shared_prim = primclex.shared_prim();
      auto xtal_tol = primclex.crystallography_tol();
      auto basis_set_name = get_clex_description(cmd).bset;

      if(!primclex.has_basis_set_specs(basis_set_name)) {
        auto bspecs_path = primclex.dir().bspecs(basis_set_name);
        throw std::runtime_error {"'bspecs.json' file not found at: " + bspecs_path.string()};
      }

      // TODO: update this function with ClusterSpecs

      auto basis_set_specs = primclex.basis_set_specs(basis_set_name);
      CLUSTER_PERIODICITY_TYPE clex_periodicity_type = basis_set_specs.contains("local_bpsecs") ?
                                                       CLUSTER_PERIODICITY_TYPE::LOCAL : CLUSTER_PERIODICITY_TYPE::PRIM_PERIODIC;

      if(clex_periodicity_type == CLUSTER_PERIODICITY_TYPE::PRIM_PERIODIC) {
        std::vector<PrimPeriodicIntegralClusterOrbit> orbits;
        PrimPeriodicSymCompare<IntegralCluster> sym_compare {shared_prim, xtal_tol};
        primclex.orbits(
          basis_set_name,
          std::back_inserter(orbits),
          sym_compare);

        print_bset(cmd, orbits);
      }
      else {
        std::vector<LocalIntegralClusterOrbit> orbits;
        LocalSymCompare<IntegralCluster> sym_compare {shared_prim, xtal_tol};
        primclex.orbits(
          basis_set_name,
          std::back_inserter(orbits),
          sym_compare);

        print_bset(cmd, orbits);
      }
    }

    /// TODO: Move this
    //
    // jsonParser &local_diff_trans_cspecs_to_json(
    //   std::string phenomenal_name,
    //   const PrimClex &primclex,
    //   jsonParser &json) {
    //
    //   const auto &shared_prim = primclex.shared_prim();
    //
    //   const auto &const_db = primclex.db<PrimPeriodicDiffTransOrbit>();
    //   auto dtorbit_it = const_db.find(phenomenal_name);
    //   if(dtorbit_it == const_db.end()) {
    //     throw std::runtime_error("Error: diff_trans '" + phenomenal_name + "' does not exist.");
    //   }
    //   auto prototype = dtorbit_it->prototype();
    //
    //   typedef PrimPeriodicSymCompare<Kinetics::DiffusionTransformation> DTSymCompareType;
    //   DTSymCompareType dt_sym_compare {shared_prim, shared_prim->lattice().tol()};
    //
    //   SymGroup generating_group = make_invariant_subgroup(
    //                                 prototype, shared_prim->factor_group(), dt_sym_compare);
    //
    //   to_json(generating_group, json["generating_group"]);
    //   to_json(prototype.cluster(), json["phenomenal"]);
    //   return json;
    // }

  }


  const std::string BsetCommand::name = "bset";

  BsetCommand::BsetCommand(const CommandArgs &_args, Completer::BsetOption &_opt) :
    APICommand<Completer::BsetOption>(_args, _opt) {}

  BsetCommand::~BsetCommand() {}

  int BsetCommand::vm_count_check() const {
    if(!in_project()) {
      err_log().error("No casm project found");
      err_log() << std::endl;
      return ERR_NO_PROJ;
    }

    if(vm().count("update")) {
      return 0;
    }
    else if(vm().count("orbits") || vm().count("clusters") || vm().count("functions")) {
      return 0;
    }
    else {
      err_log() << "Error in 'casm bset'. Choose --update or one or more of --orbits, --clusters, ----functions." << std::endl;
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
    log().indent() << "Generate and inspect cluster basis functions. A bspecs.json file should be available at\n";

    log().increase_indent();
    log().indent() << "$ROOT/basis_set/$current_bset/bspecs.json\n";
    log().decrease_indent();

    log().indent() << std::endl;
    log().indent() << "For a description of the bspecs.json file use one of:" << std::endl;
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
      if(vm.count("update")) {
        bset_impl::update_bset(*this);
      }
      else if(vm.count("orbits") || vm.count("clusters") || vm.count("functions")) {
        bset_impl::print_bset(*this);
      }
      else {
        err_log() << "Error in 'casm bset'. Choose --update or one or more of --orbits, --clusters, ----functions." << std::endl;
        return ERR_INVALID_ARG;
      }
    }
    catch(CASM::runtime_error &e) {
      err_log() << e.what() << std::endl;
      return e.code();
    }
    catch(std::exception &e) {
      err_log() << e.what() << std::endl;
      return ERR_UNKNOWN;
    }

    return 0;
  }

}
