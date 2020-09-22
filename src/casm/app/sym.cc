#include "casm/app/sym.hh"
#include "casm/app/io/json_io_impl.hh"
#include "casm/app/sym/dof_space_analysis.hh"
#include "casm/app/sym/json_io.hh"
#include "casm/app/sym/symmetrize.hh"
#include "casm/app/sym/write_prim_symmetry.hh"
#include "casm/clex/PrimClex.hh"

namespace CASM {

  namespace Completer {
    SymOption::SymOption(): EnumOptionBase("sym") {}

    void SymOption::initialize() {
      bool required = false;

      add_help_suboption();
      add_coordtype_suboption();
      add_selection_suboption("NONE");
      add_confignames_suboption();
      //add_scelnames_suboption();
      add_dofs_suboption();
      add_settings_suboption(required);
      add_input_suboption(required);

      m_desc.add_options()
      ("lattice-point-group", "Pretty print lattice point group")
      ("factor-group", "Pretty print factor group")
      ("crystal-point-group", "Pretty print crystal point group")
      ("calc-wedge", "Perform calculation of irreducible wedge (may significantly slow down analysis)")
      //("no-directions", "Skip calculation of high-symmetry direction and irreducible wedge (for faster evaluation)")
      ("tol", po::value<double>(&m_tol)->default_value(1.0e-5), "Tolerance (in Angstr.) used for symmetrization (default 1e-5)")
      ("symmetrize", po::value<fs::path>(&m_poscar_path)->value_name(ArgHandler::path()), "symmetrize a POSCAR specified by path to a given tolerance");

      return;
    }
  }
}

namespace CASM {
  // ///////////////////////////////////////
  // 'sym' function for casm
  //    (add an 'if-else' statement in casm.cpp to call this)

  const std::string SymCommand::name = "sym";

  SymCommand::SymCommand(const CommandArgs &_args, Completer::SymOption &_opt) :
    APICommand<Completer::SymOption>(_args, _opt) {}


  int SymCommand::vm_count_check() const {
    if(!vm().count("symmetrize") && !in_project()) {
      help();
      err_log().error("No casm project found");
      err_log() << std::endl;
      return ERR_NO_PROJ;
    }

    if(vm().count("settings") + vm().count("input") == 2) {
      help();
      err_log() << "Error in 'casm sym'. The options --settings or --input may not both be chosen." << std::endl;
      return ERR_INVALID_ARG;
    }

    return 0;
  }

  int SymCommand::help() const {
    log() << "\n";
    log() << opt().desc() << std::endl;

    return 0;
  }

  int SymCommand::desc() const {
    log() << "\n";
    log() << opt().desc() << std::endl;
    log() << "DESCRIPTION" << std::endl;
    log() << "    Display symmetry group information.\n";

    log() << write_prim_symmetry_desc();
    log() << symmetrize_desc();
    log() << dof_space_analysis_desc();
    return 0;
  }

  int SymCommand::run() const {

    jsonParser json_options = make_json_input(opt());

    jsonParser cli_options_as_json;
    to_json(opt(), cli_options_as_json);

    if(vm().count("symmetrize")) {
      symmetrize(*this, json_options, cli_options_as_json);
    }
    else if(vm().count("dof-space-analysis")) {
      dof_space_analysis(*this, json_options, cli_options_as_json);
    }
    else {
      write_prim_symmetry(*this, json_options, cli_options_as_json);
    }
    return 0;
  }

}
