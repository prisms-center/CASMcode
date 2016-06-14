#include<cstring>

#include "casm/CASM_global_definitions.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/app/DirectoryStructure.hh"
#include "casm/app/casm_functions.hh"
#include "casm/app/AppIO.hh"
#include "casm/completer/Handlers.hh"

namespace CASM {

  void display(std::ostream &sout, const CompositionAxes &opt) {

    if(opt.standard.size()) {
      std::cout << "Standard composition axes:\n\n";
      display_composition_axes(sout, opt.standard);
    }

    if(opt.custom.size()) {
      std::cout << "Custom composition axes:\n\n";
      display_composition_axes(sout, opt.custom);
    }

    std::cout << "\n";
    if(opt.has_current_axes) {
      std::cout << "Currently selected composition axes: " << opt.curr_key << "\n";

      std::cout << "\n";
      std::cout << "Parametric composition:\n";
      display_comp(sout, opt.curr, 2);

      std::cout << "\n";
      std::cout << "Composition:\n";
      display_comp_n(sout, opt.curr, 2);

      std::cout << "\n";
      std::cout << "Parametric chemical potentials:\n";
      display_param_chem_pot(sout, opt.curr, 2);
    }


  }

  namespace Completer {

    CompositionOption::CompositionOption(): OptionHandlerBase("composition") {};

    void CompositionOption::initialize() {
      add_help_suboption();
      m_desc.add_options()
      ("display,d", "Display composition axes file")
      ("calc,c", "Calculate the standard composition axes")
      ("select,s", po::value<std::string>(&m_axis_choice_str), "Select composition axes");

      return;
    }

    const std::string &CompositionOption::axis_choice_str() const {
      return m_axis_choice_str;
    }


  }

  // ///////////////////////////////////////
  // 'composition' function for casm
  //    (add an 'if-else' statement in casm.cpp to call this)

  int composition_command(const CommandArgs &args) {
    po::variables_map vm;

    Completer::CompositionOption comp_opt;
    //("update,u", "Update composition and references based on current 'composition_axes.json' file");

    try {
      po::store(po::parse_command_line(args.argc, args.argv, comp_opt.desc()), vm);

      bool call_help = false;

      //quit out if there are no arguments
      if(!vm.count("help")) {
        if(vm.count("calc") + vm.count("select") + vm.count("display") + vm.count("update") != 1) {
          std::cout << "Error in 'casm composition'. You need to use either --calc, --select, --display, or --update." << std::endl;
          call_help = true;
        }
      }

      /** --help option
       */
      if(vm.count("help") || call_help) {
        std::cout << std::endl;
        std::cout << comp_opt.desc() << std::endl;

        std::cout << "DESCRIPTION" << std::endl;
        std::cout << "    Setup the composition axes.\n";
        std::cout << "    - expects a PRIM file in the project root directory \n";
        std::cout << "    - custom composition axes can be added to the 'composition_axes.json' file \n";
        std::cout << "      and then compositions and references updated with 'casm composition --update'\n";
        std::cout << "      (word of caution: compatibility with PRIM is currently not checked for.)\n";
        std::cout << " \n";
        std::cout << "    Examples:\n";
        std::cout << "      casm composition --calc \n";
        std::cout << "      - Calculate standard composition axes, but does not select any. \n";
        std::cout << "\n";
        std::cout << "      casm composition --display \n";
        std::cout << "      - Display the possible composition axes\n";
        std::cout << "      - Possible axes are stored in the 'composition_axes.json' file.\n";
        std::cout << "\n";
        std::cout << "      casm composition --select 0 \n";
        std::cout << "      - Select the composition axes by key name.\n";
        std::cout << "      - Updates and writes compositions and reference properties.\n";
        std::cout << "\n";
        /*std::cout << "      casm composition --update \n";
        std::cout << "      - Updates and writes compositions and reference properties based.\n";
        std::cout << "        on the contents of the 'composition_axes.json' file.\n";
        std::cout << "\n";*/
        if(call_help)
          return ERR_INVALID_ARG;

        return 0;
      }

      po::notify(vm);


    }
    catch(po::error &e) {
      std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
      std::cerr << comp_opt.desc() << std::endl;
      return ERR_INVALID_ARG;
    }
    catch(std::exception &e) {
      std::cerr << "Unhandled Exception reached the top of main: "
                << e.what() << ", application will now exit" << std::endl;
      return ERR_UNKNOWN;

    }

    const fs::path &root = args.root;
    if(root.empty()) {
      args.err_log.error("No casm project found");
      args.err_log << std::endl;
      return ERR_NO_PROJ;
    }

    // If 'args.primclex', use that, else construct PrimClex in 'uniq_primclex'
    // Then whichever exists, store reference in 'primclex'
    std::unique_ptr<PrimClex> uniq_primclex;
    PrimClex &primclex = make_primclex_if_not(args, uniq_primclex);

    const DirectoryStructure &dir = primclex.dir();
    std::string calctype = primclex.settings().calctype();
    std::string ref = primclex.settings().ref();

    CompositionAxes opt;
    if(fs::exists(dir.composition_axes(calctype, ref))) {
      opt.read(dir.composition_axes(calctype, ref));
    }

    if(vm.count("display")) {

      std::cout << "\n***************************\n\n";

      display(std::cout, opt);

      std::cout << "\n\n";

      if(opt.err_code) {
        std::cout << opt.err_message << std::endl;
      }
      else if(opt.standard.size() && !opt.has_current_axes) {
        std::cout << "Please use 'casm composition --select' to choose your composition axes.\n\n";
      }
      else if(!opt.standard.size() && !opt.has_current_axes) {
        std::cout << "Please use 'casm composition --calc' to calculate standard composition axes.\n\n";
      }

      return 0;
    }
    else {

      if(vm.count("calc")) {

        std::cout << "\n***************************\n\n";

        std::cout << "Using the PRIM to enumerate standard composition axes for this space.\n\n";

        if(opt.standard.size()) {
          std::cout << "Overwriting existing standard composition axes.\n\n";
        }

        opt.standard.clear();
        std::vector<CompositionConverter> v;
        standard_composition_axes(primclex.get_prim(), std::back_inserter(v));
        for(int i = 0; i < v.size(); i++) {
          opt.standard[std::to_string(i)] = v[i];
        }

        display(std::cout, opt);

        std::cout << "\n\n";


        if(!fs::exists(dir.ref_dir(calctype, ref))) {
          fs::create_directories(dir.ref_dir(calctype, ref));
        }

        opt.write(dir.composition_axes(calctype, ref));

        std::cout << "Wrote: " << dir.composition_axes(calctype, ref) << "\n\n";

        if(!opt.has_current_axes) {
          std::cout << "Please use 'casm composition --select' to choose your composition axes.\n\n";
        }
      }

      if(vm.count("select")) {

        std::cout << "\n***************************\n\n";

        if(opt.standard.size() + opt.custom.size() == 0) {

          std::cout << "Error: No composition axes found.\n\n";

          std::cout << "Please use 'casm composition --calc' to calculate standard composition axes,\n" <<
                    "or add custom composition axes to " << dir.composition_axes(calctype, ref) << "\n\n";

          return ERR_MISSING_DEPENDS;

        }

        if(opt.standard.find(comp_opt.axis_choice_str()) != opt.standard.end() &&
           opt.custom.find(comp_opt.axis_choice_str()) != opt.custom.end()) {

          std::cout << "Error: The selected composition axes '" << comp_opt.axis_choice_str() << "' can be \n" <<
                    "found in both the standard and custom compostion axes. Please  \n" <<
                    "edit the custom composition axes to remove this ambiguity.     \n\n" <<

                    "File: " << dir.composition_axes(calctype, ref) << "\n\n";

          return ERR_INVALID_INPUT_FILE;

        }
        else if(opt.standard.find(comp_opt.axis_choice_str()) == opt.standard.end() &&
                opt.custom.find(comp_opt.axis_choice_str()) == opt.custom.end()) {

          std::cout << "Error: The selected composition axes '" << comp_opt.axis_choice_str() << "' can not \n" <<
                    "be found in either the standard or custom compostion axes. Please\n" <<
                    "re-select composition axes.                                     \n\n";

          return ERR_INVALID_INPUT_FILE;

        }

        opt.select(comp_opt.axis_choice_str());

        display(std::cout, opt);

        std::cout << "\n\n";

        opt.write(dir.composition_axes(calctype, ref));

        std::cout << "Wrote: " << dir.composition_axes(calctype, ref) << "\n\n";

      }

      return 0;
    }
  }

}

