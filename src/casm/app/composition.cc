#include "casm/app/composition.hh"

#include <boost/filesystem.hpp>
#include <cstring>

#include "casm/app/DirectoryStructure.hh"
#include "casm/app/ProjectSettings.hh"
#include "casm/app/casm_functions.hh"
#include "casm/casm_io/Log.hh"
#include "casm/clex/CompositionAxes_impl.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/io/file/CompositionAxes_file_io.hh"
#include "casm/completer/Handlers.hh"
#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/global/definitions.hh"

namespace CASM {

void display(std::ostream &sout, const CompositionAxes &comp_axes) {
  if (comp_axes.all_axes.size()) {
    sout << "Possible composition axes:\n\n";
    display_composition_axes(sout, comp_axes.all_axes);
  }
  sout << "\n";

  if (comp_axes.has_current_axes()) {
    sout << "Currently selected composition axes: " << comp_axes.curr_key
         << "\n";

    sout << "\n";
    sout << "Parametric composition:\n";
    display_comp(sout, comp_axes.curr, 2);

    sout << "\n";
    sout << "Composition:\n";
    display_comp_n(sout, comp_axes.curr, 2);

    sout << "\n";
    sout << "Parametric chemical potentials:\n";
    display_param_chem_pot(sout, comp_axes.curr, 2);
  }
}

namespace Completer {

CompositionOption::CompositionOption() : OptionHandlerBase("composition"){};

void CompositionOption::initialize() {
  add_help_suboption();
  m_desc.add_options()("display,d", "Display composition axes file")(
      "calc,c", "Calculate the standard composition axes")(
      "select,s", po::value<std::string>(&m_axis_choice_str),
      "Select composition axes");

  return;
}

const std::string &CompositionOption::axis_choice_str() const {
  return m_axis_choice_str;
}

}  // namespace Completer

// ///////////////////////////////////////
// 'composition' function for casm
//    (add an 'if-else' statement in casm.cpp to call this)

int composition_command(const CommandArgs &args) {
  po::variables_map vm;

  Completer::CompositionOption comp_opt;
  //("update,u", "Update composition and references based on current
  //'composition_axes.json' file");

  try {
    po::store(po::parse_command_line(args.argc(), args.argv(), comp_opt.desc()),
              vm);

    bool call_help = false;

    // quit out if there are no arguments
    if (!vm.count("help") && !vm.count("desc")) {
      if (vm.count("calc") + vm.count("select") + vm.count("display") +
              vm.count("update") !=
          1) {
        log() << "Error in 'casm composition'. You need to use either --calc, "
                 "--select, --display, or --update."
              << std::endl;
        call_help = true;
      }
    }

    /** --help option
     */
    if (vm.count("help") || call_help) {
      log() << std::endl;
      log() << comp_opt.desc() << std::endl;

      return 0;
    }

    if (vm.count("desc")) {
      log() << std::endl;
      log() << comp_opt.desc() << std::endl;

      log() << "DESCRIPTION" << std::endl;
      log() << "    Setup the composition axes.\n";
      log() << "    - expects a PRIM file in the project root directory \n";
      log() << "    - custom composition axes can be added to the "
               "'composition_axes.json' file \n";
      log() << "      and then compositions and references updated with 'casm "
               "composition --update'\n";
      log() << "      (word of caution: compatibility with PRIM is currently "
               "not checked for.)\n";
      log() << " \n";
      log() << "    Examples:\n";
      log() << "      casm composition --calc \n";
      log() << "      - Calculate standard composition axes, but does not "
               "select any. \n";
      log() << "\n";
      log() << "      casm composition --display \n";
      log() << "      - Display the possible composition axes\n";
      log() << "      - Possible axes are stored in the "
               "'composition_axes.json' file.\n";
      log() << "\n";
      log() << "      casm composition --select 0 \n";
      log() << "      - Select the composition axes by key name.\n";
      log() << "      - Updates and writes compositions and reference "
               "properties.\n";
      log() << "\n";
      /*log() << "      casm composition --update \n";
      log() << "      - Updates and writes compositions and reference properties
      based.\n"; log() << "        on the contents of the
      'composition_axes.json' file.\n"; log() << "\n";*/
      if (call_help) return ERR_INVALID_ARG;

      return 0;
    }

    po::notify(vm);

  } catch (po::error &e) {
    err_log() << "ERROR: " << e.what() << std::endl << std::endl;
    err_log() << comp_opt.desc() << std::endl;
    return ERR_INVALID_ARG;
  } catch (std::exception &e) {
    err_log() << "Unhandled Exception reached the top of main: " << e.what()
              << ", application will now exit" << std::endl;
    return ERR_UNKNOWN;
  }

  const fs::path &root = args.root;
  if (root.empty()) {
    err_log().error("No casm project found");
    err_log() << std::endl;
    return ERR_NO_PROJ;
  }

  // If 'args.primclex', use that, else construct PrimClex in 'uniq_primclex'
  // Then whichever exists, store reference in 'primclex'
  std::unique_ptr<PrimClex> uniq_primclex;
  PrimClex &primclex = make_primclex_if_not(args, uniq_primclex);

  const DirectoryStructure &dir = primclex.dir();
  fs::path comp_axes_path = dir.composition_axes();

  CompositionAxes comp_axes;
  if (fs::exists(comp_axes_path)) {
    comp_axes = read_composition_axes(comp_axes_path);
  }

  if (vm.count("display")) {
    log() << "\n***************************\n\n";

    display(log(), comp_axes);

    log() << "\n\n";

    if (comp_axes.err_code) {
      log() << comp_axes.err_message << std::endl;
    } else if (comp_axes.all_axes.size() && !comp_axes.has_current_axes()) {
      log() << "Please use 'casm composition --select' to choose your "
               "composition axes.\n\n";
    } else if (!comp_axes.all_axes.size()) {
      log() << "Please use 'casm composition --calc' to calculate standard "
               "composition axes.\n\n";
    }

    return 0;
  } else {
    if (vm.count("calc")) {
      log() << "\n***************************\n\n";

      log() << "Using the PRIM to enumerate standard composition axes for this "
               "space.\n\n";

      if (comp_axes.enumerated.size()) {
        log() << "Overwriting existing standard composition axes.\n\n";
      }

      comp_axes.erase_enumerated();

      std::vector<CompositionConverter> v;
      standard_composition_axes(xtal::allowed_molecule_names(primclex.prim()),
                                std::back_inserter(v));
      comp_axes.insert_enumerated(v.begin(), v.end());

      display(log(), comp_axes);

      log() << "\n\n";

      write_composition_axes(comp_axes_path, comp_axes);

      log() << "Wrote: " << comp_axes_path << "\n\n";

      if (!comp_axes.has_current_axes()) {
        log() << "Please use 'casm composition --select' to choose your "
                 "composition axes.\n\n";
      }
    }

    if (vm.count("select")) {
      log() << "\n***************************\n\n";

      if (comp_axes.all_axes.empty()) {
        log() << "Error: No composition axes found.\n\n";

        log() << "Please use 'casm composition --calc' to calculate standard "
                 "composition axes,\n"
              << "or add custom composition axes to " << comp_axes_path
              << "\n\n";

        return ERR_MISSING_DEPENDS;
      }

      if (comp_axes.all_axes.count(comp_opt.axis_choice_str()) == 0) {
        log() << "Error: The selected composition axes '"
              << comp_opt.axis_choice_str() << "' can not \n"
              << "be found in either the standard or custom compostion axes. "
                 "Please\n"
              << "re-select composition axes.                                  "
                 "   \n\n";

        return ERR_INVALID_INPUT_FILE;
      }

      comp_axes.select(comp_opt.axis_choice_str());

      display(log(), comp_axes);

      log() << "\n\n";

      write_composition_axes(comp_axes_path, comp_axes);

      log() << "Wrote: " << comp_axes_path << "\n\n";

      if (args.primclex) {
        args.primclex->refresh(false, true, true, false);
      };
    }

    return 0;
  }
}

}  // namespace CASM
