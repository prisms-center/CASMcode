#include <cstring>

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include "casm/app/AppIO.hh"
#include "casm/app/casm_functions.hh"
#include "casm/app/DirectoryStructure.hh"
#include "casm/app/HamiltonianModules.hh"
#include "casm/app/ProjectBuilder.hh"
#include "casm/app/ProjectSettings.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/clex/Configuration.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/SimpleStructureTools.hh"
#include "casm/completer/Handlers.hh"
#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/BasicStructureTools.hh"
#include "casm/crystallography/CanonicalForm.hh"
#include "casm/crystallography/DoFSet.hh"
#include "casm/crystallography/SimpleStructure.hh"
#include "casm/crystallography/SimpleStructureTools.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/global/definitions.hh"
#include "casm/misc/CASM_Eigen_math.hh"

namespace CASM {

  const std::string subproject_opt = "sub";

  const std::string write_prim_opt = "write-prim";

  const std::string relaxed_opt = "relaxed";

  const std::string molecule_opt = "as-molecules";

  const std::string include_va_opt = "include-va";

  // returns {error message, new file extension, new structure}
  std::tuple<std::string, std::string, BasicStructure > standardize_prim(BasicStructure const &prim, bool force) {
    /// Check if PRIM is primitive
    xtal::BasicStructure true_prim = xtal::make_primitive(prim);
    std::string err;
    if(true_prim.basis().size() < prim.basis().size()) {
      //err
      err +=
        "       The structure in the provided file is not primitive. Some CASM\n"
        "       features cannot be used with a non-primitive starting structure.\n\n";

      if(!force) {
        // current is_primitive
        Lattice lat_niggli = xtal::canonical::equivalent(true_prim.lattice());
        true_prim.set_lattice(lat_niggli, CART);
        return {err, ".primitive.json", true_prim};
      }
    }

    /// Check that the PRIM is in reduced form:
    Lattice niggli_lat = xtal::canonical::equivalent(prim.lattice());

    bool is_standard_niggli = almost_equal(niggli_lat.lat_column_mat(), prim.lattice().lat_column_mat());

    if(!is_standard_niggli) {
      err +=
        "       The structure in the provided file is not the Niggli cell in the CASM\n"
        "       standard orientation.\n\n";

      if(!force) {
        xtal::BasicStructure tmp(prim);
        tmp.set_lattice(niggli_lat, CART);
        return {err, ".niggli.json", tmp};
      }

    }

    //Check if the lattice is right handed, and if not print PRIM.right_handed.json
    if(!prim.lattice().is_right_handed()) {
      err +=
        "       The lattice vectors of the provided structure do not form a right-handed\n"
        "       coordinate system. Some electronic-structure codes will not accept this\n"
        "       input.\n\n";

      if(!force) {

        xtal::BasicStructure tmp(prim);
        tmp.set_lattice(Lattice(prim.lattice()).make_right_handed(), CART);
        tmp.within();
        return {err, ".right_handed.json", tmp};
      }
    }
    return {err, "", prim};
  }

  namespace Completer {
    InitOption::InitOption(): OptionHandlerBase("init") {}

    void InitOption::initialize() {
      add_help_suboption();
      add_prim_path_suboption("prim.json");
      add_file_path_suboption();
      add_configlist_suboption("NONE");
      add_confignames_suboption();
      add_coordtype_suboption();
      add_dofs_suboption();
      m_desc.add_options()
      (subproject_opt.c_str(), "Initialize a project in sub-directory of an existing project. After initialization, commands executed below the sub-directory will act on the sub-directory project; commmands executed above the sub-directory act on the original project.")
      (write_prim_opt.c_str(), "Create prim.json file for specified configuration(s). Each prim.json is written to the training_data directory of the corresponding configuration.")
      (relaxed_opt.c_str(), "Utilize relaxed coordinates for writing configuration-specific prim.json files.")
      (molecule_opt.c_str(), "Keep multi-atom species as molecules. By default, multi-atom species are split into constituent atoms.")
      (include_va_opt.c_str(), "Print sites that can only be vacancies; otherwise these sites are excluded.")
      ("force,f", "Force using a non-reduced, non-primitive, or left-handed PRIM.");
      return;
    }
  }

  // ///////////////////////////////////////
  // 'init' function for casm
  //    (add an 'if-else' statement in casm.cpp to call this)

  int init_command(const CommandArgs &args) {

    std::string name;
    po::variables_map vm;
    HamiltonianModules modules;

    /// Set command line options using boost program_options
    Completer::InitOption init_opt;

    try {
      po::store(po::parse_command_line(args.argc(), args.argv(), init_opt.desc()), vm); // can throw

      /** --help option
       */
      if(vm.count("help")) {
        log() << "\n";
        log() << init_opt.desc() << std::endl;

        return 0;
      }

      if(vm.count("desc")) {
        log() << "\n";
        log() << init_opt.desc() << std::endl;

        log() << "DESCRIPTION                                                \n" <<
              "    Initialize a new CASM project in the current directory.\n" <<
              "    - Expects a prim.json file in the current directory    \n" <<
              "    - If not found, looks for a PRIM file in the current   \n" <<
              "      directory and creates prim.json.                     \n\n";

        return 0;
      }

      po::notify(vm); // throws on error, so do after help in case
      // there are any problems
    }
    catch(po::error &e) {
      err_log() << "ERROR: " << e.what() << std::endl << std::endl;
      err_log() << init_opt.desc() << std::endl;
      return ERR_INVALID_ARG;
    }
    catch(std::exception &e) {
      err_log() << "Unhandled Exception reached the top of main: "
                << e.what() << ", application will now exit" << std::endl;
      return ERR_UNKNOWN;

    }


    fs::path root = fs::current_path();
    if(!init_opt.file_path().empty()) {
      root = init_opt.file_path();
    }
    else if(!args.root.empty()) {
      root = args.root;
    }

    // Going to do conversion:
    if(vm["config"].defaulted() && init_opt.config_strs().empty()) {
      DirectoryStructure dir(root);
      BasicStructure prim;
      std::string err_msg;
      std::string extension;
      // Read PRIM or prim.json:
      try {
        prim = read_prim(init_opt.prim_path(), TOL, &modules);
        std::tie(err_msg, extension, prim) = standardize_prim(prim, vm.count("force"));
      }
      catch(std::runtime_error &e) {
        err_log() << e.what() << std::endl;
        err_log() << "To initialize your project anyway, use the --force option." << std::endl;
        return ERR_INVALID_INPUT_FILE;
      }

      // Add the new DoFs before doing anything else
      std::map<DoFKey, xtal::DoFSet> new_global_dofs = prim.global_dofs();
      std::map<DoFKey, xtal::SiteDoFSet> new_site_dofs;
      for(std::string const &doftype : init_opt.dof_strs()) {
        AnisoValTraits ttraits(doftype);
        if(ttraits.global() && !new_global_dofs.count(doftype))
          new_global_dofs.emplace(doftype, ttraits);
        else
          new_site_dofs.emplace(doftype, ttraits);
      }
      prim.set_global_dofs(new_global_dofs);

      for(xtal::Site &site : prim.set_basis()) {
        std::map<DoFKey, xtal::SiteDoFSet> site_dofs = site.dofs();
        site_dofs.insert(new_site_dofs.begin(), new_site_dofs.end());
        site.set_dofs(std::move(site_dofs));
      }

      // Check error message to see if PRIM did not meet standards; report if so
      if(!err_msg.empty()) {
        if(vm.count("force")) {
          err_log()
              << "WARNING: " << init_opt.prim_path() << " failed the following check(s): \n"
              << err_msg
              << "Continuing due to usage of '--force' modifier.\n\n";
        }
        else {
          std::string new_path = init_opt.prim_path().string();
          std::string tpath = new_path.substr(0, new_path.find(".json"));
          new_path = tpath.substr(0, new_path.find(".JSON"));
          new_path += extension;
          err_log()
              << "Validation ERROR for " << init_opt.prim_path() << ":\n"
              << err_msg
              << "Standardizing structure and writing to JSON file: " << new_path << "\n";
          write_prim(prim, new_path, init_opt.coordtype_enum(), vm.count(include_va_opt));
          return ERR_INVALID_INPUT_FILE;
        }
      }
      if(vm.count(write_prim_opt)) {
        jsonParser json_prim;
        write_prim(prim, json_prim, init_opt.coordtype_enum(), vm.count(include_va_opt));
        log() << json_prim << std::endl;
      }
      else {
        log() << "\n***************************\n" << std::endl;
        // Actually going to initialize:
        root = fs::current_path();
        if(!init_opt.file_path().empty()) {
          root = init_opt.file_path();
        }
        fs::path existing = find_casmroot(root);
        if(vm.count(subproject_opt)) {
          if(existing == root) {
            log() << "Directory '" << root << "' is already the head directory of a casm project." << std::endl;
            return ERR_OTHER_PROJ;
          }
          if(existing.empty()) {
            log() << "Cannot create sub-directory project at '" << root
                  << "' because no existing project was found at a higher level." << std::endl
                  << "To initialize a top-level project, try again without the --" << subproject_opt << "flag." << std::endl;
            return ERR_OTHER_PROJ;

          }
        }
        else if(!existing.empty()) {
          log() << "Already in a casm project. To create a project at " << root
                << ", try again using the --" << subproject_opt << " flag." << std::endl;
          return ERR_OTHER_PROJ;
        }//End new control block

        try {

          //std::string orig_prim_title = prim.title;
          if(prim.title().empty() || !is_valid_project_name(prim.title())) {
            log() << "Please enter a short title for this project.\n";
            log() << "  Use something suitable as a prefix for files specific to this project, such as 'ZrO' or 'TiAl'.\n\n";

            std::string ttitle;
            log() << "Title: ";
            std::cin >> ttitle;
            log() << "\n\n";
            prim.set_title(ttitle);
          }
          log() << "Initializing CASM project '" << prim.title() << "'" << std::endl;
          auto project_settings = make_default_project_settings(prim, prim.title(), root);
          build_project(project_settings, Structure {prim});
        }
        catch(std::runtime_error &e) {
          err_log() << "ERROR: Could not build CASM project.\n";
          err_log() << e.what() << std::endl;
          return ERR_INVALID_INPUT_FILE;
        }

        log() << "  DONE" << std::endl;
        log() << std::endl;

      }
    }
    else if(vm.count(write_prim_opt)) {
      log() << "\n***************************\n" << std::endl;

      // Start from config in existing project:
      // If 'args.primclex', use that, else construct PrimClex in 'uniq_primclex'
      // Then whichever exists, store reference in 'primclex'

      std::unique_ptr<PrimClex> uniq_primclex;
      PrimClex &primclex = make_primclex_if_not(args, uniq_primclex);

      DB::Selection<Configuration> config_select(primclex.db<Configuration>(), init_opt.selection_path());
      for(std::string const &config : init_opt.config_strs()) {
        config_select.data()[config] = true;
      }
      for(const auto &config : config_select.selected()) {
        SimpleStructure::SpeciesMode species_mode = SimpleStructure::SpeciesMode::ATOM;
        if(vm.count(molecule_opt))
          species_mode = SimpleStructure::SpeciesMode::MOL;
        SimpleStructure sstruc = make_simple_structure(config, {}, vm.count(relaxed_opt));
        BasicStructure new_prim = xtal::make_basic_structure(sstruc, init_opt.dof_strs(), species_mode);
        fs::path config_dir = primclex.dir().configuration_dir(config.name());
        try {
          fs::create_directories(config_dir);
        }
        catch(const fs::filesystem_error &ex) {
          err_log() << "Error making directory " << config_dir << std::endl;
          err_log() << ex.what() << std::endl;

        }

        std::string prim_title = config.name();
        auto it = prim_title.begin();
        for(; it != prim_title.end(); ++it) {
          if((*it) == '/') {
            (*it) = 'c';
            break;
          }
        }

        prim_title.insert(it, '_');

        new_prim.set_title(prim_title);

        write_prim(new_prim, config_dir / "prim.json", init_opt.coordtype_enum(), vm.count(include_va_opt));

        log() << "Wrote file " << (config_dir / "prim.json") << std::endl;
      }
    }


    return 0;
  };

}
