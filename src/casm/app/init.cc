#include <cstring>

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include "casm/global/definitions.hh"
#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/app/DirectoryStructure.hh"
#include "casm/app/HamiltonianModules.hh"
#include "casm/app/ProjectBuilder.hh"
#include "casm/app/casm_functions.hh"
#include "casm/app/AppIO.hh"
#include "casm/completer/Handlers.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/Configuration.hh"
#include "casm/crystallography/CanonicalForm.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/crystallography/BasicStructure_impl.hh"
#include "casm/crystallography/SimpleStructure.hh"
#include "casm/crystallography/SimpleStructureTools.hh"

namespace CASM {

  const std::string subproject_opt = "sub";

  const std::string write_prim_opt = "write-prim";

  const std::string relaxed_opt = "relaxed";

  const std::string molecule_opt = "as-molecules";

  // returns {error message, new file extension, new structure}
  std::tuple<std::string, std::string, BasicStructure<Site> > standardize_prim(BasicStructure<Site> const &prim, bool force) {
    /// Check if PRIM is primitive
    Structure true_prim;
    std::string err;
    true_prim.set_title(prim.title());
    if(!prim.is_primitive(true_prim)) {
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
        Structure tmp(prim);
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

        Structure tmp(prim);
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
        args.log() << "\n";
        args.log() << init_opt.desc() << std::endl;

        return 0;
      }

      if(vm.count("desc")) {
        args.log() << "\n";
        args.log() << init_opt.desc() << std::endl;

        args.log() << "DESCRIPTION                                                \n" <<
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
      args.err_log() << "ERROR: " << e.what() << std::endl << std::endl;
      args.err_log() << init_opt.desc() << std::endl;
      return ERR_INVALID_ARG;
    }
    catch(std::exception &e) {
      args.err_log() << "Unhandled Exception reached the top of main: "
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

    args.log() << "\n***************************\n" << std::endl;
    // Going to do conversion:
    if(vm["config"].defaulted() || !vm.count(write_prim_opt)) {
      DirectoryStructure dir(root);
      BasicStructure<Site> prim;
      std::string err_msg;
      std::string extension;
      // Read PRIM or prim.json:
      try {
        prim = Structure(read_prim(init_opt.prim_path(), TOL, &modules));

        std::tie(err_msg, extension, prim) = standardize_prim(prim, vm.count("force"));
      }
      catch(std::runtime_error &e) {
        args.err_log() << e.what() << std::endl;

        return ERR_INVALID_INPUT_FILE;
      }

      // Add the new DoFs before doing anything else
      for(std::string const &doftype : init_opt.dof_strs()) {
        prim.add_dof(DoFSet::make_default(doftype));
      }

      // Check error message to see if PRIM did not meet standards; report if so
      if(!err_msg.empty()) {
        if(vm.count("force")) {
          args.err_log()
              << "WARNING: " << init_opt.prim_path() << " failed the following check(s): \n"
              << err_msg
              << "Continuing due to usage of '--force' modifier.\n\n";
        }
        else {
          std::string new_path = init_opt.prim_path().string();
          std::string tpath = new_path.substr(0, new_path.find(".json"));
          new_path = tpath.substr(0, new_path.find(".JSON"));
          new_path += extension;
          args.err_log()
              << "Validation ERROR for " << init_opt.prim_path() << ":\n"
              << err_msg
              << "Standardizing structure and writing to JSON file: " << new_path << "\n";
          write_prim(prim, new_path, init_opt.coordtype_enum());
          return ERR_INVALID_INPUT_FILE;
        }
      }
      if(vm.count(write_prim_opt)) {
        jsonParser json_prim;
        write_prim(prim, json_prim, init_opt.coordtype_enum());
        args.log() << json_prim << std::endl;
      }
      else {
        // Actually going to initialize:
        fs::path existing = find_casmroot(root);
        if(vm.count(subproject_opt)) {
          if(existing == root) {
            args.log() << "Directory '" << root << "' is already the head directory of a casm project." << std::endl;
            return ERR_OTHER_PROJ;
          }
          if(!root.empty()) {
            args.log() << "No existing project found. Cannot create sub-directory project at '" << root << "'." << std::endl;
            return ERR_OTHER_PROJ;

          }
        }
        else if(!existing.empty()) {
          args.log() << "Already in a casm project." << std::endl;
          return ERR_OTHER_PROJ;
        }//End new control block

        try {

          //std::string orig_prim_title = prim.title;
          if(prim.title().empty() || !ProjectBuilder::valid_title(prim.title())) {
            args.log() << "Please enter a short title for this project.\n";
            args.log() << "  Use something suitable as a prefix for files specific to this project, such as 'ZrO' or 'TiAl'.\n\n";

            std::string ttitle;
            args.log() << "Title: ";
            std::cin >> ttitle;
            args.log() << "\n\n";
            prim.set_title(ttitle);
          }
          args.log() << "Initializing CASM project '" << prim.title() << "'" << std::endl;
          ProjectBuilder builder(prim, root, prim.title(), "formation_energy");
          builder.build();
        }
        catch(std::runtime_error &e) {
          args.err_log() << "ERROR: Could not build CASM project.\n";
          args.err_log() << e.what() << std::endl;
          return ERR_INVALID_INPUT_FILE;
        }

        args.log() << "  DONE" << std::endl;
        args.log() << std::endl;

      }
    }
    else if(vm.count(write_prim_opt)) {
      // Start from config in existing project:
      // If 'args.primclex', use that, else construct PrimClex in 'uniq_primclex'
      // Then whichever exists, store reference in 'primclex'

      std::unique_ptr<PrimClex> uniq_primclex;
      PrimClex &primclex = make_primclex_if_not(args, uniq_primclex);

      DB::Selection<Configuration> config_select(primclex.db<Configuration>(), init_opt.selection_path());
      for(const auto &config : config_select.selected()) {
        SimpleStructure::SpeciesMode species_mode = SimpleStructure::SpeciesMode::ATOM;
        if(vm.count(molecule_opt))
          species_mode = SimpleStructure::SpeciesMode::MOL;
        SimpleStructure sstruc = xtal::make_simple_structure(config, {}, vm.count(relaxed_opt));
        BasicStructure<Site> new_prim = make_basic_structure(sstruc, init_opt.dof_strs(), species_mode);
        fs::path config_dir = primclex.dir().configuration_dir(config.name());
        try {
          fs::create_directories(config_dir);
        }
        catch(const fs::filesystem_error &ex) {
          args.err_log() << "Error making directory " << config_dir << std::endl;
          args.err_log() << ex.what() << std::endl;

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

        write_prim(new_prim, config_dir / "prim.json", init_opt.coordtype_enum());

        args.log() << "Wrote file " << (config_dir / "prim.json") << std::endl;
      }
    }


    return 0;
  };

}
