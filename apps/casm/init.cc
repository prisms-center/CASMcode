#include "init.hh"

#include <cstring>

#include <casm/core>
#include "casm/app/DirectoryStructure.hh"
#include "casm/app/ProjectBuilder.hh"

namespace CASM {

  // ///////////////////////////////////////
  // 'init' function for casm
  //    (add an 'if-else' statement in casm.cpp to call this)

  int init_command(int argc, char *argv[]) {

    std::string name;
    po::variables_map vm;

    try {

      /// Set command line options using boost program_options
      po::options_description desc("'casm init' usage");
      desc.add_options()
      ("help,h", "Write help documentation")
      ("force,f", "Force using a non-reduced, non-primitive, or left-handed PRIM");

      try {
        po::store(po::parse_command_line(argc, argv, desc), vm); // can throw

        /** --help option
        */
        if(vm.count("help")) {
          std::cout << "\n";
          std::cout << desc << std::endl;

          std::cout << "DESCRIPTION                                                \n" <<
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
        std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
        std::cerr << desc << std::endl;
        return ERR_INVALID_ARG;
      }
    }
    catch(std::exception &e) {
      std::cerr << "Unhandled Exception reached the top of main: "
                << e.what() << ", application will now exit" << std::endl;
      return ERR_UNKNOWN;

    }

    if(fs::current_path() == find_casmroot(fs::current_path())) {
      std::cout << "Already in a casm project." << std::endl;
      return ERR_OTHER_PROJ;
    }

    std::cout << "\n***************************\n" << std::endl;

    fs::path root = fs::current_path();

    DirectoryStructure dir(root);
    Structure prim;


    // if prim.json does not exist, try to read PRIM and create prim.json
    if(!fs::is_regular_file(dir.prim())) {

      if(!fs::is_regular_file(dir.PRIM())) {
        std::cout << "Error in 'casm init': Neither 'prim.json' nor 'PRIM' found.\n\n";

        std::cout << "Run 'casm format --prim' for the format of the 'prim.json' file.\n\n";

        std::cout << "For step by step help use: 'casm status -n'\n\n";

        return ERR_MISSING_INPUT_FILE;
      }

      try {
        fs::ifstream poscar_prim;
        poscar_prim.open(dir.PRIM());
        prim.read(poscar_prim);
        poscar_prim.close();
      }
      catch(std::runtime_error &e) {

        std::cerr << "ERROR: No prim.json exists. PRIM exists, but it could not be read.\n";
        std::cerr << e.what() << std::endl;
        return ERR_INVALID_INPUT_FILE;
      }


      std::string poscar_prim_title = prim.title;
      std::cout << "Converting 'PRIM' to 'prim.json'.\n\n" << std::endl;

      std::cout << "Please enter a short title for this project.\n";
      std::cout << "  Use something suitable as a prefix for files specific to this project, such as 'ZrO' or 'TiAl'.\n\n";

      std::cout << "Title: ";
      std::cin >> prim.title;
      std::cout << "\n\n";

      jsonParser json;
      write_prim(prim, json, FRAC);
      json["description"] = poscar_prim_title;
      fs::ofstream primfile(dir.prim());
      json.print(primfile);
      primfile.close();

    }

    jsonParser prim_json;

    try {

      prim_json = jsonParser(dir.prim());

      prim = Structure(read_prim(prim_json));
    }
    catch(std::runtime_error &e) {
      std::cerr << e.what() << std::endl;

      return ERR_INVALID_INPUT_FILE;
    }

    /// Check if PRIM is primitive
    BasicStructure<Site> true_prim;
    true_prim.title = prim.title;
    if(!prim.is_primitive(true_prim)) {
      if(!vm.count("force")) {
        std::cerr << "ERROR: The structure in the prim.json file is not primitive. Writing the most       \n"
                  << "       primitive structure to file 'prim.true.json'.\n\n";

        Structure tmp(true_prim);
        Lattice lat_niggli = niggli(true_prim.lattice(), tmp.point_group(), TOL);
        tmp.set_lattice(lat_niggli, CART);

        fs::ofstream primfile(root / "prim.true.json");
        jsonParser json;
        write_prim(tmp, json, FRAC);
        json["description"] = prim_json["description"];
        json.print(primfile);
        primfile.close();

        std::cerr << "If you want to use the current prim.json anyway, re-run with the --force option. Some\n"
                  << "CASM features cannot be used with a non-primitive starting structure.\n";
        return ERR_INVALID_INPUT_FILE;
      }
      else {
        std::cerr << "WARNING: The structure in the prim.json file is not primitive. Continuing anyway    \n"
                  << "         because the --force option is on.\n\n";
      }
    }

    /// Check that the PRIM is in reduced form:
    Lattice niggli_lat = niggli(prim.lattice(), prim.point_group(), TOL);

    bool is_standard_niggli = niggli_lat.lat_column_mat().is_equal(prim.lattice().lat_column_mat());

    if(!is_standard_niggli) {
      if(!vm.count("force")) {
        if(!is_standard_niggli) {
          std::cerr << "ERROR: The structure in the prim.json file is not the niggli cell in the CASM standard\n" <<
                    "       orientation. Writing the suggested structure to 'prim.niggli.json'.\n\n";
        }

        Structure tmp(true_prim);
        Lattice lat_niggli = niggli(true_prim.lattice(), tmp.point_group(), TOL);
        tmp.set_lattice(lat_niggli, CART);

        fs::ofstream primfile(root / "prim.niggli.json");
        jsonParser json;
        write_prim(tmp, json, FRAC);
        json["description"] = prim_json["description"];
        json.print(primfile);
        primfile.close();

        if(!is_standard_niggli) {
          std::cerr << "If you want to use the current prim.json anyway, re-run with the --force option. Some\n" <<
                    "CASM features cannot be used with a non-reduced starting structure.\n";
          return ERR_INVALID_INPUT_FILE;
        }
      }
      else {
        std::cerr << "WARNING: The structure in the prim.json file is not the standard orientation Niggli\n"
                  << "         cell. Continuing anyway because the --force option is on.\n\n";
        //return 1;
      }

    }

    //Check if the lattice is right handed, and if not print PRIM.right_handed.json
    if(!prim.lattice().is_right_handed()) {
      if(!vm.count("force")) {
        std::cerr << "ERROR: The structure in prim.json is not right-handed. Some electronic-"
                  << "structure codes will not accept this input. If you would like to "
                  << "keep this PRIM, re-run with the --force option. Writing the "
                  << "right-handed structure to PRIM.right_handed.json" << std::endl;

        prim.set_lattice(Lattice(prim.lattice()).make_right_handed(), CART);
        prim.within();

        fs::ofstream primfile(root / "prim.right_handed.json");
        jsonParser json;
        write_prim(prim, json, FRAC);
        json["description"] = prim_json["description"];
        json.print(primfile);
        primfile.close();
        return ERR_INVALID_INPUT_FILE;
      }
      else {
        std::cerr << "WARNING: The structure in the prim.json file is not right-handed. Continuing anyway    \n"
                  << "         because the --force option is on.\n\n";
      }
    }

    std::cout << "Initializing CASM project '" << prim.title << "'" << std::endl;

    ProjectBuilder builder(root, prim.title, "formation_energy");
    builder.build();

    std::cout << "  DONE" << std::endl;
    std::cout << std::endl;

    return 0;
  };

}

