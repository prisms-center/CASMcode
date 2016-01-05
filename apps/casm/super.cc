#include "super.hh"

#include <cstring>

#include "casm_functions.hh"
#include "casm/CASM_classes.hh"

namespace CASM {


  // ///////////////////////////////////////
  // 'super' function for casm
  //    (add an 'if-else' statement in casm.cpp to call this)

  int super_command(int argc, char *argv[]) {

    //casm enum [—supercell min max] [—config supercell ] [—hopconfigs hop.background]
    //- enumerate supercells and configs and hop local configurations

    std::string scelname, configname, unitscelname;
    fs::path tmatfile, structfile;
    fs::path abs_tmatfile, abs_structfile;
    double tol;
    COORD_TYPE coordtype;
    po::variables_map vm;

    try {

      /// Set command line options using boost program_options
      po::options_description desc("'casm super' usage");
      desc.add_options()
      ("help,h", "Write help documentation")
      ("transf_mat", po::value<fs::path>(&tmatfile), "File containing a 3x3 transformation matrix used to create a supercell.")
      ("get_transf_mat", "If it exists, find the transformation matrix.")
      ("structure", po::value<fs::path>(&structfile), "File with structure (POSCAR type) to use.")
      ("configname", po::value<std::string>(&configname), "Name of a configuration. For ex. \"SCEL4_2_2_1_0_0_0/4\".")
      ("scelname", po::value<std::string>(&scelname), "Name of supercell. For ex. \"SCEL4_2_2_1_0_0_0\".")
      ("unitcell", po::value<std::string>(&unitscelname), "Name of supercell to use as unit cell. For ex. \"SCEL2_2_1_1_0_0_0\".")
      ("add-canonical,a", "Will add the generated super configuration in it's canonical form in the equivalent niggli supercell.")
      ("vasp5", "Print using VASP5 style (include atom name line)")
      ("tol", po::value<double>(&tol)->default_value(CASM::TOL), "Tolerance used for checking symmetry")
      ("coord", po::value<COORD_TYPE>(&coordtype)->default_value(CASM::CART), "Coord mode: FRAC=0, or CART=1");

      try {
        po::store(po::parse_command_line(argc, argv, desc), vm); // can throw

        if(!vm.count("help")) {
          if(vm.count("transf_mat") + vm.count("get_transf_mat") != 1) {
            std::cerr << "Error in 'casm super'. Either the --transf_mat or --get_transf_mat option must be included." << std::endl;
            return 1;
          }
        }

        /** --help option
        */
        if(vm.count("help")) {
          std::cout << "\n";
          std::cout << desc << std::endl;

          std::cout << "DESCRIPTION" << std::endl;
          std::cout << "                                                                      \n" <<
                    "  casm super --transf_mat T                                           \n" <<
                    "  - Print super lattice of the PRIM lattice                           \n" <<
                    "                                                                      \n" <<
                    "  casm super --structure POSCAR --transf_mat T                        \n" <<
                    "  - Print superstructure of a POSCAR                                  \n" <<
                    "                                                                      \n" <<
                    "  casm super --configname configname --transf_mat T                   \n" <<
                    "  - Print superstructure of a configuration                           \n" <<
                    "                                                                      \n" <<
                    "  casm super --structure POSCAR --unitcell scelname --get_transf_mat  \n" <<
                    "  - Check if POSCAR lattice is a supercell of unit cell lattice and   \n" <<
                    "    if so print the transformation matrix                             \n" <<
                    "  - Uses primitive cell for unitcell if none given                    \n" <<
                    "                                                                      \n" <<
                    "  casm super --scelname scelname --unitcell scelname --get_transf_mat\n" <<
                    "  - Check if configuration lattice is a supercell of unit cell lattice.\n" <<
                    "    and print the transformation matrix                               \n" <<
                    "  - Uses primitive cell for unitcell if none given                    \n\n";

          return 0;
        }

        po::notify(vm); // throws on error, so do after help in case
        // there are any problems

      }
      catch(po::error &e) {
        std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
        std::cerr << desc << std::endl;
        return 1;
      }
    }
    catch(std::exception &e) {
      std::cerr << "Unhandled Exception reached the top of main: "
                << e.what() << ", application will now exit" << std::endl;
      return 1;

    }


    // -- no casm project necessary for super cell of a POSCAR -------

    // want absolute paths
    abs_tmatfile = fs::absolute(tmatfile);
    abs_structfile = fs::absolute(structfile);

    if(vm.count("structure") && vm.count("transf_mat")) {

      if(!fs::exists(abs_structfile)) {
        std::cout << "ERROR: " << abs_tmatfile << " not found." << std::endl;
        return 1;
      }
      BasicStructure<Site> unitcell(abs_structfile);

      // -- read transf matrix ---
      if(!fs::exists(abs_tmatfile)) {
        std::cout << "ERROR: " << abs_tmatfile << " not found." << std::endl;
        return 1;
      }
      Matrix3<int> Tm;
      fs::ifstream file(abs_tmatfile);
      file >> Tm;
      file.close();

      auto super = unitcell.create_superstruc(make_supercell(unitcell.lattice(), Tm));

      if(vm.count("vasp5")) {
        super.print5(std::cout);
      }
      else {
        super.print(std::cout);
      }

      return 0;
    }


    COORD_MODE C(coordtype);

    fs::path orig = fs::current_path();
    fs::path root = find_casmroot(fs::current_path());
    if(root.empty()) {
      std::cout << "Error: No casm project found." << std::endl;
      return 1;
    }
    fs::current_path(root);


    std::cout << "\n***************************\n" << std::endl;

    // initialize primclex
    std::cout << "Initialize primclex: " << root << std::endl << std::endl;
    PrimClex primclex(root, std::cout);
    std::cout << "  DONE." << std::endl << std::endl;

    if(vm.count("transf_mat")) {

      Matrix3<int> Tm;
      if(!fs::exists(abs_tmatfile)) {
        std::cout << "ERROR: " << abs_tmatfile << " not found." << std::endl;
        return 1;
      }
      fs::ifstream file(abs_tmatfile);
      file >> Tm;
      file.close();

      Eigen::Matrix3i T = Tm;

      std::cout << "Read transformation matrix, T: \n" << T << "\n\n";


      // super structure
      if(vm.count("configname")) {

        std::stringstream ss;
        const Configuration &con = primclex.configuration(configname);
        const Supercell &scel = con.get_supercell();
        con.print(ss, FRAC);

        std::istringstream iss(ss.str());
        BasicStructure<Site> unit;
        unit.read(iss);

        std::cout << "Unit structure:";
        std::cout << "\n------\n";
        if(vm.count("vasp5")) {
          unit.print5(std::cout);
        }
        else {
          unit.print(std::cout);
        }
        std::cout << "\n------\n";
        std::cout << "\n\n";


        BasicStructure<Site> super = unit.create_superstruc(make_supercell(unit.lattice(), T));
        super.title = std::string("Supercell of ") + con.name();

        std::cout << "Super structure:";
        std::cout << "\n------\n";
        if(vm.count("vasp5")) {
          super.print5(std::cout);
        }
        else {
          super.print(std::cout);
        }
        std::cout << "\n------\n";
        
        if(vm.count("add-canonical")) {
          
          int map_opt = ConfigMapper::none;
          double tol = TOL;
          double vol_tol = 0.25;
          double lattice_weight = 0.5;
          ConfigMapper configmapper(primclex, lattice_weight, vol_tol, map_opt, tol);

          std::string imported_name;
          Eigen::Matrix3d cart_op;
          std::vector<Index> best_assignment;
          jsonParser fullrelax_data;
          if(configmapper.import_structure_occupation(super, 
                                                      imported_name, 
                                                      fullrelax_data, 
                                                      best_assignment, 
                                                      cart_op, 
                                                      true)) {
            std::cout << "  The configuration was imported successfully as " 
                      << imported_name << std::endl << std::endl;
            
          }
          else {
            std::cout << "  The configuration was mapped onto pre-existing equivalent structure " 
                      << imported_name << std::endl << std::endl;
          }
          
          jsonParser json_src;
          json_src["supercell_of"] = configname;
          primclex.configuration(imported_name).push_back_source(json_src);
          
          //Update directories
          std::cout << "  Writing SCEL..." << std::endl;
          primclex.print_supercells();
          std::cout << "  Writing config_list..." << std::endl << std::endl;
          primclex.write_config_list();
          std::cout << "  DONE" << std::endl << std::endl;

        }

        return 0;

      }
      // super lattice of prim lattice
      else {

        BasicStructure<Site> unit = primclex.get_prim();
        SymGroup pg = Structure(unit).point_group();

        Eigen::Matrix3d U = unit.lattice().lat_column_mat();
        Eigen::Matrix3d S = U * T.cast<double>();
        Eigen::Matrix3i H_canon;
        Eigen::Matrix3d op_canon;

        Eigen::Matrix3d S_niggli = niggli(Lattice(S), pg, tol).lat_column_mat();
        Eigen::Matrix3i T_niggli = iround(U.inverse() * S_niggli);
        Eigen::Matrix3i H_niggli = hermite_normal_form(T_niggli).first;

        std::stringstream s_name;
        s_name << "SCEL" << H_niggli(0, 0)*H_niggli(1, 1)*H_niggli(2, 2) << "_"
               << H_niggli(0, 0) << "_" << H_niggli(1, 1) << "_" << H_niggli(2, 2) << "_"
               << H_niggli(1, 2) << "_" << H_niggli(0, 2) << "_" << H_niggli(0, 1);

        // S = U*T;
        // S_canon = op_canon*S = U*T', where H_canon*V = U.inv*op_canon*U*T = T'

        std::cout << "--- Lattices as column vector matrices ---\n\n";

        std::cout << "Prim lattice, U:\n" << U << "\n\n";

        std::cout << "Super lattice, S = U*T:\n" << S << "\n\n";

        std::cout << "This is equivalent to '" << s_name.str() << "', the equivalent super lattice \n" <<
                  "in the standard orientation niggli cell, S_niggli:\n" << S_niggli << "\n\n";

        std::cout << "The transformation matrix (S_niggli = U*T) for '" << s_name.str() << "' is:\n" << T_niggli << "\n\n";


        std::cout << "--- Lattices as row vector matrices ---\n\n";

        std::cout << "Prim lattice:\n" << U.transpose() << "\n\n";

        std::cout << "Super lattice:\n" << S.transpose() << "\n\n";

        std::cout << "This is equivalent to '" << s_name.str() << "', the equivalent super lattice \n" <<
                  "in the standard orientation niggli cell:\n" << S_niggli.transpose() << "\n\n";

        return 0;
      }

    }

    if(vm.count("get_transf_mat")) {

      Matrix3<double> tmat;

      Lattice unit_lat = primclex.get_prim().lattice();

      if(vm.count("unitcell")) {
        unit_lat = primclex.get_supercell(unitscelname).get_real_super_lattice();
      }

      Lattice super_lat;

      if(vm.count("structure")) {
        super_lat = BasicStructure<Site>(abs_structfile).lattice();
      }
      else if(vm.count("scelname")) {
        super_lat = primclex.get_supercell(scelname).get_real_super_lattice();
      }
      else {
        std::cout << "Error in 'casm super --get_transf_mat'. No --structure or --scelname given." << std::endl << std::endl;
        return 1;
      }

      std::cout << "--- Lattices as column vector matrices ---\n\n";

      std::cout << "Unit lattice, U:\n" << unit_lat.lat_column_mat() << "\n\n";

      std::cout << "Super lattice, S:\n" << super_lat.lat_column_mat() << "\n\n";

      // see if super_lat is a supercell of unitlat
      // S == U*T
      Matrix3<double> T = unit_lat.lat_column_mat().inverse() * super_lat.lat_column_mat();

      if(T.is_integer() && !T.is_zero()) {
        std::cout << "The super lattice is a supercell of the unit lattice.\n\n";

        std::cout << "The transformation matrix, T, where S = U*T, is: \n" << iround(T) << "\n\n";
      }
      else {
        std::cout << "The super lattice is NOT a supercell of the unit lattice.\n\n";

        std::cout << "The transformation matrix, T, where S = U*T, is: \n" << T << "\n\n";
      }

      return 0;
    }

    return 0;
  };

}


