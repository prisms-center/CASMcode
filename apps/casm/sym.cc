#include "sym.hh"

#include<cstring>

#include "casm/CASM_classes.hh"
#include "casm_functions.hh"

namespace CASM {

  // ///////////////////////////////////////
  // 'sym' function for casm
  //    (add an 'if-else' statement in casm.cpp to call this)

  int sym_command(const CommandArgs& args) {
    std::string name;
    COORD_TYPE coordtype;
    po::variables_map vm;

    try {

      /// Set command line options using boost program_options
      po::options_description desc("'casm sym' usage");
      desc.add_options()
      ("help,h", "Write help documentation")
      ("lattice-point-group", "Pretty print lattice point group")
      ("factor-group", "Pretty print factor group")
      ("crystal-point-group", "Pretty print crystal point group")
      ("coord", po::value<COORD_TYPE>(&coordtype)->default_value(CASM::CART), "Coord mode: FRAC=0, or CART=1");

      try {
        po::store(po::parse_command_line(args.argc, args.argv, desc), vm); // can throw

        /** --help option
        */
        if(vm.count("help")) {
          std::cout << "\n";
          std::cout << desc << std::endl;

          std::cout << "DESCRIPTION" << std::endl;
          std::cout << "    Display symmetry group information.\n";

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

    COORD_MODE C(coordtype);

    fs::path &root = args.root;
    
    DirectoryStructure dir(root);
    ProjectSettings set(root);
    Structure prim(read_prim(dir.prim()));

    std::cout << "Generating lattice point group. " << std::endl << std::endl;
    SymGroup prim_pg;
    prim.lattice().generate_point_group(prim_pg, set.tol());
    prim_pg.character_table();


    std::cout << "  Lattice point group size: " << prim_pg.size() << std::endl;
    std::cout << "  Lattice point group is: " << prim_pg.get_name() << std::endl << std::endl;

    std::cout << "Generating factor group. " << std::endl << std::endl;

    prim.generate_factor_group(set.tol());
    prim.set_site_internals();

    std::cout << "  Factor group size: " << prim.factor_group().size() << std::endl;

    std::cout << "  Crystal point group is: " << prim.point_group().get_name() << std::endl;


    if(vm.count("lattice-point-group")) {
      std::cout << "\n***************************\n" << std::endl;
      std::cout << "Lattice point group:\n\n" << std::endl;
      prim_pg.print(std::cout, coordtype);
    }

    if(vm.count("factor-group")) {
      std::cout << "\n***************************\n" << std::endl;
      std::cout << "Factor group:\n\n" << std::endl;
      prim.factor_group().print(std::cout, coordtype);
    }

    if(vm.count("crystal-point-group")) {
      std::cout << "\n***************************\n" << std::endl;
      std::cout << "Crystal point group:\n\n" << std::endl;
      prim.point_group().print(std::cout, coordtype);
    }

    // Write symmetry info files
    set.new_symmetry_dir();

    // Write lattice point group
    {
      fs::ofstream outfile;
      jsonParser json;
      outfile.open(dir.lattice_point_group());
      write_symgroup(prim_pg, json);
      json.print(outfile);
      outfile.close();
    }

    // Write factor group
    {
      fs::ofstream outfile;
      jsonParser json;
      outfile.open(dir.factor_group());
      write_symgroup(prim.factor_group(), json);
      json.print(outfile);
      outfile.close();
    }

    // Write crystal point group
    {
      fs::ofstream outfile;
      jsonParser json;
      outfile.open(dir.crystal_point_group());
      write_symgroup(prim.point_group(), json);
      json.print(outfile);
      outfile.close();
    }

    std::cout << std::endl;

    return 0;

  };

}


