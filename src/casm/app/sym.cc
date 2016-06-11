#include "casm/CASM_classes.hh"
#include "casm/app/casm_functions.hh"

#include "casm/completer/Handlers.hh"

namespace CASM {

  namespace Completer {
    SymOption::SymOption(): OptionHandlerBase("sym") {}

    void SymOption::initialize() {
      add_help_suboption();
      add_coordtype_suboption();
      m_desc.add_options()
      ("lattice-point-group", "Pretty print lattice point group")
      ("factor-group", "Pretty print factor group")
      ("crystal-point-group", "Pretty print crystal point group");

      return;
    }
  }


  // ///////////////////////////////////////
  // 'sym' function for casm
  //    (add an 'if-else' statement in casm.cpp to call this)

  int sym_command(const CommandArgs &args) {
    //std::string name;
    COORD_TYPE coordtype;
    po::variables_map vm;

    /// Set command line options using boost program_options
    Completer::SymOption sym_opt;
    try {
      po::store(po::parse_command_line(args.argc, args.argv, sym_opt.desc()), vm); // can throw
      coordtype = sym_opt.coordtype_enum();

      /** --help option
      */
      if(vm.count("help")) {
        std::cout << "\n";
        std::cout << sym_opt.desc() << std::endl;

        std::cout << "DESCRIPTION" << std::endl;
        std::cout << "    Display symmetry group information.\n";

        return 0;
      }

      po::notify(vm); // throws on error, so do after help in case
      // there are any problems
    }
    catch(po::error &e) {
      std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
      std::cerr << sym_opt.desc() << std::endl;
      return ERR_INVALID_ARG;
    }
    catch(std::exception &e) {
      std::cerr << "Unhandled Exception reached the top of main: "
                << e.what() << ", application will now exit" << std::endl;
      return ERR_UNKNOWN;

    }

    COORD_MODE C(coordtype);

    const fs::path &root = args.root;
    if(root.empty()) {
      args.err_log.error("No casm project found");
      args.err_log << std::endl;
      return ERR_NO_PROJ;
    }

    DirectoryStructure dir(root);
    ProjectSettings set(root);
    Structure prim(read_prim(dir.prim()));

    std::cout << "Generating lattice point group. " << std::endl << std::endl;
    SymGroup prim_pg;
    prim.lattice().generate_point_group(prim_pg, set.crystallography_tol());
    prim_pg.character_table();


    std::cout << "  Lattice point group size: " << prim_pg.size() << std::endl;
    std::cout << "  Lattice point group is: " << prim_pg.get_name() << std::endl << std::endl;

    std::cout << "Generating factor group. " << std::endl << std::endl;

    prim.generate_factor_group(set.crystallography_tol());
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


