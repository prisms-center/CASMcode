#include <boost/filesystem/fstream.hpp>
#include "casm/crystallography/CoordinateSystems.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/crystallography/SimpleStructureTools.hh"
#include "casm/crystallography/SymTools.hh"
#include "casm/app/ProjectSettings.hh"
#include "casm/app/DirectoryStructure.hh"
#include "casm/app/AppIO.hh"
#include "casm/app/casm_functions.hh"
#include "casm/crystallography/io/VaspIO.hh"


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
      ("crystal-point-group", "Pretty print crystal point group")
      ("tol", po::value<double>(&m_tol)->default_value(1.0e-5), "tolerance to use with symmetrize in Angstroms default (1e-5)")
      ("symmetrize", po::value<fs::path>(&m_poscar_path)->value_name(ArgHandler::path()), "symmetrize a POSCAR specified by path to a given tolerance");

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
      po::store(po::parse_command_line(args.argc(), args.argv(), sym_opt.desc()), vm); // can throw

      /** --help option
      */
      if(vm.count("help")) {
        args.log() << "\n";
        args.log() << sym_opt.desc() << std::endl;

        return 0;
      }

      if(vm.count("desc")) {
        args.log() << "\n";
        args.log() << sym_opt.desc() << std::endl;
        args.log() << "DESCRIPTION" << std::endl;
        args.log() << "    Display symmetry group information.\n";

        return 0;
      }



      po::notify(vm); // throws on error, so do after help in case
      // there are any problems
    }
    catch(po::error &e) {
      args.err_log() << "ERROR: " << e.what() << std::endl << std::endl;
      args.err_log() << sym_opt.desc() << std::endl;
      return ERR_INVALID_ARG;
    }
    catch(std::exception &e) {
      args.err_log() << "Unhandled Exception reached the top of main: "
                     << e.what() << ", application will now exit" << std::endl;
      return ERR_UNKNOWN;

    }

    coordtype = sym_opt.coordtype_enum();
    COORD_MODE C(coordtype);

    const fs::path &root = args.root;
    if(root.empty()) {
      args.err_log().error("No casm project found");
      args.err_log() << std::endl;
      return ERR_NO_PROJ;
    }

    DirectoryStructure dir(root);
    ProjectSettings set(root);
    Structure prim(read_prim(dir.prim(), set.crystallography_tol(), &(set.hamiltonian_modules())));

    args.log() << "Generating lattice point group. " << std::endl << std::endl;
    SymGroup prim_pg(SymGroup::lattice_point_group(prim.lattice()));
    prim_pg.character_table();


    args.log() << "  Lattice point group size: " << prim_pg.size() << std::endl;
    args.log() << "  Lattice point group is: " << prim_pg.get_name() << std::endl << std::endl;

    args.log() << "Generating factor group. " << std::endl << std::endl;

    prim.generate_factor_group();
    prim.set_site_internals();

    args.log() << "  Factor group size: " << prim.factor_group().size() << std::endl;

    args.log() << "  Crystal point group is: " << prim.point_group().get_name() << std::endl;


    if(vm.count("lattice-point-group")) {
      args.log() << "\n***************************\n" << std::endl;
      args.log() << "Lattice point group:\n\n" << std::endl;
      prim_pg.print(args.log(), coordtype);
    }

    if(vm.count("factor-group")) {
      args.log() << "\n***************************\n" << std::endl;
      args.log() << "Factor group:\n\n" << std::endl;
      prim.factor_group().print(args.log(), coordtype);
    }

    if(vm.count("crystal-point-group")) {
      args.log() << "\n***************************\n" << std::endl;
      args.log() << "Crystal point group:\n\n" << std::endl;
      prim.point_group().print(args.log(), coordtype);
    }

    if(vm.count("symmetrize")) {
      fs::path poscar_path = sym_opt.m_poscar_path;
      double tol = sym_opt.m_tol;
      args.log() << "\n***************************\n" << std::endl;
      args.log() << "Symmetrizing: " << poscar_path << std::endl;
      args.log() << "with tolerance: " << tol << std::endl;
      Structure struc(poscar_path);
      Structure tprim;
      if(!struc.is_primitive(tprim)) {
        struc = tprim;
      }
      int biggest = struc.factor_group().size();
      Structure tmp = struc;
      // a) symmetrize the lattice vectors
      Lattice lat = tmp.lattice();
      lat = xtal::symmetrize(lat, tol);
      lat.set_tol(tol);

      tmp.set_lattice(lat, FRAC);

      tmp.factor_group();
      // b) find factor group with same tolerance
      tmp.fg_converge(tol);
      // c) symmetrize the basis sites
      SymGroup g = tmp.factor_group();
      tmp = xtal::symmetrize(tmp, g);

      //TODO: Why are we doing this twice?
      g = tmp.factor_group();
      tmp = xtal::symmetrize(tmp, g);
      if(tmp.factor_group().is_group(tol) && (tmp.factor_group().size() > biggest)) {
        struc = tmp;
      }
      if(!struc.is_primitive(tprim)) {
        struc = tprim;
      }
      fs::ofstream file_i;
      fs::path POSCARpath_i = "POSCAR_sym";
      file_i.open(POSCARpath_i);
      VaspIO::PrintPOSCAR p_i(make_simple_structure(struc), struc.title());
      p_i.print(file_i);
      file_i.close();
      return 0;

    }
    coordtype = sym_opt.coordtype_enum();


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

    args.log() << std::endl;

    return 0;

  };

}
