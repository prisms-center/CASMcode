#include <boost/filesystem/fstream.hpp>
#include <ostream>
#include "casm/crystallography/BasicStructureTools.hh"
#include "casm/crystallography/CoordinateSystems.hh"
#include "casm/crystallography/Lattice.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/crystallography/SimpleStructureTools.hh"
#include "casm/crystallography/SymTools.hh"
#include "casm/app/ProjectSettings.hh"
#include "casm/app/DirectoryStructure.hh"
#include "casm/app/AppIO.hh"
#include "casm/app/casm_functions.hh"
#include "casm/crystallography/io/VaspIO.hh"


#include "casm/completer/Handlers.hh"
#include "casm/symmetry/SymGroup.hh"

namespace {
  using namespace CASM;
  void print_factor_group_convergence(const Structure &struc, double small_tol, double large_tol, double increment, std::ostream &print_stream) {
    std::vector<double> tols;
    std::vector<bool> is_group;
    std::vector<int> num_ops, num_enforced_ops;
    std::vector<std::string> name;

    xtal::Lattice lattice = struc.lattice();

    double orig_tol = lattice.tol();
    for(double i = small_tol; i < large_tol; i += increment) {
      tols.push_back(i);
      lattice.set_tol(i);

      xtal::SymOpVector factor_group_operations = xtal::make_factor_group(struc.structure());
      CASM::SymGroup factor_group = adapter::Adapter<SymGroup, xtal::SymOpVector>()(factor_group_operations, lattice);

      factor_group.get_multi_table();
      num_ops.push_back(factor_group.size());
      is_group.push_back(factor_group.is_group(i));
      factor_group.enforce_group(i);
      num_enforced_ops.push_back(factor_group.size());
      factor_group.character_table();
      name.push_back(factor_group.get_name());
    }
    lattice.set_tol(orig_tol);

    for(Index i = 0; i < tols.size(); i++) {
      std::cout << tols[i] << "\t" << num_ops[i] << "\t" << is_group[i] << "\t" << num_enforced_ops[i] << "\t name: " << name[i] << "\n";
    }

    return;
  }
}

namespace CASM {

  namespace Completer {
    SymOption::SymOption(): OptionHandlerBase("sym") {}

    void SymOption::initialize() {
      add_help_suboption();
      add_coordtype_suboption();
      add_configlist_suboption();
      add_confignames_suboption();
      add_scelnames_suboption();
      add_dofs_suboption();

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

    ProjectSettings set = open_project_settings(root);
    DirectoryStructure const &dir = set.dir();
    Structure prim(read_prim(dir.prim(), set.crystallography_tol(), &(set.hamiltonian_modules())));

    args.log() << "Generating lattice point group. " << std::endl << std::endl;
    SymGroup prim_pg(SymGroup::lattice_point_group(prim.lattice()));
    prim_pg.character_table();


    args.log() << "  Lattice point group size: " << prim_pg.size() << std::endl;
    args.log() << "  Lattice point group is: " << prim_pg.get_name() << std::endl << std::endl;

    args.log() << "Generating factor group. " << std::endl << std::endl;

    /* prim.generate_factor_group(); */
    /* prim.set_site_internals(); */

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
      struc = Structure(xtal::make_primitive(struc));

      int biggest = struc.factor_group().size();
      BasicStructure basic_tmp = struc;
      // a) symmetrize the lattice vectors
      Lattice lat = basic_tmp.lattice();
      lat = xtal::symmetrize(lat, tol);
      lat.set_tol(tol);
      basic_tmp.set_lattice(lat, FRAC);

      Structure tmp(basic_tmp);

      tmp.factor_group();
      // b) find factor group with same tolerance
      ::print_factor_group_convergence(tmp, tmp.structure().lattice().tol(), tol, (tol - tmp.structure().lattice().tol()) / 10.0, std::cout);
      // c) symmetrize the basis sites
      SymGroup g = tmp.factor_group();
      tmp = xtal::symmetrize(tmp, g);

      //TODO: Why are we doing this twice?
      g = tmp.factor_group();
      tmp = xtal::symmetrize(tmp, g);
      if(tmp.factor_group().is_group(tol) && (tmp.factor_group().size() > biggest)) {
        struc = tmp;
      }
      struc = Structure(xtal::make_primitive(struc));
      fs::ofstream file_i;
      fs::path POSCARpath_i = "POSCAR_sym";
      file_i.open(POSCARpath_i);
      VaspIO::PrintPOSCAR p_i(xtal::make_simple_structure(struc), struc.structure().title());
      p_i.print(file_i);
      file_i.close();
      return 0;

    }
    coordtype = sym_opt.coordtype_enum();


    // Write symmetry info files
    dir.new_symmetry_dir();

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
