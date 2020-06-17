#include <boost/filesystem/fstream.hpp>
#include <ostream>
#include "casm/crystallography/BasicStructureTools.hh"
#include "casm/crystallography/CoordinateSystems.hh"
#include "casm/crystallography/Lattice.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/crystallography/SimpleStructureTools.hh"
#include "casm/crystallography/SymTools.hh"

#include "casm/symmetry/json_io.hh"
#include "casm/symmetry/SymRepTools.hh"

#include "casm/enumerator/Enumerator.hh"
#include "casm/enumerator/DoFSpace.hh"
#include "casm/enumerator/io/json/DoFSpace.hh"

#include "casm/clex/PrimClex.hh"
#include "casm/app/ProjectSettings.hh"
#include "casm/app/DirectoryStructure.hh"
#include "casm/app/AppIO.hh"
#include "casm/app/casm_functions.hh"
#include "casm/app/sym.hh"
#include "casm/crystallography/io/VaspIO.hh"


#include "casm/completer/Handlers.hh"
#include "casm/symmetry/SymGroup.hh"

namespace Local {
  using namespace CASM;
  static void _print_factor_group_convergence(const Structure &struc, double small_tol, double large_tol, double increment, std::ostream &print_stream) {
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
      name.push_back(factor_group.get_name());
    }
    lattice.set_tol(orig_tol);

    for(Index i = 0; i < tols.size(); i++) {
      std::cout << tols[i] << "\t" << num_ops[i] << "\t" << is_group[i] << "\t" << num_enforced_ops[i] << "\t name: " << name[i] << "\n";
    }

    return;
  }

  static void _write_config_symmetry_files(ConfigEnumInput const &config, fs::path sym_dir) {
    // Write lattice point group
    {
      SymGroup lattice_pg(SymGroup::lattice_point_group(config.config().ideal_lattice()));
      fs::ofstream outfile;
      jsonParser json;
      outfile.open(sym_dir / "lattice_point_group.json");
      write_symgroup(lattice_pg, json);
      json.print(outfile);
      outfile.close();
    }

    // Write factor group
    {
      fs::ofstream outfile;
      jsonParser json;
      outfile.open(sym_dir / "factor_group.json");
      write_symgroup(make_sym_group(config.group().begin(),
                                    config.group().end(),
                                    config.config().ideal_lattice()), json);
      json.print(outfile);
      outfile.close();
    }

    // Write crystal point group
    {
      fs::ofstream outfile;
      jsonParser json;
      outfile.open(sym_dir / "crystal_point_group.json");
      write_symgroup(make_point_group(config.group().begin(),
                                      config.group().end(),
                                      config.config().ideal_lattice()), json);
      json.print(outfile);
      outfile.close();
    }
  }
}

namespace CASM {

  namespace Completer {
    SymOption::SymOption(): EnumOptionBase("sym") {}

    void SymOption::initialize() {
      bool required = false;

      add_help_suboption();
      add_coordtype_suboption();
      add_selection_suboption("NONE");
      add_confignames_suboption();
      add_scelnames_suboption();
      add_dofs_suboption();
      add_settings_suboption(required);
      add_input_suboption(required);

      m_desc.add_options()
      ("lattice-point-group", "Pretty print lattice point group")
      ("factor-group", "Pretty print factor group")
      ("crystal-point-group", "Pretty print crystal point group")
      ("calc-wedge", "Perform calculation of irreducible wedge (may significantly slow down analysis)")
      //("no-directions", "Skip calculation of high-symmetry direction and irreducible wedge (for faster evaluation)")
      ("tol", po::value<double>(&m_tol)->default_value(1.0e-5), "Tolerance (in Angstr.) used for symmetrization (default 1e-5)")
      ("symmetrize", po::value<fs::path>(&m_poscar_path)->value_name(ArgHandler::path()), "symmetrize a POSCAR specified by path to a given tolerance");

      return;
    }
  }
}

namespace Local {
  static bool _dof_analysis(CASM::Completer::SymOption const &opt) {
    return opt.selection_path() != "NONE" || opt.config_strs().size() || opt.supercell_strs().size() || opt.dof_strs().size();
  }

  static bool _symmetrize(CASM::po::variables_map const &vm) {
    return vm.count("symmetrize");
  }

}

namespace CASM {
  // ///////////////////////////////////////
  // 'sym' function for casm
  //    (add an 'if-else' statement in casm.cpp to call this)

  const std::string SymCommand::name = "sym";

  SymCommand::SymCommand(const CommandArgs &_args, Completer::SymOption &_opt) :
    APICommand<Completer::SymOption>(_args, _opt) {}


  int SymCommand::vm_count_check() const {
    if(!Local::_symmetrize(vm()) && !in_project()) {
      help();
      err_log().error("No casm project found");
      err_log() << std::endl;
      return ERR_NO_PROJ;
    }

    if(vm().count("settings") + vm().count("input") == 2) {
      help();
      err_log() << "Error in 'casm sym'. The options --settings or --input may not both be chosen." << std::endl;
      return ERR_INVALID_ARG;
    }

    return 0;
  }

  int SymCommand::help() const {
    log() << "\n";
    log() << opt().desc() << std::endl;

    return 0;
  }

  int SymCommand::desc() const {
    log() << "\n";
    log() << opt().desc() << std::endl;
    log() << "DESCRIPTION" << std::endl;
    log() << "    Display symmetry group information.\n";
    return 0;
  }

  int SymCommand::run() const {

    if(Local::_symmetrize(vm())) {
      fs::path poscar_path = opt().poscar_path();
      double tol = opt().tol();
      log() << "\n***************************\n" << std::endl;
      log() << "Symmetrizing: " << poscar_path << std::endl;
      log() << "with tolerance: " << tol << std::endl;
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
      Local::_print_factor_group_convergence(tmp, tmp.structure().lattice().tol(), tol, (tol - tmp.structure().lattice().tol()) / 10.0, std::cout);
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

    //std::string name;
    COORD_TYPE coordtype;

    coordtype = opt().coordtype_enum();
    COORD_MODE C(coordtype);

    Structure const &prim = primclex().prim();

    SymGroup lattice_pg(SymGroup::lattice_point_group(prim.lattice()));

    if(!Local::_dof_analysis(opt())) {
      log() << "  Lattice point group size: " << lattice_pg.size() << std::endl;
      log() << "  Lattice point group is: " << lattice_pg.get_name() << std::endl << std::endl;


      log() << "  Factor group size: " << prim.factor_group().size() << std::endl;

      log() << "  Crystal point group is: " << prim.point_group().get_name() << std::endl;


      if(vm().count("lattice-point-group")) {
        log() << "\n***************************\n" << std::endl;
        log() << "Lattice point group:\n\n" << std::endl;
        lattice_pg.print(log(), coordtype);
      }

      if(vm().count("factor-group")) {
        log() << "\n***************************\n" << std::endl;
        log() << "Factor group:\n\n" << std::endl;
        prim.factor_group().print(log(), coordtype);
      }

      if(vm().count("crystal-point-group")) {
        log() << "\n***************************\n" << std::endl;
        log() << "Crystal point group:\n\n" << std::endl;
        prim.point_group().print(log(), coordtype);
      }
    }

    coordtype = opt().coordtype_enum();


    // Write symmetry info files
    primclex().settings().new_symmetry_dir();

    // Write lattice point group
    {
      fs::ofstream outfile;
      jsonParser json;
      outfile.open(primclex().dir().lattice_point_group());
      write_symgroup(lattice_pg, json);
      json.print(outfile);
      outfile.close();
    }

    // Write factor group
    {
      fs::ofstream outfile;
      jsonParser json;
      outfile.open(primclex().dir().factor_group());
      write_symgroup(prim.factor_group(), json);
      json.print(outfile);
      outfile.close();
    }

    // Write crystal point group
    {
      fs::ofstream outfile;
      jsonParser json;
      outfile.open(primclex().dir().crystal_point_group());
      write_symgroup(prim.point_group(), json);
      json.print(outfile);
      outfile.close();
    }

    // Perform DoF Analysis for specified degrees of freedom (DoFs)
    if(Local::_dof_analysis(opt())) {
      jsonParser kwargs;
      if(vm().count("settings")) {
        kwargs = jsonParser {opt().settings_path()};
      }
      else if(vm().count("input")) {
        kwargs = jsonParser::parse(opt().input_str());
      }
      else {
        kwargs = jsonParser::parse(std::string("{\"supercells\" : {\"min\" : 0, \"max\" : 0}}"));
      }

      // Create supercells, configurations, and/or local environments
      std::vector<ConfigEnumInput> configs = make_enumerator_input_configs(primclex(), kwargs, opt(), nullptr);

      // For each enumeration envrionment, perform analysis and write files.
      for(ConfigEnumInput const &config : configs) {
        fs::path sym_dir = primclex().dir().symmetry_dir(config.name());
        fs::create_directories(sym_dir);

        Local::_write_config_symmetry_files(config, sym_dir);

        std::vector<DoFKey> dofs = opt().dof_strs();
        if(dofs.empty()) {
          dofs = all_local_dof_types(primclex().prim());
          for(DoFKey const &dof : global_dof_types(primclex().prim())) {
            dofs.push_back(dof);
          }
        }
        for(DoFKey const &dof : dofs) {
          DoFSpace dspace(config, dof);
          jsonParser report;
          std::string filename = "dof_analysis_" + dof + ".json";
          if(fs::is_regular_file(sym_dir / filename)) {
            report.read(sym_dir / filename);
          }

          report = vector_space_sym_report(dspace,
                                           vm().count("calc-wedge"));

          to_json(dspace, report);

          report.write(sym_dir / filename);

        }

      }

    }
    log() << std::endl;

    return 0;

  };

}
