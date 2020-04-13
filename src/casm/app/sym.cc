#include <boost/filesystem/fstream.hpp>
#include "casm/crystallography/CoordinateSystems.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/crystallography/SimpleStructureTools.hh"
#include "casm/crystallography/SymTools.hh"
#include "casm/symmetry/json_io.hh"
#include "casm/symmetry/SymRepTools.hh"

#include "casm/enumerator/Enumerator.hh"
#include "casm/enumerator/DoFSpace.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/app/ProjectSettings.hh"
#include "casm/app/DirectoryStructure.hh"
#include "casm/app/AppIO.hh"
#include "casm/app/casm_functions.hh"
#include "casm/app/sym.hh"
#include "casm/crystallography/io/VaspIO.hh"


#include "casm/completer/Handlers.hh"

namespace CASM {

  namespace Completer {
    SymOption::SymOption(): EnumOptionBase("sym") {}

    void SymOption::initialize() {
      bool required = false;

      add_help_suboption();
      add_coordtype_suboption();
      add_selection_suboption();
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
    return opt.selection_path().size() || opt.config_strs().size() || opt.supercell_strs().size() || opt.dof_strs().size();
  }

  static bool _prim_analysis(CASM::po::variables_map const &vm) {
    return vm.count("lattice-point-group") || vm.count("factor-group") || vm.count("crystal-point-group");

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
    if(!vm().count("symmetrize") && !in_project()) {
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
    //std::string name;
    COORD_TYPE coordtype;

    coordtype = opt().coordtype_enum();
    COORD_MODE C(coordtype);

    Structure const &prim = primclex().prim();

    SymGroup lattice_pg(SymGroup::lattice_point_group(prim.lattice()));
    lattice_pg.character_table();

    if(Local::_prim_analysis(vm())) {
      log() << "  Lattice point group size: " << lattice_pg.size() << std::endl;
      log() << "  Lattice point group is: " << lattice_pg.get_name() << std::endl << std::endl;

      log() << "Generating factor group. " << std::endl << std::endl;

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

    if(vm().count("symmetrize")) {
      fs::path poscar_path = opt().poscar_path();
      double tol = opt().tol();
      log() << "\n***************************\n" << std::endl;
      log() << "Symmetrizing: " << poscar_path << std::endl;
      log() << "with tolerance: " << tol << std::endl;
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

      std::vector<ConfigEnumInput> configs = make_enumerator_input_configs(primclex(), kwargs, opt(), nullptr);

      for(ConfigEnumInput const &config : configs) {
        std::cout << "WORKING ON CONFIG: " << config.name() << std::endl;
        jsonParser report;
        fs::path sym_dir = primclex().dir().symmetry_dir(config.name());

        if(fs::is_regular_file(sym_dir / "report.json"))
          report.read(sym_dir / "report.json");

        std::vector<DoFKey> dofs = opt().dof_strs();
        if(dofs.empty()) {
          dofs = all_local_dof_types(primclex().prim());
          for(DoFKey const &dof : global_dof_types(primclex().prim())) {
            dofs.push_back(dof);
          }
        }
        for(DoFKey const &dof : dofs) {
          DoFSpace dspace(config, dof);
          report[dof] = vector_space_sym_report(dspace,
                                                vm().count("calc-wedge"));
        }

        fs::create_directories(sym_dir);
        report.write(sym_dir / "report.json");
      }

    }
    log() << std::endl;

    return 0;

  };

}
