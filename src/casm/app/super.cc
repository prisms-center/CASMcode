#include <boost/filesystem/fstream.hpp>
#include <fstream>
#include "casm/crystallography/BasicStructureTools.hh"
#include "casm/external/Eigen/Core"
#include "casm/external/Eigen/src/Core/Matrix.h"
#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/CoordinateSystems.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/crystallography/SimpleStructureTools.hh"
#include "casm/crystallography/SymType.hh"
#include "casm/crystallography/SuperlatticeEnumerator.hh"
#include "casm/crystallography/CanonicalForm.hh"
#include "casm/crystallography/Adapter.hh"
#include "casm/crystallography/SymTools.hh"
#include "casm/clex/FillSupercell.hh"
#include "casm/clex/Supercell_impl.hh"
#include "casm/clex/Configuration_impl.hh"
#include "casm/clex/ConfigMapping.hh"
#include "casm/clex/SimpleStructureTools.hh"
#include "casm/database/Selection_impl.hh"
#include "casm/database/DatabaseTypes.hh"
#include "casm/database/ScelDatabase.hh"
#include "casm/database/ConfigDatabase.hh"
#include "casm/crystallography/io/VaspIO.hh"
#include "casm/app/casm_functions.hh"

#include "casm/completer/Handlers.hh"

namespace CASM {

  namespace Completer {
    SuperOption::SuperOption(): OptionHandlerBase("super") {}

    const std::vector<fs::path> &SuperOption::transf_mat_paths() const {
      return m_transf_mat_paths;
    }

    const fs::path &SuperOption::struct_path() const {
      return m_struct_path;
    }

    const std::string &SuperOption::unit_scel_str() const {
      return m_unit_scel_str;
    }

    Index SuperOption::min_vol() const {
      return m_min_vol;
    }

    double SuperOption::tolerance() const {
      return m_tolerance;
    }

    void SuperOption::initialize() {
      add_help_suboption();
      add_confignames_suboption();
      add_scelnames_suboption();
      add_configlists_nodefault_suboption();
      add_coordtype_suboption();

      m_desc.add_options()
      ("transf-mat",
       po::value<std::vector<fs::path> >(&m_transf_mat_paths)->multitoken()->value_name(ArgHandler::path()),
       "1 or more files containing a 3x3 transformation matrix used to create a supercell.")

      ("get-transf-mat",
       "If it exists, find the transformation matrix.")

      ("structure",
       po::value<fs::path>(&m_struct_path)->value_name(ArgHandler::path()),
       "File with structure (POSCAR type) to use.")

      ("unitcell",
       po::value<std::string>(&m_unit_scel_str)->value_name(ArgHandler::supercell()),
       "Name of supercell to use as unit cell. For ex. 'SCEL2_2_1_1_0_0_0'.")

      ("duper",
       "Construct the superduperlattice, the minimum supercell of all input supercells "
       "and configurations.")

      ("fixed-orientation",
       "When constructing the superduperlattice, do not consider other symmetrically "
       "equivalent orientations.")

      ("min-volume",
       po::value<Index>(&m_min_vol),
       "Transforms the transformation matrix, T -> T', where T' = T*M, such that "
       "(T').determininant() >= V. This has the effect that a supercell has a "
       "particular volume.")

      ("fixed-shape",
       "Used with --min-volume to enforce that T' = T*m*I, where I is the identity "
       "matrix, and m is a scalar. This has the effect of preserving the shape "
       "of the resulting supercell, but increasing the volume.")

      ("verbose",
       "When used with --duper, show how the input lattices are transformed "
       "to tile the superduperlattice.")

      ("add-canonical,a", "Will add the generated super configuration in it's "
       "canonical form in the equivalent niggli supercell.")

      ("vasp5",
       "Print using VASP5 style (include atom name line)")

      ("tol",
       po::value<double>(&m_tolerance)->default_value(CASM::TOL),
       "Tolerance used for checking symmetry");

      return;
    }
  }


  // ///////////////////////////////////////
  // 'super' function for casm
  //    (add an 'if-else' statement in casm.cpp to call this)

  int super_command(const CommandArgs &args) {

    //casm enum [—supercell min max] [—config supercell ] [—hopconfigs hop.background]
    //- enumerate supercells and configs and hop local configurations

    std::vector<std::string> scelname, configname;
    std::string unitscelname;
    fs::path structfile;
    std::vector<fs::path> tmatfile, abs_tmatfile, config_path;
    fs::path abs_structfile;
    Index min_vol;
    double tol;
    COORD_TYPE coordtype;
    po::variables_map vm;

    /// Set command line options using boost program_options
    Completer::SuperOption super_opt;

    try {
      po::store(po::parse_command_line(args.argc(), args.argv(), super_opt.desc()), vm); // can throw

      if(!vm.count("help") && !vm.count("desc")) {
        if(!vm.count("duper")) {
          if(vm.count("transf-mat") + vm.count("get-transf-mat") != 1) {
            err_log() << "Error in 'casm super'. Only one of --transf-mat or --get-transf-mat may be chosen." << std::endl;
            return ERR_INVALID_ARG;
          }
          if(configname.size() > 1 || scelname.size() > 1 || tmatfile.size() > 1) {
            err_log() << "ERROR: more than one --confignames, --scelnames, or --transf-mat argument "
                      "is only allowed for option --duper" << std::endl;
            return ERR_INVALID_ARG;
          }
          if(config_path.size() > 0) {
            err_log() << "ERROR: the --configs option is only allowed with option --duper" << std::endl;
            return ERR_INVALID_ARG;
          }
        }
      }

      /** --help option
      */
      if(vm.count("help")) {
        log() << "\n";
        log() << super_opt.desc() << std::endl;

        return 0;
      }

      if(vm.count("desc")) {
        log() << "\n";
        log() << super_opt.desc() << std::endl;
        log() << "DESCRIPTION" << std::endl;
        log() << "                                                                      \n" <<
              "  casm super --transf-mat T                                           \n" <<
              "  - Print super lattice of the PRIM lattice                           \n" <<
              "                                                                      \n" <<
              "  casm super --structure POSCAR --transf-mat T                        \n" <<
              "  - Print superstructure of a POSCAR                                  \n" <<
              "                                                                      \n" <<
              "  casm super --confignames configname --transf-mat T                   \n" <<
              "  - Print superstructure of a configuration                           \n" <<
              "                                                                      \n" <<
              "  casm super --structure POSCAR --unitcell scelname --get-transf-mat  \n" <<
              "  - Check if POSCAR lattice is a supercell of unit cell lattice and   \n" <<
              "    if so print the transformation matrix                             \n" <<
              "  - Uses primitive cell for unitcell if none given                    \n" <<
              "                                                                      \n" <<
              "  casm super --scelnames scelname --unitcell scelname --get-transf-mat\n" <<
              "  - Check if configuration lattice is a supercell of unit cell lattice.\n" <<
              "    and print the transformation matrix                               \n" <<
              "  - Uses primitive cell for unitcell if none given                    \n\n" <<

              "  casm super --duper --scelnames scel1 [scel2 ...] --confignames con1 [con2 ...]\n"
              "    --configs [mylist ...] --transf-mat M1 [M2 ...]                    \n" <<
              "  - Makes the superduperlattice of the lattices of all inputs            \n" <<
              "  - Using '--configs' with no arguments is equivalent to '--configs MASTER',\n" <<
              "    which uses the master config list                                 \n" <<
              "  - Default applies prim point group ops to try to find minimum volume\n" <<
              "    superduperlattice, disable with '--fixed-orientation'                \n\n";

        return 0;
      }

      po::notify(vm); // throws on error, so do after help in case
      // there are any problems

      scelname = super_opt.supercell_strs();
      configname = super_opt.config_strs();
      unitscelname = super_opt.unit_scel_str();
      structfile = super_opt.struct_path();
      tmatfile = super_opt.transf_mat_paths();
      config_path = super_opt.selection_paths();
      min_vol = super_opt.min_vol();
      tol = super_opt.tolerance();
      coordtype = super_opt.coordtype_enum();
    }
    catch(po::error &e) {
      err_log() << "ERROR: " << e.what() << std::endl << std::endl;
      err_log() << super_opt.desc() << std::endl;
      return 1;
    }
    catch(std::exception &e) {
      err_log() << "Unhandled Exception reached the top of main: "
                << e.what() << ", application will now exit" << std::endl;
      return 1;

    }

    xtal::COORD_MODE C(coordtype);

    // lambda for printing
    auto print = [&](const BasicStructure & struc) {
      VaspIO::PrintPOSCAR printer(make_simple_structure(struc), struc.title());

      if(vm.count("vasp5")) {
        printer.set_atom_names_on();
      }
      else {
        printer.set_atom_names_off();
      }
      printer.set_coord_mode(coordtype);
      printer.print(log());
    };



    // -- no casm project necessary for super cell of a POSCAR -------

    // want absolute paths
    for(auto &&file : tmatfile) {
      abs_tmatfile.push_back(fs::absolute(file));
    }
    abs_structfile = fs::absolute(structfile);

    if(vm.count("structure") && vm.count("transf-mat")) {

      if(!fs::exists(abs_structfile)) {
        log() << "ERROR: " << abs_tmatfile[0] << " not found." << std::endl;
        return 1;
      }
      std::ifstream abs_structfile_stream(abs_structfile.string());
      BasicStructure unitcell = BasicStructure::from_poscar_stream(abs_structfile_stream);

      // -- read transf matrix ---
      if(!fs::exists(abs_tmatfile[0])) {
        log() << "ERROR: " << abs_tmatfile[0] << " not found." << std::endl;
        return 1;
      }
      Eigen::Matrix3i Tm;
      fs::ifstream file(abs_tmatfile[0]);
      file >> Tm;
      file.close();

      BasicStructure super = xtal::make_superstructure(unitcell, Tm);
      super.set_title(std::string("Supercell of ") + unitcell.title());

      print(super);

      return 0;
    }


    const fs::path &root = args.root;
    if(root.empty()) {
      err_log().error("No casm project found");
      err_log() << std::endl;
      return ERR_NO_PROJ;
    }

    // If 'args.primclex', use that, else construct PrimClex in 'uniq_primclex'
    // Then whichever exists, store reference in 'primclex'
    std::unique_ptr<PrimClex> uniq_primclex;
    PrimClex &primclex = make_primclex_if_not(args, uniq_primclex);

    if(vm.count("duper")) {

      // collect all the Lattice to make the superduperlattice of
      std::map<std::string, Lattice> lat;
      std::map<std::string, Lattice> config_lat;

      // collect lattices by constructing from transformation matrices
      if(vm.count("transf-mat")) {
        for(auto it = abs_tmatfile.begin(); it != abs_tmatfile.end(); ++it) {
          Eigen::Matrix<int, 3, 3, Eigen::RowMajor> T;
          fs::ifstream file(*it);
          for(int i = 0; i < 9; i++) {
            file >> T.data()[i];
          }
          file.close();
          lat[it->string()] = make_superlattice(primclex.prim().lattice(), T);
        }
      }

      // collect supercells from --scelnames
      if(vm.count("scelnames")) {
        for(auto it = scelname.begin(); it != scelname.end(); ++it) {
          lat[*it] = primclex.db<Supercell>().find(*it)->lattice();
        }
      }

      // collect configs from --confignames
      if(vm.count("confignames")) {
        for(auto it = configname.begin(); it != configname.end(); ++it) {
          config_lat[*it] = lat[*it] = primclex.db<Configuration>().find(*it)->supercell().lattice();
        }
      }

      // collect configs from lists via --configs
      if(vm.count("configs")) {

        // MASTER config list if '--configs' only
        if(config_path.size() == 0) {
          config_path.push_back("MASTER");
        }

        for(const auto &p : config_path) {
          DB::Selection<Configuration> selection(primclex, p);
          for(const auto &config : selection.selected()) {
            config_lat[config.name()] = lat[config.name()] = config.supercell().lattice();
          }
        }

      }

      std::vector<Lattice> lat_only;
      for(auto it = lat.begin(); it != lat.end(); ++it) {
        lat_only.push_back(it->second);
      }

      // create superduperlattice
      auto begin = primclex.prim().point_group().begin();
      auto end = primclex.prim().point_group().end();
      if(vm.count("fixed-orientation")) {
        end = begin;
      }
      Lattice superduper = xtal::make_equivalent_superduperlattice(lat_only.begin(), lat_only.end(), begin, end);

      /// enforce a minimum volume
      if(vm.count("min-volume")) {

        log() << "  Enforcing minimum volume: " << min_vol;
        if(vm.count("fixed-shape")) {
          log() << " (with fixed shape)";
        }
        log() << "\n\n";

        auto prim_lat = primclex.prim().lattice();
        const SymGroup &pg = primclex.prim().point_group();
        auto T = is_superlattice(superduper, prim_lat, TOL).second;

        log() << "  superduperlattice lattice: \n" << superduper.lat_column_mat() << "\n\n";

        log() << "    Initial transformation matrix:\n" << iround(T)
              << "\n    (volume = " << iround(T).cast<double>().determinant() << ")\n\n";

        auto M = enforce_min_volume(pg.begin(), pg.end(), prim_lat, iround(T),  min_vol, vm.count("fixed-shape"));

        superduper = xtal::canonical::equivalent(make_superlattice(superduper, M), pg, TOL);

        auto S = is_superlattice(superduper, prim_lat, TOL).second;

        log() << "  superduperlattice lattice: \n" << superduper.lat_column_mat() << "\n\n";

        log() << "    Transformation matrix, after enforcing mininum volume:\n"
              << iround(S) << "\n    (volume = " << iround(S).cast<double>().determinant() << ")\n\n";

      }

      const Supercell &superduper_scel = *Supercell(&primclex, superduper).insert().first;
      superduper = superduper_scel.lattice();

      log() << "--- Lattices as column vector matrices ---\n\n";

      log() << "  superduperlattice: " << superduper_scel.name() << "\n\n";

      log() << "  superduperlattice lattice: \n" << superduper.lat_column_mat() << "\n\n";

      log() << "  Transformation matrix, relative the primitive cell:\n";
      log() << iround(is_superlattice(superduper, primclex.prim().lattice(), TOL).second) << "\n\n";

      if(vm.count("verbose")) {
        log() << "Transformation matrices: \n";
        for(auto it = lat.begin(); it != lat.end(); ++it) {
          log() << "--- \n";
          log() << "  Unit: " << it->first << ":\n"
                << it->second.lat_column_mat() << "\n\n";

          auto res = xtal::is_equivalent_superlattice(superduper, it->second, begin, end, TOL);
          log() << "  Superduper = (op*unit) * T\n\nop:\n";
          log() << res.first->matrix() << "\n\n";
          log() << "  T:\n";
          log() << iround(res.second) << "\n\n";

        }
        log() << "--- \n";
      }

      log() << "Write supercell database..." << std::endl;
      primclex.db<Supercell>().commit();
      log() << "  DONE" << std::endl << std::endl;


      if(vm.count("add-canonical")) {
        log() << "Add super configurations:\n";
        for(auto it = config_lat.begin(); it != config_lat.end(); ++it) {
          FillSupercell f {superduper_scel};
          auto insert_res = f(*primclex.db<Configuration>().find(it->first)).insert();
          if(insert_res.insert_canonical) {
            log() << "  " << it->first << "  ->  " << insert_res.canonical_it->name() << "\n";
          }
        }
        log() << "\n";
        log() << "Writing configuration database..." << std::endl;
        primclex.db<Configuration>().commit();
        log() << "  DONE" << std::endl;
      }

      return 0;

    }
    else if(vm.count("transf-mat")) {

      Eigen::Matrix3i T;
      if(!fs::exists(abs_tmatfile[0])) {
        log() << "ERROR: " << abs_tmatfile[0] << " not found." << std::endl;
        return 1;
      }
      fs::ifstream file(abs_tmatfile[0]);
      file >> T;
      file.close();

      log() << "Read transformation matrix, T: \n" << T << "\n\n";

      /// enforce a minimum volume
      if(vm.count("min-volume")) {

        if(!vm.count("fixed-shape")) {
          log() << "  Enforcing minimum volume: \n";
          log() << "    Finding T' = T*M, such that (T').determinant() >= " << min_vol;
        }
        else {
          log() << "  Enforcing minimum volume (with fixed shape): \n";
          log() << "    Finding T' = T*m*I, such that (T').determinant() >= " << min_vol;
        }
        log() << "\n\n";

        auto prim_lat = primclex.prim().lattice();
        const SymGroup &pg = primclex.prim().point_group();

        log() << "    Initial transformation matrix:\n" << T
              << "\n    (volume = " << T.cast<double>().determinant() << ")\n\n";

        auto M = enforce_min_volume(
                   pg.begin(),
                   pg.end(),
                   primclex.prim().lattice(),
                   T,
                   min_vol,
                   vm.count("fixed-shape"));

        Eigen::Matrix3i super_lat_matrix = T * M;
        Lattice niggli_lat = xtal::canonical::equivalent(make_superlattice(prim_lat, super_lat_matrix), pg, TOL);
        T = iround(is_superlattice(niggli_lat, prim_lat, TOL).second);

        log() << "    Transformation matrix, after enforcing mininum volume:\n"
              << T << "\n    (volume = " << iround(T).cast<double>().determinant() << ")\n\n";
      }


      // super lattice
      if(vm.count("scelnames")) {

        const Supercell &scel = *primclex.db<Supercell>().find(scelname[0]);

        log() << "  Unit cell: " << scelname[0] << "\n\n";

        log() << "  Unit cell lattice: \n" << scel.lattice().lat_column_mat() << "\n\n";

        Lattice super_lat = make_superlattice(scel.lattice(), T);
        const Supercell &super_scel = *Supercell(&primclex, super_lat).insert().first;

        log() << "  Add supercell: " << super_scel.name() << "\n\n";

        log() << "  Supercell lattice: \n" << super_scel.lattice().lat_column_mat() << "\n\n";

        log() << "  Transformation matrix: \n" << super_scel.transf_mat() << "\n\n";

        log() << "Write supercell database..." << std::endl;
        primclex.db<Supercell>().commit();
        log() << "  DONE" << std::endl << std::endl;


      }
      // super structure
      else if(vm.count("confignames")) {

        std::stringstream ss;
        const Configuration &con = *primclex.db<Configuration>().find(configname[0]);

        VaspIO::PrintPOSCAR p(make_simple_structure(con), con.name());
        p.sort();
        p.print(ss);

        std::istringstream iss(ss.str());
        BasicStructure unit = BasicStructure::from_poscar_stream(iss);

        log() << "Unit structure:";
        log() << "\n------\n";
        print(unit);
        log() << "\n------\n";
        log() << "\n\n";

        //TODO: Why is there so much code copied and pasted in all the casm commands?
        BasicStructure super = xtal::make_superstructure(unit, T);
        super.set_title(std::string("Supercell of ") + con.name());

        log() << "Super structure:";
        log() << "\n------\n";
        print(super);
        log() << "\n------\n";

        if(vm.count("add-canonical")) {

          double tol = TOL;

          ConfigMapper configmapper(primclex, ConfigMapping::Settings(), tol);

          auto map_res = configmapper.import_structure(make_simple_structure(super));

          if(map_res.success()) {
            auto insert_res = (map_res.maps.begin()->second).config.insert();
            Configuration imported_config = *insert_res.canonical_it;

            if(insert_res.insert_canonical) {
              log() << "  The configuration was imported successfully as "
                    << imported_config.name() << std::endl << std::endl;

            }
            else {
              log() << "  The configuration was mapped onto pre-existing equivalent structure "
                    << imported_config.name() << std::endl << std::endl;
            }
            jsonParser json_src;
            json_src["supercell_of"] = configname[0];
            imported_config.push_back_source(json_src);
            primclex.db<Configuration>().update(imported_config);
          }
          else {
            log() << "  Failure to map configuration." << std::endl << std::endl;
          }

          //Update directories
          log() << "Write supercell database..." << std::endl;
          primclex.db<Supercell>().commit();
          log() << "  DONE" << std::endl << std::endl;

          log() << "Writing configuration database..." << std::endl;
          primclex.db<Configuration>().commit();
          log() << "  DONE" << std::endl;

        }

        return 0;

      }
      // super lattice of prim lattice
      else {

        BasicStructure unit = primclex.prim();
        SymGroup pg = Structure(unit).point_group();

        Eigen::Matrix3d U = unit.lattice().lat_column_mat();
        Eigen::Matrix3d S = U * T.cast<double>();
        Eigen::Matrix3i H_canon;
        Eigen::Matrix3d op_canon;

        Eigen::Matrix3d S_niggli = xtal::canonical::equivalent(Lattice(S), pg, tol).lat_column_mat();
        Eigen::Matrix3i T_niggli = iround(U.inverse() * S_niggli);
        Eigen::Matrix3i H_niggli = hermite_normal_form(T_niggli).first;

        std::stringstream s_name;
        s_name << "SCEL" << H_niggli(0, 0)*H_niggli(1, 1)*H_niggli(2, 2) << "_"
               << H_niggli(0, 0) << "_" << H_niggli(1, 1) << "_" << H_niggli(2, 2) << "_"
               << H_niggli(1, 2) << "_" << H_niggli(0, 2) << "_" << H_niggli(0, 1);

        // S = U*T;
        // S_canon = op_canon*S = U*T', where H_canon*V = U.inv*op_canon*U*T = T'

        log() << "--- Lattices as column vector matrices ---\n\n";

        log() << "Prim lattice, U:\n" << U << "\n\n";

        log() << "Super lattice, S = U*T:\n" << S << "\n\n";

        log() << "This is equivalent to '" << s_name.str() << "', the equivalent super lattice \n" <<
              "in the standard orientation niggli cell, S_niggli:\n" << S_niggli << "\n\n";

        log() << "The transformation matrix (S_niggli = U*T) for '" << s_name.str() << "' is:\n" << T_niggli << "\n\n";


        log() << "--- Lattices as row vector matrices ---\n\n";

        log() << "Prim lattice:\n" << U.transpose() << "\n\n";

        log() << "Super lattice:\n" << S.transpose() << "\n\n";

        log() << "This is equivalent to '" << s_name.str() << "', the equivalent super lattice \n" <<
              "in the standard orientation niggli cell:\n" << S_niggli.transpose() << "\n\n";

        return 0;
      }

    }

    if(vm.count("get-transf-mat")) {

      Lattice unit_lat = primclex.prim().lattice();

      if(vm.count("unitcell")) {
        unit_lat = primclex.db<Supercell>().find(unitscelname)->lattice();
      }

      Lattice super_lat;

      if(vm.count("structure")) {
        std::ifstream abs_structfile_stream(abs_structfile.string());
        super_lat = BasicStructure::from_poscar_stream(abs_structfile_stream).lattice();
      }
      else if(vm.count("scelnames")) {
        super_lat = primclex.db<Supercell>().find(unitscelname)->lattice();
      }
      else {
        log() << "Error in 'casm super --get-transf-mat'. No --structure or --scelnames given." << std::endl << std::endl;
        return 1;
      }

      log() << "--- Lattices as column vector matrices ---\n\n";

      log() << "Unit lattice, U:\n" << unit_lat.lat_column_mat() << "\n\n";

      log() << "Super lattice, S:\n" << super_lat.lat_column_mat() << "\n\n";

      // see if super_lat is a supercell of unitlat
      // S == U*T
      Eigen::Matrix3d T;
      bool is_superlat;
      std::tie(is_superlat, T) = xtal::is_superlattice(super_lat, unit_lat, super_lat.tol());
      if(is_superlat) {
        log() << "The super lattice is a supercell of the unit lattice.\n\n";

        log() << "The transformation matrix, T, where S = U*T, is: \n" << iround(T) << "\n\n";
      }
      else {
        log() << "The super lattice is NOT a supercell of the unit lattice.\n\n";

        log() << "The transformation matrix, T, where S = U*T, is: \n" << T << "\n\n";
      }

      return 0;
    }

    return 0;
  };

}
