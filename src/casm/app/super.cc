#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/CoordinateSystems.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/ConfigMapping.hh"
#include "casm/clex/ConfigSelection.hh"
#include "casm/crystallography/SupercellEnumerator.hh"
#include "casm/casm_io/VaspIO.hh"
#include "casm/app/casm_functions.hh"
#include "casm/crystallography/Niggli.hh"

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
       "Construct the superdupercell, the minimum supercell of all input supercells "
       "and configurations.")

      ("fixed-orientation",
       "When constructing the superdupercell, do not consider other symmetrically "
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
       "to tile the superdupercell.")

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
      po::store(po::parse_command_line(args.argc, args.argv, super_opt.desc()), vm); // can throw

      if(!vm.count("help") && !vm.count("desc")) {
        if(!vm.count("duper")) {
          if(vm.count("transf-mat") + vm.count("get-transf-mat") != 1) {
            args.err_log << "Error in 'casm super'. Only one of --transf-mat or --get-transf-mat may be chosen." << std::endl;
            return ERR_INVALID_ARG;
          }
          if(configname.size() > 1 || scelname.size() > 1 || tmatfile.size() > 1) {
            args.err_log << "ERROR: more than one --confignames, --scelnames, or --transf-mat argument "
                         "is only allowed for option --duper" << std::endl;
            return ERR_INVALID_ARG;
          }
          if(config_path.size() > 0) {
            args.err_log << "ERROR: the --configs option is only allowed with option --duper" << std::endl;
            return ERR_INVALID_ARG;
          }
        }
      }

      /** --help option
      */
      if(vm.count("help")) {
        args.log << "\n";
        args.log << super_opt.desc() << std::endl;

        return 0;
      }

      if(vm.count("desc")) {
        args.log << "\n";
        args.log << super_opt.desc() << std::endl;
        args.log << "DESCRIPTION" << std::endl;
        args.log << "                                                                      \n" <<
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
                 "  - Makes the superdupercell of the lattices of all inputs            \n" <<
                 "  - Using '--configs' with no arguments is equivalent to '--configs MASTER',\n" <<
                 "    which uses the master config list                                 \n" <<
                 "  - Default applies prim point group ops to try to find minimum volume\n" <<
                 "    superdupercell, disable with '--fixed-orientation'                \n\n";

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
      args.err_log << "ERROR: " << e.what() << std::endl << std::endl;
      args.err_log << super_opt.desc() << std::endl;
      return 1;
    }
    catch(std::exception &e) {
      args.err_log << "Unhandled Exception reached the top of main: "
                   << e.what() << ", application will now exit" << std::endl;
      return 1;

    }

    COORD_MODE C(coordtype);

    // lambda for printing
    auto print = [&](const BasicStructure<Site> &struc) {
      VaspIO::PrintPOSCAR printer(struc);

      if(vm.count("vasp5")) {
        printer.set_atom_names_on();
      }
      else {
        printer.set_atom_names_off();
      }
      printer.set_coord_mode(coordtype);
      printer.print(args.log);
    };



    // -- no casm project necessary for super cell of a POSCAR -------

    // want absolute paths
    for(auto && file : tmatfile) {
      abs_tmatfile.push_back(fs::absolute(file));
    }
    abs_structfile = fs::absolute(structfile);

    if(vm.count("structure") && vm.count("transf-mat")) {

      if(!fs::exists(abs_structfile)) {
        args.log << "ERROR: " << abs_tmatfile[0] << " not found." << std::endl;
        return 1;
      }
      BasicStructure<Site> unitcell(abs_structfile);

      // -- read transf matrix ---
      if(!fs::exists(abs_tmatfile[0])) {
        args.log << "ERROR: " << abs_tmatfile[0] << " not found." << std::endl;
        return 1;
      }
      Eigen::Matrix3i Tm;
      fs::ifstream file(abs_tmatfile[0]);
      file >> Tm;
      file.close();

      auto super = unitcell.create_superstruc(make_supercell(unitcell.lattice(), Tm));
      super.title = std::string("Supercell of ") + unitcell.title;

      print(super);

      return 0;
    }


    const fs::path &root = args.root;
    if(root.empty()) {
      args.err_log.error("No casm project found");
      args.err_log << std::endl;
      return ERR_NO_PROJ;
    }

    // If 'args.primclex', use that, else construct PrimClex in 'uniq_primclex'
    // Then whichever exists, store reference in 'primclex'
    std::unique_ptr<PrimClex> uniq_primclex;
    PrimClex &primclex = make_primclex_if_not(args, uniq_primclex);

    if(vm.count("duper")) {

      // collect all the Lattice to make the superdupercell of
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
          lat[it->string()] = make_supercell(primclex.get_prim().lattice(), T);
        }
      }

      // collect supercells from --scelnames
      if(vm.count("scelnames")) {
        for(auto it = scelname.begin(); it != scelname.end(); ++it) {
          lat[*it] = primclex.get_supercell(*it).get_real_super_lattice();
        }
      }

      // collect configs from --confignames
      if(vm.count("confignames")) {
        for(auto it = configname.begin(); it != configname.end(); ++it) {
          config_lat[*it] = lat[*it] = primclex.configuration(*it).get_supercell().get_real_super_lattice();
        }
      }

      // collect configs from lists via --configs
      if(vm.count("configs")) {

        // MASTER config list if '--configs' only
        if(config_path.size() == 0) {
          ConstConfigSelection selection(primclex);
          for(auto it = selection.selected_config_begin(); it != selection.selected_config_end(); ++it) {
            config_lat[it.name()] = lat[it.name()] = it->get_supercell().get_real_super_lattice();
          }
        }
        // all input config list if '--configs X Y ...'
        else {
          for(auto c_it = config_path.begin(); c_it != config_path.end(); ++c_it) {
            if(c_it->string() == "MASTER") {
              ConstConfigSelection selection(primclex);
              for(auto it = selection.selected_config_begin(); it != selection.selected_config_end(); ++it) {
                config_lat[it.name()] = lat[it.name()] = it->get_supercell().get_real_super_lattice();
              }
            }
            else {
              ConstConfigSelection selection(primclex, fs::absolute(*c_it));
              for(auto it = selection.selected_config_begin(); it != selection.selected_config_end(); ++it) {
                config_lat[it.name()] = lat[it.name()] = it->get_supercell().get_real_super_lattice();
              }
            }
          }
        }

      }

      std::vector<Lattice> lat_only;
      for(auto it = lat.begin(); it != lat.end(); ++it) {
        lat_only.push_back(it->second);
      }

      // create superdupercell
      auto begin = primclex.get_prim().point_group().begin();
      auto end = primclex.get_prim().point_group().end();
      if(vm.count("fixed-orientation")) {
        end = begin;
      }
      Lattice superduper = superdupercell(lat_only.begin(), lat_only.end(), begin, end);

      /// enforce a minimum volume
      if(vm.count("min-volume")) {

        args.log << "  Enforcing minimum volume: " << min_vol;
        if(vm.count("fixed-shape")) {
          args.log << " (with fixed shape)";
        }
        args.log << "\n\n";

        auto prim_lat = primclex.get_prim().lattice();
        const SymGroup &pg = primclex.get_prim().point_group();
        auto T = is_supercell(superduper, prim_lat, TOL).second;

        args.log << "  Superdupercell lattice: \n" << superduper.lat_column_mat() << "\n\n";

        args.log << "    Initial transformation matrix:\n" << T
                 << "\n    (volume = " << T.cast<double>().determinant() << ")\n\n";

        auto M = enforce_min_volume(prim_lat, T, pg, min_vol, vm.count("fixed-shape"));

        superduper = canonical_equivalent_lattice(make_supercell(superduper, M), pg, TOL);

        auto S = is_supercell(superduper, prim_lat, TOL).second;

        args.log << "  Superdupercell lattice: \n" << superduper.lat_column_mat() << "\n\n";

        args.log << "    Transformation matrix, after enforcing mininum volume:\n"
                 << S << "\n    (volume = " << S.cast<double>().determinant() << ")\n\n";

      }

      Index index = primclex.add_supercell(superduper);
      Supercell &superduper_scel = primclex.get_supercell(index);

      args.log << "--- Lattices as column vector matrices ---\n\n";

      args.log << "  Superdupercell: " << primclex.get_supercell(index).get_name() << "\n\n";

      args.log << "  Superdupercell lattice: \n" << superduper.lat_column_mat() << "\n\n";

      args.log << "  Transformation matrix, relative the primitive cell:\n";
      args.log << is_supercell(superduper, primclex.get_prim().lattice(), TOL).second << "\n\n";

      if(vm.count("verbose")) {
        args.log << "Transformation matrices: \n";
        for(auto it = lat.begin(); it != lat.end(); ++it) {
          args.log << "--- \n";
          args.log << "  Unit: " << it->first << ":\n"
                   << it->second.lat_column_mat() << "\n\n";

          auto res = is_supercell(superduper, it->second, begin, end, TOL);
          args.log << "  Superduper = (op*unit) * T\n\nop:\n";
          args.log << res.first->matrix() << "\n\n";
          args.log << "  T:\n";
          args.log << res.second << "\n\n";

        }
        args.log << "--- \n";
      }

      args.log << "  Writing SCEL..." << std::endl;
      primclex.print_supercells();
      args.log << "  DONE\n";


      if(vm.count("add-canonical")) {
        args.log << "Add super configurations:\n";
        for(auto it = config_lat.begin(); it != config_lat.end(); ++it) {
          auto res = is_supercell(superduper, it->second, begin, end, TOL);
          FillSupercell f(superduper_scel, *res.first);
          Index config_index;
          Supercell::permute_const_iterator permute_it;
          bool result = superduper_scel.add_config(
                          f(primclex.configuration(it->first)),
                          config_index,
                          permute_it);
          if(result) {
            args.log << "  " << it->first << "  ->  " << superduper_scel.get_config(config_index).name() << "\n";
          }
        }
        args.log << "\n";
        args.log << "  Writing config_list..." << std::endl << std::endl;
        primclex.write_config_list();
        args.log << "  DONE\n";
      }

      return 0;

    }
    else if(vm.count("transf-mat")) {

      Eigen::Matrix3i T;
      if(!fs::exists(abs_tmatfile[0])) {
        args.log << "ERROR: " << abs_tmatfile[0] << " not found." << std::endl;
        return 1;
      }
      fs::ifstream file(abs_tmatfile[0]);
      file >> T;
      file.close();

      args.log << "Read transformation matrix, T: \n" << T << "\n\n";

      /// enforce a minimum volume
      if(vm.count("min-volume")) {

        if(!vm.count("fixed-shape")) {
          args.log << "  Enforcing minimum volume: \n";
          args.log << "    Finding T' = T*M, such that (T').determinant() >= " << min_vol;
        }
        else {
          args.log << "  Enforcing minimum volume (with fixed shape): \n";
          args.log << "    Finding T' = T*m*I, such that (T').determinant() >= " << min_vol;
        }
        args.log << "\n\n";

        auto prim_lat = primclex.get_prim().lattice();
        const SymGroup &pg = primclex.get_prim().point_group();

        args.log << "    Initial transformation matrix:\n" << T
                 << "\n    (volume = " << T.cast<double>().determinant() << ")\n\n";

        auto M = enforce_min_volume(
                   primclex.get_prim().lattice(),
                   T,
                   pg,
                   min_vol,
                   vm.count("fixed-shape"));

        Lattice niggli_lat = canonical_equivalent_lattice(make_supercell(prim_lat, T * M), pg, TOL);
        T = is_supercell(niggli_lat, prim_lat, TOL).second;

        args.log << "    Transformation matrix, after enforcing mininum volume:\n"
                 << T << "\n    (volume = " << T.cast<double>().determinant() << ")\n\n";
      }


      // super lattice
      if(vm.count("scelnames")) {

        Supercell &scel = primclex.get_supercell(scelname[0]);

        args.log << "  Unit cell: " << scelname[0] << "\n\n";

        args.log << "  Unit cell lattice: \n" << scel.get_real_super_lattice().lat_column_mat() << "\n\n";

        Lattice super_lat = make_supercell(scel.get_real_super_lattice(), T);
        Index index = primclex.add_supercell(super_lat);
        Supercell &super_scel = primclex.get_supercell(index);

        args.log << "  Add supercell: " << super_scel.get_name() << "\n\n";

        args.log << "  Supercell lattice: \n" << super_scel.get_real_super_lattice().lat_column_mat() << "\n\n";

        args.log << "  Transformation matrix: \n" << super_scel.get_transf_mat() << "\n\n";

        args.log << "  Writing SCEL..." << std::endl;
        primclex.print_supercells();
        args.log << "  DONE\n";

      }
      // super structure
      else if(vm.count("confignames")) {

        std::stringstream ss;
        const Configuration &con = primclex.configuration(configname[0]);

        VaspIO::PrintPOSCAR p(con);
        p.sort();
        p.print(ss);

        std::istringstream iss(ss.str());
        BasicStructure<Site> unit;
        unit.read(iss);

        args.log << "Unit structure:";
        args.log << "\n------\n";
        print(unit);
        args.log << "\n------\n";
        args.log << "\n\n";


        BasicStructure<Site> super = unit.create_superstruc(make_supercell(unit.lattice(), T));
        super.title = std::string("Supercell of ") + con.name();

        args.log << "Super structure:";
        args.log << "\n------\n";
        print(super);
        args.log << "\n------\n";

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
            args.log << "  The configuration was imported successfully as "
                     << imported_name << std::endl << std::endl;

          }
          else {
            args.log << "  The configuration was mapped onto pre-existing equivalent structure "
                     << imported_name << std::endl << std::endl;
          }

          jsonParser json_src;
          json_src["supercell_of"] = configname[0];
          primclex.configuration(imported_name).push_back_source(json_src);

          //Update directories
          args.log << "  Writing SCEL..." << std::endl;
          primclex.print_supercells();
          args.log << "  Writing config_list..." << std::endl << std::endl;
          primclex.write_config_list();
          args.log << "  DONE" << std::endl << std::endl;

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

        Eigen::Matrix3d S_niggli = canonical_equivalent_lattice(Lattice(S), pg, tol).lat_column_mat();
        Eigen::Matrix3i T_niggli = iround(U.inverse() * S_niggli);
        Eigen::Matrix3i H_niggli = hermite_normal_form(T_niggli).first;

        std::stringstream s_name;
        s_name << "SCEL" << H_niggli(0, 0)*H_niggli(1, 1)*H_niggli(2, 2) << "_"
               << H_niggli(0, 0) << "_" << H_niggli(1, 1) << "_" << H_niggli(2, 2) << "_"
               << H_niggli(1, 2) << "_" << H_niggli(0, 2) << "_" << H_niggli(0, 1);

        // S = U*T;
        // S_canon = op_canon*S = U*T', where H_canon*V = U.inv*op_canon*U*T = T'

        args.log << "--- Lattices as column vector matrices ---\n\n";

        args.log << "Prim lattice, U:\n" << U << "\n\n";

        args.log << "Super lattice, S = U*T:\n" << S << "\n\n";

        args.log << "This is equivalent to '" << s_name.str() << "', the equivalent super lattice \n" <<
                 "in the standard orientation niggli cell, S_niggli:\n" << S_niggli << "\n\n";

        args.log << "The transformation matrix (S_niggli = U*T) for '" << s_name.str() << "' is:\n" << T_niggli << "\n\n";


        args.log << "--- Lattices as row vector matrices ---\n\n";

        args.log << "Prim lattice:\n" << U.transpose() << "\n\n";

        args.log << "Super lattice:\n" << S.transpose() << "\n\n";

        args.log << "This is equivalent to '" << s_name.str() << "', the equivalent super lattice \n" <<
                 "in the standard orientation niggli cell:\n" << S_niggli.transpose() << "\n\n";

        return 0;
      }

    }

    if(vm.count("get-transf-mat")) {

      Lattice unit_lat = primclex.get_prim().lattice();

      if(vm.count("unitcell")) {
        unit_lat = primclex.get_supercell(unitscelname).get_real_super_lattice();
      }

      Lattice super_lat;

      if(vm.count("structure")) {
        super_lat = BasicStructure<Site>(abs_structfile).lattice();
      }
      else if(vm.count("scelnames")) {
        super_lat = primclex.get_supercell(scelname[0]).get_real_super_lattice();
      }
      else {
        args.log << "Error in 'casm super --get-transf-mat'. No --structure or --scelnames given." << std::endl << std::endl;
        return 1;
      }

      args.log << "--- Lattices as column vector matrices ---\n\n";

      args.log << "Unit lattice, U:\n" << unit_lat.lat_column_mat() << "\n\n";

      args.log << "Super lattice, S:\n" << super_lat.lat_column_mat() << "\n\n";

      // see if super_lat is a supercell of unitlat
      // S == U*T
      Eigen::Matrix3d T = unit_lat.lat_column_mat().inverse() * super_lat.lat_column_mat();

      if(is_integer(T, TOL) && !almost_zero(T, TOL)) {
        args.log << "The super lattice is a supercell of the unit lattice.\n\n";

        args.log << "The transformation matrix, T, where S = U*T, is: \n" << iround(T) << "\n\n";
      }
      else {
        args.log << "The super lattice is NOT a supercell of the unit lattice.\n\n";

        args.log << "The transformation matrix, T, where S = U*T, is: \n" << T << "\n\n";
      }

      return 0;
    }

    return 0;
  };

}


