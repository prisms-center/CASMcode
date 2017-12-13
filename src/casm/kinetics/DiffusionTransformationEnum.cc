#include "casm/crystallography/Site.hh"
#include "casm/clusterography/ClusterOrbits.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/app/AppIO.hh"
#include "casm/database/DiffTransOrbitDatabase.hh"

#include "casm/kinetics/DiffusionTransformationEnum_impl.hh"
#include "casm/clusterography/ClusterSymCompare_impl.hh"
#include "casm/app/AppIO_impl.hh"

extern "C" {
  CASM::EnumInterfaceBase *make_DiffusionTransformationEnum_interface() {
    return new CASM::EnumInterface<CASM::Kinetics::DiffusionTransformationEnum>();
  }
}

namespace CASM {

  namespace Kinetics {

    namespace {

      Log &operator<<(Log &out, const IntegralCluster &clust) {
        SitesPrinter printer;
        printer.print(clust, out);
        out << std::endl;
        return out;
      }
    }

    DiffTransEnumParser::DiffTransEnumParser(
      const PrimClex &_primclex,
      jsonParser &_input,
      fs::path _path,
      bool _required) :
      DiffTransEnumParser(_primclex, _input, Completer::EnumOption(), _path, _required) {}

    DiffTransEnumParser::DiffTransEnumParser(
      const PrimClex &_primclex,
      jsonParser &_input,
      const Completer::EnumOption &_enum_opt,
      fs::path _path,
      bool _required) :
      EnumInputParser(_primclex, _input, _enum_opt, _path, _required) {

      if(exists()) {
        auto _relpath = Relpath(_path);

        // check "cspecs"
        PrimPeriodicSymCompare<IntegralCluster> sym_compare(_primclex);
        this->kwargs["cspecs"] = m_cspecs_parser =
                                   std::make_shared<PrimPeriodicClustersByMaxLength>(
                                     _primclex, _primclex.prim().factor_group(), sym_compare, input, _relpath("cspecs"), true);

        // check "require" and "exclude"
        this->kwargs["require"] = m_require =
                                    std::make_shared<SpeciesSetParser>(
                                      _primclex, ALLOWED_SPECIES_TYPES::ALL, "require", input, _relpath("require"), true);

        this->kwargs["exclude"] = m_exclude =
                                    std::make_shared<SpeciesSetParser>(
                                      _primclex, ALLOWED_SPECIES_TYPES::ALL, "exclude", input, _relpath("exclude"), true);

        // check values for optional dry_run, coordinate_mode, orbit_print_mode
        optional<bool>("dry_run");
        optional<COORD_TYPE>(traits<COORD_TYPE>::name);
        optional<ORBIT_PRINT_MODE>(traits<ORBIT_PRINT_MODE>::name);

        warn_unnecessary(expected());
      }
    }

    std::set<std::string> DiffTransEnumParser::required_species() const {
      return m_require->values();
    }

    std::set<std::string> DiffTransEnumParser::excluded_species() const {
      return m_exclude->values();
    }

    const PrimPeriodicClustersByMaxLength &DiffTransEnumParser::cspecs() const {
      return *m_cspecs_parser;
    }

    std::set<std::string> DiffTransEnumParser::expected() {
      std::set<std::string> res = EnumInputParser::expected();
      res.insert({"cspecs", "require", "exclude"});
      return res;
    }


    /// \brief Construct with an IntegralCluster
    DiffusionTransformationEnum::DiffusionTransformationEnum(const IntegralCluster &clust) :
      m_cluster(clust),
      m_current(new DiffusionTransformation(clust.prim())) {

      // initialize to/from counter
      _init_occ_counter();

      // initialize the specie trajectory for the first valid set of to/from occupation values
      m_from_loc = _init_from_loc(m_occ_counter());
      m_to_loc = _init_to_loc(m_occ_counter());

      // set initial DiffTrans
      _set_current();

      if(!m_current->is_valid()) {
        increment();
      }

      if(!m_occ_counter.valid()) {
        _invalidate();
      }
      else {
        this->_initialize(&(*m_current));
        _set_step(0);
      }
    }

    const std::string DiffusionTransformationEnum::enumerator_name = "DiffusionTransformationEnum";
    const std::string DiffusionTransformationEnum::interface_help =
      std::string("DiffusionTransformationEnum: \n\n")
      + PrimPeriodicClustersByMaxLength::cspecs_help
      + SpeciesSetParser::require_all_help
      + SpeciesSetParser::exclude_all_help
      + EnumInputParser::standard_help +

      "  Example:\n"
      "  {\n"
      "   \"require\":[\"Va\"],\n"
      "   \"exclude\":[],\n"
      "   \"cspecs\":{\n"
      "      \"orbit_branch_specs\" : { \n"
      "       \"2\" : {\"max_length\" : 5.01},\n"
      "       \"3\" : {\"max_length\" : 5.01}\n"
      "      }\n"
      "    }\n"
      "  }\n\n"
      ;

    //    template<typename T, typename...Args>
    //    T get_else(const jsonParser &json, const std::string &key, const T &default_value, Args &&... args) {
    //      auto it = json.find(key);
    //      if(it != json.end()) {
    //        return it->get<T>(std::forward<Args>(args)...);
    //      }
    //      else {
    //        return default_value;
    //      }
    //    }

    /// Implements increment
    void DiffusionTransformationEnum::increment() {

      // get the next valid DiffTrans
      // to do so, get the next valid species trajectory
      do {

        // by taking a permutation of possible 'to' species position
        bool valid_perm = std::next_permutation(m_to_loc.begin(), m_to_loc.end());

        // if no more possible species trajectory,
        if(!valid_perm) {

          // get next valid from/to occupation values
          do {
            m_occ_counter++;
            _update_current_occ_transform();
          }
          while(m_occ_counter.valid() && !m_current->is_valid_occ_transform());

          // if no more possible from/to occupation values, return
          if(!m_occ_counter.valid()) {
            _invalidate();
            return;
          }
          else {
            m_from_loc = _init_from_loc(m_occ_counter());
            m_to_loc = _init_to_loc(m_occ_counter());
            _update_current_occ_transform();
            _set_current_loc();
          }
        }
        _update_current_to_loc();
      }
      while(!m_current->is_valid_species_traj());

      _increment_step();
    }

    /// Implements run
    template<typename DatabaseType>
    int DiffusionTransformationEnum::run(
      DiffTransEnumParser &parser,
      DatabaseType &db) {

      const PrimClex &primclex = parser.primclex();
      Log &log = primclex.log();
      COORD_TYPE coord_mode = parser.coord_mode();
      ORBIT_PRINT_MODE orbit_print_mode = parser.orbit_print_mode();
      std::string dry_run_msg = parser.dry_run_msg();
      std::set<std::string> require = parser.required_species();
      std::set<std::string> exclude = parser.excluded_species();
      std::string lead = dry_run_msg;

      log.begin(enumerator_name);
      Index Ninit = db.size();
      log << lead << "# diffusion transformations in this project: " << Ninit << "\n" << std::endl;

      log.increase_indent();
      log << std::endl;

      if(!parser.valid()) {
        parser.print_errors(log);
        log << std::endl << parser.report() << std::endl;
        log.decrease_indent();
        return 1;
      }
      if(parser.all_warnings().size()) {
        parser.print_warnings(log);
        log << std::endl << parser.report() << std::endl;
      }

      Printer<Kinetics::DiffusionTransformation> dt_printer(6, '\n', coord_mode);

      log.begin<Log::verbose>("Calculate cluster orbits");
      std::vector<PrimPeriodicIntegralClusterOrbit> orbits;
      auto end = make_prim_periodic_orbits(
                   primclex.prim(),
                   parser.cspecs().self, // TODO
                   alloy_sites_filter,
                   primclex.crystallography_tol(),
                   std::back_inserter(orbits),
                   primclex.log());
      print_clust(orbits.begin(), orbits.end(), log, orbit_print_mode, coord_mode);
      log << std::endl;

      log.begin<Log::verbose>("Calculate diff_trans orbits");
      std::vector< PrimPeriodicDiffTransOrbit > diff_trans_orbits;
      auto end2 = make_prim_periodic_diff_trans_orbits(
                    orbits.begin(),
                    orbits.end(),
                    primclex.crystallography_tol(),
                    std::back_inserter(diff_trans_orbits),
                    &primclex);
      print_clust(diff_trans_orbits.begin(), diff_trans_orbits.end(), log, orbit_print_mode, coord_mode);
      log << std::endl;

      log.begin<Log::verbose>("Check diff_trans orbits");
      log << lead << "COORD_MODE = " << coord_mode << std::endl << std::endl;
      lead = log.indent_str() + dry_run_msg;
      for(auto &diff_trans_orbit : diff_trans_orbits) {
        log << lead << "Checking: \n";
        log.increase_indent();
        dt_printer.print(diff_trans_orbit.prototype(), log);
        auto species_count = diff_trans_orbit.prototype().species_count();

        if(!includes_all(species_count, require.begin(), require.end())) {
          log << lead << "- Missing required species, do not insert" << std::endl;
          log.decrease_indent();
          log << std::endl;
          continue;
        }
        if(!excludes_all(species_count, exclude.begin(), exclude.end())) {
          log << lead << "- Includes excluded species, do not insert" << std::endl;
          log.decrease_indent();
          log << std::endl;
          continue;
        }

        //insert current into database
        auto res = db.insert(diff_trans_orbit);
        if(res.second) {
          log << lead << "- Inserted as: " << res.first->name() << std::endl;
        }
        else {
          log << lead << "- Already exists: " << res.first->name() << std::endl;
        }
        log << std::endl;
        log.decrease_indent();

      }

      log << lead << "  DONE." << std::endl << std::endl;
      log.decrease_indent();
      lead = log.indent_str() + dry_run_msg;

      Index Nfinal = db.size();

      log << lead << "# new diffusion transformations: " << Nfinal - Ninit << "\n";
      log << lead << "# diffusion transformations in this project: " << Nfinal << "\n" << std::endl;

      return 0;
    }

    /// Implements run using any set-like container
    template<typename DatabaseType>
    int DiffusionTransformationEnum::run(
      const PrimClex &primclex,
      const jsonParser &_kwargs,
      const Completer::EnumOption &enum_opt,
      DatabaseType &db) {

      jsonParser input {_kwargs};
      DiffTransEnumParser parser(primclex, input, enum_opt, fs::path(), true);
      return DiffusionTransformationEnum::run(parser, db);
    }

    /// Implements run using default database  (and commits)
    int DiffusionTransformationEnum::run(
      const PrimClex &primclex,
      const jsonParser &_kwargs,
      const Completer::EnumOption &enum_opt) {

      jsonParser input {_kwargs};
      DiffTransEnumParser parser(primclex, input, enum_opt, fs::path(), true);
      auto &db = primclex.db<PrimPeriodicDiffTransOrbit>();
      int res = DiffusionTransformationEnum::run(parser, db);
      if(res) {
        return res;
      }

      if(!parser.dry_run()) {
        primclex.log() << "Writing diffusion transformation database..." << std::endl;
        db.commit();
        primclex.log() << "  DONE" << std::endl;
      }

      return 0;
    }

    template int DiffusionTransformationEnum::run<std::set<PrimPeriodicDiffTransOrbit>>(
      const PrimClex &,
      const jsonParser &,
      const Completer::EnumOption &,
      std::set<PrimPeriodicDiffTransOrbit> &);

    // -- Unique -------------------

    const Structure &DiffusionTransformationEnum::prim() const {
      return m_cluster.prim();
    }

    const IntegralCluster &DiffusionTransformationEnum::cluster() const {
      return m_cluster;
    }

    /// \brief The occ_counter contains the from/to occupation values for each site
    ///
    /// - Layout is: [from values | to values ]
    ///
    void DiffusionTransformationEnum::_init_occ_counter() {
      Index N = cluster().size();
      std::vector<Index> init_occ(N * 2, 0);
      std::vector<Index> final_occ(N * 2);
      for(int i = 0; i < N; i++) {
        final_occ[i] = cluster()[i].site().site_occupant().size() - 1;
        final_occ[i + N] = final_occ[i];
      }
      std::vector<Index> incr(N * 2, 1);

      m_occ_counter = Counter<std::vector<Index> >(init_occ, final_occ, incr);
    }

    /// \brief Returns container of 'from' species locations
    std::vector<SpeciesLocation> DiffusionTransformationEnum::_init_from_loc(const std::vector<Index> &occ_values) {
      return _init_loc(occ_values, 0);
    }

    /// \brief Returns container of 'to' species locations
    std::vector<SpeciesLocation> DiffusionTransformationEnum::_init_to_loc(const std::vector<Index> &occ_values) {
      return _init_loc(occ_values, cluster().size());
    }

    /// \brief Returns container of 'from' or 'to' species locations
    ///
    /// - offset == 0 for 'from', N for 'to' species locations
    ///
    std::vector<SpeciesLocation> DiffusionTransformationEnum::_init_loc(const std::vector<Index> &occ_values, Index offset) {

      Index N = cluster().size();
      std::vector<SpeciesLocation> loc;
      // for each 'from' occupant
      for(Index i = 0; i < N; ++i) {
        Index occ = occ_values[i + offset];
        UnitCellCoord uccoord = cluster()[i];
        Index mol_size = uccoord.site().site_occupant()[occ].size();
        // for each species
        for(Index j = 0; j < mol_size; ++j) {
          loc.emplace_back(uccoord, occ, j);
        }
      }
      return loc;
    }

    /// \brief Uses m_cluster, m_occ_counter, m_from_loc, and m_to_loc to set m_current
    void DiffusionTransformationEnum::_set_current() {
      m_current->occ_transform().clear();
      m_current->occ_transform().reserve(cluster().size());
      for(const auto &uccoord : cluster()) {
        m_current->occ_transform().emplace_back(uccoord, 0, 0);
      }
      _update_current_occ_transform();
      _set_current_loc();
      _update_current_to_loc();
    }

    void DiffusionTransformationEnum::_update_current_occ_transform() {
      auto N = cluster().size();
      Index i = 0;
      for(auto &t : m_current->occ_transform()) {
        t.from_value = m_occ_counter()[i];
        t.to_value = m_occ_counter()[i + N];
        ++i;
      }
    }

    void DiffusionTransformationEnum::_set_current_loc() {
      m_current->species_traj().clear();
      m_current->species_traj().reserve(m_from_loc.size());
      for(const auto &t : m_from_loc) {
        m_current->species_traj().emplace_back(t, t);
      }
    }

    void DiffusionTransformationEnum::_update_current_to_loc() {
      auto it = m_to_loc.begin();
      for(auto &t : m_current->species_traj()) {
        t.to = *it++;
      }
    }

#define  PRIM_PERIODIC_DIFF_TRANS_ORBITS_INST(INSERTER, CLUSTER_IT) \
    \
    template INSERTER make_prim_periodic_diff_trans_orbits<INSERTER, CLUSTER_IT>( \
      CLUSTER_IT begin, \
      CLUSTER_IT end, \
      double xtal_tol, \
      INSERTER result, \
      const PrimClex*); \
    \

#define _VECTOR_IT(ORBIT) std::vector<ORBIT>::iterator
#define _VECTOR_INSERTER(ORBIT) std::back_insert_iterator<std::vector<ORBIT> >

#define  PRIM_PERIODIC_DIFF_TRANS_ORBITS_VECTOR_INST(DIFF_TRANS_ORBIT, CLUSTER_ORBIT) \
      PRIM_PERIODIC_DIFF_TRANS_ORBITS_INST( \
        _VECTOR_INSERTER(DIFF_TRANS_ORBIT), \
        _VECTOR_IT(CLUSTER_ORBIT))

    PRIM_PERIODIC_DIFF_TRANS_ORBITS_VECTOR_INST(
      PrimPeriodicDiffTransOrbit,
      PrimPeriodicIntegralClusterOrbit)

  }
}
