#include "casm/kinetics/DiffTransConfigEnumOccPerturbations.hh"

#include "casm/kinetics/DiffTransConfiguration_impl.hh"
#include "casm/symmetry/ConfigSubOrbits_impl.hh"
#include "casm/clusterography/ClusterOrbits_impl.hh"
#include "casm/clex/Supercell.hh"
#include "casm/clex/FilteredConfigIterator.hh"

#include "casm/app/ProjectSettings.hh"
#include "casm/app/AppIO_impl.hh"
#include "casm/app/QueryHandler_impl.hh"
#include "casm/database/Selection_impl.hh"
#include "casm/database/ScelDatabase.hh"
#include "casm/database/ConfigDatabase.hh"
#include "casm/database/DiffTransConfigDatabase.hh"
#include "casm/database/DiffTransOrbitDatabase.hh"


extern "C" {
  CASM::EnumInterfaceBase *make_DiffTransConfigEnumOccPerturbations_interface() {
    return new CASM::EnumInterface<CASM::Kinetics::DiffTransConfigEnumOccPerturbations>();
  }
}


namespace CASM {

  namespace Kinetics {

    DiffTransConfigEnumOccPerturbations::DiffTransConfigEnumOccPerturbations(
      const Configuration &background_config,
      const PrimPeriodicDiffTransOrbit &diff_trans_orbit,
      const jsonParser &local_cspecs) :
      m_background_config(background_config),
      m_diff_trans_orbit(diff_trans_orbit),
      m_local_cspecs(local_cspecs),
      m_scel_sym_compare(_supercell()),
      m_include_unperturbed(true),
      m_skip_subclusters(true),
      m_curr(OccPerturbation(_prim())) {

      this->_initialize();

      // initialize data
      _init_base();
      _init_local_orbits();
      _init_perturbations_data();

      // check if initial perturb is valid or not
      // (currently null perturbation is valid, so I think this should always be valid)
      _update_current_perturb();
      if(!check_increment()) {
        increment();
      }
      else {
        _set_current(m_curr.perturb);
      }

      if(valid()) {
        this->_set_step(0);
        m_current->set_source(this->source(step()));
      }
    }

    const std::string DiffTransConfigEnumOccPerturbations::enumerator_name = "DiffTransConfigEnumOccPerturbations";

    std::string DiffTransConfigEnumOccPerturbations::interface_help() {
      return "DiffTransConfigEnumOccPerturbations: \n\n"

             "  orbit_names: JSON array of strings \n"
             "    The names of diffusion transformation orbits whose local environments\n"
             "    will be perturbed.\n\n"

             "  orbit_selection: string \n"
             "    The names of a selection diffusion transformation orbits whose local \n"
             "    environments will be perturbed.\n\n"

             "  background_scel: string (optional, default=None)\n"
             "    The name of the supercell in which DiffTransConfiguration should be \n"
             "    created. If given, the supercell must be a supercell that can contain \n"
             "    all of the configurations specified by \"background_confignames\" and \n"
             "    \"background_selection\". Prim point group operations are allowed to \n"
             "    transform configurations so that they fit into the supercell. Newly \n"
             "    created configurations in the background supercell will be added to \n"
             "    the configuration database.\n\n"

             "  background_confignames: JSON array of string \n "
             "    Names of background configurations to be perturbed. Should be names of \n"
             "    Configurations that exist in your CASM project.\n\n"

             "  background_selection: string\n "
             "    The name of a selection of background configurations to be perturbed.\n\n"

             "  local_cspecs: JSON object (optional,default= {}) \n "
             "    Specify the clusters that should be use for perturbations. The string \n"
             "    \"local_cspecs\" should be a local cspecs style initialization as used \n "
             "    in 'casm bset' local basis function enumeration. This option takes \n"
             "    precedence over the following option.\n\n"

             "  local_cspecs_filepath: string (optional,default=\"\") \n "
             "    Indicate the local cspecs file that specifies the clusters around \n"
             "    the transformation that should be perturbed. The string \n"
             "    \"local_cspecs_filepath\" should be the file path to a JSON file \n"
             "    containing the local cspecs.\n\n"

             "  Example:\n"
             "  {\n"
             "    \"orbit_names\": [\n"
             "      \"diff_trans/0\",\n"
             "      \"diff_trans/1\"\n"
             "    ],\n"
             "    \"orbit_selection\": \"low_barrier_diff_trans\",\n"
             "    \"background_configs\": [\n"
             "      \"SCEL8_2_2_2_0_0_0/2\",\n"
             "      \"SCEL8_2_2_2_0_0_0/13\"\n"
             "     ],\n"
             "    \"background_selection\": \"groundstates\",\n"
             "    \"local_cspecs\": {\n"
             "      \"orbit_branch_specs\" : { \n"
             "        \"1\" : {\"cutoff_radius\" : 6.0},\n"
             "        \"2\" : {\"max_length\" : 6.01,\"cutoff_radius\" : 6.0},\n"
             "        \"3\" : {\"max_length\" : 4.01,\"cutoff_radius\" : 5.0}\n"
             "      }\n"
             "    }\n"
             "  }\n\n";
    }


    namespace {

      /// --- Implementation for 'run' ---

      /// Take input selection of background configurations, make each fill the
      /// background supercell, and return updated selection with just configurations
      /// in the background supercell selected.
      ///
      /// - Will saved newly created configurations in database
      /// - Throws if any will not fill supercell
      DB::Selection<Configuration> _fill_background_scel(
        const DB::Selection<Configuration> &background_sel,
        const Supercell &background_scel) {

        const PrimClex &primclex = background_sel.primclex();
        Log &log = primclex.log();
        std::string scelname = background_scel.name();
        std::string enum_name = DiffTransConfigEnumOccPerturbations::enumerator_name;

        log << "Checking if configurations will tile supercell " << scelname << std::endl;

        // Construct a selection to store background configurations in
        DB::Selection<Configuration> tmp_sel(primclex, "EMPTY");

        // Check if all requested configurations will fill the background scel
        auto begin = primclex.prim().factor_group().begin();
        auto end = primclex.prim().factor_group().end();
        auto tol = primclex.crystallography_tol();
        for(const auto &config : background_sel.selected()) {

          // Check if configuration can tile background scel
          auto res = is_supercell(background_scel.lattice(), config.ideal_lattice(), begin, end, tol);
          if(res.first == end) {
            std::string msg = "Error in " + enum_name + ": Configuration "
                              + config.name() + " cannot tile supercell " + scelname;
            throw std::runtime_error(msg);
          }

          // Insert configuration in background scel
          Configuration background_config = config.fill_supercell(background_scel, *res.first);
          auto insert_res = background_config.insert();

          log << config.name() << " fills the supercell as " << insert_res.canonical_it.name() << std::endl;

          // Store in selection
          tmp_sel.data()[insert_res.canonical_it.name()] = true;
        }
        log << "Writing supercell database..." << std::endl;
        primclex.db<Supercell>().commit();
        log << "  DONE" << std::endl;

        log << "Writing configuration database..." << std::endl;
        primclex.db<Configuration>().commit();
        log << "  DONE" << std::endl;

        return tmp_sel;
      }

      /// Parse input to get local_cspecs JSON
      jsonParser _parse_local_cspecs(
        const jsonParser &kwargs) {

        if(kwargs.contains("local_cspecs")) {
          // nothing necessary
          return kwargs["local_cspecs"].get<jsonParser>();
        }
        else if(kwargs.contains("local_cspecs_filepath")) {
          //look for filepath given and load it
          fs::path local_cspecs_filepath = kwargs["local_cspecs_filepath"].get<fs::path>();
          if(!fs::exists(local_cspecs_filepath)) {
            std::string msg = "Error in "
                              + DiffTransConfigEnumOccPerturbations::enumerator_name + ": "
                              "\"local_cspecs_filepath\" file: " + local_cspecs_filepath.string()
                              + "does not exist.";
            throw std::runtime_error(msg);
          }
          return jsonParser{local_cspecs_filepath};
        }

        std::string msg = "Error in " + DiffTransConfigEnumOccPerturbations::enumerator_name + ": One of "
                          "\"local_cspecs\" or \"local_cspecs_filepath\" must be given.";
        throw std::runtime_error(msg);
      }

      /// Check if local clusters overlap
      void _check_overlap(
        const PrimClex &primclex,
        const Configuration &bg_config,
        const PrimPeriodicDiffTransOrbit &dtorbit,
        const jsonParser &local_cspecs) {

        /// check if supercell is big enough for local cspecs here
        /// give warning if not

        std::vector<ScelPeriodicOrbit<IntegralCluster>> local_orbits;
        std::vector<PermuteIterator> diff_trans_g {
          dtorbit.prototype().invariant_subgroup(bg_config.supercell())};
        SymGroup diff_trans_sym_g { make_sym_group(diff_trans_g) };
        ScelPeriodicSymCompare<IntegralCluster> scel_sym_compare {bg_config.supercell()};

        make_local_orbits(
          dtorbit.prototype(),
          diff_trans_sym_g,
          scel_sym_compare,
          local_cspecs,
          alloy_sites_filter,
          primclex.crystallography_tol(),
          std::back_inserter(local_orbits),
          null_log());

        if(has_local_neighborhood_overlap(local_orbits, bg_config.supercell())) {
          std::string msg = "Warning in " +
                            DiffTransConfigEnumOccPerturbations::enumerator_name + ": Choice of background "
                            "configuration " + bg_config.name() + " results in an overlap in the local "
                            "clusters of " + dtorbit.name() + " with their periodic images. Consider "
                            "choosing a larger background configuration or smaller set of local clusters.";
          primclex.log() << msg << std::endl;
        }
      }

      /// Construct and execute enumerator
      template<typename DatabaseType>
      void _enumerate(
        const PrimClex &primclex,
        const Configuration &bg_config,
        const PrimPeriodicDiffTransOrbit &dtorbit,
        const jsonParser &local_cspecs,
        std::vector<std::string> filter_expr,
        DatabaseType &db,
        bool dry_run) {

        std::string dry_run_msg = CASM::dry_run_msg(dry_run);
        primclex.log() << dry_run_msg << "\tUsing " << dtorbit.name() << "... " << std::flush;
        Index Ninit_spec = db.size();

        DiffTransConfigEnumOccPerturbations enumerator(bg_config, dtorbit, local_cspecs);
        auto begin = enumerator.begin();
        auto end = enumerator.end();
        if(!filter_expr.empty()) {
          try {
            auto fbegin = filter_begin(
                            begin,
                            end,
                            filter_expr,
                            primclex.settings().query_handler<DiffTransConfiguration>().dict());
            auto fend = filter_end(enumerator.end());
            db.insert(fbegin, fend);
          }
          catch(std::exception &e) {
            std::string msg = "Cannot filter " + traits<DiffTransConfiguration>::name
                              + " using the expression provided: " + e.what();
            throw std::runtime_error(msg);
          }
        }
        else {
          db.insert(begin, end);
        }

        Index Nfinal_spec = db.size();
        primclex.log() << dry_run_msg << "Found " << Nfinal_spec - Ninit_spec
                       << " new " << traits<DiffTransConfiguration>::short_name << std::endl;
      }

    }

    /// Enumerate DiffTransConfigEnumOccPerturbations into any std::set-like database
    template<typename DatabaseType>
    int DiffTransConfigEnumOccPerturbations::run(
      const PrimClex &primclex,
      const jsonParser &kwargs,
      const Completer::EnumOption &enum_opt,
      DatabaseType &db) {

      Log &log = primclex.log();

      bool dry_run = CASM::dry_run(kwargs, enum_opt);
      std::string dry_run_msg = CASM::dry_run_msg(dry_run);

      // Validate and construct input
      DB::Selection<PrimPeriodicDiffTransOrbit> dtorbit_sel = make_selection<PrimPeriodicDiffTransOrbit>(
                                                                primclex, kwargs, "orbit_names", "orbit_selection", enumerator_name, OnError::THROW);

      DB::Selection<Configuration> background_sel = make_selection<Configuration>(
                                                      primclex, kwargs, "background_confignames", "background_selection", enumerator_name, OnError::THROW);

      // If requested, fill all configurations into background supercell
      if(kwargs.contains("background_scel")) {
        // Get requested background scel
        std::string scelname = kwargs["background_scel"].get<std::string>();
        const auto &background_scel = make_supercell(primclex, scelname);
        background_sel = _fill_background_scel(background_sel, background_scel);
      }

      jsonParser local_cspecs = _parse_local_cspecs(kwargs);

      std::vector<std::string> filter_expr = make_enumerator_filter_expr(kwargs, enum_opt);
      std::string type_name = traits<DiffTransConfiguration>::name;

      Index Ninit = db.size();
      log << dry_run_msg << "# " << type_name << " in this project: " << Ninit << "\n" << std::endl;

      log.begin(enumerator_name);
      for(auto bg_config : background_sel.selected()) {
        auto prim_config = bg_config.primitive().in_canonical_supercell();
        log << dry_run_msg << "Searching in " << bg_config.name()
            << " (primitive = " << prim_config.name() << ") ..." << std::endl;

        for(const auto &dtorbit : dtorbit_sel.selected()) {
          _check_overlap(primclex, bg_config, dtorbit, local_cspecs);
          _enumerate(primclex, bg_config, dtorbit, local_cspecs, filter_expr, db, dry_run);
        }
      }

      log << dry_run_msg << "  DONE." << std::endl << std::endl;

      Index Nfinal = db.size();

      log << dry_run_msg << "# new " << type_name << ": " << Nfinal - Ninit << "\n";
      log << dry_run_msg << "# " << type_name << " in this project: " << Nfinal << "\n" << std::endl;

      return 0;
    }

    /// Enumerate DiffTransConfigEnumOccPerturbations into the project database
    int DiffTransConfigEnumOccPerturbations::run(
      const PrimClex &primclex,
      const jsonParser &kwargs,
      const Completer::EnumOption &enum_opt) {

      auto &db = primclex.db<DiffTransConfiguration>();
      std::string type_name = traits<DiffTransConfiguration>::name;

      return DiffTransConfigEnumOccPerturbations::run(
               primclex,
               kwargs,
               enum_opt,
               primclex.db<DiffTransConfiguration>());

      if(!CASM::dry_run(kwargs, enum_opt)) {
        primclex.log() << "Writing " << type_name << " database..." << std::endl;
        db.commit();
        primclex.log() << "  DONE" << std::endl;
      }
      return 0;
    }

    template int DiffTransConfigEnumOccPerturbations::run<std::set<DiffTransConfiguration>>(
      const PrimClex &primclex,
      const jsonParser &kwargs,
      const Completer::EnumOption &,
      std::set<DiffTransConfiguration> &);

    ///------------------------------------Internal functions-------------------------------///

    void DiffTransConfigEnumOccPerturbations::partial_increment(bool complete_perturb) {
      // increment occ_counter once
      // if no more occupations, increment m_local_orbit_it
      // if no more local orbits, increment m_base_it
      // if no more base diff trans, _invalidate

      ++m_occ_counter;
      ++m_occ_counter_index;
      if(m_occ_counter.valid()) {
        _update_current_perturb(complete_perturb);
      }
      else if(++m_local_orbit_it != m_local_orbit.end()) {
        _init_perturbations_data();
        _update_current_perturb(complete_perturb);
      }
      else if(++m_base_it != m_base.end()) {
        _init_local_orbits();
        _init_perturbations_data();
        _update_current_perturb(complete_perturb);
      }
      else {
        // at this point m_base_it == m_base.end(): pass
      }
      return;
    }

    /// \returns True if current perturb is valid and canonical or enumerator is complete
    bool DiffTransConfigEnumOccPerturbations::check_increment() {

      bool valid = m_curr.is_not_subcluster && m_curr.is_canonical;
      if(!valid && m_base_it != m_base.end()) {
        return false;
      }

      if(valid) {
        // find canonical from_config set current
        this->_increment_step();
        _set_current(m_curr.perturb);
      }
      else {
        // if no more base diff trans, we're done
        _invalidate();
      }
      return true;
    }

    void DiffTransConfigEnumOccPerturbations::increment() {

      // increment occ_counter until a canonical perturb is found
      // if no more occupations, increment m_local_orbit_it
      // if no more local orbits, increment m_base_it
      // if no more base diff trans, _invalidate

      // current Perturbation & 'is_valid'

      do {
        partial_increment();
      }
      while(!check_increment());
    };

    double DiffTransConfigEnumOccPerturbations::_tol() const {
      return m_background_config.primclex().crystallography_tol();
    }

    /// Supercell
    const Structure &DiffTransConfigEnumOccPerturbations::_prim() const {
      return m_background_config.prim();
    }

    /// Base DiffTransConfig supercell
    const Supercell &DiffTransConfigEnumOccPerturbations::_supercell() const {
      return m_background_config.supercell();
    }

    void DiffTransConfigEnumOccPerturbations::_init_base() {

      // Make suborbit generating DiffusionTransformations
      std::vector<DiffusionTransformation> subprototypes;
      make_suborbit_generators(m_diff_trans_orbit, m_background_config, std::back_inserter(subprototypes));

      // Put them in canonical form, and store similarly transformed
      //   background config. Base constructor will also generate diff trans
      //   invariant group.
      ScelCanonicalGenerator<DiffusionTransformation> gen(_supercell());
      for(const auto &diff_trans : subprototypes) {
        auto canonical_diff_trans = gen(diff_trans);
        m_base.emplace_back(
          canonical_diff_trans.sorted(),
          copy_apply(gen.to_canonical(), m_background_config));
        if(!canonical_diff_trans.sorted().is_canonical(_supercell())) {
          throw std::runtime_error("ScelCanonicalGenerator did not work to make diff_trans canonical");
        }
      }

      m_base_it = m_base.begin();
    }

    /// Constructor for data structure holding base diff trans configuration in
    /// their canonical form
    ///
    /// \param _diff_trans DiffTrans in canonical form for the supercell
    /// \param _config Background config transformed appropriately for _diff_trans
    DiffTransConfigEnumOccPerturbations::Base::Base(
      const DiffusionTransformation &_diff_trans,
      const Configuration &_config) :
      diff_trans(_diff_trans),
      config(make_attachable(diff_trans, _config)),
      diff_trans_g(diff_trans.invariant_subgroup(config.supercell())),
      generating_g(config.invariant_subgroup(diff_trans_g.begin(), diff_trans_g.end())),
      generating_sym_g(make_sym_group(generating_g)) {
    }

    const std::vector<DiffTransConfigEnumOccPerturbations::Base> &
    DiffTransConfigEnumOccPerturbations::base() const {
      return m_base;
    }

    Index DiffTransConfigEnumOccPerturbations::base_index() const {
      return std::distance(
               m_base.begin(),
               std::vector<Base>::const_iterator(m_base_it));
    }

    const std::vector<ScelPeriodicOrbit<IntegralCluster>> &
    DiffTransConfigEnumOccPerturbations::local_orbit() const {
      return m_local_orbit;
    }

    Index DiffTransConfigEnumOccPerturbations::local_orbit_index() const {
      return std::distance(
               m_local_orbit.begin(),
               std::vector<ScelPeriodicOrbit<IntegralCluster>>::const_iterator(m_local_orbit_it));
    }

    Index DiffTransConfigEnumOccPerturbations::occ_counter_index() const {
      return m_occ_counter_index;
    }

    /// Return the current OccPerturbation (non-canonical) and whether it is valid
    ///
    /// - If not valid, the OccPerturbation itself may not be complete
    const DiffTransConfigEnumOccPerturbations::CurrentPerturbation &
    DiffTransConfigEnumOccPerturbations::current_perturb() const {
      return m_curr;
    }

    /// Generate local orbits for current base diff trans
    void DiffTransConfigEnumOccPerturbations::_init_local_orbits() {

      /// Make all local orbits
      ///
      /// - local orbit generating group is set of operations that leave
      ///   diff_trans and config invariant
      std::vector<ScelPeriodicOrbit<IntegralCluster>> _tmp;
      make_local_orbits(
        m_base_it->diff_trans,
        m_base_it->generating_sym_g,
        m_scel_sym_compare,
        m_local_cspecs,
        alloy_sites_filter,
        _tol(),
        std::back_inserter(_tmp),
        null_log());

      /// Exclude orbits that would alter the hop due small supercell size

      // get list of linear indices of hopping sites
      std::set<Index> hop_cluster_indices;
      for(auto &traj : m_base_it->diff_trans.species_traj()) {
        auto res = hop_cluster_indices.insert(_supercell().linear_index(traj.from.uccoord));
        if(!res.second) {
          std::string msg = "Error in " + enumerator_name + ": Background "
                            "configuration is too small to contain the requested diffusion transformation.\n";
          auto &err_log = _supercell().primclex().err_log();
          err_log << "Background config: " << m_background_config.name() << std::endl;
          err_log << "Diff trans: \n" << m_diff_trans_orbit.prototype() << std::endl;

          throw std::runtime_error(msg);
        }
      }

      // lambda function returns true if uccoord is not in hop cluster
      auto uccoord_does_not_overlap = [&](const UnitCellCoord & uccoord) {
        Index index = this->_supercell().linear_index(uccoord);
        return hop_cluster_indices.find(index) == hop_cluster_indices.end();
      };

      // lambda function returns true if no overlap between hop cluster and local cluster
      auto orbit_does_not_overlap = [&](const ScelPeriodicOrbit<IntegralCluster> &test) {
        const auto &proto = test.prototype();
        return std::all_of(proto.begin(), proto.end(), uccoord_does_not_overlap);
      };
      m_local_orbit.clear();
      std::copy_if(_tmp.begin(), _tmp.end(), std::back_inserter(m_local_orbit), orbit_does_not_overlap);
      m_local_orbit_it = m_local_orbit.begin();


      /* print_clust(
         m_local_orbit.begin(),
         m_local_orbit.end(),
         std::cout,
         PrototypePrinter<IntegralCluster>());
       */
    }

    /// Generate the 'from_value' for the perturbation,
    ///   the 'to_value' counter for the perturbation,
    ///   and the local orbit prototype invariant subgroup
    ///   (w/ respect to base diff trans invariant group)
    void DiffTransConfigEnumOccPerturbations::_init_perturbations_data() {

      const Configuration &from_config = m_base_it->config;
      const IntegralCluster &proto = m_local_orbit_it->prototype();

      /// Set 'm_from_value'
      m_from_value = Eigen::VectorXi::Zero(proto.size());
      for(int i = 0; i < proto.size(); ++i) {
        m_from_value(i) = from_config.occ(_supercell().linear_index(proto[i]));
      }

      /// Construct counter
      Eigen::VectorXi max_count(proto.size());
      for(int i = 0; i < proto.size(); ++i) {
        max_count(i) = proto.prim().basis()[proto[i].sublat()].site_occupant().size() - 1;
      }
      m_occ_counter = EigenCounter<Eigen::VectorXi>(
                        Eigen::VectorXi::Zero(proto.size()),
                        max_count,
                        Eigen::VectorXi::Constant(proto.size(), 1));
      m_occ_counter_index = 0;
    }

    void DiffTransConfigEnumOccPerturbations::_update_current_perturb(bool complete_perturb) {

      m_curr.perturb.elements().clear();
      const auto &proto = m_scel_sym_compare.prepare(m_local_orbit_it->prototype());

      // check for null perturbation
      if(proto.size() == 0 && m_include_unperturbed) {
        m_curr.is_not_subcluster = true;
        m_curr.is_canonical = true;
      }

      // generate next perturbation consistent with from_config
      m_curr.is_not_subcluster = true;
      for(int i = 0; i < proto.size(); ++i) {
        // check for subcluster perturbation
        if(m_from_value(i) == m_occ_counter()[i] && m_skip_subclusters) {
          m_curr.is_not_subcluster = false;
          m_curr.is_canonical = false;
          if(!complete_perturb) {
            return;
          }
        }
        m_curr.perturb.elements().emplace_back(proto[i], m_from_value(i), m_occ_counter()[i]);
      }

      /*auto begin = m_base_it->generating_g.begin();
      auto end = m_base_it->generating_g.end();
      // return perturbation and check if canonical wrt local orbit sub group
      if(m_curr.is_not_subcluster) {
        m_curr.is_canonical = m_curr.perturb.is_canonical(_supercell(), begin, end);
      }
      }

      /// Applies current perturbation to m_base_config and stores result in m_current
      void DiffTransConfigEnumOccPerturbations::_set_current(const OccPerturbation &perturb) {

      // generate perturbed from_config
      Configuration perturbed_from_config {m_base_it->config};
      perturb.apply_to(perturbed_from_config);
      // here we use the group of operations that leaves the diff_trans invariant,
      // and not necessarily the from_config, to find the canonical perturbed from_config

      auto to_canonical = perturbed_from_config.to_canonical(
                            m_base_it->diff_trans_g.begin(),
                            m_base_it->diff_trans_g.end());
      /*PermuteIterator to_canonical = *(m_base_it->diff_trans_g.begin());
      Configuration greatest = perturbed_from_config;
      for(auto it = m_base_it->diff_trans_g.begin(); it != m_base_it->diff_trans_g.end(); ++it) {
        DiffTransConfiguration tmp(make_attachable(m_base_it->diff_trans, copy_apply(*it, perturbed_from_config)), m_base_it->diff_trans);
        if(copy_apply(*it, perturbed_from_config) > greatest && tmp.has_valid_from_occ()) {
          greatest = copy_apply(*it, perturbed_from_config);
          to_canonical = *it;
        }
      }*/

      /// construct canonical DiffTransConfiguration as m_current
      m_current = notstd::make_cloneable<DiffTransConfiguration>(
                    make_attachable(m_base_it->diff_trans, copy_apply(to_canonical, perturbed_from_config)),
                    m_base_it->diff_trans);
      if(!m_current->has_valid_from_occ()) {
        throw std::runtime_error("Make attachable didn't seem to work in DTCEOP");
      }
      m_current->set_orbit_name("test");
      m_current->set_orbit_name(m_diff_trans_orbit.name());
      m_current->set_suborbit_ind(std::distance(m_base.begin(), m_base_it));
      m_current->set_bg_configname(m_background_config.primitive().name());
      m_current->set_source(this->source(step()));
      this->_set_current_ptr(&(*m_current));
      if(m_current->is_dud()) {
        throw std::runtime_error("Error in DTEOCP is dud");
      }


      // --- debug check ---
      // m_current should be in canonical form at this point
      if(!m_current->is_canonical()) {
        if(m_current->is_dud()) {
          throw std::runtime_error("Error in DTEOCP is dud and not canonical");
        }
        *m_current = m_current->canonical_form();
        if(m_current->is_dud()) {
          throw std::runtime_error("Error in DTEOCP is dud and canonical now");
        }
        // diff_trans is canonical wrt all supercell permutations
        const auto &dt = m_current->diff_trans();
        if(!dt.is_canonical(_supercell())) {
          throw std::runtime_error("Error in DiffTransConfigEnumOccPerturbations: diff trans not canonical");
        }

        // from_config is canonical wrt all supercell permutations that leave diff_trans invariant
        /*if(!m_current->from_config().is_canonical(m_base_it->diff_trans_g.begin(), m_base_it->diff_trans_g.end())) {
          throw std::runtime_error("Error in DiffTransConfigEnumOccPerturbations: from_config not canonical");
        }
        throw std::runtime_error("Error in DiffTransConfigEnumOccPerturbations: not canonical, for unknown reason");*/
      }

    }

    bool has_local_bubble_overlap(std::vector<LocalOrbit<IntegralCluster>> &local_orbits, const Supercell &scel) {
      std::set<int> present;
      std::set<UnitCellCoord> coords;
      for(auto &orbit : local_orbits) {
        for(auto &cluster : orbit) {
          for(int i = 0; i < cluster.size(); ++i) {
            coords.insert(cluster[i]);
          }
        }
      }
      for(auto &coord : coords) {
        if(!present.insert(scel.linear_index(coord)).second) {
          return true;
        }
      }
      //If no set insertion collision then no problems
      return false;
    }

    std::vector<Supercell> viable_supercells(std::vector<LocalOrbit<IntegralCluster>> &local_orbits, std::vector<Supercell> scel_options) {
      std::vector<Supercell> results;
      for(auto &scel : scel_options) {
        if(!has_local_bubble_overlap(local_orbits, scel)) {
          results.push_back(scel);
        }
      }
      return results;
    }

  }
}
