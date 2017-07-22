#include "casm/kinetics/DiffTransConfigEnumPerturbations.hh"

#include "casm/symmetry/ConfigSubOrbits_impl.hh"
#include "casm/symmetry/ScelOrbitGeneration_impl.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/Supercell.hh"
#include "casm/clex/Configuration.hh"
#include "casm/clex/ConfigCompare.hh"
#include "casm/clusterography/ClusterOrbits.hh"
#include "casm/clusterography/ClusterOrbits_impl.hh"
#include "casm/kinetics/DiffusionTransformation.hh"
#include "casm/kinetics/DiffTransEnumEquivalents.hh"

#include "casm/app/ProjectSettings.hh"
#include "casm/app/AppIO.hh"
#include "casm/app/AppIO_impl.hh"
#include "casm/app/QueryHandler.hh"
#include "casm/database/Selection.hh"
#include "casm/database/ConfigDatabase.hh"
#include "casm/database/DiffTransConfigDatabase.hh"
#include "casm/database/DiffTransOrbitDatabase.hh"


extern "C" {
  CASM::EnumInterfaceBase *make_DiffTransConfigEnumPerturbations_interface() {
    return new CASM::EnumInterface<CASM::Kinetics::DiffTransConfigEnumPerturbations>();
  }
}


namespace CASM {

  namespace Kinetics {

    Perturbation::Perturbation() {}

    Perturbation::Perturbation(std::set<OccupationTransformation> &from_set) {
      for(const OccupationTransformation &item : from_set) {
        this->insert(item);
      }
    }

    Perturbation &Perturbation::apply_sym(const SymOp &op) {
      std::set<OccupationTransformation> tmp;
      for(const OccupationTransformation &occ_trans : *this) {
        tmp.insert(copy_apply(op, occ_trans));
      }
      this->clear();
      for(const OccupationTransformation &item : tmp) {
        this->insert(item);
      }
      return *this;
    }

    Configuration &Perturbation::apply_to(Configuration &config) const {
      for(auto it = this->begin(); it != this->end(); ++it) {
        it->apply_to(config);
      }
      return config;
    }

    template<typename PermuteIteratorIt>
    bool Perturbation::is_canonical(PermuteIteratorIt begin, PermuteIteratorIt end) const {

      for(auto op = begin; op != end; ++op) {

        // loop over OccupationTransformation to perform lexicographical check
        for(auto it = this->begin(); it != this->end(); ++it) {
          auto tmp = copy_apply(op->sym_op(), *it);
          if(tmp > *it) {
            return false;
          }
        }
      }
      return true;
    }


    DiffTransConfigEnumPerturbations::DiffTransConfigEnumPerturbations(
      const Configuration &background_config,
      const PrimPeriodicDiffTransOrbit &diff_trans_orbit,
      const jsonParser &local_cspecs) :
      m_background_config(background_config),
      m_diff_trans_orbit(diff_trans_orbit),
      m_local_cspecs(local_cspecs),
      m_skip_subclusters(true) {

      std::cout << "begin constructor" << std::endl;
      this->_initialize();

      // initialize data
      _init_base();
      _init_local_orbits();
      _init_perturbations_data();

      std::cout << "check initial perturb" << std::endl;
      // check if initial perturb is valid or not
      auto res = _current_perturb();
      if(!res.second) {
        std::cout << "initial not valid" << std::endl;
        increment();
      }
      else {
        std::cout << "initial valid" << std::endl;
        _set_current(res.first);
      }

      std::cout << "initialize" << std::endl;
      if(valid()) {
        this->_set_step(0);
        m_current->set_source(this->source(step()));
      }

      std::cout << "end constructor" << std::endl;

    }

    const std::string DiffTransConfigEnumPerturbations::enumerator_name = "DiffTransConfigEnumPerturbations";
    const std::string DiffTransConfigEnumPerturbations::interface_help =
      "DiffTransConfigEnumPerturbations: \n\n"

      "  orbits: JSON array of strings \n"
      "    Indicate which diffusion transformation orbits are of interest. The \n"
      "    JSON array \"orbits\" should be the names of the orbits of interest.\n\n"

      "  background_configs: JSON array of strings \n "
      "    Indicate which configurations will be the background structures for the transformations to occur in.\n"
      "    The JSON array of strings \"background_configs\" should be names of Configurations\n"
      "    that exist in your CASM project.\n\n"

      "  local_cspecs: JSON object (optional,default= {}) \n "
      "    Specify the clusters around the transformation that should be \n"
      "    perturbed. The string \"local_cspecs\" should be a local cspecs style\n"
      "    initialization used in casm bset enumeration.\n "
      "    This option takes precedence over the following option.\n\n"

      "  local_cspecs_filepath: string (optional,default=\"\") \n "
      "    Indicate the local cspecs file that specifies the clusters around \n"
      "    the transformation that should be perturbed. The string \n"
      "    \"local_cspecs_filepath\" should be the file path to a JSON file \n"
      "    containing the local cspecs.\n\n"

      "  Example:\n"
      "  {\n"
      "   \"orbits\":[\"diff_trans/0\",\"diff_trans/1\"],\n"
      "   \"background_configs\":[\"SCEL8_2_2_2_0_0_0/2\"],\n"
      "    \"local_cspecs\":{\n"
      "       \"basis_functions\" : {\n"
      "        \"site_basis_functions\" : \"occupation\"\n"
      "      },\n"
      "      \"orbit_branch_specs\" : { \n"
      "       \"1\" : {\"cutoff_radius\" : 6.0},\n"
      "       \"2\" : {\"max_length\" : 6.01,\"cutoff_radius\" : 6.0},\n"
      "       \"3\" : {\"max_length\" : 4.01,\"cutoff_radius\" : 5.0}\n"
      "      }\n"
      "    }\n"
      "  }\n\n"
      ;

    void DiffTransConfigEnumPerturbations::increment() {

      std::cout << "begin increment" << std::endl;
      // increment occ_counter until a canonical perturb is found
      // if no more occupations, increment m_local_orbit_it
      // if no more local orbits, increment m_base_it
      // if no more base diff trans, _invalidate

      do {

        do {

          std::cout << "  here 0" << std::endl;
          // get next occupation on the current local orbit prototype
          while(++m_occ_counter) {

            std::cout << "  here 1" << std::endl;
            // std::pair<Perturbation, valid> = _current_perturb();
            auto res = _current_perturb();

            std::cout << "  here 2" << std::endl;
            if(res.second) {

              std::cout << "  here 3" << std::endl;
              // find canonical from_config set current
              this->_increment_step();
              std::cout << "  here 4" << std::endl;
              _set_current(res.first);
              std::cout << "  end increment (5)" << std::endl;
              return;
            }
          }

          std::cout << "  here 6" << std::endl;
          // if no more occupation, increment m_local_orbit_it
          if(++m_local_orbit_it != m_local_orbit.end()) {
            std::cout << "  here 7" << std::endl;
            _init_perturbations_data();
          }
          std::cout << "  here 8" << std::endl;

        }
        while(m_local_orbit_it != m_local_orbit.end());

        std::cout << "  here 9" << std::endl;
        // if no more local orbit, increment m_base_it
        if(++m_base_it != m_base.end()) {
          std::cout << "  here 10" << std::endl;
          _init_local_orbits();
        }
        std::cout << "  here 11" << std::endl;

      }
      while(m_base_it != m_base.end());

      std::cout << "  here 12" << std::endl;
      // if no more base diff trans, we're done
      _invalidate();

      std::cout << "  end increment (13)" << std::endl;
    };

    int DiffTransConfigEnumPerturbations::run(const PrimClex &primclex, const jsonParser &_kwargs, const Completer::EnumOption &enum_opt) {

      jsonParser orbitnames;
      if(!_kwargs.get_if(orbitnames, "orbits")) {
        std::cerr << "DiffTransConfigEnumPerturbations currently has no default and requires a correct JSON with a orbits tag within it" << std::endl;
        std::cerr << "Core dump will occur because cannot find proper input" << std::endl;
      }

      jsonParser confignames;
      if(!_kwargs.get_if(confignames, "background_configs")) {
        std::cerr << "DiffTransConfigEnumPerturbations currently has no default and requires a correct JSON with a background_configs tag within it" << std::endl;
        std::cerr << "Core dump will occur because cannot find proper input" << std::endl;
      }

      jsonParser kwargs;
      jsonParser local_cspecs;
      if(_kwargs.get_if(kwargs, "local_cspecs")) {
        local_cspecs = kwargs;
      }
      else if(_kwargs.get_if(kwargs, "local_cspecs_filepath")) {
        //look for filepath given and load it
        fs::path local_cspecs_path = _kwargs.get<std::string>();
        jsonParser tmp {local_cspecs_path};
        local_cspecs = tmp;
      }
      else {
        //look for default local cspecs file.
        fs::path local_cspecs_path = "default/path/to/local_cspecs.json";
        jsonParser tmp {local_cspecs_path};
        local_cspecs = tmp;
      }


      Log &log = primclex.log();
      auto &db_diff_trans_configs = primclex.db<DiffTransConfiguration>();

      Index Ninit = db_diff_trans_configs.size();
      log << "# DiffTransConfiguration in this project: " << Ninit << "\n" << std::endl;

      log.begin(enumerator_name);
      for(const auto &configname_parser : confignames) {
        std::string configname = configname_parser.get<std::string>();
        log << "Searching in " << configname << "..." << std::endl;

        Configuration bg_config = *primclex.db<Configuration>().find(configname);
        for(const auto &orbitname_parser : orbitnames) {
          std::string orbitname = orbitname_parser.get<std::string>();
          log << "\tUsing " << orbitname << "..." << std::flush;
          Index Ninit_spec = db_diff_trans_configs.size() ;

          PrimPeriodicDiffTransOrbit dtorbit = *primclex.db<PrimPeriodicDiffTransOrbit>().find(orbitname);
          /// check if configuration is big enough for local cspecs here
          /// give warning if not

          std::vector<LocalIntegralClusterOrbit> local_orbits;
          make_local_orbits(
            dtorbit.prototype(),
            local_cspecs,
            alloy_sites_filter,
            primclex.crystallography_tol(),
            std::back_inserter(local_orbits),
            null_log());
          if(has_local_bubble_overlap(local_orbits, bg_config.supercell())) {
            log << "WARNING!!! CHOICE OF BACKGROUND CONFIGURATION " << configname <<
                "\nRESULTS IN AN OVERLAP IN THE LOCAL CLUSTERS OF " << orbitname <<
                "\nWITH ITS PERIODIC IMAGES. CONSIDER CHOOSING \n" <<
                "A LARGER BACKGROUND CONFIGURATION." << std::endl;
          }
          DiffTransConfigEnumPerturbations enumerator(bg_config, dtorbit, local_cspecs);
          auto enum_it = enumerator.begin();
          std::vector<std::string> filter_expr = make_enumerator_filter_expr(_kwargs, enum_opt);
          if(!filter_expr.empty()) {
            try {
              DataFormatter<DiffTransConfiguration> filter = primclex.settings().query_handler<DiffTransConfiguration>().dict().parse(filter_expr);
              while(enum_it != enumerator.end()) {
                ValueDataStream<bool> _stream;
                _stream << filter(*enum_it);
                if(_stream.value()) {
                  db_diff_trans_configs.insert(enum_it->canonical_form());
                }
                ++enum_it;
              }
            }
            catch(std::exception &e) {
              primclex.err_log() << "Cannot filter difftransconfigs using the expression provided: \n" << e.what() << "\nExiting...\n";
              return ERR_INVALID_ARG;
            }
          }
          else {
            primclex.db<DiffTransConfiguration>().insert(enumerator.begin(), enumerator.end());
          }

          Index Nfinal_spec = db_diff_trans_configs.size();
          log << "Found " << Nfinal_spec - Ninit_spec << " new difftransconfigs" << std::endl;
        }
      }

      log << "  DONE." << std::endl << std::endl;

      Index Nfinal = db_diff_trans_configs.size();

      log << "# new DiffTransConfiguration: " << Nfinal - Ninit << "\n";
      log << "# DiffTransConfiguration in this project: " << Nfinal << "\n" << std::endl;

      log << "Writing DiffTransConfiguration database..." << std::endl;
      db_diff_trans_configs.commit();
      log << "  DONE" << std::endl;
      return 0;
    }

    ///------------------------------------Internal functions-------------------------------///

    double DiffTransConfigEnumPerturbations::_tol() const {
      return m_background_config.primclex().crystallography_tol();
    }

    /// Base DiffTransConfig supercell
    const Supercell &DiffTransConfigEnumPerturbations::_supercell() const {
      return m_background_config.supercell();
    }

    void DiffTransConfigEnumPerturbations::_init_base() {

      std::cout << "  _init_base 0" << std::endl;

      // Make suborbit generating DiffusionTransformations
      std::vector<DiffusionTransformation> subprototypes;
      make_suborbit_generators(m_diff_trans_orbit, m_background_config, std::back_inserter(subprototypes));
      std::cout << "  _init_base 1" << std::endl;

      // Put them in canonical form, and store similarly transformed
      //   background config. Base constructor will also generate diff trans
      //   invariant group.
      ScelCanonicalGenerator<DiffusionTransformation> gen(_supercell());
      for(const auto &diff_trans : subprototypes) {
        auto canonical_diff_trans = gen(diff_trans);
        m_base.emplace_back(
          canonical_diff_trans,
          copy_apply(gen.to_canonical(), m_background_config));
      }

      std::cout << "  _init_base 2" << std::endl;
      m_base_it = m_base.begin();

      std::cout << "  _init_base 3" << std::endl;

    }

    /// Constructor for data structure holding base diff trans configuration in
    /// their canonical form
    ///
    /// \param _diff_trans DiffTrans in canonical form for the supercell
    /// \param _config Background config transformed appropriately for _diff_trans
    DiffTransConfigEnumPerturbations::Base::Base(
      const DiffusionTransformation &_diff_trans,
      const Configuration &_config) :
      diff_trans(_diff_trans),
      config(make_attachable(diff_trans, _config)),
      diff_trans_g(make_invariant_subgroup(diff_trans, config.supercell())),
      diff_trans_sym_g(make_sym_group(diff_trans_g)) {}

    /// Generate local orbits for current base diff trans
    void DiffTransConfigEnumPerturbations::_init_local_orbits() {

      std::cout << "  _init_local_orbits 0" << std::endl;
      /// Make all local orbits
      std::vector<LocalIntegralClusterOrbit> _tmp;
      make_local_orbits(
        m_base_it->diff_trans,
        m_local_cspecs,
        alloy_sites_filter,
        _tol(),
        std::back_inserter(_tmp),
        null_log(),
        m_base_it->diff_trans_sym_g);
      std::cout << "  _tmp.size(): " << _tmp.size() << std::endl;

      /// Exclude orbits that would alter the hop due small supercell size

      std::cout << "  _init_local_orbits 1" << std::endl;
      // get list of linear indices of hopping sites
      std::set<Index> hop_cluster_indices;
      for(auto &traj : m_base_it->diff_trans.specie_traj()) {
        hop_cluster_indices.insert(_supercell().linear_index(traj.from.uccoord));
      }

      std::cout << "  _init_local_orbits 2" << std::endl;
      // lambda function returns true if uccoord is not in hop cluster
      auto uccoord_does_not_overlap = [&](const UnitCellCoord & uccoord) {
        Index index = this->_supercell().linear_index(uccoord);
        return hop_cluster_indices.find(index) == hop_cluster_indices.end();
      };

      std::cout << "  _init_local_orbits 3" << std::endl;
      // lambda function returns true if no overlap between hop cluster and local cluster
      auto orbit_does_not_overlap = [&](const LocalIntegralClusterOrbit & test) {
        const auto &proto = test.prototype();
        return std::all_of(proto.begin(), proto.end(), uccoord_does_not_overlap);
      };

      std::cout << "  _init_local_orbits 4" << std::endl;
      m_local_orbit.clear();
      std::copy_if(_tmp.begin(), _tmp.end(), std::back_inserter(m_local_orbit), orbit_does_not_overlap);
      std::cout << "  m_local_orbit.size(): " << m_local_orbit.size() << std::endl;
      m_local_orbit_it = m_local_orbit.begin();
      std::cout << "  _init_local_orbits 5" << std::endl;
    }

    /// Generate the 'from_value' for the perturbation,
    ///   the 'to_value' counter for the perturbation,
    ///   and the local orbit prototype invariant subgroup
    ///   (w/ respect to base diff trans invariant group)
    void DiffTransConfigEnumPerturbations::_init_perturbations_data() {

      std::cout << "  _init_perturbations_data 0" << std::endl;
      const Configuration &from_config = m_base_it->config;
      std::cout << "  _init_perturbations_data 0b" << std::endl;
      std::cout << "  m_local_orbit.size(): " << m_local_orbit.size() << std::endl;
      const IntegralCluster &proto = m_local_orbit_it->prototype();

      std::cout << "  _init_perturbations_data 1" << std::endl;
      /// Set 'm_from_value'
      m_from_value = Eigen::VectorXi::Zero(proto.size());
      for(int i = 0; i < proto.size(); ++i) {
        m_from_value(i) = from_config.occ(_supercell().linear_index(proto[i]));
      }

      std::cout << "  _init_perturbations_data 2" << std::endl;
      /// Construct counter
      Eigen::VectorXi max_count(proto.size());
      for(int i = 0; i < proto.size(); ++i) {
        max_count(i) = proto.prim().basis[proto[i].sublat()].site_occupant().size();
      }
      std::cout << "  _init_perturbations_data 3" << std::endl;
      m_occ_counter = EigenCounter<Eigen::VectorXi>(
                        Eigen::VectorXi::Zero(proto.size()),
                        max_count,
                        Eigen::VectorXi::Constant(proto.size(), 1));

      std::cout << "  _init_perturbations_data 4" << std::endl;
      /// Construct local orbit prototype invariant group
      ///   (w/ respect to base diff trans invariant group)
      m_local_orbit_sub_g = make_invariant_subgroup(
                              proto,
                              _supercell(),
                              m_base_it->diff_trans_g.begin(),
                              m_base_it->diff_trans_g.end());
      std::cout << "  _init_perturbations_data 5" << std::endl;
    }

    std::pair<Perturbation, bool> DiffTransConfigEnumPerturbations::_current_perturb() const {

      std::cout << "  _current_perturb 0" << std::endl;
      Perturbation perturb;
      const auto &proto = m_local_orbit_it->prototype();
      bool is_subcluster = false;
      for(int i = 0; i < proto.size(); ++i) {
        if(m_from_value(i) == m_occ_counter()[i] && m_skip_subclusters) {
          is_subcluster = true;
          break;
        }
        perturb.emplace(proto[i], m_from_value(i), m_occ_counter()[i]);
      }

      std::cout << "  _current_perturb 1" << std::endl;
      bool is_valid = !(is_subcluster && m_skip_subclusters)
                      && perturb.is_canonical(m_local_orbit_sub_g.begin(), m_local_orbit_sub_g.end());

      std::cout << "  _current_perturb 2" << std::endl;
      // if canonical perturbation (w/ local orbit sub group)
      return std::make_pair(perturb, is_valid);
    }

    /// Applies current perturbation to m_base_config and stores result in m_current
    void DiffTransConfigEnumPerturbations::_set_current(const Perturbation &perturb) {

      std::cout << "  _set_current 0" << std::endl;
      // generate perturbed from_config
      Configuration perturbed_from_config {m_base_it->config};
      perturb.apply_to(perturbed_from_config);

      std::cout << "  _set_current 1" << std::endl;
      // check perturbed_from_config < test*perturbed_from_config, to find max
      ConfigCompare compare(perturbed_from_config, _tol());
      auto max = m_base_it->diff_trans_g[0] * m_local_orbit_sub_g[0];
      for(const auto &op_i : m_local_orbit_sub_g) {
        for(const auto &op_j : m_base_it->diff_trans_g) {
          //auto test_config = op_j*op_i*(from_config + perturbation)
          auto test = op_j * op_i;
          if(compare(max, test)) {
            max = test;
          }
        }
      }

      std::cout << "  _set_current 2" << std::endl;
      /// construct canonical DiffTransConfiguration as m_current
      m_current = notstd::make_cloneable<DiffTransConfiguration>(
                    copy_apply(max, perturbed_from_config),
                    m_base_it->diff_trans);
      m_current->set_orbit_name(m_diff_trans_orbit.name());
      m_current->set_source(this->source(step()));
      this->_set_current_ptr(&(*m_current));

      std::cout << "  _set_current 3" << std::endl;
      // --- debug check ---
      // m_current should be in canonical form at this point
      if(!m_current->is_canonical()) {

        throw std::runtime_error("Error in DiffTransConfigEnumPerturbations: not canonical");
      }
      std::cout << "  _set_current 4" << std::endl;

    }

    bool has_local_bubble_overlap(std::vector<LocalIntegralClusterOrbit> &local_orbits, const Supercell &scel) {
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

    std::vector<Supercell> viable_supercells(std::vector<LocalIntegralClusterOrbit> &local_orbits, std::vector<Supercell> scel_options) {
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
