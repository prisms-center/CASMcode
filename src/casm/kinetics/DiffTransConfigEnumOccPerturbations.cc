#include "casm/kinetics/DiffTransConfigEnumOccPerturbations.hh"

#include "casm/clex/PrimClex.hh"
#include "casm/clex/Supercell.hh"
#include "casm/clex/Configuration.hh"
#include "casm/clex/ConfigCompare.hh"
#include "casm/clex/ConfigIsEquivalent.hh"
#include "casm/kinetics/DiffusionTransformation.hh"

#include "casm/app/ProjectSettings.hh"
#include "casm/app/AppIO.hh"
#include "casm/app/QueryHandler.hh"
#include "casm/database/Selection.hh"
#include "casm/database/ConfigDatabase.hh"
#include "casm/database/DiffTransConfigDatabase.hh"
#include "casm/database/DiffTransOrbitDatabase.hh"

#include "casm/app/AppIO_impl.hh"
#include "casm/symmetry/ConfigSubOrbits_impl.hh"
#include "casm/symmetry/ScelOrbitGeneration_impl.hh"
#include "casm/clusterography/ClusterOrbits_impl.hh"
#include "casm/clusterography/ClusterSymCompare_impl.hh"


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
      m_scel_sym_compare(_supercell().prim_grid(), _tol()),
      m_include_unperturbed(true),
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

    const std::string DiffTransConfigEnumOccPerturbations::enumerator_name = "DiffTransConfigEnumOccPerturbations";
    const std::string DiffTransConfigEnumOccPerturbations::interface_help =
      "DiffTransConfigEnumOccPerturbations: \n\n"

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
      "  }\n\n";

    void DiffTransConfigEnumOccPerturbations::increment() {

      std::cout << "begin increment" << std::endl;
      // increment occ_counter until a canonical perturb is found
      // if no more occupations, increment m_local_orbit_it
      // if no more local orbits, increment m_base_it
      // if no more base diff trans, _invalidate

      // current Perturbation & 'is_valid'
      std::pair<OccPerturbation, bool> curr_perturb {OccPerturbation(_prim()), false};

      do {
        std::cout << "-- increment iteration" << std::endl;
        ++m_occ_counter;
        if(m_occ_counter.valid()) {
          std::cout << "  next perturbation" << std::endl;
          curr_perturb = _current_perturb();
        }
        else if(++m_local_orbit_it != m_local_orbit.end()) {
          std::cout << "  next local orbit: " << std::distance(m_local_orbit.begin(), m_local_orbit_it)
                    << "/" << m_local_orbit.size() << std::endl;
          _init_perturbations_data();
          curr_perturb = _current_perturb();
        }
        else if(++m_base_it != m_base.end()) {
          std::cout << "  next base: " << std::distance(m_base.begin(), m_base_it)
                    << "/" << m_base.size() << std::endl;
          _init_local_orbits();
          _init_perturbations_data();
          curr_perturb = _current_perturb();
        }
        std::cout << "  curr_perturb.second: " << curr_perturb.second << std::endl;
        std::cout << "  m_base_it != m_base.end(): " << (m_base_it != m_base.end()) << std::endl;
        std::cout << "  check: " << (!curr_perturb.second && m_base_it != m_base.end()) << std::endl;
      }
      while(!curr_perturb.second && m_base_it != m_base.end());


      if(curr_perturb.second) {
        std::cout << "  incr: is valid perturb" << std::endl;
        // find canonical from_config set current
        this->_increment_step();
        _set_current(curr_perturb.first);
      }
      else {
        std::cout << "  invalidate" << std::endl;
        // if no more base diff trans, we're done
        _invalidate();
      }
      std::cout << "  end increment" << std::endl;
    };

    int DiffTransConfigEnumOccPerturbations::run(const PrimClex &primclex, const jsonParser &_kwargs, const Completer::EnumOption &enum_opt) {

      jsonParser orbitnames;
      if(!_kwargs.get_if(orbitnames, "orbits")) {
        std::cerr << "DiffTransConfigEnumOccPerturbations currently has no default and requires a correct JSON with a orbits tag within it" << std::endl;
        std::cerr << "Core dump will occur because cannot find proper input" << std::endl;
      }

      jsonParser confignames;
      if(!_kwargs.get_if(confignames, "background_configs")) {
        std::cerr << "DiffTransConfigEnumOccPerturbations currently has no default and requires a correct JSON with a background_configs tag within it" << std::endl;
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

          std::vector<LocalOrbit<IntegralCluster>> local_orbits;
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
          DiffTransConfigEnumOccPerturbations enumerator(bg_config, dtorbit, local_cspecs);
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

      std::cout << "  begin _init_base" << std::endl;

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
          canonical_diff_trans,
          copy_apply(gen.to_canonical(), m_background_config));
      }

      std::cout << "  m_base.size(): " << m_base.size() << std::endl;
      std::cout << "Base: " << std::endl;
      Index index = 0;
      for(const auto &base : m_base) {
        std::cout << "  Index: " << index << std::endl;
        std::cout << "    DiffTrans: \n" << base.diff_trans << std::endl;
        ++index;
      }
      m_base_it = m_base.begin();

      std::cout << "  end _init_base" << std::endl;
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
      diff_trans_g(make_invariant_subgroup(diff_trans, config.supercell())),
      diff_trans_sym_g(make_sym_group(diff_trans_g)) {}

    /// Generate local orbits for current base diff trans
    void DiffTransConfigEnumOccPerturbations::_init_local_orbits() {

      std::cout << "  begin _init_local_orbits" << std::endl;
      /// Make all local orbits
      std::vector<LocalOrbit<IntegralCluster>> _tmp;
      make_local_orbits(
        m_base_it->diff_trans,
        m_local_cspecs,
        alloy_sites_filter,
        _tol(),
        std::back_inserter(_tmp),
        null_log(),
        m_base_it->diff_trans_sym_g);

      /// Exclude orbits that would alter the hop due small supercell size

      // get list of linear indices of hopping sites
      std::set<Index> hop_cluster_indices;
      for(auto &traj : m_base_it->diff_trans.specie_traj()) {
        hop_cluster_indices.insert(_supercell().linear_index(traj.from.uccoord));
      }

      // lambda function returns true if uccoord is not in hop cluster
      auto uccoord_does_not_overlap = [&](const UnitCellCoord & uccoord) {
        Index index = this->_supercell().linear_index(uccoord);
        return hop_cluster_indices.find(index) == hop_cluster_indices.end();
      };

      // lambda function returns true if no overlap between hop cluster and local cluster
      auto orbit_does_not_overlap = [&](const LocalOrbit<IntegralCluster> &test) {
        const auto &proto = test.prototype();
        return std::all_of(proto.begin(), proto.end(), uccoord_does_not_overlap);
      };

      m_local_orbit.clear();
      std::copy_if(_tmp.begin(), _tmp.end(), std::back_inserter(m_local_orbit), orbit_does_not_overlap);
      std::cout << "  m_local_orbit.size(): " << m_local_orbit.size() << std::endl;
      m_local_orbit_it = m_local_orbit.begin();

      std::cout << "  m_local_orbits: \n" << std::endl;
      print_clust(
        m_local_orbit.begin(),
        m_local_orbit.end(),
        std::cout,
        PrototypePrinter<IntegralCluster>());
      std::cout << "  end _init_local_orbits" << std::endl;
    }

    /// Generate the 'from_value' for the perturbation,
    ///   the 'to_value' counter for the perturbation,
    ///   and the local orbit prototype invariant subgroup
    ///   (w/ respect to base diff trans invariant group)
    void DiffTransConfigEnumOccPerturbations::_init_perturbations_data() {

      std::cout << "  begin _init_perturbations_data" << std::endl;
      const Configuration &from_config = m_base_it->config;
      const IntegralCluster &proto = m_local_orbit_it->prototype();

      std::cout << "    orbit index: " << std::distance(m_local_orbit.begin(), m_local_orbit_it)
                << "  branch: " << proto.size() << std::endl;

      /// Set 'm_from_value'
      m_from_value = Eigen::VectorXi::Zero(proto.size());
      for(int i = 0; i < proto.size(); ++i) {
        m_from_value(i) = from_config.occ(_supercell().linear_index(proto[i]));
      }

      /// Construct counter
      Eigen::VectorXi max_count(proto.size());
      for(int i = 0; i < proto.size(); ++i) {
        max_count(i) = proto.prim().basis[proto[i].sublat()].site_occupant().size() - 1;
      }
      m_occ_counter = EigenCounter<Eigen::VectorXi>(
                        Eigen::VectorXi::Zero(proto.size()),
                        max_count,
                        Eigen::VectorXi::Constant(proto.size(), 1));

      /// Construct local orbit prototype invariant group
      ///   (w/ respect to base diff trans invariant group)
      m_local_orbit_sub_g = make_invariant_subgroup(
                              proto,
                              _supercell(),
                              m_base_it->diff_trans_g.begin(),
                              m_base_it->diff_trans_g.end());

      std::cout << "  end _init_perturbations_data" << std::endl;
    }

    std::pair<OccPerturbation, bool> DiffTransConfigEnumOccPerturbations::_current_perturb() const {

      OccPerturbation perturb {_prim()};
      const auto &proto = m_scel_sym_compare.prepare(m_local_orbit_it->prototype());

      std::cout << "  begin _current_perturb" << std::endl;
      std::cout << "  from: " << m_from_value.transpose() << std::endl;
      std::cout << "  to: " << m_occ_counter().transpose() << std::endl;

      // check for null perturbation
      if(proto.size() == 0 && m_include_unperturbed) {
        std::cout << "  perturb: (size: " << perturb.size() << ")\n" << perturb << std::endl;
        return std::make_pair(perturb, true);
      }

      // generate next perturbation consistent with from_config
      for(int i = 0; i < proto.size(); ++i) {
        // check for subcluster perturbation
        if(m_from_value(i) == m_occ_counter()[i] && m_skip_subclusters) {
          std::cout << "  not valid Perturbation: is subcluster" << std::endl;
          return std::make_pair(perturb, false);
        }
        perturb.elements().emplace_back(proto[i], m_from_value(i), m_occ_counter()[i]);
      }

      auto begin = m_local_orbit_sub_g.begin();
      auto end = m_local_orbit_sub_g.end();

      std::cout << "  perturb: (size: " << perturb.size() << ")\n" << perturb << std::endl;
      std::cout << "  is_canonical: " << perturb.is_canonical(_supercell(), begin, end) << std::endl;

      // return perturbation and check if canonical wrt local orbit sub group
      return std::make_pair(perturb, perturb.is_canonical(_supercell(), begin, end));
    }

    /// Applies current perturbation to m_base_config and stores result in m_current
    void DiffTransConfigEnumOccPerturbations::_set_current(const OccPerturbation &perturb) {

      std::cout << "  begin _set_current" << std::endl;
      // generate perturbed from_config
      Configuration perturbed_from_config {m_base_it->config};
      perturb.apply_to(perturbed_from_config);

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

      /// construct canonical DiffTransConfiguration as m_current
      m_current = notstd::make_cloneable<DiffTransConfiguration>(
                    copy_apply(max, perturbed_from_config),
                    m_base_it->diff_trans);
      m_current->set_orbit_name(m_diff_trans_orbit.name());
      m_current->set_source(this->source(step()));
      this->_set_current_ptr(&(*m_current));

      // --- debug check ---
      // m_current should be in canonical form at this point
      if(!m_current->is_canonical()) {

        throw std::runtime_error("Error in DiffTransConfigEnumOccPerturbations: not canonical");
      }
      std::cout << "  end _set_current" << std::endl;

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
