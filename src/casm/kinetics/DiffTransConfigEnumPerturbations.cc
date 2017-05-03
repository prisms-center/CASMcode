#include "casm/kinetics/DiffusionTransformation.hh"
#include "casm/clex/Configuration.hh"
#include "casm/kinetics/DiffTransEnumEquivalents.hh"
#include "casm/kinetics/DiffTransConfigEnumPerturbations.hh"
#include "casm/clex/Supercell.hh"
#include "casm/clusterography/ClusterOrbits.hh"
#include "casm/clusterography/ClusterOrbits_impl.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/app/AppIO.hh"
#include "casm/app/AppIO_impl.hh"

extern "C" {
  CASM::EnumInterfaceBase *make_DiffTransConfigEnumPerturbations_interface() {
    return new CASM::EnumInterface<CASM::Kinetics::DiffTransConfigEnumPerturbations>();
  }
}

namespace CASM {

  namespace Kinetics {

    class Perturbation;

    DiffTransConfigEnumPerturbations::DiffTransConfigEnumPerturbations(
      const Configuration &background_config,
      const PrimPeriodicDiffTransOrbit &diff_trans_orbit, // or const DiffusionTransformation &diff_trans
      const jsonParser &local_bspecs // or iterators over IntegralClusters
    ) :
      m_background_config(background_config), m_diff_trans_orbit(diff_trans_orbit),
      m_local_bspecs(local_bspecs), m_base_config(background_config, diff_trans_orbit.prototype()) {

      //Enumerate unique diffusion transformation set for this background config
      _init_unique_difftrans();


      if(!m_unique_difftrans.size()) {
        _invalidate();
      }

      //Pick the first base DiffTransConfiguration
      _init_base_config();

      //Initialize perturbation set for m_base_dtc
      _init_perturbations();

      if(!m_perturbations.size()) {
        _invalidate();
      }

      //Set current DiffTransConfiguration
      _set_current();


      if(false) {
        increment();
      }
      else {
        this-> _initialize(&(*m_current));
        _set_step(0);
      }
    }

    const std::string DiffTransConfigEnumPerturbations::enumerator_name = "DiffTransConfigEnumPerturbations";
    const std::string DiffTransConfigEnumPerturbations::interface_help =
      "DiffTransConfigEnumPerturbations: \n\n"

      "  orbits: JSON array of strings \n"
      "    Indicate which diffusion transformation orbits are of interest. The \n"
      "    JSON array \"orbits\" should be the names of the orbits of interest.\n"
      "              \n"
      "  background_configs: JSON array of strings \n "
      "    Indicate which configurations will be the background structures for the transformations to occur in.\n"
      "    The JSON array of strings \"background_configs\" should be names of Configurations\n"
      "    that exist in your CASM project.\n\n"
      ""
      "  local_bspecs: JSON object (optional,default= {}) \n "
      "    Indicate the local bspecs hat indicate the clusters around \n"
      "    the transformation that should be perturbed. The string \"local_bspecs\" should be a local bspecs style\n"
      "    initialization used in casm bset enumeration.\n "
      "    This option takes precedence over the following option.\n\n"
      ""
      ""
      "  local_bspecs_filepath: string (optional,default=\"\") \n "
      "    Indicate the local bspecs file that indicates the clusters around \n"
      "    the transformation that should be perturbed. The string \"local_bspecs_filepath\" should be the file path\n"
      "    to the JSON file containing the local bspecs.\n "
      ""
      "  Example:\n"
      "  {\n"
      "   \"orbits\":[\"DT0\",\"DT1\"],\n"
      "   \"background_configs\":[\"SCEL8_2_2_2_0_0_0/2\"],\n"
      "    \"local_bspecs\":{\n"
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
      if(m_perturb_it == m_perturbations.end()) {
        _increment_base_config();
        _init_perturbations();
        m_perturb_it = m_perturbations.begin();
      }

      _set_current();
      m_perturb_it++;

    };

    int DiffTransConfigEnumPerturbations::run(PrimClex &primclex, const jsonParser &_kwargs, const Completer::EnumOption &enum_opt) {

      jsonParser orbitnames;
      if(!_kwargs.get_if(orbitnames, "orbits")) {
        std::cerr << "DiffTransConfigEnumPerturbations currently has no default and requires a correct JSON with a orbits tag within it" << std::endl;
        std::cerr << "Core dump will occur because cannot find proper input" << std::endl;
      }

      jsonParser confignames;
      if(!_kwargs.get_if(orbitnames, "background_configs")) {
        std::cerr << "DiffTransConfigEnumPerturbations currently has no default and requires a correct JSON with a background_configs tag within it" << std::endl;
        std::cerr << "Core dump will occur because cannot find proper input" << std::endl;
      }

      jsonParser kwargs;
      jsonParser local_bspecs;
      if(_kwargs.get_if(kwargs, "local_bspecs")) {
        jsonParser local_bspecs {kwargs};

      }
      else if(_kwargs.get_if(kwargs, "local_bspecs_filepath")) {
        //look for filepath given and load it
        jsonParser local_bspecs {kwargs};
      }
      else {
        //look for default local bspecs file.
        // Find default path string here
        //jsonParser local_bspecs {kwargs.get()};
      }
      for(const auto &configname : confignames) {
        for(const auto &orbitname : orbitnames) {
          /// PrimPeriodicDiffTransOrbit dtorbit = select(orbitname);
          /// Configuration bg_config = select(configname);
          /// check if configuration is big enough for local bspecs here
          /// give warning if not
          /// DiffTransConfigEnumPerturbations enumerator(bg_config,dtorbit,local_bspecs)
          /// while (enumerator.is_valid()){
          ///   if (enumerator.current() passes filter){
          ///     enumerator.store()
          ///   }
          ///
          ///   ++enumerator;
          ///}
        }
      }
      return 0;
    }

    ///------------------------------------Internal functions-------------------------------///

    ///sets the first m_base_config as the orbit prototype placed in m_background_config
    void DiffTransConfigEnumPerturbations::_init_base_config() {
      Configuration bg_config = make_attachable(*(m_unique_difftrans.begin()), m_background_config);
      DiffTransConfiguration tmp(bg_config, *(m_unique_difftrans.begin()));
      m_base_config = tmp.canonical_form();
      ++m_unique_difftrans_it;
      return;
    }

    /// Uses DiffTransEnumEquivalents to initialize set of the unique diffusion transformations in
    /// this configuration
    void DiffTransConfigEnumPerturbations::_init_unique_difftrans() {
      /// Need to find sub prototypes
      /* Inefficient but correct
        for (auto &e : m_diff_trans_orbit){
          subprototypes.push_back(e);
        }
      */
      std::vector<DiffusionTransformation> subprototypes;
      subprototypes.push_back(m_diff_trans_orbit.prototype());
      for(auto it = subprototypes.begin(); it != subprototypes.end(); ++it) {
        auto begin = m_background_config.supercell().permute_begin();
        auto end = m_background_config.supercell().permute_end();
        DiffTransEnumEquivalents diff_trans_unique(*it, begin, end, m_background_config);
        m_unique_difftrans.insert(diff_trans_unique.begin(), diff_trans_unique.end());
      }
      return;
    }

    /// Uses m_local_bspecs to initialize set of the possible perturbations in
    /// this configuration
    void DiffTransConfigEnumPerturbations::_init_perturbations() {
      m_perturbations.clear();
      const SymGroup &scel_grp = m_base_config.from_config().supercell().factor_group();
      Kinetics::ScelPeriodicDiffTransSymCompare dt_sym_compare(m_base_config.from_config().supercell().prim_grid(),
                                                               m_base_config.from_config().primclex().crystallography_tol());
      SymGroup generating_group = invariant_subgroup(m_base_config.diff_trans(), scel_grp, dt_sym_compare);
      std::vector<LocalIntegralClusterOrbit> local_orbits;
      make_local_orbits(
        m_base_config.diff_trans(),
        m_local_bspecs,
        alloy_sites_filter,
        m_base_config.from_config().primclex().crystallography_tol(),
        std::back_inserter(local_orbits),
        m_base_config.from_config().primclex().log(),
        generating_group);

      for(const auto &orbit : local_orbits) {
        Eigen::VectorXi max_count(orbit.prototype().size());
        for(int i = 0; i < max_count.size(); ++i) {
          max_count(i) = m_base_config.from_config().supercell().max_allowed_occupation()
                         [m_base_config.from_config().supercell().linear_index(orbit.prototype()[i])];
        }
        EigenCounter<Eigen::VectorXi> occ_counter(Eigen::VectorXi::Zero(orbit.prototype().size()),
                                                  max_count, Eigen::VectorXi::Constant(orbit.prototype().size(), 1));
        do {
          std::set<OccupationTransformation> my_occ_transforms;
          for(int i = 0; i < orbit.prototype().size(); ++i) {
            my_occ_transforms.emplace(orbit.prototype()[i],
                                      m_base_config.from_config().occ(m_base_config.from_config().supercell().linear_index(orbit.prototype()[i])),
                                      occ_counter()[i]);
          }
          Perturbation proto_perturb(my_occ_transforms);
          Perturbation max_perturbation = proto_perturb;
          for(auto &op : generating_group) {
            if(copy_apply(op, proto_perturb) > max_perturbation) {
              max_perturbation = copy_apply(op, proto_perturb);
            }
          }
          m_perturbations.insert(max_perturbation);
        }
        while(++occ_counter);

      }
      m_perturb_it = m_perturbations.begin();
      return;
    }

    /// Applies current perturbation to m_base_config and stores result in m_current
    void DiffTransConfigEnumPerturbations::_set_current() {
      /// apply_perturbation(perturb,m_base_dtc);
      Configuration tmp {m_base_config.from_config()};
      for(const auto &occ_trans : *m_perturb_it) {
        occ_trans.apply_to(tmp);
      }
      DiffTransConfiguration ret_dtc(tmp, m_base_config.diff_trans());
      *m_current = ret_dtc;
      return;
    }

    /// Moves to next unique diffusion transformation and places in m_background_config
    void DiffTransConfigEnumPerturbations::_increment_base_config() {
      if(m_unique_difftrans_it != m_unique_difftrans.end()) {
        Configuration bg_config = make_attachable(*m_unique_difftrans_it, m_background_config);
        DiffTransConfiguration tmp(bg_config, *m_unique_difftrans_it);
        m_base_config = tmp.canonical_form();
        ++m_unique_difftrans_it;
      }
      else {
        _invalidate();
      }
      return;
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
