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

    DiffTransConfigEnumPerturbations::DiffTransConfigEnumPerturbations(
      const Configuration &background_config,
      const PrimPeriodicDiffTransOrbit &diff_trans_orbit, // or const DiffusionTransformation &diff_trans
      const jsonParser &local_bspecs // or iterators over IntegralClusters
    ) :
      m_background_config(background_config), m_diff_trans_orbit(diff_trans_orbit), m_local_bspecs(local_bspecs), m_base_dtc(background_config, diff_trans_orbit.prototype()) {

      //Enumerate unique diffusion transformation set for this background config
      _init_unique_dts();

      //Pick the first base DiffTransConfiguration
      _init_base_dtc();

      //Initialize generic perturbation set
      _init_perturbations();

      //Set current DiffTransConfiguration
      _set_current();
      /*
      if (!m_current->is_valid()){
        increment();
      }

      if (!m_perturbations.size()){
        _invalidate();
      }

      else{
        this-> _initialize(&(*m_current));
        _set_step(0);
      }*/
    }

    const std::string DiffTransConfigEnumPerturbations::enumerator_name = "DiffTransConfigEnumPerturbations";
    const std::string DiffTransConfigEnumPerturbations::interface_help =
      "DiffTransConfigEnumPerturbations: \n\n"

      "  orbits: JSON array of strings \n"
      "    Indicate which diffusion transformation orbits are of interest. The \n"
      "    JSON array \"orbits\" should be the names of the orbits of interest.\n"
      "              \n\n"
      "  background_configs: JSON array of strings \n "
      "    Indicate which configurations will be the background structures for the transformations to occur in.\n"
      "    The JSON array of strings \"background_configs\" should be names of Configurations\n"
      "    that exist in your CASM project.\n"
      ""
      "  local_bspecs: string (optional,default="") \n "
      "    Indicate the local bspecs file that indicates the clusters around \n"
      "    the transformation that should be perturbed. The string \"local_bspecs\" should be the file path\n"
      "    to the JSON file containing the local bspecs.\n "
      ""
      ;

    void DiffTransConfigEnumPerturbations::increment() {
      ///if (m_perturb_it == m_perturbations.end()){
      /// _increment_base_dtc();
      /// _init_perturbations();
      /// m_perturb_it = m_perturbations.begin();
      ///}

      ///_set_current();
      ///m_perturb_it++;

    };

    int DiffTransConfigEnumPerturbations::run(PrimClex &primclex, const jsonParser &kwargs, const Completer::EnumOption &enum_opt) {
      return 0;
    }

    ///------------------------------------Internal functions-------------------------------///

    ///sets the first m_base_dtc as the orbit prototype placed in m_background_config
    void DiffTransConfigEnumPerturbations::_init_base_dtc() {
      Configuration bg_config = make_attachable(*(m_unique_dts.begin()), m_background_config);
      DiffTransConfiguration tmp(bg_config, *(m_unique_dts.begin()));
      m_base_dtc = tmp.canonical_form();
      ++m_unique_dts_it;
      return;
    }

    /// Uses DiffTransEnumEquivalents to initialize set of the unique diffusion transformations in
    /// this configuration
    void DiffTransConfigEnumPerturbations::_init_unique_dts() {
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
        m_unique_dts.insert(diff_trans_unique.begin(), diff_trans_unique.end());
      }
      return;
    }

    /// Uses m_local_bspecs to initialize set of the possible perturbations in
    /// this configuration
    void DiffTransConfigEnumPerturbations::_init_perturbations() {
      std::vector<LocalIntegralClusterOrbit> local_orbits;
      make_local_orbits(
        m_base_dtc.diff_trans(),
        m_local_bspecs,
        alloy_sites_filter,
        m_base_dtc.from_config().primclex().crystallography_tol(),
        std::back_inserter(local_orbits),
        m_base_dtc.from_config().primclex().log());
      /* change clusters to perturbations vector somehow */
      return;
    }

    /// Applies current perturbation to m_base_dtc and stores result in m_current
    void DiffTransConfigEnumPerturbations::_set_current() {
      /// apply_perturbation(perturb,m_base_dtc);
      return;
    }

    /// Moves to next unique diffusion transformation and places in m_background_config
    void DiffTransConfigEnumPerturbations::_increment_base_dtc() {
      if(m_unique_dts_it != m_unique_dts.end()) {
        Configuration bg_config = make_attachable(*m_unique_dts_it, m_background_config);
        DiffTransConfiguration tmp(bg_config, *m_unique_dts_it);
        m_base_dtc = tmp.canonical_form();
        ++m_unique_dts_it;
      }
      else {
        _invalidate();
      }
      return;
    }

  }
}
