#ifndef CASM_DiffTransConfigEnumOccPerturbations
#define CASM_DiffTransConfigEnumOccPerturbations

#include <iterator>
#include "casm/container/InputEnumerator.hh"
#include "casm/container/Counter.hh"
#include "casm/symmetry/Orbit.hh"
#include "casm/symmetry/ScelOrbitGeneration.hh"
#include "casm/kinetics/DiffusionTransformation.hh"
#include "casm/kinetics/OccPerturbation.hh"
#include "casm/kinetics/DiffTransConfiguration.hh"
#include "casm/clusterography/ClusterSymCompare.hh"

namespace CASM {

  namespace Kinetics {

    /// \brief Enumerate DiffTransConfiguration for a particular DiffusionTransformation,
    ///        set of local clusters, and a particular initial Configuration
    ///
    /// Note:
    /// - DiffTransConfiguration output are unique and canonical
    ///
    /// Algorithm:
    /// - Make suborbit generators of the specified diff_trans_orbit due to
    ///   symmetry breaking by the background_config
    /// - Find the operation that transforms the suborbit generators into their
    ///   canonical form in the supercell. Transform both the suborbit generators
    ///   and the background config accordingly, and save each pair. These are
    ///   the 'base' DiffTransConfigurations that will be perturbed. The subgroup
    ///   of the supercell factor group that leaves the 'base' diff trans invariant
    ///   is the 'diff_trans_g'.
    /// - For each 'base' diff trans, generate local orbits using the base config
    ///   & diff trans invariant group. Note: local orbits are actually
    ///   ScelPeriodicOrbit<IntegralCluster> because we take supercell symmetry
    ///   into account. This means there may be more perturbations then the local
    ///   basis functions may indicate.
    /// - For each local orbit, enumerate perturbations. Only include perturbations
    ///   that modify every site in the cluster to avoid repeats.
    /// - The subgroup of the 'diff_trans_g' that leaves the local orbit protype
    ///   invariant is the 'local_orbit_sub_g'. Use the 'local_orbit_sub_g' to
    ///   identify unique perturbations.
    /// - Use the 'local_orbit_sub_g' and 'diff_trans_g' to find the canonical
    ///   form of the 'from_config'.
    /// - Enumerate DiffTransConfiguration using the base diff trans and the
    ///   canonical 'from_config'.
    class DiffTransConfigEnumOccPerturbations : public InputEnumeratorBase<DiffTransConfiguration> {

      // -- Required members -------------------
      //class OccPerturbation;
    public:

      /// \brief Construct with an IntegralCluster
      DiffTransConfigEnumOccPerturbations(
        const Configuration &background_config,
        const PrimPeriodicDiffTransOrbit &diff_trans_orbit, // or const DiffusionTransformation &diff_trans
        const jsonParser &local_cspecs // or iterators over IntegralClusters
      );

      std::string name() const override {
        return enumerator_name;
      }

      static const std::string enumerator_name;
      static const std::string interface_help;

      /// Enumerate DiffTransConfigEnumOccPerturbations into the project database
      static int run(
        const PrimClex &primclex,
        const jsonParser &_kwargs,
        const Completer::EnumOption &enum_opt);

      /// Enumerate DiffTransConfigEnumOccPerturbations into any std::set-like database
      template<typename DatabaseType>
      static int run(
        const PrimClex &primclex,
        const jsonParser &kwargs,
        const Completer::EnumOption &enum_opt,
        DatabaseType &db);

      /// Stores base DiffTrans and Config to enable construction of DiffTransConfig
      /// in canonical form
      struct Base {

        Base(const DiffusionTransformation &_diff_trans, const Configuration &_config);

        /// sub-orbit prototype of diff_trans_orbit in the canonical position in _supercell()
        DiffusionTransformation diff_trans;

        /// m_background_config transformed in the same manner as diff_trans
        Configuration config;

        /// The subgroup of supercell permutations leaving diff_trans invariant
        std::vector<PermuteIterator> diff_trans_g;

        /// The subgroup of supercell permutations leaving config and diff_trans invariant
        /// - Generating group for local clusters / perturbations
        std::vector<PermuteIterator> generating_g;

        /// Alternative representation of generating_g
        SymGroup generating_sym_g;
      };

      /// Stores the current OccPertubation (non-canonical) and its validity
      struct CurrentPerturbation {
        CurrentPerturbation(const OccPerturbation &_perturb) :
          perturb(_perturb),
          is_not_subcluster(false),
          is_canonical(false) {}

        OccPerturbation perturb;
        bool is_not_subcluster;
        bool is_canonical;
      };

      const std::vector<Base> &base() const;

      Index base_index() const;

      const std::vector<ScelPeriodicOrbit<IntegralCluster>> &local_orbit() const;

      Index local_orbit_index() const;

      Index occ_counter_index() const;

      void partial_increment(bool complete_perturb = false);

      bool check_increment();

      /// Return the current OccPerturbation (non-canonical) and whether it is valid
      const CurrentPerturbation &current_perturb() const;

    private:

      /// Implements increment: generate the next DiffTransConfiguration
      void increment() override;

      /// Crystallography tolerance
      double _tol() const;

      /// Supercell
      const Structure &_prim() const;

      /// Supercell
      const Supercell &_supercell() const;

      /// Determines unique canonical diff trans for given background config
      void _init_base();

      /// Generate local orbits for current base diff trans
      void _init_local_orbits();

      /// Generate the 'from_value' for the perturbation,
      ///   the 'to_value' counter for the perturbation,
      ///   and the local orbit prototype invariant subgroup
      ///   (w/ respect to base diff trans invariant group)
      void _init_perturbations_data();

      /// Generate the current OccPerturbation (non-canonical) and whether it is valid
      ///
      /// - If not valid, the OccPerturbation itself may not be complete
      void _update_current_perturb(bool complete_perturb = false);

      /// Apply perturbation, find canonical 'from_config' and set current DiffTransConfiguration
      void _set_current(const OccPerturbation &perturb);


      /// The background configuration which will be perturbed
      const Configuration m_background_config;

      /// The DiffusingTransformation to be applied in all symmetrically unique
      /// places in the background config
      const PrimPeriodicDiffTransOrbit m_diff_trans_orbit;

      /// Include the base DiffTransConfiguration in the output
      bool m_include_unperturbed;

      /// Avoid repeating perturbations
      bool m_skip_subclusters;

      /// Vector of bases to be perturbed
      /// - Each represents the DiffusionTransformation in unique local
      ///   environment of the background config, with the DiffTrans in its
      ///   canonical form inside the supercell
      std::vector<Base> m_base;

      /// Current base to be perturbed
      std::vector<Base>::iterator m_base_it;

      /// Local orbit specs
      const jsonParser m_local_cspecs;

      /// SymCompare to translate local orbits into the supercell
      ScelPeriodicSymCompare<IntegralCluster> m_scel_sym_compare;

      /// Local orbits for the current base diff trans
      std::vector<ScelPeriodicOrbit<IntegralCluster>> m_local_orbit;

      /// The current local orbit
      std::vector<ScelPeriodicOrbit<IntegralCluster>>::iterator m_local_orbit_it;

      /// The subgroup of the current diff trans group that leaves the local orbit prototype invariant
      /// - i.e. The subgroup of m_base_it->diff_trans_g that leaves
      ///        m_local_orbit_it->prototype() invariant
      std::vector<PermuteIterator> m_local_orbit_sub_g;

      /// The current 'from_value' for perturbations
      /// - determined from current base config and current local orbit
      Eigen::VectorXi m_from_value;

      /// Counts over 'to_value' for the perturbation
      EigenCounter<Eigen::VectorXi> m_occ_counter;

      /// Counts over m_occ_counter
      Index m_occ_counter_index;

      /// The base DiffTransConfiguration with the current perturbation applied
      notstd::cloneable_ptr<DiffTransConfiguration> m_current;

      /// Current OccPerturbation and its validity
      CurrentPerturbation m_curr;

    };

  }
}
#endif
