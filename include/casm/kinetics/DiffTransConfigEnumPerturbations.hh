#ifndef CASM_DiffTransConfigEnumPerturbations
#define CASM_DiffTransConfigEnumPerturbations

#include "casm/container/InputEnumerator.hh"
#include "casm/container/Counter.hh"
#include "casm/kinetics/DiffusionTransformation.hh"
#include "casm/symmetry/OrbitGeneration.hh"
#include "casm/kinetics/DiffTransEnumEquivalents.hh"
#include "casm/kinetics/DiffTransConfiguration.hh"

namespace CASM {

  namespace Kinetics {

    class Perturbation : public std::set<OccupationTransformation> {

    public:

      Perturbation();

      Perturbation(std::set<OccupationTransformation> &from_set);

      Perturbation &apply_sym(const SymOp &op);

      Configuration &apply_to(Configuration &config) const;

      template<typename PermuteIteratorIt>
      bool is_canonical(PermuteIteratorIt begin, PermuteIteratorIt end) const;
    };


    /// \brief Enumerate DiffTransConfiguration for a particular DiffusionTransformation,
    ///        set of local clusters, and a particular initial Configuration
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
    /// - For each 'base' diff trans, generate local orbits
    /// - For each local orbit, enumerate perturbations. Only include perturbations
    ///   that modify every site in the cluster to avoid repeats.
    /// - The subgroup of the 'diff_trans_g' that leaves the local orbit protype
    ///   invariant is the 'local_orbit_sub_g'. Use the 'local_orbit_sub_g' to
    ///   identify unique perturbations.
    /// - Use the 'local_orbit_sub_g' and 'diff_trans_g' to find the canonical
    ///   form of the 'from_config'.
    /// - Enumerate DiffTransConfiguration using the base diff trans and the
    ///   canonical 'from_config'.
    class DiffTransConfigEnumPerturbations : public InputEnumeratorBase<DiffTransConfiguration> {

      // -- Required members -------------------
      //class Perturbation;
    public:

      /// \brief Construct with an IntegralCluster
      DiffTransConfigEnumPerturbations(
        const Configuration &background_config,
        const PrimPeriodicDiffTransOrbit &diff_trans_orbit, // or const DiffusionTransformation &diff_trans
        const jsonParser &local_cspecs // or iterators over IntegralClusters
      );

      std::string name() const override {
        return enumerator_name;
      }

      static const std::string enumerator_name;
      static const std::string interface_help;

      static int run(const PrimClex &primclex, const jsonParser &_kwargs, const Completer::EnumOption &enum_opt);

    private:

      /// Implements increment: generate the next DiffTransConfiguration
      void increment() override;

      /// Crystallography tolerance
      double _tol() const;

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

      /// Generate the current Perturbation (non-canonical) and whether it is valid
      std::pair<Perturbation, bool> _current_perturb() const;

      /// Apply perturbation, find canonical 'from_config' and set current DiffTransConfiguration
      void _set_current(const Perturbation &perturb);


      /// The background configuration which will be perturbed
      const Configuration m_background_config;

      /// The DiffusingTransformation to be applied in all symmetrically unique
      /// places in the background config
      const PrimPeriodicDiffTransOrbit m_diff_trans_orbit;

      /// Avoid repeating perturbations
      bool m_skip_subclusters;

    public:
      /// Stores base DiffTrans and Config to enable construction of DiffTransConfig
      /// in canonical form
      struct Base {

        Base(const DiffusionTransformation &_diff_trans, const Configuration &_config);

        /// sub-orbit prototype of diff_trans_orbit in the canonical position in _supercell()
        DiffusionTransformation diff_trans;

        /// m_background_config transformed in the same manner as diff_trans
        Configuration config;

        /// The subgroup of _supercell() leaving diff_trans invariant
        std::vector<PermuteIterator> diff_trans_g;

        /// Alternative representation of diff_trans_g
        SymGroup diff_trans_sym_g;
      };

    private:
      /// Vector of bases to be perturbed
      /// - Each represents the DiffusionTransformation in unique local
      ///   environment of the background config, with the DiffTrans in its
      ///   canonical form inside the supercell
      std::vector<Base> m_base;

      /// Current base to be perturbed
      std::vector<Base>::iterator m_base_it;

      /// Local orbit specs
      const jsonParser m_local_cspecs;

      /// Local orbits for the current base diff trans
      std::vector<LocalIntegralClusterOrbit> m_local_orbit;

      /// The current local orbit
      std::vector<LocalIntegralClusterOrbit>::iterator m_local_orbit_it;

      /// The subgroup of the current diff trans group that leaves the local orbit prototype invariant
      /// - i.e. The subgroup of m_base_it->diff_trans_g that leaves
      ///        m_local_orbit_it->prototype() invariant
      std::vector<PermuteIterator> m_local_orbit_sub_g;

      /// The current 'from_value' for perturbations
      /// - determined from current base config and current local orbit
      Eigen::VectorXi m_from_value;

      /// Counts over 'to_value' for the perturbation
      EigenCounter<Eigen::VectorXi> m_occ_counter;

      /// The base DiffTransConfiguration with the current perturbation applied
      notstd::cloneable_ptr<DiffTransConfiguration> m_current;
    };

    bool has_local_bubble_overlap(std::vector<LocalIntegralClusterOrbit> &local_orbits, const Supercell &scel);

    std::vector<Supercell> viable_supercells(std::vector<LocalIntegralClusterOrbit> &local_orbits, std::vector<Supercell>);


  }
}
#endif
