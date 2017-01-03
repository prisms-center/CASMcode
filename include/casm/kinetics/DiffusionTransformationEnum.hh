#ifndef CASM_DiffusionTransformationEnum
#define CASM_DiffusionTransformationEnum

#include "casm/container/Counter.hh"
#include "casm/container/InputEnumerator.hh"
#include "casm/kinetics/DiffusionTransformation.hh"
#include "casm/symmetry/OrbitGeneration.hh"
#include "casm/misc/cloneable_ptr.hh"

namespace CASM {

  namespace Kinetics {

    /// \brief Enumerate DiffusionTransformation for a particular IntegralCluster
    ///
    /// - InputEnumerator
    /// - Outputs all valid DiffusionTransformation for a particular
    ///   IntegralCluster. Results may include duplicates, non-canonical, unsorted forms.
    /// -
    /// - To get unique orbits, see for example make_prim_periodic_diffusion_transformation_orbits
    ///
    class DiffusionTransformationEnum : public InputEnumeratorBase<DiffusionTransformation> {

      // -- Required members -------------------

    public:

      /// \brief Construct with an IntegralCluster
      DiffusionTransformationEnum(const IntegralCluster &clust);

      std::string name() const override {
        return enumerator_name;
      }

      static const std::string enumerator_name;

    private:


      /// Implements increment
      void increment() override;


      // -- Unique -------------------

      const Structure &prim() const;

      const IntegralCluster &cluster() const;

      /// \brief The occ_counter contains the from/to occupation values for each site
      void _init_occ_counter();

      /// \brief Returns container of 'from' specie locations
      std::vector<SpecieLocation> _init_from_loc(const std::vector<Index> &occ_values);

      /// \brief Returns container of 'to' specie locations
      std::vector<SpecieLocation> _init_to_loc(const std::vector<Index> &occ_values);

      /// \brief Returns container of 'from' or 'to' specie locations
      std::vector<SpecieLocation> _init_loc(const std::vector<Index> &occ_values, Index offset);

      /// \brief Uses m_cluster, m_occ_counter, m_from_loc, and m_to_loc to set m_current
      void _set_current();
      void _update_current_occ_transform();
      void _set_current_loc();
      void _update_current_to_loc();


      Counter<std::vector<Index> > m_occ_counter;
      std::vector<SpecieLocation> m_from_loc;
      std::vector<SpecieLocation> m_to_loc;

      IntegralCluster m_cluster;
      notstd::cloneable_ptr<DiffusionTransformation> m_current;
    };

    template<typename OrbitOutputIterator, typename IntegralClusterOrbitInputIterator>
    OrbitOutputIterator make_prim_periodic_diff_trans_orbits(
      IntegralClusterOrbitInputIterator begin,
      IntegralClusterOrbitInputIterator end,
      double xtal_tol,
      OrbitOutputIterator result);
  }
}

#endif
