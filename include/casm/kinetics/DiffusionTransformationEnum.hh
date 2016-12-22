#ifndef CASM_DiffusionTransformationEnum
#define CASM_DiffusionTransformationEnum

#include "casm/container/Counter.hh"
#include "casm/container/InputEnumerator.hh"
#include "casm/symmetry/SymCompare.hh"
#include "casm/kinetics/DoFTransformation.hh"
#include "casm/misc/cloneable_ptr.hh"

namespace CASM {

  namespace Kinetics {

    /// \brief Enumerate DiffusionTransformation for a particular IntegralCluster
    ///
    /// - InputEnumerator
    /// - Outputs all unique, canonical DiffusionTransformation for a particular
    ///   IntegralCluster and generating group
    class DiffusionTransformationEnum : public InputEnumeratorBase<DiffusionTransformation> {

      // -- Required members -------------------

    public:

      /// \brief Construct with an IntegralCluster, prim factor group, and SymCompare object
      template<typename SymCompareType>
      DiffusionTransformationEnum(const IntegralCluster &clust, const SymCompareType &sym_compare) :
        DiffusionTransformationEnum(clust, invariant_subgroup(clust, clust.prim().factor_group(), sym_compare)) {}

      /// \brief Construct with an IntegralCluster, and generating_grp
      DiffusionTransformationEnum(const IntegralCluster &clust, const SymGroup &generating_grp);

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
      const SymGroup *m_generating_grp;
      notstd::cloneable_ptr<DiffusionTransformation> m_current;
    };
  }
}

#endif
