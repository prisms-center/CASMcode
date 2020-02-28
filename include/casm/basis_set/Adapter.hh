#ifndef BASIS_SET_ADAPTER_HH
#define BASIS_SET_ADAPTER_HH

#include "casm/basis_set/DoF.hh"
#include "casm/basis_set/DoFSet.hh"
#include "casm/crystallography/DoFSet.hh"
#include "casm/basis_set/OccupationDoFTraits.hh"

namespace CASM {
  namespace adapter {
    template <typename ToType, typename FromType>
    struct Adapter;

    /// Convert the molecular degrees of freedom in the simple vector form (stored this way in Site)
    /// to the specialized OccupantDoF type, which can be used to construct basis functions.
    //TODO: Squash OccupantDoF and DiscreteDoF into a single DiscreteDoF class that does NOT depend on Molecule (i.e. DiscreteType)
    template <typename DiscreteType>
    struct Adapter<OccupantDoF<DiscreteType>, std::vector<DiscreteType>> {
      OccupantDoF<DiscreteType> operator()(const std::vector<xtal::Molecule> &adaptable, int dof_id) {
        OccupantDoF<DiscreteType> dof(DoFType::occupation().val_traits(), "s", adaptable);
        dof.set_ID(dof_id);
        return dof;
      }
    };

    /// Convert the xtal implementation of SiteDoFSet into the CASM::DoFSet that is used to generate
    /// basis functions elsewhere in CASM
    template<>
    struct Adapter<CASM::DoFSet, xtal::SiteDoFSet> {
      CASM::DoFSet operator()(const xtal::SiteDoFSet &adaptable, int dof_id) {
        CASM::DoFSet adapted_dofset(adaptable.traits(), adaptable.excluded_occs());
        adapted_dofset.set_ID(dof_id);
        return adapted_dofset;
      }
    };
  }
}



#endif
