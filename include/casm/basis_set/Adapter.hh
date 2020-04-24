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
    /// The int parameter dof_id is used to set the ID of the OccupantDoF, and is most likely
    /// the index of the relevant Site within the Structure basis.
    //TODO: Squash OccupantDoF and DiscreteDoF into a single DiscreteDoF class that does NOT depend on Molecule (i.e. DiscreteType)
    template <typename DiscreteType>
    struct Adapter<OccupantDoF<DiscreteType>, std::vector<DiscreteType>> {
      OccupantDoF<DiscreteType> operator()(const std::vector<xtal::Molecule> &adaptable, int dof_id) {
        OccupantDoF<DiscreteType> dof(DoFType::occupation().val_traits(), "s", adaptable);
        dof.set_ID(dof_id);
        return dof;
      }
    };

    /// Convert the xtal implementation of DoFSet into the CASM::DoFSet that is used to generate
    /// basis functions elsewhere in CASM.
    /// IDs of each component are set to values 0-n-1, and are then locked to avoid future changes to them.
    template<>
    struct Adapter<CASM::DoFSet, xtal::DoFSet> {
      CASM::DoFSet operator()(const xtal::DoFSet &adaptable) {
        CASM::DoFSet adapted_dofset(adaptable.traits(),
                                    adaptable.component_names(),
                                    adaptable.basis());
        adapted_dofset.set_sequential_IDs();
        adapted_dofset.lock_IDs();
        return adapted_dofset;
      }
    };

    /// Convert the xtal implementation of SiteDoFSet into the CASM::DoFSet that is used to generate
    /// basis functions elsewhere in CASM.
    /// The int parameter dof_id is used to set the ID of the DoFSet, and is most likely
    /// the index of the relevant Site within the Structure basis.
    template<>
    struct Adapter<CASM::DoFSet, xtal::SiteDoFSet> {
      CASM::DoFSet operator()(const xtal::SiteDoFSet &adaptable, int dof_id) {
        CASM::DoFSet adapted_dofset(adaptable.traits(),
                                    adaptable.component_names(),
                                    adaptable.basis(),
                                    adaptable.excluded_occupants());
        adapted_dofset.set_ID(dof_id);
        return adapted_dofset;
      }
    };
  }
}



#endif
