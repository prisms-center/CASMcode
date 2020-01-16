#ifndef CASM_OccupationTransformationTraits
#define CASM_OccupationTransformationTraits

#include "casm/global/definitions.hh"
#include "casm/symmetry/OrbitDecl.hh"

namespace CASM {

  template <typename Base>
  class CopyApplyWithPrim_crtp;

  // --- OccupationTransformation ---
  //
  // should rename:
  //   OccupationTransformation -> OccTransElement
  //   OccPerturbation -> OccupationTransformation

  namespace Kinetics {
    class OccupationTransformation;
  }
  typedef Kinetics::OccupationTransformation OccupationTransformation;

  namespace xtal {
    class UnitCellCoord;
  }

  class OccPerturbation;
  class OccPerturbationInvariants;

  /// Traits necessary for SymCompare
  template<>
  struct traits<OccupationTransformation> {
    template <typename Base>
    using copy_apply_crtp_type = CopyApplyWithPrim_crtp<Base>;
  };

  template<>
  struct traits<OccPerturbation> {
    typedef OccupationTransformation Element;
    typedef OccPerturbationInvariants InvariantsType;
    typedef Index size_type;
    static xtal::UnitCellCoord position(const OccPerturbation &perturb);
    template <typename Base>
    using copy_apply_crtp_type = CopyApplyWithPrim_crtp<Base>;
  };

  typedef PrimPeriodicSymCompare<OccPerturbation> PrimPeriodicOccPerturbSymCompare;
  typedef ScelPeriodicSymCompare<OccPerturbation> ScelPeriodicOccPerturbSymCompare;

  typedef PrimPeriodicOrbit<OccPerturbation> PrimPeriodicOccPerturbOrbit;
  typedef ScelPeriodicOrbit<OccPerturbation> ScelPeriodicOccPerturbOrbit;
}

#endif
