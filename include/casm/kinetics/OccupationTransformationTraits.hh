#ifndef CASM_OccupationTransformationTraits
#define CASM_OccupationTransformationTraits

#include "casm/symmetry/OrbitDecl.hh"

namespace CASM {

  // --- OccupationTransformation ---
  //
  // should rename:
  //   OccupationTransformation -> OccTransElement
  //   OccPerturbation -> OccupationTransformation

  namespace Kinetics {
    class OccupationTransformation;
  }
  typedef Kinetics::OccupationTransformation OccupationTransformation;

  class OccPerturbation;
  class OccPerturbationInvariants;
  class UnitCellCoord;

  /// Traits necessary for SymCompare
  template<>
  struct traits<OccPerturbation> {
    typedef OccupationTransformation Element;
    typedef OccPerturbationInvariants InvariantsType;
    typedef Index size_type;
    static UnitCellCoord position(const OccPerturbation &perturb);
  };

  typedef PrimPeriodicSymCompare<OccPerturbation> PrimPeriodicOccPerturbSymCompare;
  typedef ScelPeriodicSymCompare<OccPerturbation> ScelPeriodicOccPerturbSymCompare;

  typedef PrimPeriodicOrbit<OccPerturbation> PrimPeriodicOccPerturbOrbit;
  typedef ScelPeriodicOrbit<OccPerturbation> ScelPeriodicOccPerturbOrbit;



}

#endif
