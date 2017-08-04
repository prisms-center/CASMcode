#ifndef CASM_DiffusionTransformationTraits
#define CASM_DiffusionTransformationTraits

#include "casm/symmetry/OrbitDecl.hh"
#include <iostream>

namespace CASM {
  class UnitCellCoord;

  namespace Kinetics {
    class DiffusionTransformation;
    class DiffTransInvariants;
  }
  typedef PrimPeriodicSymCompare<Kinetics::DiffusionTransformation> PrimPeriodicDiffTransSymCompare;
  typedef ScelPeriodicSymCompare<Kinetics::DiffusionTransformation> ScelPeriodicDiffTransSymCompare;

  /// Traits necessary for SymCompare
  template<>
  struct traits<Kinetics::DiffusionTransformation> {

    typedef typename Kinetics::DiffusionTransformation MostDerived;
    typedef typename Kinetics::DiffTransInvariants InvariantsType;
    static UnitCellCoord position(const Kinetics::DiffusionTransformation &diff_trans);
  };

  /// Specialization gives required functions for storing PrimPeriodicDiffTransOrbit in a Database
  template<>
  struct OrbitTraits<Kinetics::DiffusionTransformation, PrimPeriodicDiffTransSymCompare> {
    typedef Kinetics::DiffusionTransformation Element;
    typedef PrimPeriodicDiffTransSymCompare SymCompareType;
    typedef DatabaseTypeOrbit<Element, SymCompareType> OrbitType;

    static void write_pos(const OrbitType &orbit, std::ostream &sout);
    static std::string generate_name_impl(const OrbitType &orbit);
  };

  typedef OrbitTraits<Kinetics::DiffusionTransformation, PrimPeriodicDiffTransSymCompare> PrimPeriodicDiffTransOrbitTraits;

  typedef PrimPeriodicOrbit<Kinetics::DiffusionTransformation> PrimPeriodicDiffTransOrbit;
  typedef ScelPeriodicOrbit<Kinetics::DiffusionTransformation> ScelPeriodicDiffTransOrbit;


}

#endif
