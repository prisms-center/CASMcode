#ifndef CASM_DiffusionTransformationTraits
#define CASM_DiffusionTransformationTraits

#include "casm/symmetry/OrbitDecl.hh"
#include <iostream>

namespace CASM {
  namespace xtal {
    class UnitCellCoord;
  }

  namespace sym {
    class CopyApplyDefault_f;
  }

  namespace Kinetics {
    class DiffusionTransformation;
    class DiffTransInvariants;
    typedef PrimPeriodicSymCompare<Kinetics::DiffusionTransformation> PrimPeriodicDiffTransSymCompare;
    typedef ScelPeriodicSymCompare<Kinetics::DiffusionTransformation> ScelPeriodicDiffTransSymCompare;
  }

  typedef Kinetics::PrimPeriodicDiffTransSymCompare PrimPeriodicDiffTransSymCompare;
  typedef Kinetics::ScelPeriodicDiffTransSymCompare ScelPeriodicDiffTransSymCompare;

  template <typename Base>
  class CopyApplyDefault;

  /// Traits necessary for SymCompare
  template<>
  struct traits<Kinetics::DiffusionTransformation> {

    typedef typename Kinetics::DiffusionTransformation MostDerived;
    typedef typename Kinetics::DiffTransInvariants InvariantsType;
    static xtal::UnitCellCoord position(const Kinetics::DiffusionTransformation &diff_trans);
    static std::string name;
    template<typename Base>
    using CopyApplyType = CopyApplyDefault<Base>;
    typedef sym::CopyApplyDefault_f copy_apply_f_type;
  };

  /// Specialization gives required for DatabaseTypeOrbit
  template<>
  struct OrbitTraits<PrimPeriodicDiffTransSymCompare> {
    typedef Kinetics::DiffusionTransformation Element;
    typedef PrimPeriodicDiffTransSymCompare SymCompareType;
    typedef DatabaseTypeOrbit<SymCompareType> OrbitType;

    static void write_pos(const OrbitType &orbit, std::ostream &sout);
    static std::string generate_name_impl(const OrbitType &orbit);
  };

  typedef OrbitTraits<PrimPeriodicDiffTransSymCompare> PrimPeriodicDiffTransOrbitTraits;

  namespace Kinetics {
    typedef PrimPeriodicOrbit<Kinetics::DiffusionTransformation> PrimPeriodicDiffTransOrbit;
    typedef ScelPeriodicOrbit<Kinetics::DiffusionTransformation> ScelPeriodicDiffTransOrbit;
  }
  typedef Kinetics::PrimPeriodicDiffTransOrbit PrimPeriodicDiffTransOrbit;
  typedef Kinetics::ScelPeriodicDiffTransOrbit ScelPeriodicDiffTransOrbit;

  /// Traits necessary for database types
  template<>
  struct traits<PrimPeriodicDiffTransOrbit> {
    static const std::string name;
    static const std::string short_name;
    static const std::string orbit_type_name;
    static bool name_compare(std::string A, std::string B);
  };
}

#endif
