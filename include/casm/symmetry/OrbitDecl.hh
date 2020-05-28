#ifndef CASM_OrbitDecl
#define CASM_OrbitDecl

//#include "casm/symmetry/SymCompare.hh"
namespace CASM {

  template<typename Derived> struct traits;

  template<typename Derived>
  class SymCompare;

  // /// \brief Template class to be specialized for comparisons with aperiodic symmetry
  // template<typename Element, typename U = void>
  // class AperiodicSymCompare;
  //
  // template<typename Element, typename U = void>
  // using LocalSymCompare = AperiodicSymCompare<Element, U>;
  //
  // /// \brief Template class to be specialized for comparisons with periodic symmetry
  // /// of the primitive lattice
  // template<typename Element, typename U = void>
  // class PrimPeriodicSymCompare;
  //
  // /// \brief Template class to be specialized for comparisons with periodic symmetry
  // /// of the supercell lattice
  // template<typename Element, typename U = void>
  // class ScelPeriodicSymCompare;
  //
  // /// \brief Template class to be specialized for comparisons with periodic symmetry
  // /// of the supercell lattice, for clusters that do not extend outside the supercell
  // template<typename Element, typename U = void>
  // class WithinScelSymCompare;

  /// \brief Template class to be specialized for comparisons of vector directions
  template<typename Element, typename SymApply>
  class DirectionSymCompare;

  /// \brief Template class to be specialized for comparisons of vector subspaces
  template<typename Element, typename SymApply>
  class SubspaceSymCompare;


  template<typename _SymCompareType>
  class Orbit;

  // template<typename _SymCompareType>
  // class GenericOrbit;
  //
  // template<typename _SymCompareType>
  // class DatabaseTypeOrbit;

  // TODO: remove
  // // -- Previously from SymCompare.hh --
  // /// \brief Traits class for AperiodicSymCompare
  // template<typename _Element>
  // struct traits<AperiodicSymCompare<_Element>> {
  //   typedef _Element Element;
  //   typedef AperiodicSymCompare<Element> MostDerived;
  // };
  //
  // /// \brief Traits class for PrimPeriodicSymCompare
  // template<typename _Element>
  // struct traits<PrimPeriodicSymCompare<_Element>> {
  //   typedef _Element Element;
  //   typedef PrimPeriodicSymCompare<Element> MostDerived;
  // };
  //
  // /// \brief Traits class for ScelPeriodicSymCompare
  // template<typename _Element>
  // struct traits<ScelPeriodicSymCompare<_Element>> {
  //   typedef _Element Element;
  //   typedef ScelPeriodicSymCompare<Element> MostDerived;
  // };
  //
  // /// \brief Traits class for WithinScelPeriodicSymCompare
  // template<typename _Element>
  // struct traits<WithinScelSymCompare<_Element>> {
  //   typedef _Element Element;
  //   typedef WithinScelSymCompare<Element> MostDerived;
  // };

  /// \brief Traits class for DirectionSymCompare
  template<typename _Element, typename _SymApply>
  struct traits<DirectionSymCompare<_Element, _SymApply>> {
    typedef _Element Element;
    typedef _SymApply SymApply;
    typedef DirectionSymCompare<_Element, _SymApply> MostDerived;
  };

  /// \brief Traits class for SubspaceSymCompare
  template<typename _Element, typename _SymApply>
  struct traits<SubspaceSymCompare<_Element, _SymApply>> {
    typedef _Element Element;
    typedef _SymApply SymApply;
    typedef SubspaceSymCompare<_Element, _SymApply> MostDerived;
  };
  //\-- from SymCompare.hh --

  // /// \brief OrbitTraits can be specialized for Orbit types that will be stored in a database
  // ///
  // /// Example specialized OrbitTraits:
  // /// \code
  // /// template<>
  // /// struct OrbitTraits<DiffusionTransformation, PrimPeriodicSymCompare<DiffusionTransformation>> {
  // ///   typedef DiffusionTransformation Element;
  // ///   typedef PrimPeriodicSymCompare<DiffusionTransformation> SymCompareType;
  // ///   typedef DatabaseTypeOrbit<SymCompareType> OrbitType;
  // ///
  // ///   static write_pos(const OrbitType& orbit);
  // ///   static std::string generate_name_impl(const OrbitType& orbit);
  // /// };
  // /// \endcode
  // ///
  // template <typename _SymCompareType>
  // struct OrbitTraits {
  //   using Element = typename traits<_SymCompareType>::Element;
  //   using SymCompareType = _SymCompareType;
  //   using OrbitType = GenericOrbit<_SymCompareType>;
  // };

  /// Alias to select between GenericOrbit and DatabaseTypeOrbit by checking
  ///   OrbitTraits<_SymCompareType>::IsDatabaseType
  // template <typename _SymCompareType>
  // using Orbit = typename traits<_SymCompareType>::OrbitType;


  // template<typename Element>
  // using AperiodicOrbit = Orbit<AperiodicSymCompare<Element>>;
  //
  // template<typename Element>
  // using LocalOrbit = AperiodicOrbit<Element>;
  //
  // template<typename Element>
  // using ScelPeriodicOrbit = Orbit<ScelPeriodicSymCompare<Element>>;
  //
  // template<typename Element>
  // using PrimPeriodicOrbit = Orbit<PrimPeriodicSymCompare<Element>>;
  //
  // template<typename Element>
  // using WithinScelOrbit = Orbit<WithinScelSymCompare<Element>>;

}

#endif
