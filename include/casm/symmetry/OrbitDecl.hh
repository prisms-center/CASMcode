#ifndef CASM_OrbitDecl
#define CASM_OrbitDecl

namespace CASM {

  template<typename Derived> struct traits;

  template<typename Derived>
  class SymCompare;

  /// \brief Template class to be specialized for comparisons with aperiodic symmetry
  template<typename Element, typename U = void>
  class AperiodicSymCompare;

  template<typename Element, typename U = void>
  using LocalSymCompare = AperiodicSymCompare<Element, U>;

  /// \brief Template class to be specialized for comparisons with periodic symmetry
  /// of the primitive lattice
  template<typename Element, typename U = void>
  class PrimPeriodicSymCompare;

  /// \brief Template class to be specialized for comparisons with periodic symmetry
  /// of the supercell lattice
  template<typename Element, typename U = void>
  class ScelPeriodicSymCompare;

  template<typename _Element, typename _SymCompareType>
  class GenericOrbit;

  template<typename _Element, typename _SymCompareType>
  class DatabaseTypeOrbit;

  /// \brief OrbitTraits can be specialized for Orbit types that will be stored in a database
  ///
  /// Example specialized OrbitTraits:
  /// \code
  /// template<>
  /// struct OrbitTraits<DiffusionTransformation, PrimPeriodicSymCompare<DiffusionTransformation>> {
  ///   typedef DiffusionTransformation Element;
  ///   typedef PrimPeriodicSymCompare<DiffusionTransformation> SymCompareType;
  ///   typedef DatabaseTypeOrbit<Element, SymCompareType> OrbitType;
  ///
  ///   static write_pos(const OrbitType& orbit);
  ///   static std::string generate_name_impl(const OrbitType& orbit);
  /// };
  /// \endcode
  ///
  template <
    typename _Element,
    typename _SymCompareType >
  struct OrbitTraits {
    typedef _Element Element;
    typedef _SymCompareType SymCompareType;
    typedef GenericOrbit<_Element, _SymCompareType> OrbitType;
  };

  /// Alias to select between GenericOrbit and DatabaseTypeOrbit by checking
  ///   OrbitTraits<_Element, _SymCompareType>::IsDatabaseType
  template <
    typename _Element,
    typename _SymCompareType >
  using Orbit = typename OrbitTraits<_Element, _SymCompareType>::OrbitType;


  template<typename Element>
  using AperiodicOrbit = Orbit<Element, AperiodicSymCompare<Element>>;

  template<typename Element>
  using LocalOrbit = AperiodicOrbit<Element>;

  template<typename Element>
  using ScelPeriodicOrbit = Orbit<Element, ScelPeriodicSymCompare<Element>>;

  template<typename Element>
  using PrimPeriodicOrbit = Orbit<Element, PrimPeriodicSymCompare<Element>>;

}

#endif
