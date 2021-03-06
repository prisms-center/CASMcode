#ifndef CASM_symmetry_SimpleOrbit
#define CASM_symmetry_SimpleOrbit

#include "casm/container/multivector.hh"

namespace CASM {

class SymOp;
class SymGroup;

/// Generate an orbit of unique Element generated by a group
///
/// SimpleOrbit is similar to Orbit, but does not require use of SymOp and
/// SymGroup and does not use an existing multiplication table. This allows it
/// to be used with large groups (for example, groups that include translations
/// within a supercell). Therefore it also does not maintain the generating
/// group or equivalence map and element sorting is not based on the
/// equivalence map.
///
/// Use of SimpleOrbit requires implementation of a SymCompareType that
/// describes how to apply group elements to the orbit elements and how to
/// compare orbit elements. The struct `VectorSymCompare` in `symmetry/
/// Symmetrizer.cc` provides an example implementation.
///
/// \ingroup OrbitGeneration
///
template <typename _SymCompareType>
class SimpleOrbit {
 public:
  using size_type = unsigned int;
  using SymCompareType = _SymCompareType;
  using Element = typename _SymCompareType::Element;
  using InvariantsType = typename _SymCompareType::InvariantsType;
  using const_iterator = typename std::vector<Element>::const_iterator;

  /// Construct a SimpleOrbit using a range of SymOp or PermuteIterator
  template <typename GroupIterator>
  SimpleOrbit(typename SymCompareType::Element const &_generating_element,
              GroupIterator _group_begin, GroupIterator _group_end,
              SymCompareType const &_sym_compare);

  const_iterator begin() const { return m_element.cbegin(); }

  const_iterator end() const { return m_element.cend(); }

  const_iterator cbegin() const { return m_element.cbegin(); }

  const_iterator cend() const { return m_element.cend(); }

  size_type size() const { return m_element.size(); }

  /// \brief Identical to element(0)
  Element const &prototype() const { return m_element[0]; }

  /// \brief Return Element at index, without bounds checking
  ///
  /// - May not be prepared
  Element const &operator[](size_type index) const { return element(index); }

  /// \brief Equivalent to operator[](size_type index) const
  ///
  /// - May not be prepared
  Element const &element(size_type index) const { return m_element[index]; }

  /// \brief const Access vector of Element
  ///
  /// - May not be prepared
  std::vector<Element> const &elements() const { return m_element; }

  /// \brief Not implemented: Will throw. Exists only for compatibility with
  /// orbit printing.
  const multivector<SymOp>::X<2> &equivalence_map() const {
    throw std::runtime_error("No equivalence map for SimpleOrbit");
  }

  /// \brief Not implemented: Will throw. Exists only for compatibility with
  /// orbit printing.
  const SymGroup &generating_group() const {
    throw std::runtime_error("No generating group for SimpleOrbit");
  }

  /// \brief Return the SymCompare functor reference
  ///
  /// - implements symmetry properties of this orbit
  SymCompareType const &sym_compare() const { return m_sym_compare; }

  InvariantsType const &invariants() const { return m_invariants; }

  /// Compare orbits, using SymCompareType::inter_orbit_compare
  bool operator<(const SimpleOrbit &B) const;

 private:
  SymCompareType m_sym_compare;
  InvariantsType m_invariants;
  std::vector<Element> m_element;
};
}  // namespace CASM

#endif
