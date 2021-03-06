#ifndef CASM_OrbitGeneration_impl
#define CASM_OrbitGeneration_impl

#include "casm/symmetry/OrbitGeneration.hh"
#include "casm/symmetry/PermuteIterator.hh"

namespace CASM {

template <typename _OrbitType>
OrbitGenerators<_OrbitType>::OrbitGenerators(const SymGroup &_group,
                                             const SymCompareType &_sym_compare)
    : group(_group),
      sym_compare(_sym_compare),
      m_element_compare(sym_compare),
      m_generate_canonical(group, sym_compare),
      elements(m_element_compare) {}

/// \brief Try inserting an element, after generating the canonical form
template <typename _OrbitType>
std::pair<typename OrbitGeneratorSet<_OrbitType>::iterator, bool>
OrbitGenerators<_OrbitType>::insert(const Element &test) {
  return elements.insert(m_generate_canonical(test));
}

/// \brief Try inserting an element, assuming it is in canonical form
template <typename _OrbitType>
std::pair<typename OrbitGeneratorSet<_OrbitType>::iterator, bool>
OrbitGenerators<_OrbitType>::insert_canonical(const Element &test) {
  return elements.insert(test);
}

/// \brief Construct Orbit from all generating elements
template <typename _OrbitType>
template <typename OrbitOutputIterator>
OrbitOutputIterator OrbitGenerators<_OrbitType>::make_orbits(
    OrbitOutputIterator result) {
  for (const auto &e : elements) {
    *result++ = OrbitType(e, group, sym_compare);
  }
  return result;
}

// /// \brief Construct Orbit from all generating elements, including PrimClex
// /// pointer
// template <typename _OrbitType>
// template <typename OrbitOutputIterator>
// OrbitOutputIterator OrbitGenerators<_OrbitType>::make_orbits(
//     OrbitOutputIterator result, const PrimClex &primclex) {
//   for (const auto &e : elements) {
//     *result++ = OrbitType(e, group, sym_compare, &primclex);
//   }
//   return result;
// }

template <typename SymCompareType>
std::vector<Orbit<SymCompareType>> generate_orbits(
    std::vector<typename SymCompareType::Element> const &generating_elements,
    SymGroup const &generating_group, SymCompareType const &sym_compare) {
  typedef Orbit<SymCompareType> OrbitType;
  OrbitGenerators<OrbitType> generators(generating_group, sym_compare);
  for (auto const &el : generating_elements) {
    generators.insert(el);
  }
  std::vector<Orbit<SymCompareType>> orbits;
  orbits.reserve(generating_elements.size());
  generators.make_orbits(std::back_inserter(orbits));
  return orbits;
}

// --- template<typename _OrbitType> class OrbitGeneratorCompare ---

template <typename _OrbitType>
OrbitGeneratorCompare<_OrbitType>::OrbitGeneratorCompare(
    const SymCompareType &_sym_compare)
    : sym_compare(_sym_compare) {}

template <typename _OrbitType>
bool OrbitGeneratorCompare<_OrbitType>::operator()(const Element &A,
                                                   const Element &B) const {
  auto A_invariants = sym_compare.make_invariants(A);
  auto B_invariants = sym_compare.make_invariants(B);
  return sym_compare.inter_orbit_compare(A, A_invariants, B, B_invariants);
}

// --- template<typename _OrbitType> class CanonicalGenerator ---

template <typename _OrbitType>
CanonicalGenerator<_OrbitType>::CanonicalGenerator(
    const std::vector<SymOp> &_generating_group,
    const SymCompareType &_sym_compare)
    : generating_group(_generating_group),
      sym_compare(_sym_compare),
      m_to_canonical(nullptr) {}

/// \brief Applies symmetry to return an equivalent Element in a canonical form
template <typename _OrbitType>
typename CanonicalGenerator<_OrbitType>::Element
CanonicalGenerator<_OrbitType>::operator()(const Element &e) const {
  Element result = sym_compare.prepare(e);
  for (const auto &op : generating_group) {
    auto test = sym_compare.prepare(sym_compare.copy_apply(op, e));
    if (sym_compare.compare(result, test)) {
      result = test;
      m_to_canonical = &op;
    }
  }
  return result;
}

/// \brief After using call operator, this can be checked
template <typename _OrbitType>
const SymOp &CanonicalGenerator<_OrbitType>::to_canonical() const {
  return *m_to_canonical;
}

/// \brief After using call operator, this can be checked
template <typename _OrbitType>
const SymOp &CanonicalGenerator<_OrbitType>::from_canonical() const {
  return to_canonical().inverse();
}

// // --- template<typename _SymCompareType> class PermuteCanonicalGenerator ---
//
// template<typename _SymCompareType>
// PermuteCanonicalGenerator<_SymCompareType>::PermuteCanonicalGenerator(
//   PermuteIterator _permute_begin,
//   PermuteIterator _permute_end,
//   SymCompareType const &_sym_compare) :
//   permute_begin(_permute_begin),
//   permute_end(_permute_end),
//   sym_compare(_sym_compare),
//   m_to_canonical(nullptr) {}
//
// /// \brief Applies symmetry to return an equivalent Element in a canonical
// form template<typename _SymCompareType> typename
// PermuteCanonicalGenerator<_SymCompareType>::Element
// PermuteCanonicalGenerator<_SymCompareType>::operator()(Element const &e)
// const {
//   Element result = sym_compare.prepare(e);
//   for(auto permute_it = permute_begin; permute_it != permute_end;
//   ++permute_it) {
//     auto test =
//     sym_compare.prepare(sym_compare.copy_apply(permute_it.sym_op(), e));
//     if(sym_compare.compare(result, test)) {
//       result = test;
//       m_to_canonical = &permute_it;
//     }
//   }
//   return result;
// }
//
// /// \brief After using call operator, this can be checked
// template<typename _SymCompareType>
// PermuteIterator const
// &PermuteCanonicalGenerator<_SymCompareType>::to_canonical() const {
//   return *m_to_canonical;
// }
//
// /// \brief After using call operator, this can be checked
// template<typename _SymCompareType>
// PermuteIterator const
// &PermuteCanonicalGenerator<_SymCompareType>::from_canonical() const {
//   return to_canonical().inverse();
// }

// --- template<typename _OrbitType> class IsCanonical<_OrbitType> ---

template <typename _OrbitType>
IsCanonical<_OrbitType>::IsCanonical(
    const std::vector<SymOp> &_generating_group,
    const SymCompareType &_sym_compare)
    : generating_group(_generating_group), sym_compare(_sym_compare) {}

/// \brief Applies symmetry to check if any Element is greater than e
template <typename _OrbitType>
bool IsCanonical<_OrbitType>::operator()(const Element &e) const {
  return std::none_of(
      generating_group.begin(), generating_group.end(), [&](const SymOp &op) {
        auto test = sym_compare.prepare(sym_compare.copy_apply(op, e));
        return sym_compare.compare(e, test);
      });
}

}  // namespace CASM

#endif
