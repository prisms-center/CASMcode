#ifndef CASM_HasCanonicalForm_impl
#define CASM_HasCanonicalForm_impl

#include "casm/clex/HasCanonicalForm.hh"
#include "casm/clex/Supercell.hh"
#include "casm/crystallography/CanonicalForm.hh"
#include "casm/crystallography/SymTools.hh"
#include "casm/symmetry/InvariantSubgroup_impl.hh"
#include "casm/symmetry/OrbitGeneration.hh"
#include "casm/symmetry/PermuteIterator.hh"
#include "casm/symmetry/ScelOrbitGeneration.hh"
#include "casm/symmetry/SymGroup.hh"
#include "casm/symmetry/SymOp.hh"

namespace CASM {

// --- template<typename _Base> class CanonicalForm ---

template <typename _Base>
template <typename SymCompareType>
bool CanonicalForm<_Base>::is_canonical(
    const std::vector<SymOp> &g, const SymCompareType &sym_compare) const {
  IsCanonical<Orbit<SymCompareType>> f(g, sym_compare);
  return f(derived());
}

template <typename _Base>
template <typename SymCompareType>
typename CanonicalForm<_Base>::MostDerived CanonicalForm<_Base>::canonical_form(
    const std::vector<SymOp> &g, const SymCompareType &sym_compare) const {
  CanonicalGenerator<Orbit<SymCompareType>> f(g, sym_compare);
  return f(derived());
}

template <typename _Base>
template <typename SymCompareType>
bool CanonicalForm<_Base>::is_sym_equivalent(
    const MostDerived &other, const std::vector<SymOp> &g,
    const SymCompareType &sym_compare) const {
  CanonicalGenerator<Orbit<SymCompareType>> f(g, sym_compare);
  return sym_compare.equal(f(derived()), f(other));
}

template <typename _Base>
template <typename ObjIterator, typename SymCompareType>
ObjIterator CanonicalForm<_Base>::find_sym_equivalent(
    ObjIterator begin, ObjIterator end, const std::vector<SymOp> &g,
    const SymCompareType &sym_compare) const {
  CanonicalGenerator<Orbit<SymCompareType>> f(g, sym_compare);
  auto canon = f(derived());
  auto is_sym_equiv = [&](const MostDerived &test) {
    return sym_compare.equal(canon, f(test));
  };
  return std::find_if(begin, end, is_sym_equiv);
}

template <typename _Base>
template <typename SymCompareType>
SymOp CanonicalForm<_Base>::to_canonical(
    const std::vector<SymOp> &g, const SymCompareType &sym_compare) const {
  CanonicalGenerator<Orbit<SymCompareType>> f(g, sym_compare);
  f(derived());
  return f.to_canonical();
}

template <typename _Base>
template <typename SymCompareType>
SymOp CanonicalForm<_Base>::from_canonical(
    const std::vector<SymOp> &g, const SymCompareType &sym_compare) const {
  return to_canonical(g, sym_compare).inverse();
}

template <typename _Base>
template <typename SymCompareType>
SymGroup CanonicalForm<_Base>::invariant_subgroup(
    const SymGroup &super_grp, const SymCompareType &sym_compare) const {
  return make_invariant_subgroup(derived(), super_grp, sym_compare);
}

template <typename _Base>
bool CanonicalForm<_Base>::is_canonical(const Supercell &scel) const {
  return is_canonical(scel, scel.sym_info().permute_begin(),
                      scel.sym_info().permute_end());
}

template <typename _Base>
typename CanonicalForm<_Base>::MostDerived CanonicalForm<_Base>::canonical_form(
    const Supercell &scel) const {
  return canonical_form(scel, scel.sym_info().permute_begin(),
                        scel.sym_info().permute_end());
}

/// True if this and B have same canonical form
template <typename _Base>
bool CanonicalForm<_Base>::is_sym_equivalent(const MostDerived &B,
                                             const Supercell &scel) const {
  return is_sym_equivalent(scel, scel.sym_info().permute_begin(),
                           scel.sym_info().permute_end());
}

template <typename _Base>
SymOp CanonicalForm<_Base>::to_canonical(const Supercell &scel) const {
  return to_canonical(scel, scel.sym_info().permute_begin(),
                      scel.sym_info().permute_end());
}

template <typename _Base>
SymOp CanonicalForm<_Base>::from_canonical(const Supercell &scel) const {
  return from_canonical(scel, scel.sym_info().permute_begin(),
                        scel.sym_info().permute_end());
}

template <typename _Base>
std::vector<PermuteIterator> CanonicalForm<_Base>::invariant_subgroup(
    const Supercell &scel) const {
  return invariant_subgroup(scel, scel.sym_info().permute_begin(),
                            scel.sym_info().permute_end());
}

template <typename _Base>
template <typename PermuteIteratorIt>
bool CanonicalForm<_Base>::is_canonical(const Supercell &scel,
                                        PermuteIteratorIt begin,
                                        PermuteIteratorIt end) const {
  ScelIsCanonical<MostDerived> f(scel);
  return f(derived(), begin, end);
}

template <typename _Base>
template <typename PermuteIteratorIt>
typename CanonicalForm<_Base>::MostDerived CanonicalForm<_Base>::canonical_form(
    const Supercell &scel, PermuteIteratorIt begin,
    PermuteIteratorIt end) const {
  ScelCanonicalGenerator<MostDerived> f(scel);
  return f(derived());
}

/// True if this and B have same canonical form
template <typename _Base>
template <typename PermuteIteratorIt>
bool CanonicalForm<_Base>::is_sym_equivalent(const MostDerived &B,
                                             const Supercell &scel,
                                             PermuteIteratorIt begin,
                                             PermuteIteratorIt end) const {
  ScelCanonicalGenerator<MostDerived> f(scel);
  return f.sym_compare.equal(f(derived(), begin, end), f(B, begin, end));
}

/// Find element that has the same canonical form, with respect to a subgroup
template <typename _Base>
template <typename ObjIterator, typename PermuteIteratorIt>
ObjIterator CanonicalForm<_Base>::find_sym_equivalent(
    ObjIterator obj_begin, ObjIterator obj_end, const Supercell &scel,
    PermuteIteratorIt begin, PermuteIteratorIt end) const {
  ScelCanonicalGenerator<MostDerived> f(scel);
  auto canon = f(derived(), begin, end);
  auto is_sym_equiv = [&](const MostDerived &test) {
    return f.sym_compare.equal(canon, f(test, begin, end));
  };
  return std::find_if(obj_begin, obj_end, is_sym_equiv);
}

template <typename _Base>
template <typename PermuteIteratorIt>
SymOp CanonicalForm<_Base>::to_canonical(const Supercell &scel,
                                         PermuteIteratorIt begin,
                                         PermuteIteratorIt end) const {
  ScelCanonicalGenerator<MostDerived> f(scel);
  f(derived());
  return f.to_canonical();
}

template <typename _Base>
template <typename PermuteIteratorIt>
SymOp CanonicalForm<_Base>::from_canonical(const Supercell &scel,
                                           PermuteIteratorIt begin,
                                           PermuteIteratorIt end) const {
  // Please use -Werror=return-type in your compiler flags
  return to_canonical(scel, begin, end).inverse();
}

template <typename _Base>
template <typename PermuteIteratorIt>
std::vector<PermuteIterator> CanonicalForm<_Base>::invariant_subgroup(
    const Supercell &scel, PermuteIteratorIt begin,
    PermuteIteratorIt end) const {
  return make_invariant_subgroup(derived(), scel, begin, end);
}

// --- template<typename Base> class ConfigCanonicalForm<Base>

template <typename Base>
bool ConfigCanonicalForm<Base>::is_sym_equivalent(const MostDerived &B) const {
  return this->canonical_form() == B.canonical_form();
}

template <typename Base>
template <typename ConfigIterator>
ConfigIterator ConfigCanonicalForm<Base>::find_sym_equivalent(
    const MostDerived &B, ConfigIterator obj_begin,
    ConfigIterator obj_end) const {
  auto canon = this->canonical_form();
  auto is_sym_equiv = [&](const MostDerived &test) {
    return canon == test.canonical_form();
  };
  return std::find_if(obj_begin, obj_end, is_sym_equiv);
}

template <typename Base>
bool ConfigCanonicalForm<Base>::is_canonical() const {
  return is_canonical(derived().supercell().sym_info().permute_begin(),
                      derived().supercell().sym_info().permute_end());
}

template <typename Base>
typename ConfigCanonicalForm<Base>::MostDerived
ConfigCanonicalForm<Base>::canonical_form() const {
  return canonical_form(derived().supercell().sym_info().permute_begin(),
                        derived().supercell().sym_info().permute_end());
}

template <typename Base>
PermuteIterator ConfigCanonicalForm<Base>::to_canonical() const {
  return to_canonical(derived().supercell().sym_info().permute_begin(),
                      derived().supercell().sym_info().permute_end());
}

template <typename Base>
PermuteIterator ConfigCanonicalForm<Base>::from_canonical() const {
  return from_canonical(derived().supercell().sym_info().permute_begin(),
                        derived().supercell().sym_info().permute_end());
}

template <typename Base>
std::vector<PermuteIterator> ConfigCanonicalForm<Base>::invariant_subgroup()
    const {
  return invariant_subgroup(derived().supercell().sym_info().permute_begin(),
                            derived().supercell().sym_info().permute_end());
}

template <typename Base>
template <typename PermuteIteratorIt>
bool ConfigCanonicalForm<Base>::is_canonical(PermuteIteratorIt begin,
                                             PermuteIteratorIt end) const {
  return std::none_of(begin, end, derived().less());
}

/// True if this and B have same canonical form
template <typename Base>
template <typename PermuteIteratorIt>
bool ConfigCanonicalForm<Base>::is_sym_equivalent(const MostDerived &B,
                                                  PermuteIteratorIt begin,
                                                  PermuteIteratorIt end) const {
  return this->canonical_form(begin, end) == B.canonical_form(begin, end);
}

template <typename Base>
template <typename ConfigIterator, typename PermuteIteratorIt>
ConfigIterator ConfigCanonicalForm<Base>::find_sym_equivalent(
    ConfigIterator obj_begin, ConfigIterator obj_end, PermuteIteratorIt begin,
    PermuteIteratorIt end) const {
  auto canon = this->canonical_form(begin, end);
  auto is_sym_equiv = [&](const MostDerived &test) {
    return canon == test.canonical_form(begin, end);
  };
  return std::find_if(obj_begin, obj_end, is_sym_equiv);
}

template <typename Base>
template <typename PermuteIteratorIt>
typename ConfigCanonicalForm<Base>::MostDerived
ConfigCanonicalForm<Base>::canonical_form(PermuteIteratorIt begin,
                                          PermuteIteratorIt end) const {
  return copy_apply(to_canonical(begin, end), derived());
}

template <typename Base>
template <typename PermuteIteratorIt>
PermuteIterator ConfigCanonicalForm<Base>::to_canonical(
    PermuteIteratorIt begin, PermuteIteratorIt end) const {
  return *std::max_element(begin, end, derived().less());
}

template <typename Base>
template <typename PermuteIteratorIt>
PermuteIterator ConfigCanonicalForm<Base>::from_canonical(
    PermuteIteratorIt begin, PermuteIteratorIt end) const {
  // simplest version: use the inverse of the first element that results in the
  // canonical form return to_canonical(begin, end).inverse();

  // alternate version: the lowest index element that transforms canonical form
  // to this
  auto less = derived().less();
  auto _to_canonical = begin;
  auto _from_canonical = _to_canonical.inverse();
  for (auto it = begin; it < end; ++it) {
    if (less(_to_canonical, it)) {
      _to_canonical = it;
      _from_canonical = _to_canonical.inverse();
    }
    // other permutations that result in canonical config may have a lower index
    // inverse
    else if (!less(it, _to_canonical)) {
      auto it_inv = it.inverse();
      if (it_inv < _from_canonical) {
        _from_canonical = it_inv;
      }
    }
  }
  return _from_canonical;
}

template <typename Base>
template <typename PermuteIteratorIt>
std::vector<PermuteIterator> ConfigCanonicalForm<Base>::invariant_subgroup(
    PermuteIteratorIt begin, PermuteIteratorIt end) const {
  // std::cout << "INSIDE invariant_subgroup()" << std::endl;
  std::vector<PermuteIterator> sub_grp;
  // std::cout << "begin: " << begin->factor_group_index() << std::endl;
  std::copy_if(begin, end, std::back_inserter(sub_grp), derived().equal_to());
  // std::cout << "Complete!\n";
  return sub_grp;
}

// --- template<typename Base> class SupercellCanonicalForm

template <typename Base>
bool SupercellCanonicalForm<Base>::is_canonical() const {
  return xtal::canonical::check(derived().lattice(),
                                derived().prim().point_group());
}

template <typename Base>
SymOp SupercellCanonicalForm<Base>::to_canonical() const {
  auto to_canonical_ix = xtal::canonical::operation_index(
      derived().lattice(), derived().prim().point_group());
  return derived().prim().point_group()[to_canonical_ix];
}

template <typename Base>
SymOp SupercellCanonicalForm<Base>::from_canonical() const {
  return to_canonical().inverse();
}

template <typename Base>
Lattice SupercellCanonicalForm<Base>::canonical_lattice() const {
  return xtal::canonical::equivalent(derived().lattice(),
                                     derived().prim().point_group(),
                                     derived().crystallography_tol());
}

template <typename Base>
const Supercell &SupercellCanonicalForm<Base>::canonical_form() const {
  if (!m_canonical) {
    m_canonical = &*derived().insert().first;
  }
  return *m_canonical;
}

/// \brief Construct the subgroup of permutations that leaves a Supercell
/// unchanged
///
/// \param scel_B Supercell associated with the supergroup [begin, end)
/// \param begin,end Range of PermuteIterator describing the supergroup
///
/// - 'this' Supercell should be a supercell of (or the same as) scel_B
template <typename Base>
template <typename PermuteIteratorIt>
std::vector<PermuteIterator> SupercellCanonicalForm<Base>::invariant_subgroup(
    const Supercell &scel_B, PermuteIteratorIt begin, PermuteIteratorIt end) {
  return make_invariant_subgroup(derived(), scel_B, begin, end);
}

}  // namespace CASM

#endif
