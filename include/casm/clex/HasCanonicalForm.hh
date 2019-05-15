#ifndef CASM_HasCanonicalForm
#define CASM_HasCanonicalForm

#include <vector>

namespace CASM {
  class SymGroup;
  class SymOp;
  class PermuteIterator;
  class Supercell;

  /// \brief Implements canonical form finding when using SymCompare
  ///
  template<typename Base>
  class CanonicalForm: public Base {
  public:
    typedef typename Base::MostDerived MostDerived;
    using Base::derived;


    /// Check if canonical
    template<typename SymCompareType>
    bool is_canonical(
      const std::vector<SymOp> &g,
      const SymCompareType &sym_compare) const;

    /// Return canonical form
    template<typename SymCompareType>
    MostDerived canonical_form(
      const std::vector<SymOp> &g,
      const SymCompareType &sym_compare) const;

    /// Check if two elements have the same canonical form
    template<typename SymCompareType>
    bool is_sym_equivalent(
      const MostDerived &other,
      const std::vector<SymOp> &g,
      const SymCompareType &sym_compare) const;

    /// Find element that has the same canonical form
    template<typename ObjIterator, typename SymCompareType>
    ObjIterator find_sym_equivalent(
      ObjIterator begin,
      ObjIterator end,
      const std::vector<SymOp> &g,
      const SymCompareType &sym_compare) const;

    /// Return op that transforms this into canonical form
    template<typename SymCompareType>
    SymOp to_canonical(
      const std::vector<SymOp> &g,
      const SymCompareType &sym_compare) const;

    /// Return op that transforms the canonical form into this
    template<typename SymCompareType>
    SymOp from_canonical(
      const std::vector<SymOp> &g,
      const SymCompareType &sym_compare) const;

    /// Return subgroup of super_grp that leaves this invariant
    template<typename SymCompareType>
    SymGroup invariant_subgroup(
      const SymGroup &super_grp,
      const SymCompareType &sym_compare) const;


    /// Check if canonical in Supercell
    bool is_canonical(const Supercell &scel) const;

    /// Return canonical form in Supercell
    MostDerived canonical_form(const Supercell &scel) const;

    /// True if this and B have same canonical form in Supercell
    bool is_sym_equivalent(
      const MostDerived &B,
      const Supercell &scel) const;

    /// Return operation in Supercell that transforms this into canonical form
    SymOp to_canonical(const Supercell &scel) const;

    /// Return operation in Supercell that transforms canonical form into this
    SymOp from_canonical(const Supercell &scel) const;

    /// Return subgroup of Supercell permutations that leave this invariant
    std::vector<PermuteIterator> invariant_subgroup(const Supercell &scel) const;


    /// Check if canonical in Supercell, with respect to a subgroup of Supercell permutations
    template<typename PermuteIteratorIt>
    bool is_canonical(
      const Supercell &scel,
      PermuteIteratorIt begin,
      PermuteIteratorIt end) const;

    /// Return canonical form in Supercell, with respect to a subgroup of Supercell permutations
    template<typename PermuteIteratorIt>
    MostDerived canonical_form(
      const Supercell &scel,
      PermuteIteratorIt begin,
      PermuteIteratorIt end) const;

    /// True if this and B have same canonical form, with respect to a subgroup
    /// of Supercell permutations
    template<typename PermuteIteratorIt>
    bool is_sym_equivalent(
      const MostDerived &B,
      const Supercell &scel,
      PermuteIteratorIt begin,
      PermuteIteratorIt end) const;

    /// Find element that has the same canonical form, with respect to a subgroup
    template<typename ObjIterator, typename PermuteIteratorIt>
    ObjIterator find_sym_equivalent(
      ObjIterator obj_begin,
      ObjIterator obj_end,
      const Supercell &scel,
      PermuteIteratorIt begin,
      PermuteIteratorIt end) const;

    /// Return operation in Supercell that transforms this into canonical form,
    /// with respect to a subgroup of Supercell permutations
    template<typename PermuteIteratorIt>
    SymOp to_canonical(
      const Supercell &scel,
      PermuteIteratorIt begin,
      PermuteIteratorIt end) const;

    /// Return operation in Supercell that transforms canonical form into this,
    /// with respect to a subgroup of Supercell permutations
    template<typename PermuteIteratorIt>
    SymOp from_canonical(
      const Supercell &scel,
      PermuteIteratorIt begin,
      PermuteIteratorIt end) const;

    /// Return subgroup of Supercell permutations that leave this invariant,
    /// with respect to a subgroup of Supercell permutations
    template<typename PermuteIteratorIt>
    std::vector<PermuteIterator> invariant_subgroup(
      const Supercell &scel,
      PermuteIteratorIt begin,
      PermuteIteratorIt end) const;

  };

  /// \brief Implements canonical form finding for Configuration and DiffTransConfiguration
  ///
  /// Requires MostDerived implements:
  /// - Functor less() const;
  ///   - Functor is a `ConfigCompare` like functor that can compare A*config < B*other
  /// - SuperGroupPermuteIteratorType default_permute_begin() const;
  /// - SuperGroupPermuteIteratorType default_permute_end() const;
  ///   - where SuperGroupPermuteIteratorType may be PermuteIterator or
  ///     std::vector<PermuteIterator>::const_iterator
  ///
  template<typename Base>
  class ConfigCanonicalForm: public Base {
  public:
    typedef typename Base::MostDerived MostDerived;
    using Base::derived;


    /// True if this and B have same canonical form
    bool is_sym_equivalent(const MostDerived &B) const;

    /// Find Config that has same canonical form
    template<typename ConfigIterator>
    ConfigIterator find_sym_equivalent(
      const MostDerived &B,
      ConfigIterator obj_begin,
      ConfigIterator obj_end) const;

    bool is_canonical() const;

    MostDerived canonical_form() const;

    PermuteIterator to_canonical() const;

    PermuteIterator from_canonical() const;

    std::vector<PermuteIterator> invariant_subgroup() const;


    /// True if this and B have same canonical form
    template<typename PermuteIteratorIt>
    bool is_sym_equivalent(const MostDerived &B, PermuteIteratorIt begin, PermuteIteratorIt end) const;

    /// True if this and B have same canonical form
    template<typename ConfigIterator, typename PermuteIteratorIt>
    ConfigIterator find_sym_equivalent(
      ConfigIterator obj_begin,
      ConfigIterator obj_end,
      PermuteIteratorIt begin,
      PermuteIteratorIt end) const;

    template<typename PermuteIteratorIt>
    bool is_canonical(PermuteIteratorIt begin, PermuteIteratorIt end) const;

    template<typename PermuteIteratorIt>
    MostDerived canonical_form(PermuteIteratorIt begin, PermuteIteratorIt end) const;

    template<typename PermuteIteratorIt>
    PermuteIterator to_canonical(PermuteIteratorIt begin, PermuteIteratorIt end) const;

    template<typename PermuteIteratorIt>
    PermuteIterator from_canonical(PermuteIteratorIt begin, PermuteIteratorIt end) const;

    template<typename PermuteIteratorIt>
    std::vector<PermuteIterator> invariant_subgroup(PermuteIteratorIt begin, PermuteIteratorIt end) const;

    // --- Required in MostDerived:

    //Functor less() const;
    //const Supercell& supercell() const;

  };

  /// Supercell canonical form finding is a special case that returns references
  ///   to Supercell in the Database<Supercell>
  template<typename Base>
  class SupercellCanonicalForm : public Base {
  public:

    using Base::derived;

    SupercellCanonicalForm() : m_canonical(nullptr) {}

    /// True if lattice().is_canonical(g), where g is the prim().point_group()
    bool is_canonical() const;

    /// The first SymOp in the derived().prim().point_group() for which:
    ///   derived().lattice().is_equivalent(op, canonical_form().lattice())) == true,
    SymOp to_canonical() const;

    /// The inverse of to_canonical()
    SymOp from_canonical() const;

    /// \brief Return canonical equivalent lattice
    Lattice canonical_lattice() const;

    /// \brief Return canonical equivalent Supercell
    ///
    /// - Will be inserted in Database if necessary
    const Supercell &canonical_form() const;

    /// \brief Construct the subgroup of permutations that leaves a Supercell unchanged
    template<typename PermuteIteratorIt>
    std::vector<PermuteIterator> invariant_subgroup(
      const Supercell &scel_B,
      PermuteIteratorIt begin,
      PermuteIteratorIt end);

  private:

    /// Store a pointer to the canonical equivalent Supercell
    mutable const Supercell *m_canonical;

  };

}

#endif
