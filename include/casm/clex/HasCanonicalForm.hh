#ifndef CASM_HasCanonicalForm
#define CASM_HasCanonicalForm

namespace CASM {

  class SymOp;
  class SymGroup;
  class PermuteIterator;

  /// \brief Implements canonical form finding when using SymCompare
  ///
  template<typename Base>
  class CanonicalForm: public Base {
  public:
    typedef typename Base::MostDerived MostDerived;
    using Base::derived;

    template<typename SymCompareType>
    bool is_canonical(
      const SymGroup &g,
      const SymCompareType &sym_compare) const;

    template<typename SymCompareType>
    MostDerived canonical_form(
      const SymGroup &g,
      const SymCompareType &sym_compare) const;

    template<typename SymCompareType>
    bool is_equivalent(
      const MostDerived &other,
      const SymGroup &g,
      const SymCompareType &sym_compare) const;

    template<typename SymCompareType>
    SymOp to_canonical(
      const SymGroup &g,
      const SymCompareType &sym_compare) const;

    template<typename SymCompareType>
    SymOp from_canonical(
      const SymGroup &g,
      const SymCompareType &sym_compare) const;

    template<typename SymCompareType>
    SymGroup invariant_subgroup(
      const SymGroup &super_grp,
      const SymCompareType &sym_compare) const;


    template<typename PermuteIteratorIt>
    bool is_canonical(
      const Supercell &scel,
      PermuteIteratorIt begin,
      PermuteIteratorIt end) const;

    template<typename PermuteIteratorIt>
    MostDerived canonical_form(
      const Supercell &scel,
      PermuteIteratorIt begin,
      PermuteIteratorIt end) const;

    /// True if this and B have same canonical form
    template<typename PermuteIteratorIt>
    bool is_equivalent(
      const MostDerived &B,
      const Supercell &scel,
      PermuteIteratorIt begin,
      PermuteIteratorIt end) const;

    template<typename PermuteIteratorIt>
    SymOp to_canonical(
      const Supercell &scel,
      PermuteIteratorIt begin,
      PermuteIteratorIt end) const;

    template<typename PermuteIteratorIt>
    SymOp from_canonical(
      const Supercell &scel,
      PermuteIteratorIt begin,
      PermuteIteratorIt end) const;

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
    bool is_equivalent(const MostDerived &B) const;

    bool is_canonical() const;

    MostDerived canonical_form() const;

    PermuteIterator to_canonical() const;

    PermuteIterator from_canonical() const;

    std::vector<PermuteIterator> invariant_subgroup() const;


    /// True if this and B have same canonical form
    template<typename PermuteIteratorIt>
    bool is_equivalent(const MostDerived &B, PermuteIteratorIt begin, PermuteIteratorIt end) const;

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

    bool is_canonical() const;

    SymOp to_canonical() const;

    SymOp from_canonical() const;

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
