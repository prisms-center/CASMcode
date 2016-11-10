#ifndef CASM_ScelEnum
#define CASM_ScelEnum

#include "casm/misc/cloneable_ptr.hh"
#include "casm/crystallography/SupercellEnumerator.hh"
#include "casm/container/InputEnumerator.hh"
#include "casm/container/RandomAccessEnumerator.hh"
#include "casm/clex/Supercell.hh"

namespace CASM {

  ENUMERATOR_VARIABLECONST_TRAITS(ScelEnumByNameT)

  template<bool IsConst = true>
  class ScelEnumByNameT : public RandomAccessEnumeratorBase<Supercell, IsConst> {

    // -- Required members -------------------

  public:

    using typename RandomAccessEnumeratorBase<Supercell, IsConst>::step_type;

    /// \brief Construct with PrimClex and ScelEnumProps settings
    ScelEnumByNameT(PrimClex &primclex, std::initializer_list<std::string> scelnames);

    /// \brief Construct with PrimClex and ScelEnumProps settings
    template<typename ScelNameIterator>
    ScelEnumByNameT(PrimClex &primclex, ScelNameIterator begin, ScelNameIterator end);

    /// \brief Construct with PrimClex and array of supercell names
    ScelEnumByNameT(PrimClex &primclex, const jsonParser &input);

    ENUMERATOR_MEMBERS(ScelEnumByNameT)

  private:

    /// Implements at_step
    Supercell *at_step(step_type n) override;

    // -- Unique -------------------

    void _init();

    PrimClex *m_primclex;

    std::vector<Supercell *> m_scelptr;
  };

  typedef ScelEnumByNameT<true> ConstScelEnumByName;
  typedef ScelEnumByNameT<false> ScelEnumByName;



  ENUMERATOR_VARIABLECONST_TRAITS(ScelEnumByPropsT)

  template<bool IsConst = true>
  class ScelEnumByPropsT : public InputEnumeratorBase<Supercell, IsConst> {

  public:

    /// \brief Construct with PrimClex and ScelEnumProps settings
    ScelEnumByPropsT(PrimClex &primclex, const ScelEnumProps &enum_props);

    /// \brief Construct with PrimClex and ScelEnumProps JSON settings
    ScelEnumByPropsT(PrimClex &primclex, const jsonParser &input);

    ScelEnumByPropsT(const ScelEnumByPropsT &) = delete;
    ScelEnumByPropsT &operator=(const ScelEnumByPropsT &) = delete;


    ENUMERATOR_MEMBERS(ScelEnumByPropsT)

  private:

    /// Implements increment over supercells
    void increment() override;

    PrimClex *m_primclex;

    std::unique_ptr<SupercellEnumerator<Lattice> > m_lattice_enum;
    SupercellEnumerator<Lattice>::const_iterator m_lat_it;
    SupercellEnumerator<Lattice>::const_iterator m_lat_end;
  };

  typedef ScelEnumByPropsT<true> ConstScelEnumByProps;
  typedef ScelEnumByPropsT<false> ScelEnumByProps;



  ENUMERATOR_INTERFACE_VARIABLECONST_TRAITS(ScelEnumT)

  /// \brief Enumerate over Supercell
  ///
  /// - Unified Interface for ScelEnumByName and ScelEnumByProps
  ///
  /// \ingroup ScelEnum
  ///
  template<bool IsConst = true>
  class ScelEnumT : public InputEnumeratorBase<Supercell, IsConst> {

  public:

    /// \brief Construct with PrimClex and JSON settings
    ScelEnumT(PrimClex &primclex, const jsonParser &input);

    ScelEnumT(const ScelEnumT &) = delete;
    ScelEnumT &operator=(const ScelEnumT &) = delete;


    ENUMERATOR_MEMBERS(ScelEnumT)

  private:

    /// Implements increment over all occupations
    void increment() override;

    InputEnumIterator<Supercell, false> m_it;
    InputEnumIterator<Supercell, false> m_end;
    InputEnumerator<Supercell, false> m_enum;
  };

  typedef ScelEnumT<true> ConstScelEnum;
  typedef ScelEnumT<false> ScelEnum;


}

#endif
