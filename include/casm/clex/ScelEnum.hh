#ifndef CASM_ScelEnum
#define CASM_ScelEnum

#include "casm/misc/cloneable_ptr.hh"
#include "casm/crystallography/SupercellEnumerator.hh"
#include "casm/container/InputEnumerator.hh"
#include "casm/container/RandomAccessEnumerator.hh"
#include "casm/clex/Supercell.hh"

ENUMERATOR_VARIABLECONST_TRAITS(ScelEnumByNameT)

namespace CASM {

  /** \defgroup ScelEnum
   *
   *  \ingroup Enumerator
   *
   *  @{
  */


  /// \brief Enumerate over Supercell
  ///
  /// - Specify Supercell by providing a list of names of Supercell already
  ///   included in PrimClex
  ///
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

  /// \relates ScelEnumByNameT
  typedef ScelEnumByNameT<true> ConstScelEnumByName;

  /// \relates ScelEnumByNameT
  typedef ScelEnumByNameT<false> ScelEnumByName;
}


ENUMERATOR_VARIABLECONST_TRAITS(ScelEnumByPropsT)

namespace CASM {

  /// \brief Enumerate over Supercell
  ///
  /// - Specify Supercell using ScelEnumProps (min/max volume, dirs, unit_cell)
  /// - Enumerated Supercell are canonical, included in PrimClex
  ///
  template<bool IsConst = true>
  class ScelEnumByPropsT : public InputEnumeratorBase<Supercell, IsConst> {

  public:

    /// \brief Construct with PrimClex and ScelEnumProps settings
    ScelEnumByPropsT(PrimClex &primclex, const ScelEnumProps &enum_props, bool existing_only = false);

    /// \brief Construct with PrimClex and ScelEnumProps JSON settings
    ScelEnumByPropsT(PrimClex &primclex, const jsonParser &input);

    ScelEnumByPropsT(const ScelEnumByPropsT &) = delete;
    ScelEnumByPropsT &operator=(const ScelEnumByPropsT &) = delete;


    ENUMERATOR_MEMBERS(ScelEnumByPropsT)

  private:

    /// Check for existing supercells
    bool _include(const Lattice &lat) const;

    /// Implements increment over supercells
    void increment() override;

    PrimClex *m_primclex;

    std::unique_ptr<SupercellEnumerator<Lattice> > m_lattice_enum;
    SupercellEnumerator<Lattice>::const_iterator m_lat_it;
    SupercellEnumerator<Lattice>::const_iterator m_lat_end;

    bool m_existing_only;
  };

  /// \relates ScelEnumByPropsT
  typedef ScelEnumByPropsT<true> ConstScelEnumByProps;

  /// \relates ScelEnumByPropsT
  typedef ScelEnumByPropsT<false> ScelEnumByProps;

}


ENUMERATOR_INTERFACE_VARIABLECONST_TRAITS(ScelEnumT)

namespace CASM {

  /// \brief Enumerate over Supercell
  ///
  /// - Provides a unified Interface for ScelEnumByName and ScelEnumByProps
  /// - Enumerated Supercell are canonical, included in PrimClex
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

  /// \relates ScelEnumT
  typedef ScelEnumT<true> ConstScelEnum;

  /// \relates ScelEnumT
  typedef ScelEnumT<false> ScelEnum;

  /** @}*/
}

#endif
