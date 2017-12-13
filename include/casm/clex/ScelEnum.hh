#ifndef CASM_ScelEnum
#define CASM_ScelEnum

#include "casm/misc/cloneable_ptr.hh"
#include "casm/crystallography/SupercellEnumerator.hh"
#include "casm/enumerator/InputEnumerator.hh"
#include "casm/enumerator/RandomAccessEnumerator.hh"
#include "casm/clex/Supercell.hh"

/** \defgroup ScelEnumGroup Supercell Enumerators
 *
 *  \ingroup Enumerator
 *  \ingroup Supercell
 *  \brief Enumerates Supercell
 *  @{
*/

extern "C" {
  CASM::EnumInterfaceBase *make_ScelEnum_interface();
}

namespace CASM {


  /// \brief Enumerate over Supercell
  ///
  /// - Specify Supercell by providing a list of names of Supercell already
  ///   included in PrimClex
  ///
  class ScelEnumByName : public RandomAccessEnumeratorBase<Supercell> {

    // -- Required members -------------------

  public:

    using typename RandomAccessEnumeratorBase<Supercell>::step_type;

    /// \brief Construct with PrimClex and ScelEnumProps settings
    ScelEnumByName(const PrimClex &primclex, std::initializer_list<std::string> scelnames);

    /// \brief Construct with PrimClex and ScelEnumProps settings
    template<typename ScelNameIterator>
    ScelEnumByName(const PrimClex &primclex, ScelNameIterator begin, ScelNameIterator end);

    /// \brief Construct with PrimClex and array of supercell names
    ScelEnumByName(const PrimClex &primclex, const jsonParser &input);

    std::string name() const override {
      return enumerator_name;
    }

    static const std::string enumerator_name;


  private:

    /// Implements at_step
    const Supercell *at_step(step_type n) override;

    // -- Unique -------------------

    void _init();

    const PrimClex *m_primclex;

    std::vector<const Supercell *> m_scelptr;
  };


  /// \brief Enumerate over Supercell
  ///
  /// - Specify Supercell using ScelEnumProps (min/max volume, dirs, unit_cell)
  /// - Enumerated Supercell are canonical, included in PrimClex
  /// - References are invalidated after incrementing an iterator
  ///
  class ScelEnumByProps : public InputEnumeratorBase<Supercell> {

  public:

    /// \brief Construct with PrimClex and ScelEnumProps settings
    ScelEnumByProps(const PrimClex &primclex, const ScelEnumProps &enum_props, bool existing_only = false);

    /// \brief Construct with PrimClex and ScelEnumProps JSON settings
    ScelEnumByProps(const PrimClex &primclex, const jsonParser &input);

    ScelEnumByProps(const ScelEnumByProps &) = delete;
    ScelEnumByProps &operator=(const ScelEnumByProps &) = delete;


    std::string name() const override {
      return enumerator_name;
    }

    static const std::string enumerator_name;


  private:

    /// Check for existing supercells
    bool _include(const Lattice &lat) const;

    /// Implements increment over supercells
    void increment() override;

    const PrimClex *m_primclex;

    std::unique_ptr<SupercellEnumerator<Lattice> > m_lattice_enum;
    SupercellEnumerator<Lattice>::const_iterator m_lat_it;
    SupercellEnumerator<Lattice>::const_iterator m_lat_end;

    bool m_existing_only;
  };


  /// \brief Enumerate over Supercell
  ///
  /// - Provides a unified Interface for ScelEnumByName and ScelEnumByProps
  /// - Enumerated Supercell are canonical, included in PrimClex
  /// - If ScelEnumByProps, references are invalidated after incrementing an
  ///   iterator
  ///
  class ScelEnum : public InputEnumeratorBase<Supercell> {

  public:

    /// \brief Construct with PrimClex and JSON settings
    ScelEnum(const PrimClex &primclex, const jsonParser &input);

    ScelEnum(const ScelEnum &) = delete;
    ScelEnum &operator=(const ScelEnum &) = delete;


    std::string name() const override {
      return enumerator_name;
    }

    static const std::string enumerator_name;
    static std::string interface_help();
    static int run(const PrimClex &primclex, const jsonParser &kwargs, const Completer::EnumOption &enum_opt);


  private:

    /// Implements increment over all occupations
    void increment() override;

    InputEnumIterator<Supercell> m_it;
    InputEnumIterator<Supercell> m_end;
    InputEnumerator<Supercell> m_enum;
  };

}

/** @}*/

#endif
