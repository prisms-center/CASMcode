#ifndef CASM_ScelEnum
#define CASM_ScelEnum

#include "casm/misc/cloneable_ptr.hh"
#include "casm/crystallography/SuperlatticeEnumerator.hh"
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

    std::string name() const override;

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
  /// - Enumerated Supercell are not canonical
  /// - References are invalidated after incrementing an iterator
  ///
  class ScelEnumByProps : public InputEnumeratorBase<Supercell> {

  public:

    /// \brief Construct with shared prim Structure and ScelEnumProps settings
    ScelEnumByProps(std::shared_ptr<const Structure> &shared_prim, const ScelEnumProps &enum_props);

    /// \brief Construct with PrimClex and ScelEnumProps settings
    ScelEnumByProps(const PrimClex &primclex, const ScelEnumProps &enum_props, bool existing_only = false);

    /// \brief Construct with PrimClex and ScelEnumProps JSON settings
    ScelEnumByProps(const PrimClex &primclex, const jsonParser &input);

    ScelEnumByProps(const ScelEnumByProps &) = delete;
    ScelEnumByProps &operator=(const ScelEnumByProps &) = delete;
    ~ScelEnumByProps() {}


    std::string name() const override;

    static const std::string enumerator_name;


  private:

    /// Check for existing supercells
    bool _include(const Lattice &lat) const;

    /// Implements increment over supercells
    void increment() override;

    std::shared_ptr<Structure const> m_shared_prim;
    notstd::cloneable_ptr<Supercell> m_current;
    PrimClex const *m_primclex;

    std::unique_ptr<SuperlatticeEnumerator > m_lattice_enum;
    SuperlatticeEnumerator::const_iterator m_lat_it;
    SuperlatticeEnumerator::const_iterator m_lat_end;

    bool m_existing_only;
  };


  /// \brief Enumerate over Supercell
  ///
  /// - Provides a unified Interface for ScelEnumByName and ScelEnumByProps
  /// - Enumerated Supercell are not
  /// - If ScelEnumByProps, references are invalidated after incrementing an
  ///   iterator
  ///
  class ScelEnum : public InputEnumeratorBase<Supercell> {

  public:

    /// \brief Construct with PrimClex and JSON settings
    ScelEnum(const PrimClex &primclex, const jsonParser &input);

    ScelEnum(const ScelEnum &) = delete;
    ScelEnum &operator=(const ScelEnum &) = delete;


    std::string name() const override;

    static const std::string enumerator_name;
    static std::string interface_help();
    static int run(const PrimClex &primclex, const jsonParser &kwargs, const Completer::EnumOption &enum_opt, EnumeratorMap const *interface_map);


  private:

    /// Implements increment over all occupations
    void increment() override;

    InputEnumIterator<Supercell> m_it;
    InputEnumIterator<Supercell> m_end;
    InputEnumerator<Supercell> m_enum;
  };

  /// \brief Read unit cell transformation matrix from JSON input
  ///
  /// Input json with one of the following forms:
  /// - \code
  ///   {
  ///     "unit_cell" : [
  ///       [X, X, X],
  ///       [X, X, X],
  ///       [X, X, X]
  ///     ]
  ///   }
  ///   \endcode
  /// - \code
  ///   { "unit_cell" : "SCEL..."}
  ///   \endcode
  /// - If input is null or does not contain "unit_cell", returns identity matrix.
  ///
  Eigen::Matrix3i make_unit_cell(const PrimClex &primclex, const jsonParser &json);

  /// \brief Make a ScelEnumProps object from JSON input
  ///
  /// - min: int (default=1)
  /// - max: int (default=max existing supercell size)
  /// - dirs: string, (default="abc")
  /// - unit_cell: 3x3 matrix of int, or string (default=identity matrix)
  /// \code
  /// {
  ///   "min" : 1,
  ///   "max" : 5,
  ///   "dirs" : "abc",
  ///   "unit_cell" : [
  ///     [0, 0, 0],
  ///     [0, 0, 0],
  ///     [0, 0, 0]
  ///   ],
  ///   "unit_cell" : "SCEL...",
  /// }
  /// \endcode
  ///
  ScelEnumProps make_scel_enum_props(const PrimClex &primclex, const jsonParser &input);

  template <>
  struct jsonConstructor<xtal::ScelEnumProps> {
    static xtal::ScelEnumProps from_json(const jsonParser &json, const PrimClex &primclex);
  };


}

/** @}*/

#endif
