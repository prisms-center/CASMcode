#ifndef CASM_ScelEnum
#define CASM_ScelEnum

#include "casm/clex/Supercell.hh"
#include "casm/crystallography/SuperlatticeEnumerator.hh"
#include "casm/enumerator/InputEnumerator.hh"
#include "casm/misc/cloneable_ptr.hh"

/** \defgroup ScelEnumGroup Supercell Enumerators
 *
 *  \ingroup Enumerator
 *  \ingroup Supercell
 *  \brief Enumerates Supercell
 *  @{
 */

namespace CASM {

class ScelEnumByProps;

/// Always true for ScelEnumByProps
template <>
bool is_guaranteed_for_database_insert(ScelEnumByProps const &enumerator);

/// Enumerate symmetrically unique Supercell
///
/// - Specify Supercell using xtal::ScelEnumProps (min/max volume, dirs,
/// unit_cell)
/// - Enumerated Supercell are canonical
/// - Enumerated Supercell do not have a `PrimClex const *`
/// - References are invalidated after incrementing an iterator
///
class ScelEnumByProps : public InputEnumeratorBase<Supercell> {
 public:
  /// \brief Construct with shared prim Structure and ScelEnumProps settings
  ScelEnumByProps(std::shared_ptr<const Structure> const &shared_prim,
                  const xtal::ScelEnumProps &enum_props);

  ScelEnumByProps(const ScelEnumByProps &) = delete;
  ScelEnumByProps &operator=(const ScelEnumByProps &) = delete;
  ~ScelEnumByProps() {}

  std::string name() const override;

  static const std::string enumerator_name;

 private:
  /// Implements increment over supercells
  void increment() override;

  std::shared_ptr<Structure const> m_shared_prim;
  notstd::cloneable_ptr<Supercell> m_current;

  std::unique_ptr<xtal::SuperlatticeEnumerator> m_lattice_enum;
  xtal::SuperlatticeEnumerator::const_iterator m_lat_it;
  xtal::SuperlatticeEnumerator::const_iterator m_lat_end;
};

}  // namespace CASM

/** @}*/

#endif
