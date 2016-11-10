#ifndef CASM_ConfigEnumAllOccupations
#define CASM_ConfigEnumAllOccupations

#include "casm/container/Counter.hh"
#include "casm/container/InputEnumerator.hh"
#include "casm/clex/Configuration.hh"
#include "casm/misc/cloneable_ptr.hh"

namespace CASM {

  ENUMERATOR_INTERFACE_TRAITS(ConfigEnumAllOccupations)

  /// \brief Enumerate over all possible occupations in a particular Supercell
  ///
  /// \ingroup ConfigEnum
  ///
  class ConfigEnumAllOccupations : public InputEnumeratorBase<Configuration> {

    // -- Required members -------------------

  public:

    /// \brief Construct with a Supercell, using all permutations
    ConfigEnumAllOccupations(Supercell &_scel);

    ENUMERATOR_MEMBERS(ConfigEnumAllOccupations)

  private:

    /// Implements increment
    void increment() override;


    // -- Unique -------------------

    /// Returns true if current() is primitive and canonical
    bool _check_current() const;

    Counter<Array<int> > m_counter;
    notstd::cloneable_ptr<Configuration> m_current;
  };

}

#endif
