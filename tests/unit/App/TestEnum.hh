#ifndef CASM_TestEnum
#define CASM_TestEnum

#include "casm/container/Counter.hh"
#include "casm/container/InputEnumerator.hh"
#include "casm/clex/Configuration.hh"
#include "casm/misc/cloneable_ptr.hh"

ENUMERATOR_INTERFACE_TRAITS(TestEnum)

namespace CASM {

  /// \brief Enumerate over all possible occupations in a particular Supercell
  ///
  /// \ingroup ConfigEnum
  ///
  class TestEnum : public InputEnumeratorBase<Configuration> {

    // -- Required members -------------------

  public:

    /// \brief Construct with a Supercell, using all permutations
    TestEnum(Supercell &_scel);

    ENUMERATOR_MEMBERS(TestEnum)

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
