#ifndef CASM_ConfigEnumAllOccupations
#define CASM_ConfigEnumAllOccupations

#include "casm/container/Counter.hh"
#include "casm/container/InputEnumerator.hh"
#include "casm/clex/Configuration.hh"
#include "casm/misc/cloneable_ptr.hh"

extern "C" {
  CASM::EnumInterfaceBase *make_ConfigEnumAllOccupations_interface();
}

namespace CASM {

  /** \defgroup ConfigEnumGroup Configuration Enumerators
   *  \ingroup Configuration
   *  \ingroup Enumerator
   *  \brief Enumerates Configuration
   *  @{
  */

  /// \brief Enumerate over all possible occupations in a particular Supercell
  ///
  class ConfigEnumAllOccupations : public InputEnumeratorBase<Configuration> {

    // -- Required members -------------------

  public:

    /// \brief Construct with a Supercell, using all permutations
    ConfigEnumAllOccupations(Supercell &_scel);

    std::string name() const override {
      return enumerator_name;
    }

    static const std::string enumerator_name;
    static const std::string interface_help;
    static int run(PrimClex &primclex, const jsonParser &kwargs, const Completer::EnumOption &enum_opt);

  private:


    /// Implements increment
    void increment() override;


    // -- Unique -------------------

    /// Returns true if current() is primitive and canonical
    bool _check_current() const;

    Counter<Array<int> > m_counter;
    notstd::cloneable_ptr<Configuration> m_current;
  };

  /** @}*/
}

#endif
