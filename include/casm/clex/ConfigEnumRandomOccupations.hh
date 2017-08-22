#ifndef CASM_ConfigEnumRandomOccupations
#define CASM_ConfigEnumRandomOccupations

#include "casm/container/Counter.hh"
#include "casm/container/InputEnumerator.hh"
#include "casm/clex/Configuration.hh"
#include "casm/misc/cloneable_ptr.hh"

extern "C" {
  CASM::EnumInterfaceBase *make_ConfigEnumRandomOccupations_interface();
}

class MTRand;

namespace CASM {

  /** \defgroup ConfigEnumGroup Configuration Enumerators
   *  \ingroup Configuration
   *  \ingroup Enumerator
   *  \brief Enumerates Configuration
   *  @{
  */

  /// \brief Enumerate n random occupations in a particular Supercell
  ///
  class ConfigEnumRandomOccupations : public InputEnumeratorBase<Configuration> {

    // -- Required members -------------------

  public:

    /// \brief Construct with a Supercell, using all permutations
    ConfigEnumRandomOccupations(
      Supercell &_scel,
      Index _n_config,
      MTRand &_mtrand);

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

    void randomize();

    Index m_n_config;
    MTRand *m_mtrand;
    Array<int> m_max_allowed;
    notstd::cloneable_ptr<Configuration> m_current;
  };

  /** @}*/
}

#endif
