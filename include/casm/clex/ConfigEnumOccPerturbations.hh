/*
#ifndef CASM_ConfigEnumOccPerturbations
#define CASM_ConfigEnumOccPerturbations

#include "casm/container/Counter.hh"
#include "casm/enumerator/InputEnumerator.hh"
#include "casm/clex/Configuration.hh"
#include "casm/misc/cloneable_ptr.hh"

extern "C" {
  CASM::EnumInterfaceBase *make_ConfigEnumOccPerturbations_interface();
}

namespace CASM {

  /** \defgroup ConfigEnumGroup Configuration Enumerators
   *  \ingroup Configuration
   *  \ingroup Enumerator
   *  \brief Enumerates Configuration
   *  @{
  */
/*
  /// \brief Enumerate occupation perturbations about a specified Configuration
  ///
  class ConfigEnumOccPerturbations : public InputEnumeratorBase<Configuration> {

    // -- Required members -------------------

  public:

    ConfigEnumOccPerturbations(const Configuration &_background_config);

    std::string name() const override {
      return enumerator_name;
    }

    static const std::string enumerator_name;
    static std::string interface_help();

    static int run(const PrimClex &primclex, const jsonParser &kwargs, const Completer::EnumOption &enum_opt);

  private:


    /// Implements increment
    void increment() override;


    // -- Unique -------------------

  };

  /** @}*/
/*
}

#endif
*/
