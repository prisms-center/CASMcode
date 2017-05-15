#ifndef CASM_TestEnum
#define CASM_TestEnum

#include "casm/container/Counter.hh"
#include "casm/container/InputEnumerator.hh"
#include "casm/clex/Configuration.hh"
#include "casm/misc/cloneable_ptr.hh"

extern "C" {
  CASM::EnumInterfaceBase *make_TestEnum_interface();
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
  class TestEnum : public InputEnumeratorBase<Configuration> {

    // -- Required members -------------------

  public:

    /// \brief Construct with a Supercell, using all permutations
    TestEnum(const Supercell &_scel);

    std::string name() const override {
      return enumerator_name;
    }

    static const std::string enumerator_name;
    static const std::string interface_help;

    static int run(const PrimClex &primclex, const jsonParser &kwargs, const Completer::EnumOption &enum_opt);

    template<typename ScelIterator>
    static int run(
      const PrimClex &primclex,
      ScelIterator begin,
      ScelIterator end,
      const std::vector<std::string> &filter_expr = {});

  private:


    /// Implements increment
    void increment() override;


    // -- Unique -------------------

    /// Returns true if current() is primitive and canonical
    bool _check_current() const;

    Counter<std::vector<int> > m_counter;
    notstd::cloneable_ptr<Configuration> m_current;
  };

  /** @}*/
}

#endif
