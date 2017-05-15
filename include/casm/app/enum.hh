#ifndef CASM_enum
#define CASM_enum

#include "casm/app/APICommand.hh"
#include "casm/completer/Handlers.hh"

namespace CASM {
  namespace Completer {

    /// Options set for `casm enum`.
    class EnumOption : public OptionHandlerBase {

    public:

      EnumOption();

      using OptionHandlerBase::settings_path;
      using OptionHandlerBase::input_str;
      using OptionHandlerBase::supercell_strs;

      const std::vector<std::string> &desc_vec() const {
        return m_desc_vec;
      }

      std::string method() const {
        return m_method;
      }

      int min_volume() const {
        return m_min_volume;
      }

      int max_volume() const {
        return m_max_volume;
      }

      bool all_existing() const {
        return m_all_existing;
      }

      const std::vector<std::string> &filter_strs() const {
        return m_filter_strs;
      }

    private:

      void initialize() override;

      std::vector<std::string> m_desc_vec;

      std::string m_method;
      int m_min_volume;
      int m_max_volume;
      bool m_all_existing;
      std::vector<std::string> m_filter_strs;

    };

  }
}

namespace CASM {

  typedef InterfaceMap<Completer::EnumOption> EnumeratorMap;

  /// 'casm enum' implementation
  class EnumCommand : public APICommand<Completer::EnumOption> {

  public:

    static const std::string name;

    EnumCommand(const CommandArgs &_args, Completer::EnumOption &_opt);

    int vm_count_check() const override;

    int help() const override;

    int desc() const override;

    int run() const override;

    // -- custom --

    const EnumeratorMap &enumerators() const;

    void print_names(std::ostream &sout, const EnumeratorMap &enumerators) const;

  private:

    mutable std::unique_ptr<EnumeratorMap> m_standard_enumerators;
    mutable const EnumeratorMap *m_enumerators;
  };
}

#endif
