#ifndef CASM_EnumInputParser
#define CASM_EnumInputParser

#include "casm/app/AppIO_impl.hh"
#include "casm/casm_io/InputParser_impl.hh"
#include "casm/clex/HasPrimClex.hh"
#include "casm/completer/Handlers.hh"
#include "casm/enumerator/Enumerator_impl.hh"
#include "casm/misc/CRTPBase.hh"

namespace CASM {
  class EnumInputParser : public InputParser, public HasPrimClex<CRTPBase<EnumInputParser>> {

  public:

    static const std::string standard_help;

    EnumInputParser(
      const PrimClex &_primclex,
      jsonParser &_input,
      const Completer::EnumOption &_enum_opt,
      fs::path _path,
      bool _required);

    bool dry_run() const;

    std::string dry_run_msg() const;

    COORD_TYPE coord_mode() const;

    ORBIT_PRINT_MODE orbit_print_mode() const;

    int verbosity() const;

    std::vector<std::string> filter_expr() const;

    const PrimClex &primclex() const;

  protected:
    static std::set<std::string> expected();

  private:
    const PrimClex &m_primclex;
    Completer::EnumOption m_enum_opt;
    bool m_dry_run;
    COORD_TYPE m_coord_mode;
    ORBIT_PRINT_MODE m_orbit_print_mode;
    int m_verbosity;
    std::vector<std::string> m_filter_expr;

  };
}

#endif
