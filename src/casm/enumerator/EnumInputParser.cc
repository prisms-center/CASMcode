#include "casm/enumerator/EnumInputParser.hh"

namespace CASM {

  std::string EnumInputParser::standard_help() {
    return InputParser::dry_run_help()
           + InputParser::verbosity_help()
           + InputParser::coordinate_mode_help();
  }


  EnumInputParser::EnumInputParser(
    const PrimClex &_primclex,
    jsonParser &_input,
    const Completer::EnumOption &_enum_opt,
    fs::path _path,
    bool _required) :
    InputParser(_input, _path, _required),
    m_primclex(_primclex),
    m_enum_opt(_enum_opt),
    m_dry_run(false),
    m_coord_type(FRAC),
    m_verbosity(Log::standard) {

    m_enum_opt.desc();

    if(exists()) {

      m_dry_run = parse_dry_run(m_enum_opt);
      m_coord_type = parse_coord_type(m_enum_opt);
      m_verbosity = parse_verbosity(m_enum_opt);
      m_filter_expr = parse_filter_expr(m_enum_opt);
      // warn_unnecessary should be done in derived
    }
  }

  bool EnumInputParser::dry_run() const {
    return m_dry_run;
  }

  std::string EnumInputParser::dry_run_msg() const {
    return CASM::dry_run_msg(dry_run());
  }

  COORD_TYPE EnumInputParser::coord_type() const {
    return m_coord_type;
  }

  int EnumInputParser::verbosity() const {
    return m_verbosity;
  }

  std::vector<std::string> EnumInputParser::filter_expr() const {
    return m_filter_expr;
  }

  std::set<std::string> EnumInputParser::expected() {
    return std::set<std::string>({"dry_run", "coordinate_mode", "verbosity", "filter"});
  }

  const PrimClex &EnumInputParser::primclex() const {
    return m_primclex;
  }

}
