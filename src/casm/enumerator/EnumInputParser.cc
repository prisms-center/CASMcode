#include "casm/enumerator/EnumInputParser.hh"

namespace CASM {

  const std::string EnumInputParser::standard_help =
    InputParser::dry_run_help
    + InputParser::coordinate_mode_help
    + InputParser::orbit_print_mode_help
    + InputParser::verbosity_help;


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
    m_coord_mode(FRAC),
    m_orbit_print_mode(ORBIT_PRINT_MODE::PROTO),
    m_verbosity(Log::standard) {

    if(exists()) {
      m_dry_run = parse_dry_run(m_enum_opt);
      m_coord_mode = parse_coord_type(m_enum_opt);
      m_orbit_print_mode = parse_orbit_print_mode(m_enum_opt);
      m_verbosity = parse_verbosity(m_enum_opt);
      // warn_unnecessary should be done in derived
    }
  }

  bool EnumInputParser::dry_run() const {
    return m_dry_run;
  }

  std::string EnumInputParser::dry_run_msg() const {
    return CASM::dry_run_msg(dry_run());
  }

  COORD_TYPE EnumInputParser::coordinate_mode() const {
    return m_coord_mode;
  }

  ORBIT_PRINT_MODE EnumInputParser::orbit_print_mode() const {
    return m_orbit_print_mode;
  }

  int EnumInputParser::verbosity() const {
    return m_verbosity;
  }

  std::set<std::string> expected() const {
    return std::set<std::string>({"dry_run", "coordinate_mode", "verbosity", "orbit_print_mode"});
  }

}
