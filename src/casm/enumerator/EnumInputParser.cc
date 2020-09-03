#include "casm/completer/help.hh"
#include "casm/enumerator/EnumInputParser.hh"

namespace CASM {

  std::string SymInfoOptionsParser::brief_help() {
    return
      "  sym_info_opt: JSON object (optional, see 'casm format --sym-info')\n"
      "    Options controlling the printing of symmetry operations.\n\n";
  }

  std::string SymInfoOptionsParser::print_matrix_tau_help() {
    return
      "  print_matrix_tau: bool (optional, default=false)\n"
      "    Print the symmetry operations that map the prototype to each orbit element.\n\n";
  }

  std::string SymInfoOptionsParser::standard_help() {
    return coordinate_mode_help()
           + prec_help("printing coordinates of symmetry operations", 7)
           + print_matrix_tau_help();
  }

  const SymInfoOptions &SymInfoOptionsParser::sym_info_opt() const {
    return m_sym_info_opt;
  }

  std::set<std::string> SymInfoOptionsParser::expected() {
    return std::set<std::string>({"tol", "coordinate_mode", "prec", "print_matrix_tau"});
  }


  std::string OrbitPrinterOptionsParser::print_coordinates_help() {
    return
      "  print_coordinates: bool (optional, default=true)\n"
      "    Print coordinates of orbit elements.\n\n";
  }

  std::string OrbitPrinterOptionsParser::print_equivalence_map_help() {
    return
      "  print_equivalence_map: bool (optional, default=false)\n"
      "    Print the symmetry operations that map the prototype to each orbit element.\n\n";
  }

  std::string OrbitPrinterOptionsParser::print_invariant_group_help() {
    return
      "  print_invariant_group: bool (optional, default=false)\n"
      "    Print the symmetry operations that leave each orbit element invariant.\n\n";
  }

  std::string OrbitPrinterOptionsParser::brief_help() {
    return
      "  orbit_printer_opt: JSON object (optional, see 'casm format --orbit-printer')\n"
      "    Options controlling the printing of orbits.\n\n";
  }

  std::string OrbitPrinterOptionsParser::standard_help() {
    return indent_space_help()
           + prec_help("printing coordinates of orbit elements", 7)
           + coordinate_mode_help()
           + orbit_print_mode_help()
           + SymInfoOptionsParser::brief_help()
           + print_coordinates_help()
           + print_equivalence_map_help()
           + print_invariant_group_help();
  }

  const OrbitPrinterOptions &OrbitPrinterOptionsParser::orbit_printer_opt() const {
    return m_orbit_printer_opt;
  }

  std::set<std::string> OrbitPrinterOptionsParser::expected() {
    return std::set<std::string>({"indent_space", "prec", traits<COORD_TYPE>::name,
                                  traits<ORBIT_PRINT_MODE>::name, "print_coordinates", "print_equivalence_map",
                                  "print_invariant_group", "sym_info_opt"
                                 });
  }


  std::string EnumInputParser::standard_help() {
    return dry_run_help()
           + verbosity_help()
           + coordinate_mode_help();
  }


  EnumInputParser::EnumInputParser(
    const PrimClex &_primclex,
    jsonParser &_input,
    const Completer::EnumOption &_enum_opt,
    fs::path _path,
    bool _required) :
    InputParser<std::nullptr_t>(_input, _path, _required),
    m_primclex(_primclex),
    m_enum_opt(_enum_opt),
    m_dry_run(false),
    m_coord_type(FRAC),
    m_verbosity(Log::standard) {

    m_enum_opt.desc();

    if(exists()) {

      m_dry_run = parse_dry_run(*this, m_enum_opt);
      m_coord_type = parse_coord_type(*this, m_enum_opt);
      m_verbosity = parse_verbosity(*this, m_enum_opt);
      m_filter_expr = parse_filter_expr(*this, m_enum_opt);
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
