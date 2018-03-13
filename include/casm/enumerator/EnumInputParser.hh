#ifndef CASM_EnumInputParser
#define CASM_EnumInputParser

#include "casm/casm_io/InputParser_impl.hh"
#include "casm/app/AppIO_impl.hh"
#include "casm/clex/HasPrimClex.hh"
#include "casm/completer/Handlers.hh"
#include "casm/enumerator/Enumerator_impl.hh"
#include "casm/misc/CRTPBase.hh"

namespace CASM {

  class SymInfoOptionsParser : public InputParser {
  public:
    static std::string brief_help() {
      return
        "  sym_info_opt: JSON object (optional, see 'casm format --sym-info')\n"
        "    Options controlling the printing of symmetry operations.\n\n";
    }

    static std::string print_matrix_tau_help() {
      return
        "  print_matrix_tau: bool (optional, default=false)\n"
        "    Print the symmetry operations that map the prototype to each orbit element.\n\n";
    }

    static std::string standard_help() {
      return InputParser::coordinate_mode_help()
             + InputParser::prec_help("printing coordinates of symmetry operations", 7)
             + print_matrix_tau_help();
    }

    template<typename CompleterOptionType>
    SymInfoOptionsParser(
      const PrimClex &_primclex,
      jsonParser &_input,
      const CompleterOptionType &_opt,
      fs::path _path,
      bool _required) :
      InputParser(_input, _path, _required)  {
      _opt.desc();

      if(exists()) {
        m_sym_info_opt = jsonConstructor<SymInfoOptions>::from_json(self);

        // check for unexpected values
        warn_unnecessary(expected());
      }
      else {
        m_sym_info_opt = SymInfoOptions();
      }
      m_sym_info_opt.tol = _primclex.crystallography_tol();
      m_sym_info_opt.coord_type = parse_coord_type(_opt);
      m_sym_info_opt.prec = optional_else<int>("prec", 7);
    }

    const SymInfoOptions &sym_info_opt() const {
      return m_sym_info_opt;
    }

  protected:
    static std::set<std::string> expected() {
      return std::set<std::string>({"tol", "coordinate_mode", "prec", "print_matrix_tau"});
    }

  private:
    SymInfoOptions m_sym_info_opt;

  };

  class OrbitPrinterOptionsParser : public InputParser {
  public:

    static std::string print_coordinates_help() {
      return
        "  print_coordinates: bool (optional, default=true)\n"
        "    Print coordinates of orbit elements.\n\n";
    }

    static std::string print_equivalence_map_help() {
      return
        "  print_equivalence_map: bool (optional, default=false)\n"
        "    Print the symmetry operations that map the prototype to each orbit element.\n\n";
    }

    static std::string print_invariant_grp_help() {
      return
        "  print_invariant_grp: bool (optional, default=false)\n"
        "    Print the symmetry operations that leave each orbit element invariant.\n\n";
    }

    static std::string brief_help() {
      return
        "  orbit_printer_opt: JSON object (optional, see 'casm format --orbit-printer')\n"
        "    Options controlling the printing of orbits.\n\n";
    }

    static std::string standard_help() {
      return InputParser::indent_space_help()
             + InputParser::prec_help("printing coordinates of orbit elements", 7)
             + InputParser::coordinate_mode_help()
             + InputParser::orbit_print_mode_help()
             + SymInfoOptionsParser::brief_help()
             + print_coordinates_help()
             + print_equivalence_map_help()
             + print_invariant_grp_help();
    }

    template<typename CompleterOptionType>
    OrbitPrinterOptionsParser(
      const PrimClex &_primclex,
      jsonParser &_input,
      const CompleterOptionType &_opt,
      fs::path _path,
      bool _required) :
      InputParser(_input, _path, _required)  {
      _opt.desc();

      if(exists()) {
        m_orbit_printer_opt = jsonConstructor<OrbitPrinterOptions>::from_json(self);

        // check for unexpected values
        warn_unnecessary(expected());
      }
      else {
        m_orbit_printer_opt = OrbitPrinterOptions();
      }
      m_orbit_printer_opt.coord_type = parse_coord_type(_opt);
      m_orbit_printer_opt.prec = optional_else<int>("prec", 7);

      this->kwargs["sym_info_opt"] = m_sym_info_opt_parser =
                                       std::make_shared<SymInfoOptionsParser>(_primclex, input, _opt, relpath("sym_info_opt"), false);
      m_orbit_printer_opt.sym_info_opt = m_sym_info_opt_parser->sym_info_opt();
    }

    const OrbitPrinterOptions &orbit_printer_opt() const {
      return m_orbit_printer_opt;
    }

  protected:
    static std::set<std::string> expected() {
      return std::set<std::string>({"indent_space", "prec", traits<COORD_TYPE>::name,
                                    traits<ORBIT_PRINT_MODE>::name, "print_coordinates", "print_equivalence_map",
                                    "print_invariant_grp", "sym_info_opt"
                                   });
    }

  private:
    OrbitPrinterOptions m_orbit_printer_opt;
    std::shared_ptr<SymInfoOptionsParser> m_sym_info_opt_parser;
  };

  class EnumInputParser : public InputParser, public HasPrimClex<CRTPBase<EnumInputParser>> {

  public:

    static std::string standard_help();

    EnumInputParser(
      const PrimClex &_primclex,
      jsonParser &_input,
      const Completer::EnumOption &_enum_opt,
      fs::path _path,
      bool _required);

    bool dry_run() const;

    std::string dry_run_msg() const;

    COORD_TYPE coord_type() const;

    int verbosity() const;

    std::vector<std::string> filter_expr() const;

    const PrimClex &primclex() const;

  protected:
    static std::set<std::string> expected();

  private:
    const PrimClex &m_primclex;
    Completer::EnumOption m_enum_opt;
    bool m_dry_run;
    COORD_TYPE m_coord_type;
    int m_verbosity;
    std::vector<std::string> m_filter_expr;
  };
}

#include "casm/clex/HasPrimClex_impl.hh"

#endif
