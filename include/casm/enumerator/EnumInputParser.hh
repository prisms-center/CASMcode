#ifndef CASM_EnumInputParser
#define CASM_EnumInputParser

#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/app/AppIO_impl.hh"
#include "casm/clex/HasPrimClex.hh"
#include "casm/completer/Handlers.hh"
#include "casm/enumerator/Enumerator_impl.hh"
#include "casm/global/enum/stream_io.hh"
#include "casm/misc/CRTPBase.hh"

namespace CASM {

  // --- Optional, check CLI and JSON options ---

  template<typename OptHandlerType>
  int parse_verbosity(KwargsParser &parser, const OptHandlerType &opt) {
    std::string verbosity_str = opt.vm().count("verbosity") ?
                                opt.verbosity_str() :
                                parser.optional_else<std::string>("verbosity", "standard");
    auto val = Log::verbosity_level(verbosity_str);
    if(val.first) {
      return val.second;
    }
    else {
      parser.error.insert(Log::invalid_verbosity_msg(verbosity_str));
      return Log::standard;
    }
    return 0;
  }

  template<typename OptHandlerType>
  bool parse_dry_run(KwargsParser &parser, const OptHandlerType &opt) {
    return opt.vm().count("dry-run") ?
           true :
           parser.optional_else<bool>("dry_run", false);
  }

  template<typename OptHandlerType>
  COORD_TYPE parse_coord_type(KwargsParser &parser, const OptHandlerType &opt) {
    return opt.vm().count("coord") ?
           opt.coordtype_enum() :
           parser.optional_else<COORD_TYPE>(traits<COORD_TYPE>::name, COORD_TYPE::FRAC);
  }

  template<typename OptHandlerType>
  ORBIT_PRINT_MODE parse_orbit_print_mode(KwargsParser &parser, const OptHandlerType &opt) {
    return parser.optional_else<ORBIT_PRINT_MODE>(traits<ORBIT_PRINT_MODE>::name, ORBIT_PRINT_MODE::PROTO);
  }

  template<typename OptHandlerType>
  std::vector<std::string> parse_filter_expr(KwargsParser &parser, const OptHandlerType &opt) {
    if(opt.vm().count("filter")) {
      return opt.filter_strs();
    }
    else {
      auto ptr = parser.optional<std::string>("filter");
      if(ptr) {
        return std::vector<std::string>({*ptr});
      }
      else {
        return std::vector<std::string>();
      }
    }
  }

  class SymInfoOptionsParser : public InputParser<std::nullptr_t> {
  public:
    static std::string brief_help();

    static std::string print_matrix_tau_help();

    static std::string standard_help();

    template<typename CompleterOptionType>
    SymInfoOptionsParser(
      const PrimClex &_primclex,
      jsonParser &_input,
      const CompleterOptionType &_opt,
      fs::path _path,
      bool _required) :
      InputParser<std::nullptr_t>(_input, _path, _required)  {
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
      m_sym_info_opt.coord_type = parse_coord_type(*this, _opt);
      m_sym_info_opt.prec = optional_else<int>("prec", 7);

    }

    const SymInfoOptions &sym_info_opt() const;

  protected:
    static std::set<std::string> expected();

  private:
    SymInfoOptions m_sym_info_opt;

  };

  class OrbitPrinterOptionsParser : public InputParser<std::nullptr_t> {
  public:

    static std::string print_coordinates_help();

    static std::string print_equivalence_map_help();

    static std::string print_invariant_grp_help();

    static std::string brief_help();

    static std::string standard_help();

    template<typename CompleterOptionType>
    OrbitPrinterOptionsParser(
      const PrimClex &_primclex,
      jsonParser &_input,
      const CompleterOptionType &_opt,
      fs::path _path,
      bool _required) :
      InputParser<std::nullptr_t>(_input, _path, _required)  {
      _opt.desc();

      if(exists()) {
        m_orbit_printer_opt = jsonConstructor<OrbitPrinterOptions>::from_json(self);

        // check for unexpected values
        warn_unnecessary(expected());
      }
      else {
        m_orbit_printer_opt = OrbitPrinterOptions();
      }
      m_orbit_printer_opt.coord_type = parse_coord_type(*this, _opt);
      m_orbit_printer_opt.prec = optional_else<int>("prec", 7);

      // this->kwargs["sym_info_opt"] = m_sym_info_opt_parser =
      //                                  std::make_shared<SymInfoOptionsParser>(_primclex, input, _opt, relpath("sym_info_opt"), false);

      m_sym_info_opt_parser = std::make_shared<SymInfoOptionsParser>(_primclex, input, _opt, relpath("sym_info_opt"), false);
      this->insert(m_sym_info_opt_parser->path, m_sym_info_opt_parser);

      m_orbit_printer_opt.sym_info_opt = m_sym_info_opt_parser->sym_info_opt();

    }

    const OrbitPrinterOptions &orbit_printer_opt() const;

  protected:
    static std::set<std::string> expected();

  private:
    OrbitPrinterOptions m_orbit_printer_opt;
    std::shared_ptr<SymInfoOptionsParser> m_sym_info_opt_parser;
  };

  class EnumInputParser : public InputParser<std::nullptr_t>, public HasPrimClex<CRTPBase<EnumInputParser>> {

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
