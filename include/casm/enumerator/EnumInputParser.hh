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

  class SymInfoOptionsParser : public InputParser {
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

      std::cout << "end SymInfoOptionsParser: " << _path << std::endl;
      std::cout << "opt.vm().count(\"coord\"): " << _opt.vm().count("coord") << std::endl;
      if(_opt.vm().count("coord")) {
        std::cout << "_opt.coordtype_enum(): " << _opt.coordtype_enum() << std::endl;
      }
      std::cout << "optinal_else: " << optional_else<COORD_TYPE>(traits<COORD_TYPE>::name, COORD_TYPE::FRAC) << std::endl;
      std::cout << "m_sym_info_opt.coord_type: " << to_string(m_sym_info_opt.coord_type) << std::endl;
    }

    const SymInfoOptions &sym_info_opt() const;

  protected:
    static std::set<std::string> expected();

  private:
    SymInfoOptions m_sym_info_opt;

  };

  class OrbitPrinterOptionsParser : public InputParser {
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

      std::cout << "OrbitPrinterOptionsParser, before parsing sym_info_opt: " << _path << std::endl;
      std::cout << "opt.vm().count(\"coord\"): " << _opt.vm().count("coord") << std::endl;
      if(_opt.vm().count("coord")) {
        std::cout << "_opt.coordtype_enum(): " << _opt.coordtype_enum() << std::endl;
      }
      std::cout << "optinal_else: " << optional_else<COORD_TYPE>(traits<COORD_TYPE>::name, COORD_TYPE::FRAC) << std::endl;
      std::cout << "orbit_printer_opt().coord_type: " << to_string(orbit_printer_opt().coord_type) << std::endl;
      std::cout << "orbit_printer_opt().sym_info_opt.coord_type: " << to_string(orbit_printer_opt().sym_info_opt.coord_type) << std::endl;

      this->kwargs["sym_info_opt"] = m_sym_info_opt_parser =
                                       std::make_shared<SymInfoOptionsParser>(_primclex, input, _opt, relpath("sym_info_opt"), false);
      m_orbit_printer_opt.sym_info_opt = m_sym_info_opt_parser->sym_info_opt();

      std::cout << "end OrbitPrinterOptionsParser: " << _path << std::endl;
      std::cout << "orbit_printer_opt().coord_type: " << to_string(orbit_printer_opt().coord_type) << std::endl;
      std::cout << "orbit_printer_opt().sym_info_opt.coord_type: " << to_string(orbit_printer_opt().sym_info_opt.coord_type) << std::endl;
    }

    const OrbitPrinterOptions &orbit_printer_opt() const;

  protected:
    static std::set<std::string> expected();

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
