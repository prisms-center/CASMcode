#include "casm/app/CLIParse.hh"
#include <boost/program_options.hpp>
#include "casm/completer/Handlers.hh"

namespace CASM {

  /// Non-owning
  CLIParse::CLIParse(int _argc, char **_argv, const Logging &_logging) :
    Logging(_logging),
    m_argc(_argc),
    m_argv(_argv),
    m_parse_result(0),
    m_free_p(false) {
  }

  /// Owning
  CLIParse::CLIParse(std::string _args, const Logging &_logging) :
    Logging(_logging),
    m_free_p(false) {

    // parse _args -> argc, argv
    m_parse_result = wordexp(_args.c_str(), &m_p, 0);
    if(m_parse_result) {
      err_log() << "Error parsing query: '" << _args << "'" << std::endl;
      err_log() << "wordexp() error: " << m_parse_result << std::endl;
      switch(m_parse_result) {
      case 1: {
        err_log() << "Check for illegal unescaped characters: |, &, ;, <, >, (, ), {, }" << std::endl;
        break;
      }
      default: {
        err_log() << "Check 'man wordexp' for error code meaning" << std::endl;
      }
      }
      return;
    }

    m_free_p = true;
    m_argc = m_p.we_wordc;
    m_argv = m_p.we_wordv;
  }

  CLIParse::~CLIParse() {
    if(m_free_p) {
      wordfree(&m_p);
    }
  }

  /// Take CLI args string, 'casm X ...', and use boost::program_options to parse into Options
  void parse_args(Completer::OptionHandlerBase &opt, std::string args, const Logging &_logging) {
    CLIParse parse(args, _logging);
    po::store(po::parse_command_line(parse.argc(), parse.argv(), opt.desc()), opt.vm());
    // gets default values
    po::notify(opt.vm());
  }
}
