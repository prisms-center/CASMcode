#ifndef CASM_CLIParse
#define CASM_CLIParse

#include <wordexp.h>
#include <string>

namespace CASM {

  class CLIParse {
  public:

    /// Non-owning
    CLIParse(int _argc, char **_argv);

    /// Owning
    CLIParse(std::string _args);

    ~CLIParse();

    int argc() const {
      return m_argc;
    }

    char **argv() const {
      return m_argv;
    }

    int parse_result() const {
      return m_parse_result;
    }

  private:

    int m_argc;
    char **m_argv;

    /// stores error codes when attempting to parse std::string _args -> argc, argv
    int m_parse_result;

    bool m_free_p;
    wordexp_t m_p;
  };

  namespace Completer {
    class OptionHandlerBase;
  }

  /// Take CLI args string, 'casm X ...', and use boost::program_options to parse into Options
  void parse_args(Completer::OptionHandlerBase &opt, std::string args);
}

#endif
