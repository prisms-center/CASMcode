#ifndef Popen_HH
#define Popen_HH

#include <iostream>
#include <string>
#include <functional>
#include <limits.h>

namespace CASM {

  /// \brief Remember how to use popen
  class Popen {

  public:

    /// \brief Construct a Popen object
    Popen(std::function<void(FILE *)> _popen_handler = Popen::default_popen_handler,
          std::function<void(int)> _pclose_handler = Popen::default_pclose_handler,
          bool _combine_stdout_stderr = true);

    /// \brief Execute popen for a given command
    void popen(std::string _command);

    /// \brief Returns the stdout resulting from the last popen call
    std::string gets() const;

    /// \brief Print the last command executed and the resulting stdout
    void print(std::ostream &sout) const;

    /// \brief Returns pclose(fp)/256
    int exit_code() const;

    /// \brief Default popen error handler throws std::runtime_error if popen failed
    static void default_popen_handler(FILE *fp);

    /// \brief Default pclose error handler throws std::runtime_error if pclose failed
    static void default_pclose_handler(int status);

  private:

    std::string m_command;
    std::string m_stdout;
    int m_pclose_result;
    std::function<void(FILE *)> m_popen_handler;
    std::function<void(int)> m_pclose_handler;
    bool m_combine_stdout_stderr;
  };

}

#endif

