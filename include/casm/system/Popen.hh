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
    ///
    /// Allows custom popen and pclose error handlers
    Popen(std::function<void(FILE *)> _popen_handler = Popen::default_popen_handler,
          std::function<void(int)> _pclose_handler = Popen::default_pclose_handler) :
      m_popen_handler(_popen_handler),
      m_pclose_handler(_pclose_handler) {}

    /// \brief Execute popen for a given command
    ///
    /// - Uses provide error handlers
    /// - Stores popen stdout, which can be accessed using gets
    void popen(std::string _command) {

      m_command = _command;

      FILE *fp;
      char path[PATH_MAX];

      fp = ::popen(m_command.c_str(), "r");
      m_popen_handler(fp);

      m_stdout = "";
      while(fgets(path, PATH_MAX, fp) != NULL) {
        m_stdout += path;
      }

      m_pclose_handler(pclose(fp));
    }

    /// \brief Returns the stdout resulting from the last popen call
    std::string gets() const {
      return m_stdout;
    }

    /// \brief Print the last command executed and the resulting stdout
    void print(std::ostream &sout) const {
      sout << m_command << "\n";
      sout << m_stdout;
    }

    /// \brief Default popen error handler throws std::runtime_error if popen failed
    static void default_popen_handler(FILE *fp) {
      if(fp == NULL) {
        throw std::runtime_error("Error: popen failed.");
      }
    }

    /// \brief Default pclose error handler throws std::runtime_error if pclose failed
    static void default_pclose_handler(int status) {
      if(status == -1) {
        throw std::runtime_error("Error: pclose failed.");
      }
    }

  private:

    std::string m_command;
    std::string m_stdout;
    std::function<void(FILE *)> m_popen_handler;
    std::function<void(int)> m_pclose_handler;
  };

}

#endif

