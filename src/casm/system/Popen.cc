#include "casm/system/Popen.hh"
#include <stdexcept>

namespace CASM {

  /// \brief Construct a Popen object
  Popen::Popen(
    std::function<void(FILE *)> _popen_handler,
    std::function<void(int)> _pclose_handler,
    bool _combine_stdout_stderr) :
    m_popen_handler(_popen_handler),
    m_pclose_handler(_pclose_handler),
    m_combine_stdout_stderr(_combine_stdout_stderr) {}

  /// \brief Execute popen for a given command
  ///
  /// - Uses provide error handlers
  /// - Stores popen stdout, which can be accessed using gets
  void Popen::popen(std::string _command) {

    m_command = _command;

    if(m_combine_stdout_stderr) {
      m_command += " 2>&1";
    }

    FILE *fp;
    char path[PATH_MAX];

    fp = ::popen(m_command.c_str(), "r");
    m_popen_handler(fp);

    m_stdout = "";
    while(fgets(path, PATH_MAX, fp) != nullptr) {
      m_stdout += path;
    }

    m_pclose_result = pclose(fp);

    m_pclose_handler(m_pclose_result);
  }

  /// \brief Returns the stdout resulting from the last popen call
  std::string Popen::gets() const {
    return m_stdout;
  }

  /// \brief Print the last command executed and the resulting stdout
  void Popen::print(std::ostream &sout) const {
    sout << m_command << "\n";
    sout << m_stdout;
  }

  /// \brief Returns pclose(fp)/256
  int Popen::exit_code() const {
    return m_pclose_result / 256;
  }

  /// \brief Default popen error handler throws std::runtime_error if popen failed
  void Popen::default_popen_handler(FILE *fp) {
    if(fp == nullptr) {
      throw std::runtime_error("Error: popen failed.");
    }
  }

  /// \brief Default pclose error handler throws std::runtime_error if pclose failed
  void Popen::default_pclose_handler(int status) {
    if(status == -1) {
      throw std::runtime_error("Error: pclose failed.");
    }
  }

}

