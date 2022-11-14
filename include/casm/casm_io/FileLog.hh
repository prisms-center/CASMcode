#ifndef CASM_FileLog
#define CASM_FileLog

#include <fstream>

#include "casm/casm_io/Log.hh"

namespace CASM {

class FileLog : public Log {
 public:
  /// \brief Construct a FileLog
  ///
  /// \param verbosity The amount to be printed
  ///
  /// For verbosity:
  /// - 0: print nothing
  /// - 10: print all standard output
  /// - 100: print all possible output
  FileLog(std::ofstream &&other, int _verbosity = standard,
          bool _show_clock = false)
      : Log(std::cout, _verbosity, _show_clock), m_fs(std::move(other)) {
    reset(m_fs);
  }

  std::ofstream &fs() { return m_fs; };

  const std::ofstream &fs() const { return m_fs; };

 private:
  std::ofstream m_fs;
};

}  // namespace CASM

#endif
