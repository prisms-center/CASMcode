#ifndef CASM_global_errors
#define CASM_global_errors

#include <string>

/// \brief Main CASM namespace
namespace CASM {

class libcasm_runtime_error : public std::runtime_error {
 public:
  libcasm_runtime_error(std::string _what) : std::runtime_error(_what) {}

  virtual ~libcasm_runtime_error() {}

 private:
  std::string m_what;
};

}  // namespace CASM

#endif
