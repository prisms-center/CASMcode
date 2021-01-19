#ifndef CASM_Validator
#define CASM_Validator

#include <set>
#include <string>

namespace CASM {

/// Data structure to hold error and warning messages
struct Validator {
  std::set<std::string> error;
  std::set<std::string> warning;

  void clear() {
    error.clear();
    warning.clear();
  }

  Validator &insert(const Validator &other) {
    error.insert(other.error.begin(), other.error.end());
    warning.insert(other.warning.begin(), other.warning.end());
    return *this;
  }

  bool valid() const { return !error.size(); }
};

}  // namespace CASM

#endif
