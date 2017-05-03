#ifndef CASM_ConfigurationTraits
#define CASM_ConfigurationTraits

#include <string>

namespace CASM {
  template<typename T> struct traits;

  class Configuration;

  template<>
  struct traits<Configuration> {
    static const std::string name;
    static const std::string short_name;
  };
}

#endif
