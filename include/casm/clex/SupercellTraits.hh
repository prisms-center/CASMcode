#ifndef CASM_SupercellTraits
#define CASM_SupercellTraits

#include <string>

namespace CASM {
  template<typename T> struct traits;

  class Supercell;

  template<>
  struct traits<Supercell> {
    static const std::string name;
    static const std::string short_name;
    static bool name_compare(std::string A, std::string B);
  };
}

#endif
