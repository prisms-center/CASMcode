#ifndef CASM_ConfigurationTraits
#define CASM_ConfigurationTraits

#include <string>
#include <vector>
namespace CASM {
template <typename T>
struct traits;

class Configuration;

template <>
struct traits<Configuration> {
  static const std::string name;
  static const std::string short_name;
  static bool name_compare(std::string A, std::string B);
};
}  // namespace CASM

#endif
