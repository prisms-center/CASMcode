#ifndef CASM_DiffTransConfigurationTraits
#define CASM_DiffTransConfigurationTraits

#include <string>
#include <vector>
#include "casm/kinetics/DiffTransConfiguration.hh"
namespace CASM {
  template<typename T> struct traits;

  namespace Kinetics {
    class DiffTransConfiguration;
  }

  template<>
  struct traits<CASM::Kinetics::DiffTransConfiguration> {
    static const std::string name;
    static const std::string short_name;
    static bool name_compare(std::string A, std::string B);
  };
}

#endif