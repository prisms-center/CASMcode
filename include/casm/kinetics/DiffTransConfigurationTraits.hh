#ifndef CASM_DiffTransConfigurationTraits
#define CASM_DiffTransConfigurationTraits

#include <string>

//Traits Classes are used for conventional naming schemes i.e. how to refer to objects through selection/queries
namespace CASM {
  template<typename T> struct traits;
  template<typename ConfigType, typename IsEqualImpl> class GenericConfigCompare;

  namespace Kinetics {
    class DiffTransConfigIsEqualFast;
    class DiffTransConfigIsEqualSimple;
    class DiffTransConfiguration;

    using DiffTransConfigIsEqual = DiffTransConfigIsEqualSimple;
    using DiffTransConfigCompare = GenericConfigCompare<DiffTransConfiguration, DiffTransConfigIsEqual>;
  }

  template<>
  struct traits<CASM::Kinetics::DiffTransConfiguration> {
    static const std::string name;
    static const std::string short_name;
    static bool name_compare(std::string A, std::string B);
  };
}

#endif
