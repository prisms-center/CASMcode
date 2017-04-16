#include "casm/app/DirectoryStructure.hh"
#include "casm/clex/Configuration.hh"

namespace CASM {

  const std::set<std::string> &config_types() {
    static std::set<std::string> _config_types {
      QueryTraits<Configuration>::name
    };
    return _config_types;
  };

  const std::set<std::string> &config_types_short() {
    static std::set<std::string> _config_types_short {
      QueryTraits<Configuration>::short_name
    };
    return _config_types_short;
  };

}
