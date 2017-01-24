#include "TestConfigname.hh"
#include "casm/clex/Configuration.hh"

extern "C" {
  CASM::BaseDatumFormatter<CASM::Configuration> *make_TestConfigname_formatter() {
    return CASM::ConfigIO::test_configname().clone().release();
  }
}

namespace CASM {

  namespace ConfigIO {

    GenericConfigFormatter<std::string> test_configname() {
      auto lambda = [](const Configuration & config)->std::string {
        return config.name();
      };
      return GenericConfigFormatter<std::string>(
               "test_configname",
               "Configuration name, in the form 'SCEL#_#_#_#_#_#_#/#'",
               lambda);
    }

  }
}
