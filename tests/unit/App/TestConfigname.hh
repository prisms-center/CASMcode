#ifndef CASM_TestConfigname
#define CASM_TestConfigname

#include "casm/clex/ConfigIO.hh"

extern "C" {
  CASM::BaseDatumFormatter<CASM::Configuration> *make_TestConfigname_formatter();
}

namespace CASM {

  namespace ConfigIO {

    ConfigIO::GenericConfigFormatter<std::string> test_configname();

  }

}

#endif
