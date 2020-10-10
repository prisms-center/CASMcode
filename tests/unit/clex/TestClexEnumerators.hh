#ifndef CASMtest_TestClexEnumerators
#define CASMtest_TestClexEnumerators

#include "casm/app/CLIParse.hh"
#include "casm/app/enum/EnumInterface.hh"
#include "casm/app/io/json_io_impl.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/completer/Handlers.hh"

namespace test {

  template<typename EnumInterfaceType>
  void run_enum_interface(std::string cli_str,
                          CASM::PrimClex &primclex,
                          jsonParser const &json_options = jsonParser()) {
    CASM::Completer::EnumOption opt;
    parse_args(opt, cli_str);
    CASM::jsonParser cli_options_as_json {opt};

    EnumInterfaceType interface;
    interface.run(primclex, json_options, cli_options_as_json);
  }
}

#endif
