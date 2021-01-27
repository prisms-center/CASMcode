#ifndef CASMtest_TestEnumeratorInterface
#define CASMtest_TestEnumeratorInterface

#include "casm/app/CLIParse.hh"
#include "casm/app/enum/EnumInterface.hh"
#include "casm/app/io/json_io_impl.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/completer/Handlers.hh"

namespace test {

// read `json_options` from --input string or --settings file
template <typename EnumInterfaceType>
void run_enum_interface(std::string cli_str, CASM::PrimClex &primclex) {
  CASM::Completer::EnumOption opt;
  CASM::parse_args(opt, cli_str);
  CASM::jsonParser json_options = CASM::make_json_input(opt);
  CASM::jsonParser cli_options_as_json{opt};

  EnumInterfaceType interface;
  interface.run(primclex, json_options, cli_options_as_json);
}

// pass in `json_options` and ignore --input string or --settings file
template <typename EnumInterfaceType>
void run_enum_interface(std::string cli_str, CASM::PrimClex &primclex,
                        CASM::jsonParser const &json_options) {
  CASM::Completer::EnumOption opt;
  CASM::parse_args(opt, cli_str);
  CASM::jsonParser cli_options_as_json{opt};

  EnumInterfaceType interface;
  interface.run(primclex, json_options, cli_options_as_json);
}
}  // namespace test

#endif
