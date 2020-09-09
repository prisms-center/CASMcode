#ifndef CASM_app_enum_json_io
#define CASM_app_enum_json_io

namespace CASM {

  template<typename T> class InputParser;
  struct InsertAlgorithmsOptions;

  void parse(InputParser<InsertAlgorithmsOptions> &parser);

}

#endif
