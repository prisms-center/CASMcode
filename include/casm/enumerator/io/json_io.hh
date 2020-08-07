#ifndef CASM_enumerator_json_io
#define CASM_enumerator_json_io

#include <vector>

namespace CASM {

  template<typename T> class InputParser;
  struct InsertAlgorithmsOptions;

  void parse(InputParser<InsertAlgorithmsOptions> &parser);

}

#endif
