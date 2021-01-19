#ifndef CASM_ECIContainer_json_io
#define CASM_ECIContainer_json_io

namespace CASM {

template <typename T>
class InputParser;
class ECIContainer;

///  Make ECIContainer from JSON (eci.json file)
void parse(InputParser<ECIContainer> &parser);
}  // namespace CASM

#endif
