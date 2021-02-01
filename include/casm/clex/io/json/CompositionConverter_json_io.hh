#ifndef CASM_clex_CompositionConverter_json_io
#define CASM_clex_CompositionConverter_json_io

namespace CASM {

class CompositionConverter;
class jsonParser;

/// \brief Serialize CompositionConverter to JSON
jsonParser &to_json(const CompositionConverter &f, jsonParser &json);

/// \brief Deserialize CompositionConverter from JSON
void from_json(CompositionConverter &f, const jsonParser &json);

}  // namespace CASM

#endif
