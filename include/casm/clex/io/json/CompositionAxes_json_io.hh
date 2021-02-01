#ifndef CASM_clex_CompositionAxes_json_io
#define CASM_clex_CompositionAxes_json_io

namespace CASM {

struct CompositionAxes;
class jsonParser;

void from_json(CompositionAxes &composition_axes, const jsonParser &json);

jsonParser &to_json(CompositionAxes const &composition_axes, jsonParser &json);

}  // namespace CASM

#endif
