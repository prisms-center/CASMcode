#ifndef CASM_sym_json_io
#define CASM_sym_json_io

namespace CASM {

namespace Completer {
class SymOption;
}
class jsonParser;

/// Convert `casm sym` CLI input to JSON
jsonParser &to_json(const Completer::SymOption &sym_opt, jsonParser &json);

}  // namespace CASM

#endif
