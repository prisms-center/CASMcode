#ifndef CASM_symmetry_SymInfo_json_io
#define CASM_symmetry_SymInfo_json_io

#include "casm/casm_io/enum/json_io.hh"
#include "casm/global/definitions.hh"
#include "casm/symmetry/SymInfo.hh"
#include "casm/symmetry/io/stream/SymInfo_stream_io.hh"

namespace CASM {

struct SymInfo;
struct SymInfoOptions;
template <typename T>
struct jsonConstructor;
class jsonParser;

ENUM_JSON_IO_DECL(symmetry_type)

jsonParser &to_json(const SymInfoOptions &opt, jsonParser &json);

/// \brief Read from JSON
void from_json(SymInfoOptions &opt, const jsonParser &json);

template <>
struct jsonConstructor<SymInfoOptions> {
  static SymInfoOptions from_json(const jsonParser &json);
};

/// \brief Adds to existing JSON object
void to_json(const SymInfo &info, jsonParser &json);

}  // namespace CASM

#endif
