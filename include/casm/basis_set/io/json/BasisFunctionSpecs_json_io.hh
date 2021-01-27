#ifndef CASM_BasisFunctionSpecs_json_io
#define CASM_BasisFunctionSpecs_json_io

#include "casm/basis_set/io/io_traits.hh"
#include "casm/casm_io/enum/io_traits.hh"
#include "casm/casm_io/enum/json_io.hh"
#include "casm/casm_io/json/InputParser.hh"
#include "casm/misc/ParsingDictionary.hh"

namespace CASM {

namespace DoFType {
class Traits;
}
struct BasisFunctionSpecs;
class Structure;
enum class PARAM_PACK_TYPE;

ENUM_JSON_IO_DECL(PARAM_PACK_TYPE)

// see source file for JSON format documentation
void parse(InputParser<BasisFunctionSpecs> &parser, Structure const &prim,
           ParsingDictionary<DoFType::Traits> const *dof_dict);

jsonParser &to_json(const BasisFunctionSpecs &basis_function_specs,
                    jsonParser &json, Structure const &prim,
                    ParsingDictionary<DoFType::Traits> const *dof_dict);

}  // namespace CASM

#endif
