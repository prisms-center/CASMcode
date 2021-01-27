#ifndef CASM_ClexBasisSpecs_json_io
#define CASM_ClexBasisSpecs_json_io

#include <utility>

namespace CASM {

namespace DoFType {
class Traits;
}

template <typename T>
class InputParser;
template <typename T>
class ParsingDictionary;

struct ClexBasisSpecs;
class Structure;
class jsonParser;

/// Parse ClexBasisSpecs from JSON & validate (bspecs.json input)
///
/// Expects:
/// \code
/// {
///   "basis_function_specs": <BasisFunctionSpecs JSON>,
///   "cluster_specs": <ClusterSpecs JSON>
/// }
/// \endcode
///
/// Note:
/// - For BasisFunctionSpecs JSON format, see
/// src/casm/basis_set/io/json/BasisFunctionSpecs_json_io.cc
/// - For ClusterSpecs JSON format, see
/// src/casm/clusterography/io/json/ClusterSpecs_json_io.cc
///
void parse(InputParser<ClexBasisSpecs> &parser,
           std::shared_ptr<Structure const> shared_prim,
           ParsingDictionary<DoFType::Traits> const *dof_dict);

/// Write ClexBasisSpecs to JSON (bspecs.json output)
///
/// Example:
/// \code
/// {
///    "basis_function_specs": <BasisFunctionSpecs JSON>,
///    "cluster_specs": <ClusterSpecs JSON>
/// }
/// \endcode
///
jsonParser &to_json(ClexBasisSpecs const &basis_set_specs, jsonParser &json,
                    Structure const &prim,
                    ParsingDictionary<DoFType::Traits> const *dof_dict);

}  // namespace CASM

#endif
