#include "casm/basis_set/io/json/BasisFunctionSpecs_json_io.hh"
#include "casm/basis_set/BasisFunctionSpecs.hh"
#include "casm/basis_set/DoFTraits.hh"
#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/InputParser.hh"
#include "casm/crystallography/Structure.hh"

#include "casm/casm_io/json/InputParser_impl.hh"

namespace CASM {

  ENUM_JSON_IO_DEF(PARAM_PACK_TYPE)

  /// \brief Parse BasisFunctionSpecs from JSON & validate
  ///
  /// Allowed JSON attributes:
  ///     dofs: array of string (optional, default= all prim dof keys)
  ///         An array of string of dof type names that should be used to construct basis functions.
  ///         The default value is all DoF types included in the prim.
  ///     max_poly_order: int (optional, default=-1)
  ///         Specify maximum polynomial order to construct. Default value of -1 specifies context-
  ///         dependent default behaviour.
  ///     param_pack_type: string (optional, default="DEFAULT")
  ///         Specifies the ParamPack type to use when writing the Clexulator. Options are "DEFAULT"
  ///         or "DIFF", which uses fadbad automatic differentiating.
  ///     site_basis_functions: string or array (required if "occ" dof included)
  ///         Must be one of "chebychev" or "occupation" or an array specifying sublat compositions.
  ///         Example sublat composition specification:
  ///             [
  ///               {
  ///                 "sublat_indices": [0, 1],
  ///                 "composition": {"A": 0.25, "B": 0.75}
  ///               },
  ///             {
  ///               "sublat_indices": [2, 3],
  ///               "composition": {"A": 0.75, "B": 0.25}
  ///             }
  ///           ]
  /// \endcode
  void parse(
    InputParser<BasisFunctionSpecs> &parser,
    Structure const &prim,
    ParsingDictionary<DoFType::Traits> const *dof_dict) {

    // get all_dof_keys, to use as default value for "dofs", and for parsing DoFSpecs
    auto all_dof_keys = all_dof_types(prim);

    // read generic JSON into a BasisFunctionSpecs object
    parser.value = std::unique_ptr<BasisFunctionSpecs>();
    auto &basis_function_specs = *parser.value;
    parser.optional_else(basis_function_specs.dof_keys, "dofs", all_dof_keys);
    parser.optional_else(basis_function_specs.max_poly_order, "max_poly_order", Index {-1});
    parser.optional_else(basis_function_specs.parampack_type, "parampack_type", PARAM_PACK_TYPE::DEFAULT);

    // parse & validate specialized DoFType specs from JSON here...
    for(const auto &key : all_dof_keys) {
      dof_dict->lookup(key).parse_dof_specs(parser, prim);
    }
  }

  jsonParser &to_json(
    const BasisFunctionSpecs &basis_function_specs,
    jsonParser &json,
    Structure const &prim,
    ParsingDictionary<DoFType::Traits> const *dof_dict) {

    json = jsonParser::object();
    json["dofs"] = basis_function_specs.dof_keys;
    json["max_poly_order"] = basis_function_specs.max_poly_order;
    json["parampack_type"] = basis_function_specs.parampack_type;

    // parse & validate specialized DoFType specs from JSON here...
    auto all_dof_keys = all_dof_types(prim);
    for(const auto &key : all_dof_keys) {
      dof_dict->lookup(key).dof_specs_to_json(basis_function_specs, json, prim);
    }

    return json;
  }

}
