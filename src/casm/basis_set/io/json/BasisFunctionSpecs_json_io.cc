#include "casm/basis_set/io/json/BasisFunctionSpecs_json_io.hh"
#include "casm/basis_set/BasisFunctionSpecs.hh"
#include "casm/basis_set/DoFTraits.hh"
#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/InputParser.hh"
#include "casm/crystallography/Structure.hh"

#include "casm/casm_io/json/InputParser_impl.hh"

namespace CASM {

  ENUM_JSON_IO_DEF(PARAM_PACK_TYPE)

  // TODO: JCT, check this, particularly the continuous DoF options documentation

  /// Parse BasisFunctionSpecs from JSON & validate
  ///
  /// Allowed JSON attributes:
  ///     dofs: array of string (optional, default= all prim dof keys)
  ///         An array of string of dof type names that should be used to construct basis functions.
  ///         The default value is all DoF types included in the prim.
  ///     max_poly_order: int (optional, default=-1)
  ///         Specify maximum polynomial order to construct for cluster basis functions. By default
  ///         polynomials of order up to cluster size are created. If max_poly_order is greater than
  ///         the cluster size, higher order polynomials up to that maximum value will also be
  ///         included.
  ///     param_pack_type: string (optional, default="default")
  ///         Specifies the ParamPack type to use when writing the Clexulator. Options are "default"
  ///         or "diff", which enables fadbad automatic differentiating.
  ///     dof_specs: object (required for some dofs)
  ///         Provides DoF-particular specifications for constructing basis functions. Not all DoF
  ///         types require their own DoFSpecs. See documentation for a particular DoF type to
  ///         determine if it is required.
  ///
  ///         For "occupation": (required if "occupation" dof included)
  ///             site_basis_functions: string or array (required)
  ///                 Must be one of:
  ///                 - "chebychev": For basis functions generated about the random alloy.
  ///                 - "occupation": For basis functions generated about the order alloy generated
  ///                   by the first occupant listed for every sublattice in the prim structure.
  ///                 - An array specifying sublat compositions, for "composition" basis functions
  ///                   generated about an average composition speficified for each sublattice.
  ///                   Example sublattice composition specification, for a prim structure with four
  ///                   sublattices and two allowed occupants ("A" and "B") on each sublattice:
  ///                       [
  ///                         { // composition on sublattices 0 and 1, as listed in prim
  ///                           "sublat_indices": [0, 1],
  ///                           "composition": {"A": 0.25, "B": 0.75}
  ///                         },
  ///                         { // composition on sublattices 2 and 3, as listed in prim
  ///                           "sublat_indices": [2, 3],
  ///                           "composition": {"A": 0.75, "B": 0.25}
  ///                         }
  ///                       ]
  ///
  ///         For "<flavor>magspin": (optional if "<flavor>magspin" dof included)
  ///             max_poly_order: int (optional, default=-1)
  ///                 Specifies the maximum polynomial order for site basis functions.
  ///
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
