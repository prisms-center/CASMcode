#include "casm/basis_set/io/json/BasisFunctionSpecs_json_io.hh"

#include "casm/basis_set/BasisFunctionSpecs.hh"
#include "casm/basis_set/DoFTraits.hh"
#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/InputParser.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/crystallography/Structure.hh"

namespace CASM {

ENUM_JSON_IO_DEF(PARAM_PACK_TYPE)

// TODO: JCT, check this, particularly the continuous DoF options documentation

/// Parse BasisFunctionSpecs from JSON & validate
///
/// Allowed JSON attributes:
///     dofs: array of string (optional, default= all prim dof keys)
///         An array of string of dof type names that should be used to
///         construct basis functions. The default value is all DoF types
///         included in the prim.
///     global_max_poly_order: int (optional, default=-1)
///         See `orbit_branch_max_poly_order` documentation.
///     orbit_branch_max_poly_order: object (optional, default={})
///         By default, for a given cluster orbit, polynomials of order up to
///         the cluster size are created. Higher order polynomials can be
///         requested either on a per-orbit-branch or global basis. The most
///         specific level specified is used. Orbit branches are specified using
///         the string value of the cluster size as a key.
///
///         Example:
///             "orbit_branch_max_poly_order": {
///                 "4": 7  // use maximum polynomial order == 7, for orbits of
///                 cluster size == 4
///             },
///             "global_max_poly_order": 3, // use max(3, cluster size), for all
///             other orbits
///
///     param_pack_type: string (optional, default="default")
///         Specifies the ParamPack type to use when writing the Clexulator.
///         Options are "default" or "diff", which enables fadbad automatic
///         differentiating.
///     dof_specs: object (required for some dofs)
///         Provides DoF-particular specifications for constructing basis
///         functions. Not all DoF types require their own DoFSpecs. See
///         documentation for a particular DoF type to determine if it is
///         required.
///
///         For "occ": (required if occupation dof included)
///             site_basis_functions: string or array (required)
///                 Must be one of:
///                 - "chebychev": For basis functions generated about the
///                 random alloy.
///                 - "occupation": For basis functions generated about the
///                 order alloy generated
///                   by the first occupant listed for every sublattice in the
///                   prim structure.
///                 - An array specifying sublat compositions, for "composition"
///                 basis functions
///                   generated about an average composition speficified for
///                   each sublattice.
///
///                   Example sublattice composition specification, for a prim
///                   structure with four sublattices and two allowed occupants
///                   ("A" and "B") on each sublattice:
///
///                       [
///                         { // composition on sublattices 0 and 1, as listed
///                         in prim
///                           "sublat_indices": [0, 1],
///                           "composition": {"A": 0.25, "B": 0.75}
///                         },
///                         { // composition on sublattices 2 and 3, as listed
///                         in prim
///                           "sublat_indices": [2, 3],
///                           "composition": {"A": 0.75, "B": 0.25}
///                         }
///                       ]
///
///         For "<flavor>magspin": (optional if "<flavor>magspin" dof included)
///             max_poly_order: int (optional, default=-1)
///                 Specifies the maximum polynomial order for site basis
///                 functions.
///
void parse(InputParser<BasisFunctionSpecs> &parser, Structure const &prim,
           ParsingDictionary<DoFType::Traits> const *dof_dict) {
  // read generic JSON into a BasisFunctionSpecs object
  parser.value = notstd::make_unique<BasisFunctionSpecs>();
  auto &basis_function_specs = *parser.value;

  // parse "dofs", using all dof types in prim as default value
  auto all_dof_keys = all_dof_types(prim.structure());
  parser.optional_else(basis_function_specs.dof_keys, "dofs", all_dof_keys);

  // parse "global_max_poly_order"
  parser.optional_else(basis_function_specs.global_max_poly_order,
                       "global_max_poly_order", Index{-1});

  // parse "orbit_branch_max_poly_order"
  if (parser.self.contains("orbit_branch_max_poly_order")) {
    std::stringstream err_msg;
    err_msg << "Error: could not read 'orbit_branch_max_poly_order': Expected "
               "an object with "
               "pairs of \"<branch>\":<max_poly_order>, where <branch> is the "
               "number of sites in a cluster "
               "(as a string), and <max_poly_order> is an integer.";

    auto const &json = parser.self["orbit_branch_max_poly_order"];

    try {
      for (auto it = json.begin(); it != json.end(); ++it) {
        Index branch = std::stoi(it.name());
        Index max_poly_order = it->get<Index>();
        parser.value->orbit_branch_max_poly_order[branch] = max_poly_order;
      }
    } catch (std::exception &e) {
      parser.error.insert(err_msg.str());
    }
  }

  // parse "include_functions", "exclude_functions"
  parser.optional(basis_function_specs.include_functions, "include_functions");
  parser.optional(basis_function_specs.exclude_functions, "exclude_functions");

  // parse "param_pack_type"
  parser.optional_else(basis_function_specs.param_pack_type, "param_pack_type",
                       PARAM_PACK_TYPE::DEFAULT);

  // parse & validate "dof_specs" for each DoF type here...
  for (const auto &key : all_dof_keys) {
    dof_dict->lookup(key).parse_dof_specs(parser, prim);
  }
}

jsonParser &to_json(const BasisFunctionSpecs &basis_function_specs,
                    jsonParser &json, Structure const &prim,
                    ParsingDictionary<DoFType::Traits> const *dof_dict) {
  json = jsonParser::object();
  json["dofs"] = basis_function_specs.dof_keys;
  json["global_max_poly_order"] = basis_function_specs.global_max_poly_order;
  for (auto const &pair : basis_function_specs.orbit_branch_max_poly_order) {
    json["orbit_branch_max_poly_order"][std::to_string(pair.first)] =
        pair.second;
  }
  json["param_pack_type"] = basis_function_specs.param_pack_type;
  if (basis_function_specs.include_functions.size()) {
    json["include_functions"] = basis_function_specs.include_functions;
  }
  if (basis_function_specs.exclude_functions.size()) {
    json["exclude_functions"] = basis_function_specs.exclude_functions;
  }

  // parse & validate specialized DoFType specs from JSON here...
  auto all_dof_keys = all_dof_types(prim.structure());
  for (const auto &key : all_dof_keys) {
    dof_dict->lookup(key).dof_specs_to_json(basis_function_specs, json, prim);
  }

  return json;
}

}  // namespace CASM
