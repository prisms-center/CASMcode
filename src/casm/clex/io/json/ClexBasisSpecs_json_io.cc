#include "casm/basis_set/io/json/BasisFunctionSpecs_json_io.hh"
#include "casm/basis_set/BasisFunctionSpecs.hh"
#include "casm/casm_io/json/InputParser.hh"
#include "casm/clex/io/json/ClexBasisSpecs_json_io.hh"
#include "casm/clex/ClexBasisSpecs.hh"
#include "casm/clusterography/io/json/ClusterSpecs_json_io.hh"
#include "casm/clusterography/ClusterSpecs.hh"
#include "casm/misc/ParsingDictionary.hh"

#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clusterography/ClusterSpecs_impl.hh"


namespace CASM {

  void parse(
    InputParser<ClexBasisSpecs> &parser,
    std::shared_ptr<Structure const> shared_prim,
    ParsingDictionary<DoFType::Traits> const *dof_dict) {

    auto bfuncs_subparser = parser.subparse<BasisFunctionSpecs>("basis_function_specs", *shared_prim, dof_dict);
    auto cspecs_subparser = parser.subparse<ClusterSpecs>("cluster_specs", shared_prim);
    if(parser.valid()) {
      parser.value = notstd::make_unique<ClexBasisSpecs>(
                       *bfuncs_subparser->value,
                       std::move(cspecs_subparser->value));
    }
  }

  jsonParser &to_json(
    ClexBasisSpecs const &basis_set_specs,
    jsonParser &json,
    Structure const &prim,
    ParsingDictionary<DoFType::Traits> const *dof_dict) {

    to_json(basis_set_specs.basis_function_specs, json["basis_function_specs"], prim, dof_dict);
    to_json(*basis_set_specs.cluster_specs, json["cluster_specs"]);
    return json;
  }
}
