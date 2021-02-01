#include "casm/clex/io/stream/ClexBasis_stream_io.hh"

#include "casm/app/DirectoryStructure.hh"
#include "casm/casm_io/Log.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/clex/ClexBasisSpecs.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clusterography/io/OrbitPrinter_impl.hh"
#include "casm/crystallography/Structure.hh"

namespace CASM {

ClexBasisFunctionPrinter::ClexBasisFunctionPrinter(
    Log &_log, std::shared_ptr<Structure const> _shared_prim,
    ClexBasisSpecs const &_basis_set_specs)
    : m_shared_prim(_shared_prim),
      m_basis_set_specs(_basis_set_specs),
      m_log(_log) {}

/// Pretty-print basis functions -- generate, then print
void print_basis_functions(Log &log,
                           std::shared_ptr<Structure const> const &shared_prim,
                           ClexBasisSpecs const &basis_set_specs) {
  ClexBasisFunctionPrinter printer{log, shared_prim, basis_set_specs};
  for_all_orbits(*basis_set_specs.cluster_specs, log, printer);
}

/// Pretty-print basis functions -- read clusters from existing clust.json file
void print_basis_functions(Log &log, PrimClex const &primclex,
                           std::string basis_set_name) {
  auto const &shared_prim = primclex.shared_prim();
  auto const &basis_set_specs = primclex.basis_set_specs(basis_set_name);
  auto const &cluster_specs = *basis_set_specs.cluster_specs;

  std::vector<IntegralCluster> prototypes;
  jsonParser clust_json{primclex.dir().clust(basis_set_name)};
  read_clust(std::back_inserter(prototypes), clust_json, *shared_prim);

  ClexBasisFunctionPrinter printer{log, shared_prim, basis_set_specs};
  for_all_orbits(cluster_specs, prototypes, printer);
}

}  // namespace CASM
