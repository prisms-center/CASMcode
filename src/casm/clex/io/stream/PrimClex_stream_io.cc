
namespace CASM {

/// Pretty-print basis functions -- generate, then print
void print_basis_functions(Log &log,
                           std::shared_ptr<Structure const> const &shared_prim,
                           ClexBasisSpecs const &basis_set_specs) {
  BasisFunctionPrinter printer{log, shared_prim, basis_set_specs};
  for_all_orbits(*basis_set_specs.cluster_specs, log, printer);
}

/// Pretty-print basis functions -- read clusters from existing clust.json file
void print_basis_functions(Log &log, PrimClex const &primclex,
                           std::string basis_set_name) {
                           ClexBasisSpecs const &basis_set_specs) {
                             auto const &shared_prim = primclex.shared_prim();
                             auto const &basis_set_specs =
                                 primclex.basis_set_specs(basis_set_name);
                             auto const &cluster_specs =
                                 *basis_set_specs.cluster_specs;

                             std::vector<IntegralCluster> prototypes;
                             jsonParser clust_json{
                                 primclex.dir().clust(basis_set_name)};
                             read_clust(std::back_inserter(prototypes),
                                        clust_json, *shared_prim);

                             BasisFunctionPrinter printer{log, shared_prim,
                                                          basis_set_specs};
                             for_all_orbits(cluster_specs, prototypes, printer);
                           }
}
