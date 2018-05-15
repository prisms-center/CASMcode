#include "casm/clex/OrbitFunctionTraits.hh"
#include "casm/basis_set/BasisSet.hh"
#include "casm/clusterography/IntegralCluster_impl.hh"

namespace CASM {
  BasisSet InvariantPolyBasisBuilder::build_proto(IntegralCluster const &_prototype,
                                                  SymGroup const &_generating_group,
                                                  std::vector<BasisSet const *> const &_arg_bases,
                                                  Index max_poly_order,
                                                  Index min_poly_order) const {

    BasisSet res(name());
    max_poly_order = valid_index(max_poly_order) ? max_poly_order : _prototype.size();

    Array<BasisSet const *> arg_bases_array(_arg_bases.begin(), _arg_bases.end());
    for(Index i = min_poly_order; i <= max_poly_order; i++) {
      BasisSet tres(name());
      tres.construct_invariant_polynomials(arg_bases_array, _generating_group, i, 1);
      res.append(tres);
    }

    return res;
  }


}
