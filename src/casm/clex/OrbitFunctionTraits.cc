#include "casm/clex/OrbitFunctionTraits.hh"
#include "casm/basis_set/BasisSet.hh"
#include "casm/clusterography/IntegralCluster_impl.hh"
#include "casm/casm_io/stream_io/container.hh"
#include <ostream>

namespace CASM {
  BasisSet InvariantPolyBasisBuilder::build_proto(IntegralCluster const &_prototype,
                                                  SymGroup const &_generating_group,
                                                  std::vector<BasisSet const *> const &_arg_bases,
                                                  Index max_poly_order,
                                                  Index min_poly_order) const {
    std::cout << "Begin build_proto " << std::endl;

    BasisSet res(name());
    std::cout << "At 1" << std::endl;
    std::cout << "_arg_bases.size(): " << _arg_bases.size() <<  std::endl;
    max_poly_order = valid_index(max_poly_order) ? max_poly_order : _prototype.size();

    std::cout << "at 2, arg_bases dof_IDs:" << std::endl;
    for(Index i = 0; i < _arg_bases.size(); i++) {
      std::cout << "i = " << i << ": " << _arg_bases[i]->dof_IDs() << std::endl;
    }

    Array<BasisSet const *> arg_bases_array(_arg_bases.begin(), _arg_bases.end());

    for(Index i = min_poly_order; i <= max_poly_order; i++) {
      std::cout << "poly_order " << i << std::endl;
      BasisSet tres(name());
      std::cout << "construct poly " << std::endl;
      tres.construct_invariant_polynomials(arg_bases_array, _generating_group, i, 1);
      std::cout << "End construct poly " << std::endl;
      res.append(tres);
    }

    std::cout << "end build_proto" << std::endl;
    return res;
  }


}
