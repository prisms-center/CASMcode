#include "casm/basis_set/StrainDoFTraits.hh"
#include "casm/basis_set/FunctionVisitor.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/symmetry/SymOp.hh"
#include "casm/clusterography/IntegralCluster.hh"
#include "casm/clusterography/ClusterSymCompare_impl.hh"
#include "casm/clusterography/ClusterOrbits_impl.hh"
#include "casm/clex/ClexBasis.hh"
#include "casm/clex/NeighborList.hh"

namespace CASM {
  namespace DoF_impl {

    /// \brief Generate a symmetry representation for the supporting vector space
    Eigen::MatrixXd StrainDoFTraits::symop_to_matrix(SymOp const &op) const {
      Eigen::MatrixXd result(6, 6);
      auto const &S = op.matrix();
      result <<
             S(1, 1)*S(1, 1), S(1, 2)*S(1, 2), S(1, 3)*S(1, 3), sqrt(2)*S(1, 2)*S(1, 3), sqrt(2)*S(1, 3)*S(1, 1), sqrt(2)*S(1, 1)*S(1, 2),
             S(2, 1)*S(2, 1), S(2, 2)*S(2, 2), S(2, 3)*S(2, 3), sqrt(2)*S(2, 2)*S(2, 3), sqrt(2)*S(2, 3)*S(2, 1), sqrt(2)*S(2, 1)*S(2, 2),
             S(3, 1)*S(3, 1), S(3, 2)*S(3, 2), S(3, 3)*S(3, 3), sqrt(2)*S(3, 2)*S(3, 3), sqrt(2)*S(3, 3)*S(3, 1), sqrt(2)*S(3, 1)*S(3, 2),
             sqrt(2)*S(2, 1)*S(3, 1), sqrt(2)*S(2, 2)*S(3, 2), sqrt(2)*S(2, 3)*S(3, 3), S(2, 2)*S(3, 3) + S(2, 3)*S(3, 2), S(2, 1)*S(3, 3) + S(2, 3)*S(3, 1), S(2, 2)*S(3, 1) + S(2, 1)*S(3, 2),
             sqrt(2)*S(3, 1)*S(1, 1), sqrt(2)*S(3, 2)*S(1, 2), sqrt(2)*S(3, 3)*S(1, 3), S(3, 2)*S(1, 3) + S(3, 3)*S(1, 2), S(3, 1)*S(1, 3) + S(3, 3)*S(1, 1), S(3, 2)*S(1, 1) + S(3, 1)*S(1, 2),
             sqrt(2)*S(1, 1)*S(2, 1), sqrt(2)*S(1, 2)*S(2, 2), sqrt(2)*S(1, 3)*S(2, 3), S(1, 2)*S(2, 3) + S(1, 3)*S(2, 2), S(1, 1)*S(2, 3) + S(1, 3)*S(2, 1), S(1, 2)*S(2, 1) + S(1, 1)*S(2, 2);
      return result;
    }

  }

  namespace DoFType {
    DoF_impl::StrainDoFTraits GLstrain() {
      return DoF_impl::StrainDoFTraits("GL");
    }

    DoF_impl::StrainDoFTraits EAstrain() {
      return DoF_impl::StrainDoFTraits("EA");
    }

    DoF_impl::StrainDoFTraits Hstrain() {
      return DoF_impl::StrainDoFTraits("H");
    }

  }
}
