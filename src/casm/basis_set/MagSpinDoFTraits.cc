#include "casm/basis_set/MagSpinDoFTraits.hh"
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
    Eigen::MatrixXd MagSpinDoFTraits::symop_to_matrix(SymOp const &op) const {
      return op.matrix().determinant() * op.matrix();
    }


  }
}
