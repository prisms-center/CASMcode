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
             S(0, 0)*S(0, 0), S(0, 1)*S(0, 1), S(0, 2)*S(0, 2), sqrt(2)*S(0, 1)*S(0, 2), sqrt(2)*S(0, 2)*S(0, 0), sqrt(2)*S(0, 0)*S(0, 1),
             S(1, 0)*S(1, 0), S(1, 1)*S(1, 1), S(1, 2)*S(1, 2), sqrt(2)*S(1, 1)*S(1, 2), sqrt(2)*S(1, 2)*S(1, 0), sqrt(2)*S(1, 0)*S(1, 1),
             S(2, 0)*S(2, 0), S(2, 1)*S(2, 1), S(2, 2)*S(2, 2), sqrt(2)*S(2, 1)*S(2, 2), sqrt(2)*S(2, 2)*S(2, 0), sqrt(2)*S(2, 0)*S(2, 1),
             sqrt(2)*S(1, 0)*S(2, 0), sqrt(2)*S(1, 1)*S(2, 1), sqrt(2)*S(1, 2)*S(2, 2), S(1, 1)*S(2, 2) + S(1, 2)*S(2, 1), S(1, 0)*S(2, 2) + S(1, 2)*S(2, 0), S(1, 1)*S(2, 0) + S(1, 0)*S(2, 1),
             sqrt(2)*S(2, 0)*S(0, 0), sqrt(2)*S(2, 1)*S(0, 1), sqrt(2)*S(2, 2)*S(0, 2), S(2, 1)*S(0, 2) + S(2, 2)*S(0, 1), S(2, 0)*S(0, 2) + S(2, 2)*S(0, 0), S(2, 1)*S(0, 0) + S(2, 0)*S(0, 1),
             sqrt(2)*S(0, 0)*S(1, 0), sqrt(2)*S(0, 1)*S(1, 1), sqrt(2)*S(0, 2)*S(1, 2), S(0, 1)*S(1, 2) + S(0, 2)*S(1, 1), S(0, 0)*S(1, 2) + S(0, 2)*S(1, 0), S(0, 1)*S(1, 0) + S(0, 0)*S(1, 1);
      return result;
    }



    /// \brief Construct the site basis (if DOF_MODE is LOCAL) for a DoF, given its site
    std::vector<BasisSet> StrainDoFTraits::construct_site_bases(Structure const &_prim,
                                                                std::vector<Orbit<IntegralCluster, PrimPeriodicSymCompare<IntegralCluster> > > &_asym_unit,
                                                                jsonParser const &_bspecs) const {


      //std::cout << "Using " << func_type << " site basis functions." << std::endl << std::endl;
      if(_prim.global_dofs().find(type_name()) == _prim.global_dofs().end())
        return std::vector<BasisSet>();

      std::vector<BasisSet> result(1);
      result[0].set_variable_basis(_prim.global_dof(type_name()));

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
