#include "casm/symmetry/SymOp.hh"
#include "casm/basis_set/StrainDoFTraits.hh"
#include "casm/basis_set/FunctionVisitor.hh"
#include "casm/basis_set/BasisSet.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/crystallography/SimpleStructure.hh"
#include "casm/crystallography/Site.hh"
#include "casm/clusterography/IntegralCluster.hh"
#include "casm/clusterography/ClusterSymCompare_impl.hh"
#include "casm/clusterography/ClusterOrbits_impl.hh"
#include "casm/clex/ClexBasis.hh"
#include "casm/clex/ConfigDoF.hh"
#include "casm/clex/NeighborList.hh"
#include "casm/strain/StrainConverter.hh"

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




    /// \brief Apply DoF values for this DoF to _struc
    void StrainDoFTraits::apply_dof(ConfigDoF const &_dof, BasicStructure<Site> const &_reference, SimpleStructure &_struc) const {
      Eigen::VectorXd unrolled_metric = _dof.global_dof(type_name()).standard_values();
      StrainConverter c(m_metric);
      Eigen::Matrix3d F = c.unrolled_strain_metric_to_F(unrolled_metric);

      _struc.lat_column_mat = F * _struc.lat_column_mat;
      _struc.mol_info.coords = F * _struc.mol_info.coords;
      _struc.atom_info.coords = F * _struc.atom_info.coords;
      _struc.global_dofs[type_name()]["deformation"] = F;
      to_json_array(unrolled_metric, _struc.global_dofs[type_name()]["value"]);
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
