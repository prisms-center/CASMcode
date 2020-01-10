#include "casm/basis_set/StrainDoFTraits.hh"
#include "casm/basis_set/FunctionVisitor.hh"
#include "casm/basis_set/BasisSet.hh"
#include "casm/crystallography/SymType.hh"
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

    /// \brief Construct the site basis (if DOF_MODE is LOCAL) for a DoF, given its site
    std::vector<BasisSet> StrainDoFTraits::construct_site_bases(Structure const &_prim,
                                                                std::vector<Orbit<PrimPeriodicSymCompare<IntegralCluster> > > &_asym_unit,
                                                                jsonParser const &_bspecs) const {


      //std::cout << "Using " << func_type << " site basis functions." << std::endl << std::endl;
      if(_prim.global_dofs().find(name()) == _prim.global_dofs().end())
        return std::vector<BasisSet>();

      std::vector<BasisSet> result(1);
      result[0].set_variable_basis(_prim.global_dof(name()));

      return result;
    }




    /// \brief Apply DoF values for this DoF to _struc
    void StrainDoFTraits::apply_dof(ConfigDoF const &_dof, BasicStructure const &_reference, SimpleStructure &_struc) const {
      Eigen::VectorXd unrolled_metric = _dof.global_dof(name()).standard_values();
      StrainConverter c(m_metric);
      Eigen::Matrix3d F = c.unrolled_strain_metric_to_F(unrolled_metric);

      _struc.lat_column_mat = F * _struc.lat_column_mat;
      _struc.mol_info.coords = F * _struc.mol_info.coords;
      _struc.atom_info.coords = F * _struc.atom_info.coords;
    }

    jsonParser StrainDoFTraits::dof_to_json(ConfigDoF const &_dof, BasicStructure const &_reference) const {
      Eigen::VectorXd unrolled_metric = _dof.global_dof(name()).standard_values();
      StrainConverter c(m_metric);
      Eigen::Matrix3d F = c.unrolled_strain_metric_to_F(unrolled_metric);

      jsonParser json;
      json["deformation"] = F;
      to_json_array(unrolled_metric, json["value"]);
      return json;
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
