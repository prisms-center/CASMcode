#include "casm/symmetry/SymOp.hh"
#include "casm/symmetry/Orbit_impl.hh"
#include "casm/basis_set/DisplacementDoFTraits.hh"
#include "casm/basis_set/FunctionVisitor.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/crystallography/SimpleStructure.hh"
#include "casm/clusterography/IntegralCluster.hh"
#include "casm/clusterography/ClusterSymCompare_impl.hh"
#include "casm/clusterography/ClusterOrbits_impl.hh"
#include "casm/clex/ClexBasis.hh"
#include "casm/clex/Configuration.hh"
#include "casm/clex/NeighborList.hh"

namespace CASM {
  namespace DoF_impl {
    DoFType::BasicTraits *DisplacementDoFTraits::_clone() const {
      return new DisplacementDoFTraits(*this);
    }

    /// \brief Generate a symmetry representation for the supporting vector space
    Eigen::MatrixXd DisplacementDoFTraits::symop_to_matrix(SymOp const &op) const {
      return op.matrix();
    }

    /// \brief Construct the site basis (if DOF_MODE is LOCAL) for a DoF, given its site
    std::vector<BasisSet> DisplacementDoFTraits::construct_site_bases(Structure const &_prim,
                                                                      std::vector<Orbit<PrimPeriodicSymCompare<IntegralCluster> > > &_asym_unit,
                                                                      jsonParser const &_bspecs) const {
      std::vector<BasisSet> result(_prim.basis().size());

      //std::cout << "Using " << func_type << " site basis functions." << std::endl << std::endl;

      for(Index b = 0; b < _prim.basis().size(); b++) {

        if(!_prim.basis()[b].has_dof(type_name()))
          continue;

        result[b].set_variable_basis(_prim.basis()[b].dof(type_name()));
        //std::cout << "+:+:+:+Created variable set for site " << b << ", size " << result[b].size() << "\n";
      }
      return result;
    }

    //  Apply DoF values for this DoF to _struc
    void DisplacementDoFTraits::apply_dof(ConfigDoF const &_dof, BasicStructure<Site> const &_reference, SimpleStructure &_struc) const {
      _struc.mol_info.coords += _dof.local_dof(type_name()).standard_values();
      _struc.mol_info.dofs[type_name()]["value"] = _dof.local_dof(type_name()).standard_values().transpose();
      // Any considerations here for Selective Dynamics?
    }

  }

  namespace DoFType {
    DoF_impl::DisplacementDoFTraits displacement() {
      return DoF_impl::DisplacementDoFTraits();
    }
  }

}
