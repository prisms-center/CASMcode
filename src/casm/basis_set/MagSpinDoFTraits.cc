#include "casm/basis_set/MagSpinDoFTraits.hh"
#include "casm/basis_set/FunctionVisitor.hh"
#include "casm/crystallography/SymType.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/clusterography/IntegralCluster.hh"
#include "casm/clusterography/ClusterSymCompare_impl.hh"
#include "casm/clusterography/ClusterOrbits_impl.hh"
#include "casm/clex/ClexBasis.hh"
#include "casm/clex/NeighborList.hh"

namespace CASM {
  namespace DoF_impl {
    /// \brief Construct the site basis (if DOF_MODE is LOCAL) for a DoF, given its site
    std::vector<BasisSet> MagSpinDoFTraits::construct_site_bases(Structure const &_prim,
                                                                 std::vector<Orbit<PrimPeriodicSymCompare<IntegralCluster> > > &_asym_unit,
                                                                 jsonParser const &_bspecs) const {
      std::vector<BasisSet> result(_prim.basis().size());

      //std::cout << "Using " << func_type << " site basis functions." << std::endl << std::endl;
      Index order = 0;
      if(_bspecs.contains("max_poly_order")) {
        order = _bspecs["max_poly_order"].get<Index>();
      }
      for(Index b = 0; b < _prim.basis().size(); b++) {

        if(!_prim.basis()[b].has_dof(name()))
          continue;
        BasisSet tresult;
        tresult.set_variable_basis(_prim.basis()[b].dof(name()));
        Array<BasisSet const *> tsubs(1, &tresult);
        result[b].construct_harmonic_polynomials(tsubs, 2, 1, false);
        result[b].get_symmetry_representation(_prim.factor_group());

        result[b].set_name(name() + "_site_func");
        //std::cout << "+:+:+:+Created variable set for site " << b << ", size " << result[b].size() << "\n";
      }
      return result;
    }

  }

  namespace DoFType {

    DoF_impl::MagSpinDoFTraits magspin() {
      return DoF_impl::MagSpinDoFTraits();
    }
  }

}
