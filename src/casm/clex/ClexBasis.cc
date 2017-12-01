#include "casm/clex/ClexBasis_impl.hh"
#include "casm/clusterography/ClusterDecl.hh"
#include "casm/clusterography/ClusterOrbits_impl.hh"
#include "casm/clusterography/ClusterSymCompare_impl.hh"
#include "casm/container/algorithm.hh"
#include "casm/basis_set/DoFTraits.hh"
#include "casm/basis_set/FunctionVisitor.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/app/AppIO.hh" // necessary for write_prim() function

#include "casm/CASM_global_Eigen.hh"

namespace CASM {

  template void ClexBasis::generate<std::vector<AperiodicIntegralClusterOrbit>::iterator>(
    std::vector<AperiodicIntegralClusterOrbit>::iterator,
    std::vector<AperiodicIntegralClusterOrbit>::iterator,
    jsonParser const &,
    Index);

  template void ClexBasis::generate<std::vector<PrimPeriodicIntegralClusterOrbit>::iterator>(
    std::vector<PrimPeriodicIntegralClusterOrbit>::iterator,
    std::vector<PrimPeriodicIntegralClusterOrbit>::iterator,
    jsonParser const &,
    Index);

  template
  BasisSet ClexBasis::_construct_prototype_basis<AperiodicIntegralClusterOrbit>(
    AperiodicIntegralClusterOrbit const &_orbit,
    std::vector<DoFKey> const &local_keys,
    std::vector<DoFKey> const &global_keys,
    Index max_poly_order) const;

  ClexBasis::ClexBasis(Structure const &_prim) :
    m_prim_ptr(&_prim) {

  }

  //********************************************************************

  Index print_clust_basis(std::ostream &stream,
                          BasisSet _clust_basis,
                          IntegralCluster const &_prototype,
                          Index func_ind,
                          int space,
                          char delim) {
    for(Index np = 0; np < _prototype.size(); np++) {

      stream << std::string(space, ' ');

      stream.setf(std::ios::showpoint, std::ios_base::fixed);
      stream.precision(5);
      stream.width(9);
      _prototype.coordinate(np).print(stream);
      stream << "  basis_index: " << _prototype[np].sublat() << "  clust_index: " << np << " ";
      if(delim)
        stream << delim;
    }
    stream << "\n"
           << "            Basis Functions:\n";

    _clust_basis.set_dof_IDs(sequence<Index>(0, _prototype.size() - 1));
    _clust_basis.accept(OccFuncLabeler("\\phi_%b_%f(s_%n)"));
    Index i;
    for(i = 0; i < _clust_basis.size(); i++) {
      stream << "              \\Phi_" << func_ind + i << " = " << _clust_basis[i]->tex_formula() << std::endl;
    }
    return _clust_basis.size();
  }


  //************************************************************

  Structure const &ClexBasis::prim() const {
    return *m_prim_ptr;
  }

  /// \brief Total number of basis sites in primitive cell
  Index ClexBasis::n_sublat() const {
    throw std::runtime_error("ClexBasis::n_sublat is not implemented");
  }

  /// \brief Total number of basis functions
  Index ClexBasis::n_functions() const {
    throw std::runtime_error("ClexBasis::n_functions is not implemented");
  }


  /*
  void print_eci_in(ClexBasis const &_clex_basis,
                    std::ostream &out,
                    std::vector<Orbit<IntegralCluster> > const &_tree) const {
    if(_tree.index.size() != _tree.size())
      get_index();
    if(_tree.subcluster.size() != _tree.size())
      _tree.get_hierarchy();

    out << std::left
        << std::setw(8) << "label"
        << std::setw(8) << "weight"
        << std::setw(8) << "mult"
        << std::setw(8) << "size"
        << std::setw(12) << "length"
        << std::setw(8) << "hierarchy" << std::endl;


    int functioncount = 0;
    for(Index i = 0; i < _tree.size(); i++) {
      for(Index j = 0; j < _tree.size(i); j++) {

        for(Index k = 0; k < m_bset_tree[_tree.orbit(i, j).index()].size(); k++, functioncount++) {

          out << std::left
              << std::setw(8) << functioncount
              << std::setw(8) << 0
              << std::setw(8) << _tree.orbit(i, j).size()
              << std::setw(8) << _tree.prototype(i, j).size()
              << std::setw(12) << _tree.prototype(i, j).max_length();

          // print hierarchy
          out << std::left << std::setw(8) << 0;
          for(Index l = 0; l < _tree.subcluster[ _tree.index[i][j]].size(); l++) {
            out << std::left
                << std::setw(8) << _tree.subcluster[ _tree.index[i][j] ][l];
          }
          out << '\n' << std::flush;
        }

      }
    }

    //std::cout << "finish print_eci_in" << std::endl;

  }
  */
  //********************************************************************

  void ClexBasis::_populate_site_bases(Structure const &_prim) {
    std::vector<Orbit<IntegralCluster, PrimPeriodicSymCompare<IntegralCluster> > > asym_unit;
    std::ostream nullstream(0);
    make_prim_periodic_asymmetric_unit(_prim,
                                       CASM_TMP::ConstantFunctor<bool>(true),
                                       TOL,
                                       std::back_inserter(asym_unit),
                                       nullstream);
    for(DoFKey const &dof_type : ClexBasis_impl::extract_dof_types(_prim))
      m_site_bases[dof_type] = DoFType::traits(dof_type).construct_site_bases(_prim, asym_unit, bspecs());
  }

  //********************************************************************


  namespace ClexBasis_impl {
    std::vector<ClexBasis::DoFKey> extract_dof_types(Structure const &_prim) {
      throw std::runtime_error("ClexBasis_impl::extract_dof_types() is not yet implemented!!\n");
      return std::vector<ClexBasis::DoFKey>();
    }

    BasisSet construct_clust_dof_basis(IntegralCluster const &_clust, std::vector<BasisSet const *> const &site_dof_sets) {
      throw std::runtime_error("ClexBasis_impl::construct_clust_dof_basis() needs to be re-implemented!\n");
      BasisSet result;
      // UPDATE to replace _clust.clust_group(), _clust.nlist_inds(), etc
      /*
      result.set_dof_IDs(_clust.nlist_inds());
      std::vector<SymGroupRep const *> subspace_reps;
      for(BasisSet const *site_bset_ptr : site_dof_sets) {
        if(site_bset_ptr) {
          result.append(*site_bset_ptr);
          subspace_reps.push_back(SymGroupRep::RemoteHandle(_clust.clust_group(),
                                                            site_bset_ptr->basis_symrep_ID()).rep_ptr());
        }
        else {
          subspace_reps.push_back(SymGroupRep::RemoteHandle(_clust.clust_group(),
                                                            SymGroupRepID::identity(0)).rep_ptr());
        }
      }
      result.set_basis_symrep_ID(permuted_direct_sum_rep(*(_clust.permute_rep().rep_ptr()),
                                                               subspace_reps).add_copy_to_master());
      */
      return result;


    }

  }

}
