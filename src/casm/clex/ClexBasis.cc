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

  ClexBasis::ClexBasis(Structure const &_prim, jsonParser const &_bspecs) :
    m_prim_ptr(&_prim),
    m_basis_builder(std::unique_ptr<ClexBasisBuilder>(new InvariantPolyBasisBuilder())),
    m_bspecs(_bspecs) {

    _populate_site_bases();

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
    return prim().basis.size();
  }

  /// \brief Total number of cluster orbits
  Index ClexBasis::n_orbits() const {
    return m_bset_tree.size();
  }

  /// \brief Total number of basis functions
  Index ClexBasis::n_functions() const {
    Index nf = 0;
    for(auto const &bo : m_bset_tree) {
      nf += bo.size();
    }
    return nf;
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

  void ClexBasis::_populate_site_bases() {
    std::vector<Orbit<IntegralCluster, PrimPeriodicSymCompare<IntegralCluster> > > asym_unit;
    std::ostream nullstream(0);
    make_prim_periodic_asymmetric_unit(prim(),
                                       CASM_TMP::ConstantFunctor<bool>(true),
                                       TOL,
                                       std::back_inserter(asym_unit),
                                       nullstream);

    for(DoFKey const &dof_type : prim().local_dof_types())
      m_site_bases[dof_type] = DoFType::traits(dof_type).construct_site_bases(prim(), asym_unit, bspecs());
  }

  //********************************************************************


  namespace ClexBasis_impl {
  }

}
