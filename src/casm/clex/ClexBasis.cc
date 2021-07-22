#include "casm/basis_set/DoFTraits.hh"
#include "casm/basis_set/FunctionVisitor.hh"
#include "casm/clex/ClexBasis_impl.hh"
#include "casm/clusterography/ClusterDecl.hh"
#include "casm/clusterography/ClusterOrbits_impl.hh"
#include "casm/clusterography/ClusterSymCompare_impl.hh"
#include "casm/container/algorithm.hh"
#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/crystallography/io/BasicStructureIO.hh"
#include "casm/global/eigen.hh"

namespace CASM {

template void
    ClexBasis::generate<std::vector<AperiodicIntegralClusterOrbit>::iterator>(
        std::vector<AperiodicIntegralClusterOrbit>::iterator,
        std::vector<AperiodicIntegralClusterOrbit>::iterator);

template void ClexBasis::generate<
    std::vector<PrimPeriodicIntegralClusterOrbit>::iterator>(
    std::vector<PrimPeriodicIntegralClusterOrbit>::iterator,
    std::vector<PrimPeriodicIntegralClusterOrbit>::iterator);

template BasisSet
ClexBasis::_construct_prototype_basis<AperiodicIntegralClusterOrbit>(
    AperiodicIntegralClusterOrbit const &_orbit,
    std::vector<DoFKey> const &local_keys,
    std::vector<DoFKey> const &global_keys, Index max_poly_order) const;

ClexBasis::ClexBasis(PrimType_ptr _prim_ptr,
                     ClexBasisSpecs const &_basis_set_specs,
                     ParsingDictionary<DoFType::Traits> const *_dof_dict)
    : m_prim_ptr(_prim_ptr),
      m_basis_set_specs(_basis_set_specs),
      m_dof_dict(_dof_dict),
      m_basis_builder(std::unique_ptr<ClexBasisBuilder>(
          new InvariantPolyBasisBuilder("invariant_poly"))) {
  _populate_site_bases();
}

//********************************************************************

Index print_clust_basis(std::ostream &stream, BasisSet _clust_basis,
                        IntegralCluster const &_prototype, Index func_ind,
                        int space, char delim) {
  for (Index np = 0; np < _prototype.size(); np++) {
    stream << std::string(space, ' ');

    stream.setf(std::ios::showpoint, std::ios_base::fixed);
    stream.precision(5);
    stream.width(9);
    _prototype.coordinate(np).print(stream);
    if (delim) stream << delim;
  }
  stream << "\n"
         << "            Basis Functions:\n";

  _clust_basis.set_dof_IDs(sequence<Index>(0, _prototype.size() - 1));
  _clust_basis.accept(OccFuncLabeler("\\phi_{%b,%f}(s_{%n})"));
  Index i;
  for (i = 0; i < _clust_basis.size(); i++) {
    stream << "              \\Phi_{" << func_ind + i
           << "} = " << _clust_basis[i]->tex_formula() << std::endl;
  }
  return _clust_basis.size();
}

//************************************************************

const ClexBasis::PrimType &ClexBasis::prim() const { return *m_prim_ptr; }

ClexBasis::PrimType_ptr ClexBasis::shared_prim() const { return m_prim_ptr; }

ClexBasisSpecs const &ClexBasis::basis_set_specs() const {
  return m_basis_set_specs;
}

ParsingDictionary<DoFType::Traits> const *ClexBasis::dof_dict() const {
  return m_dof_dict;
}

DoFType::Traits const &ClexBasis::lookup_dof_type_traits(
    std::string const &key) const {
  // return DoFType::traits(key);
  return this->dof_dict()->lookup(key);
}

/// \brief Total number of basis sites in primitive cell
Index ClexBasis::n_sublat() const { return prim().basis().size(); }

/// \brief Total number of cluster orbits
Index ClexBasis::n_orbits() const { return m_bset_tree.size(); }

/// \brief Total number of basis functions
Index ClexBasis::n_functions() const {
  Index nf = 0;

  for (auto const &bo : m_bset_tree) {
    if (bo.size()) nf += bo[0].size();
  }
  return nf;
}

//********************************************************************

void ClexBasis::_populate_site_bases() {
  std::vector<Orbit<PrimPeriodicSymCompare<IntegralCluster> > > asym_unit;
  std::ostream nullstream(0);
  make_prim_periodic_asymmetric_unit(shared_prim(),
                                     CASM_TMP::ConstantFunctor<bool>(true), TOL,
                                     std::back_inserter(asym_unit), nullstream);

  for (DoFKey const &dof_type : xtal::all_local_dof_types(prim())) {
    m_site_bases[dof_type] =
        lookup_dof_type_traits(dof_type).construct_site_bases(
            prim(), asym_unit, basis_set_specs().basis_function_specs);
  }
  for (DoFKey const &dof_type : xtal::global_dof_types(prim())) {
    std::vector<BasisSet> tbasis =
        lookup_dof_type_traits(dof_type).construct_site_bases(
            prim(), asym_unit, basis_set_specs().basis_function_specs);
    if (tbasis.empty()) {
      throw std::runtime_error(
          "In ClexBasis::_populate_site_bases(), unable to lookup global DoF "
          "type " +
          dof_type);
    }
    m_global_bases[dof_type] = tbasis;
  }
}

//********************************************************************

namespace ClexBasis_impl {}

namespace sym {
ClexBasis copy_apply(const SymOp &op, ClexBasis clex_basis) {
  ClexBasis new_clex_basis = clex_basis;
  return new_clex_basis;
}
} // namespace sym

}  // namespace CASM
