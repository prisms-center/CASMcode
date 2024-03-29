#ifndef CASM_ClexBasis_impl
#define CASM_ClexBasis_impl

#include "casm/basis_set/DoFTraits.hh"
#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/clex/ClexBasis.hh"
#include "casm/clex/NeighborList.hh"
#include "casm/clex/OrbitFunctionTraits.hh"
#include "casm/clusterography/ClusterSpecs_impl.hh"
#include "casm/clusterography/IntegralCluster.hh"
#include "casm/container/algorithm.hh"
#include "casm/crystallography/UnitCellCoord.hh"
#include "casm/symmetry/SymGroup.hh"
#include "casm/symmetry/SymGroupRep.hh"
#include "casm/symmetry/SymRepTools.hh"

namespace CASM {

/// \brief The extended equivalence map maps a prototype to equivalent
///     clusters in the current and equivalent orbitrees
///
/// Usage:
/// \code
/// IntegralCluster equiv_cluster = sym_compare.copy_apply(
///     extended_equivalence_map[phenom_index][cluster_index][symop_index],
///     prototype)
/// \endcode
/// Where `equiv_cluster` is the cluster_index-th cluster around the
/// phenom_index-th equivalent phenomenal cluster, for all symop_index.
///
/// For periodic ClexBasis, phenom_index = 0 only, and
/// extended_equivalence_map[0] == orbit.equivalence_map()
template <typename OrbitType>
SymGroupRepID make_canonization_rep(
    OrbitType const &orbit,
    multivector<SymOp>::X<3> const &extended_equivalence_map) {
  auto const &master_group = orbit.generating_group()[0].master_group();
  Index total_size = 0;
  for (auto const &table : extended_equivalence_map) {
    for (auto const &row : table) {
      total_size += row.size();
    }
  }

  if (master_group.size() != total_size) {
    std::stringstream msg;
    msg << "Error in make_canonization_rep: "
        << "extended_equivalence_map total size  (" << total_size
        << ") != master group size (" << master_group.size() << ")";
    throw libcasm_runtime_error(msg.str());
  }

  if (orbit.size() == 0) {
    return SymGroupRepID::identity(0);
  }

  SymGroupRepID canonization_rep_ID = master_group.allocate_representation();

  for (Index k = 0; k < extended_equivalence_map[0][0].size(); k++) {
    IntegralCluster equiv_cluster = orbit.sym_compare().copy_apply(
        extended_equivalence_map[0][0][k], orbit.prototype());

    std::unique_ptr<SymOpRepresentation> new_rep =
        orbit.sym_compare().canonical_transform(equiv_cluster)->inverse();

    for (Index i = 0; i < extended_equivalence_map.size(); i++) {
      for (Index j = 0; j < extended_equivalence_map[i].size(); j++) {
        extended_equivalence_map[i][j][k].set_rep(canonization_rep_ID,
                                                  *new_rep);
      }
    }
  }
  return canonization_rep_ID;
}

/// allow basis_function_specs to specify max_poly_order by branch, orbit, etc.
template <typename OrbitType>
Index _orbit_max_poly_order(OrbitType const &orbit,
                            BasisFunctionSpecs const &basis_function_specs) {
  auto it = basis_function_specs.orbit_branch_max_poly_order.find(
      orbit.prototype().size());
  if (it != basis_function_specs.orbit_branch_max_poly_order.end()) {
    return it->second;
  }
  return basis_function_specs.global_max_poly_order;
}

template <typename OrbitIteratorType>
void ClexBasis::generate(OrbitIteratorType _orbit_begin,
                         OrbitIteratorType _orbit_end) {
  std::vector<DoFKey> dof_keys =
      basis_set_specs().basis_function_specs.dof_keys;
  std::vector<DoFKey> global_keys;
  std::vector<DoFKey> local_keys;

  if (dof_keys.empty()) {
    for (auto const &dof : m_global_bases) {
      global_keys.push_back(dof.first);
    }
    for (auto const &dof : m_site_bases) {
      local_keys.push_back(dof.first);
    }
  } else if (dof_keys.size() > 1 || dof_keys[0] != "none") {
    // separate local_args from global_args
    for (DoFKey const &key : dof_keys) {
      if (m_global_bases.find(key) != m_global_bases.end()) {
        global_keys.push_back(key);
      } else if (m_site_bases.find(key) != m_site_bases.end()) {
        local_keys.push_back(key);
      } else {
        assert(0);
        std::stringstream ss;
        ss << "Attempting to build Clex basis set, but missing degree of "
              "freedom '"
           << key << "'\n";
        throw std::runtime_error(ss.str());
      }
    }
  }

  if (basis_set_specs().cluster_specs->periodicity_type() ==
      CLUSTER_PERIODICITY_TYPE::LOCAL) {
    m_bset_tree_equivalence_map = make_equivalents_generating_ops(
        shared_prim(),
        basis_set_specs().cluster_specs->get_phenomenal_cluster(),
        basis_set_specs().cluster_specs->get_generating_group());
  } else {
    m_bset_tree_equivalence_map.push_back(SymOp());
  }
  Index n_equiv = m_bset_tree_equivalence_map.size();

  m_bset_tree.resize(n_equiv);
  for (Index i = 0; i < n_equiv; ++i) {
    m_bset_tree[i].resize(std::distance(_orbit_begin, _orbit_end));
  }

  Index orbit_index = 0;
  for (; _orbit_begin != _orbit_end; ++_orbit_begin, ++orbit_index) {
    auto const &orbit = *_orbit_begin;

    // extended equivalence map
    multivector<SymOp>::X<3> extended_equivalence_map =
        make_extended_equivalence_map(orbit.equivalence_map(),
                                      m_bset_tree_equivalence_map);
    // make canonization rep
    SymGroupRepID canonization_rep_ID =
        make_canonization_rep(orbit, extended_equivalence_map);

    // construct prototype cluster basis functions
    Index max_poly_order =
        _orbit_max_poly_order(orbit, basis_set_specs().basis_function_specs);
    BasisSet prototype_basis_set = _construct_prototype_basis(
        orbit, local_keys, global_keys, max_poly_order, canonization_rep_ID);

    // construct equivalent cluster basis functions
    for (Index i = 0; i < n_equiv; ++i) {
      m_bset_tree[i][orbit_index].reserve(orbit.size());
    }
    for (Index i = 0; i < n_equiv; ++i) {
      for (Index j = 0; j < orbit.size(); j++) {
        m_bset_tree[i][orbit_index].push_back(
            extended_equivalence_map[i][j][0] * prototype_basis_set);
      }
    }
  }
}

//*******************************************************************************************
// @param local_args[i][j] is BasisSet for i'th DoFspace at j'th site of cluster
//                                site  0            site 1            site 2
// DoFKey 0: displacement       {x0,y0,z0}       {x1, y1, z1}       {x2,y2,z2}
// DoFKey 1: configuration       {pA0,pB0}            {}            {pA2, pB2}
//
// Step 1:  Get the kroenecker product of cluster permutation with DoF symrep
//
//  permutation  |    kronecker prod  |    DoF Symrep (e.g., x--y displacement)
//   [ 0  1 ]              v                     [cos -sin]
//   [ 1  0 ]             XkX                    [sin  cos]
//
//       [x0]      [  0    0   cos -sin ]   [x0]
//       [y0]      [  0    0   sin  cos ]   [y0]
//   S * [x1]  =   [ cos -sin   0    0  ]   [x1]
//       [x2]      [ sin  cos   0    0  ]   [x2]
//
// ----------------------------------------------------------------------
//
// Step 2: mix-in @param global_args to get all_argsets
//
// GLOBAL ARGS
//
// strain             { e1, e2, e3, e4, e5, e6}
// composition        { comp_a, comp_b }
//
// arg_subsets =   [ {x0,y0,z0,x1,y1,z1,x2,y2,z2},
//                   {pA0,pB0,pA2,pB2},
//                   {e1,e2,e3,e4,e5,e6},
//                   {comp_a,comp_b}]
//
template <typename OrbitType>
BasisSet ClexBasis::_construct_prototype_basis(
    OrbitType const &_orbit, std::vector<DoFKey> const &local_keys,
    std::vector<DoFKey> const &global_keys, Index max_poly_order,
    SymGroupRepID canonization_rep_ID) const {
  // std::cout<<"In IntegralCluster::generate_clust_basis, the size of this
  // cluster is:"<<size()<<std::endl; std::cout<<"valid_index evaluates
  // to:"<<valid_index(max_poly_order)<<std::endl;

  // Default polynomial order is cluster size
  // if(!valid_index(max_poly_order))
  max_poly_order = max(max_poly_order, Index(_orbit.prototype().size()));

  // std::cout<<"Max_poly_order "<<max_poly_order<<std::endl;

  // record pointers to global dof arguments
  std::vector<BasisSet const *> arg_subsets;
  for (DoFKey const &key : global_keys) {
    auto find_it = m_global_bases.find(key);
    if (find_it != m_global_bases.end())
      arg_subsets.push_back(&(find_it->second[0]));
    else
      throw std::runtime_error(
          "Unable to construct basis sets. No known global DoF: " + key + "\n");
  }

  // copy local site bases to a temporary location where we can alter their DoF
  // IDs
  std::vector<BasisSet> all_local;
  all_local.reserve(local_keys.size());

  // Loop over dof's
  for (DoFKey const &key : local_keys) {
    // Make copies of local arguments to ensure that they are distinguishable by
    // their DoF_IDs i.e., make copies in 'tlocal' and reset the DoF_IDs to
    // {0,1,2,etc...}
    auto find_it = m_site_bases.find(key);
    if (find_it == m_site_bases.end())
      throw std::runtime_error(
          "Unable to construct basis sets. No known local DoF: " + key + "\n");

    std::vector<BasisSet> const &arg_vec(find_it->second);
    // std::cout << "skdsa Working on DoF " << key << " size() : " <<
    // find_it->second.size() << "\n";
    std::vector<BasisSet> tlocal;
    tlocal.reserve(_orbit.prototype().size());
    // Loop over sites
    for (Index i = 0; i < _orbit.prototype().size(); i++) {
      tlocal.push_back(arg_vec[_orbit.prototype()[i].sublattice()]);
      tlocal.back().set_dof_IDs(std::vector<Index>(1, i));
    }

    all_local.push_back(ClexBasis_impl::construct_proto_dof_basis(
        _orbit, BasisSet::ArgList(tlocal), canonization_rep_ID));

    if (all_local.back().size()) arg_subsets.push_back(&(all_local.back()));
  }
  SymGroup clust_group(_orbit.equivalence_map(0).first,
                       _orbit.equivalence_map(0).second);
  return m_basis_builder->build_proto(_orbit.prototype(), clust_group,
                                      arg_subsets, max_poly_order, 1);
}

namespace ClexBasis_impl {

/// Tool for hiding type of orbits
template <typename FunctorOfClexBasisAndOrbits>
struct ForClexBasisAndOrbits {
  ForClexBasisAndOrbits(std::shared_ptr<Structure const> const &_shared_prim,
                        ClexBasisSpecs const &_basis_set_specs,
                        FunctorOfClexBasisAndOrbits _f)
      : shared_prim(_shared_prim), basis_set_specs(_basis_set_specs), f(_f) {}

  std::shared_ptr<Structure const> shared_prim;
  ClexBasisSpecs const &basis_set_specs;
  FunctorOfClexBasisAndOrbits f;

  /// Gets called by `for_all_orbits`
  template <typename OrbitVecType>
  void operator()(OrbitVecType const &orbits) {
    // construct ClexBasis
    ParsingDictionary<DoFType::Traits> const *dof_dict =
        &DoFType::traits_dict();
    ClexBasis clex_basis{shared_prim, basis_set_specs, dof_dict};
    clex_basis.generate(orbits.begin(), orbits.end());
    f(clex_basis, orbits);
  }
};

}  // namespace ClexBasis_impl

/// \brief Act on ClexBasis and orbits generated by ClexBasisSpecs
///
/// \param shared_prim Shared prim structure
/// \param clex_basis_specs Describes how to generate cluster basis functions
/// \param f Functor that implements:
///
///     template <typename OrbitVecType>
///     void operator()(ClexBasis const &clex_basis, OrbitVecType const
///          &orbits) const;
///
/// This method will construct orbits as specified by
/// `clex_basis_specs.cluster_specs`, use those orbits to construct a
/// ClexBasis, and then call `f(clex_basis, orbits)`. This is a way to
/// use templated methods of OrbitVecType and ClexBasis.
///
template <typename FunctorOfClexBasisAndOrbits>
void for_clex_basis_and_orbits(
    std::shared_ptr<Structure const> const &shared_prim,
    ClexBasisSpecs const &clex_basis_specs, std::ostream &status,
    FunctorOfClexBasisAndOrbits f) {
  ClexBasis_impl::ForClexBasisAndOrbits g{shared_prim, clex_basis_specs, f};
  for_all_orbits(*clex_basis_specs.cluster_specs, status, g);
}

namespace ClexBasis_impl {
template <typename OrbitType>
BasisSet construct_proto_dof_basis(OrbitType const &_orbit,
                                   BasisSet::ArgList const &site_dof_sets,
                                   SymGroupRepID canonization_rep_ID) {
  SymGroup clust_group(_orbit.equivalence_map(0).first,
                       _orbit.equivalence_map(0).second);

  std::vector<SymGroupRep const *> subspace_reps;
  for (BasisSet const *site_bset_ptr : site_dof_sets) {
    if (site_bset_ptr) {
      subspace_reps.push_back(SymGroupRep::RemoteHandle(
                                  clust_group, site_bset_ptr->basis_symrep_ID())
                                  .rep_ptr());
    } else {
      subspace_reps.push_back(
          SymGroupRep::RemoteHandle(clust_group, SymGroupRepID::identity(0))
              .rep_ptr());
    }
  }
  SymGroupRep const *permute_rep =
      SymGroupRep::RemoteHandle(clust_group, canonization_rep_ID).rep_ptr();

  BasisSet result = direct_sum(site_dof_sets);
  result.set_basis_symrep_ID(clust_group.master_group().add_representation(
      permuted_direct_sum_rep(*(permute_rep), subspace_reps)));

  return result;
}

}  // namespace ClexBasis_impl
}  // namespace CASM
#endif
