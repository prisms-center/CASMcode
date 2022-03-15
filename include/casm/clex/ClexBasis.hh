#ifndef CLEXBASIS_HH
#define CLEXBASIS_HH

#include <string>

#include "casm/basis_set/BasisSet.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/clex/ClexBasisSpecs.hh"
#include "casm/clex/OrbitFunctionTraits.hh"
#include "casm/clusterography/ClusterDecl.hh"
#include "casm/crystallography/DoFDecl.hh"
#include "casm/global/enum.hh"
#include "casm/misc/ParsingDictionary.hh"

namespace CASM {
namespace xtal {
class Site;
class BasicStructure;
class UnitCell;
}  // namespace xtal
using xtal::BasicStructure;
using xtal::Site;
using xtal::UnitCell;

namespace clexulator {
class PrimNeighborList;
}
using clexulator::PrimNeighborList;

namespace DoFType {
class Traits;
}
class ClexBasisBuilder;

class ClexBasis {
 public:
  typedef std::vector<BasisSet> BSetOrbit;
  typedef std::vector<BSetOrbit>::const_iterator BSetOrbitIterator;

  typedef Structure PrimType;
  typedef std::shared_ptr<const Structure> PrimType_ptr;

  /// Constructor
  ///
  /// \param _prim_ptr Shared pointer to prim Structure
  /// \param _basis_function_specs Parameters specifying how to construct basis
  /// functions. See
  ///        BasisFunctionSpecs for documentation of options.
  ///
  /// The constructor constructs site bases, so that they are ready to be used
  /// when the `generate` method is called with a range of cluster orbits to
  /// generate cluster expasion basis functions.
  ClexBasis(PrimType_ptr _prim_ptr, ClexBasisSpecs const &_basis_set_specs,
            ParsingDictionary<DoFType::Traits> const *_dof_dict);

  PrimType const &prim() const;

  PrimType_ptr shared_prim() const;

  ClexBasisSpecs const &basis_set_specs() const;

  ParsingDictionary<DoFType::Traits> const *dof_dict() const;

  DoFType::Traits const &lookup_dof_type_traits(std::string const &key) const;

  /// \brief Total number of basis sites in primitive cell
  Index n_sublat() const;

  /// \brief Total number of BasisSet orbits
  Index n_orbits() const;

  /// \brief Total number of basis functions
  Index n_functions() const;

  /// \brief The number of equivalent bset tree
  ///
  /// Notes:
  /// - For periodic ClexBasis: n_equivalent_tree == 1
  /// - For local ClexBasis: n_equivalent_tree >= 1
  Index n_equivalent_tree() const { return m_bset_tree.size(); }

  /// \brief An equivalence map for bset tree
  ///
  /// Bset tree are equivalent according to:
  /// - bset_tree[i] == bset_tree_equivalence_map[i] * bset_tree[0]
  std::vector<SymOp> const &bset_tree_equivalence_map() const {
    return m_bset_tree_equivalence_map;
  }

  /// \brief Set the current equivalent bset tree
  ///
  /// Usage:
  /// - This is used for printing equivalent local clexulator
  void set_current_tree(Index tree_index) { m_tree_index = tree_index; }

  /// \brief Const access of clust basis of orbit @param orbit_ind and
  /// equivalent cluster @param equiv_ind
  BasisSet const &clust_basis(Index orbit_ind, Index equiv_ind) const {
    return m_bset_tree[m_tree_index][orbit_ind][equiv_ind];
  }

  /// \brief Const access of BSetOrbit of orbit @param orbit_ind
  BSetOrbit const &bset_orbit(Index orbit_ind) const {
    return m_bset_tree[m_tree_index][orbit_ind];
  }

  /// \brief Const iterator to first BasisSet orbit
  BSetOrbitIterator begin() const { return m_bset_tree[m_tree_index].cbegin(); }

  /// \brief Const iterator to first BasisSet orbit
  BSetOrbitIterator cbegin() const {
    return m_bset_tree[m_tree_index].cbegin();
  }

  /// \brief Const past-the-end iterator for BasisSet orbits
  BSetOrbitIterator end() const { return m_bset_tree[m_tree_index].cbegin(); }

  /// \brief Const past-the-end iterator for BasisSet orbits
  BSetOrbitIterator cend() const { return m_bset_tree[m_tree_index].cbegin(); }

  /// \brief Const access to dictionary of all site BasisSets
  std::map<DoFKey, std::vector<BasisSet>> const &site_bases() const {
    return m_site_bases;
  }

  /// \brief Const access to dictionary of all global BasisSets
  std::map<DoFKey, std::vector<BasisSet>> const &global_bases() const {
    return m_global_bases;
  }

  /// \brief generate clust_basis for all equivalent clusters in @param
  /// _orbitree
  ///
  /// \param _begin, _end A range of cluster orbits for which to generate
  /// cluster basis functions
  template <typename OrbitIterType>
  void generate(OrbitIterType _begin, OrbitIterType _end);

 private:
  template <typename OrbitType>
  BasisSet _construct_prototype_basis(OrbitType const &_orbit,
                                      std::vector<DoFKey> const &local_keys,
                                      std::vector<DoFKey> const &global_keys,
                                      Index max_poly_order,
                                      SymGroupRepID canonization_rep_ID) const;

  /// \brief Performs heavy lifting for populating site bases in m_site_bases
  void _populate_site_bases();

  PrimType_ptr m_prim_ptr;

  ClexBasisSpecs m_basis_set_specs;

  ParsingDictionary<DoFType::Traits> const *m_dof_dict;

  /// \brief pointer to class that constructs cluster functions
  notstd::cloneable_ptr<ClexBasisBuilder> m_basis_builder;

  /// \brief An equivalence map for bset_tree
  ///
  /// Bset tree are equivalent according to:
  /// - m_bset_tree[i] == m_bset_tree_equivalence_map[i] * m_bset_tree[0]
  std::vector<SymOp> m_bset_tree_equivalence_map;

  /// \brief Collections of all cluster BasisSets, one per cluster orbit (may
  /// be > 1 for local ClexBasis)
  std::vector<std::vector<BSetOrbit>> m_bset_tree;

  /// \brief Determines which bset_tree information is returned about
  ///
  /// Note:
  /// - This is used as `m_bset_tree[tree_index]` to give one of the equivalent
  ///   bset_tree (i.e. an Orbitree / cluster expansion)
  /// - Iterate, varying this via `set_tree_index`, to print out equivalent
  ///   local clexulator
  Index m_tree_index;

  /// \brief Dictionary of all site BasisSets, initialized on construction
  /// m_site_basis[DOF][b] gives site basis functions for 'DOF' at site 'b' of
  /// prim
  std::map<DoFKey, std::vector<BasisSet>> m_site_bases;

  /// \brief Dictionary of all global BasisSets, initialized
  std::map<DoFKey, std::vector<BasisSet>> m_global_bases;
};

/// \brief Act on ClexBasis and orbits generated by ClexBasisSpecs
template <typename FunctorOfClexBasisAndOrbits>
void for_clex_basis_and_orbits(
    std::shared_ptr<Structure const> const &shared_prim,
    ClexBasisSpecs const &clex_basis_specs, std::ostream &status,
    FunctorOfClexBasisAndOrbits f);

/// Print cluster with basis_index and nlist_index (from 0 to size()-1),
/// followed by cluster basis functions Functions are labeled \Phi_{i}, starting
/// from i = @param begin_ind Returns the number of functions that were printed
Index print_clust_basis(std::ostream &stream, BasisSet _clust_basis,
                        IntegralCluster const &_prototype, Index func_ind = 0,
                        int space = 18, char delim = '\n');

template <typename OrbitType>
void print_proto_clust_funcs(ClexBasis const &_clex_basis, std::ostream &out,
                             BasicStructure const &_prim,
                             std::vector<OrbitType> const &_tree);

namespace ClexBasis_impl {
std::vector<DoFKey> extract_dof_types(Structure const &_prim);

template <typename OrbitType>
BasisSet construct_proto_dof_basis(OrbitType const &_orbit,
                                   BasisSet::ArgList const &site_dof_sets,
                                   SymGroupRepID canonization_rep_ID);

// template<typename UCCIterType, typename IntegralClusterSymCompareType>
// std::map<UnitCellCoord, std::set<UnitCellCoord> > unique_ucc(UCCIterType
// begin, UCCIterType end, IntegralClusterSymCompareType const &sym_compare);
}  // namespace ClexBasis_impl

/// \brief Copy and apply SymOp to a Site.
namespace sym {
ClexBasis copy_apply(const SymOp &op, ClexBasis clex_basis);
}  // namespace sym
}  // namespace CASM

#endif
