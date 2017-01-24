#ifndef CLEXBASIS_HH
#define CLEXBASIS_HH

#include <string>
#include <vector>
#include "casm/container/multivector.hh"
#include "casm/basis_set/BasisSet.hh"

namespace CASM {

  class Structure;
  template <typename ClusterType> class CoordCluster;
  class UnitCellCoord;
  typedef CoordCluster<UnitCellCoord> IntegralCluster;

  class ClexBasis {

  public:

    typedef std::vector<BasisSet> BSetOrbit;
    typedef std::string DoFType;
    typedef std::vector<BSetOrbit>::const_iterator BSetOrbitIterator;

    /// \brief Initialize from Structure, in order to get Site DoF and global DoF info
    ClexBasis(Structure const &_prim) {}

    /// \brief Total number of BasisSet orbits
    Index n_orbits() const;

    /// \brief Total number of basis functions
    Index n_funcions() const;

    /// \brief Const access of clust basis of orbit @no and equivalent cluster @ne
    BasisSet const &clust_basis(Index no,
                                Index ne) const;

    /// \brief Const iterator to first BasisSet orbit
    BSetOrbitIterator begin() const {
      return m_bset_tree.cbegin();
    }

    /// \brief Const iterator to first BasisSet orbit
    BSetOrbitIterator cbegin() const {
      return m_bset_tree.cbegin();
    }

    /// \brief Const past-the-end iterator for BasisSet orbits
    BSetOrbitIterator end() const {
      return m_bset_tree.cbegin();
    }

    /// \brief Const past-the-end iterator for BasisSet orbits
    BSetOrbitIterator cend() const {
      return m_bset_tree.cbegin();
    }

    /// \brief generate clust_basis for all equivalent clusters in @param _orbitree
    template<typename ClusterOrbitIterator>
    void generate(ClusterOrbitIterator begin,
                  ClusterOrbitIterator end,
                  const jsonParser &bspecs,
                  std::vector<DoFType> const &dof_keys,
                  Index max_poly_order = -1);

  private:
    /// \brief Performs heavy lifting for populating site bases in m_site_bases
    void _populate_site_bases(Structure const &_prim);

    /// \brief Collection of all cluster BasisSets, one per cluster orbit
    std::vector<BSetOrbit> m_bset_tree;

    /// \brief Dictionary of all site BasisSets, initialized on construction
    std::map<DoFType, std::vector<BasisSet> > m_site_bases;

    /// \brief Dictionary of all global BasisSets, initialized
    std::map<DoFType, BasisSet> m_global_bases;
  };


  /// Print cluster with basis_index and nlist_index (from 0 to size()-1), followed by cluster basis functions
  /// Functions are labeled \Phi_{i}, starting from i = @param begin_ind
  void print_clust_basis(ClexBasis const &_basis_set,
                         std::ostream &stream,
                         Index begin_ind = 0,
                         int space = 18,
                         char delim = 0,
                         COORD_TYPE mode = COORD_DEFAULT);

  /// returns std::vector of std::string, each of which is
  std::vector<std::string> orbit_function_cpp_strings(ClexBasis const &_basis_set,
                                                      std::vector<FunctionVisitor *> const &labelers);

  /// nlist_index is the index into the nlist for the site the flower centers on
  std::vector<std::string> flower_function_cpp_strings(ClexBasis const &_basis_set,
                                                       std::vector<FunctionVisitor *> const &labelers,
                                                       Index nlist_index);

  /// b_index is the basis site index, f_index is the index of the configurational site basis function in Site::occupant_basis
  /// nlist_index is the index into the nlist for the site the flower centers on
  std::vector<std::string> delta_occfunc_flower_function_cpp_strings(ClexBasis const &_basis_set,
                                                                     BasisSet site_basis,
                                                                     const std::vector<FunctionVisitor *> &labelers,
                                                                     Index nlist_index,
                                                                     Index b_index,
                                                                     Index f_index);


  template<typename ClusterOrbitIterator>
  void ClexBasis::generate(ClusterOrbitIterator begin,
                           ClusterOrbitIterator end,
                           const jsonParser &bspecs,
                           std::vector<DoFType> const &dof_keys,
                           Index _max_poly_order /*= -1*/) {
    /*
    std::vector<DoFType> global_keys;
    std::vector<DoFType> local_keys;
    //separate local_args from global_args
    for(DoFType const &key : dof_keys) {
      if(m_global_bases.find(key) != m_global_bases.end()) {
        global_args.push_back(key);
      }
      else if(m_local_bases.find(key) != m_local_bases.end()) {
        local_keys.push_back
      }
      else {
        throw std::runtime_error(std::string("Attempting to build Clex basis set, but missing degree of freedom \"") + key + "\n");
      }
    }
    m_bset_tree.resize(std::distance(begin,end));
    Index l = 0;
    for(auto it=begin; it!=end; ++it) {
      m_bset_tree[l].reserve(it->size());
      m_bset_tree[l].push_back(generate_clust_basis(it->prototype(), local_keys, global_keys));
      for(Index k = 1; k < it->size(); k++) {
        m_bset_tree[l].push_back(_orbitree.orbit(i, j).equivalence_map[k][0]*m_bset_tree[l][0]);
      }
    }
    */
  }

  namespace ClexBasis_impl {
    void generate_clust_basis(multivector<BasisSet const *>::X<2> const &local_args,
                              std::vector<BasisSet const *> const &global_args,
                              Index max_poly_order = -1);


    BasisSet construct_clust_dof_basis(IntegralCluster const &_clust,
                                       std::vector<BasisSet const *> const &site_dof_sets);
  }
}

#endif
