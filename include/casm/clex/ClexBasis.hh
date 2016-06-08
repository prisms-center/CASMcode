#ifndef CLEXBASIS_HH
#define CLEXBASIS_HH

#include <string>
#include "casm/basis_set/BasisSet.hh"

class SiteOrbitree;

namespace CASM {
  class BasisBuilder;

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

    /// \brief Const access of clust basis of orbit @param orbit_ind and equivalent cluster @param equiv_ind
    BasisSet const &clust_basis(Index orbit_ind,
                                Index equiv_ind) const;

    /// \brief Const access of BSetOrbit of orbit @param orbit_ind
    BSetOrbit const &orbit_basis(Index orbit_ind) const;

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
    template<typename OrbitIterType>
    void generate(OrbitIterType _begin,
                  OrbitIterType _end,
                  std::vector<DoFType> const &dof_keys,
                  Index max_poly_order = -1);

  private:
    template<typename OrbitType>
    BasisSet _construct_prototype_basis(OrbitType const &_orbit,
                                        std::vector<DoFType> const &local_keys,
                                        std::vector<DoFType> const &global_keys,
                                        Index max_poly_order) const;

    /// \brief Performs heavy lifting for populating site bases in m_site_bases
    void _populate_site_bases(Structure const &_prim);

    notstd::cloneable_ptr<BasisBuilder> m_builder;

    /// \brief Collection of all cluster BasisSets, one per cluster orbit
    std::vector<BSetOrbit> m_bset_tree;

    /// \brief Dictionary of all site BasisSets, initialized on construction
    std::map<DoFType, std::vector<BasisSet> > m_site_bases;

    /// \brief Dictionary of all global BasisSets, initialized
    std::map<DoFType, BasisSet> m_global_bases;
  };


  class BasisBuilder {
  public:
    virtual ~BasisBuilder() {}

    virtual void prepare(Structure const &_prim) {

    }

    virtual std::vector<ClexBasis::DoFType> filter_dof_types(std::vector<ClexBasis::DoFType> const &_dof_types) {
      return _dof_types;
    }

    virtual void pre_generate() {

    }

    std::unique_ptr<BasisBuilder> clone()const {
      return std::unique_ptr<BasisBuilder>(_clone());
    }
  private:
    virtual BasisBuilder *_clone()const = 0;

  };

  /// Print cluster with basis_index and nlist_index (from 0 to size()-1), followed by cluster basis functions
  /// Functions are labeled \Phi_{i}, starting from i = @param begin_ind
  void print_clust_basis(ClexBasis const &_basis_set,
                         PrimNeighborList &_nlist,
                         std::ostream &_stream,
                         Index begin_ind = 0,
                         int space = 18,
                         char delim = 0,
                         COORD_TYPE mode = COORD_DEFAULT);

  /// returns std::vector of std::string, each of which is
  template<typename OrbitType>
  std::vector<std::string> orbit_function_cpp_strings(ClexBasis::BSetOrbit _bset_orbit, // used as temporary
                                                      OrbitType const &_clust_orbit,
                                                      PrimNeighborList &_nlist,
                                                      std::vector<FunctionVisitor *> const &labelers);

  /// nlist_index is the index into the nlist for the site the flower centers on
  template<typename OrbitType>
  std::vector<std::string> flower_function_cpp_strings(ClexBasis::BSetOrbit _bset_orbit, // used as temporary
                                                       OrbitType const &_clust_orbit,
                                                       PrimNeighborList &_nlist,
                                                       std::vector<FunctionVisitor *> const &labelers,
                                                       Index sublat_index);

  /// b_index is the basis site index, f_index is the index of the configurational site basis function in Site::occupant_basis
  /// nlist_index is the index into the nlist for the site the flower centers on
  template<typename OrbitType>
  std::map< UnitCell, std::vector< std::string > > delta_occfunc_flower_function_cpp_strings(ClexBasis::BSetOrbit _bset_orbit, // used as temporary
      OrbitType const &_clust_orbit,
      PrimNeighborList &_nlist,
      BasisSet site_basis,
      const std::vector<FunctionVisitor *> &labelers,
      Index nlist_index,
      Index b_index,
      Index f_index);


  namespace ClexBasis_impl {
    void generate_clust_basis(multivector<BasisSet const *>::X<2> const &local_args,
                              std::vector<BasisSet const *> const &global_args,
                              Index max_poly_order = -1);


    BasisSet construct_clust_dof_basis(SiteCluster const &_clust,
                                       std::vector<BasisSet const *> const &site_dof_sets);
  }
}
#include "casm/clex/ClexBasis_impl.hh"
#endif
