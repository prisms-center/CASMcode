#ifndef ORBIT_HH
#define ORBIT_HH

#include <iostream>

#include "casm/clusterography/Cluster.hh"

namespace CASM {

  template<typename CoordType> class BasicStructure;
  class Site;

  template<typename ClustType>
  class GenericOrbit;

  class SiteCluster;
  typedef GenericOrbit<SiteCluster> SiteOrbit;

  class HopCluster;
  typedef GenericOrbit<HopCluster> HopOrbit;

  class FunctionVisitor;

  class BasisSet;

  //GenericOrbit class definition begins here:
  /**GenericOrbit is the set of all clusters that are
     equivalent with respect to a SymGroup and have at least one point in the unit cell.
     If two clusters can be mapped onto each other by translation
     only one of them will be included in orbit (typically) **/
  template<typename ClustType>
  class GenericOrbit : public Array<ClustType > {
  private:
    ///Linear index of this orbit, when many orbits are stored in a complicated structure (e.g., Orbitree)
    mutable Index index;
    Index permute_rep_ID, coord_rep_ID;
    ///Pointers to all orbits whose clusters are subclusters of clusters in this orbit
    Array<GenericOrbit *> sub_cluster;

  public:
    using Array<ClustType > :: size;
    using Array<ClustType > :: at;
    using Array<ClustType > :: clear;

    ///One example of the clusters contained in this orbit
    ///All clusters in the Orbit are equivalent to 'prototype' by symmetry
    ClustType  prototype;

    ///equivalence_map[i][j] is a symmetry operation that maps 'prototype' onto cluster 'i'
    /// 'j' indexes ALL the symmetry operations that perform this map
    Array< Array<SymOp> > equivalence_map;

    GenericOrbit(const ClustType  &init_prototype) : index(-1), permute_rep_ID(-1), coord_rep_ID(-1), prototype(init_prototype) {};

    ///calls set_lattice on prototype, and all equivalent clusters
    void set_lattice(const Lattice &new_home, COORD_TYPE mode);

    ///Take prototype cluster, find all of the equivalent clusters (fills orbit array)
    ///only translates if the periodicity flag is on
    ///Also, fill equivalence_map
    ///example:  GenericOrbit<SiteCluster> my_orbit(my_cluster); //specifies prototype
    ///          my_orbit.get_equivalent(my_point_group);      //maps protype onto all its equivalents
    void get_equivalent(const SymGroup &sym_group, double tol);

    ///Apply symmetry to prototype and all the clusters in the orbit
    GenericOrbit &apply_sym(const SymOp &op);

    ///go through all cluster in orbit array and see if test_clust is among them,
    /// to within a lattice translation
    ///only translates if the periodicity flag is on
    bool contains(const ClustType  &test_clust, double tol) const;
    /// check if test_site is contained in any of the equivalent clusters,
    /// to within a lattice translation
    /// only translates if the periodicity flag is on
    bool contains(const typename ClustType::WhichCoordType &test_site, double tol) const;

    /// Same as contains, but returns index of equivalent cluster that maps onto test_clust
    Index find(const ClustType  &test_clust, double tol) const;
    /// Same as contains, but returns index of equivalent cluster that maps onto test_clust
    /// translation that maps equivalent cluster onto test_clust is stored in 'trans'
    Index find(const ClustType  &test_clust, Coordinate &trans, double tol) const;

    /// calls collect_basis_info on all clusters in orbit
    void collect_basis_info(const Array<typename ClustType::WhichCoordType> &basis, const Coordinate &shift);
    /// calls collect_basis_info on all clusters in orbit, with respect to a translation
    void collect_basis_info(const Array<typename ClustType::WhichCoordType> &basis);

    /// get permutation representation of every operation in equivalence_map to describe how operations permute site order
    SymGroupRep const *get_full_permutation_representation();

    /// get symmetry representation of every operation in equivalence_map to describe how operations map coordinates at site 'i' of prototype
    /// onto coordinates of site 'j' of equivalent cluster
    SymGroupRep const *get_full_coord_representation();

    /// return max_length of clusters in Orbit
    double max_length() const {
      return prototype.max_length();
    }
    /// return min_length of clusters in Orbit
    double min_length() const {
      return prototype.min_length();
    }

    //Not implemented
    void read(std::istream &stream, int num_sites, COORD_TYPE mode);
    /// reads in an Orbit
    void read(std::istream &stream, COORD_TYPE mode, const SymGroup &sym_group); //Modified by Ivy


    /// returns Array of std::string, each of which is
    ReturnArray<std::string> orbit_function_cpp_strings(const Array<FunctionVisitor *> &labelers);

    /// nlist_index is the index into the nlist for the site the flower centers on
    ReturnArray<std::string> flower_function_cpp_strings(const Array<FunctionVisitor *> &labelers, Index nlist_index);

    /// b_index is the basis site index, f_index is the index of the configurational site basis function in Site::occupant_basis
    /// nlist_index is the index into the nlist for the site the flower centers on
    ReturnArray<std::string> delta_occfunc_flower_function_cpp_strings(BasisSet site_basis, const Array<FunctionVisitor *> &labelers, Index nlist_index, Index b_index, Index f_index);

    //Access index, which is private
    Index get_index() const {
      return this->index;
    };

    //Allows you to set the index since it's private
    void set_index(Index ind) const { //Added by Ivy 05/05/2013 -- This may be a temporary fix depending on how we want to treat index
      this->index = ind;
    };

    jsonParser &to_json(jsonParser &json) const;

    /// Assumes the prototype lattice is already set
    void from_json(const jsonParser &json);

  };

  // Some template weirdness is not allowing this to compile on flux.
  //   So I'm specializing for SiteCluster in ClusterFunctions.hh/.cc
  /*
  template<typename ClustType>
  jsonParser &to_json(const GenericOrbit<ClustType> &orbit, jsonParser &json) {
    return orbit.to_json(json);
  }

  /// Assumes the prototype lattice is already set
  template<typename ClustType>
  void from_json(GenericOrbit<ClustType> &orbit, const jsonParser &json) {
    orbit.from_json(json);
  }
  */


};

#include "casm/clusterography/Orbit_impl.hh"

#endif

