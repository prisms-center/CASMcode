
#ifndef ORBITBRANCH_HH
#define ORBITBRANCH_HH

#include <iostream>
#include <fstream>

#include "casm/clusterography/Orbit.hh"
#include "casm/clusterography/SiteCluster.hh"

namespace CASM {
  class Lattice;

  class Supercell;

  template<typename ClustType>
  class GenericOrbitBranch;

  template<typename ClustType>
  class GenericOrbitree;

  typedef GenericOrbitBranch<SiteCluster> SiteOrbitBranch;
  typedef GenericOrbitree<SiteCluster> SiteOrbitree;

  //GenericOrbitBranch Class definition begins here:
  /// GenericOrbitBranch contains a set of GericOrbits that belong together (decided by user)
  /// with the constraint that all Orbits in an OrbitBranch have clusters with the same number of points
  /// i.e, an OrbitBranch contains only pairs, only triplets, only quadruplets, etc.
  template<typename ClustType>
  class GenericOrbitBranch : public Array< GenericOrbit<ClustType> > {
    ///The number of points for all Orbits in this OrbitBranch
    Index m_num_sites;
  public:

    using Array< GenericOrbit<ClustType> > :: size;
    using Array< GenericOrbit<ClustType> > :: at;
    using Array< GenericOrbit<ClustType> > :: back;

    /// Pivot is used for 'local' OrbitBranch; pivot is the cluster around which
    /// all clusters of the OrbitBranch are centered (i.e., the symmetry properties of pivot)
    ClustType pivot;

    ///'index' is index of each orbit in an ordered linear array
    Array<Index> index;

    GenericOrbitBranch(const Lattice &init_home) : m_num_sites(-1), pivot(init_home) { };
    GenericOrbitBranch(const Lattice &init_home, Index tnum_sites) : m_num_sites(tnum_sites), pivot(init_home) { };
    GenericOrbitBranch(const ClustType &init_pivot) : m_num_sites(-1), pivot(init_pivot) { };
    GenericOrbitBranch(const ClustType &init_pivot, Index tnum_sites) : m_num_sites(tnum_sites), pivot(init_pivot) { };

    void clear();

    ///Method to access orbits
    GenericOrbit<ClustType> &orbit(Index no);
    const GenericOrbit<ClustType> &orbit(Index no) const;

    ///How many points are allowed in clusters belonging to this OrbitBranch
    Index num_sites() const {
      return m_num_sites;
    };

    ///Method to access prototypes
    ClustType &prototype(Index no);
    const ClustType &prototype(Index no) const;

    ///Method to access equivalent clusters of Orbit 'no'
    const ClustType &equiv(Index no, Index ne) const;
    ClustType &equiv(Index no, Index ne);

    ///How many equivalent clusters are int orbit 'no'
    Index size(Index no) const;

    /// Calls set_lattice on all orbits of OrbitBranch
    void set_lattice(const Lattice &new_lat, COORD_TYPE mode);

    void set_pivot(const ClustType &new_pivot);

    /// Adds new orbit to OrbitBranch, but only if it has clusters of the correct size
    /// given by num_site(); if num_site() is undefined, push_back uses new_orbit to define its value
    void push_back(const GenericOrbit<ClustType> &new_orbit);

    void print(std::ostream &stream, COORD_TYPE mode = FRAC);

    ///Sorts all of the orbits in OrbitBranch by max_length
    void sort(); //Done - Alex

    // Alex did these
    ///If cluster is contained in OrbitBranch, return linear index of Orbit that contains it;
    ///else, return number of orbits in orbitree
    Index find(const ClustType &test_clust, double tol)const;

    ///If cluster exists in current OrbitBranch, return true
    bool contains(const ClustType &test_clust, double tol)const;

    ///apply_sym to everything in this OrbitBranch (i.e, pivot and all Orbits)
    GenericOrbitBranch &apply_sym(const SymOp &op);

    void generate_asymmetric_unit(const Array<typename ClustType::WhichCoordType > &basis, const SymGroup &factor_group, double tol);

    ///Finds Orbits of Clusters for which 'pivot' is a subcluster and adds them to 'flowerbranch'
    /// e.g., my_big_branch.extract_orbits_including(my_little_cluster, my_flower_branch);
    bool extract_orbits_including(const ClustType &pivot, GenericOrbitBranch &flowerbranch, double tol) const;

    //Translation operators need to be defined
    GenericOrbitBranch &operator+=(const Coordinate &shift) {
      return (*this);
    };
    GenericOrbitBranch &operator-=(const Coordinate &shift) {
      return (*this);
    };

    jsonParser &to_json(jsonParser &json) const;

    /// Assumes the pivot lattice is already set
    void from_json(const jsonParser &json);

  };

  // Some template weirdness is not allowing this to compile on flux.
  //   So I'm specializing for SiteCluster in ClusterFunctions.hh/.cc
  /*
  template<typename ClustType>
  jsonParser &to_json(const GenericOrbitBranch<ClustType> &branch, jsonParser &json) {
    return branch.to_json(json);
  }

  /// Assumes the prototype lattice is already set
  template<typename ClustType>
  void from_json(GenericOrbitBranch<ClustType> &branch, const jsonParser &json) {
    branch.from_json(json);
  }
  */
};

#include "casm/clusterography/OrbitBranch_impl.hh"

#endif
