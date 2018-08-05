#ifndef CLUSTER_HH
#define CLUSTER_HH

#include <iostream>
#include <complex>
#include <cmath>

#include "casm/crystallography/Lattice.hh"
#include "casm/crystallography/Coordinate.hh"
#include "casm/symmetry/SymGroupRep.hh"

namespace CASM {

  class Lattice;
  template < typename CoordType >
  class BasicStructure;

  template < typename ClustType >
  class GenericOrbitree;

  template < typename CoordType >
  class GenericCluster;

  //Begin Cluster class definition:
  /// A cluster is a group of points (CoordType) in the crystal
  /// A CoordType is any type that derives from Coordinate
  template < typename CoordType >
  class GenericCluster : public Array<CoordType> {
  public:
    using Array<CoordType> :: size;
    using Array<CoordType> :: at;
    using Array<CoordType> :: back;
    using Array<CoordType> :: operator[];
    using Array<CoordType> :: contains;
    using Array<CoordType> :: find;

    typedef CoordType WhichCoordType;

    //GenericCluster() : m_min_length(0.0), m_max_length(0.0), clust_group(LOCAL) {}
    GenericCluster(const Lattice &init_home);

    const Lattice &home() const {
      return *m_lat_ptr;
    }

    const SymGroup &clust_group() const {
      return m_clust_group;
    }

    const SymGroupRep::RemoteHandle &permute_rep() const {
      return m_permute_rep;
    }

    void set_clust_group(const Array<SymOp> &new_group) {
      _clust_group() = new_group;
    }

    void set_permute_rep(SymGroupRepID perm_rep_ID) {
      if(clust_group().size() == 0) {
        std::cerr << "CRITICAL ERROR: In GenericCluster<>::set_permute_rep(), Cluster::clust_group() has not yet been initialized!\n"
                  << "                Groups must always be initialized before their representations. Exiting...\n";
        assert(0);
        exit(1);
      }
      _permute_rep().set_rep(clust_group(), perm_rep_ID);
    }

    void set_lattice(const Lattice &new_home, COORD_TYPE mode);

    void push_back(const CoordType &new_coord);

    ///Translate entire cluster so that point at(pivot_ind) is inside unit cell
    void within(Index pivot_ind = 0);

    /// \brief Translate entire cluster so that point at(pivot_ind) is inside unit cell
    ///  @param Coordinate trans contains translation used to map within
    void within(Index pivot_ind, Coordinate &trans);

    ///Map every point of cluster inside unit cell
    void all_within();

    void update_data_members(const BasicStructure<CoordType> &ref_struc);

    /// \brief Checks to see if cluster is "compact" relative to (Lattice cell)
    /// in other words, period images of the cluster points are farther away than
    /// the points themselves. Returns true if cluster is not "compact"  **/
    bool image_check(const Lattice &cell, int nV = 0) const;

    /// permute sites of the cluster, and everything that depends on the site order
    GenericCluster &permute(const Array<Index> &perm);
    GenericCluster &permute(const Permutation &perm);

    /// apply symmetry to all points of the cluster
    GenericCluster &apply_sym(const SymOp &op);
    /// apply symmetry to all points of the cluster without translation
    GenericCluster &apply_sym_no_trans(const SymOp &op);

    /// Finds the sub_group of super_group that is the point group of the cluster
    void generate_clust_group(const SymGroup &super_group, std::vector<Permutation> *perm_array_ptr = nullptr, double tol = TOL);

    /// Finds the Permutation corresponding to each element of clust_group
    std::vector<Permutation> clust_group_permutations(double tol) const;

    /// gets max_length and min_length
    void calc_properties();  //Alex do this
    void calc_properties(GenericCluster<CoordType> phenom_clust);

    double max_length() const {
      return m_max_length;
    }

    double min_length() const {
      return m_min_length;
    }

    ///Returns the geometric center of "mass" of a cluster (treats all sites as having equal mass)
    Coordinate geometric_center() const;

    ///Performs preparatory steps on prototype before doing Orbit::get_equivalent()
    void prepare_prototype() {}

    /// is test_cluster a subcluster of (*this)
    bool contains(const GenericCluster &test_cluster) const;

    /// Like Array<CoordType>::contains(), but takes periodicity mode into account
    bool contains_periodic(const CoordType &test_coord, double tol) const;

    /// \brief is test_cluster a subcluster of (*this), and how do the indices map
    /// points of test_cluster
    Index find(const CoordType &test_elem, double tol) const;

    /// \brief is test_cluster a subcluster of (*this), and how do the indices map
    /// 'index' is populated with the indices of (*this) that correspond to the
    /// points of test_cluster
    bool find(const GenericCluster &test_cluster, Array<Index> &index, double tol) const;

    /// if pivot is a sub_cluster, return true and translate (*this) by a lattice translation
    /// so that the points of 'pivot' are coincident with subcluster points in (*this)
    bool map_onto_subcluster(const GenericCluster &pivot, double tol = TOL);
    bool map_onto_subcluster(const GenericCluster &pivot, int num_maps, double tol = TOL);

    /// Figure out which basis atoms in basis correspond to the points in cluster (*this)
    void collect_basis_info(const Array<CoordType> &basis);
    /// Figure out which basis atoms in basis correspond to the points in cluster (*this)
    /// when cluster is translated by 'shift'
    void collect_basis_info(const Array<CoordType> &basis, const Coordinate &shift);


    void read(std::istream &stream, int num_sites, COORD_TYPE mode, bool SD_is_on);
    /// Reads the cluster
    void read(std::istream &stream, COORD_TYPE mode); //Modified by Ivy

    void print(std::ostream &stream, char delim = '\n', COORD_TYPE mode = COORD_DEFAULT) const;
    void print_shifted(std::ostream &stream, const Coordinate &shift, char delim = '\n', COORD_TYPE mode = COORD_DEFAULT) const;
    void print_sites(std::ostream &stream, int space, char delim = '\n', COORD_TYPE mode = COORD_DEFAULT) const;
    void print_basis_info(std::ostream &stream, int space, char delim = '\n', COORD_TYPE mode = COORD_DEFAULT) const;
    void print_decorated_sites(std::ostream &stream, int space, char delim = '\n', COORD_TYPE mode = COORD_DEFAULT) const;

    /// adds unique points of 'RHS' to (*this)
    void merge(const GenericCluster &RHS);
    /// Adds new point to cluster, but only if it is unique
    void merge(const CoordType &RHS);

    /// in=place translation of a cluster
    GenericCluster &operator+=(const Coordinate &RHS);
    GenericCluster &operator-=(const Coordinate &RHS);

    /// create translated cluster
    GenericCluster operator+(const GenericCluster &RHS);

    /// are two clusters identical, to within a permutation
    bool operator==(const GenericCluster &RHS) const;

    /// are two clusters identical, to within permutation and translation
    bool is_equivalent(const GenericCluster &test_clust) const;

    /// \brief are two clusters identical, to within permutation and translation
    /// translation that maps clusters is stored in 'trans'
    bool is_equivalent(const GenericCluster &test_clust, Coordinate &trans) const;

    /// if is_equivalent(test_clust) is true, return true and map (*this) onto test_clust
    bool map_onto(const GenericCluster &test_clust, double tol);

    /// if is_equivalent(test_clust) is true, return true and map (*this) onto test_clust
    /// translation that maps clusters is stored in 'trans'
    bool map_onto(const GenericCluster &test_clust, Coordinate &trans, double tol);

    /// Write GenericCluster to json. Does not write lattice.
    jsonParser &to_json(jsonParser &json) const;
    void from_json(const jsonParser &json);

  protected:
    SymGroup &_clust_group() {
      return m_clust_group;
    }

    SymGroupRep::RemoteHandle &_permute_rep() {
      return m_permute_rep;
    }

  private:
    Lattice const *m_lat_ptr;
    double m_min_length, m_max_length;
    SymGroup m_clust_group;
    SymGroupRep::RemoteHandle m_permute_rep;
  };

  ///create symmetry-transformed cluster
  template <typename CoordType>
  GenericCluster<CoordType> operator*(const SymOp &LHS, const GenericCluster<CoordType> &RHS);

  ///create translated cluster
  template <typename CoordType>
  GenericCluster<CoordType> operator+(const GenericCluster<CoordType> &LHS, const Coordinate &RHS);

  ///create translated cluster
  template <typename CoordType>
  GenericCluster<CoordType> operator-(const GenericCluster<CoordType> &LHS, const Coordinate &RHS);

  template <typename CoordType>
  bool almost_equal(const GenericCluster<CoordType> &LHS, const GenericCluster<CoordType> &RHS, double tol);



  /// Write GenericCluster to json. Does not write lattice.
  template<typename CoordType>
  jsonParser &GenericCluster<CoordType>::to_json(jsonParser &json) const {

    json.put_obj();

    // members not included:
    //   Lattice const *m_lat_ptr;


    // inherits: public Array<CoordType>
    const Array<CoordType> &coord_array_ref = *this;
    // std::cout << "coord_array_ref:" << coord_array_ref.size() << std::endl;
    // std::cout << coord_array_ref << std::endl;
    json["coordtype"] = coord_array_ref;

    // members included:

    // double m_min_length, m_max_length;
    json["min_length"] = m_min_length;
    json["max_length"] = m_max_length;

    //Turning off the printing of cluster symmetry
    //TODO: Come up with some way of enabling this
    // // SymGroup clust_group;
    // if(clust_group.size() > 0)
    //   json["clust_group"] = clust_group;

    return json;
  }

  // Read from json. Lattice must be set previously. Assumes CoordType::CoordType(Lattice) exists.
  template<typename CoordType>
  void GenericCluster<CoordType>::from_json(const jsonParser &json) {
    try {

      // members not read:
      //   Lattice const *m_lat_ptr;

      // inherits: public Array<CoordType>

      //std::cout<<"Number of coordinates: "<<json["coordtype"].size()<<std::endl;
      //this->resize(json["coordtype"].size(), coord);
      for(int i = 0; i < json["coordtype"].size(); i++) {
        CoordType coord(home());
        CASM::from_json(coord, json["coordtype"][i]);
        (*this).push_back(coord);
      }
      //std::cout<<"Read in the coordinates"<<std::endl;

      // double m_min_length, m_max_length;
      //   Read these, because if local clusters they were generated from phenom_clust
      CASM::from_json(m_min_length, json["min_length"]);
      CASM::from_json(m_max_length, json["max_length"]);
      //std::cout<<"Read in the min and max lengths"<<std::endl;


      // // SymGroup clust_group;
      // //This is a hack to initialize the appropriate SymOp
      // Coordinate temp_coord(Vector3<double>(0.0,0.0,0.0),home());
      // SymOp temp_symOp(temp_coord);
      // if(clust_group.size()!=0)
      //   std::cerr<<"WARNING in Cluster::from_json. Your clust_group "
      //            <<"is not empty. Clearing anyways"<<std::endl;
      // clust_group.clear();
      // clust_group.push_back(temp_symOp);
      // //std::cout<<"Reading the clust_group"<<std::endl;
      // if(json.contains("clust_group"))
      //   CASM::from_json(clust_group, json["clust_group"]);

    }
    catch(...) {
      /// re-throw exceptions
      throw;
    }
  }


  template<typename T>
  jsonParser &to_json(const GenericCluster<T> &clust, jsonParser &json) {
    return clust.to_json(json);
  }

  template<typename T>
  void from_json(GenericCluster<T> &clust, const jsonParser &json) {
    clust.from_json(json);
  }

};

#include "casm/clusterography/Cluster_impl.hh"

#endif
