#ifndef CLUSTER_HH
#define CLUSTER_HH

#include <iostream>
#include <complex>
#include <cmath>

#include "casm/symmetry/SymGroupRep.hh"
#include "casm/container/Tensor.hh"
#include "casm/crystallography/Site.hh"

namespace CASM {

  class Lattice;
  template < typename CoordType >
  class BasicStructure;
  class Structure;

  template < typename ClustType >
  class GenericOrbitree;

  template < typename CoordType >
  class GenericCluster;

  class ContinuousDoF;
  class DiscreteDoF;

  //class Site;

  //Begin Cluster class definition:
  /** A cluster is a group of points (CoordType) in the crystal
      A CoordType is any type that derives from Coordinate **/
  template < typename CoordType >
  class GenericCluster : public Array<CoordType> {
  private:
    Lattice const *home;
    double min_length_val, max_length_val;
    //DynamicVector<ClusterDoF> cluster_variables;

    //s2s_vec[i][j] describes vector pointing from site 'i' to site 'j'
    Array< Array<Coordinate> > s2s_vec;
    Array< Array<Coordinate> > s2s_norm_vec;

  public:
    using Array<CoordType> :: size;
    using Array<CoordType> :: at;
    using Array<CoordType> :: back;
    using Array<CoordType> :: operator[];
    using Array<CoordType> :: contains;
    using Array<CoordType> :: find;

    typedef CoordType WhichCoordType;

    Array<ContinuousDoF> clust_DoF;
    Index DoF_rep;

    SymGroup clust_group;
    SymGroupRep permute_group;

    TensorBasis<double> tensor_basis;

    Tensor<double> eci;

    //GenericCluster() : min_length_val(0.0), max_length_val(0.0), clust_group(LOCAL) {};
    GenericCluster(const Lattice &init_home);

    void set_lattice(const Lattice &new_home, COORD_TYPE mode);
    const Lattice &get_home() const {
      return *home;
    };

    virtual void push_back(const CoordType &new_coord);

    ///Translate entire cluster so that point at(pivot_ind) is inside unit cell
    void within(Index pivot_ind = 0);
    /**Translate entire cluster so that point at(pivot_ind) is inside unit cell
       Coordinate trans contains translation used to map within**/
    void within(Index pivot_ind, Coordinate &trans);

    ///Map every point of cluster inside unit cell
    void all_within();

    void update_data_members(const BasicStructure<CoordType> &ref_struc);

    /**Checks to see if cluster is "compact" relative to (Lattice cell)
       in other words, period images of the cluster points are farther away than
       the points themselves. Returns true if cluster is not "compact"  **/
    bool image_check(const Lattice &cell, int nV = 0) const;

    ReturnArray<int> SuperScreener(Array<Lattice> &s_cells, const Array<SymOp> &symoplist);

    /// permute sites of the cluster, and everything that depends on the site order
    GenericCluster &permute(const Array<Index> &perm);
    GenericCluster &permute(const Permutation &perm);

    ///apply symmetry to all points of the cluster
    GenericCluster &apply_sym(const SymOp &op);
    ///apply symmetry to all points of the cluster without translation
    GenericCluster &apply_sym_no_trans(const SymOp &op);

    ///Finds the sub_group of super_group that is the point group of the cluster
    void get_clust_group(const SymGroup &super_group);   //calc_clust_symmetry in old code

    ///Uses clust_group to populate permute_group (i.e., how each operation permutes the cluster points)
    void get_permute_group();

    /// Uses clust_group to population tensor_basis; tensors are of rank 'nrank'
    void get_tensor_basis(Index nrank);

    /// gets phase factor for the vector from site i to site j, at k-point 'k'
    std::complex<double> get_phase(const Coordinate &k, int i = 0, int j = 1);

    /// populate s2s_vec
    void get_s2s_vec();
    Coordinate s2s(int i, int j) const {
      return s2s_vec[i][j];
    };

    /// gets max_length and min_length
    void calc_properties();  //Alex do this
    void calc_properties(GenericCluster<CoordType> phenom_clust);

    ///Call CoordType::update on each of the sites
    void update();
    double max_length() const {
      return max_length_val;
    }
    double min_length() const {
      return min_length_val;
    }

    ///Returns the geometric center of "mass" of a cluster (treats all sites as having equal mass)
    Coordinate geometric_center() const;

    ///Performs preparatory steps on prototype before doing Orbit::get_equivalent()
    void prepare_prototype() {};

    //Anna do this
    /// is test_cluster a subcluster of (*this)
    bool contains(const GenericCluster &test_cluster) const;

    // Like Array<CoordType>::contains(), but takes periodicity mode into account
    bool contains_periodic(const CoordType &test_coord) const;

    /** is test_cluster a subcluster of (*this), and how do the indices map
    'index' is populated with the indices of (*this) that correspond to the
    points of test_cluster **/
    bool find(const GenericCluster &test_cluster, Array<Index> &index) const;

    /** if pivot is a sub_cluster, return true and translate (*this) by a lattice translation
    so that the points of 'pivot' are coincident with subcluster points in (*this) **/
    bool map_onto_subcluster(const GenericCluster &pivot);
    bool map_onto_subcluster(const GenericCluster &pivot, int num_maps);

    //Only use for cluster of sites -- get overlaps of orbitals of site 0 with orbitals of site 1
    // This is tight-binding voodoo
    void get_overlap_tensor();

    ///Figure out which basis atoms in basis correspond to the points in cluster (*this)
    void collect_basis_info(const Array<CoordType> &basis);
    /**Figure out which basis atoms in basis correspond to the points in cluster (*this)
       when cluster is translated by 'shift' **/
    void collect_basis_info(const Array<CoordType> &basis, const Coordinate &shift);


    void read(std::istream &stream, int num_sites, COORD_TYPE mode, bool SD_is_on);
    /// Reads the cluster and maybe the tensor basis and eci tensor
    void read(std::istream &stream, COORD_TYPE mode, bool read_tensors); //Modified by Ivy


    //Alex do this
    void print(std::ostream &stream, char delim = 0, COORD_TYPE mode = COORD_DEFAULT) const;
    void print_shifted(std::ostream &stream, const Coordinate &shift, char delim = 0, COORD_TYPE mode = COORD_DEFAULT) const;
    void print_sites(std::ostream &stream, int space, char delim = 0, COORD_TYPE mode = COORD_DEFAULT) const;
    void print_basis_info(std::ostream &stream, int space, char delim = 0, COORD_TYPE mode = COORD_DEFAULT) const;
    void print_decorated_sites(std::ostream &stream, int space, char delim = 0, COORD_TYPE mode = COORD_DEFAULT) const;

    ///adds unique points of 'RHS' to (*this)
    void merge(const GenericCluster &RHS);
    ///Adds new point to cluster, but only if it is unique
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
    /** are two clusters identical, to within permutation and translation
    translation that maps clusters is stored in 'trans' **/
    bool is_equivalent(const GenericCluster &test_clust, Coordinate &trans) const;

    /// if is_equivalent(test_clust) is true, return true and map (*this) onto test_clust
    bool map_onto(const GenericCluster &test_clust);

    /** if is_equivalent(test_clust) is true, return true and map (*this) onto test_clust
    translation that maps clusters is stored in 'trans' **/
    bool map_onto(const GenericCluster &test_clust, Coordinate &trans);

    /// Generate an orbitree that contains all the subclusters of this cluster
    //    GenericOrbitree< GenericCluster<CoordType> > enumerate_subclusters(const Structure &prim, bool verbose) const;

    /// Write GenericCluster to json. Does not write lattice.
    jsonParser &to_json(jsonParser &json) const;
    void from_json(const jsonParser &json);
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


  /// Write GenericCluster to json. Does not write lattice.
  template<typename CoordType>
  jsonParser &GenericCluster<CoordType>::to_json(jsonParser &json) const {

    json.put_obj();

    // members not included:
    //   Lattice const *home;
    //   Array< Array<Coordinate> > s2s_vec;
    //   Array< Array<Coordinate> > s2s_norm_vec;


    // inherits: public Array<CoordType>
    const Array<CoordType> &coord_array_ref = *this;
    // std::cout << "coord_array_ref:" << coord_array_ref.size() << std::endl;
    // std::cout << coord_array_ref << std::endl;
    json["coordtype"] = coord_array_ref;

    // members included:

    // double min_length_val, max_length_val;
    json["min_length"] = min_length_val;
    json["max_length"] = max_length_val;

    // CASM::Array<CASM::DoF *> clust_DoF;
    if(clust_DoF.size() > 0) {
      json["clust_DoF"] = clust_DoF;
    }

    // int DoF_rep;
    json["DoF_rep"] = DoF_rep;

    //Turning off the printing of cluster symmetry
    //TODO: Come up with some way of enabling this
    // // SymGroup clust_group;
    // if(clust_group.size() > 0)
    //   json["clust_group"] = clust_group;

    // // SymGroupRep permute_group;
    // if(permute_group.size() > 0)
    //   json["permute_group"] = permute_group;

    // TensorBasis<double> tensor_basis;
    if(tensor_basis.size() > 0)
      json["tensor_basis"] = tensor_basis;

    // Tensor<double> eci;
    if(eci.size() > 0)
      json["eci"] = eci;

    return json;
  }

  // Read from json. Lattice must be set previously. Assumes CoordType::CoordType(Lattice) exists.
  template<typename CoordType>
  void GenericCluster<CoordType>::from_json(const jsonParser &json) {
    try {

      // members not read:
      //   Lattice const *home;
      //   Array< Array<Coordinate> > s2s_vec;
      //   Array< Array<Coordinate> > s2s_norm_vec;

      // inherits: public Array<CoordType>

      //std::cout<<"Number of coordinates: "<<json["coordtype"].size()<<std::endl;
      //this->resize(json["coordtype"].size(), coord);
      for(int i = 0; i < json["coordtype"].size(); i++) {
        CoordType coord(*home);
        CASM::from_json(coord, json["coordtype"][i]);
        (*this).push_back(coord);
      }
      //std::cout<<"Read in the coordinates"<<std::endl;

      // double min_length_val, max_length_val;
      //   Read these, because if local clusters they were generated from phenom_clust
      CASM::from_json(min_length_val, json["min_length"]);
      CASM::from_json(max_length_val, json["max_length"]);
      //std::cout<<"Read in the min and max lengths"<<std::endl;


      clust_DoF.clear();
      //std::cout<<"Reading the clust_DoF"<<std::endl;
      if(json.contains("clust_DoF")) {
        clust_DoF.resize(json["clust_DoF"].size());
        for(int i = 0; i < json["clust_DoF"].size(); i++) {

          // This allocates a new object to clust_DoF[i].
          //   It needs a Lattice in case it is a MoleculeOccupant
          CASM::from_json(clust_DoF[i], json["clust_DoF"][i]);
        }
      }

      //std::cout<<"Reading the DoF_rep"<<std::endl;
      // int DoF_rep;
      CASM::from_json(DoF_rep, json["DoF_rep"]);

      // // SymGroup clust_group;
      // //This is a hack to initialize the appropriate SymOp
      // Coordinate temp_coord(Vector3<double>(0.0,0.0,0.0),*home);
      // SymOp temp_symOp(temp_coord);
      // if(clust_group.size()!=0)
      //   std::cerr<<"WARNING in Cluster::from_json. Your clust_group "
      //            <<"is not empty. Clearing anyways"<<std::endl;
      // clust_group.clear();
      // clust_group.push_back(temp_symOp);
      // //std::cout<<"Reading the clust_group"<<std::endl;
      // if(json.contains("clust_group"))
      //   CASM::from_json(clust_group, json["clust_group"]);

      // //std::cout<<"Reading the permute_group"<<std::endl;
      // // SymGroupRep permute_group;
      // if(json.contains("permute_group"))
      //   CASM::from_json(permute_group, json["permute_group"], *home);

      //std::cout<<"Reading the tensor_basis"<<std::endl;
      // TensorBasis<double> tensor_basis;
      if(json.contains("tensor_basis"))
        CASM::from_json(tensor_basis, json["tensor_basis"]);

      //std::cout<<"Reading the eci"<<std::endl;
      // Tensor<double> eci;
      if(json.contains("eci"))
        CASM::from_json(eci, json["eci"]);

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
