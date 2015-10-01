#ifndef LATTICE_HH
#define LATTICE_HH

#include <iostream>
#include <cmath>

#include "casm/container/Array.hh"
#include "casm/container/LinearAlgebra.hh"
#include "casm/container/Counter.hh"

namespace CASM {
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  class SymGroup;
  class SymOp;

  class Lattice {
  private:
    Array<Vector3<double> > vecs;

    mutable Array<Vector3< double > > voronoi_table;

    //Coordinate Conversion Matrices
    //0 is fractional to cartesion; 1 is cartesian to fractional
    //use FRAC and CART globals to index
    //e.g., coord[CART]=coord_trans[FRAC]*coord[FRAC]; //converts frac to cart
    //Word to the wise: coord_trans[FRAC] is the matrix with columns equal to the lattice vectors
    Matrix3< double > coord_trans_mat[2];

    //int periodicity_dim;   //dimension of periodicity
    //int periodicity_axis;  //index of lattice vector that is non-periodic (2d) or periodic (1d)
  public:
    //Lattice properties
    //lengths -> lattice parameters
    //angles  -> angle between lattice vectors angle[i]=angle(get((i+1)%3), get((i+2)%3))
    Vector3<double> lengths, angles;


    //Lattice();

    Lattice(const Vector3<double> &vec1, const Vector3<double> &vec2,
            const Vector3<double> &vec3);

    ///Construct Lattice from a matrix of lattice vectors, where lattice vectors are columns
    ///(e.g., lat_mat is equivalent to coord_trans[FRAC])
    Lattice(const Matrix3<double> &lat_mat = Matrix3<double>::identity());
    Lattice(const Eigen::Matrix3d &lat_mat);

    Lattice(const Lattice &RHS);

    Lattice &operator=(const Lattice &RHS);
    const Vector3<double> &operator[](Index i)const {
      return vecs[i];
    };

    Lattice scaled_lattice(double scale);

    void calc_conversions();
    void calc_properties();

    double vol()const {
      return lat_column_mat().determinant();
    }

    //const acces of coord_trans Matrices. Should only be accessed using coord_trans(CART) or coord_trans(FRAC)
    const Matrix3< double > &coord_trans(int i)const {
      return coord_trans_mat[i];
    }

    ///Access coord_trans[FRAC] in an intuitive way.
    const Matrix3< double > &lat_column_mat() const {
      return coord_trans_mat[FRAC];
    };

    ///Access coord_trans[FRAC] in an intuitive way.
    const Matrix3< double > &inv_lat_column_mat() const {
      return coord_trans_mat[CART];
    };

    //Calculates the kpoint mesh for a supercell lattice given the kpoint mesh
    //for the primitive
    Array<int> calc_kpoints(Array<int> prim_kpoints, Lattice prim_lat);

    ///Return reciprocal lattice
    Lattice get_reciprocal() const;  //Anna did this

    /**pg_tol can be increased to find point group of lattice vectors
       that are slightly distorted due to numerical noise **/
    void generate_point_group(SymGroup &point_group, double pg_tol = TOL) const;

    void find_invariant_subgroup(const SymGroup &super_group, SymGroup &sub_group, double pg_tol = TOL) const;

    void generate_supercells(Array<Lattice> &supercell, const SymGroup &effective_pg, int max_prim_vol, int min_prim_vol = 1) const; //Donghee did this, ARN100113
    //void generate_supercells(Array<Lattice> &supercell, const MasterSymGroup &factor_group, int max_prim_vol, int min_prim_vol)const;

    template <typename T>
    Lattice make_supercell(const Matrix3<T> &trans_mat) const;

    ///Find the lattice vectors which give most compact unit cell
    Lattice get_reduced_cell() const;  //Alex did this

    /// populate voronoi information.
    //It should get called in any method of Lattice that tries to use an empty voronoi table
    void generate_voronoi_table() const;
    void print_voronoi_table(std::ostream &stream) const;

    /// Radius of largest sphere that totally fits inside the voronoi cell
    double min_voronoi_radius() const;

    /// returns voronoi vector that maximizes dot product with 'pos'
    Vector3<double> max_voronoi_vector(const Vector3<double> &pos) const;

    /// return number of voronoi cell faces that 'pos' is on, within +/- TOL
    /// 0 indicates that 'pos' is within the voronoi cell, -1 indicates that it is outside
    int voronoi_number(const Vector3<double> &pos) const;

    /**Gives a scaling of the lattice vectors such that after scaling,
       The eight parallelipipeds touching the origin enclose a sphere
       of given radius**/
    Vector3<int> enclose_sphere(double radius) const;

    template<typename CoordType, typename CoordType2>
    Array<CoordType> gridstruc_build(double max_radius, double min_radius, Array<CoordType> basis, CoordType2 lat_point); //Anirudh

    void read(std::istream &stream);
    void print(std::ostream &stream) const;

    void symmetrize(double _tol = TOL);

    //template <class T>
    //Lattice &operator*=(const Matrix3<T> &RHS);

    //Lattice &operator*=(const Eigen::Matrix3d &RHS);

    ///Are two lattices the same, even if they have different lattice vectors
    bool is_equivalent(const Lattice &RHS) const;

    ///Are lattice vectors identical for two lattices
    bool operator==(const Lattice &RHS) const;

    ///Matrix that relates two lattices (e.g., strain or slat)
    //Matrix3<double> operator/(const Lattice &RHS);

    //John G 121212
    ///Checks if lattice is a supercell of tile, acting on multiplication matrix. Check is performed applying operations from symlist
    bool is_supercell_of(const Lattice &tile, Matrix3<double> &multimat, double _tol=TOL) const;
    bool is_supercell_of(const Lattice &tile, const Array<SymOp> &symlist, Matrix3<double> &multimat, double _tol = TOL) const;

    ///Checks if lattice is a supercell of tile, applying operations from symlist
    bool is_supercell_of(const Lattice &tile, double _tol = TOL) const;
    bool is_supercell_of(const Lattice &tile, const Array<SymOp> &symlist, double _tol=TOL) const;

    /// Finds 'new_scel' equivalent to '*this' and 'new_prim' equivalent to 'prim', such that 'new_prim' perfectly tiles 'new_scel'
    /// Returns true if tesselation cannot be found.
    bool find_tessellation(Lattice &prim, Lattice &new_scel, Lattice &new_prim) const;

    ///Return a lattice with diagonal matrix that fits around starting lattice
    Lattice box(const Lattice &prim, const Lattice &scel, bool verbose = false) const;

    ///Finds smallest supercell of given list of supercells that has the same angles as the primitive cell. The value boxed defaults to 0 to make calculations faster
    void superduper_size_me(const Lattice &prim, const Array<Lattice> &supercells);

    ///Flip c vector if it's on the wrong side of a-b plane -- return (*this)
    Lattice &make_right_handed();
    ///Check if the lattice is right handed
    bool is_right_handed() const;
    ///Given a normal vector, a Vector3 containing the miller indeces for the lattice is generated
    Vector3< int > get_millers(Vector3< double > plane_normal, double tolerance = TOL) const;
    ///Generates a lattice with vectors a and b parallel to the plane described by the miller indeces

    Lattice get_lattice_in_plane(Vector3< int > millers, int max_vol = 20) const; //John G 121030

    Array<double> pg_converge(double large_tol);
    void pg_converge(double small_tol, double large_tol, double increment);

    ///Returns mirror image of lattice
    //Lattice get_reflection(bool override = 0) const;
    //\John G 121212

    //John G 011013
    void symmetrize(const SymGroup &relaxed_points);
    void symmetrize(const double &tolerance);

  };

  // write Lattice in json as array of vectors
  jsonParser &to_json(const Lattice &lat, jsonParser &json);
  void from_json(Lattice &lat, const jsonParser &json);


  /* never write a Matrix*Lattice operator, PLEASE

     template <class T>
     Lattice operator*(const Matrix3<T> &LHS, const Lattice &RHS);

     Lattice operator*(const Eigen::Matrix3d &LHS, const Lattice &RHS);
  */

  namespace niggli_impl {

    /// \brief Returns an equivalent Lattice in Niggli form, but without setting standard orientation
    Lattice _niggli(const Lattice &lat, double tol);
  }

  /// \brief Returns an equivalent Lattice in Niggli form with a standard orientation
  Lattice niggli(const Lattice &lat, const SymGroup &point_grp, double tol);

  /// \brief Rotate the Lattice to a standard orientation using allowed point group operations
  Lattice standard_orientation(const Lattice &lat, const SymGroup &point_grp, double tol);

  /// \brief Returns the volume of a Lattice
  double volume(const Lattice &lat);

  /// \brief Returns a super Lattice
  Lattice make_supercell(const Lattice &lat, const Eigen::Matrix3i &transf_mat);

  std::istream &operator>>(std::istream &in, const Lattice &lattice_in);

  //Implementation of template methods must live in *.hh file for proper compilation:

  //********************************************************************
  /**
   * This function generates a grid of points between max_radius and
   * min_radius. Additionally, it also fills up the points with a basis
   */
  //********************************************************************

  template<typename CoordType, typename CoordType2>
  Array<CoordType> Lattice::gridstruc_build(double max_radius, double min_radius, Array<CoordType> basis, CoordType2 lat_point) {
    Vector3<int> dim;
    dim = enclose_sphere(max_radius);
    Counter<Vector3<int> > grid_count(-dim, dim, Vector3<int>(1));
    double min_dist, dist;
    Array<CoordType> gridstruc;
    Vector3<int> temp;

    do {
      lat_point(FRAC) = grid_count();

      for(Index i = 0; i < basis.size(); i++) {
        CoordType tatom(basis[i] + lat_point);
        //get distance to closest basis site in the unit cell at the origin

        min_dist = 1e20;
        for(Index j = 0; j < basis.size(); j++) {
          dist = tatom.dist(basis[j]);
          if(dist < min_dist)
            min_dist = dist;
        }
        if(min_dist < min_radius) {
          continue;
        }
        if(min_dist < max_radius) {
          gridstruc.push_back(tatom);
          //          std::cout<<"tatom"<<tatom<<"\t Min Dist"<<min_dist<<"\n";
        }
      }
    }
    while(++grid_count);

    return gridstruc;
  }

  //********************************************************************
  // A column of trans_mat specifies a lattice vector of the supercell in terms of the
  // lattice vectors of (*this) lattice.
  template <typename T>
  Lattice Lattice::make_supercell(const Matrix3<T> &trans_mat) const {
    return Lattice(lat_column_mat() * trans_mat);
  }

  /// \brief Returns the volume of a Lattice
  ///
  /// \returns volume of the Lattice
  ///
  /// \param lat a \ref Lattice
  ///
  inline double volume(const Lattice &lat) {
    return lat[0].dot(lat[1].cross(lat[2]));
  }

  /// \brief Returns a super Lattice
  inline Lattice make_supercell(const Lattice &lat, const Eigen::Matrix3i &transf_mat) {
    return Lattice(Eigen::Matrix3d(lat.lat_column_mat()) * transf_mat.cast<double>());
  }

}
#endif
