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
  public:
    typedef Eigen::Matrix3d::ConstColXpr LatVec;

    Lattice(const Eigen::Vector3d &vec1, const Eigen::Vector3d &vec2,
            const Eigen::Vector3d &vec3);

    ///Construct Lattice from a matrix of lattice vectors, where lattice vectors are columns
    ///(e.g., lat_mat is equivalent to coord_trans[FRAC])
    Lattice(const Eigen::Matrix3d &lat_mat = Eigen::Matrix3d::Identity());

    /// \brief Construct FCC primitive cell of unit volume
    static Lattice fcc();

    /// \brief Construct BCC primitive cell of unit volume
    static Lattice bcc();

    /// \brief Construct simple cubic primitive cell of unit volume
    static Lattice cubic();

    /// \brief Construct cubic primitive cell of unit volume
    static Lattice hexagonal();

    LatVec operator[](Index i)const {
      return m_lat_mat.col(i);
    }

    Lattice scaled_lattice(double scale) const;

    double angle(Index i)const;

    double length(Index i)const;

    double vol()const {
      return lat_column_mat().determinant();
    }

    ///Access coord_trans[FRAC] in an intuitive way.
    const Eigen::Matrix3d &lat_column_mat() const {
      return m_lat_mat;
    }

    ///Access coord_trans[CART] in an intuitive way.
    const Eigen::Matrix3d &inv_lat_column_mat() const {
      return m_inv_lat_mat;
    }

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
    Lattice make_supercell(const Eigen::Matrix<T, 3, 3> &trans_mat) const;

    ///Find the lattice vectors which give most compact unit cell
    Lattice get_reduced_cell() const;  //Alex did this

    /// populate voronoi information.
    //It should get called in any method of Lattice that tries to use an empty voronoi table
    void generate_voronoi_table() const;
    void print_voronoi_table(std::ostream &stream) const;

    /// Radius of largest sphere that totally fits inside the voronoi cell
    double min_voronoi_radius() const;

    /// returns voronoi vector that maximizes dot product with 'pos'
    Eigen::Vector3d max_voronoi_vector(const Eigen::Vector3d &pos) const;

    /// return number of voronoi cell faces that 'pos' is on, within +/- TOL
    /// 0 indicates that 'pos' is within the voronoi cell, -1 indicates that it is outside
    int voronoi_number(const Eigen::Vector3d &pos) const;

    /**Gives a scaling of the lattice vectors such that after scaling,
       The eight parallelipipeds touching the origin enclose a sphere
       of given radius**/
    Eigen::Vector3i enclose_sphere(double radius) const;

    template<typename CoordType, typename CoordType2>
    Array<CoordType> gridstruc_build(double max_radius, double min_radius, Array<CoordType> basis, CoordType2 lat_point); //Anirudh

    void read(std::istream &stream);
    void print(std::ostream &stream) const;

    void symmetrize(double _tol = TOL);

    //template <class T>
    //Lattice &operator*=(const Matrix3<T> &RHS);

    //Lattice &operator*=(const Eigen::Matrix3d &RHS);

    ///Are two lattices the same, even if they have different lattice vectors
    bool is_equivalent(const Lattice &RHS, double tol) const;

    ///Are lattice vectors identical for two lattices
    bool operator==(const Lattice &RHS) const;

    ///Matrix that relates two lattices (e.g., strain or slat)
    //Eigen::Matrix3d operator/(const Lattice &RHS);

    //John G 121212
    ///Checks if lattice is a supercell of tile, acting on multiplication matrix. Check is performed applying operations from symlist
    bool is_supercell_of(const Lattice &tile, Eigen::Matrix3d &multimat, double _tol = TOL) const;
    bool is_supercell_of(const Lattice &tile, const Array<SymOp> &symlist, Eigen::Matrix3d &multimat, double _tol = TOL) const;

    ///Checks if lattice is a supercell of tile, applying operations from symlist
    bool is_supercell_of(const Lattice &tile, double _tol = TOL) const;
    bool is_supercell_of(const Lattice &tile, const Array<SymOp> &symlist, double _tol = TOL) const;

    ///Return a lattice with diagonal matrix that fits around starting lattice
    Lattice box(const Lattice &prim, const Lattice &scel, bool verbose = false) const;

    ///Flip c vector if it's on the wrong side of a-b plane -- return (*this)
    Lattice &make_right_handed();
    ///Check if the lattice is right handed
    bool is_right_handed() const;
    ///Given a normal vector, a Vector3 containing the miller indeces for the lattice is generated
    Eigen::Vector3i get_millers(Eigen::Vector3d plane_normal, double tolerance = TOL) const;
    ///Generates a lattice with vectors a and b parallel to the plane described by the miller indeces

    Lattice get_lattice_in_plane(Eigen::Vector3i millers, int max_vol = 20) const; //John G 121030

    Array<double> pg_converge(double large_tol);
    void pg_converge(double small_tol, double large_tol, double increment);

    ///Returns mirror image of lattice
    //Lattice get_reflection(bool override = 0) const;
    //\John G 121212

    //John G 011013
    void symmetrize(const SymGroup &relaxed_points);
    void symmetrize(const double &tolerance);

  private:

    mutable std::vector<Eigen::Vector3d> voronoi_table;

    //Coordinate Conversion Matrices
    //0 is fractional to cartesion; 1 is cartesian to fractional
    //use FRAC and CART globals to index
    //e.g., coord[CART]=coord_trans[FRAC]*coord[FRAC]; //converts frac to cart
    //Word to the wise: coord_trans[FRAC] is the matrix with columns equal to the lattice vectors
    Eigen::Matrix3d m_lat_mat, m_inv_lat_mat;

    //int periodicity_dim;   //dimension of periodicity
    //int periodicity_axis;  //index of lattice vector that is non-periodic (2d) or periodic (1d)

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

  /// \brief Apply SymOp to a Lattice
  Lattice &apply(const SymOp &op, Lattice &lat);

  /// \brief Copy and apply SymOp to a Lattice
  Lattice copy_apply(const SymOp &op, const Lattice &lat);

  /// \brief Returns a super Lattice
  Lattice make_supercell(const Lattice &lat, const Eigen::Matrix3i &transf_mat);

  /// Check if scel is a supercell of unitcell unit and some integer transformation matrix T
  std::pair<bool, Eigen::MatrixXi> is_supercell(const Lattice &scel, const Lattice &unit, double tol);

  /// Check if there is a symmetry operation, op, and transformation matrix T,
  ///   such that scel is a supercell of the result of applying op to unit
  template<typename Object, typename OpIterator>
  std::pair<OpIterator, Eigen::MatrixXi> is_supercell(
    const Object &scel,
    const Object &unit,
    OpIterator begin,
    OpIterator end,
    double tol);

  std::istream &operator>>(std::istream &in, const Lattice &lattice_in);

  ///\brief returns Lattice that is smallest possible supercell of both input Lattice
  Lattice superdupercell(const Lattice &lat1, const Lattice &lat2);

  ///\brief returns Lattice that is smallest possible supercell of all input Lattice
  template<typename LatIterator, typename SymOpIterator>
  Lattice superdupercell(LatIterator begin,
                         LatIterator end,
                         SymOpIterator op_begin = SymOpIterator(),
                         SymOpIterator op_end = SymOpIterator());


  //********************************************************************
  /**
   * This function generates a grid of points between max_radius and
   * min_radius. Additionally, it also fills up the points with a basis
   */
  //********************************************************************

  template<typename CoordType, typename CoordType2>
  Array<CoordType> Lattice::gridstruc_build(double max_radius, double min_radius, Array<CoordType> basis, CoordType2 lat_point) {
    Eigen::Vector3i dim;
    dim = enclose_sphere(max_radius);
    EigenCounter<Eigen::Vector3i > grid_count(-dim, dim, Eigen::Vector3i(1));
    double min_dist, dist;
    Array<CoordType> gridstruc;
    Eigen::Vector3i temp;

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


  ///\brief returns Lattice that is smallest possible supercell of all input Lattice
  ///
  /// If SymOpIterator are provided they are applied to each Lattice in an attempt
  /// to find the smallest possible superdupercell of all symmetrically transformed Lattice
  template<typename LatIterator, typename SymOpIterator>
  Lattice superdupercell(LatIterator begin,
                         LatIterator end,
                         SymOpIterator op_begin,
                         SymOpIterator op_end) {

    Lattice best = *begin;
    for(auto it = ++begin; it != end; ++it) {
      Lattice tmp_best = superdupercell(best, *it);
      for(auto op_it = op_begin; op_it != op_end; ++op_it) {
        Lattice test = superdupercell(best, copy_apply(*op_it, *it));
        if(std::abs(volume(test)) < std::abs(volume(tmp_best))) {
          tmp_best = test;
        }
      }
      best = tmp_best;
    }
    return best;
  }

  //********************************************************************
  // A column of trans_mat specifies a lattice vector of the supercell in terms of the
  // lattice vectors of (*this) lattice.
  template <typename T>
  Lattice Lattice::make_supercell(const Eigen::Matrix<T, 3, 3> &trans_mat) const {
    return Lattice(lat_column_mat() * trans_mat);
  }

  /// \brief Returns the volume of a Lattice
  ///
  /// \returns volume of the Lattice
  ///
  /// \param lat a \ref Lattice
  ///
  inline double volume(const Lattice &lat) {
    return lat.lat_column_mat().determinant();
  }

  /// \brief Returns a super Lattice
  inline Lattice make_supercell(const Lattice &lat, const Eigen::Matrix3i &transf_mat) {
    return Lattice(Eigen::Matrix3d(lat.lat_column_mat()) * transf_mat.cast<double>());
  }

  /// Check if there is a symmetry operation, op, and transformation matrix T,
  ///   such that scel is a supercell of the result of applying op to unit
  ///
  /// \returns pair corresponding to first successful op and T, or with op=end if not successful
  template<typename Object, typename OpIterator>
  std::pair<OpIterator, Eigen::MatrixXi> is_supercell(
    const Object &scel,
    const Object &unit,
    OpIterator begin,
    OpIterator end,
    double tol) {

    std::pair<bool, Eigen::MatrixXi> res;
    for(auto it = begin; it != end; ++it) {
      res = is_supercell(scel, copy_apply(*it, unit), tol);
      if(res.first) {
        return std::make_pair(it, res.second);
      }
    }
    return std::make_pair(end, res.second);
  }

}
#endif
