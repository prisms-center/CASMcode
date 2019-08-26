#ifndef LATTICE_HH
#define LATTICE_HH

#include <iostream>
#include <cmath>
#include <vector>

#include "casm/CASM_global_Eigen.hh"
#include "casm/misc/Comparisons.hh"
#include "casm/CASM_global_definitions.hh"
#include "casm/crystallography/LatticeCanonicalForm.hh"

namespace CASM {
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  class SymOp;
  class ScelEnumProps;
  class jsonParser;

  /** \defgroup Crystallography
   *
   *  \brief Relates to crystallography
   */

  /** \defgroup Lattice
   *  \ingroup Crystallography
   *  \brief Relates to Lattice
   *
   *  @{
   */

  class Lattice : public LatticeCanonicalForm<Comparisons<CRTPBase<Lattice>>> {
  public:
    typedef Eigen::Matrix3d::ColXpr LatVec;
    typedef Eigen::Matrix3d::ConstColXpr ConstLatVec;

    Lattice(Eigen::Ref<const Eigen::Vector3d> const &vec1,
            Eigen::Ref<const Eigen::Vector3d> const &vec2,
            Eigen::Ref<const Eigen::Vector3d> const &vec3,
            double xtal_tol = TOL,
            bool force = false);

    ///Construct Lattice from a matrix of lattice vectors, where lattice vectors are columns
    ///(e.g., lat_mat is equivalent to coord_trans[FRAC])
    Lattice(Eigen::Ref<const Eigen::Matrix3d> const &lat_mat = Eigen::Matrix3d::Identity(),
            double xtal_tol = TOL, bool force = false);

    /// \brief Construct FCC primitive cell of unit volume
    static Lattice fcc(double tol = TOL);

    /// \brief Construct BCC primitive cell of unit volume
    static Lattice bcc(double tol = TOL);

    /// \brief Construct simple cubic primitive cell of unit volume
    static Lattice cubic(double tol = TOL);

    /// \brief Construct cubic primitive cell of unit volume
    static Lattice hexagonal(double tol = TOL);

    static std::vector<Eigen::Matrix3d> const &skew_transforms();

    /// \brief Get i'th lattice vector as column expression
    LatVec operator[](Index i) {
      return m_lat_mat.col(i);
    }

    /// \brief Get i'th lattice vector as column expression
    ConstLatVec operator[](Index i)const {
      return m_lat_mat.col(i);
    }

    std::tuple<LatVec, LatVec, LatVec> vectors() {
      return std::make_tuple(m_lat_mat.col(0), m_lat_mat.col(1), m_lat_mat.col(2));
    }

    std::tuple<ConstLatVec, ConstLatVec, ConstLatVec> vectors() const {
      return std::make_tuple(m_lat_mat.col(0), m_lat_mat.col(1), m_lat_mat.col(2));
    }

    /// \brief Return scaled copy of this lattice (Note: Volume will be scaled by scale^3)
    Lattice scaled_lattice(double scale) const;

    /// \brief Return angle between lattice vectors (*this)[(i+1)%3] and (*this)[(i+2)%3], in degrees
    double angle(Index i)const;

    /// \brief Return length of i'th lattice vector
    double length(Index i)const;

    /// \brief Return *signed* volume of this lattice
    double vol()const {
      return lat_column_mat().determinant();
    }

    /// \brief Return voronoi table, which specifies outward-pointing normals of Lattice Voronoi cell.
    /// outward-pointing normals are given as rows of the matrix, and defined such that
    /// if (voronoi_table()*coord.cart()).maxCoeff()>1, then 'coord' is outside of the voronoi cell
    Eigen::MatrixXd const &voronoi_table() const {
      if(!m_voronoi_table.size()) {
        _generate_voronoi_table();
      }
      return m_voronoi_table;
    }

    /// \brief Radius of the largest sphere that fits wholly within the voronoi cell
    double inner_voronoi_radius() const {
      voronoi_table();
      return m_inner_voronoi_radius;
    }

    /// \brief 3x3 matrix with lattice vectors as its columne
    const Eigen::Matrix3d &lat_column_mat() const {
      return m_lat_mat;
    }

    /// \brief Inverse of Lattice::lat_column_mat()
    /// It is the transformation matrix 'C2F', such that
    ///    f = C2F * c
    /// where 'c' is Cartesian coordinate/vector and 'f'
    /// and 'f' is fractional coordinate/vector
    const Eigen::Matrix3d &inv_lat_column_mat() const {
      return m_inv_lat_mat;
    }

    /// Calculates the kpoint mesh for a supercell lattice given the kpoint mesh
    /// for the primitive lattice
    std::vector<int> calc_kpoints(std::vector<int> prim_kpoints, Lattice prim_lat);

    /// \brief Return reciprocal lattice
    Lattice reciprocal() const;

    /// \brief Return boxiness factor directly proportional to volume/SA ratio
    double boxiness() const;


    /// \brief Populate \param supercell with symmetrically distinct supercells of this lattice
    /// Superlattices are enumerated with volumes \param min_prim_vol <= volume <= \param max_prim_vol
    /// \param effective_pg is a group that should either be equivalent to the full point group of this lattice
    /// or be a subgroup of that full point group
    void generate_supercells(std::vector<Lattice> &supercell, const std::vector<SymOp> &effective_pg, const ScelEnumProps &enum_props) const;

    /// \brief make a supercell of this lattice.
    /// Equivalent to Lattice(lat_column_mat()*trans_mat)
    template <typename T>
    Lattice make_supercell(const Eigen::Matrix<T, 3, 3> &trans_mat) const;

    /// \brief Find the lattice vectors which give most compact unit cell
    /// Compactness is measured by how close lat_column_mat().transpose()*lat_column_mat() is to a diagonal matrix
    Lattice reduced_cell() const;

    /// \brief Find the lattice vectors which give most compact unit cell
    /// Compactness is measured by how close lat_column_mat().transpose()*lat_column_mat() is to a diagonal matrix
    Lattice reduced_cell2() const;

    void print_voronoi_table(std::ostream &stream) const;

    /// Radius of largest sphere that totally fits inside the voronoi cell
    double min_voronoi_radius() const;

    /// Given cartesian vector 'pos', returns double v_dist and populates 'lattice_trans'.
    /// v_dist is the fractional number of half lattice-planes along 'lattice_trans' between 'pos' and the lattice plane that passes through the origin.
    /// lattice_trans is a lattice translation that is normal to a face of the voronoi cell and translating pos to 'pos-lattice_trans'
    /// will yield a translationally-equivalent vector between the -1 and +1 half-plane.
    /// v_dist is maximized over all faces of the voronoi cell.
    double max_voronoi_measure(const Eigen::Vector3d &pos, Eigen::Vector3d &lattice_trans) const;

    /// return number of voronoi cell faces that 'pos' is on, within +/- TOL
    ///  0 indicates that 'pos' is within the voronoi cell, and the origin is the nearest lattice site
    /// -1 indicates that 'pos' is outside the voronoi cell, and there is a lattice site closer than the origin
    /// Values of 1<=n<=7 indicate that there n lattice sites equally as close as the origin
    int voronoi_number(const Eigen::Vector3d &pos) const;

    /// Gives a scaling of the lattice vectors such that after scaling,
    /// The eight parallelipipeds touching the origin enclose a sphere of given radius
    Eigen::Vector3i enclose_sphere(double radius) const;

    /// Make a grid of lattice sites such that min_radius <= distance <= max_radius from \param lat_point
    // ***This should live somewhere else
    template<typename CoordType, typename CoordType2>
    std::vector<CoordType> gridstruc_build(double max_radius, double min_radius, std::vector<CoordType> basis, CoordType2 lat_point);

    void read(std::istream &stream);
    void print(std::ostream &stream, int _prec = 8) const;

    /// \brief Compare two Lattice
    bool operator<(const Lattice &RHS) const;

    ///Matrix that relates two lattices (e.g., strain or slat)
    //Eigen::Matrix3d operator/(const Lattice &RHS);

    //John G 121212
    ///Checks if lattice is a supercell of tile, acting on multiplication matrix. Check is performed applying operations from symlist
    bool is_supercell_of(const Lattice &tile, Eigen::Matrix3d &multimat) const;
    bool is_supercell_of(const Lattice &tile, const std::vector<SymOp> &symlist, Eigen::Matrix3d &multimat) const;

    ///Checks if lattice is a supercell of tile, applying operations from symlist
    bool is_supercell_of(const Lattice &tile) const;
    bool is_supercell_of(const Lattice &tile, const std::vector<SymOp> &symlist) const;

    ///Return a lattice with diagonal matrix that fits around starting lattice
    Lattice box(const Lattice &prim, const Lattice &scel, bool verbose = false) const;

    ///Flip c vector if it's on the wrong side of a-b plane -- return (*this)
    Lattice &make_right_handed();

    ///Check if the lattice is right handed
    bool is_right_handed() const;

    ///Given a normal vector, a Vector3 containing the miller indeces for the lattice is generated
    Eigen::Vector3i millers(Eigen::Vector3d plane_normal) const;

    ///Generates a lattice with vectors a and b parallel to the plane described by the miller indeces
    Lattice lattice_in_plane(Eigen::Vector3i millers, int max_vol = 20) const; //John G 121030

    std::vector<double> pg_converge(double large_tol);
    void pg_converge(double small_tol, double large_tol, double increment);

    /// \brief Symmetrized copy of this lattice, having symmetry of group \param _pg
    Lattice symmetrized(const std::vector<SymOp> &_pg) const;

    /// \brief Symmetrized copy of this lattice, having symmetry of group formed by fractional-coordinate symops \param _pg
    Lattice symmetrized(const std::vector<Eigen::Matrix3i> &_pg) const;

    /// \brief Symmetrized copy of this lattice having symmetry of point group calculated based on tolerance \param _tol
    Lattice symmetrized(double _tol) const;

    double tol() const {
      return m_tol;
    }

    void set_tol(double _tol) {
      m_tol = _tol;
    }

  private:

    friend Comparisons<Lattice>;

    /// Are lattice vectors identical for two lattices, within TOL
    bool _eq(const Lattice &RHS) const;

    /// \brief populate voronoi information.
    void _generate_voronoi_table() const;

    mutable double m_inner_voronoi_radius;
    mutable Eigen::MatrixXd m_voronoi_table;

    //Coordinate Conversion Matrices
    //0 is fractional to cartesion; 1 is cartesian to fractional
    //use FRAC and CART globals to index
    //e.g., coord[CART]=coord_trans[FRAC]*coord[FRAC]; //converts frac to cart
    //Word to the wise: coord_trans[FRAC] is the matrix with columns equal to the lattice vectors
    Eigen::Matrix3d m_lat_mat, m_inv_lat_mat;

    double m_tol;
  };


  // write Lattice in json as array of vectors
  jsonParser &to_json(const Lattice &lat, jsonParser &json);
  void from_json(Lattice &lat, const jsonParser &json, double xtal_tol);


  /* never write a Matrix*Lattice operator, PLEASE

     template <class T>
     Lattice operator*(const Matrix3<T> &LHS, const Lattice &RHS);

     Lattice operator*(const Eigen::Matrix3d &LHS, const Lattice &RHS);
  */

  /// \brief Populate \param point_group with the point group of this lattice
  /// \param point_group should be empty
  /// \param pg_tol can be increased to find point group of lattice vectors
  /// that are slightly distorted due to numerical noise
  std::vector<SymOp> calc_point_group(Lattice const &_lat);

  /// \brief Populate \param point_group with the point group of this lattice
  /// \param point_group should be empty
  /// \param pg_tol can be increased to find point group of lattice vectors
  /// that are slightly distorted due to numerical noise
  std::vector<SymOp> calc_point_group(Lattice const &_lat, double _tol);


  /// \brief Returns the volume of a Lattice
  double volume(const Lattice &lat);

  /// \brief Apply SymOp to a Lattice
  Lattice &apply(const SymOp &op, Lattice &lat);

  /// \brief Copy and apply SymOp to a Lattice
  Lattice copy_apply(const SymOp &op, const Lattice &lat);

  /// \brief Returns a super Lattice
  Lattice make_supercell(const Lattice &lat, const Eigen::Matrix3i &transf_mat);

  /// Check if scel is a supercell of unitcell unit and some integer transformation matrix T
  std::pair<bool, Eigen::Matrix3d> is_supercell(const Lattice &scel, const Lattice &unit, double tol);

  /// Check if there is a symmetry operation, op, and transformation matrix T,
  ///   such that scel is a supercell of the result of applying op to unit
  template<typename Object, typename OpIterator>
  std::pair<OpIterator, Eigen::Matrix3d> is_supercell(
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

  Lattice replace_vector(const Lattice &lat, const Eigen::Vector3d &new_vector, double tol);




  //********************************************************************
  ///\brief Returns 'frac_mat' which is transformation of 'cart_mat'
  /// if
  ///    cart_vec_after = cart_mat*cart_vec
  /// then
  ///    frac_vec_after = frac_mat*frac_vec
  /// where cart_vec = lat.lat_column_mat()*frac_vec
  /// and cart_vec_after = lat.lat_column_mat()*frac_vec_after
  inline
  Eigen::Matrix3d cart2frac(const Eigen::Ref<const Eigen::Matrix3d> &cart_mat, const Lattice &lat) {
    return lat.inv_lat_column_mat() * cart_mat * lat.lat_column_mat();
  }

  //********************************************************************
  ///\brief Returns 'cart_mat' which is transformation of 'frac_mat'
  /// if
  ///    cart_vec_after = cart_mat*cart_vec
  /// then
  ///    frac_vec_after = frac_mat*frac_vec
  /// where cart_vec = lat.lat_column_mat()*frac_vec
  /// and cart_vec_after = lat.lat_column_mat()*frac_vec_after
  inline
  Eigen::Matrix3d frac2cart(const Eigen::Ref<const Eigen::Matrix3d> &frac_mat, const Lattice &lat) {
    return lat.lat_column_mat() * frac_mat * lat.inv_lat_column_mat();
  }

  //********************************************************************
  // A column of trans_mat specifies a lattice vector of the supercell in terms of the
  // lattice vectors of (*this) lattice.
  template <typename T>
  Lattice Lattice::make_supercell(const Eigen::Matrix<T, 3, 3> &trans_mat) const {
    return Lattice(lat_column_mat() * trans_mat, tol());
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
    return Lattice(Eigen::Matrix3d(lat.lat_column_mat()) * transf_mat.cast<double>(), lat.tol());
  }

  /** @} */

}
#endif
