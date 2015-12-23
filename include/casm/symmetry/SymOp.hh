#ifndef SYMOP_HH
#define SYMOP_HH

#include <iostream>
#include <cmath>

#include "casm/crystallography/Coordinate.hh"
#include "casm/symmetry/SymOpRepresentation.hh"

namespace CASM {
  class MasterSymGroup;
  class Lattice;

  ///\brief SymOp is the Coordinate representation of a symmetry operation
  /// it keeps fraction (FRAC) and Cartesian (CART) information about how
  /// a symetry operation transforms 3D spatial coordinates
  class SymOp : public SymOpRepresentation {
  private:

    /*** Inherited from SymOpRepresentation ***/
    /*
    ///Index into MasterSymGroup that specifies the operation
    Index op_index, rep_ID;

    /// Pointer to the MasterSymGroup where prototype of this SymOp lives
    MasterSymGroup const *master_group;
    */

    ///Lattice in terms of  which the fractional
    ///coordinates are defined
    Lattice const *home;

    ///matrix representation of symettry operation in
    ///Cartesian (symmetry_mat[CART])
    ///Fractional (symmetry_mat[FRAC])
    mutable Matrix3<double> symmetry_mat[2];
    ///Keep track of flags to specify whether
    ///Cartesian/Fractional symmetry matrix is up to date
    mutable bool is_current[2];

    ///translation vector that is applied to
    ///a point after matrix transformation
    mutable Coordinate tau_vec;

    ///eigenvec describes the orientation of the symmetry-generating element
    ///(i.e., mirror/glide plane normal or vector parallel to rotation/screw axis)
    mutable Coordinate eigenvec;

    ///location is a point in the lattice that describes the location
    ///of the symmetry-generating element (e.g, mirror plane, rotation axis, inversion point)
    ///location is related to tau_vec in that
    /// new_coord=symmetry_mat*(old_coord-location) + location
    mutable Coordinate location;

    ///For rotation/rotoinversion, 'rotation_angle' is the corresponding angle (in degrees)
    mutable double rotation_angle;

    ///shift vector for screw or glide, perpendicular to eigenvector
    mutable Coordinate screw_glide_shift;

    /// This stores the mapping error associated with this SymOp, which will depend on the tolerances
    ///  you choose when you attempt to generate symmetry groups.
    double map_error;

    int mode_ind() const {
      return (int) COORD_MODE::CHECK();
    };
  public:

    SymOp(const Coordinate &trans) :
      home(trans.get_home()),
      tau_vec(trans),
      eigenvec(Vector3<double>(0, 0, 0), *(trans.get_home()), CART),
      location(eigenvec),
      screw_glide_shift(eigenvec),
      map_error(TOL) {
      symmetry = identity_op;
      op_index = 0;
      symmetry_mat[CART] = Matrix3<double>::identity();
      symmetry_mat[FRAC] = Matrix3<double>::identity();
      is_current[CART] = true;
      is_current[FRAC] = true;
    };

    ///Create new SymOp from Matrix3 in specified mode, but must also specify the lattice
    ///Assumes translation is zero
    SymOp(const Matrix3<double> &tsym_mat, const Lattice &t_home,
          COORD_TYPE mode, double _map_error = TOL) :
      home(&t_home),
      tau_vec(Vector3<double>(0, 0, 0), t_home, mode),
      eigenvec(tau_vec),
      location(tau_vec),
      screw_glide_shift(tau_vec),
      map_error(_map_error) {
      symmetry_mat[mode] = tsym_mat;
      is_current[mode] = true;
      is_current[!mode] = false;
      if(mode == FRAC) {
        calc(CART);
      }
      else {
        calc(FRAC);
      }

    };

    ///Create new SymOp from Matrix3 and tau translation in specified mode, but must also specify the lattice
    SymOp(const Matrix3<double> &tsym_mat, const Vector3<double> &t_tau, const Lattice &t_home,
          COORD_TYPE mode, double _map_error = TOL) :
      home(&t_home),
      tau_vec(t_tau, t_home, mode),
      eigenvec(Vector3<double>(0, 0, 0), t_home, mode),
      location(eigenvec),
      screw_glide_shift(eigenvec),
      map_error(_map_error) {
      symmetry_mat[mode] = tsym_mat;
      is_current[mode] = true;
      is_current[!mode] = false;
      if(mode == FRAC)
        calc(CART);
      else
        calc(FRAC);
    };


    /// Const access of elements of the symmetry matrix
    double operator()(int i, int j, COORD_TYPE mode = COORD_DEFAULT) const;

    /// Const access of the translation vector, 'tau', in specified COORD_MODE
    const Vector3<double> &tau(COORD_TYPE mode) const; //AAB
    /// Const access of the translation vector, 'tau', in current COORD_MODE
    const Coordinate &tau() const; //AAB

    /// Const access of elements of the translation vector, 'tau', in specified COORD_MODE
    double tau(int i, COORD_TYPE mode = COORD_DEFAULT) const ;

    /// Const access of entire symmetry matrix in specified COORD_MODE
    const Matrix3<double> &get_matrix(COORD_TYPE mode) const; //AAB
    /// Const access of entire symmetry matrix in the default COORD_MODE
    const Matrix3<double> &get_matrix() const; //AAB

    /// Const access of eigenvector in the specified COORD_MODE
    const Vector3<double> &get_eigenvec(COORD_TYPE mode) const; //AAB
    /// Const access of eigenvector in the current COORD_MODE
    const Vector3<double> &get_eigenvec() const; //AAB

    /// Const access of screw/glide shift vector in the specified COORD_MODE
    const Vector3<double> &get_screw_glide_shift(COORD_TYPE mode) const; //Donghee
    /// Const access of screw/glide shift vector in the current COORD_MODE
    const Coordinate &get_screw_glide_shift() const; //Donghee

    /// Const access of location Coordinate in the specified COORD_MODE
    const Vector3<double> &get_location(COORD_TYPE mode) const;
    /// Const access of location Coordinate in the current COORD_MODE
    const Coordinate &get_location() const;

    /// set master_group and op_index for the SymOp
    /// Performs simple checks to ensure that (*this) is compatible with new_group[new_index]
    void set_index(const MasterSymGroup &new_group, Index new_index);

    /// set representation for SymOp corresponding to rep_ID
    void set_rep(Index rep_ID, const SymOpRepresentation &op_rep) const;

    /// get pointer to matrix representation corresponding to rep_ID
    Eigen::MatrixXd const *get_matrix_rep(Index rep_ID) const;

    /// get pointer to permutation representation corresponding to rep_ID
    Permutation const *get_permutation_rep(Index rep_ID) const;

    /// get pointer to BasisPermute representation corresponding to rep_ID
    SymBasisPermute const* get_basis_permute_rep(Index rep_ID) const;

    /// get array of pointers to matrix representations for representations corresponding to rep_IDs
    Array<Eigen::MatrixXd const * > get_matrix_reps(Array<Index> rep_IDs) const;

    /// Allows access to the map_error.
    const double &get_map_error() const;
    void set_map_error(const double &value);

    //Routines that relate symmetry options:
    ///get a new SymOp that is equivalent to subsequent application of both SymOps
    SymOp operator*(const SymOp &RHS) const;

    /// Cartesian translation of the origin of the symmetry operation by Coordinate RHS
    SymOp &operator+=(const Coordinate &RHS);
    SymOp &operator-=(const Coordinate &RHS);

    /// get the inverse of this SymOp
    SymOp inverse() const;

    /// Get copy of the SymOp without translation
    SymOp no_trans() const;

    /// Brings the shift vector within...
    void within();

    /// Check equality of SymOps, (matrix and translation). Does not necessarily return true for translational equivalents
    bool operator==(const SymOp &RHS) const;
    /// Check equality of SymOps, with specified tolerance, (matrix and translation). Returns true for translational equivalents
    bool compare(const SymOp &RHS, double eq_tol = TOL) const;

    /// Const access of the symmetry type
    symmetry_type type() const {
      if(symmetry == invalid_op)
        get_sym_type();
      return symmetry;
    }

    /// Calculate the matrix in the specified COORD_MODE
    bool calc(COORD_TYPE mode = COORD_DEFAULT) const;

    /// Change home lattice
    /// Preserves Cartesian representation and changes fractional
    /// WARNING: Never add a default value for 'mode'
    void set_lattice(const Lattice &new_lat, COORD_TYPE mode);

    /// Const access of the home lattice
    const Lattice &get_home() const;


    /// Performs conjugation of this SymOp with SymOp 'op'
    /// In other words, transforms coordinates of this SymOp by
    /// transformation specified by 'op':  op.inverse()*(*this)*op
    SymOp &apply_sym(const SymOp &op);

    /// Use 'symmetry_mat' and 'tau_vec' to populate 'location'
    void find_location() const;

    /// What does this do?
    mutable double plane_XYZ_equ[4], axis_XYZ_equ[6]; //Donghee
    void get_plane_axis(COORD_TYPE coord_mode) const;  //Donghee
    void print_plane_axis(std::ostream &stream, COORD_TYPE coord_mode) const; //Donghee
    void printII_plane_axis(std::ostream &stream, COORD_TYPE coord_mode) const;  //Donghee

    /// Use 'symmetry_mat' and 'tau_vec' to populate 'symmetry'
    void get_sym_type() const;

    /// Auxiliary routines to get_sym_type()
    void glide_check() const;
    void screw_check() const;
    void mirror_check() const;
    void calc_rotation_angle(Matrix3<double> mat, double det, double trace) const;

    /// const access of 'rotation_angle'
    double get_rotation_angle() const;

    /// calculate and return character of this SymOp
    double get_character() const {
      return get_matrix(CART).trace();
    };

    /// Print this SymOP as header, followed by symmetry_mat and tau_vec (aligned vertically)
    ///COORD_DEFAULT is default argument, this means default to the global mode
    void print(std::ostream &stream, COORD_TYPE coord_mode = COORD_DEFAULT) const;

    ///Prints abridged description of SymOp, including its type, angle, and eigenvector
    ///does NOT print out the entire symmetry operation matrix
    void print_short(std::ostream &stream, COORD_TYPE coord_mode = COORD_DEFAULT) const;

    /// Read in a SymOp in the same format as it is printed by 'print(...)'
    //  NOT YET IMPLEMENTED
    void read(std::istream &stream, COORD_TYPE coord_mode = COORD_DEFAULT) {};

    /// Makes sure that all properties of SymOp are up to date in the specified COORD_MODE
    /// THIS HAS NOT BEEN IMPLEMENTED --- we should discuss specifics
    void update(COORD_TYPE coord_mode);


    /// Methods for checking symmetry type.
    bool is_identity() const;
    bool is_mirror() const;
    bool is_glide() const;
    bool is_rotation() const;
    bool is_screw() const;
    bool is_inversion() const;
    bool is_rotoinversion() const;
    bool is_invalid() const;

    /// Return pointer to a new copy of this SymOp
    SymOpRepresentation *copy() const {
      return new SymOp(*this);
    };

    jsonParser &to_json(jsonParser &json) const;

    void from_json(const jsonParser &json);
  };

  jsonParser &to_json(const SymOp &sym, jsonParser &json);
  void from_json(SymOp &sym, const jsonParser &json);

  //std::istream& operator>> (std::istream &stream, SymOp &op){op.read(stream); return stream;};
  //std::ostream& operator<< (std::ostream &stream, SymOp &op){op.print(stream); return stream;};

}
#endif
