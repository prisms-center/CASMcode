#ifndef SYMMATRIXXD_HH
#define SYMMATRIXXD_HH

#include <iostream>
#include <cmath>

#include "casm/external/Eigen/Dense"

#include "casm/symmetry/SymOpRepresentation.hh"

namespace CASM {
  class MasterSymGroup;
  class Lattice;
  class jsonParser;

  ///\brief  Generalized symmetry matrix representation for arbitrary dimension
  /// Can be used to describe application of symmetry to N-dimensional vector spaces
  /// Use for 3-dimensional transformations if they do not describe coordinate transformations
  class SymMatrixXd : public SymOpRepresentation {
  private:
    Eigen::MatrixXd mat;
  public:
    /// fixes alignment of mat
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    SymMatrixXd(const Eigen::MatrixXd &init_mat) : mat(init_mat) {};
    SymOpRepresentation *copy() const {
      return new SymMatrixXd(*this);
    };
    Eigen::MatrixXd const *get_MatrixXd() const {
      return &mat;
    };
    double get_character() const {
      return mat.trace();
    };
    jsonParser &to_json(jsonParser &json) const;

    void from_json(const jsonParser &json);
  };

  jsonParser &to_json(const SymMatrixXd &sym, jsonParser &json);
  void from_json(SymMatrixXd &sym, const jsonParser &json);

}
#endif
