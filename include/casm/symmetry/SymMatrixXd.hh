#ifndef SYMMATRIXXD_HH
#define SYMMATRIXXD_HH

#include <iostream>
#include <cmath>

#include "casm/external/Eigen/Dense"

#include "casm/symmetry/SymOpRepresentation.hh"

namespace CASM {
  namespace xtal {
    class Lattice;
  }
  using xtal::Lattice;

  class MasterSymGroup;

  /** \ingroup Symmetry
   *  @{
   */

  ///\brief  Generalized symmetry matrix representation for arbitrary dimension
  /// Can be used to describe application of symmetry to N-dimensional vector spaces
  /// Use for 3-dimensional transformations if they do not describe coordinate transformations
  class SymMatrixXd : public SymOpRepresentation {
  private:
    Eigen::MatrixXd mat;
  public:
    SymMatrixXd(const Eigen::MatrixXd &init_mat) : mat(init_mat) {};

    SymOpRepresentation *copy() const override {
      return new SymMatrixXd(*this);
    }

    Eigen::MatrixXd const *MatrixXd() const override {
      return &mat;
    }

    double character() const override {
      return mat.trace();
    }

  };

  /** @} */
}
#endif
