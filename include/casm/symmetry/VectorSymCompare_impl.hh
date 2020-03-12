#ifndef CASM_VectorSymCompare_impl
#define CASM_VectorSymCompare_impl

#include "casm/symmetry/VectorSymCompare.hh"
#include "casm/symmetry/SymMatrixXd.hh"
#include "casm/misc/CASM_Eigen_math.hh"


namespace Eigen {
  /// \brief Prepare an element for comparison
  ///
  /// - Attempts to find a sparse set of spanning vectors, and sort them so that subspace matrix is nearly upper triangular (if possible)
  /// - Also enures that first nonzero element of each row (if there is one) is positive
  template<typename Derived>
  representation_prepare_impl(Eigen::MatrixBase<Derived> const &obj, double _tol) const {
    //std::cout << "reduced_column_echelon of \n" << obj << "\n is\n" <<reduced_column_echelon(obj, this->tol()) << "\n";
    Eigen::MatrixBase<Derived> result = Eigen::MatrixBase<Derived>(reduced_column_echelon(obj, _tol).householderQr().householderQ()).leftCols(obj.cols());
    Index col = 0;
    for(Index row = 0; row < result.rows(); ++row) {
      Index i = 0;
      for(i = col; i < result.cols(); ++i) {
        if(!almost_zero(result(row, i), _tol)) {
          if(result(row, i) < 0)
            result.col(i) *= -1;
          ++col;
          break;
        }
      }
    }
    return result;
  }
}

namespace CASM {

  /// \brief Return tolerance
  template<typename Base>
  double EigenSymCompare<Base>::tol() const {
    return m_tol;
  }


  /// \brief Orders 'prepared' elements in the same orbit
  ///
  /// - Returns 'true' to indicate A < B
  /// - Equivalence is indicated by \code !compare(A,B) && !compare(B,A) \endcode
  /// - Assumes elements are 'prepared' before being compared
  /// Implementation:
  /// - First compares by number of sites in cluster
  /// - Then compare all displacements, from longest to shortest
  template<typename Base>
  bool EigenSymCompare<Base>::invariants_compare_impl(const Element &A, const Element &B) const {
    if(A.cols() == B.cols()) {
      return CASM::compare(A.norm(), B.norm(), tol());
    }
    return A.cols() < B.cols();
  }

  /// \brief Compares 'prepared' elements
  ///
  /// - Returns 'true' to indicate B < A
  /// - Equivalence is indicated by \code !compare(A,B) && !compare(B,A) \endcode
  /// - Assumes elements are 'prepared' before being compared
  template<typename Base>
  bool EigenSymCompare<Base>::compare_impl(const Element &A, const Element &B) const {
    //std::cout << "Comparing \n   " << A.transpose() << "\n   " << B.transpose() << "\n\n";
    bool result = colmajor_lex_compare(B, A, tol());
    //std::cout << "RESULT: " << result << "\n\n";
    return result;
  }


  /// \brief Returns transformation that takes 'obj' to its prepared (canonical) form
  ///
  /// - For now returns pointer to SymPermutation object that encodes permutation due to sorting elements
  template<typename Base>
  std::unique_ptr<SymOpRepresentation> EigenSymCompare<Base>::canonical_transform_impl(Element const &obj)const {
    return std::unique_ptr<SymOpRepresentation>(new SymMatrixXd(obj.colPivHouseholderQr().solve(derived().representation_prepare_impl(obj))));
  }


  /// \brief Prepare an element for comparison
  ///
  /// - Attempts to find a sparse set of spanning vectors, and sort them so that subspace matrix is nearly upper triangular (if possible)
  /// - Also enures that first nonzero element of each row (if there is one) is positive
  template<typename Element, typename SymApply>
  Element SubspaceSymCompare<Element, SymApply>::
  representation_prepare_impl(Element obj) const {
    return representation_prepare_impl(obj, this->tol());
  }


}

#endif
