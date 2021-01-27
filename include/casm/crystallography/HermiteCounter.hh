#ifndef HERMITECOUNTER_HH
#define HERMITECOUNTER_HH

#include "casm/container/Counter.hh"
#include "casm/external/Eigen/Dense"
#include "casm/global/definitions.hh"

namespace CASM {
namespace xtal {
/**
 * Given the dimensions of a square matrix and its determinant,
 * HermiteCounter will cycle through every possible matrix that
 * maintains it's Hermite normal form:
 *  -Upper triangular matrix
 *  -Determinant remains constant
 *  -row values to the right of the diagonal will always be smaller than
 *  the value of the diagonal
 *
 * In addition, this class is limited to SQUARE matrices, and NONZERO
 * determinants. The idea is to use it for supcercell enumerations,
 * where these conditions are always met.
 *
 * For a determinant det, the initial value of the counter will be
 * a n x n identity matrix H with H(0,0)=det.
 * The final position will be a n x n identity matrix with H(n-1,n-1)=det.
 * Once the final position for a particular determinant is reached,
 * the counter starts over with the next integer determinant value.
 *
 * There are two main steps in the counter:
 *  -Incrementing the diagonal of the matrix such that its product remains
 *  equal to the determinant
 *  -Incrementing the upper triangular values such that they never exceed
 *  the diagonal
 *
 * The diagonal increments are achieved by working only with two adjacent
 * values at a time, distributing factors to the next diagonal element
 * only if it equals 1. If not, the next adjacent pair is selected.
 */

class HermiteCounter {
 public:
  typedef Eigen::VectorXi::Scalar value_type;
  typedef CASM::Index Index;

  /// \brief constructor to satisfy iterator requirements. Do not recommend.
  // HermiteCounter() {};

  /// \brief constructor given the desired determinant and square matrix
  /// dimensions
  HermiteCounter(int init_determinant, int init_dim);

  // You probably will never need these. They're just here for testing more than
  // anything. Either way, they're safe to call.
  Index position() const;
  Eigen::VectorXi diagonal() const;
  // value_type low() const;
  // value_type high() const;

  /// \brief Get the current matrix the counter is on
  Eigen::MatrixXi current() const;

  /// \brief Get the current determinant
  value_type determinant() const;

  /// \brief Get the dimensions of *this
  Index dim() const;

  /// \brief reset the counter to the first iteration of the current determinant
  void reset_current();

  /// \brief Skip the remaining iterations and start at the next determinant
  /// value
  void next_determinant();

  /// \brief Reset the diagonal to the specified determinant and set the other
  /// elements to zero
  void jump_to_determinant(value_type new_det);

  /// \brief Jump to the next available HNF matrix.
  HermiteCounter &operator++();

  /// \brief Get the current matrix the counter is on
  Eigen::MatrixXi operator()() const;

 private:
  /// \brief Keeps track of the current diagonal element that needs to be
  /// factored
  Index m_pos;

  /// \brief Vector holding diagonal element values
  Eigen::VectorXi m_diagonal;

  /// \brief unrolled vector of the upper triangle (does not include diagonal
  /// elements)
  EigenVectorXiCounter m_upper_tri;

  /// \brief Go to the next values of diagonal elements that keep the same
  /// determinant
  Index _increment_diagonal();
};

namespace HermiteCounter_impl {
/// \brief Find the next factor of the specified position and share with next
/// element. Use attempt as starting point.
HermiteCounter::Index _spill_factor(Eigen::VectorXi &diag,
                                    HermiteCounter::Index position,
                                    HermiteCounter::value_type attempt);

/// \brief Spill the next factor of the specified element with its neighbor, and
/// return new position
HermiteCounter::Index next_spill_position(Eigen::VectorXi &diag,
                                          HermiteCounter::Index position);

/// \brief Determine the number of elements in the upper triangular matrix
/// (excluding diagonal)
HermiteCounter::Index upper_size(HermiteCounter::Index init_dim);

/// \brief Create a counter for the elements above the diagonal based on the
/// current diagonal value
EigenVectorXiCounter _upper_tri_counter(const Eigen::VectorXi &current_diag);

/// \brief Assemble a matrix diagonal and unrolled upper triangle values into a
/// matrix
Eigen::MatrixXi _zip_matrix(const Eigen::VectorXi &current_diag,
                            const Eigen::VectorXi &current_upper_tri);

/// \brief Expand a n x n Hermite normal matrix into a m x m one (e.g. for 2D
/// supercells)
Eigen::MatrixXi _expand_dims_old(const Eigen::MatrixXi &hermit_mat,
                                 const Eigen::VectorXi &active_dims);

/// \brief Expand a n x n Hermite normal matrix (H) into a m x m one through a m
/// x m generating matrix (G) (e.g. for arbitrary 2D supercells)
Eigen::MatrixXi _expand_dims(const Eigen::MatrixXi &H,
                             const Eigen::MatrixXi &G);

/// \brief Unroll a Hermit normal form square matrix into a vector such that
/// it's canonical form is easy to compare
Eigen::VectorXi _canonical_unroll(const Eigen::MatrixXi &hermit_mat);

/// \brief Compare two integer matrices and see which one is lexicographically
/// greatest. Returns true if H0<H1
bool _canonical_compare(const Eigen::MatrixXi &H0, const Eigen::MatrixXi &H1);
}  // namespace HermiteCounter_impl

}  // namespace xtal
}  // namespace CASM

#endif
