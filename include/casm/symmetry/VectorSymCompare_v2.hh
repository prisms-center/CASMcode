#ifndef CASM_symmetry_VectorSymCompare_v2
#define CASM_symmetry_VectorSymCompare_v2

#include <set>

#include "casm/container/multivector.hh"
#include "casm/external/Eigen/Core"
#include "casm/global/definitions.hh"

namespace CASM {

namespace SymRepTools_v2 {

typedef std::vector<Eigen::MatrixXd> MatrixRep;
typedef std::set<Index> GroupIndices;
typedef std::set<GroupIndices> GroupIndicesOrbit;
typedef std::vector<GroupIndicesOrbit> GroupIndicesOrbitVector;

/// Used for generating SimpleOrbit of high symmetry direction vectors
class VectorInvariants {
 public:
  VectorInvariants(Eigen::VectorXcd const &vector);

  double cols() const;
  double norm() const;

 private:
  double m_cols;
  double m_norm;
};

}  // namespace SymRepTools_v2

bool almost_equal(SymRepTools_v2::VectorInvariants const &A_invariants,
                  SymRepTools_v2::VectorInvariants const &B_invariants,
                  double tol);

bool compare(SymRepTools_v2::VectorInvariants const &A_invariants,
             SymRepTools_v2::VectorInvariants const &B_invariants, double tol);

namespace SymRepTools_v2 {

/// Used for constructing SimpleOrbit of high symmetry direction vectors
struct VectorSymCompare {
  VectorSymCompare(MatrixRep const &matrix_rep, double tol);

  typedef Index SymOpRepType;
  typedef Eigen::VectorXcd Element;
  typedef VectorInvariants InvariantsType;

  // lexicographical comparison (reversed, uses vector_B < vector_A)
  bool compare(Eigen::VectorXcd const &vector_A,
               Eigen::VectorXcd const &vector_B) const;

  // return VectorInvariants
  VectorInvariants make_invariants(Eigen::VectorXcd const &vector) const;

  // apply matrix rep to vector
  Eigen::VectorXcd copy_apply(Index const &op_index,
                              Eigen::VectorXcd vector) const;

  // no change needed to prepare for comparison
  Eigen::VectorXcd prepare(Eigen::VectorXcd vector) const;

  // compare orbits by first comparing invariants, then orbit prototypes
  bool inter_orbit_compare(Eigen::VectorXcd const &A_prototype,
                           VectorInvariants const &A_invariants,
                           Eigen::VectorXcd const &B_prototype,
                           VectorInvariants const &B_invariants) const;

 private:
  MatrixRep const &m_matrix_rep;
  double m_tol;
};

/// Vector space preparation for comparison
Eigen::MatrixXcd vector_space_prepare(Eigen::MatrixXcd const &vector_space,
                                      double tol);

/// Vector space preparation for comparison
Eigen::MatrixXd vector_space_prepare(Eigen::MatrixXd const &vector_space,
                                     double tol);

}  // namespace SymRepTools_v2

}  // namespace CASM

#endif
