#include "casm/crystallography/Strain.hh"

#include "casm/external/Eigen/Core"

namespace CASM {
namespace strain {

template <>
Eigen::Matrix3d deformation_tensor_to_metric<METRIC::GREEN_LAGRANGE>(
    const Eigen::Ref<const Eigen::Matrix3d> &deformation_tensor) {
  const Eigen::Matrix3d &F = deformation_tensor;
  return 0.5 * (F.transpose() * F - Eigen::MatrixXd::Identity(3, 3));
}

template <>
Eigen::Matrix3d metric_to_deformation_tensor<METRIC::GREEN_LAGRANGE>(
    const Eigen::Ref<const Eigen::Matrix3d> &metric_tensor) {
  const Eigen::Matrix3d &E = metric_tensor;
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(
      2 * E + Eigen::MatrixXd::Identity(3, 3));
  return es.operatorSqrt();
}

//***********************************************************//

template <>
Eigen::Matrix3d deformation_tensor_to_metric<METRIC::BIOT>(
    const Eigen::Ref<const Eigen::Matrix3d> &deformation_tensor) {
  const Eigen::Matrix3d &F = deformation_tensor;
  return (right_stretch_tensor(F) - Eigen::MatrixXd::Identity(3, 3));
}

template <>
Eigen::Matrix3d metric_to_deformation_tensor<METRIC::BIOT>(
    const Eigen::Ref<const Eigen::Matrix3d> &metric_tensor) {
  const Eigen::Matrix3d &B = metric_tensor;
  return B + Eigen::MatrixXd::Identity(3, 3);
}

//***********************************************************//

template <>
Eigen::Matrix3d deformation_tensor_to_metric<METRIC::HENCKY>(
    const Eigen::Ref<const Eigen::Matrix3d> &deformation_tensor) {
  const Eigen::Matrix3d &F = deformation_tensor;
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(F.transpose() * F);
  return es.eigenvectors() *
         es.eigenvalues().array().log().matrix().asDiagonal() *
         es.eigenvectors().inverse() / 2.0;
}

template <>
Eigen::Matrix3d metric_to_deformation_tensor<METRIC::HENCKY>(
    const Eigen::Ref<const Eigen::Matrix3d> &metric_tensor) {
  const Eigen::Matrix3d &H = metric_tensor;
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(H);
  return es.eigenvectors() *
         es.eigenvalues().array().exp().matrix().asDiagonal() *
         es.eigenvectors().inverse();
}

//***********************************************************//

template <>
Eigen::Matrix3d deformation_tensor_to_metric<METRIC::EULER_ALMANSI>(
    const Eigen::Ref<const Eigen::Matrix3d> &deformation_tensor) {
  const Eigen::Matrix3d &F = deformation_tensor;
  return 0.5 *
         (Eigen::MatrixXd::Identity(3, 3) - (F * F.transpose()).inverse());
}

template <>
Eigen::Matrix3d metric_to_deformation_tensor<METRIC::EULER_ALMANSI>(
    const Eigen::Ref<const Eigen::Matrix3d> &metric_tensor) {
  const Eigen::Matrix3d &A = metric_tensor;
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(
      Eigen::MatrixXd::Identity(3, 3) - 2 * A);
  return es.operatorInverseSqrt();
}

//***********************************************************//

/* template <> */
/* Eigen::Matrix3d deformation_tensor_to_metric<METRIC::DISP_GRAD>(const
 * Eigen::Ref<const Eigen::Matrix3d>& deformation_tensor) */
/* { */
/*     const Eigen::Matrix3d& F=deformation_tensor; */
/* } */

/* template <> */
/* Eigen::Matrix3d metric_to_deformation_tensor<METRIC::DISP_GRAD>(const
 * Eigen::Ref<const Eigen::Matrix3d>& metric_tensor) */
/* { */
/* } */

//***********************************************************//

Eigen::Matrix3d metric_tensor(
    const Eigen::Ref<const Eigen::Matrix3d> &deformation_tensor) {
  const Eigen::Matrix3d &F = deformation_tensor;
  return F.transpose() * F;
}

Eigen::Matrix3d right_stretch_tensor(
    const Eigen::Ref<const Eigen::Matrix3d> &deformation_tensor) {
  const Eigen::Matrix3d &F = deformation_tensor;
  Eigen::Matrix3d C = metric_tensor(F);
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigen_solver(C);
  return eigen_solver.operatorSqrt();
}

}  // namespace strain
}  // namespace CASM
