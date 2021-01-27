#ifndef XTAL_STRAIN_HH
#define XTAL_STRAIN_HH

#include "casm/external/Eigen/Dense"
#include "casm/external/Eigen/src/Core/Matrix.h"

namespace CASM {
namespace strain {
/**
 * Enumerates different types of available strain metrics, which can be used to
 * specify conversion between deformation tensors and the desired metric:
 *
 * F: deformation tensor (F = R * U)
 * R: rotation tensor
 * U: stretch tensor
 * C: F^{T} * F
 *
 * GREEN_LAGRANGE = 1/2 * (F^{T} F - I)
 * BIOT = (U-I)
 * HENCKY = log(C)/2
 * EULER_ALMANSI = (I-(F F^{T})^(-1))/2
 */

enum class METRIC {
  GREEN_LAGRANGE,
  BIOT,
  HENCKY,
  EULER_ALMANSI,
  STRETCH,
  DISP_GRAD
};

/// Convert a deformation tensor to any metric tensor, specified by the template
/// parameter. For example, in order to turn a deformation tensor F into a Green
/// Lagrange metric tensor GL:
/// Eigen::Matrix3d GL=deformation_tensor_to_metric<METRIC::GREEN_LAGRANGE>(F);
template <METRIC>
Eigen::Matrix3d deformation_tensor_to_metric(
    const Eigen::Ref<const Eigen::Matrix3d> &deformation_tensor);

/// Convert a strain metric tensor to the deformation tensor, where the given
/// strain metric is specified with the template parameter. For example, to
/// convert a Green Lagrange strain tensor GL to the deformation tensor F:
/// Eigen::Matrix3d F=metric_to_deformation_tensor<METRIC::GREEN_LAGRANGE>(GL);
template <METRIC>
Eigen::Matrix3d metric_to_deformation_tensor(
    const Eigen::Ref<const Eigen::Matrix3d> &metric_tensor);

/// Calculates the metric tensor of the the deformation gradient as
/// F^{T}F
Eigen::Matrix3d metric_tensor(
    const Eigen::Ref<const Eigen::Matrix3d> &deformation_tensor);

/// Calculates and returns the value of U where F = R*U
Eigen::Matrix3d right_stretch_tensor(
    const Eigen::Ref<const Eigen::Matrix3d> &deformation_tensor);
}  // namespace strain
}  // namespace CASM

#endif
