#include "casm/crystallography/Lattice.hh"
#include "casm/crystallography/Strain.hh"
#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/strain/StrainConverter.hh"
#include "gtest/gtest.h"

// StrainConverter
// ---------------
//
// StrainConverter enables converting between the deformation tensor and various
// strain metrics, and unrolling the strain metric as a vector. The
// StrainConverter is constructed with "mode_name" (std::string) indicating one
// of several possible strain metrics (strain::METRIC):
// - "GL": strain::METRIC::GREEN_LAGRANGE, E = 1/2 * (F^{T} F - I)
// - "B": strain::METRIC::BIOT, E = (U-I)
// - "H": strain::METRIC::HENCKY, E = log(C)/2
// - "EA": strain::METRIC::EULER_ALMANSI, E = 0.5 * (I-(F F^{T})^(-1))
// - "U": strain::METRIC::STRETCH, E = U
// - "F": strain::METRIC::DISP_GRAD, E = F
//
// Naming conventions:
// - F: deformation tensor (F = R * U)
// - E: strain metric (calculated based on what "mode_name" / strain::METRIC is)
// - R: rotation tensor
// - U: stretch tensor
// - C: F^{T} * F

namespace {

template <CASM::strain::METRIC METRIC_TYPE>
void check_conversions(Eigen::Matrix3d const &F) {
  // Convert from F -> E
  Eigen::Matrix3d strain_metric =
      CASM::strain::deformation_tensor_to_metric<METRIC_TYPE>(F);
  std::cout << "strain_metric: \n" << strain_metric << std::endl << std::endl;

  // Convert from E -> F
  Eigen::Matrix3d deformation_tensor =
      CASM::strain::metric_to_deformation_tensor<METRIC_TYPE>(strain_metric);

  EXPECT_TRUE(almost_equal(deformation_tensor, F, CASM::TOL));
}
}  // namespace

TEST(ExampleCrystallographyStrain, StrainConverter) {
  // A cubic lattice (a, b, c)
  Eigen::Vector3d a{1.0, 0.0, 0.0};
  Eigen::Vector3d b{0.0, 1.0, 0.0};
  Eigen::Vector3d c{0.0, 0.0, 1.0};
  CASM::xtal::Lattice ideal_lattice{a, b, c};
  Eigen::Matrix3d L = ideal_lattice.lat_column_mat();

  // Deformed lattice
  Eigen::Vector3d da{0.01, 0.0, 0.0};
  Eigen::Vector3d db{0.0, 0.0, 0.0};
  Eigen::Vector3d dc{0.0, 0.0, 0.0};
  CASM::xtal::Lattice deformed_lattice{a + da, b + db, c + dc};
  Eigen::Matrix3d D = deformed_lattice.lat_column_mat();

  // Deformation tensor, F:
  //     deformed_lattice.lat_column_mat() = F * ideal_lattice.lat_column_mat()
  //     D = F * L  ->  D.t = L.t * F.t
  Eigen::Matrix3d F =
      L.transpose().fullPivHouseholderQr().solve(D.transpose()).transpose();
  std::cout << "F: \n" << F << std::endl;
  EXPECT_TRUE(almost_equal(D, F * L, CASM::TOL));

  // Example strain conversions between F and E using:
  // - CASM::strain::deformation_tensor_to_metric<METRIC_TYPE>
  // - CASM::strain::metric_to_deformation_tensor<METRIC_TYPE>
  check_conversions<CASM::strain::METRIC::GREEN_LAGRANGE>(F);
  check_conversions<CASM::strain::METRIC::BIOT>(F);
  check_conversions<CASM::strain::METRIC::HENCKY>(F);
  check_conversions<CASM::strain::METRIC::EULER_ALMANSI>(F);

  // CASM represents all property values as vectors.  For strain metric
  // properties, unroll as a vector (and convert back) using StrainConverter:
  CASM::StrainConverter strain_converter{"GL"};
  Eigen::VectorXd unrolled_strain_metric =
      strain_converter.unrolled_strain_metric(F);
  Eigen::Matrix3d F2 =
      strain_converter.unrolled_strain_metric_to_F(unrolled_strain_metric);
  EXPECT_TRUE(almost_equal(F2, F, CASM::TOL));
}
