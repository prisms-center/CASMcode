#include <casm/crystallography/Superlattice.hh>

namespace CASM {
namespace xtal {
Superlattice::Superlattice(const Lattice &tiling_unit,
                           const Lattice &superlattice)
    : m_primitive_lattice(tiling_unit),
      m_superlattice(superlattice),
      m_transformation_matrix_to_super(
          xtal::make_transformation_matrix_to_super(this->prim_lattice(),
                                                    this->superlattice(), TOL)),
      m_size(std::abs(m_transformation_matrix_to_super.determinant())) {}

Superlattice::Superlattice(const Lattice &tiling_unit,
                           const Eigen::Matrix3l &transformation_matrix)
    : m_primitive_lattice(tiling_unit),
      m_superlattice(
          xtal::make_superlattice(tiling_unit, transformation_matrix)),
      m_transformation_matrix_to_super(transformation_matrix),
      m_size(std::abs(m_transformation_matrix_to_super.determinant())) {}

Superlattice Superlattice::smooth_prim(const Lattice &tiling_unit,
                                       const Lattice &superlattice) {
  Superlattice result_superlattice(tiling_unit, superlattice);

  Eigen::Matrix3d trans_mat_inverse =
      result_superlattice.transformation_matrix_to_super()
          .cast<double>()
          .inverse();
  // By using the inverse transformatiom matrix, the new prim is guaranteed to
  // perfectly tile the old prim
  Lattice reconstructed_tiling_unit(
      result_superlattice.superlattice().lat_column_mat() * trans_mat_inverse,
      result_superlattice.prim_lattice().tol());

  // Pry right into the primate members and set the proper tiling unit
  result_superlattice.m_primitive_lattice = reconstructed_tiling_unit;
  return result_superlattice;
}

}  // namespace xtal
}  // namespace CASM
