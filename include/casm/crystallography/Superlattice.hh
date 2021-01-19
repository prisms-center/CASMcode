#ifndef XTAL_SUPERLATTICE_HH
#define XTAL_SUPERLATTICE_HH

#include <casm/crystallography/Lattice.hh>

namespace CASM {
namespace xtal {
class Superlattice {
 private:
  Lattice m_primitive_lattice;
  Lattice m_superlattice;
  /// Integer matrix that convers the primitive lattice into the superlattice
  Eigen::Matrix3l m_transformation_matrix_to_super;

  Index m_size;

 public:
  const Lattice &superlattice() const { return m_superlattice; }

  const Lattice &prim_lattice() const { return m_primitive_lattice; }

  /// The integer transformation matrix that converts the tiling unit (primitive
  /// lattice) into the superlattice
  const Eigen::Matrix3l &transformation_matrix_to_super() const {
    return m_transformation_matrix_to_super;
  }

  /// Returns the number of tiling units (primitive lattices) that fit inside
  /// the superlattice
  Index size() const { return m_size; }

  Superlattice(const Lattice &tiling_unit, const Lattice &superlattice);
  Superlattice(const Lattice &tiling_unit,
               const Eigen::Matrix3l &transformation_matrix);

  /// Constructs by taking the superlattice vectors to be the correct values,
  /// and then uses the rounded integer transformation matrix to slightly modify
  /// the primitive lattice so that it exactly transforms under the
  /// transformation matrix
  static Superlattice smooth_prim(const Lattice &tiling_unit,
                                  const Lattice &superlattice);

  // TODO: Implement the day you need it.
  /// Constructs by taking the primitive vectors to be the correct values, and
  /// then uses the rounded integer transformation matrix to slightly modify the
  /// superlattice so that the primitive lattice exactly transforms to it
  /// through the transformation matrix
  /* static Superlattice smooth_superlattice(const Lattice &tiling_unit, const
   * Lattice &superlattice); */
};

}  // namespace xtal
}  // namespace CASM

#endif
