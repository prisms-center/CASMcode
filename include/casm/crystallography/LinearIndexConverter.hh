#ifndef LINEARINDEXCONVERTER_HH
#define LINEARINDEXCONVERTER_HH

#include <unordered_map>
#include <vector>

#include "casm/crystallography/IntegralCoordinateWithin.hh"
#include "casm/crystallography/UnitCellCoord.hh"

namespace CASM {
namespace xtal {

/**
 * Converts back and forth between UnitCellCoord and its linear index,
 * where the linear index is guaranteed to preserve order based on the
 * sublattice index of the UnitCellCoord, and the Smith Normal Form
 * of the UnitCell.
 *
 * By default, the constructed index converter will always accept
 * UnitCellCoord values that fall outside of the superlattice. When this
 * happens, the given UnitCellCoord is brought into the superlattice
 * using superlattice vector translations, and the index of the resulting
 * UnitCellCoord is returned.
 */

class UnitCellCoordIndexConverter {
 public:
  typedef impl::OrderedLatticePointGenerator::matrix_type matrix_type;

  /// Initialize with the transformation that defines how to convert from the
  /// tiling unit (prim) to the superlattice, and the number of basis sites in
  /// the primitive cell
  UnitCellCoordIndexConverter(const matrix_type &transformation_matrix,
                              int basis_sites_in_prim)
      : m_linear_index_to_bijk(_make_all_ordered_bijk_values(
            impl::OrderedLatticePointGenerator(transformation_matrix),
            basis_sites_in_prim)),
        m_basis_sites_in_prim(basis_sites_in_prim),
        m_automatically_bring_bijk_within(true),
        m_bring_within_f(transformation_matrix) {
    _throw_if_bad_basis_sites_in_prim(m_basis_sites_in_prim);

    for (Index ix = 0; ix < m_linear_index_to_bijk.size(); ++ix) {
      const auto &bijk = m_linear_index_to_bijk[ix];
      m_bijk_to_linear_index[bijk] = ix;
    }

    // By default the construction will always allow receiving bijk values that
    // fall outside of the superlattice
    this->always_bring_within();
  }

  /// Prevent the index converter from bringing UnitCellCoord within the
  /// supercell when querying for the index
  void never_bring_within();

  /// Automatically bring UnitCellCoord values within the supercell when
  /// querying for the index (on by default)
  void always_bring_within();

  /// Bring the given UnitCellCoord into the superlattice using lattice
  /// translations
  UnitCellCoord bring_within(const UnitCellCoord &bijk) const;

  /// Given the linear index, retreive the corresponding UnitCellCoord
  const UnitCellCoord &operator()(Index ix) const;

  /// Given the UnitCellCoord, retreive its corresponding linear index.
  /// If applicable, brings the UnitCellCoord within the superlattice
  Index operator()(const UnitCellCoord &bijk) const;

  /// Returns the total number of sites within the superlattice
  Index total_sites() const { return m_linear_index_to_bijk.size(); }

 private:
  /// Convert from linear index to UnitCellCoord
  std::vector<UnitCellCoord> m_linear_index_to_bijk;

  /// Convert from UnitCellCoord to linear index
  std::unordered_map<UnitCellCoord, Index> m_bijk_to_linear_index;

  /// Stores a cache of UnitCellCoord values that landed outside of the
  /// superlattice
  mutable std::unordered_map<UnitCellCoord, Index>
      m_bijk_to_linear_index_outside_of_superlattice;

  /// How many blocks of "b", i.e. number of atoms in the primitive cell, as
  /// specified at construction
  int m_basis_sites_in_prim;

  /// If set to true, UnitCellCoord values will be brought within the supercell
  /// before querying for the index
  bool m_automatically_bring_bijk_within;

  /// Functor to bring UnitCellCoord values back into the superlattice
  IntegralCoordinateWithin_f m_bring_within_f;

  /// Throws exception if the specified index is out of the allowed range
  void _throw_if_incompatible_index(Index ix) const;

  /// Throws exception if the specified UnitCellCoord has a sublattice index
  /// that isn't compatible. If the state is set to not automatically bring the
  /// UnitCellCoord within the superlattice, then any UnitCellCoord outside the
  /// boundary will also trigger an exception
  void _throw_if_incompatible_bijk(const UnitCellCoord &bijk) const;

  /// Enumerates every possible UnitCellCoord and returns them in the expected
  /// order (blocks by basis site, with the Smith Normal Form order within each
  /// block
  static std::vector<UnitCellCoord> _make_all_ordered_bijk_values(
      const impl::OrderedLatticePointGenerator &make_point,
      int basis_sites_in_prim);

  /// Throws exception if the number of sites specified in the tiling unit is
  /// less than 1
  static void _throw_if_bad_basis_sites_in_prim(int basis_sites_in_prim);
};

/**
 * Converts back and forth between UnitCelland its linear index,
 * where the linear index is guaranteed to preserve order based on the
 * and the Smith Normal Form of the UnitCell.
 *
 * By default, the constructed index converter will always accept
 * UnitCell values that fall outside of the superlattice. When this
 * happens, the given UnitCellCoord is brought into the superlattice
 * using superlattice vector translations, and the index of the resulting
 * UnitCell is returned.
 */

class UnitCellIndexConverter : UnitCellCoordIndexConverter {
 public:
  typedef UnitCellCoordIndexConverter::matrix_type matrix_type;
  UnitCellIndexConverter(const matrix_type &transformation_matrix)
      : UnitCellCoordIndexConverter(transformation_matrix, 1) {}

  using UnitCellCoordIndexConverter::always_bring_within;
  using UnitCellCoordIndexConverter::never_bring_within;
  using UnitCellCoordIndexConverter::total_sites;

  /// Bring the given UnitCell into the superlattice using lattice translations
  UnitCell bring_within(const UnitCell &ijk) const {
    return UnitCellCoordIndexConverter::bring_within(UnitCellCoord(0, ijk))
        .unitcell();
  }

  /// Given the linear index, retreive the corresponding UnitCell
  const UnitCell &operator()(Index ix) const {
    return UnitCellCoordIndexConverter::operator()(ix).unitcell();
  }

  /// Given the UnitCell, retreive its corresponding linear index.
  /// If applicable, brings the UnitCell within the superlattice
  Index operator()(const UnitCell &ijk) const {
    return UnitCellCoordIndexConverter::operator()(UnitCellCoord(0, ijk));
  }

 private:
  // TODO:
  // Secretly cache things to "speed things up" if you're *really* that
  // concerned about the UnitCellCoord copies slowing you down. Don't implement
  // this until you've profiled it.
};

/**
 * These templates are just for syntactic convenience.
 * The template class LinearIndexConverter can be used for both
 * UnitCell and UnitCellCoord types, by specializing the
 * desired type as a template argument.
 */

template <typename IntegralCoordinateType>
struct _LinearIndexConverter {};

template <>
struct _LinearIndexConverter<UnitCell> {
  typedef UnitCellIndexConverter type;
};

template <>
struct _LinearIndexConverter<UnitCellCoord> {
  typedef UnitCellCoordIndexConverter type;
};

template <typename IntegralCoordinateType>
using LinearIndexConverter =
    typename _LinearIndexConverter<IntegralCoordinateType>::type;

}  // namespace xtal
}  // namespace CASM

#endif
