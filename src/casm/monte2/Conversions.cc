#include "casm/monte2/Conversions.hh"

#include "casm/clex/Configuration.hh"
#include "casm/clex/Supercell.hh"
#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/BasicStructureTools.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/misc/algorithm.hh"

namespace CASM {
namespace Monte2 {

namespace {
/// Return indices of equivalent basis sites
std::vector<Index> make_b_to_asym(const xtal::BasicStructure &struc) {
  std::vector<Index> b_to_asym(struc.basis().size());
  std::set<std::set<Index>> asym_unit = xtal::make_asymmetric_unit(struc);
  Index asym = 0;
  for (auto const &orbit : asym_unit) {
    for (Index b : orbit) {
      b_to_asym[b] = asym;
    }
    ++asym;
  }
  return b_to_asym;
}
}  // namespace

/// \brief Constructor (uses asymmetric unit determined from prim factor group)
///
/// \param prim The primitive structure
/// \param transformation_matrix_to_super Defines a supercell lattice,
///     S = P * T, where S = supercell lattice column matrix, P = prim lattice
///     column matrix, T = transformation_matrix_to_super.
///
/// This overload uses the prim factor group symmetry to determine the
/// asymmetric unit.
Conversions::Conversions(xtal::BasicStructure const &prim,
                         Eigen::Matrix3l const &transformation_matrix_to_super)
    : Conversions(prim, transformation_matrix_to_super,
                  Eigen::Matrix3l::Identity(), make_b_to_asym(prim)) {}

/// \brief Constructor (user specified asymmetric unit with reduced symmetry)
///
/// \param prim The primitive structure
/// \param transformation_matrix_to_super Defines a supercell lattice,
///     S = P * T, where S = supercell lattice column matrix, P = prim lattice
///     column matrix, T = transformation_matrix_to_super.
/// \param b_to_asym Specifies the asymmetric unit orbit index
///     corresponding to each sublattice in the prim. Asymmetric unit orbit
///     indices are distinct indices `(0, 1, ...)` indicating that sites with
///     the same index map onto each other via symmetry operations.
///
/// This overload allows specifying lower symmetry than the prim factor group
/// (but same periodicity) to determine the asymmetric unit.
///
Conversions::Conversions(xtal::BasicStructure const &prim,
                         Eigen::Matrix3l const &transformation_matrix_to_super,
                         std::vector<Index> const &b_to_asym)
    : Conversions(prim, transformation_matrix_to_super,
                  Eigen::Matrix3l::Identity(), b_to_asym) {}

/// \brief Constructor (user specified asymmetric unit with reduced
/// translational symmetry)
///
/// \param prim The primitive structure
/// \param transformation_matrix_to_super Defines a supercell lattice,
///     S = P * T, where S = supercell lattice column matrix, P = prim lattice
///     column matrix, T = transformation_matrix_to_super.
/// \param unit_transformation_matrix_to_super Defines a sub-supercell lattice,
///     U = P * T, where U = supercell lattice column matrix, P = prim lattice
///     column matrix, T_unit = unit_transformation_matrix_to_super. U must
///     tile into S (i.e. S = U * T', where T' is an integer matrix). Allows
///     specifying an asymmetric unit which does not fit in the primitive cell.
/// \param unitl_to_asym Specifies the asymmetric unit orbit index
///     corresponding to each site in the supercell U. Asymmetric unit orbit
///     indices are distinct indices `(0, 1, ...)` indicating that sites with
///     the same index map onto each other via symmetry operations.
///
/// This overload allows specifying an asymmetric unit which does not fit in
/// the primitive cell.
///
Conversions::Conversions(
    xtal::BasicStructure const &prim,
    Eigen::Matrix3l const &transformation_matrix_to_super,
    Eigen::Matrix3l const &unit_transformation_matrix_to_super,
    std::vector<Index> const &unitl_to_asym)
    : m_unit_transformation_matrix_to_super(
          unit_transformation_matrix_to_super),
      m_unitl_and_bijk_converter(unit_transformation_matrix_to_super,
                                 prim.basis().size()),
      m_transformation_matrix_to_super(transformation_matrix_to_super),
      m_l_and_bijk_converter(transformation_matrix_to_super,
                             prim.basis().size()),
      m_struc_mol(xtal::struc_molecule(prim)),
      m_struc_molname(xtal::struc_molecule_name(prim)),
      m_unitl_to_asym(unitl_to_asym) {
  // find m_Nasym
  m_Nasym =
      *std::max_element(m_unitl_to_asym.begin(), m_unitl_to_asym.end()) + 1;

  // make m_asym_to_unitl & m_asym_to_b
  Index unit_Nsites =
      unit_transformation_matrix_to_super.determinant() * prim.basis().size();
  m_asym_to_unitl.resize(m_Nasym);
  m_asym_to_b.resize(m_Nasym);
  for (Index unitl = 0; unitl < unit_Nsites; ++unitl) {
    Index asym = m_unitl_to_asym[unitl];
    m_asym_to_unitl[asym].insert(unitl);
    m_asym_to_b[asym].insert(unitl_to_b(unitl));
  }

  // make m_occ_to_species and m_species_to_occ

  // [b][occ] -> species
  auto index_converter = make_index_converter(prim, m_struc_molname);

  // [b][species] -> occ, index_converter[b].size() if not allowed
  auto index_converter_inv =
      make_index_converter_inverse(prim, m_struc_molname);

  m_occ_to_species.resize(m_Nasym);
  m_species_to_occ.resize(m_Nasym);
  for (Index asym = 0; asym < m_Nasym; ++asym) {
    Index b = *(m_asym_to_b[asym].begin());
    m_occ_to_species[asym] = index_converter[b];
    m_species_to_occ[asym] = index_converter_inv[b];
  }
}

Index Conversions::l_to_b(Index l) const {
  return m_l_and_bijk_converter(l).sublattice();
}

xtal::UnitCell Conversions::l_to_ijk(Index l) const {
  return m_l_and_bijk_converter(l).unitcell();
}

xtal::UnitCellCoord Conversions::l_to_bijk(Index l) const {
  return m_l_and_bijk_converter(l);
}

Index Conversions::l_to_unitl(Index l) const {
  return bijk_to_unitl(l_to_bijk(l));
}

Index Conversions::l_to_asym(Index l) const {
  return m_unitl_to_asym[l_to_unitl(l)];
}

Index Conversions::bijk_to_l(xtal::UnitCellCoord const &bijk) const {
  return m_l_and_bijk_converter(bijk);
}
Index Conversions::bijk_to_unitl(xtal::UnitCellCoord const &bijk) const {
  return m_unitl_and_bijk_converter(bijk);
}

Index Conversions::bijk_to_asym(xtal::UnitCellCoord const &bijk) const {
  return l_to_asym(bijk_to_l(bijk));
}

Index Conversions::unitl_to_b(Index unitl) const {
  return m_unitl_and_bijk_converter(unitl).sublattice();
}

xtal::UnitCellCoord Conversions::unitl_to_bijk(Index unitl) const {
  return m_unitl_and_bijk_converter(unitl);
}

Index Conversions::unitl_to_asym(Index unitl) const {
  return m_unitl_to_asym[unitl];
}

Index Conversions::asym_size() const { return m_Nasym; }

std::set<Index> const &Conversions::asym_to_b(Index asym) const {
  return m_asym_to_b[asym];
}

std::set<Index> const &Conversions::asym_to_unitl(Index asym) const {
  return m_asym_to_unitl[asym];
}

Index Conversions::occ_size(Index asym) const {
  return m_occ_to_species[asym].size();
}

Index Conversions::species_index(Index asym, Index occ_index) const {
  return m_occ_to_species[asym][occ_index];
}

Index Conversions::occ_index(Index asym, Index species_index) const {
  // returns occ_size(asym) if species not allowed
  return m_species_to_occ[asym][species_index];
}

bool Conversions::species_allowed(Index asym, Index species_index) const {
  return occ_index(asym, species_index) != occ_size(asym);
}

Index Conversions::species_size() const { return m_struc_mol.size(); }

Index Conversions::species_index(std::string species_name) const {
  return find_index(m_struc_molname, species_name);
}
xtal::Molecule const &Conversions::species_to_mol(Index species_index) const {
  return m_struc_mol[species_index];
}
std::string const &Conversions::species_name(Index species_index) const {
  return m_struc_molname[species_index];
}
Index Conversions::components_size(Index species_index) const {
  return species_to_mol(species_index).size();
}

}  // namespace Monte2
}  // namespace CASM
