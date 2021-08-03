#include "casm/monte2/Conversions.hh"

#include "casm/clex/Configuration.hh"
#include "casm/clex/Supercell.hh"
#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/misc/algorithm.hh"

namespace CASM {
namespace Monte2 {

Conversions::Conversions(std::shared_ptr<Supercell const> const &supercell)
    : Conversions(Configuration{supercell}, supercell) {}

Conversions::Conversions(Configuration const &unit_config,
                         std::shared_ptr<Supercell const> const &supercell)
    : m_unit_supercell(unit_config.shared_supercell()),
      m_unitl_and_bijk_converter(unit_config.shared_supercell()
                                     ->sym_info()
                                     .unitcellcoord_index_converter()),
      m_supercell(supercell),
      m_l_and_bijk_converter(
          supercell->sym_info().unitcellcoord_index_converter()),
      m_struc_mol(xtal::struc_molecule(m_supercell->prim())),
      m_struc_molname(xtal::struc_molecule_name(m_supercell->prim())) {
  // make m_unitl_to_asym, m_Nasym
  Index asym = 0;
  Index unit_Nsites = m_unit_supercell->num_sites();
  m_unitl_to_asym.resize(unit_Nsites, -1);
  std::vector<PermuteIterator> fg = unit_config.factor_group();
  for (Index unitl = 0; unitl < unit_Nsites; ++unitl) {
    if (m_unitl_to_asym[unitl] == -1) {
      for (auto it = fg.begin(); it != fg.end(); ++it) {
        // permutation defined by: after[i] = before[it->permute_ind(i)]
        m_unitl_to_asym[it->permute_ind(unitl)] = asym;
      }
      ++asym;
    }
  }
  m_Nasym = asym;

  // make m_asym_to_unitl & m_asym_to_b
  m_asym_to_unitl.resize(m_Nasym);
  m_asym_to_b.resize(m_Nasym);
  for (Index unitl = 0; unitl < unit_Nsites; ++unitl) {
    Index asym = m_unitl_to_asym[unitl];
    m_asym_to_unitl[asym].insert(unitl);
    m_asym_to_b[asym].insert(unitl_to_b(unitl));
  }

  // make m_occ_to_species and m_species_to_occ

  // [b][occ] -> species
  auto index_converter =
      make_index_converter(m_supercell->prim(), m_struc_molname);

  // [b][species] -> occ, index_converter[b].size() if not allowed
  auto index_converter_inv =
      make_index_converter_inverse(m_supercell->prim(), m_struc_molname);

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

std::shared_ptr<Supercell const> const &Conversions::unit_supercell() const {
  return m_unit_supercell;
}

std::shared_ptr<Supercell const> const &Conversions::supercell() const {
  return m_supercell;
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
