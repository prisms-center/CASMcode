#include "casm/monte_carlo/Conversions.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/Supercell.hh"
#include "casm/clex/Configuration.hh"
#include "casm/crystallography/Structure.hh"

namespace CASM {
  namespace Monte {

    namespace {
      Configuration default_prim_config(PrimClex &primclex) {
        auto &scel = primclex.get_supercell(primclex.add_supercell(primclex.get_prim().lattice()));
        Configuration res(scel);
        res.init_occupation();
        return res;
      }
    }

    Conversions::Conversions(const Supercell &mc_scel) :
      Conversions(default_prim_config(mc_scel.get_primclex()), mc_scel) {}


    Conversions::Conversions(const Configuration &unit_config, const Supercell &mc_scel) :
      m_unit_scel(&unit_config.get_supercell()),
      m_mc_scel(&mc_scel),
      m_struc_mol(m_mc_scel->get_prim().get_struc_molecule()),
      m_struc_molname(m_mc_scel->get_prim().get_struc_molecule_name()) {

      // make m_unitl_to_asym, m_Nasym
      Index asym = 0;
      Index unit_Nsites = m_unit_scel->num_sites();
      m_unitl_to_asym.resize(unit_Nsites, -1);
      std::vector<PermuteIterator> fg = unit_config.factor_group();
      for(Index unitl = 0; unitl < unit_Nsites; ++unitl) {
        if(m_unitl_to_asym[unitl] == -1) {
          for(auto it = fg.begin(); it != fg.end(); ++it) {
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
      for(Index unitl = 0; unitl < unit_Nsites; ++unitl) {
        Index asym = m_unitl_to_asym[unitl];
        m_asym_to_unitl[asym].insert(unitl);
        m_asym_to_b[asym].insert(unitl_to_b(unitl));
      }

      // make m_occ_to_species and m_species_to_occ

      // [b][occ] -> species
      auto index_converter = get_index_converter(m_mc_scel->get_prim(), m_struc_molname);

      // [b][species] -> occ, index_converter[b].size() if not allowed
      auto index_converter_inv = get_index_converter_inverse(m_mc_scel->get_prim(), m_struc_molname);

      m_occ_to_species.resize(m_Nasym);
      m_species_to_occ.resize(m_Nasym);
      for(Index asym = 0; asym < m_Nasym; ++asym) {
        Index b = *(m_asym_to_b[asym].begin());
        m_occ_to_species[asym] = index_converter[b];
        m_species_to_occ[asym] = index_converter_inv[b];
      }
    }

    Index Conversions::l_to_b(Index l) const {
      return m_mc_scel->get_b(l);
    }
    UnitCell Conversions::l_to_ijk(Index l) const {
      return m_mc_scel->prim_grid().unitcell(l % m_mc_scel->volume());
    }
    UnitCellCoord Conversions::l_to_bijk(Index l) const {
      return m_mc_scel->uccoord(l);
    }
    Index Conversions::l_to_unitl(Index l) const {
      return bijk_to_unitl(l_to_bijk(l));
    }
    Index Conversions::l_to_asym(Index l) const {
      return m_unitl_to_asym[l_to_unitl(l)];
    }

    Index Conversions::bijk_to_l(const UnitCellCoord &bijk) const {
      return m_mc_scel->find(bijk);
    }
    Index Conversions::bijk_to_unitl(const UnitCellCoord &bijk) const {
      return m_unit_scel->find(bijk);
    }
    Index Conversions::bijk_to_asym(const UnitCellCoord &bijk) const {
      return l_to_asym(bijk_to_l(bijk));
    }

    Index Conversions::unitl_to_b(Index unitl) const {
      return m_unit_scel->get_b(unitl);
    }
    UnitCellCoord Conversions::unitl_to_bijk(Index unitl) const {
      return m_unit_scel->uccoord(unitl);
    }
    Index Conversions::unitl_to_asym(Index unitl) const {
      return m_unitl_to_asym[unitl];
    }

    Index Conversions::asym_size() const {
      return m_Nasym;
    }
    const std::set<Index> &Conversions::asym_to_b(Index asym) const {
      return m_asym_to_b[asym];
    }
    const std::set<Index> &Conversions::asym_to_unitl(Index asym) const {
      return m_asym_to_unitl[asym];
    }

    const Supercell &Conversions::unit_scel() const {
      return *m_unit_scel;
    }
    const Supercell &Conversions::mc_scel() const {
      return *m_mc_scel;
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

    Index Conversions::species_size() const {
      return m_struc_mol.size();
    }
    Index Conversions::species_index(std::string species_name) const {
      return find_index(m_struc_molname, species_name);
    }
    const Molecule &Conversions::species_to_mol(Index species_index) const {
      return m_struc_mol[species_index];
    }
    const std::string &Conversions::species_name(Index species_index) const {
      return m_struc_molname[species_index];
    }
    Index Conversions::components_size(Index species_index) const {
      return species_to_mol(species_index).size();
    }

  }
}
