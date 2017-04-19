#ifndef CASM_Monte_Conversions_HH
#define CASM_Monte_Conversions_HH

#include <vector>
#include <set>
#include "casm/CASM_global_definitions.hh"

namespace CASM {

  class PrimClex;
  class Supercell;
  class Configuration;
  class UnitCell;
  class UnitCellCoord;
  class Molecule;

  namespace Monte {

    /// l: linear index into mc_scel
    /// b: prim basis site index
    /// ijk: prim unit cell indices
    /// bijk: prim basis site index + unit cell indices
    /// unitl: ref config basis site index
    /// asym: asymmetric unit index
    class Conversions {

    public:

      Conversions(const Supercell &mc_scel);

      Conversions(const Configuration &unit_config, const Supercell &mc_scel);


      Index l_to_b(Index l) const;
      UnitCell l_to_ijk(Index l) const;
      UnitCellCoord l_to_bijk(Index l) const;
      Index l_to_unitl(Index l) const;
      Index l_to_asym(Index l) const;

      Index bijk_to_l(const UnitCellCoord &bijk) const;
      Index bijk_to_unitl(const UnitCellCoord &bijk) const;
      Index bijk_to_asym(const UnitCellCoord &bijk) const;

      Index unitl_to_b(Index unitl) const;
      UnitCellCoord unitl_to_bijk(Index unitl) const;
      Index unitl_to_asym(Index unitl) const;

      Index asym_size() const;
      const std::set<Index> &asym_to_b(Index asym) const;
      const std::set<Index> &asym_to_unitl(Index asym) const;

      const Supercell &unit_scel() const;
      const Supercell &mc_scel() const;

      Index occ_size(Index asym) const;
      Index species_index(Index asym, Index occ_index) const;
      Index occ_index(Index asym, Index species_index) const;
      bool species_allowed(Index asym, Index species_index) const;

      Index species_size() const;
      Index species_index(std::string species_name) const;
      const Molecule &species_to_mol(Index species_index) const;
      const std::string &species_name(Index species_index) const;
      Index components_size(Index species_index) const;

    private:

      const Supercell *m_unit_scel;
      const Supercell *m_mc_scel;
      std::vector<Molecule> m_struc_mol;
      std::vector<std::string> m_struc_molname;

      Index m_Nasym;
      std::vector<Index> m_unitl_to_asym;
      std::vector<std::set<Index> > m_asym_to_unitl;
      std::vector<std::set<Index> > m_asym_to_b;

      /// m_occ_to_species[asym][occ_index] -> species_index
      std::vector<std::vector<Index> > m_occ_to_species;

      /// m_species_to_occ[asym][species_index] -> occ_index
      std::vector<std::vector<Index> > m_species_to_occ;
    };
  }
}

#endif
