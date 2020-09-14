#ifndef CASM_ScelDatabaseTools
#define CASM_ScelDatabaseTools

#include "casm/database/ScelDatabase.hh"

namespace CASM {

  namespace DB {

    /// Make canonical supercell and insert into supercell database
    std::pair<Database<Supercell>::iterator, bool> make_canonical_and_insert(
      PrimClex const *primclex,
      Lattice const &super_lattice,
      Database<Supercell> &supercell_db);

    /// Make canonical supercell and insert into supercell database
    std::pair<Database<Supercell>::iterator, bool> make_canonical_and_insert(
      PrimClex const *primclex,
      Eigen::Matrix3l const &transformation_matrix_to_super,
      Database<Supercell> &supercell_db);

    /// Make canonical supercell and insert into supercell database
    std::pair<Database<Supercell>::iterator, bool> make_canonical_and_insert(
      std::shared_ptr<Structure const> const &shared_prim,
      Lattice const &super_lattice,
      Database<Supercell> &supercell_db);

    /// Make canonical supercell and insert into supercell database
    std::pair<Database<Supercell>::iterator, bool> make_canonical_and_insert(
      std::shared_ptr<Structure const> const &shared_prim,
      Eigen::Matrix3l const &transformation_matrix_to_super,
      Database<Supercell> &supercell_db);
  }
}

#endif
