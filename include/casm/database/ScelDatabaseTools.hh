#ifndef CASM_ScelDatabaseTools
#define CASM_ScelDatabaseTools

#include "casm/database/ScelDatabase.hh"

namespace CASM {

  namespace DB {

    // Note: Supercell is being transitioned to no longer have a `PrimClex const *`. While this is
    // occuring there are overloads of these functions with and without a `PrimClex const *`.
    // Currently, the `PrimClex const *` is required for Configuration calculated properties, but
    // not for enumeration without calculated properties.


    /// Return const reference to canonical equivalent supercell in the supercell database
    Supercell const &canonical_supercell(
      Supercell const &supercell,
      PrimClex const *primclex,
      Database<Supercell> &supercell_db);

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


    /// Return const reference to canonical equivalent supercell in the supercell database
    Supercell const &canonical_supercell(
      Supercell const &supercell,
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
