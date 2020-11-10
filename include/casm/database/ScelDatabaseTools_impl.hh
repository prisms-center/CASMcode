#ifndef CASM_ScelDatabaseTools_impl
#define CASM_ScelDatabaseTools_impl

#include "casm/database/ScelDatabaseTools.hh"

namespace CASM {

  namespace DB {

    /// Make canonical supercell and insert into supercell database
    ///
    /// - This version checks `is_guaranteed_for_database_insert(enumerator)` and either inserts
    ///   directly or makes canonical and then inserts
    template<typename EnumeratorType>
    std::pair<Database<Supercell>::iterator, bool> make_canonical_and_insert(
      EnumeratorType const &enumerator,
      Supercell const &supercell,
      Database<Supercell> &supercell_db) {
      if(is_guaranteed_for_database_insert(enumerator)) {
        return supercell_db.insert(supercell);
      }
      else {
        return make_canonical_and_insert(supercell.shared_prim(),
                                         supercell.lattice(),
                                         supercell_db);
      }
    }
  }
}

#endif
