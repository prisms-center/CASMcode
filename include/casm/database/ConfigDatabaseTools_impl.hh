#ifndef CASM_ConfigDatabaseTools_impl
#define CASM_ConfigDatabaseTools_impl

#include "casm/database/ConfigDatabaseTools.hh"

namespace CASM {

  namespace DB {

    // Note: Supercell is being transitioned to no longer have a `PrimClex const *`. While this is
    // occuring there are overloads of these functions with and without a `PrimClex const *`.
    // Currently, the `PrimClex const *` is required for Configuration calculated properties, but
    // not for enumeration without calculated properties.

    /// Insert this configuration (in primitive & canonical form) in the database
    template<typename EnumeratorType>
    ConfigInsertResult make_canonical_and_insert(
      EnumeratorType const &enumerator,
      Configuration const &configuration,
      Database<Supercell> &supercell_db,
      Database<Configuration> &configuration_db,
      bool primitive_only) {
      if(is_guaranteed_for_database_insert(enumerator)) {
        configuration_db.insert(configuration);
      }
      else {
        make_canonical_and_insert(
          configuration,
          supercell_db,
          configuration_db,
          primitive_only);
      }
    }
  }
}

#endif
