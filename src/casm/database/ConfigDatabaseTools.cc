#include "casm/database/ConfigDatabaseTools.hh"

#include "casm/clex/Configuration.hh"
#include "casm/clex/FillSupercell.hh"
#include "casm/database/ScelDatabaseTools.hh"

namespace CASM {
namespace DB {

// /// Returns the canonical form Configuration in the canonical Supercell
// ///
// /// Note:
// /// - The canonical Supercell will be inserted in the "supercell_db"
// /// - The canonical Supercell will have a `PrimClex const*` (in the future
// this will be
// ///   deprecated)
// Configuration in_canonical_supercell(
//   Configuration const &configuration,
//   PrimClex const *primclex,
//   Database<Supercell> &supercell_db) {
//
//   Supercell const &canon_supercell_of_configuration =
//   *(make_canonical_and_insert(
//                                                           primclex,
//                                                           configuration.supercell().sym_info().supercell_lattice(),
//                                                           supercell_db).first);
//
//   return fill_supercell(configuration,
//   canon_supercell_of_configuration).canonical_form();
// }
//
// /// Insert this configuration (in primitive & canonical form) in the database
// ///
// /// \param primitive_only If true, only the primitive Configuration is
// inserted.
// ///
// /// Note:
// /// - By convention, the primitive canonical form of a configuration must
// ///   always be saved in the config list.
// /// - Supercells are inserted in the Supercell database as necessary
// /// - If "configuration" is already known to be primitive & canonical, can
// use
// ///   `configuration_db.insert(configuration)` directly.
// /// - The canonical Supercell will have a `PrimClex const*` (this will be
// deprecated) ConfigInsertResult make_canonical_and_insert(
//   Configuration const &configuration,
//   PrimClex const *primclex,
//   Database<Supercell> &supercell_db,
//   Database<Configuration> &configuration_db,
//   bool primitive_only) {
//
//   ConfigInsertResult res;
//
//   Supercell const &canon_supercell_of_configuration = canonical_supercell(
//                                                         configuration.supercell(),
//                                                         primclex,
//                                                         supercell_db);
//
//   Configuration prim_config_in_canon_supercell = in_canonical_supercell(
//                                                    configuration.primitive(),
//                                                    primclex,
//                                                    supercell_db);
//
//   std::tie(res.primitive_it, res.insert_primitive) =
//     configuration_db.insert(prim_config_in_canon_supercell);
//
//   // if the primitive supercell is the same as the equivalent canonical
//   supercell if(canon_supercell_of_configuration ==
//   prim_config_in_canon_supercell.supercell()) {
//     res.insert_canonical = res.insert_primitive;
//     res.canonical_it = res.primitive_it;
//   }
//   else {
//     if(primitive_only) {
//       res.insert_canonical = false;
//     }
//     else {
//       Configuration configuration_in_canon_supercell = fill_supercell(
//                                                          configuration,
//                                                          canon_supercell_of_configuration).canonical_form();
//
//       std::tie(res.canonical_it, res.insert_canonical) =
//         configuration_db.insert(configuration_in_canon_supercell);
//     }
//   }
//   return res;
// }

/// Returns the canonical form Configuration in the canonical Supercell
///
/// Note:
/// - The canonical Supercell will be inserted in the "supercell_db"
/// - The canonical Supercell will not have a `PrimClex const*` (this is
/// preferred if possible)
Configuration in_canonical_supercell(Configuration const &configuration,
                                     Database<Supercell> &supercell_db) {
  Supercell const &canon_supercell_of_configuration =
      *(make_canonical_and_insert(
            configuration.supercell().shared_prim(),
            configuration.supercell().sym_info().supercell_lattice(),
            supercell_db)
            .first);

  return fill_supercell(configuration, canon_supercell_of_configuration)
      .canonical_form();
}

/// Insert this configuration (in primitive & canonical form) in the database
///
/// \param primitive_only If true, only the primitive Configuration is inserted.
///
/// Note:
/// - By convention, the primitive canonical form of a configuration must
///   always be saved in the config list.
/// - Supercells are inserted in the Supercell database as necessary
/// - If "configuration" is already known to be primitive & canonical, can use
///   `configuration_db.insert(configuration)` directly.
/// - The canonical Supercell will not have a `PrimClex const*` (this is
/// preferred if possible)
ConfigInsertResult make_canonical_and_insert(
    Configuration const &configuration, Database<Supercell> &supercell_db,
    Database<Configuration> &configuration_db, bool primitive_only) {
  ConfigInsertResult res;

  Supercell const &canon_supercell_of_configuration =
      canonical_supercell(configuration.supercell(), supercell_db);

  Configuration prim_config_in_canon_supercell =
      in_canonical_supercell(configuration.primitive(), supercell_db);

  std::tie(res.primitive_it, res.insert_primitive) =
      configuration_db.insert(prim_config_in_canon_supercell);

  // if the primitive supercell is the same as the equivalent canonical
  // supercell
  if (canon_supercell_of_configuration ==
      prim_config_in_canon_supercell.supercell()) {
    res.insert_canonical = res.insert_primitive;
    res.canonical_it = res.primitive_it;
  } else {
    if (primitive_only) {
      res.insert_canonical = false;
      res.canonical_it = configuration_db.end();
    } else {
      Configuration configuration_in_canon_supercell =
          fill_supercell(configuration, canon_supercell_of_configuration)
              .canonical_form();

      std::tie(res.canonical_it, res.insert_canonical) =
          configuration_db.insert(configuration_in_canon_supercell);
    }
  }
  return res;
}
}  // namespace DB
}  // namespace CASM
