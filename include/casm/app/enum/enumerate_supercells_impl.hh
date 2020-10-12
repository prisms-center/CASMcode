#ifndef CASM_enum_enumerate_supercells_impl
#define CASM_enum_enumerate_supercells_impl

#include "casm/app/enum/EnumInterface.hh"
#include "casm/app/enum/enumerate_supercells.hh"
#include "casm/casm_io/Log.hh"
#include "casm/clex/Supercell_impl.hh"
#include "casm/database/ScelDatabase.hh"
#include "casm/database/ScelDatabaseTools_impl.hh"

namespace CASM {

  /// Enumerate supercells
  ///
  /// Note:
  /// - To avoid unnecessary supercell canonicalization, specialize the following method:
  ///   `template<> bool is_guaranteed_for_database_insert(EnumeratorType const &);`
  ///
  /// \param options See EnumerateSupercellsOptions for method options
  /// \param enumerator A supercell enumerator.
  /// \param supercell_db Will commit any new Supercell if `options.dry_run==false`.
  ///
  /// Note:
  /// - Uses CASM::log() for logging progress
  ///
  template<typename EnumeratorType>
  void enumerate_supercells(
    EnumerateSupercellsOptions const &options,
    EnumeratorType &enumerator,
    DB::Database<Supercell> &supercell_db) {

    Log &log = CASM::log();
    std::pair<DB::Database<Supercell>::iterator, bool> insert_result;
    std::string dry_run_msg = CASM::dry_run_msg(options.dry_run);

    log.set_verbosity(options.verbosity);
    log.begin<Log::standard>(options.method_name);

    for(Supercell const &supercell : enumerator) {

      /// Use while transitioning Supercell to no longer need a `PrimClex const *`
      if(!supercell.has_primclex()) {
        supercell.set_primclex(options.primclex_ptr);
      }

      if(options.filter && !options.filter(supercell)) {
        continue;
      }

      // checks `is_guaranteed_for_database_insert(enumerator)` to see if supercell
      // can be directly inserted, else makes canonical before inserting
      insert_result = make_canonical_and_insert(enumerator, supercell, supercell_db);

      if(insert_result.second) {
        log << dry_run_msg << "  Generated: " << insert_result.first->name() << "\n";
      }
      else {
        log << dry_run_msg << "  Generated: " << insert_result.first->name() << " (already existed)\n";
      }
    }
    log << dry_run_msg << "  DONE." << std::endl << std::endl;

    if(!options.dry_run) {
      log << "Committing database..." << std::endl;
      supercell_db.commit();
      log << "  DONE" << std::endl << std::endl;
    }
    log.end_section();
  }

}

#endif
