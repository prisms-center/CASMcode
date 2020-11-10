#ifndef CASM_enum_enumerate_supercells
#define CASM_enum_enumerate_supercells

#include <functional>
#include <map>
#include <string>

namespace CASM {

  class PrimClex;
  class Supercell;

  namespace DB {
    template<typename T> class Database;
  }

  /// Options for the `enumerate_supercells` function
  struct EnumerateSupercellsOptions {

    EnumerateSupercellsOptions(PrimClex const &primclex):
      primclex_ptr(&primclex) {}

    /// Method name, for printing progress
    std::string method_name;

    /// If filter(supercell)==true, keep supercell, else skip
    std::function<bool (Supercell const &)> filter;

    /// If dry_run==true, do not save results, just print to screen
    bool dry_run = false;

    /// Printing verbosity level
    int verbosity = 10;

    /// Use while transitioning Supercell to no longer need a `PrimClex const *`
    PrimClex const *primclex_ptr = nullptr;
  };

  template<typename EnumeratorType>
  void enumerate_supercells(
    EnumerateSupercellsOptions const &options,
    EnumeratorType &enumerator,
    DB::Database<Supercell> &supercell_db);

}

#endif
