#ifndef CASM_enum_enumerate_supercells
#define CASM_enum_enumerate_supercells

#include <functional>
#include <map>
#include <string>

#include "casm/casm_io/dataformatter/FormattedDataFile.hh"
#include "casm/clex/Supercell.hh"
#include "casm/database/ScelDatabase.hh"

namespace CASM {

class PrimClex;
class Supercell;

namespace DB {
template <typename T>
class Database;
template <typename T>
class DatabaseIterator;
}  // namespace DB

/// Options for the `enumerate_supercells` function
struct EnumerateSupercellsOptions {
  EnumerateSupercellsOptions(PrimClex const &primclex)
      : primclex_ptr(&primclex) {}

  /// Method name, for printing progress
  std::string method_name;

  /// If filter(supercell)==true, keep supercell, else skip
  std::function<bool(Supercell const &)> filter;

  /// If dry_run==true, do not save results, just print to screen
  bool dry_run = false;

  /// Printing verbosity level
  int verbosity = 10;

  /// Use while transitioning Supercell to no longer need a `PrimClex const *`
  PrimClex const *primclex_ptr = nullptr;

  /// Output enumerated supercells
  bool output_supercells = true;

  /// Options for construcing FormattedDataFile object
  FormattedDataFileOptions output_options;

  /// If true, include output for supercells that were filtered out
  bool output_filtered_supercells = false;
};

/// Collect information during `enumerate_configurations` function for optional
/// output
template <typename EnumeratorType>
struct ScelEnumData {
  ScelEnumData(PrimClex const &_primclex, EnumeratorType const &_enumerator,
               Supercell const &_supercell)
      : primclex(_primclex), enumerator(_enumerator), supercell(_supercell) {}

  PrimClex const &primclex;
  EnumeratorType const &enumerator;
  Supercell const &supercell;

  bool is_excluded_by_filter = false;
  std::pair<DB::Database<Supercell>::iterator, bool> insert_result;
};

template <typename EnumeratorType, typename ScelEnumDataType>
void enumerate_supercells(EnumerateSupercellsOptions const &options,
                          EnumeratorType &enumerator,
                          DB::Database<Supercell> &supercell_db,
                          DataFormatter<ScelEnumDataType> const &formatter);

}  // namespace CASM

#endif
