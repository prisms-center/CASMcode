#ifndef CASM_ConfigDatabase
#define CASM_ConfigDatabase

#include "casm/clex/Configuration.hh"
#include "casm/database/Database.hh"

namespace boost {

template <typename T>
class iterator_range;

}

namespace CASM {

namespace DB {

/// Derived ConfigDatabase must implement public methods:
/// - iterator update(const Configuration &config)
/// - boost::iterator_range<iterator> scel_range(const name_type& scelname)
/// const
///
/// Database insert methods by convention do not enforce canonical forms or
/// ensure that the primitive form of a Configuration is included also.
/// That logic is included in Configuration::insert, which is
/// the safest way to insert new Configuration in the database. But in cases
/// where it is known that a Configuration is generated in primitive,
/// canonical form, the Database insert methods may be used directly, for
/// instance as in the method insert_unique_canon_configs.
///
template <>
class Database<Configuration> : public ValDatabase<Configuration> {
 public:
  Database(const PrimClex &_primclex) : ValDatabase<Configuration>(_primclex) {}

  virtual ~Database() {}

  /// Update record
  virtual iterator update(const Configuration &config) = 0;

  /// Range of Configuration in a particular supecell
  ///
  /// - Should return range {end(), end()} if no Configuration in specified
  /// Supercell
  /// - Note: boost::iterator_range<iterator>::size is not valid for
  ///   DatabaseIterator.  Use boost::distance instead.
  virtual boost::iterator_range<iterator> scel_range(
      const std::string &scelname) const = 0;

  /// Find canonical Configuration in database by comparing DoF
  virtual iterator search(const Configuration &config) const;

  /// Number of Configuration in a particular supercell
  Index scel_range_size(const std::string &scelname) const;
};

}  // namespace DB
}  // namespace CASM

#endif
