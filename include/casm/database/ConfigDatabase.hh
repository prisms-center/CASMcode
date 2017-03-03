#ifndef CASM_ConfigDatabase
#define CASM_ConfigDatabase

#include "casm/database/Database.hh"
#include "casm/clex/Configuration.hh"

namespace boost {

  template<typename T>
  class iterator_range;

}

namespace CASM {

  namespace DB {

    /// Derived ConfigDatabase must implement public methods:
    /// - iterator update(const Configuration &config)
    /// - boost::iterator_range<iterator> scel_range(const name_type& scelname) const
    ///
    /// Database insert methods by convention do not enforce canonical forms or
    /// ensure that the primitive form of a Configuration is included also.
    /// That logic is included in Configuration::insert, which is
    /// the safest way to insert new Configuration in the database. But in cases
    /// where it is known that a Configuration is generated in primitive,
    /// canonical form, the Database insert methods may be used directly, for
    /// instance as in the method insert_unique_canon_configs.
    ///
    template<>
    class Database<Configuration> : public ValDatabase<Configuration> {

    public:

      /// Update record
      virtual iterator update(const Configuration &config);

      /// Range of Configuration in a particular supecell
      virtual boost::iterator_range<iterator> scel_range(const std::string &scelname) const;

    };

  }
}

#endif
