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
