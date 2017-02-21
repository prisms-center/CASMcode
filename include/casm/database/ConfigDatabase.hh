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
    /// - std::pair<iterator, bool> rename(const name_type& old_name, const name_type& new_name)
    /// - std::pair<iterator, bool> update(const Configuration &config)
    /// - boost::iterator_range<iterator> scel_range(const name_type& scelname) const
    ///
    template<>
    class Database<Configuration> : public ValDatabase<Configuration, std::string> {

    public:

      /// Set calc properties
      virtual iterator set_calc_properties(const Configuration &config);

      /// Range of Configuration in a particular supecell
      virtual boost::iterator_range<iterator> scel_range(const name_type &scelname) const;

    };

  }
}

#endif
