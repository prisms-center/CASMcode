#ifndef CASM_DiffTransConfigDatabase
#define CASM_DiffTransConfigDatabase

#include <utility>

#include "casm/database/Database.hh"
#include "casm/kinetics/DiffTransConfiguration.hh"

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
    class DiffTransConfigDatabase : public Database<DiffTransConfiguration> {

    public:

      /// For updating properties
      virtual std::pair<iterator, bool> update(const DiffTransConfiguration &config);

      /// Range of DiffTransConfiguration for a particular DiffTransOrbit
      virtual boost::iterator_range<iterator> orbit_range(const name_type &diff_trans_orbit_name) const;

      /// Range of DiffTransConfiguration in a particular DiffTransOrbit, in a particular supecell
      virtual boost::iterator_range<iterator> orbit_scel_range(const name_type &diff_trans_orbit_name, const name_type &scelname) const;

    };

  }
}

#endif
