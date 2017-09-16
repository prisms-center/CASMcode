#include "casm/database/ConfigDatabase.hh"
#include "casm/clex/Supercell.hh"
#include <boost/range/iterator_range.hpp>

namespace CASM {
  namespace DB {

    /// Number of Configuration in a particular supecell
    Index Database<Configuration>::scel_range_size(const std::string &scelname) const {
      return boost::distance(scel_range(scelname));
    }

    /// Find canonical Configuration in database by comparing DoF
    ///
    /// - Default implementation searches scel_range
    typename Database<Configuration>::iterator
    Database<Configuration>::search(const Configuration &config) const {
      auto range = scel_range(config.supercell().name());
      auto res = std::find(range.begin(), range.end(), config);
      if(res == range.end()) {
        return end();
      }
      return res;
    }

  }
}
