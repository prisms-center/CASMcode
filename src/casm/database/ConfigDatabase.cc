#include "casm/database/ConfigDatabase.hh"
#include <boost/range/iterator_range.hpp>

namespace CASM {
  namespace DB {

    /// Number of Configuration in a particular supecell
    Index Database<Configuration>::scel_range_size(const std::string &scelname) const {
      return boost::distance(scel_range(scelname));
    }

  }
}
