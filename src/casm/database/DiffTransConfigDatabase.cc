#include "casm/database/DiffTransConfigDatabase.hh"
#include <boost/range/iterator_range.hpp>

namespace CASM {
  namespace DB {

    /// Number of DiffTransConfiguration in a particular supercell
    Index Database<Kinetics::DiffTransConfiguration>::scel_range_size(const std::string &scelname) const {
      return boost::distance(scel_range(scelname));
    }

    /// Number of DiffTransConfiguration in a particular orbit
    Index Database<Kinetics::DiffTransConfiguration>::orbit_range_size(const std::string &diff_trans_name) const {
      return boost::distance(orbit_range(diff_trans_name));
    }

    /// Number of DiffTransConfiguration in a particular supercell within an orbit
    Index Database<Kinetics::DiffTransConfiguration>::orbit_scel_range_size(const std::string &diff_trans_name, const std::string &scelname) const {
      return boost::distance(orbit_scel_range(diff_trans_name, scelname));
    }

  }
}
