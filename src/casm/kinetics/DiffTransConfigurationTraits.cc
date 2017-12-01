#include "casm/kinetics/DiffTransConfigurationTraits.hh"

#include <string>

namespace CASM {

  const std::string traits<Kinetics::DiffTransConfiguration>::name = "DiffTransConfiguration";

  const std::string traits<Kinetics::DiffTransConfiguration>::short_name = "diff_trans_config";

  /// does lexicographical comparison
  bool traits<Kinetics::DiffTransConfiguration>::name_compare(std::string A, std::string B) {
    return A < B;
  };
}
