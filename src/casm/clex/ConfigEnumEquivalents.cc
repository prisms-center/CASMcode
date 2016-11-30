#include "casm/clex/ConfigEnumEquivalents.hh"
#include <algorithm>
#include "casm/casm_io/Log.hh"
#include "casm/clex/Configuration.hh"
#include "casm/clex/Supercell.hh"
#include "casm/clex/ConfigIsEquivalent.hh"
#include "casm/clex/ScelEnumEquivalents.hh"

namespace CASM {

  namespace {

    struct MakeConfigInvariantSubgroup {

      MakeConfigInvariantSubgroup() {}

      template<typename PermuteOutputIterator>
      PermuteOutputIterator operator()(const Configuration &config, PermuteIterator begin, PermuteIterator end, PermuteOutputIterator result) {
        ConfigIsEquivalent f(config, config.crystallography_tol());
        return std::copy_if(begin, end, result, f);
      }

    };
  }

  const std::string ConfigEnumEquivalents::enumerator_name = "ConfigEnumEquivalents";

  ConfigEnumEquivalents::ConfigEnumEquivalents(
    const Configuration &config) :
    ConfigEnumEquivalents(
      config,
      config.get_supercell().permute_begin(),
      config.get_supercell().permute_end()) {
  }

  ConfigEnumEquivalents::ConfigEnumEquivalents(
    const Configuration &config,
    PermuteIterator begin,
    PermuteIterator end) :
    EnumEquivalents<Configuration, PermuteIterator>(config.canonical_form(), begin, end, MakeConfigInvariantSubgroup()) {
  }

}
