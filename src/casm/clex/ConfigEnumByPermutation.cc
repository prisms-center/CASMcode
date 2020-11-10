#include "casm/clex/ConfigEnumByPermutation.hh"
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

  const std::string ConfigEnumByPermutation::enumerator_name = "ConfigEnumByPermutation";

  ConfigEnumByPermutation::ConfigEnumByPermutation(
    const Configuration &config) :
    ConfigEnumByPermutation(
      config,
      config.supercell().sym_info().permute_begin(),
      config.supercell().sym_info().permute_end()) {
  }

  ConfigEnumByPermutation::ConfigEnumByPermutation(
    const Configuration &config,
    PermuteIterator begin,
    PermuteIterator end) :
    EnumEquivalents<Configuration, PermuteIterator>(config.canonical_form(), begin, end, MakeConfigInvariantSubgroup()) {
  }

}
