#include "casm/kinetics/DiffTransEnumEquivalents.hh"
#include <algorithm>
#include "casm/casm_io/Log.hh"
#include "casm/clex/Configuration.hh"
#include "casm/clex/Supercell.hh"
#include "casm/clex/ConfigIsEquivalent.hh"
#include "casm/clex/ScelEnumEquivalents.hh"

namespace CASM {

  namespace {

    struct MakeDiffTransInvariantSubgroup {

      MakeDiffTransInvariantSubgroup() {}

      template<typename PermuteOutputIterator>
        ///  PermuteOutputIterator operator()([](const Configuration &bg_config){return bg_config.primitive().factor_group()}, PermuteIterator begin, PermuteIterator end, PermuteOutputIterator result) {
        /// We probably need a different function for DiffTrans
        /// ConfigIsEquivalent f(config, config.crystallography_tol());

      PermuteOutputIterator operator()(const Configuration &bg_config, PermuteIterator begin, PermuteIterator end, PermuteOutputIterator result) {
        ConfigIsEquivalent f(bg_config.primitive(), bg_config.primitive().crystallography_tol());


        return std::copy_if(begin, end, result, f);
      }

    };
  }

  const std::string DiffTransEnumEquivalents::enumerator_name = "DiffTransEnumEquivalents";

  DiffTransEnumEquivalents::DiffTransEnumEquivalents(
    const DiffusionTransformation &diff_trans,
    PermuteIterator begin,
    PermuteIterator end,
    const Configuration &bg_config) :
    EnumEquivalents<DiffusionTransformation, PermuteIterator, Configuration>(diff_trans.canonical_form(), begin, end, MakeDiffTransInvariantSubgroup()) {
  }

}
