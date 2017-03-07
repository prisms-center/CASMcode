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

      MakeDiffTransInvariantSubgroup(const Configuration &config_prim): m_config_prim(config_prim) {}

      template<typename PermuteOutputIterator>

      PermuteOutputIterator operator()(const DiffusionTransformation &diff_trans, PermuteIterator begin, PermuteIterator end, PermuteOutputIterator result) {
        ConfigIsEquivalent f(m_config_prim, m_config_prim.crystallography_tol());
        return std::copy_if(begin, end, result, f);
      }
      const Configuration &m_config_prim;
    };
  }

  const std::string DiffTransEnumEquivalents::enumerator_name = "DiffTransEnumEquivalents";

  DiffTransEnumEquivalents::DiffTransEnumEquivalents(
    const DiffusionTransformation &diff_trans,
    PermuteIterator begin,
    PermuteIterator end,
    const Configuration &bg_config_prim) :
    EnumEquivalents<DiffusionTransformation, PermuteIterator>(diff_trans, begin, end, MakeDiffTransInvariantSubgroup(bg_config_prim)) {
  }

}
