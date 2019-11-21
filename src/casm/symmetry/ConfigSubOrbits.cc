#include "casm/symmetry/ConfigSubOrbits_impl.hh"
#include "casm/clex/Configuration_impl.hh"

namespace CASM {

  MakeConfigSubOrbitGenerators::MakeConfigSubOrbitGenerators(const Configuration &_config) :
    m_config(_config),
    m_prim_config(m_config.primitive().in_canonical_supercell()),
    m_prim_config_fg(m_prim_config.factor_group()),
    m_config_subgroup(
      make_invariant_subgroup(
        m_config.supercell(),
        m_prim_config.supercell(),
        m_prim_config_fg.begin(),
        m_prim_config_fg.end())) {
    std::cout << "In constructor, still " << m_config.supercell().primclex().prim().basis().size() << std::endl;
    std::cout << "In constructor, also " << m_prim_config.supercell().primclex().prim().basis().size() << std::endl;
  }

}
