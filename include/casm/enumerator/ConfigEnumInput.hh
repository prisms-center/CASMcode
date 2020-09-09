#ifndef CASM_ConfigEnumInput
#define CASM_ConfigEnumInput

#include <string>

#include "casm/global/definitions.hh"
#include "casm/symmetry/PermuteIterator.hh"
#include "casm/clex/Configuration.hh"

namespace CASM {
  /** \addtogroup Enumerator
      @{
  */

  class ConfigEnumInput {
  public:
    ConfigEnumInput(Configuration const &_config, std::set<Index> const &_sites_selection = {}, std::string const &_name = "");
    ConfigEnumInput(Supercell const &_scel, std::set<Index> const &_sites_selection = {});

    std::string const &name() const {
      return m_name;
    }

    template<typename SitesContainer>
    void set_sites(SitesContainer const &_container) {
      m_sites_selection.clear();
      for(auto const &site : _container)
        _add_site(site);
    }

    std::set<Index> const &sites() const {
      return m_sites_selection;
    }

    Supercell const &supercell() const {
      return m_config.supercell();
    }

    Configuration const &config() const {
      return m_config;
    }

    ConfigDoF const &configdof() const {
      return m_config.configdof();
    }

    void set_group(std::vector<PermuteIterator> const &_group) {
      m_group = _group;
    }

    std::vector<PermuteIterator> const &group() const {
      if(m_group.empty())
        _generate_group();
      return m_group;
    }

  private:

    void _generate_group() const;

    void _add_site(Index b);

    void _add_site(UnitCellCoord const &_ucc);

    std::string m_name;
    std::set<Index> m_sites_selection;
    Configuration m_config;
    mutable std::vector<PermuteIterator> m_group;
  };

  /** @}*/
}

#endif
