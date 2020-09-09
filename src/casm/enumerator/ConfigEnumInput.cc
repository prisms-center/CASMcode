#include "casm/clex/Supercell.hh"
#include "casm/enumerator/ConfigEnumInput.hh"

namespace CASM {

  ConfigEnumInput::ConfigEnumInput(Configuration const &_config, std::set<Index> const &_sites_selection, std::string const &_name) :
    m_name(_name),
    m_sites_selection(_sites_selection),
    m_config(_config) {

    if(m_name.empty())
      m_name = _config.name();

    if(m_sites_selection.empty()) {
      for(Index i = 0; i < _config.size(); ++i)
        m_sites_selection.insert(i);
    }

  }

  ConfigEnumInput::ConfigEnumInput(Supercell const &_scel, std::set<Index> const &_sites_selection) :
    ConfigEnumInput(Configuration::zeros(_scel), _sites_selection, _scel.name()) {}

  void ConfigEnumInput::_generate_group() const {
    for(PermuteIterator const &perm_it : config().factor_group()) {
      bool add_it = true;

      for(Index s : sites()) {
        if(sites().count(perm_it.permute_ind(s)) == 0) {
          add_it = false;
          break;
        }
      }
      if(add_it)
        m_group.push_back(perm_it);
    }

  }

  void ConfigEnumInput::_add_site(Index b) {
    Index V = m_config.supercell().volume();
    for(Index i = b * V; i < (b + 1)*V; ++i) {
      m_sites_selection.insert(i);
    }
  }

  void ConfigEnumInput::_add_site(UnitCellCoord const &_ucc) {
    m_sites_selection.insert(this->supercell().sym_info().unitcellcoord_index_converter()(_ucc));
    return;
  }
}
