#include "casm/clex/ConfigSelection.hh"

#include "casm/clex/Supercell.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/ConfigIterator.hh"
#include <boost/algorithm/string.hpp>

namespace CASM {

  template<>
  ConfigSelection<false>::ConfigSelection(typename ConfigSelection<false>::PrimClexType &_primclex)
    : m_primclex(&_primclex) {
    m_col_headers.clear();
    for(auto it = _primclex.config_begin(); it != _primclex.config_end(); ++it) {
      m_config[it->name()] = it->selected();
    }
  }

  template<>
  ConfigSelection<true>::ConfigSelection(typename ConfigSelection<true>::PrimClexType &_primclex)
    : m_primclex(&_primclex) {
    m_col_headers.clear();
    for(auto it = _primclex.config_cbegin(); it != _primclex.config_cend(); ++it) {
      m_config[it->name()] = it->selected();
    }
  }


}

