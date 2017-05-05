#include "casm/database/PropertiesDatabase.hh"

namespace CASM {
  namespace DB {
    /// \brief Compare mapped properties 'from_A' and 'from_B', preferring self-mapped results
    bool PropertiesDatabase::Compare::operator()(
      const std::string &from_A,
      const std::string &from_B) const {

      if(from_A == from_B) {
        return false;
      }
      if(from_A == m_to) {
        return true;
      }
      if(from_B == m_to) {
        return false;
      }
      return m_score(*m_map->find_via_from(from_A)) < m_score(*m_map->find_via_from(from_B));
    }

    /// \brief Insert data
    std::pair<PropertiesDatabase::iterator, bool> PropertiesDatabase::insert(
      const MappedProperties &value) {

      // insert data
      auto res = _insert(value);
      if(!res.second) {
        return res;
      }

      // insert 'to' -> 'from' link
      auto tset = relaxed_from_all(value.to);
      tset.insert(value.from);
      _set_relaxed_from_all(value.to, tset);

      return res;
    }

    /// \brief Erase the data 'from' from_configname
    PropertiesDatabase::iterator PropertiesDatabase::erase(iterator pos) {

      auto tset = relaxed_from_all(pos->to);
      tset.erase(pos->from);
      _set_relaxed_from_all(pos->to, tset);

      return _erase(pos);
    }

  }
}

