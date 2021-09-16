#include "casm/database/PropertiesDatabase.hh"

#include "casm/clex/MappedProperties.hh"

namespace CASM {
namespace DB {
/// \brief Compare mapped properties 'origin_A' and 'origin_B'
bool PropertiesDatabase::Compare::operator()(
    const std::string &origin_A, const std::string &origin_B) const {
  if (origin_A == origin_B) {
    return false;
  }
  double score_A = m_score(*m_map->find_via_origin(origin_A));
  double score_B = m_score(*m_map->find_via_origin(origin_B));
  if (score_A == score_B) {
    return origin_A < origin_B;
  }
  return score_A < score_B;
}

/// \brief Insert data
std::pair<PropertiesDatabase::iterator, bool> PropertiesDatabase::insert(
    const MappedProperties &value) {
  // insert data
  auto res = _insert(value);
  if (!res.second) {
    return res;
  }

  // insert 'to' -> 'origin' link
  auto tset = all_origins(value.to);
  tset.insert(value.origin);
  _set_all_origins(value.to, tset);

  return res;
}

/// \brief Erase the 'origin' data element at provided iterator
PropertiesDatabase::iterator PropertiesDatabase::erase(iterator pos) {
  auto tset = all_origins(pos->to);
  tset.erase(pos->origin);
  _set_all_origins(pos->to, tset);

  return _erase(pos);
}

}  // namespace DB
}  // namespace CASM
