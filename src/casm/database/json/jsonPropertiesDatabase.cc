#include "casm/database/json/jsonPropertiesDatabase.hh"

#include <boost/filesystem.hpp>
#include "casm/casm_io/json_io/container.hh"
#include "casm/casm_io/SafeOfstream.hh"

namespace CASM {
  namespace DB {

    jsonPropertiesDatabase::jsonPropertiesDatabase(const PrimClex &_primclex, std::string calc_type, fs::path location) :
      PropertiesDatabase(_primclex),
      m_is_open(false),
      m_calc_type(calc_type),
      m_location(location) {}

    DatabaseBase &jsonPropertiesDatabase::open() {
      if(m_is_open || !fs::exists(m_location)) {
        return *this;
      }

      jsonParser json(m_location);

      from_json(m_default_score, json["default_conflict_score"]);

      {
        auto it = json["conflict_score"].begin();
        auto end = json["conflict_score"].end();
        for(; it != end; ++it) {
          set_score_method(it.name(), ScoreMappedProperties(*it));
        }
      }

      {
        MappedProperties obj;
        auto it = json["data"].begin();
        auto end = json["data"].end();
        for(; it != end; ++it) {
          from_json(obj, *it);
          insert(obj);
        }
      }

      m_is_open = true;
      return *this;
    }

    void jsonPropertiesDatabase::commit() {
      jsonParser json;

      json["data"] = m_data;
      json["default_conflict_score"] = m_default_score;

      json["conflict_score"].put_obj();
      jsonParser &j = json["conflict_score"];
      for(const auto &val : m_relaxed_from) {
        if(val.second.key_comp().score_method() != m_default_score) {
          j[val.first] = val.second.key_comp().score_method();
        }
      }

      SafeOfstream file;
      fs::create_directories(m_location.parent_path());
      file.open(m_location);
      //json.print(file.ofstream());
      int indent = 0;
      int prec = 12;
      json_spirit::write_stream((json_spirit::mValue &) json, file.ofstream(), indent, prec);
      file.close();
    }

    void jsonPropertiesDatabase::close() {
      m_data.clear();
      m_relaxed_from.clear();
      m_is_open = false;
    }

    /// \brief Begin iterator
    jsonPropertiesDatabase::iterator jsonPropertiesDatabase::begin() const {
      return _iterator(m_data.begin());
    }

    /// \brief End iterator
    jsonPropertiesDatabase::iterator jsonPropertiesDatabase::end() const {
      return _iterator(m_data.end());
    }

    jsonPropertiesDatabase::size_type jsonPropertiesDatabase::size() const {
      return m_data.size();
    }

    /// \brief Return iterator to MappedProperties that is the best mapping to specified config
    ///
    /// - Prefers self-mapped, else best scoring
    jsonPropertiesDatabase::iterator jsonPropertiesDatabase::find_via_to(std::string to_configname) const {
      auto it = m_relaxed_from.find(to_configname);
      if(it == m_relaxed_from.end()) {
        return end();
      }
      // it->second is set of all 'from' -> 'to'
      return find_via_from(*it->second.begin());
    }

    /// \brief Return iterator to MappedProperties that is from the specified config
    jsonPropertiesDatabase::iterator jsonPropertiesDatabase::find_via_from(std::string from_configname) const {
      m_key.from = from_configname;
      return _iterator(m_data.find(m_key));
    }


    /// \brief Names of all configurations that relaxed 'from'->'to'
    std::set<std::string, PropertiesDatabase::Compare>
    jsonPropertiesDatabase::relaxed_from_all(std::string to_configname) const {
      auto it = m_relaxed_from.find(to_configname);
      if(it == m_relaxed_from.end()) {
        return _make_set(to_configname, m_default_score);
      }
      else {
        return it->second;
      }
    }

    /// \brief Change the score method for a single configuration
    void jsonPropertiesDatabase::set_score_method(
      std::string to_configname,
      const ScoreMappedProperties &score) {

      auto it = m_relaxed_from.find(to_configname);
      if(it == m_relaxed_from.end()) {

        // do nothing if default score
        if(score == m_default_score) {
          return;
        }

        auto tmp = _make_set(to_configname, score);
        m_relaxed_from.insert({to_configname, tmp});
      }
      else {
        // if no change, return
        if(it->second.value_comp().score_method() == score) {
          return;
        }

        // construct new set and copy from old set
        auto tmp = _make_set(to_configname, score);
        for(const auto &from : it->second) {
          tmp.insert(from);
        }
        it->second = tmp;
      }
    }

    jsonPropertiesDatabase::iterator
    jsonPropertiesDatabase::_iterator(
      std::set<MappedProperties>::const_iterator _it) const {
      return iterator(jsonPropertiesDatabaseIterator(_it));
    }

    /// \brief Private _insert MappedProperties, without modifying 'relaxed_from'
    std::pair<jsonPropertiesDatabase::iterator, bool>
    jsonPropertiesDatabase::_insert(const MappedProperties &value) {
      auto res = m_data.insert(value);
      return std::make_pair(_iterator(res.first), res.second);
    }

    /// \brief Private _erase MappedProperties, without modifying 'relaxed_from'
    jsonPropertiesDatabase::iterator jsonPropertiesDatabase::_erase(iterator pos) {
      auto base_it = static_cast<jsonPropertiesDatabaseIterator *>(pos.get())->base();
      return _iterator(m_data.erase(base_it));
    }

    /// \brief Names of all configurations that relaxed 'from'->'to'
    void jsonPropertiesDatabase::_set_relaxed_from_all(
      std::string to_configname,
      const std::set<std::string, Compare> &_set) {

      auto it = m_relaxed_from.find(to_configname);
      if(it == m_relaxed_from.end()) {
        if(_set.size()) {
          m_relaxed_from.insert({to_configname, _set});
        }
      }
      else {
        if(!_set.size()) {
          m_relaxed_from.erase(it);
        }
        else {
          it->second = _set;
        }
      }
    }

    std::set<std::string, PropertiesDatabase::Compare>
    jsonPropertiesDatabase::_make_set(
      std::string to_configname,
      const ScoreMappedProperties &score) const {

      return std::set<std::string, Compare>(Compare(this, to_configname, score));
    }

  }
}

