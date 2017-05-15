#ifndef CASM_Database_impl
#define CASM_Database_impl

#include "casm/database/Database.hh"
#include "casm/casm_io/SafeOfstream.hh"
#include "casm/casm_io/jsonParser.hh"
#include "casm/casm_io/json_io/container.hh"
#include "casm/app/DirectoryStructure.hh"
#include "casm/clex/PrimClex.hh"


namespace CASM {
  namespace DB {

    /// For setting an alias
    template<typename ValueType>
    std::pair<typename ValDatabase<ValueType>::iterator, bool>
    ValDatabase<ValueType>::set_alias(const std::string &name_or_alias, const std::string &alias) {

      auto it = find(name_or_alias);

      // if name_or_alias not valid
      if(it == end()) {
        // - return (end, false)
        return std::make_pair(it, false);
      }
      // if name_or_alias valid, but alias already used
      else if(m_alias_to_name.find(alias) != m_alias_to_name.end()) {
        // - return (it, false)
        return std::make_pair(it, false);
      }

      // if name_or_alias valid, and alias available
      std::string name;

      // if name_or_alias is name
      if(m_alias_to_name.find(name_or_alias) == m_alias_to_name.end()) {
        name = name_or_alias;

        auto old_alias_it = m_name_to_alias.find(name);
        // if there is an old_alias, erase it
        if(old_alias_it != m_name_to_alias.end()) {
          m_alias_to_name.erase(old_alias_it->second);
        }
      }
      // if name_or_alias is old_alias
      else {
        //std::string old_alias = name_or_alias;
        name = m_alias_to_name[name_or_alias];
        m_alias_to_name.erase(name_or_alias);
      }

      // update data
      m_name_to_alias[name] = alias;
      m_alias_to_name[alias] = name;

      return std::make_pair(it, true);
    }

    /// Get name from name_or_alias
    ///
    /// - Checks if name_or_alias is a known alias
    /// - If known alias, returns associated name
    /// - If not known alias, assumed to be a name and returns name_or_alias
    template<typename ValueType>
    std::string ValDatabase<ValueType>::name(const std::string &name_or_alias) const {
      auto it = m_alias_to_name.find(name_or_alias);
      if(it == m_alias_to_name.end()) {
        return name_or_alias;
      }
      return it->second;
    }

    /// Get alias from name_or_alias
    ///
    /// - Checks if name_or_alias is a known alias
    /// - If yes, returns name_or_alias
    /// - If no, checks if name_or_alias is a name corresponding to a known alias
    ///   - If yes, returns alias
    ///   - If no, (because name_or_alias is an invalid name, or because
    ///     name_or_alias is a valid name without an alias), returns empty string
    template<typename ValueType>
    std::string ValDatabase<ValueType>::alias(const std::string &name_or_alias) const {
      auto it = m_alias_to_name.find(name_or_alias);

      // if name_or_alias is an alias
      if(it != m_alias_to_name.end()) {
        return name_or_alias;
      }

      auto name_it = m_name_to_alias.find(name_or_alias);

      // name has no known alias
      // - could be because name is valid, but no alias is set
      // - or because name is not valid
      if(name_it == m_name_to_alias.end()) {
        return std::string("");
      }

      // name has known alias
      return name_it->second;
    }

    template<typename ValueType>
    void ValDatabase<ValueType>::read_aliases() {
      fs::path p = primclex().dir().template aliases<ValueType>();
      if(fs::exists(p)) {
        jsonParser json(p);
        from_json(m_name_to_alias, json);
        for(const auto &val : m_name_to_alias) {
          m_alias_to_name[val.second] = val.first;
        }
      }
    }

    template<typename ValueType>
    void ValDatabase<ValueType>::write_aliases() {
      fs::path p = primclex().dir().template aliases<ValueType>();
      fs::create_directories(p.parent_path());

      SafeOfstream file;
      file.open(p);
      jsonParser json(m_name_to_alias);
      json.print(file.ofstream());
      file.close();
    }

  }
}

#endif
