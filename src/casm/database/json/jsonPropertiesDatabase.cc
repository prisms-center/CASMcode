#include "casm/database/jsonPropertiesDatabase.hh"

namespace CASM {

  namespace DB {

    jsonPropertiesDatabase::jsonPropertiesDatabase(const PrimClex &_primclex) :
      PropertiesDatabase(_primclex),
      m_is_open(false) {}

    PropertiesDatabaseBase &jsonPropertiesDatabase::open() override {
      if(m_is_open) {
        return *this;
      }

      jsonParser json(m_file);

      from_json(m_default_score, "default_conflict_score");

      {
        auto it = json["conflict_score"].begin();
        auto end = json["conflict_score"].end();
        for(; it != end; ++it) {
          set_score_method(it.name(), ScoreMappedProperties(*it));
        }
      }

      {
        auto it = json["data"].begin();
        auto end = json["data"].end();
        for(; it != end; ++it) {
          insert(MappedProperties(*it));
        }
      }

      m_is_open = true;
      return *this;
    }

    void jsonPropertiesDatabase::commit() override {
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
      file.open(m_file);
      json.print(file.ofstream());
      file.close();
    }

    void jsonPropertiesDatabase::close() override {
      m_data.clear();
      m_relaxed_from.clear();
      m_is_open = false;
    }

  }

}

