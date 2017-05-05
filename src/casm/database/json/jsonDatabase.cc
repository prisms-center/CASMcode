#include "casm/database/json/jsonDatabase.hh"

#include <boost/filesystem.hpp>
#include <boost/range/iterator_range.hpp>
#include "casm/app/DirectoryStructure.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/casm_io/SafeOfstream.hh"
#include "casm/casm_io/json_io/container.hh"
#include "casm/database/DatabaseHandler_impl.hh"
#include "casm/database/Database_impl.hh"
#include "casm/database/DatabaseTypeDefs.hh"

namespace CASM {

  const std::string traits<DB::jsonDB>::name = "jsonDB";

  const std::string traits<DB::jsonDB>::version = "1.0";

  namespace DB {

    namespace {
      struct InsertImpl {
        InsertImpl(DatabaseHandler &_db_handler) : db_handler(_db_handler) {}
        DatabaseHandler &db_handler;

        template<typename T>
        void eval() {
          db_handler.insert<T>(
            traits<jsonDB>::name,
            notstd::make_unique<jsonDatabase<T> >(db_handler.primclex()));
        }
      };
    }


    void jsonDB::insert(DatabaseHandler &db_handler) {
      DB::for_each_type(InsertImpl(db_handler));
    }

    jsonDatabase<Supercell>::jsonDatabase(const PrimClex &_primclex) :
      Database<Supercell>(_primclex),
      m_is_open(false) {}

    jsonDatabase<Supercell> &jsonDatabase<Supercell>::open() {

      if(m_is_open) {
        return *this;
      }

      if(fs::exists(primclex().dir().scel_list())) {
        _read_scel_list();
      }
      else if(fs::exists(primclex().dir().SCEL())) {
        _read_SCEL();
      }

      this->read_aliases();

      m_is_open = true;
      return *this;
    }

    void jsonDatabase<Supercell>::commit() {
      jsonParser json;

      for(const auto &scel : *this) {
        json[scel.name()] = scel.transf_mat();
      }

      SafeOfstream file;
      file.open(primclex().dir().scel_list());
      json.print(file.ofstream());
      file.close();

      this->write_aliases();
    }

    void jsonDatabase<Supercell>::close() {
      m_is_open = false;
      this->clear();
    }

    void jsonDatabase<Supercell>::_read_scel_list() {
      jsonParser json(primclex().dir().scel_list());

      // check json version
      if(!json.contains("version") || json["version"].get<std::string>() != traits<jsonDB>::version) {
        throw std::runtime_error(
          std::string("Error jsonDB version mismatch: found: ") +
          json["version"].get<std::string>() +
          " expected: " +
          traits<jsonDB>::version);
      }

      if(!json.is_array() || !json.contains("supercells")) {
        throw std::runtime_error(
          std::string("Error invalid format: ") + primclex().dir().scel_list().string());
      }

      auto it = json.begin();
      auto end = json.end();
      for(; it != end; ++it) {
        Eigen::Matrix3i mat;
        from_json(mat, *it);
        this->emplace(&primclex(), mat);
      }
    }

    void jsonDatabase<Supercell>::_read_SCEL() {
      // expect a file with format:
      //
      // Supercell Number: 0 Volume: 1
      // Supercell Transformation Matrix:
      //  1 0 0
      //  0 1 0
      //  0 0 1
      //
      // Supercell Number: 1 Volume: 2
      // Supercell Transformation Matrix:
      //  1 0 -1
      //  0 1 0
      //  0 0 2

      Eigen::Matrix3i mat;

      std::string s;
      fs::ifstream stream(primclex().dir().SCEL());
      while(!stream.eof()) {
        std::getline(stream, s);
        if(s[0] == 'S') {
          std::getline(stream, s);
          stream >> mat;

          this->emplace(&primclex(), mat);
        }
      }
    }


    jsonDatabase<Configuration>::jsonDatabase(const PrimClex &_primclex) :
      Database<Configuration>(_primclex),
      m_is_open(false) {}

    jsonDatabase<Configuration> &jsonDatabase<Configuration>::open() {
      if(m_is_open) {
        return *this;
      }

      if(!fs::exists(primclex().dir().config_list())) {
        m_is_open = true;
        return *this;
      }

      jsonParser json(primclex().dir().config_list());

      if(!json.is_obj() || !json.contains("supercells")) {
        throw std::runtime_error(
          std::string("Error invalid format: ") + primclex().dir().config_list().string());
      }

      // check json version
      if(!json.contains("version") || json["version"].get<std::string>() != traits<jsonDB>::version) {
        throw std::runtime_error(
          std::string("Error jsonDB version mismatch: found: ") +
          json["version"].get<std::string>() +
          " expected: " +
          traits<jsonDB>::version);
      }

      // read config list contents
      auto scel_it = json["supercells"].begin();
      auto scel_end = json["supercells"].end();

      for(; scel_it != scel_end; ++scel_it) {

        auto config_it = scel_it->begin();
        auto config_end = scel_it->end();

        const Supercell &scel = *primclex().db_handler().
                                db<Supercell>(traits<jsonDB>::name).find(scel_it.name());

        for(; config_it != config_end; ++config_it) {
          auto result = m_config_list.emplace(scel, config_it.name(), *config_it);
          _on_insert_or_emplace(result);
        }
      }

      // read next config id for each supercell
      from_json(m_config_id, json["config_id"]);

      this->read_aliases();

      m_is_open = true;
      return *this;
    }

    void jsonDatabase<Configuration>::commit() {

      fs::path config_list_path = primclex().dir().config_list();
      if(primclex().db_handler().db<Supercell>(traits<jsonDB>::name).size() == 0) {
        fs::remove(config_list_path);
        return;
      }

      jsonParser json;

      if(fs::exists(config_list_path)) {
        json.read(config_list_path);
      }
      else {
        json.put_obj();
      }

      for(const auto &config : m_config_list) {
        config.to_json(json["supercells"][config.supercell().name()][config.id()]);
      }

      json["config_id"] = m_config_id;

      SafeOfstream file;
      file.open(config_list_path);
      json.print(file.ofstream());
      file.close();

      this->write_aliases();
    }

    void jsonDatabase<Configuration>::close() {
      m_name_to_config.clear();
      m_config_list.clear();
      m_scel_range.clear();
      m_is_open = false;
    }

    jsonDatabase<Configuration>::iterator jsonDatabase<Configuration>::begin() const {
      return _iterator(m_config_list.begin());
    }

    jsonDatabase<Configuration>::iterator jsonDatabase<Configuration>::end() const {
      return _iterator(m_config_list.end());
    }

    jsonDatabase<Configuration>::size_type jsonDatabase<Configuration>::size() const {
      return m_config_list.size();
    }

    std::pair<jsonDatabase<Configuration>::iterator, bool> jsonDatabase<Configuration>::insert(const Configuration &config) {

      auto result = m_config_list.insert(config);

      return _on_insert_or_emplace(result);
    }

    std::pair<jsonDatabase<Configuration>::iterator, bool> jsonDatabase<Configuration>::insert(const Configuration &&config) {

      auto result = m_config_list.insert(std::move(config));

      return _on_insert_or_emplace(result);
    }

    jsonDatabase<Configuration>::iterator jsonDatabase<Configuration>::update(const Configuration &config) {

      ValDatabase<Configuration>::erase(config.name());
      return insert(config).first;
    }

    jsonDatabase<Configuration>::iterator jsonDatabase<Configuration>::erase(iterator pos) {

      // get m_config_list iterator
      auto base_it = static_cast<db_set_iterator *>(pos.get())->base();

      // erase name & alias
      m_name_to_config.erase(base_it->name());

      // update scel_range
      auto _scel_range_it = m_scel_range.find(base_it->supercell().name());
      if(_scel_range_it->second.first == _scel_range_it->second.second) {
        m_scel_range.erase(_scel_range_it);
      }
      else if(_scel_range_it->second.first == base_it) {
        ++(_scel_range_it->second.first);
      }
      else if(_scel_range_it->second.second == base_it) {
        --(_scel_range_it->second.second);
      }

      // erase Configuration
      return _iterator(m_config_list.erase(base_it));
    }

    jsonDatabase<Configuration>::iterator jsonDatabase<Configuration>::find(const std::string &name_or_alias) const {
      auto it = m_name_to_config.find(this->name(name_or_alias));
      if(it == m_name_to_config.end()) {
        return _iterator(m_config_list.end());
      }
      return _iterator(it->second);
    }

    /// Range of Configuration in a particular supecell
    boost::iterator_range<jsonDatabase<Configuration>::iterator>
    jsonDatabase<Configuration>::scel_range(const std::string &scelname) const {
      auto &res = m_scel_range.find(scelname)->second;
      return boost::make_iterator_range(_iterator(res.first), _iterator(std::next(res.second)));
    }

    /// Update m_name_to_config and m_scel_range after performing an insert or emplace
    std::pair<jsonDatabase<Configuration>::iterator, bool>
    jsonDatabase<Configuration>::_on_insert_or_emplace(std::pair<base_iterator, bool> &result) {

      if(result.second) {

        const Configuration &config = *result.first;

        // set the config id, and increment
        auto _config_id_it = m_config_id.find(config.supercell().name());
        if(_config_id_it == m_config_id.end()) {
          _config_id_it = m_config_id.insert(
                            std::make_pair(config.supercell().name(), 0)).first;
        }
        this->set_id(config, _config_id_it->second++);

        // update name -> config
        m_name_to_config.insert(std::make_pair(config.name(), result.first));

        // check if scel_range needs updating
        auto _scel_range_it = m_scel_range.find(config.supercell().name());

        // new supercell
        if(_scel_range_it == m_scel_range.end()) {
          m_scel_range.emplace(
            config.supercell().name(),
            std::make_pair(result.first, result.first));
        }
        // if new 'begin' of scel range
        else if(_scel_range_it->second.first == std::next(result.first)) {
          _scel_range_it->second.first = result.first;
        }
        // if new 'end' of scel range (!= past-the-last config in scel)
        else if(_scel_range_it->second.second == std::prev(result.first)) {
          _scel_range_it->second.second = result.first;
        }
      }

      return std::make_pair(_iterator(result.first), result.second);
    }

  }
}
