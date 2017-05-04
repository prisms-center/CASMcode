#include "casm/database/json/jsonDatabase.hh"

#include <boost/filesystem.hpp>
#include <boost/range/iterator_range.hpp>
#include "casm/app/DirectoryStructure.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/casm_io/SafeOfstream.hh"
#include "casm/casm_io/json_io/container.hh"
#include "casm/database/DatabaseHandler.hh"
#include "casm/database/DatabaseTypeDefs.hh"

namespace CASM {
  namespace DB {

    const std::string Traits<jsonDB>::name = "jsonDB";

    const std::string Traits<jsonDB>::version = "1.0";

    void Traits<jsonDB>::insert(DatabaseHandler &db_handler) {
      db_handler.insert<Supercell>(
        name,
        notstd::make_unique<jsonScelDatabase>(db_handler.primclex()));

      db_handler.insert<Configuration>(
        name,
        notstd::make_unique<jsonScelDatabase>(db_handler.primclex()));
    }

    jsonScelDatabase::jsonScelDatabase(const PrimClex &_primclex) :
      Database<Supercell>(_primclex),
      m_is_open(false) {}

    jsonScelDatabase &jsonScelDatabase::open() {

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

    void jsonScelDatabase::commit() {
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

    void jsonScelDatabase::close() {
      m_is_open = false;
      this->clear();
    }

    void jsonScelDatabase::_read_scel_list() {
      jsonParser json(primclex().dir().scel_list());

      // check json version
      if(!json.contains("version") || json["version"].get<std::string>() != Traits<jsonDB>::version) {
        throw std::runtime_error(
          std::string("Error jsonDB version mismatch: found: ") +
          json["version"].get<std::string>() +
          " expected: " +
          Traits<jsonDB>::version);
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

    void jsonScelDatabase::_read_SCEL() {
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


    jsonConfigDatabase::jsonConfigDatabase(const PrimClex &_primclex) :
      Database<Configuration>(_primclex),
      m_is_open(false) {}

    jsonConfigDatabase &jsonConfigDatabase::open() {
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
      if(!json.contains("version") || json["version"].get<std::string>() != Traits<jsonDB>::version) {
        throw std::runtime_error(
          std::string("Error jsonDB version mismatch: found: ") +
          json["version"].get<std::string>() +
          " expected: " +
          Traits<jsonDB>::version);
      }

      // read config list contents
      auto scel_it = json["supercells"].begin();
      auto scel_end = json["supercells"].end();

      for(; scel_it != scel_end; ++scel_it) {

        auto config_it = scel_it->begin();
        auto config_end = scel_it->end();

        const Supercell &scel = *primclex().db<Supercell>().find(scel_it.name());

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

    void jsonConfigDatabase::commit() {

      fs::path config_list_path = primclex().dir().config_list();
      if(primclex().db<Supercell>(Traits<jsonDB>::name).size() == 0) {
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

    void jsonConfigDatabase::close() {
      m_name_to_config.clear();
      m_config_list.clear();
      m_scel_range.clear();
      m_is_open = false;
    }

    jsonConfigDatabase::iterator jsonConfigDatabase::begin() const {
      return _iterator(m_config_list.begin());
    }

    jsonConfigDatabase::iterator jsonConfigDatabase::end() const {
      return _iterator(m_config_list.end());
    }

    jsonConfigDatabase::size_type jsonConfigDatabase::size() const {
      return m_config_list.size();
    }

    std::pair<jsonConfigDatabase::iterator, bool> jsonConfigDatabase::insert(const Configuration &config) {

      auto result = m_config_list.insert(config);

      return _on_insert_or_emplace(result);
    }

    std::pair<jsonConfigDatabase::iterator, bool> jsonConfigDatabase::insert(const Configuration &&config) {

      auto result = m_config_list.insert(std::move(config));

      return _on_insert_or_emplace(result);
    }

    jsonConfigDatabase::iterator jsonConfigDatabase::update(const Configuration &config) {

      ValDatabase<Configuration>::erase(config.name());
      return insert(config).first;
    }

    jsonConfigDatabase::iterator jsonConfigDatabase::erase(iterator pos) {

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

    jsonConfigDatabase::iterator jsonConfigDatabase::find(const std::string &name_or_alias) const {
      auto it = m_name_to_config.find(this->name(name_or_alias));
      if(it == m_name_to_config.end()) {
        return _iterator(m_config_list.end());
      }
      return _iterator(it->second);
    }

    /// Range of Configuration in a particular supecell
    boost::iterator_range<jsonConfigDatabase::iterator>
    jsonConfigDatabase::scel_range(const std::string &scelname) const {
      auto &res = m_scel_range.find(scelname)->second;
      return boost::make_iterator_range(_iterator(res.first), _iterator(std::next(res.second)));
    }

    /// Update m_name_to_config and m_scel_range after performing an insert or emplace
    std::pair<jsonConfigDatabase::iterator, bool>
    jsonConfigDatabase::_on_insert_or_emplace(std::pair<base_iterator, bool> &result) {

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
