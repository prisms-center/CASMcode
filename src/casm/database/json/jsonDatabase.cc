#include "casm/database/json/jsonDatabase.hh"
#include "casm/clex/Supercell.hh"
#include "casm/clex/Configuration.hh"

namespace CASM {

  namespace DB {

    const std::string Traits<jsonDB>::name = "jsonDB";

    void Traits<jsonDB>::insert(DatabaseHandler &db_handler) {
      db_handler.insert<Supercell>(name, jsonScelDatabase(db_handler.primclex()));
      db_handler.insert<Configuration>(name, jsonConfigDatabase(db_handler.primclex()));
    }

    jsonScelDatabase::jsonScelDatabase(const PrimClex &_primclex) :
      DatabaseBase(_primclex),
      m_is_open(false) {}

    jsonScelDatabase &jsonScelDatabase::open() {

      if(m_is_open) {
        return *this;
      }

      if(fs::exists(primclex.dir().scel_list())) {
        _read_scel_list();
      }
      else if(fs::exists(primclex.dir().SCEL())) {
        _read_SCEL();
      }

      m_is_open = true;
      return *this;
    }

    void jsonScelDatabase::commit() {
      jsonParser json;

      for(const auto &scel : *this) {
        json[scel.name()] = scel.transf_mat();
      }

      SafeOfstream file;
      file.open(primclex.dir().scel_list());
      json.print(file.ofstream());
      file.close();
    }

    void jsonScelDatabase::_read_scel_list() {
      jsonParser json(primclex.dir().scel_list());

      if(!json.is_array() || !json.contains("supercells")) {
        throw std::runtime_error(
          std::string("Error invalid format: ") + config_list_path.str());
      }

      auto it = json.begin();
      auto end = json.end();
      for(; it != end; ++it) {
        Eigen::Vector3i mat;
        from_json(mat, *it);
        this->emplace(&primclex(), Lattice(primclex().prim().lattice().lat_column_mat()*mat));
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

      Eigen::Matrix3d mat;

      std::string s;
      while(!stream.eof()) {
        std::getline(stream, s);
        if(s[0] == 'S') {
          std::getline(stream, s);
          stream >> mat;

          this->emplace(&primclex(), Lattice(primclex().prim().lattice().lat_column_mat()*mat));
        }
      }
    }


    /// Database format version, incremented separately from casm --version
    const std::string jsonConfigDatabase::version = "1";

    jsonConfigDatabase::jsonConfigDatabase(const PrimClex &_primclex) :
      DatabaseBase(_primclex),
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
          std::string("Error invalid format: ") + primclex().dir().config_list().str());
      }

      // check json version
      if(!json.contains("version") || json["version"].get<std::string>() != jsonConfigDatabase::version) {
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
          auto result = m_config_list.emplace(scel, *config_it, config_it.name());
          _on_insert_or_emplace(result);
        }
      }

      // read next config id for each supercell
      from_json(m_config_id, json["config_id"])

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
        config.write(json["supercells"][config.supercell().name()][config.id()]);
      }

      json["config_id"] = m_config_id;

      SafeOfstream file;
      file.open(config_list_path);
      json.print(file.ofstream());
      file.close();
    }

    void jsonConfigDatabase::close() {
      m_name_and_alias.clear();
      m_config_list.clear();
      m_scel_range.clear();
      m_is_open = false;
    }

    jsonConfigDatabase::iterator jsonConfigDatabase::begin() {
      return _iterator(m_config_list.begin());
    }

    jsonConfigDatabase::iterator jsonConfigDatabase::end() {
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

    jsonConfigDatabase::iterator update(const Configuration &config) {
      auto it = m_config_list.find(config);
      if(it != m_config_list.end()) {
        *it = config;
      }
      return _iterator(it);
    }

    jsonConfigDatabase::iterator jsonConfigDatabase::erase(iterator pos) {

      // get m_config_list iterator
      auto base_it = static_cast<jsonConfigDatabaseIterator *>(pos.get())->base();

      // erase name & alias
      m_name_and_alias.erase(base_it->name());
      if(!base_it->alias().empty()) {
        m_name_and_alias.erase(base_it->alias());
      }

      // update scel_range
      auto _scel_range_it = m_scel_range.find(base_it->supercell().name());
      if(_scel_range_it->first == _scel_range_it->second) {
        m_scel_range.erase(_scel_range_it);
      }
      else if(_scel_range_it->first == base_it) {
        ++(_scel_range_it->first);
      }
      else if(_scel_range_it->second == base_it) {
        --(_scel_range_it->second);
      }

      // erase Configuration
      return _iterator(m_config_list.erase(base_it));
    }

    jsonConfigDatabase::iterator jsonConfigDatabase::find(const name_type &name_or_alias) {
      auto it = m_name_and_alias.find(name_or_alias);
      if(it == m_name_and_alias.end()) {
        return _iterator(m_config_list.end());
      }
      return _iterator(it);
    }

    /// For setting alias, the new alias must not already exist
    std::pair<jsonConfigDatabase::iterator, bool>
    jsonConfigDatabase::set_alias(const name_type &name_or_alias, const name_type &alias) {

      // check that new alias doesn't already exist
      auto it = m_name_and_alias.find(alias);
      if(it != m_name_and_alias.end()) {
        return std::make_pair(_iterator(it), false);
      }

      // get existing pointer
      base_iterator base_it = *(m_name_and_alias.find(name_or_alias));

      // swap alias
      std::string old_alias = base_it->alias();
      base_it->set_alias(alias);

      // erase old alias
      if(!old_alias.empty()) {
        m_name_or_alias.erase(old_alias);
      }

      // insert new alias
      auto res = m_name_or_alias.insert(std::make_pair(alias, base_it))->set_alias(alias);
      return std::make_pair(_iterator(res.first), res.second);
    }

    /// Range of Configuration in a particular supecell
    boost::iterator_range<jsonConfigDatabase::iterator>
    jsonConfigDatabase::scel_range(const name_type &scelname) const {
      auto res = m_scel_range.find(scelname)->second;
      return boost::make_iterator_range(_iterator(res->first), _iterator(std::next(res->second)));
    }

    /// Update m_name_and_alias and m_scel_range after performing an insert or emplace
    std::pair<jsonConfigDatabase::iterator, bool>
    jsonConfigDatabase::_on_insert_or_emplace(const std::pair<base_iterator, bool> &result) {

      if(result.second) {

        Configuration &config = *result.first;

        // set the config id, and increment
        auto _config_id_it = m_config_id.find(config.supercell().name());
        if(_config_id_it == m_config_id.end()) {
          _config_id_it = m_config_id.insert(
                            std::make_pair(
                              config.supercell().name(),
                              0));
        }
        config.set_id(_config_id_it->second++);

        // update name & alias
        m_name_and_alias.insert(std::make_pair(config.name(), result.first));
        if(!config.alias().empty()) {
          m_name_and_alias.insert(std::make_pair(config.alias(), result.first));
        }

        // check if scel_range needs updating
        auto _scel_range_it = m_scel_range.find(config.supercell().name());
        if(_scel_range_it == m_scel_range.end()) {
          m_scel_range.insert(
            std::make_pair(
              config.supercell().name(),
              std::make_pair(result, result)));
        }
        else if(_scel_range_it->first == std::next(result)) {
          _scel_range_it->first = result;
        }
        else if(_scel_range_it->second == std::prev(result)) {
          _scel_range_it->second = result;
        }
      }

      return std::make_pair(_iterator(result.first), result.second));
    }

  }
}
