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
      //todo
    }

    void jsonScelDatabase::_read_scel_list() {
      jsonParser json(primclex.dir().scel_list());

      if(!json.is_array() || !json.contains("supercells")) {
        throw std::runtime_error(
          std::string("Error invalid format: ") + config_list_path.str());
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

      if(!fs::exists(primclex.dir().config_list())) {
        m_is_open = true;
        return *this;
      }

      jsonParser json(primclex.dir().config_list());

      if(!json.is_obj() || !json.contains("supercells")) {
        throw std::runtime_error(
          std::string("Error invalid format: ") + primclex.dir().config_list().str());
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

        for(; config_it != config_end; ++config_it) {
          auto result = m_config_list.emplace(
                          *config_it,
                          this->primclex(),
                          scel_it.name() + "/" + config_it.name());
          _on_insert_or_emplace(result);
        }
      }

      m_is_open = true;
      return *this;
    }

    void jsonConfigDatabase::commit() {

      // ensure commit of the Supercell db
      primclex().db<Supercell>(Traits<jsonDB>::name).commit();

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

    iterator jsonConfigDatabase::begin() {
      return _iterator(m_config_list.begin());
    }

    iterator jsonConfigDatabase::end() {
      return _iterator(m_config_list.end());
    }

    size_type jsonConfigDatabase::size() const {
      return m_config_list.size();
    }

    std::pair<iterator, bool> jsonConfigDatabase::insert(const Configuration &config) {

      auto result = m_config_list.insert(config);

      return _on_insert_or_emplace(result);
    }

    std::pair<iterator, bool> jsonConfigDatabase::insert(const Configuration &&config) {

      auto result = m_config_list.insert(std::move(config));

      return _on_insert_or_emplace(result);
    }

    iterator jsonConfigDatabase::erase(iterator pos) {

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

    iterator jsonConfigDatabase::find(const name_type &name_or_alias) {
      auto it = m_name_and_alias.find(name_or_alias);
      if(it == m_name_and_alias.end()) {
        return _iterator(m_config_list.end());
      }
      return _iterator(it);
    }

    /// For setting alias, the new alias must not already exist
    std::pair<iterator, bool> jsonConfigDatabase::set_alias(const name_type &name_or_alias, const name_type &alias) {

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

    /// Set calc properties
    iterator jsonConfigDatabase::set_calc_properties(const name_type &name_or_alias, const jsonParser &props) {
      auto it = m_config_list.find(config.name());
      it->set_calc_properties(props);
      return _iterator(it);
    }

    /// Range of Configuration in a particular supecell
    boost::iterator_range<iterator> jsonConfigDatabase::scel_range(const name_type &scelname) const {
      auto res = m_scel_range.find(scelname)->second;
      return boost::make_iterator_range(_iterator(res->first), _iterator(std::next(res->second)));
    }

  }
}
