#include "casm/database/json/jsonDatabase.hh"

#include <boost/filesystem.hpp>
#include <boost/range/iterator_range.hpp>

#include "casm/app/DirectoryStructure.hh"
#include "casm/app/QueryHandler_impl.hh"
#include "casm/casm_io/SafeOfstream.hh"
#include "casm/casm_io/container/json_io.hh"
#include "casm/clex/PrimClex_impl.hh"
#include "casm/clex/io/json/ConfigDoF_json_io.hh"
#include "casm/database/DatabaseHandler_impl.hh"
#include "casm/database/DatabaseTypes_impl.hh"
#include "casm/database/Database_impl.hh"
#include "casm/database/json/jsonPropertiesDatabase.hh"

// for testing:
#include "casm/casm_io/container/stream_io.hh"

namespace CASM {

const std::string traits<DB::jsonDB>::name = "jsonDB";

const std::string traits<DB::jsonDB>::version = "1.0";

namespace DB {

namespace {
struct InsertImpl {
  InsertImpl(DatabaseHandler &_db_handler) : db_handler(_db_handler) {}
  DatabaseHandler &db_handler;

  template <typename T>
  void eval() {
    db_handler.insert<T>(
        traits<jsonDB>::name,
        notstd::make_unique<jsonDatabase<T> >(db_handler.primclex()));
  }
};
}  // namespace

namespace {

struct InsertPropsImpl {
  InsertPropsImpl(DatabaseHandler &_db_handler)
      : db_handler(_db_handler),
        primclex(_db_handler.primclex()),
        dir(primclex.dir()),
        json_dir(dir.root_dir()) {}

  DatabaseHandler &db_handler;
  const PrimClex &primclex;
  const DirectoryStructure &dir;
  jsonDB::DirectoryStructure json_dir;

  template <typename T>
  void eval() {
    for (auto calc_type : dir.all_calctype()) {
      fs::path location = json_dir.props_list<T>(calc_type);
      db_handler.insert_props<T>(traits<jsonDB>::name, calc_type,
                                 notstd::make_unique<jsonPropertiesDatabase>(
                                     primclex, calc_type, location));
    }
  }
};
}  // namespace

void jsonDB::insert(DatabaseHandler &db_handler) {
  DB::for_each_type(InsertImpl(db_handler));
  if (db_handler.primclex().has_dir()) {
    DB::for_each_config_type(InsertPropsImpl(db_handler));
  }
}

jsonDB::DirectoryStructure::DirectoryStructure(const fs::path _root)
    : m_dir(_root) {}

template <typename DataObject>
fs::path jsonDB::DirectoryStructure::obj_list() const {
  return m_dir.casm_dir() / traits<jsonDB>::name /
         (traits<DataObject>::short_name + "_list.json");
}

template <typename DataObject>
fs::path jsonDB::DirectoryStructure::props_list(std::string calctype) const {
  return m_dir.casm_dir() / traits<jsonDB>::name / _calctype(calctype) /
         (traits<DataObject>::short_name + "_props.json");
}

std::string jsonDB::DirectoryStructure::_calctype(std::string calctype) const {
  return std::string("calctype.") + calctype;
}

jsonDatabase<Supercell>::jsonDatabase(const PrimClex &_primclex)
    : Database<Supercell>(_primclex), m_is_open(false) {}

jsonDatabase<Supercell> &jsonDatabase<Supercell>::open() {
  if (m_is_open) {
    return *this;
  }

  if (primclex().has_dir()) {
    jsonDB::DirectoryStructure dir(primclex().dir().root_dir());
    if (fs::exists(dir.obj_list<Supercell>())) {
      _read_scel_list();
    } else if (fs::exists(primclex().dir().SCEL())) {
      _read_SCEL();
    }
  }
  master_selection() = Selection<Supercell>(*this);
  this->read_aliases();

  m_is_open = true;
  return *this;
}

void jsonDatabase<Supercell>::commit() {
  if (!primclex().has_dir()) {
    throw std::runtime_error(
        "Error in jsonDatabase<Supercell>::commit(): CASM project has no root "
        "directory.");
  }

  jsonParser json;
  json["version"] = traits<jsonDB>::version;

  for (const auto &scel : *this) {
    json["supercells"][scel.name()] = scel.transf_mat();
  }

  jsonDB::DirectoryStructure dir(primclex().dir().root_dir());

  SafeOfstream file;
  fs::create_directories(dir.obj_list<Supercell>().parent_path());
  file.open(dir.obj_list<Supercell>());
  json.print(file.ofstream());
  file.close();

  this->write_aliases();
  auto handler = primclex().settings().query_handler<Supercell>();
  handler.set_selected(master_selection());
  master_selection().write(
      handler.dict(), primclex().dir().template master_selection<Supercell>(),
      false, false);
}

void jsonDatabase<Supercell>::close() {
  m_is_open = false;

  this->clear();
}

void jsonDatabase<Supercell>::_read_scel_list() {
  if (!primclex().has_dir()) {
    throw std::runtime_error(
        "Error in jsonDatabase<Supercell>::_read_scel_list(): CASM project has "
        "no root directory.");
  }
  jsonDB::DirectoryStructure dir(primclex().dir().root_dir());
  jsonParser json(dir.obj_list<Supercell>());

  // check json version
  if (!json.contains("version") ||
      json["version"].get<std::string>() != traits<jsonDB>::version) {
    throw std::runtime_error(
        std::string("Error jsonDB version mismatch: found: ") +
        json["version"].get<std::string>() +
        " expected: " + traits<jsonDB>::version);
  }

  if (!json.is_obj() || !json.contains("supercells")) {
    throw std::runtime_error(std::string("Error invalid format: ") +
                             dir.obj_list<Supercell>().string());
  }

  std::cout << "begin reading supercells" << std::endl;
  auto it = json["supercells"].begin();
  auto end = json["supercells"].end();
  for (; it != end; ++it) {
    Eigen::Matrix3l mat;
    from_json(mat, *it);
    std::cout << "mat: \n" << mat << std::endl;
    this->emplace(&primclex(), mat);
    std::cout << "emplaced" << std::endl;
  }
  std::cout << "end reading supercells" << std::endl;
}

void jsonDatabase<Supercell>::_read_SCEL() {
  if (!primclex().has_dir()) {
    throw std::runtime_error(
        "Error in jsonDatabase<Supercell>::_read_SCEL(): CASM project has no "
        "root directory.");
  }
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

  Eigen::Matrix3l mat;

  std::string s;
  fs::ifstream stream(primclex().dir().SCEL());
  while (!stream.eof()) {
    std::getline(stream, s);
    if (s[0] == 'S') {
      std::getline(stream, s);
      stream >> mat;

      this->emplace(&primclex(), mat);
    }
  }
}

jsonDatabase<Configuration>::jsonDatabase(const PrimClex &_primclex)
    : Database<Configuration>(_primclex), m_is_open(false) {}

jsonDatabase<Configuration> &jsonDatabase<Configuration>::open() {
  if (m_is_open) {
    return *this;
  }

  if (!primclex().has_dir()) {
    m_is_open = true;
    master_selection() = Selection<Configuration>(*this);
    return *this;
  }

  jsonDB::DirectoryStructure dir(primclex().dir().root_dir());
  fs::path config_list_path = dir.obj_list<Configuration>();

  if (!fs::exists(config_list_path)) {
    m_is_open = true;
    master_selection() = Selection<Configuration>(*this);
    return *this;
  }

  jsonParser json(config_list_path);

  if (!json.is_obj() || !json.contains("supercells")) {
    throw std::runtime_error(std::string("Error invalid format: ") +
                             config_list_path.string());
  }

  // check json version
  if (!json.contains("version") ||
      json["version"].get<std::string>() != traits<jsonDB>::version) {
    throw std::runtime_error(
        std::string("Error jsonDB version mismatch: found: ") +
        json["version"].get<std::string>() +
        " expected: " + traits<jsonDB>::version);
  }

  // read config list contents
  auto scel_it = json["supercells"].begin();
  auto scel_end = json["supercells"].end();
  bool is_new = false;

  for (; scel_it != scel_end; ++scel_it) {
    auto config_it = scel_it->begin();
    auto config_end = scel_it->end();

    const Supercell &scel = *primclex()
                                 .db_handler()
                                 .db<Supercell>(traits<jsonDB>::name)
                                 .find(scel_it.name());

    for (; config_it != config_end; ++config_it) {
      Configuration configuration{scel};
      from_json(configuration.configdof(), (*config_it)["dof"]);

      auto source_it = config_it->find("source");
      if (source_it != config_it->end()) {
        configuration.set_source(*source_it);
      }
      auto cache_it = config_it->find("cache");
      if (cache_it != config_it->end()) {
        configuration.set_initial_cache(*cache_it);
      }

      this->clear_name(configuration);
      // config_it.name() is the JSON attribute name, which is the config ID
      this->set_id(configuration, config_it.name());

      auto result = m_config_list.emplace(configuration);
      _on_insert_or_emplace(result, is_new);
    }
  }

  // read next config id for each supercell
  from_json(m_config_id, json["config_id"]);
  master_selection() = Selection<Configuration>(*this);
  this->read_aliases();

  m_is_open = true;
  return *this;
}

void jsonDatabase<Configuration>::commit() {
  if (!m_is_open) {
    throw std::runtime_error(
        "Error in jsonDatabase<Configuration>::commit(): Database not open");
  }
  if (!primclex().has_dir()) {
    throw std::runtime_error(
        "Error in jsonDatabase<Configuration>::commit(): CASM project has no "
        "root directory.");
  }

  jsonDB::DirectoryStructure dir(primclex().dir().root_dir());
  fs::path config_list_path = dir.obj_list<Configuration>();
  if (primclex().db_handler().db<Supercell>(traits<jsonDB>::name).size() == 0) {
    fs::remove(config_list_path);
    return;
  }

  jsonParser json;
  if (fs::exists(config_list_path)) {
    json.read(config_list_path);
  } else {
    json.put_obj();
  }
  json["version"] = traits<jsonDB>::version;

  json["supercells"] = jsonParser::object();
  for (const auto &config : m_config_list) {
    jsonParser &configjson =
        json["supercells"][config.supercell().name()][config.id()];
    to_json(config.configdof(), configjson["dof"]);
    to_json(config.source(), configjson["source"]);
    configjson["cache"].put_obj();
    if (config.cache_updated()) {
      to_json(config.cache(), configjson["cache"]);
    }
  }

  json["config_id"] = m_config_id;

  SafeOfstream file;
  fs::create_directories(config_list_path.parent_path());
  file.open(config_list_path);
  int indent = 0;
  int prec = 12;
  json_spirit::write_stream((json_spirit::mValue &)json, file.ofstream(),
                            indent, prec);
  file.close();

  this->write_aliases();
  auto handler = primclex().settings().query_handler<Configuration>();
  handler.set_selected(master_selection());

  bool write_json = false;
  bool only_selected = false;
  master_selection().write(
      handler.dict(),
      primclex().dir().template master_selection<Configuration>(), write_json,
      only_selected);
}

void jsonDatabase<Configuration>::close() {
  m_name_to_config.clear();
  m_config_list.clear();
  m_scel_range.clear();

  m_is_open = false;
}

jsonDatabase<Configuration>::iterator jsonDatabase<Configuration>::begin()
    const {
  return _iterator(m_config_list.begin());
}

jsonDatabase<Configuration>::iterator jsonDatabase<Configuration>::end() const {
  return _iterator(m_config_list.end());
}

jsonDatabase<Configuration>::size_type jsonDatabase<Configuration>::size()
    const {
  return m_config_list.size();
}

std::pair<jsonDatabase<Configuration>::iterator, bool>
jsonDatabase<Configuration>::insert(const Configuration &config) {
  auto result = m_config_list.insert(config);

  return _on_insert_or_emplace(result, true);
}

std::pair<jsonDatabase<Configuration>::iterator, bool>
jsonDatabase<Configuration>::insert(const Configuration &&config) {
  auto result = m_config_list.insert(std::move(config));

  return _on_insert_or_emplace(result, true);
}

jsonDatabase<Configuration>::iterator jsonDatabase<Configuration>::update(
    const Configuration &config) {
  ValDatabase<Configuration>::erase(config.name());
  return insert(config).first;
}

jsonDatabase<Configuration>::iterator jsonDatabase<Configuration>::erase(
    iterator pos) {
  // get m_config_list iterator
  auto base_it = static_cast<db_set_iterator *>(pos.get())->base();

  // erase name & alias
  m_name_to_config.erase(base_it->name());
  master_selection().data().erase(base_it->name());
  // update scel_range
  auto _scel_range_it = m_scel_range.find(base_it->supercell().name());
  if (_scel_range_it->second.first == _scel_range_it->second.second) {
    m_scel_range.erase(_scel_range_it);
  } else if (_scel_range_it->second.first == base_it) {
    ++(_scel_range_it->second.first);
  } else if (_scel_range_it->second.second == base_it) {
    --(_scel_range_it->second.second);
  }

  // erase Configuration
  return _iterator(m_config_list.erase(base_it));
}

jsonDatabase<Configuration>::iterator jsonDatabase<Configuration>::find(
    const std::string &name_or_alias) const {
  auto it = m_name_to_config.find(this->name(name_or_alias));
  if (it == m_name_to_config.end()) {
    return _iterator(m_config_list.end());
  }
  return _iterator(it->second);
}

/// Range of Configuration in a particular supecell
boost::iterator_range<jsonDatabase<Configuration>::iterator>
jsonDatabase<Configuration>::scel_range(const std::string &scelname) const {
  auto it = m_scel_range.find(scelname);
  if (it == m_scel_range.end()) {
    return boost::make_iterator_range(end(), end());
  } else {
    auto &res = it->second;
    return boost::make_iterator_range(_iterator(res.first),
                                      _iterator(std::next(res.second)));
  }
}

/// Find canonical Configuration in database by comparing DoF
///
/// \param config A Configuration in canonical form
///
/// - Find in set<Configuration>
typename jsonDatabase<Configuration>::iterator
jsonDatabase<Configuration>::search(const Configuration &config) const {
  // not clear if using m_scel_range to search on a sub-range would help...
  auto res = m_config_list.find(config);
  if (res == m_config_list.end()) {
    return end();
  }
  return _iterator(res);
}

/// Update m_name_to_config and m_scel_range after performing an insert or
/// emplace
std::pair<jsonDatabase<Configuration>::iterator, bool>
jsonDatabase<Configuration>::_on_insert_or_emplace(
    std::pair<base_iterator, bool> &result, bool is_new) {
  if (result.second) {
    const Configuration &config = *result.first;
    assert(&config.primclex() == &primclex() &&
           "jsonDatabase<Configuration>::_on_insert_or_emplace primclex does "
           "not match");

    if (is_new) {
      // set the config id, and increment
      auto _config_id_it = m_config_id.find(config.supercell().name());
      if (_config_id_it == m_config_id.end()) {
        _config_id_it =
            m_config_id.insert(std::make_pair(config.supercell().name(), 0))
                .first;
      }
      this->set_id(config, _config_id_it->second++);
    }

    // update name -> config
    m_name_to_config.insert(std::make_pair(config.name(), result.first));

    // check if scel_range needs updating
    auto _scel_range_it = m_scel_range.find(config.supercell().name());

    // new supercell
    if (_scel_range_it == m_scel_range.end()) {
      m_scel_range.emplace(config.supercell().name(),
                           std::make_pair(result.first, result.first));
    }
    // if new 'begin' of scel range
    else if (_scel_range_it->second.first == std::next(result.first)) {
      _scel_range_it->second.first = result.first;
    }
    // if new 'end' of scel range (!= past-the-last config in scel)
    else if (_scel_range_it->second.second == std::prev(result.first)) {
      _scel_range_it->second.second = result.first;
    }

    master_selection().data().emplace(config.name(), 0);
  }

  return std::make_pair(_iterator(result.first), result.second);
}
}  // namespace DB
}  // namespace CASM

// explicit template instantiations
#define INST_jsonDB(r, data, type) \
  template fs::path jsonDB::DirectoryStructure::obj_list<type>() const;

#define INST_jsonDB_config(r, data, type)                         \
  template fs::path jsonDB::DirectoryStructure::props_list<type>( \
      std::string calctype) const;

namespace CASM {
namespace DB {

BOOST_PP_SEQ_FOR_EACH(INST_jsonDB, _, CASM_DB_TYPES)
BOOST_PP_SEQ_FOR_EACH(INST_jsonDB_config, _, CASM_DB_CONFIG_TYPES)
}  // namespace DB
}  // namespace CASM
