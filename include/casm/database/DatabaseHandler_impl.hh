#ifndef CASM_DatabaseHandler_impl
#define CASM_DatabaseHandler_impl

#include <sstream>
#include "casm/CASM_global_definitions.hh"
#include "casm/database/DatabaseHandler.hh"
#include "casm/database/Database.hh"

namespace CASM {
  namespace DB {

    /// Insert a database
    template<typename T>
    void DatabaseHandler::insert(std::string db_name, std::unique_ptr<DatabaseBase> &&value) {
      m_db.emplace(
        std::make_pair(traits<T>::name, db_name),
        std::move(value));
    }

    /// Access default Database<T>
    template<typename T>
    Database<T> &DatabaseHandler::db() {
      return db<T>(m_default_db_name);
    }

    /// Access default Database<T>
    template<typename T>
    const Database<T> &DatabaseHandler::db() const {
      return db<T>(m_default_db_name);
    }

    /// Access default Database<T>
    template<typename T>
    const Database<T> &DatabaseHandler::const_db() {
      return const_db<T>(m_default_db_name);
    }


    /// Access specified Database<T>
    template<typename T>
    Database<T> &DatabaseHandler::db(std::string db_name) {
      auto res = _find<T>(db_name);
      return static_cast<Database<T>&>(res->second->open());
    }

    /// Access specified Database<T>
    template<typename T>
    const Database<T> &DatabaseHandler::db(std::string db_name) const {
      auto res = _find<T>(db_name);
      return static_cast<Database<T>&>(res->second->open());
    }

    /// Access specified Database<T>
    template<typename T>
    const Database<T> &DatabaseHandler::const_db(std::string db_name) {
      auto res = _find<T>(db_name);
      return static_cast<Database<T>&>(res->second->open());
    }

    /// Close all databases
    void DatabaseHandler::close() {
      for(auto &db : m_db) {
        db.second->close();
      }
    }

    template<typename T>
    DatabaseHandler::map_type::iterator DatabaseHandler::_find(std::string db_name) const {
      auto key = std::make_pair(
                   traits<T>::name,
                   db_name);
      auto res = m_db.find(key);
      if(res == m_db.end()) {
        _no_database_error<T>(db_name);
      }
      return res;
    }

    template<typename T>
    void DatabaseHandler::_no_database_error(std::string db_name) const {
      std::stringstream ss;
      ss << "Value: " << traits<T>::name;
      ss << "  Database: " << db_name;
      throw std::runtime_error("Requested database not found: " + ss.str());
    }

  }
}

#endif
