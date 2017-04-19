#ifndef CASM_DatabaseHandler
#define CASM_DatabaseHandler

#include "casm/database/DatabaseHandler.hh"

namespace CASM {
  namespace DB {

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
                   QueryTraits<T>::name,
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
      ss << "Value: " << QueryTraits<T>::name;
      ss << "  Database: " << db_name;
      throw std::runtime_error("Requested database not found: " + ss.str());
    }

  }
}

#endif
