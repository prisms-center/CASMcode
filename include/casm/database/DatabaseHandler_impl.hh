#ifndef CASM_DatabaseHandler_impl
#define CASM_DatabaseHandler_impl

#include <sstream>
#include "casm/global/definitions.hh"
#include "casm/database/DatabaseHandler.hh"
#include "casm/database/Database_impl.hh"
#include "casm/database/PropertiesDatabase.hh"

namespace CASM {
  namespace DB {

    /// Insert a Database
    template<typename T>
    void DatabaseHandler::insert(std::string db_name, std::unique_ptr<DatabaseBase> &&value) {
      m_db.emplace(
        std::make_pair(traits<T>::name, db_name),
        std::move(value));
    }

    /// Insert a PropertiesDatabase
    template<typename T>
    void DatabaseHandler::insert_props(std::string db_name, std::string calc_type, std::unique_ptr<PropertiesDatabase> &&value) {
      m_db_props.emplace(
        PropDBKey(traits<T>::name, db_name, calc_type),
        std::move(value));
    }

    /// Access default Database<T>
    template<typename T>
    ValDatabase<T> &DatabaseHandler::generic_db() {
      return generic_db<T>(m_default_db_name);
    }

    /// Access default Database<T>
    template<typename T>
    const ValDatabase<T> &DatabaseHandler::generic_db() const {
      return generic_db<T>(m_default_db_name);
    }

    /// Access default Database<T>
    template<typename T>
    const ValDatabase<T> &DatabaseHandler::const_generic_db() {
      return const_generic_db<T>(m_default_db_name);
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

    /// Access default PropertiesDatabase
    template<typename T>
    PropertiesDatabase &DatabaseHandler::db_props(std::string calc_type)  {
      return db_props<T>(m_default_db_name, calc_type);
    }

    /// Access default PropertiesDatabase
    template<typename T>
    const PropertiesDatabase &DatabaseHandler::db_props(std::string calc_type) const  {
      return db_props<T>(m_default_db_name, calc_type);
    }

    /// Access default PropertiesDatabase
    template<typename T>
    const PropertiesDatabase &DatabaseHandler::const_db_props(std::string calc_type)  {
      return const_db_props<T>(m_default_db_name, calc_type);
    }


    /// Access specified Database<T>
    template<typename T>
    ValDatabase<T> &DatabaseHandler::generic_db(std::string db_name) {
      auto res = _find<T>(db_name);
      return static_cast<ValDatabase<T>&>(res->second->open());
    }

    /// Access specified Database<T>
    template<typename T>
    const ValDatabase<T> &DatabaseHandler::generic_db(std::string db_name) const {
      auto res = _find<T>(db_name);
      return static_cast<ValDatabase<T>&>(res->second->open());
    }

    /// Access specified Database<T>
    template<typename T>
    const ValDatabase<T> &DatabaseHandler::const_generic_db(std::string db_name) {
      auto res = _find<T>(db_name);
      return static_cast<ValDatabase<T>&>(res->second->open());
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

    /// Access specified PropertiesDatabase
    template<typename T>
    PropertiesDatabase &DatabaseHandler::db_props(std::string db_name, std::string calc_type) {
      auto res = _find_props<T>(db_name, calc_type);
      return static_cast<PropertiesDatabase &>(res->second->open());
    }

    /// Access specified PropertiesDatabase
    template<typename T>
    const PropertiesDatabase &DatabaseHandler::db_props(std::string db_name, std::string calc_type) const {
      auto res = _find_props<T>(db_name, calc_type);
      return static_cast<PropertiesDatabase &>(res->second->open());
    }

    /// Access specified PropertiesDatabase
    template<typename T>
    const PropertiesDatabase &DatabaseHandler::const_db_props(std::string db_name, std::string calc_type) {
      auto res = _find_props<T>(db_name, calc_type);
      return static_cast<PropertiesDatabase &>(res->second->open());
    }

    template<typename T>
    DatabaseHandler::map_type::iterator DatabaseHandler::_find(std::string db_name) const {
      auto key = std::make_pair(traits<T>::name, db_name);
      auto res = m_db.find(key);
      if(res == m_db.end()) {
        _no_database_error<T>(db_name);
      }
      return res;
    }

    template<typename T>
    DatabaseHandler::props_map_type::iterator DatabaseHandler::_find_props(std::string db_name, std::string calc_type) const {
      auto key = PropDBKey(traits<T>::name, db_name, calc_type);
      auto res = m_db_props.find(key);
      if(res == m_db_props.end()) {
        _no_props_database_error<T>(db_name, calc_type);
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

    template<typename T>
    void DatabaseHandler::_no_props_database_error(std::string db_name, std::string calc_type) const {
      std::stringstream ss;
      ss << "Value: " << traits<T>::name;
      ss << "  Database: " << db_name;
      ss << "  CalcType: " << calc_type;
      throw std::runtime_error("Requested properties database not found: " + ss.str());
    }

  }
}

#endif
