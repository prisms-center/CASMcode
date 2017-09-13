#ifndef CASM_DatabaseHandler
#define CASM_DatabaseHandler

#include <utility>
#include <map>
#include <memory>
#include <string>

namespace CASM {

  class PrimClex;

  namespace DB {

    class DatabaseBase;
    template<typename T> class Database;
    template<typename T> class ValDatabase;
    class PropertiesDatabase;

    /// \brief Provides access to all databases
    ///
    /// - Holds map of <type name, db name> -> db
    /// - Lazy initialization
    /// - Does not do any checks on the operations of the databases, just
    ///   provides access
    class DatabaseHandler {

    public:

      /// Constructor
      ///
      /// - Uses PrimClex for constructing Supercell, etc.
      /// - Uses PrimClex.settings().db_type() to determine default database type
      DatabaseHandler(const PrimClex &_primclex);

      ~DatabaseHandler();

      const PrimClex &primclex() const;


      /// Access default Database<T>
      template<typename T>
      ValDatabase<T> &generic_db();

      /// Access default Database<T>
      template<typename T>
      const ValDatabase<T> &generic_db() const;

      /// Access default Database<T>
      template<typename T>
      const ValDatabase<T> &const_generic_db();


      /// Access default Database<T>
      template<typename T>
      Database<T> &db();

      /// Access default Database<T>
      template<typename T>
      const Database<T> &db() const;

      /// Access default Database<T>
      template<typename T>
      const Database<T> &const_db();


      /// Access default PropertiesDatabase
      template<typename T>
      PropertiesDatabase &db_props(std::string calc_type);

      /// Access default PropertiesDatabase
      template<typename T>
      const PropertiesDatabase &db_props(std::string calc_type) const;

      /// Access default PropertiesDatabase
      template<typename T>
      const PropertiesDatabase &const_db_props(std::string calc_type);


      /// Close all databases
      void close();


      /// Insert a Database
      template<typename T>
      void insert(std::string db_name, std::unique_ptr<DatabaseBase> &&value);

      /// Insert a PropertiesDatabase
      template<typename T>
      void insert_props(std::string db_name, std::string calc_type, std::unique_ptr<PropertiesDatabase> &&value);


      // --- Access non-default database (i.e. jsonDB, lmdbDB, mongoDB) ---
      //  - This would be used to migrate from one to the other, so it it not
      //    a common use case

      /// Access specified ValDatabase<T>
      template<typename T>
      ValDatabase<T> &generic_db(std::string db_name);

      /// Access specified ValDatabase<T>
      template<typename T>
      const ValDatabase<T> &generic_db(std::string db_name) const;

      /// Access specified ValDatabase<T>
      template<typename T>
      const ValDatabase<T> &const_generic_db(std::string db_name);


      /// Access specified Database<T>
      template<typename T>
      Database<T> &db(std::string db_name);

      /// Access specified Database<T>
      template<typename T>
      const Database<T> &db(std::string db_name) const;

      /// Access specified Database<T>
      template<typename T>
      const Database<T> &const_db(std::string db_name);


      /// Access specified PropertiesDatabase
      template<typename T>
      PropertiesDatabase &db_props(std::string db_name, std::string calc_type);

      /// Access specified PropertiesDatabase
      template<typename T>
      const PropertiesDatabase &db_props(std::string db_name, std::string calc_type) const;

      /// Access specified PropertiesDatabase
      template<typename T>
      const PropertiesDatabase &const_db_props(std::string db_name, std::string calc_type);


    private:

      // std::pair<traits<ValueType>::name, db_name> -> db
      typedef std::map <
      std::pair<std::string, std::string>,
          std::unique_ptr<DatabaseBase> > map_type;

      class PropDBKey {
      public:
        PropDBKey(std::string _config_type, std::string _db_type, std::string _calc_type) :
          m_value(_config_type, _db_type, _calc_type) {}

        std::string &config_type() {
          return std::get<0>(m_value);
        }
        const std::string &config_type() const {
          return std::get<0>(m_value);
        }

        std::string &db_type() {
          return std::get<1>(m_value);
        }
        const std::string &db_type() const {
          return std::get<1>(m_value);
        }

        std::string &calc_type() {
          return std::get<2>(m_value);
        }
        const std::string &calc_type() const {
          return std::get<2>(m_value);
        }

        bool operator<(const PropDBKey &B) const {
          return m_value < B.m_value;
        }

      private:
        std::tuple<std::string, std::string, std::string> m_value;
      };

      // std::pair<traits<ValueType>::name, db_name, calctype> -> db
      typedef std::map <
      PropDBKey,
      std::unique_ptr<PropertiesDatabase> > props_map_type;

      template<typename T>
      map_type::iterator _find(std::string db_name) const;

      template<typename T>
      props_map_type::iterator _find_props(std::string db_name, std::string calc_type) const;

      template<typename T>
      void _no_database_error(std::string db_name) const;

      template<typename T>
      void _no_props_database_error(std::string db_name, std::string calc_type) const;


      const PrimClex *m_primclex;

      std::string m_default_db_name;

      // std::pair<traits<ValueType>::name, db_name> -> db
      // mutable for lazy initialization
      mutable map_type m_db;

      // PropDBKey -> db_props
      // mutable for lazy initialization
      mutable props_map_type m_db_props;

    };

  }
}

#endif
