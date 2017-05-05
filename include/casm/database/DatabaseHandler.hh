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
      Database<T> &db();

      /// Access default Database<T>
      template<typename T>
      const Database<T> &db() const;

      /// Access default Database<T>
      template<typename T>
      const Database<T> &const_db();


      /// Access default PropertiesDatabase
      template<typename T>
      PropertiesDatabase &db_props();

      /// Access default PropertiesDatabase
      template<typename T>
      const PropertiesDatabase &db_props() const;

      /// Access default PropertiesDatabase
      template<typename T>
      const PropertiesDatabase &const_db_props();


      /// Close all databases
      void close();


      /// Insert a Database
      template<typename T>
      void insert(std::string db_name, std::unique_ptr<DatabaseBase> &&value);

      /// Insert a PropertiesDatabase
      template<typename T>
      void insert_props(std::string db_name, std::unique_ptr<PropertiesDatabase> &&value);


      // --- Access non-default database (i.e. jsonDB, lmdbDB, mongoDB) ---
      //  - This would be used to migrate from one to the other, so it it not
      //    a common use case

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
      PropertiesDatabase &db_props(std::string db_name);

      /// Access specified PropertiesDatabase
      template<typename T>
      const PropertiesDatabase &db_props(std::string db_name) const;

      /// Access specified PropertiesDatabase
      template<typename T>
      const PropertiesDatabase &const_db_props(std::string db_name);


    private:

      // std::pair<traits<ValueType>::name, db_name> -> db
      typedef std::map <
      std::pair<std::string, std::string>,
          std::unique_ptr<DatabaseBase> > map_type;

      // std::pair<traits<ValueType>::name, db_name> -> db
      typedef std::map <
      std::pair<std::string, std::string>,
          std::unique_ptr<PropertiesDatabase> > props_map_type;

      template<typename T>
      map_type::iterator _find(std::string db_name) const;

      template<typename T>
      props_map_type::iterator _find_props(std::string db_name) const;

      template<typename T>
      void _no_database_error(std::string db_name) const;

      template<typename T>
      void _no_props_database_error(std::string db_name) const;

      const PrimClex *m_primclex;

      std::string m_default_db_name;

      // std::pair<traits<ValueType>::name, db_name> -> db
      // mutable for lazy initialization
      mutable map_type m_db;

      // std::pair<traits<ValueType>::name, db_name> -> db_props
      // mutable for lazy initialization
      mutable props_map_type m_db_props;

    };

  }
}

#endif
