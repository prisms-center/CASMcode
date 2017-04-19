#ifndef CASM_DatabaseHandler
#define CASM_DatabaseHandler

#include <utility>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <stdexcept>

namespace CASM {

  class PrimClex;

  namespace DB {

    class DatabaseBase;
    template<typename T> class Database;

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

      /// Access default Database<T>
      template<typename T>
      Database<T> &db();

      /// Access default Database<T>
      template<typename T>
      const Database<T> &db() const;

      /// Access default Database<T>
      template<typename T>
      const Database<T> &const_db();


      /// Access specified Database<T>
      template<typename T>
      Database<T> &db(std::string db_name);

      /// Access specified Database<T>
      template<typename T>
      const Database<T> &db(std::string db_name) const;

      /// Access specified Database<T>
      template<typename T>
      const Database<T> &const_db(std::string db_name);

      /// Close all databases
      void close();

    private:

      // std::pair<QueryTraits<ValueType>::name, db_name> -> db
      typedef std::map <
      std::pair<std::string, std::string>,
          std::unique_ptr<DatabaseBase> > map_type;

      template<typename T>
      map_type::iterator _find(std::string db_name) const;

      template<typename T>
      void _no_database_error(std::string db_name) const;

      const PrimClex *m_primclex;

      std::string m_default_db_name;

      // std::pair<QueryTraits<ValueType>::name, db_name> -> db
      // mutable for lazy initialization
      mutable map_type m_db;

    };
  }
}

#endif
