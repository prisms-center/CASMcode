

namespace CASM {

  namespace DB {

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
      DatabaseHandler(const PrimClex &_primclex) :
        m_primclex(&_primclex),
        m_default_db_type(m_primclex->settings().db_type()) {

        jsonDB::insert(m_db);
      }

      /// Access default Database<T>
      template<typename T>
      Database<T> &db() {
        return db(m_default_db_name);
      }

      /// Access default Database<T>
      template<typename T>
      const Database<T> &db() const {
        return db(m_default_db_name);
      }

      /// Access default Database<T>
      template<typename T>
      const Database<T> &const_db() {
        return const_db(m_default_db_name);
      }


      /// Access specified Database<T>
      template<typename T>
      Database<T> &db(std::string db_name) {
        auto res = _find<T>(std::string db_name);
        return static_cast<Database<T>&>(res->second->open());
      }

      /// Access specified Database<T>
      template<typename T>
      const Database<T> &db(std::string db_name) const {
        auto res = _find<T>(db_name);
        return static_cast<Database<T>&>(res->second->open());
      }

      /// Access specified Database<T>
      template<typename T>
      const Database<T> &const_db(std::string db_name) {
        auto res = _find<T>(db_name);
        return static_cast<Database<T>&>(res->second->open());
      }

      /// Close all databases
      void close() {
        for(auto &db : m_db) {
          db->close();
        }
      }

    private:

      // std::pair<QueryTraits<ValueType>::name, db_name> -> db
      typedef std::map <
      std::make_pair<value_name_type, std::string>,
          std::unique_ptr<DatabaseBase> > map_type;

      template<typename T>
      map_type::iterator _find(std::string db_name) const {
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
      void _no_database_error(std::string db_name) {
        std::stringstream ss;
        ss << "Value: " << QueryTraits<T>::name;
        ss << "  Database: " << db_name;
        throw std::runtime_error("Requested database not found: " + ss.str());
      }

      const PrimClex *m_primclex;

      std::string m_default_db_type;

      // std::pair<QueryTraits<ValueType>::name, db_name> -> db
      // mutable for lazy initialization
      mutable map_type m_db;

    };
  }
}
