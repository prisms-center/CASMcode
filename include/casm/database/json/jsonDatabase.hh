#ifndef CASM_jsonDatabase
#define CASM_jsonDatabase

#include "casm/database/Database.hh"
#include "casm/database/ScelDatabase.hh"
#include "casm/database/ConfigDatabase.hh"

namespace CASM {

  namespace DB {

    struct jsonDB {};

    template<>
    struct Traits<jsonDB> {

      static const std::string name;

      static void insert(DatabaseHandler &db_handler);
    };

    /// ValueType must have:
    /// - std::string ValueType::name() const
    ///
    /// Derived classes must implement public methods:
    /// - void DatabaseBase& open()
    /// - void commit()
    /// - void close()
    ///
    class jsonScelDatabase : public Database<Supercell> {

    public:

      /// Database format version, incremented separately from casm --version
      static const std::string version;

      jsonScelDatabase(const PrimClex &_primclex);

      jsonScelDatabase &open() override;

      void commit() override;

      void close() override;


    private:

      void _read_scel_list();

      void _read_SCEL();

      bool m_is_open;
    };


    /// ValueType must have:
    /// - std::string ValueType::name() const
    ///
    /// Derived classes must implement public methods:
    /// - void DatabaseBase& open()
    /// - void commit()
    /// - void close()
    /// - iterator begin() const
    /// - iterator end() const
    /// - size_type size() const
    /// - std::pair<iterator, bool> insert(const ValueType &obj)
    /// - std::pair<iterator, bool> insert(const ValueType &&obj)
    /// - iterator erase(iterator pos)
    /// - iterator find(const std::string &name) const
    ///
    /// Derived ConfigDatabase must implement public methods:
    /// - iterator update(const Configuration &config)
    /// - boost::iterator_range<iterator> scel_range(const std::string& scelname) const
    ///
    class jsonConfigDatabase : public Database<Configuration> {

    public:

      /// Database format version, incremented separately from casm --version
      static const std::string version;

      jsonConfigDatabase(const PrimClex &_primclex);

      jsonConfigDatabase &open() override;

      void commit() override;

      void close() override;

      iterator begin() const override;

      iterator end() const override;

      size_type size() const override ;

      std::pair<iterator, bool> insert(const Configuration &config) override;

      std::pair<iterator, bool> insert(const Configuration &&config) override;

      iterator update(const Configuration &config) override;

      iterator erase(iterator pos) override;

      iterator find(const std::string &name_or_alias) const override;

      iterator find(const Configuration &obj) const override;

      /// Range of Configuration in a particular supecell
      boost::iterator_range<iterator> scel_range(const std::string &scelname) const override;

    private:

      friend jsonConfigDatabaseIterator;

      typedef std::set<Configuration>::iterator base_iterator;

      /// Update m_name_and_alias and m_scel_range after performing an insert or emplace
      std::pair<iterator, bool> _on_insert_or_emplace(const std::pair<base_iterator, bool> &result);

      iterator _iterator(base_iterator base_it) const {
        return iterator(DatabaseSetIterator<Configuration, jsonConfigDatabase>(base_it));
      }

      // map name -> Configuration
      std::map<std::string, base_iterator> m_name_to_config;

      // container of Configuration
      std::set<Configuration> m_config_list;

      // map of scelname -> begin and last (not end) configuration in scel
      std::map<std::string, std::pair<base_iterator, base_iterator> > m_scel_range;

      // map of scelname -> next id to assign to a new Configuration
      std::map<std::string, Index> m_config_id;
    };

  }
}

#endif
