#ifndef CASM_jsonDatabase
#define CASM_jsonDatabase

namespace CASM {

  namespace DB {

    struct jsonDB {};

    template<>
    struct Traits<jsonDB> {

      static const std::string name;

      static void insert(DatabaseHandler &db_handler);
    };

    /// ValueType must have:
    /// - NameType ValueType::name() const
    ///
    /// Derived classes must implement public methods:
    /// - void DatabaseBase& open()
    /// - void commit()
    /// - void close()
    /// - iterator begin()
    /// - iterator end()
    /// - size_type size() const
    /// - std::pair<iterator, bool> insert(const ValueType &obj)
    /// - iterator erase(iterator pos)
    /// - iterator find(const name_type &name) const
    ///
    /// Derived ScelDatabase must implement public methods:
    /// - iterator find(const Lattice &lat)
    /// - std::pair<iterator, bool> insert(const Lattice &lat)
    ///
    class jsonScelDatabase : public ScelDatabase {

    public:

      /// Database format version, incremented separately from casm --version
      static const std::string version;

      jsonScelDatabase(const PrimClex &_primclex);

      jsonScelDatabase &open() override;

      void commit() override;

    private:

      void _read_scel_list();

      void _read_SCEL();

      bool m_is_open;
    };


    /// ValueType must have:
    /// - NameType ValueType::name() const
    ///
    /// Derived classes must implement public methods:
    /// - void DatabaseBase& open()
    /// - void commit()
    /// - void close()
    /// - iterator begin()
    /// - iterator end()
    /// - size_type size() const
    /// - std::pair<iterator, bool> insert(const ValueType &obj)
    /// - std::pair<iterator, bool> insert(const ValueType &&obj)
    /// - iterator erase(iterator pos)
    /// - iterator find(const name_type &name) const
    ///
    /// Derived ConfigDatabase must implement public methods:
    /// - std::pair<iterator, bool> set_alias(const name_type& name_or_alias, const name_type& alias)
    /// - std::pair<iterator, bool> update(const Configuration &config)
    /// - boost::iterator_range<iterator> scel_range(const name_type& scelname) const
    ///
    class jsonConfigDatabase : public ConfigDatabase {

    public:

      /// Database format version, incremented separately from casm --version
      static const std::string version;

      jsonConfigDatabase(const PrimClex &_primclex);

      jsonConfigDatabase &open() override;

      void commit() override;

      void close() override;

      iterator begin() override;

      iterator end() override;

      size_type size() const override ;

      std::pair<iterator, bool> insert(const Configuration &config) override;

      std::pair<iterator, bool> insert(const Configuration &&config) override;

      iterator erase(iterator pos) override;

      iterator find(const name_type &name_or_alias) const override;

      /// For setting alias, the new alias must not already exist
      std::pair<iterator, bool> set_alias(const name_type &name_or_alias, const name_type &alias) override;

      /// Set calc properties
      iterator set_calc_properties(const name_type &name_or_alias, const jsonParser &props) override;

      /// Range of Configuration in a particular supecell
      boost::iterator_range<iterator> scel_range(const name_type &scelname) const override;

    private:

      friend jsonConfigDatabaseIterator;

      typedef std::set<Configuration>::iterator base_iterator;

      /// Update m_name_and_alias and m_scel_range after performing an insert or emplace
      std::pair<iterator, bool> _on_insert_or_emplace(const std::pair<base_iterator, bool> &result) {

        if(result.second) {

          const Configuration &config = *result->second;

          // update name & alias
          m_name_and_alias.insert(std::make_pair(config.name(), result.first));
          if(!config.alias().empty()) {
            m_name_and_alias.insert(std::make_pair(config.alias(), result.first));
          }

          // check if scel_range needs updating
          auto _scel_range_it = m_scel_range.find(config.supercell().name());
          if(_scel_range_it == m_scel_range.end()) {
            m_scel_range.insert(
              std::make_pair(
                config.supercell().name(),
                std::make_pair(result, result)));
          }
          else if(_scel_range_it->first == std::next(result)) {
            _scel_range_it->first = result;
          }
          else if(_scel_range_it->second == std::prev(result)) {
            _scel_range_it->second = result;
          }
        }

        return std::make_pair(_iterator(result.first), result.second));
      }

      // map name and alias -> Configuration
      std::map<std::string, base_iterator> m_name_and_alias;

      // container of Configuration
      std::set<Configuration> m_config_list;

      // map of scelname -> begin and last (not end) configuration in scel
      std::map<std::string, std::pair<base_iterator, base_iterator> > m_scel_range;
    };

  }
}

#endif
