
namespace CASM {

  namespace DB {

    struct jsonDB {};

    template<>
    struct Traits<jsonDB> {
      const std::string name;
      static void insert(DatabaseHandler &db_handler) {
        db_handler.insert<Supercell>(name, jsonScelDatabase(db_handler.primclex()));
        db_handler.insert<Configuration>(name, jsonConfigDatabase(db_handler.primclex()));
      }
    }

    /// Derived classes must implement private methods:
    /// - name_type name() const
    /// - bool is_end() const
    /// - void increment()
    /// - reference dereference() const
    /// - DatabaseIteratorBase *_clone() const
    ///
    class jsonConfigDatabaseIterator : public DatabaseIteratorBase<Configuration> {

    public:

      typedef std::set<Configuration>::iterator base_iterator;

      jsonConfigDatabaseIterator() {}

      jsonConfigDatabaseIterator(const jsonConfigDatabase *_db_ptr,
                                 std::map<name_type, Configuration>::iterator _it) :
        m_db_ptr(_db_ptr),
        m_it(_it) {}

      std::unique_ptr<DatabaseIteratorBase<Configuration> > clone() const {
        return std::unique_ptr<DatabaseIteratorBase<Configuration> >(this->_clone());
      }

      base_iterator base() const {
        return m_it;
      }

    private:

      name_type name() const override {
        return m_it->first;
      }

      bool is_end() const override {
        return !m_db_ptr || m_it == m_dp_ptr->m_config_list.end();
      }

      void increment() override {
        ++m_it;
      }

      reference dereference() const override {
        return m_it->second;
      }

      jsonConfigDatabaseIterator *_clone() const override {
        return new jsonConfigDatabaseIterator(*this);
      }

      const jsonConfigDatabase *m_db_ptr;
      base_iterator m_it;
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
    /// - iterator find(const name_type &name)
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

      jsonConfigDatabase(const PrimClex &_primclex) :
        DatabaseBase(_primclex),
        m_is_open(false) {}

      jsonConfigDatabase &open() override {
        if(m_is_open) {
          return *this;
        }

        if(!fs::exists(primclex.dir().config_list())) {
          m_is_open = true;
          return *this;
        }

        jsonParser json(primclex.dir().config_list());

        if(!json.is_obj() || !json.contains("supercells")) {
          throw std::runtime_error(
            std::string("Error invalid format: ") + config_list_path.str());
        }

        // check json version
        if(!json.contains("version") || json["version"].get<std::string>() != jsonConfigDatabase::version) {
          throw std::runtime_error(
            std::string("Error jsonDB version mismatch: found: ") +
            json["version"].get<std::string>() +
            " expected: " +
            Traits<jsonDB>::version;
        }

        // read config list contents
        auto scel_it = json["supercells"].begin();
        auto scel_end = json["supercells"].end();

        for(; scel_it != scel_end; ++scel_it) {

          auto config_it = scel_it->begin();
          auto config_end = scel_it->end();

          for(; config_it != config_end; ++config_it) {
            auto result = m_config_list.emplace(
                            *config_it,
                            this->primclex(),
                            scel_it.name() + "/" + config_it.name());
            _on_insert_or_emplace(result);
          }
        }

        m_is_open = true;
        return *this;
      }

      void commit() override {

        // ensure commit of the Supercell db
        primclex().db<Supercell>(Traits<jsonDB>::name).commit();

        fs::path config_list_path = primclex().dir().config_list();
        if(primclex().db<Supercell>(Traits<jsonDB>::name).size() == 0) {
          fs::remove(config_list_path);
          return;
        }

        jsonParser json;

        if(fs::exists(config_list_path)) {
          json.read(config_list_path);
        }
        else {
          json.put_obj();
        }

        for(const auto &config : m_config_list) {
          config.write(json["supercells"][config.supercell().name()][config.id()]);
        }

        SafeOfstream file;
        file.open(config_list_path);
        json.print(file.ofstream());
        file.close();
      }

      void close() override {
        m_name_and_alias.clear();
        m_config_list.clear();
        m_scel_range.clear();
        m_is_open = false;
      }

      iterator begin() override {
        return _iterator(m_config_list.begin());
      }

      iterator end() override {
        return _iterator(m_config_list.end());
      }

      size_type size() const override {
        return m_config_list.size();
      }

      std::pair<iterator, bool> insert(const Configuration &config) override {

        auto result = m_config_list.insert(config);

        return _on_insert_or_emplace(result);
      }

      iterator erase(iterator pos) override {

        // get m_config_list iterator
        auto base_it = static_cast<jsonConfigDatabaseIterator *>(pos.get())->base();

        // erase name & alias
        m_name_and_alias.erase(base_it->name());
        if(!base_it->alias().empty()) {
          m_name_and_alias.erase(base_it->alias());
        }

        // update scel_range
        auto _scel_range_it = m_scel_range.find(base_it->supercell().name());
        if(_scel_range_it->first == _scel_range_it->second) {
          m_scel_range.erase(_scel_range_it);
        }
        else if(_scel_range_it->first == base_it) {
          ++(_scel_range_it->first);
        }
        else if(_scel_range_it->second == base_it) {
          --(_scel_range_it->second);
        }

        // erase Configuration
        return _iterator(m_config_list.erase(base_it));
      }

      iterator find(const name_type &name_or_alias) override {
        auto it = m_name_and_alias.find(name_or_alias);
        if(it == m_name_and_alias.end()) {
          return _iterator(m_config_list.end());
        }
        return _iterator(*it);
      }

      /// For setting alias, the new alias must not already exist
      std::pair<iterator, bool> set_alias(const name_type &name_or_alias, const name_type &alias) override {

        // check that new alias doesn't already exist
        auto it = m_name_and_alias.find(alias);
        if(it != m_name_and_alias.end()) {
          return std::make_pair(_iterator(it), false);
        }

        // get existing pointer
        base_iterator base_it = *(m_name_and_alias.find(name_or_alias));

        // swap alias
        std::string old_alias = base_it->alias();
        base_it->set_alias(alias);

        // erase old alias
        if(!old_alias.empty()) {
          m_name_or_alias.erase(old_alias);
        }

        // insert new alias
        auto res = m_name_or_alias.insert(std::make_pair(alias, base_it))->set_alias(alias);
        return std::make_pair(_iterator(res.first), res.second);
      }

      /// For updating properties
      iterator update(const Configuration &config) override {
        auto it = m_config_list.find(config.name());
        it->second = config;
        return _iterator(it);
      }

      /// Range of Configuration in a particular supecell
      boost::iterator_range<iterator> scel_range(const name_type &scelname) const override {
        auto res = m_scel_range.find(scelname)->second;
        return boost::make_iterator_range(_iterator(res->first), _iterator(std::next(res->second)));
      }

    private:

      /// Update m_name_and_alias and m_scel_range after performing an insert or emplace
      std::pair<iterator, bool> _on_insert_or_emplace(base_iterator result) {

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

        return std::make_pair(_iterator(result.first, result.second));
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
