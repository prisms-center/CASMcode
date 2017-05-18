#ifndef CASM_jsonDatabase
#define CASM_jsonDatabase

#include "casm/app/DirectoryStructure.hh"
#include "casm/database/Database.hh"
#include "casm/database/ScelDatabase.hh"
#include "casm/database/ConfigDatabase.hh"
#include "casm/database/DiffTransOrbitDatabase.hh"

namespace CASM {

  template<typename T> struct traits;

  namespace DB {
    template<typename DataObject> class jsonDatabase;
    class DatabaseHandler;

    struct jsonDB;
  }

  template<>
  struct traits<DB::jsonDB> {

    static const std::string name;

    /// Database format version, incremented separately from casm --version
    static const std::string version;

    static void insert(DB::DatabaseHandler &db_handler);
  };

  namespace DB {

    struct jsonDB {
      static void insert(DatabaseHandler &);

      class DirectoryStructure {
      public:
        DirectoryStructure(const fs::path _root);

        template<typename DataObject>
        fs::path obj_list() const;

        template<typename DataObject>
        fs::path props_list(std::string calctype) const;

      private:
        std::string _calctype(std::string calctype) const;

        CASM::DirectoryStructure m_dir;
      };
    };


    /// ValueType must have:
    /// - std::string ValueType::name() const
    ///
    /// Derived classes must implement public methods:
    /// - void DatabaseBase& open()
    /// - void commit()
    /// - void close()
    ///
    /// This JSON Database has the following structure
    /// json["version"] contains the version (format) of this database.
    /// json["supercells"] is a JSON object that corresponds to a map in which the
    /// supercell name is the key and the value is the information of that supercell.
    template<>
    class jsonDatabase<Supercell> : public Database<Supercell> {

    public:

      jsonDatabase<Supercell>(const PrimClex &_primclex);

      jsonDatabase<Supercell> &open() override;

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
    /// This JSON Database has the following structure
    /// json["version"] contains the version (format) of this database.
    /// json["supercells"] is a JSON object that corresponds to a map in which the
    /// supercell name is the key and the value is the object that contains all
    //  the information of all the configurations that correspond to the supercell.
    /// Each configuration is indexed within this object.
    /// json["config_id"] is a map of supercell name (key) to the next index (value)to be assigned
    /// to the newly enumerated configuration within the given supercell
    template<>
    class jsonDatabase<Configuration> : public Database<Configuration> {

    public:

      jsonDatabase<Configuration>(const PrimClex &_primclex);


      jsonDatabase<Configuration> &open() override;

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

      /// Range of Configuration in a particular supecell
      ///
      /// - Should return range {end(), end()} if no Configuration in specified Supercell
      /// - Note: boost::iterator_range<iterator>::size is not valid for
      ///   DatabaseIterator.  Use boost::distance instead.
      boost::iterator_range<iterator> scel_range(const std::string &scelname) const override;

    private:

      typedef std::set<Configuration>::iterator base_iterator;
      typedef DatabaseSetIterator<Configuration, jsonDatabase<Configuration> > db_set_iterator;

      /// Update m_name_and_alias and m_scel_range after performing an insert or emplace
      std::pair<iterator, bool> _on_insert_or_emplace(std::pair<base_iterator, bool> &result, bool is_new);

      iterator _iterator(base_iterator name_it) const {
        return iterator(db_set_iterator(name_it));
      }

      bool m_is_open;

      // map name -> Configuration
      std::map<std::string, base_iterator> m_name_to_config;

      // container of Configuration
      std::set<Configuration> m_config_list;

      // map of scelname -> begin and last (not end) configuration in scel
      std::map<std::string, std::pair<base_iterator, base_iterator> > m_scel_range;

      // map of scelname -> next id to assign to a new Configuration
      std::map<std::string, Index> m_config_id;
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
    /// Derived DiffTransOrbitDatabase must implement public methods:
    /// nothing for now
    ///
    /// This JSON Database has three sub objects
    /// json["version"] contains the version (format) of this database.
    /// json["protoypes"] is a JSON object that corresponds to a map in which the
    /// orbit id is the key and the value is the information of the prototype of the orbit
    /// json["orbit_id"] is the next id to be assigned to a newly enumerated orbit.
    ///
    template<>
    class jsonDatabase<PrimPeriodicDiffTransOrbit> : public Database<PrimPeriodicDiffTransOrbit> {

    public:

      jsonDatabase<PrimPeriodicDiffTransOrbit>(const PrimClex &_primclex);


      jsonDatabase<PrimPeriodicDiffTransOrbit> &open() override;

      void commit() override;

      void close() override;

      iterator begin() const override;

      iterator end() const override;

      size_type size() const override ;

      std::pair<iterator, bool> insert(const PrimPeriodicDiffTransOrbit &orbit) override;

      std::pair<iterator, bool> insert(const PrimPeriodicDiffTransOrbit &&orbit) override;

      template<typename... Args>
      std::pair<iterator, bool> emplace(Args &&... args) {
        return _on_insert_or_emplace(m_orbit_list.emplace(std::forward<Args>(args)...));
      }

      iterator erase(iterator pos) override;

      iterator find(const std::string &name_or_alias) const override;

    private:

      typedef std::set<PrimPeriodicDiffTransOrbit>::iterator base_iterator;
      typedef DatabaseSetIterator<PrimPeriodicDiffTransOrbit, jsonDatabase<PrimPeriodicDiffTransOrbit> > db_set_iterator;

      /// Update m_name_and_alias and m_scel_range after performing an insert or emplace
      std::pair<iterator, bool> _on_insert_or_emplace(std::pair<base_iterator, bool> &result, bool is_new);

      iterator _iterator(base_iterator name_it) const {
        return iterator(db_set_iterator(name_it));
      }

      bool m_is_open;

      // map name -> PrimPeriodicDiffTransOrbit
      std::map<std::string, base_iterator> m_name_to_orbit;

      // container of PrimPeriodicDiffTransOrbit
      std::set<PrimPeriodicDiffTransOrbit> m_orbit_list;

      //next id to assign to a new PrimPeriodicDiffTransOrbit
      Index m_orbit_id;
    };

  }
}

#endif
