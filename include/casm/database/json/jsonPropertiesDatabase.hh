#ifndef CASM_jsonPropertiesDatabase
#define CASM_jsonPropertiesDatabase

#include <boost/filesystem/path.hpp>
#include "casm/database/PropertiesDatabase.hh"


namespace CASM {
  namespace DB {

    class jsonPropertiesDatabase;

    class jsonPropertiesDatabaseIterator : public PropertiesDatabaseIteratorBase {

    public:

      jsonPropertiesDatabaseIterator() {}

      std::unique_ptr<jsonPropertiesDatabaseIterator> clone() const {
        return std::unique_ptr<jsonPropertiesDatabaseIterator>(this->_clone());
      }

    private:

      typedef typename std::set<MappedProperties>::const_iterator base_iterator;
      friend jsonPropertiesDatabase;

      jsonPropertiesDatabaseIterator(base_iterator _it) :
        m_it(_it) {}

      base_iterator base() const {
        return m_it;
      }

      bool equal(const PropertiesDatabaseIteratorBase &other) const override {
        return m_it == static_cast<const jsonPropertiesDatabaseIterator &>(other).m_it;
      }

      void increment() override {
        ++m_it;
      }

      const MappedProperties &dereference() const override {
        return *m_it;
      }

      long distance_to(const PropertiesDatabaseIteratorBase &other) const override {
        return std::distance(m_it, static_cast<const jsonPropertiesDatabaseIterator &>(other).m_it);
      }

      jsonPropertiesDatabaseIterator *_clone() const override {
        return new jsonPropertiesDatabaseIterator(*this);
      }

      base_iterator m_it;
    };


    class jsonPropertiesDatabase : public PropertiesDatabase {

    public:

      jsonPropertiesDatabase(const PrimClex &_primclex, fs::path location);

      DatabaseBase &open() override;

      void commit() override;

      void close() override;

      /// \brief Begin iterator
      iterator begin() const override;

      /// \brief End iterator
      iterator end() const override;

      /// \brief Return iterator to MappedProperties that is the best mapping to specified config
      ///
      /// - Prefers self-mapped, else best scoring
      iterator find_via_to(std::string to_configname) const override;

      /// \brief Return iterator to MappedProperties that is from the specified config
      iterator find_via_from(std::string from_configname) const override;


      /// \brief Names of all configurations that relaxed 'from'->'to'
      std::set<std::string, Compare> relaxed_from_all(std::string to_configname) const override;

      /// \brief Change the score method for a single configuration
      void set_score_method(std::string to_configname, const ScoreMappedProperties &score) override;


    private:

      iterator _iterator(std::set<MappedProperties>::const_iterator _it) const;

      /// \brief Private _insert MappedProperties, without modifying 'relaxed_from'
      std::pair<iterator, bool> _insert(const MappedProperties &value) override;

      /// \brief Private _erase MappedProperties, without modifying 'relaxed_from'
      iterator _erase(iterator pos) override;

      /// \brief Names of all configurations that relaxed 'from'->'to'
      void _set_relaxed_from_all(
        std::string to_configname,
        const std::set<std::string, Compare> &_set) override;

      std::set<std::string, Compare> _make_set(
        std::string to_configname,
        const ScoreMappedProperties &score) const;


      bool m_is_open;

      fs::path m_location;

      ScoreMappedProperties m_default_score;

      // a key whose 'from' value is modified to find MappedProperties in m_data
      mutable MappedProperties m_key;

      // the MappedProperties container
      std::set<MappedProperties> m_data;

      // to -> {from, from, ...}, used to find best mapping
      std::map<std::string, std::set<std::string, Compare> > m_relaxed_from;

    };

  }
}

#endif
