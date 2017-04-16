#ifndef CASM_jsonPropertiesDatabase
#define CASM_jsonPropertiesDatabase

#include "casm/database/PropertiesDatabase.hh"

namespace CASM {

  namespace DB {

    class jsonPropertiesDatabaseIterator : public PropertiesDatabaseIteratorBase {

    public:

      jsonPropertiesDatabaseIterator() {}

      std::unique_ptr<jsonPropertiesDatabaseIterator> clone() const {
        return std::unique_ptr<jsonPropertiesDatabaseIterator>(this->_clone());
      }

    private:

      typedef typename std::set<MappedProperties>::iterator base_iterator;

      DatabaseSetIterator(base_iterator _it) :
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

      jsonPropertiesDatabase(const PrimClex &_primclex);

      PropertiesDatabaseBase &open() override;

      void commit() override;

      void close() override;

      /// \brief Begin iterator
      iterator begin() const override {
        return _iterator(m_data.begin());
      }

      /// \brief End iterator
      iterator end() const override {
        return _iterator(m_data.end());
      }

      /// \brief Return iterator to MappedProperties that is the best mapping to specified config
      ///
      /// - Prefers self-mapped, else best scoring
      iterator find_via_to(std::string to_configname) const override {
        auto it = m_relaxed_from.find(to_configname);
        if(it == m_relaxed_from.end()) {
          return end();
        }
        // it->second is set of all 'from' -> 'to'
        return find_from(*it->second.begin());
      }

      /// \brief Return iterator to MappedProperties that is from the specified config
      iterator find_via_from(std::string from_configname) const override {
        m_key.from = from_configname;
        return m_data.find(m_key);
      }


      /// \brief Names of all configurations that relaxed 'from'->'to'
      std::set<std::string, Compare> relaxed_from_all(std::string to_configname) const {
        return m_relaxed_from.find(to_configname)->second;
      }

      /// \brief Change the score method for a single configuration
      void set_score_method(std::string to_configname, const ScoreMappedProperties &score) {

        auto it = m_relaxed_from.find(to_configname);
        if(it == m_relaxed_from.end()) {

          // do nothing if default score
          if(score == m_default_score) {
            return;
          }

          auto tmp = _make_set(to_configname, score);
          m_relaxed.insert({to_configname, tmp});
        }
        else {
          // if no change, return
          if(it->second.value_comp().score_method() == score) {
            return;
          }

          // construct new set and copy from old set
          auto tmp = _make_set(to_configname, score);
          for(const auto &from : it->second) {
            tmp.insert(from);
          }
          it->second = tmp;
        }
      }


    private:

      iterator _iterator(std::set<MappedProperties>::const_iterator _it) const {
        return iterator(_it);
      }

      /// \brief Private _insert MappedProperties, without modifying 'relaxed_from'
      std::pair<iterator, bool> _insert(const MappedProperties &value) override {
        auto res = m_data.insert(value);
        return std::make_pair(_iterator(res.first), res.second);
      }

      /// \brief Private _erase MappedProperties, without modifying 'relaxed_from'
      iterator _erase(iterator pos) override {
        auto base_it = static_cast<jsonPropertiesDatabaseIterator *>(pos.get())->base();
        return _iterator(m_data.erase(base_it));
      }

      /// \brief Names of all configurations that relaxed 'from'->'to'
      void _set_relaxed_from_all(
        std::string to_configname,
        const std::set<std::string, Compare> &_set) const override {

        auto it = m_relaxed_from.find(to_configname);
        if(it == m_relaxed_from.end()) {
          if(_set.size()) {
            m_relaxed_from.insert({to_configname, _set});
          }
        }
        else {
          if(!_set.size()) {
            m_relaxed_from.erase(it);
          }
          else {
            it->second = _set;
          }
        }
      }

      std::set<std::string, Compare> _make_set(
        std::string to_configname,
        const ScoreMappedProperties &score) const {

        return std::set<std::string, Compare>(Compare(this, to_configname, score));
      }


      bool m_is_open;

      ScoreMappedProperties m_default_score;

      // a key whose 'from' value is modified to find MappedProperties in m_data
      MappedProperties m_key;

      // the MappedProperties container
      std::set<MappedProperties> m_data;

      // to -> {from, from, ...}, used to find best mapping
      std::map<std::string, std::set<std::string, Compare> > m_relaxed_from;

    };

  }
}

#endif
