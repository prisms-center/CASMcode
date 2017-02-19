#ifndef CASM_ScelDatabase
#define CASM_ScelDatabase

namespace {

  namespace DB {

    /// Supercell references and pointers must remain valid, therefore the
    /// ScelDatabase implementation should typically be the same for different
    /// backends, differing only for open, commit, and close. This class may
    /// be inherited by those classes to implement all other required methods.
    ///
    /// Derived ScelDatabase must implement public methods:
    /// - void DatabaseBase& open()
    /// - void commit()
    /// - void close()
    ///
    class ScelDatabase : public Database<Supercell> {

    public:

      iterator begin() override {
        return _iterator(m_scel_list.begin());
      }

      return iterator end() override {
        return _iterator(m_scel_list.begin());
      }

      size_type size() const override {
        return m_scel_list.size()
      }

      std::pair<iterator, bool> insert(const ValueType &obj) override {
        return _on_insert_or_emplace(m_scel_list.insert(obj));
      }

      std::pair<iterator, bool> insert(const ValueType &&obj) override {
        return _on_insert_or_emplace(m_scel_list.insert(std::move(obj)));
      }

      iterator erase(iterator pos) override {
        return _iterator(m_scel_list.erase(pos));
      }

      iterator find(const name_type &name_or_alias) override {
        return _iterator(m_scel_list.find(
      }


                     protected:

                       typedef DatabaseSetIterator<Supercell, ScelDatabase> base_iterator;

             template<typename... Args>
      std::pair<iterator, bool> emplace(Args &&... args) {
        return _on_insert_or_emplace(m_scel_list.emplace(std::forward<Args>(args)));
      }

      std::pair<iterator, bool> _on_insert_or_emplace(const std::pair<base_iterator, bool> &result) {

        if(result.second) {

          const value_type &obj = *result->second;

          // update name & alias
          m_name_and_alias.insert(std::make_pair(obj.name(), result.first));
          if(!obj.alias().empty()) {
            m_name_and_alias.insert(std::make_pair(obj.alias(), result.first));
          }

        }

        return std::make_pair(_iterator(result.first), result.second));
      }

      void clear() {
        m_name_or_alias.clear();
        m_scel_list.clear();
      }

      std::map<std::string, base_iterator> m_name_or_alias;
      std::set<Supercell> m_scel_list;

    };
  }
}

#endif
