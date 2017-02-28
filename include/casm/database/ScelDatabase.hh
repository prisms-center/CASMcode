#ifndef CASM_ScelDatabase
#define CASM_ScelDatabase

#include <map>
#include <set>
#include "casm/database/Database.hh"
#include "casm/database/DatabaseSetIterator.hh"
#include "casm/clex/Supercell.hh"

namespace CASM {

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
    template<>
    class Database<Supercell> : public ValDatabase<Supercell> {

    public:

      iterator begin() const override {
        return _iterator(m_scel_list.begin());
      }

      iterator end() const override {
        return _iterator(m_scel_list.begin());
      }

      size_type size() const override {
        return m_scel_list.size();
      }

      std::pair<iterator, bool> insert(const Supercell &obj) override {
        return _on_insert_or_emplace(m_scel_list.insert(obj));
      }

      std::pair<iterator, bool> insert(const Supercell &&obj) override {
        return _on_insert_or_emplace(m_scel_list.insert(std::move(obj)));
      }

      template<typename... Args>
      std::pair<iterator, bool> emplace(Args &&... args) {
        return _on_insert_or_emplace(m_scel_list.emplace(std::forward<Args>(args)...));
      }

      iterator erase(iterator pos) override {
        typedef DatabaseSetIterator<Supercell, Database<Supercell> > it_type;
        base_iterator base_it = static_cast<it_type *>(pos.get())->base();
        return _iterator(m_scel_list.erase(base_it));
      }

      iterator find(const std::string &name_or_alias) const override {
        return _iterator(m_name_or_alias.find(name_or_alias)->second);
      }

      iterator find(const Supercell &obj) const override {
        return _iterator(m_scel_list.find(obj));
      }

    protected:

      typedef std::set<Supercell>::iterator base_iterator;

      std::pair<iterator, bool> _on_insert_or_emplace(const std::pair<base_iterator, bool> &result) {

        if(result.second) {

          const value_type &obj = *result.first;

          // update name & alias
          m_name_or_alias.insert(std::make_pair(obj.name(), result.first));
          if(!obj.alias().empty()) {
            m_name_or_alias.insert(std::make_pair(obj.alias(), result.first));
          }

        }

        return std::make_pair(_iterator(result.first), result.second);
      }

      void clear() {
        m_name_or_alias.clear();
        m_scel_list.clear();
      }

      iterator _iterator(base_iterator base_it) const {
        return iterator(DatabaseSetIterator<Supercell, Database<Supercell> >(base_it));
      }

      std::map<std::string, base_iterator> m_name_or_alias;
      std::set<Supercell> m_scel_list;

    };
  }
}

#endif
