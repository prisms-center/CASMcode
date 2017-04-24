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
    /// backends, differing only for open, commit, and close. This class should
    /// be inherited by those classes to implement all other required methods.
    ///
    /// Database insert and emplace methods by convention do not enforce
    /// canonical forms. That logic is included in Supercell::insert, which is
    /// the safest way to insert new Supercell in the database. But in cases
    /// where it is known that a Supercell is generated in canonical form, the
    /// Database insert and emplace methods may be used directly.
    ///
    /// Derived ScelDatabase must implement public methods:
    /// - void DatabaseBase& open()
    /// - void commit()
    /// - void close()
    ///
    template<>
    class Database<Supercell> : public ValDatabase<Supercell> {

    public:

      Database(const PrimClex &_primclex) :
        ValDatabase<Supercell>(_primclex) {}

      virtual ~Database() {}


      iterator begin() const override;

      iterator end() const override;

      size_type size() const override;

      std::pair<iterator, bool> insert(const Supercell &obj) override;

      std::pair<iterator, bool> insert(const Supercell &&obj) override;

      template<typename... Args>
      std::pair<iterator, bool> emplace(Args &&... args) {
        return _on_insert_or_emplace(m_scel_list.emplace(std::forward<Args>(args)...));
      }

      iterator erase(iterator pos) override;

      iterator find(const std::string &name_or_alias) const override;

      iterator find(const Supercell &obj) const override;

    protected:

      typedef std::set<Supercell>::iterator base_iterator;

      std::pair<iterator, bool> _on_insert_or_emplace(const std::pair<base_iterator, bool> &result);

      void clear();

      iterator _iterator(base_iterator base_it) const;

      std::map<std::string, base_iterator> m_name_to_scel;
      std::set<Supercell> m_scel_list;

    };
  }
}

#endif
