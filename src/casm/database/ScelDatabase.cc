#include "casm/database/ScelDatabase.hh"
#include "casm/clex/PrimClex.hh"

namespace CASM {
  namespace DB {

    Database<Supercell>::iterator Database<Supercell>::begin() const {
      return _iterator(m_scel_list.begin());
    }

    Database<Supercell>::iterator Database<Supercell>::end() const {
      return _iterator(m_scel_list.end());
    }

    Database<Supercell>::size_type Database<Supercell>::size() const {
      return m_scel_list.size();
    }

    std::pair<Database<Supercell>::iterator, bool>
    Database<Supercell>::insert(const Supercell &obj) {
      return _on_insert_or_emplace(m_scel_list.insert(obj));
    }

    std::pair<Database<Supercell>::iterator, bool>
    Database<Supercell>::insert(const Supercell &&obj) {
      return _on_insert_or_emplace(m_scel_list.insert(std::move(obj)));
    }

    Database<Supercell>::iterator Database<Supercell>::erase(iterator pos) {
      typedef DatabaseSetIterator<Supercell, Database<Supercell> > it_type;
      base_iterator base_it = static_cast<it_type *>(pos.get())->base();
      m_name_to_scel.erase(base_it->name());
      master_selection().data().erase(base_it->name());
      return _iterator(m_scel_list.erase(base_it));
    }

    Database<Supercell>::iterator
    Database<Supercell>::find(const std::string &name_or_alias) const {
      std::string name = this->name(name_or_alias);
      auto it = m_name_to_scel.find(name);
      if(it == m_name_to_scel.end()) {
        return _iterator(m_scel_list.end());
      }
      else {
        return _iterator(it->second);
      }
    }

    std::pair<Database<Supercell>::iterator, bool>
    Database<Supercell>::_on_insert_or_emplace(
      const std::pair<base_iterator, bool> &result) {

      if(result.second) {

        const value_type &obj = *result.first;

        // update
        m_name_to_scel.insert(std::make_pair(obj.name(), result.first));
        master_selection().data().emplace(obj.name(), false);
      }

      return std::make_pair(_iterator(result.first), result.second);
    }

    void Database<Supercell>::clear() {
      m_name_to_scel.clear();
      m_scel_list.clear();
    }

    Database<Supercell>::iterator
    Database<Supercell>::_iterator(base_iterator base_it) const {
      return iterator(DatabaseSetIterator<Supercell, Database<Supercell> >(base_it));
    }

    /// Make canonical supercell and insert into supercell database
    std::pair<iterator, bool> make_canonical_and_insert(
      PrimClex const *primclex,
      Lattice const &super_lattice,
      Database<Supercell> &supercell_db) {

      auto const &prim = primclex->prim();
      auto const &pg = prim.point_group();
      double xtal_tol = prim.lattice().tol();
      xtal::Lattice canonical_lattice = xtal::canonical::equivalent(super_lattice, pg, xtal_tol);
      return supercell_db.emplace(primclex, canonical_lattice);
    }

    /// Make canonical supercell and insert into supercell database
    std::pair<iterator, bool> make_canonical_and_insert(
      PrimClex const *primclex,
      Eigen::Matrix3l const &transformation_matrix_to_super,
      Database<Supercell> &supercell_db) {
      auto const &prim = primclex->prim();
      xtal::Lattice super_lattice = make_superlattice(prim.lattice(), transformation_matrix_to_super);
      return make_canonical_and_insert(primclex, super_lattice, supercell_db);
    }

    /// Make canonical supercell and insert into supercell database
    std::pair<iterator, bool> make_canonical_and_insert(
      std::shared_ptr<Structure const> const &shared_prim,
      Lattice const &super_lattice,
      Database<Supercell> &supercell_db) {

      auto const &pg = shared_prim->point_group();
      double xtal_tol = shared_prim->lattice().tol();
      xtal::Lattice canonical_lattice = xtal::canonical::equivalent(super_lattice, pg, xtal_tol);
      return supercell_db.emplace(shared_prim, canonical_lattice);
    }

    /// Make canonical supercell and insert into supercell database
    std::pair<iterator, bool> make_canonical_and_insert(
      std::shared_ptr<Structure const> const &shared_prim,
      Eigen::Matrix3l const &transformation_matrix_to_super,
      Database<Supercell> &supercell_db) {

      xtal::Lattice super_lattice = make_superlattice(shared_prim->lattice(), transformation_matrix_to_super);
      return make_canonical_and_insert(shared_prim, super_lattice, supercell_db);
    }

  }
}
