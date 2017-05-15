#ifndef CASM_ScelEnum_impl
#define CASM_ScelEnum_impl

#include "casm/clex/ScelEnum.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/database/ScelDatabase.hh"

namespace CASM {

  /// \brief Construct with PrimClex and a range of Supercell names
  ///
  /// \param primclex A PrimClex for which to enumerate Supercells
  /// \param begin,end A range of names of Supercells to enumerate
  ///
  template<typename ScelNameIterator>
  ScelEnumByName::ScelEnumByName(
    const PrimClex &primclex,
    ScelNameIterator begin,
    ScelNameIterator end) :
    RandomAccessEnumeratorBase<Supercell>(),
    m_primclex(&primclex) {

    for(auto it = begin; it != end; ++it) {
      m_scelptr.push_back(&*m_primclex->db<Supercell>().find(*it));
    }
    _init();
  }

}

#endif
