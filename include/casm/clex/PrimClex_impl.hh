#ifndef CASM_PrimClex_impl
#define CASM_PrimClex_impl

#include "casm/casm_io/jsonParser.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/app/DirectoryStructure.hh"
#include "casm/app/ClexDescription.hh"
#include "casm/app/ProjectSettings.hh"
#include "casm/crystallography/Structure.hh"

namespace CASM {

  template<typename OrbitOutputIterator, typename SymCompareType>
  OrbitOutputIterator PrimClex::orbits(
    const ClexDescription &key,
    OrbitOutputIterator result,
    const SymCompareType &sym_compare) const {

    return read_clust(
             result,
             jsonParser(dir().clust(key.bset)),
             prim(),
             prim().factor_group(),
             sym_compare,
             settings().crystallography_tol());
  }

}

#endif
