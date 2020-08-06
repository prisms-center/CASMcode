#ifndef CASM_PrimClex_impl
#define CASM_PrimClex_impl

#include "casm/app/ClexDescription.hh"
#include "casm/app/DirectoryStructure.hh"
#include "casm/app/ProjectSettings.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/crystallography/Structure.hh"

namespace CASM {

  template <typename OrbitOutputIterator, typename SymCompareType>
  OrbitOutputIterator
  PrimClex::orbits(std::string const &basis_set_name, OrbitOutputIterator result, SymCompareType const &sym_compare) const {

    return read_clust(result, jsonParser {dir().clust(basis_set_name)}, prim(), prim().factor_group(),
                      sym_compare, settings().crystallography_tol());
  }

} // namespace CASM

#endif
