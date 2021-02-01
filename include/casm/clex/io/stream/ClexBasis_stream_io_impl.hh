#ifndef CASM_clex_ClexBasis_stream_io_impl
#define CASM_clex_ClexBasis_stream_io_impl

#include "casm/basis_set/DoFTraits.hh"
#include "casm/casm_io/Log.hh"
#include "casm/clex/ClexBasis.hh"
#include "casm/clex/io/ProtoFuncsPrinter_impl.hh"
#include "casm/clex/io/stream/ClexBasis_stream_io.hh"
#include "casm/clusterography/io/OrbitPrinter_impl.hh"
#include "casm/crystallography/Structure.hh"

namespace CASM {

template <typename OrbitVecType>
void ClexBasisFunctionPrinter::operator()(OrbitVecType const &orbits) const {
  ParsingDictionary<DoFType::Traits> const *dof_dict = &DoFType::traits_dict();
  ClexBasis clex_basis{m_shared_prim, m_basis_set_specs, dof_dict};
  clex_basis.generate(orbits.begin(), orbits.end());

  print_site_basis_funcs(m_shared_prim, clex_basis, m_log);
  ProtoFuncsPrinter funcs_printer{clex_basis,
                                  m_shared_prim->shared_structure()};
  print_clust(orbits.begin(), orbits.end(), m_log, funcs_printer);
}

}  // namespace CASM

#endif
