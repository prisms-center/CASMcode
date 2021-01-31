#ifndef CASM_clex_io_stream_BasisFunctionPrinter_impl
#define CASM_clex_io_stream_BasisFunctionPrinter_impl

#include "casm/casm_io/Log.hh"
#include "casm/clex/ClexBasis.hh"
#include "casm/clex/io/OrbitPrinter.hh"
#include "casm/clex/io/stream/BasisFunctionPrinter.hh"
#include "casm/clex/io/stream/ProtoFuncsPrinter.hh"
#include "casm/crystallography/Structure.hh"

namespace CASM {

BasisFunctionPrinter::BasisFunctionPrinter(
    Log &_log, std::shared_ptr<Structure const> _shared_prim,
    ClexBasisSpecs const &_basis_set_specs)
    : m_shared_prim(_shared_prim),
      m_basis_set_specs(_basis_set_specs),
      m_log(_log) {}

template <typename OrbitVecType>
void BasisFunctionPrinter::operator()(OrbitVecType const &orbits) const {
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
