#ifndef CASM_clex_ClexBasis_stream_io_impl
#define CASM_clex_ClexBasis_stream_io_impl

#include "casm/basis_set/DoFTraits.hh"
#include "casm/casm_io/Log.hh"
#include "casm/clex/ClexBasis_impl.hh"
#include "casm/clex/io/ProtoFuncsPrinter_impl.hh"
#include "casm/clex/io/stream/ClexBasis_stream_io.hh"
#include "casm/clusterography/io/OrbitPrinter_impl.hh"
#include "casm/crystallography/Structure.hh"

namespace CASM {

template <typename OrbitVecType>
void ClexBasisFunctionPrinter::operator()(OrbitVecType const &orbits) const {
  if (m_json) {
    this->to_json(orbits);
  } else {
    this->to_stream(orbits);
  }
}

template <typename OrbitVecType>
void ClexBasisFunctionPrinter::to_stream(OrbitVecType const &orbits) const {
  ParsingDictionary<DoFType::Traits> const *dof_dict = &DoFType::traits_dict();
  ClexBasis clex_basis{m_shared_prim, m_basis_set_specs, dof_dict};
  clex_basis.generate(orbits.begin(), orbits.end());

  if (m_align) {
    print_aligned_site_basis_funcs(m_shared_prim, clex_basis, m_log,
                                   m_orbit_printer_opt.indent_space,
                                   m_orbit_printer_opt.coord_type);
  } else {
    print_site_basis_funcs(m_shared_prim, clex_basis, m_log,
                           m_orbit_printer_opt.indent_space,
                           m_orbit_printer_opt.coord_type);
  }
  ProtoFuncsPrinter funcs_printer{clex_basis, m_shared_prim->shared_structure(),
                                  m_align, m_orbit_printer_opt};
  m_log << "Prototype cluster functions: \n\n";
  print_clust(orbits.begin(), orbits.end(), m_log, funcs_printer);
}

template <typename OrbitVecType>
void ClexBasisFunctionPrinter::to_json(OrbitVecType const &orbits) const {
  ParsingDictionary<DoFType::Traits> const *dof_dict = &DoFType::traits_dict();
  ClexBasis clex_basis{m_shared_prim, m_basis_set_specs, dof_dict};
  clex_basis.generate(orbits.begin(), orbits.end());

  jsonParser json;
  write_site_basis_funcs(clex_basis.shared_prim(), clex_basis, json);
  bool align = false;
  ProtoFuncsPrinter funcs_printer{clex_basis, m_shared_prim->shared_structure(),
                                  align, m_orbit_printer_opt};
  write_clust(orbits.begin(), orbits.end(), json, funcs_printer);
  m_log << json << std::endl;
}

}  // namespace CASM

#endif
