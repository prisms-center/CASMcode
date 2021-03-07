#ifndef CASM_clex_ClexBasis_stream_io
#define CASM_clex_ClexBasis_stream_io

#include <memory>

#include "casm/clex/ClexBasisSpecs.hh"
#include "casm/clusterography/io/OrbitPrinter.hh"

namespace CASM {

class Log;
class PrimClex;
class Structure;

/// This functor implements a template method so that basis functions can be
/// printed for any orbit type (supported by `for_all_orbits`), as determined
/// at runtime
class ClexBasisFunctionPrinter {
 public:
  ClexBasisFunctionPrinter(Log &_log,
                           std::shared_ptr<Structure const> _shared_prim,
                           ClexBasisSpecs const &_basis_set_specs, bool align,
                           OrbitPrinterOptions _opt = OrbitPrinterOptions());

  template <typename OrbitVecType>
  void operator()(OrbitVecType const &orbits) const;

 private:
  std::shared_ptr<Structure const> m_shared_prim;
  ClexBasisSpecs m_basis_set_specs;
  Log &m_log;
  bool m_align;
  OrbitPrinterOptions m_orbit_printer_opt;
};

/// Pretty-print basis functions -- generate, then print
void print_basis_functions(Log &log,
                           std::shared_ptr<Structure const> const &shared_prim,
                           ClexBasisSpecs const &basis_set_specs, bool align,
                           OrbitPrinterOptions opt = OrbitPrinterOptions());

/// Pretty-print basis functions -- read clusters from existing clust.json file
void print_basis_functions(Log &log, PrimClex const &primclex,
                           std::string basis_set_name, bool align,
                           OrbitPrinterOptions opt = OrbitPrinterOptions());

}  // namespace CASM

#endif
