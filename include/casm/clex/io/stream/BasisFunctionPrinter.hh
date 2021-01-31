#ifndef CASM_clex_io_stream_BasisFunctionPrinter
#define CASM_clex_io_stream_BasisFunctionPrinter

#include <memory>

#include "casm/clex/ClexBasisSpecs.hh"

namespace CASM {

class Log;
class Structure;

/// This functor implements a template method so that basis functions can be
/// printed for any orbit type, as determined at runtime from the JSON file
/// parameters
class BasisFunctionPrinter {
 public:
  BasisFunctionPrinter(Log &_log, std::shared_ptr<Structure const> _shared_prim,
                       ClexBasisSpecs const &_basis_set_specs);

  template <typename OrbitVecType>
  void operator()(OrbitVecType const &orbits) const;

 private:
  std::shared_ptr<Structure const> m_shared_prim;
  ClexBasisSpecs m_basis_set_specs;
  Log &m_log;
};

}  // namespace CASM

#endif
