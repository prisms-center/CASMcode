#ifndef CASM_clex_io_ProtoFuncsPrinter
#define CASM_clex_io_ProtoFuncsPrinter

#include <vector>

#include "casm/basis_set/FunctionVisitor.hh"
#include "casm/clusterography/io/OrbitPrinter.hh"
#include "casm/global/definitions.hh"

namespace CASM {

namespace xtal {
class BasicStructure;
}

class ClexBasis;
class Structure;
class jsonParser;

/// \brief Print Orbit<SymCompareType> & ClexBasis, including prototypes and
/// prototype basis functions
struct ProtoFuncsPrinter : public SitesPrinter {
  typedef xtal::BasicStructure PrimType;
  typedef std::shared_ptr<const PrimType> PrimType_ptr;

  ClexBasis const &clex_basis;

  std::vector<SubExpressionLabeler> labelers;

  ProtoFuncsPrinter(ClexBasis const &_clex_basis, PrimType_ptr prim_ptr,
                    bool align,
                    OrbitPrinterOptions const &_opt = OrbitPrinterOptions());

  /// \brief Print to JSON
  ///
  /// Note: for 'read_clust' to work, "prototype" must be written
  template <typename OrbitType>
  void operator()(const OrbitType &orbit, Log &out, Index orbit_index,
                  Index Norbits) const;

  template <typename OrbitType>
  jsonParser &to_json(const OrbitType &orbit, jsonParser &json,
                      Index orbit_index, Index Norbits) const;

 private:
  PrimType_ptr prim_ptr;
  bool m_align;
};

/// Print prototype cluster sites as a tex tabular
void print_tex_tabular_cluster_sites(Log &out, IntegralCluster const &cluster,
                                     xtal::BasicStructure const &prim,
                                     COORD_TYPE mode);

/// Print site basis functions, as for 'casm bset --functions'
void print_site_basis_funcs(std::shared_ptr<const Structure> prim_ptr,
                            ClexBasis const &clex_basis, Log &out,
                            Index indent_space = 6, COORD_TYPE mode = FRAC);

/// Print aligned site basis functions, as for 'casm bset --functions --align'
void print_aligned_site_basis_funcs(std::shared_ptr<const Structure> prim_ptr,
                                    ClexBasis const &clex_basis, Log &out,
                                    Index indent_space = 6,
                                    COORD_TYPE mode = FRAC);

void write_site_basis_funcs(std::shared_ptr<const Structure> prim_ptr,
                            ClexBasis const &clex_basis, jsonParser &json);

}  // namespace CASM

#endif
