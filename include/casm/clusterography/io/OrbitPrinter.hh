#ifndef CASM_cluterography_io_OrbitPrinter
#define CASM_cluterography_io_OrbitPrinter

#include "casm/casm_io/Log.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/clusterography/ClusterDecl.hh"
#include "casm/global/definitions.hh"
#include "casm/global/enum.hh"
#include "casm/symmetry/SymInfo.hh"

namespace CASM {

class IntegralCluster;
class Structure;

enum class ORBIT_PRINT_MODE { PROTO, FULL };

ENUM_IO_DECL(ORBIT_PRINT_MODE)
ENUM_JSON_IO_DECL(ORBIT_PRINT_MODE)
ENUM_TRAITS(ORBIT_PRINT_MODE)

struct OrbitPrinterOptions {
  int indent_space = 6;
  char delim = '\n';
  int prec = 7;
  COORD_TYPE coord_type = FRAC;
  ORBIT_PRINT_MODE orbit_print_mode = ORBIT_PRINT_MODE::PROTO;
  SymInfoOptions sym_info_opt;
  bool print_coordinates = true;
  bool print_equivalence_map = false;
  bool print_invariant_group = false;
};

jsonParser &to_json(const OrbitPrinterOptions &opt, jsonParser &json);

/// \brief Read from JSON
void from_json(OrbitPrinterOptions &opt, const jsonParser &json);

template <>
struct jsonConstructor<OrbitPrinterOptions> {
  static OrbitPrinterOptions from_json(const jsonParser &json);
};

struct PrinterBase {
  OrbitPrinterOptions opt;

  PrinterBase(const OrbitPrinterOptions &_opt = OrbitPrinterOptions());

  void coord_type(Log &out);

  void increase_indent(Log &out) const {
    out.increase_indent_spaces(opt.indent_space);
  }
  void decrease_indent(Log &out) const {
    out.decrease_indent_spaces(opt.indent_space);
  }

  template <typename OrbitType>
  void print_equivalence_map(const OrbitType &orbit, Index equiv_index,
                             Log &out) const;

  template <typename OrbitType>
  void print_equivalence_map(const OrbitType &orbit, Index equiv_index,
                             jsonParser &json) const;

  template <typename OrbitType>
  void print_equivalence_map(const OrbitType &orbit, Log &out) const;

  template <typename OrbitType, typename Element>
  void print_invariant_group(const OrbitType &orbit, const Element &element,
                             Log &out) const;

  template <typename OrbitType, typename Element>
  void print_invariant_group(const OrbitType &orbit, const Element &element,
                             jsonParser &json) const;
};

template <typename _Element>
struct Printer : public PrinterBase {
  typedef _Element Element;

  Printer(const OrbitPrinterOptions &_opt = OrbitPrinterOptions())
      : PrinterBase(_opt) {}

  void print(const Element &element, Log &out) const {
    xtal::COORD_MODE printer_mode(opt.coord_type);
    out << element;
  }
};

template <>
struct Printer<IntegralCluster> : public PrinterBase {
  typedef IntegralCluster Element;
  static const std::string element_name;

  Printer(const OrbitPrinterOptions &_opt = OrbitPrinterOptions())
      : PrinterBase(_opt) {}

  void print(const Element &element, Log &out) const;
};

typedef Printer<IntegralCluster> SitesPrinter;

template <typename _Element, ORBIT_PRINT_MODE>
struct OrbitPrinter {};

/// \brief Print Orbit<SymCompareType>, including only prototypes
template <typename _Element>
struct OrbitPrinter<_Element, ORBIT_PRINT_MODE::PROTO>
    : public Printer<_Element> {
  using Printer<_Element>::element_name;
  using Printer<_Element>::print;

  OrbitPrinter(const OrbitPrinterOptions &_opt = OrbitPrinterOptions())
      : Printer<_Element>(_opt) {}

  template <typename OrbitType>
  void operator()(const OrbitType &orbit, Log &out, Index orbit_index,
                  Index Norbits) const;

  template <typename OrbitType>
  jsonParser &to_json(const OrbitType &orbit, jsonParser &json,
                      Index orbit_index, Index Norbits) const;
};

template <typename _Element>
using PrototypePrinter = OrbitPrinter<_Element, ORBIT_PRINT_MODE::PROTO>;

typedef PrototypePrinter<IntegralCluster> ProtoSitesPrinter;

/// \brief Print Orbit<SymCompareType>, including all equivalents
template <typename _Element>
struct OrbitPrinter<_Element, ORBIT_PRINT_MODE::FULL>
    : public Printer<_Element> {
  using Printer<_Element>::element_name;
  using Printer<_Element>::print;

  OrbitPrinter(const OrbitPrinterOptions &_opt = OrbitPrinterOptions())
      : Printer<_Element>(_opt) {}

  template <typename OrbitType>
  void operator()(const OrbitType &orbit, Log &out, Index orbit_index,
                  Index Norbits) const;

  template <typename OrbitType>
  jsonParser &to_json(const OrbitType &orbit, jsonParser &json,
                      Index orbit_index, Index Norbits) const;
};

template <typename _Element>
using FullOrbitPrinter = OrbitPrinter<_Element, ORBIT_PRINT_MODE::FULL>;

typedef FullOrbitPrinter<IntegralCluster> FullSitesPrinter;

/// \brief Print IntegralCluster orbits
template <typename ClusterOrbitIterator, typename OrbitPrinter>
void print_clust(ClusterOrbitIterator begin, ClusterOrbitIterator end, Log &out,
                 OrbitPrinter printer);

/// \brief Print IntegralCluster orbits
template <typename ClusterOrbitIterator>
void print_clust(ClusterOrbitIterator begin, ClusterOrbitIterator end, Log &out,
                 const OrbitPrinterOptions &opt = OrbitPrinterOptions());

// ---------- clust.json IO
// ------------------------------------------------------------------

/// \brief Read JSON containing Orbit<SymCompareType> prototypes
template <typename ClusterOutputIterator, typename SymCompareType>
ClusterOutputIterator read_clust(ClusterOutputIterator result,
                                 const jsonParser &json, const Structure &prim,
                                 const SymGroup &generating_grp,
                                 const SymCompareType &sym_compare);

/// \brief Read JSON containing IntegralCluster prototypes, as IntegralCluster
template <typename ClusterOutputIterator>
ClusterOutputIterator read_clust(ClusterOutputIterator result,
                                 const jsonParser &json, const Structure &prim);

/// \brief Write Orbit<SymCompareType> to JSON
template <typename ClusterOrbitIterator, typename Printer>
jsonParser &write_clust(ClusterOrbitIterator begin, ClusterOrbitIterator end,
                        jsonParser &json, Printer printer);

/// \brief Write Orbit<SymCompareType> to JSON
template <typename ClusterOrbitIterator>
jsonParser &write_clust(ClusterOrbitIterator begin, ClusterOrbitIterator end,
                        jsonParser &json,
                        const OrbitPrinterOptions &opt = OrbitPrinterOptions());

/// \brief Write Orbit<SymCompareType> to JSON, including 'bspecs'
template <typename ClusterOrbitIterator, typename Printer>
jsonParser &write_clust(ClusterOrbitIterator begin, ClusterOrbitIterator end,
                        jsonParser &json, Printer printer,
                        const jsonParser &bspecs);

}  // namespace CASM

#endif
