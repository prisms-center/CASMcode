#ifndef CASM_cluterography_io_OrbitPrinter_impl
#define CASM_cluterography_io_OrbitPrinter_impl

#include "casm/casm_io/container/json_io.hh"
#include "casm/clusterography/ClusterOrbits_impl.hh"
#include "casm/clusterography/io/OrbitPrinter.hh"
#include "casm/clusterography/io/json/IntegralCluster_json_io.hh"
#include "casm/crystallography/io/BasicStructureIO.hh"
#include "casm/symmetry/InvariantSubgroup_impl.hh"
#include "casm/symmetry/SymTools.hh"

namespace CASM {

template <typename OrbitPrinter, typename Element>
void print_coordinates(OrbitPrinter &printer, const Element &element,
                       Log &out) {
  out << out.indent_str() << "Coordinates:" << std::endl;
  printer.increase_indent(out);
  printer.print(element, out);
  printer.decrease_indent(out);
}

template <typename OrbitType>
void PrinterBase::print_equivalence_map(const OrbitType &orbit,
                                        Index equiv_index, Log &out) const {
  out << out.indent_str() << "Equivalence map: " << std::endl;
  this->increase_indent(out);
  int j = 0;
  SymInfoOptions topt = this->opt.sym_info_opt;
  for (const auto &op : orbit.equivalence_map()[equiv_index]) {
    out << out.indent_str() << j << ": (" << op.index() << ") "
        << brief_description(op, orbit.prototype().prim().lattice(), topt)
        << std::endl;
    ++j;
  }
  this->decrease_indent(out);
}

template <typename OrbitType>
void PrinterBase::print_equivalence_map(const OrbitType &orbit,
                                        Index equiv_index,
                                        jsonParser &json) const {
  SymInfoOptions topt = this->opt.sym_info_opt;

  std::vector<Index> master_group_indices;
  for (const auto &op : orbit.equivalence_map()[equiv_index]) {
    master_group_indices.push_back(op.index());
  }
  json["equivalence_map"] = master_group_indices;

  std::vector<std::string> desc;
  for (const auto &op : orbit.equivalence_map()[equiv_index]) {
    desc.push_back(
        brief_description(op, orbit.prototype().prim().lattice(), topt));
  }
  json["equivalence_map_descriptions"] = desc;
  json["equivalence_map_descriptions"].set_force_column();
}

template <typename OrbitType>
void PrinterBase::print_equivalence_map(const OrbitType &orbit,
                                        Log &out) const {
  out << out.indent_str() << "Orbit equivalence map: " << std::endl;
  this->increase_indent(out);
  SymInfoOptions topt = this->opt.sym_info_opt;
  for (int i = 0; i < orbit.size(); ++i) {
    out << out.indent_str() << "Element: " << i << std::endl;
    this->increase_indent(out);
    int j = 0;
    for (const auto &op : orbit.equivalence_map()[i]) {
      out << out.indent_str() << j << ": (" << op.index() << ") "
          << brief_description(op, orbit.prototype().prim().lattice(), topt)
          << std::endl;
      ++j;
    }
    this->decrease_indent(out);
  }
  this->decrease_indent(out);
}

template <typename OrbitType, typename Element>
void PrinterBase::print_invariant_group(const OrbitType &orbit,
                                        const Element &element,
                                        Log &out) const {
  out << out.indent_str() << "Invariant group:" << std::endl;
  SymGroup invariant_group = make_invariant_subgroup(
      element, orbit.generating_group(), orbit.sym_compare());
  this->increase_indent(out);
  SymInfoOptions topt = this->opt.sym_info_opt;
  brief_description(out, invariant_group, orbit.prototype().prim().lattice(),
                    topt);
  this->decrease_indent(out);
}

template <typename OrbitType, typename Element>
void PrinterBase::print_invariant_group(const OrbitType &orbit,
                                        const Element &element,
                                        jsonParser &json) const {
  SymInfoOptions topt = this->opt.sym_info_opt;
  SymGroup invariant_group = make_invariant_subgroup(
      element, orbit.generating_group(), orbit.sym_compare());

  std::vector<Index> master_group_indices;
  for (const auto &op : invariant_group) {
    master_group_indices.push_back(op.index());
  }
  json["invariant_group"] = master_group_indices;

  std::vector<std::string> desc;
  for (const auto &op : invariant_group) {
    desc.push_back(
        brief_description(op, orbit.prototype().prim().lattice(), topt));
  }
  json["invariant_group_descriptions"] = desc;
  json["invariant_group_descriptions"].set_force_column();
}

// --- OrbitPrinter templates ---

template <typename _Element>
template <typename OrbitType>
void OrbitPrinter<_Element, ORBIT_PRINT_MODE::PROTO>::operator()(
    const OrbitType &orbit, Log &out, Index orbit_index, Index Norbits) const {
  out << out.indent_str() << "Prototype"
      << " of " << orbit.size() << " Equivalent " << element_name
      << " in Orbit " << orbit_index << std::endl;
  this->increase_indent(out);

  if (this->opt.print_coordinates) {
    print_coordinates(*this, orbit.prototype(), out);
  }
  if (this->opt.print_invariant_group) {
    this->print_invariant_group(orbit, orbit.prototype(), out);
  }
  this->decrease_indent(out);
}

/// \brief Print to JSON
///
/// Note: for 'read_clust' to work, "prototype" must be written
template <typename _Element>
template <typename OrbitType>
jsonParser &OrbitPrinter<_Element, ORBIT_PRINT_MODE::PROTO>::to_json(
    const OrbitType &orbit, jsonParser &json, Index orbit_index,
    Index Norbits) const {
  json.put_obj();
  json["linear_orbit_index"] = orbit_index;
  json["prototype"] = orbit.prototype();
  if (this->opt.print_invariant_group) {
    this->print_invariant_group(orbit, orbit.prototype(), json["prototype"]);
  }

  return json;
}

template <typename _Element>
template <typename OrbitType>
void OrbitPrinter<_Element, ORBIT_PRINT_MODE::FULL>::operator()(
    const OrbitType &orbit, Log &out, Index orbit_index, Index Norbits) const {
  for (Index equiv_index = 0; equiv_index != orbit.size(); ++equiv_index) {
    out << out.indent_str() << equiv_index << " of " << orbit.size()
        << " Equivalent " << element_name << " in Orbit " << orbit_index
        << std::endl;
    this->increase_indent(out);

    if (this->opt.print_coordinates) {
      print_coordinates(*this, orbit[equiv_index], out);
    }
    if (this->opt.print_invariant_group) {
      this->print_invariant_group(orbit, orbit[equiv_index], out);
    }
    if (this->opt.print_equivalence_map) {
      this->print_equivalence_map(orbit, equiv_index, out);
    }
    this->decrease_indent(out);
  }
}

/// \brief Print to JSON
///
/// Note: for 'read_clust' to work, "prototype" must be written
template <typename _Element>
template <typename OrbitType>
jsonParser &OrbitPrinter<_Element, ORBIT_PRINT_MODE::FULL>::to_json(
    const OrbitType &orbit, jsonParser &json, Index orbit_index,
    Index Norbits) const {
  json.put_obj();
  json["linear_orbit_index"] = orbit_index;
  json["prototype"] = orbit.prototype();
  if (this->opt.print_invariant_group) {
    this->print_invariant_group(orbit, orbit.prototype(), json["prototype"]);
  }
  json["elements"].put_array(orbit.begin(), orbit.end());
  for (Index equiv_index = 0; equiv_index != orbit.size(); ++equiv_index) {
    if (this->opt.print_invariant_group) {
      this->print_invariant_group(orbit, orbit[equiv_index],
                                  json["elements"][equiv_index]);
    }
    if (this->opt.print_equivalence_map) {
      this->print_equivalence_map(orbit, equiv_index,
                                  json["elements"][equiv_index]);
    }
  }
  return json;
}

/// \brief Print IntegralCluster orbits
///
/// \param begin,end Range of Orbit<SymCompareType>
/// \param out output stream
/// \param mode Coordinate output mode
/// \param printer A functor to control printing for the orbit
///
/// Printer is expected to have:
/// - \code std::string Printer::indent(); \endcode
/// - \code void Printer::coord_type(Log& out); \endcode
/// - \code void Printer::operator()(const Orbit<SymCompareType>& orbit, Log&
/// out, Index orbit_index, Index Norbits); \endcode
///
template <typename ClusterOrbitIterator, typename OrbitPrinter>
void print_clust(ClusterOrbitIterator begin, ClusterOrbitIterator end, Log &out,
                 OrbitPrinter printer) {
  printer.coord_type(out);

  out.ostream().flags(std::ios::showpoint | std::ios::fixed | std::ios::left);
  out.ostream().precision(5);

  Index branch = -1;
  Index orbit_index = 0;
  Index Norbits = std::distance(begin, end);

  for (auto it = begin; it != end; ++it) {
    if (it->prototype().size() != branch) {
      branch = it->prototype().size();
      out << out.indent_str() << "** Branch " << branch << " ** " << std::endl;
    }
    printer.increase_indent(out);
    double min_length = 0.0;
    double max_length = 0.0;
    if (it->invariants().displacement().size()) {
      min_length = it->invariants().displacement().front();
      max_length = it->invariants().displacement().back();
    }
    out << out.indent_str() << "** " << orbit_index << " of " << Norbits
        << " Orbits **"
        << "  Points: " << it->prototype().size() << "  Mult: " << it->size()
        << "  MinLength: " << min_length << "  MaxLength: " << max_length
        << std::endl;
    printer.increase_indent(out);
    printer(*it, out, orbit_index, Norbits);
    out << std::endl;
    printer.decrease_indent(out);
    printer.decrease_indent(out);
    ++orbit_index;
  }
}

/// \brief Print IntegralCluster orbits
template <typename ClusterOrbitIterator>
void print_clust(ClusterOrbitIterator begin, ClusterOrbitIterator end, Log &out,
                 const OrbitPrinterOptions &opt) {
  // typedef typename ClusterOrbitIterator::container_type container_type;
  // typedef typename container_type::value_type orbit_type;
  typedef typename std::iterator_traits<ClusterOrbitIterator>::value_type
      orbit_type;
  typedef typename orbit_type::Element Element;

  if (opt.orbit_print_mode == ORBIT_PRINT_MODE::PROTO) {
    OrbitPrinter<Element, ORBIT_PRINT_MODE::PROTO> orbit_printer(opt);
    print_clust(begin, end, out, orbit_printer);
  } else if (opt.orbit_print_mode == ORBIT_PRINT_MODE::FULL) {
    OrbitPrinter<Element, ORBIT_PRINT_MODE::FULL> orbit_printer(opt);
    print_clust(begin, end, out, orbit_printer);
  }
}

// ---------- clust.json IO
// ------------------------------------------------------------------

/// \brief Read JSON containing Orbit<SymCompareType> prototypes
///
/// - Uses 'prim', 'generating_grp', and 'sym_compare' to generate orbits from
///   prototypes read from the JSON
/// - Ignores "prim" and "bspecs" info in the JSON
///
template <typename ClusterOutputIterator, typename SymCompareType>
ClusterOutputIterator read_clust(ClusterOutputIterator result,
                                 const jsonParser &json, const Structure &prim,
                                 const SymGroup &generating_grp,
                                 const SymCompareType &sym_compare) {
  typedef Orbit<SymCompareType> orbit_type;

  for (const auto &j : json["orbits"]) {
    *result++ = orbit_type(j["prototype"].get<IntegralCluster>(prim),
                           generating_grp, sym_compare);
  }
  return result;
}

/// \brief Read JSON containing IntegralCluster prototypes, as
/// Orbit<SymCompareType>
///
/// - Uses 'prim' to generate IntegralCluster from prototypes read from the JSON
/// - Ignores "prim" and "bspecs" info in the JSON
///
template <typename ClusterOutputIterator>
ClusterOutputIterator read_clust(ClusterOutputIterator result,
                                 const jsonParser &json,
                                 const Structure &prim) {
  for (const auto &j : json["orbits"]) {
    *result++ = j["prototype"].get<IntegralCluster>(prim);
  }
  return result;
}

/// \brief Write Orbit<SymCompareType> to JSON
///
/// Format:
/// \code
/// {
///   "orbits": [
///     {
///       "prototype" : {
///         "min_length" : number,
///         "max_length" : number,
///         "sites" : (JSON array of UnitCellCoord)
///       }
///     },
///     ... for each orbit ...
///   ],
///   "prim" : (JSON object, the contents of prim.json)
/// }
/// \endcode
///
template <typename ClusterOrbitIterator, typename Printer>
jsonParser &write_clust(ClusterOrbitIterator begin, ClusterOrbitIterator end,
                        jsonParser &json, Printer printer) {
  Index Norbits = std::distance(begin, end);
  json["orbits"] = jsonParser::array(Norbits, jsonParser::object());
  Index orbit_index = 0;
  for (auto it = begin; it != end; ++it, ++orbit_index) {
    printer.to_json(*it, json["orbits"][orbit_index], orbit_index, Norbits);
  }
  write_prim(begin->prototype().prim(), json["prim"], FRAC);
  return json;
}

/// \brief Write Orbit<SymCompareType> to JSON, including 'bspecs'
template <typename ClusterOrbitIterator>
jsonParser &write_clust(ClusterOrbitIterator begin, ClusterOrbitIterator end,
                        jsonParser &json, const OrbitPrinterOptions &opt) {
  // typedef typename ClusterOrbitIterator::container_type container_type;
  // typedef typename container_type::value_type orbit_type;
  typedef typename std::iterator_traits<ClusterOrbitIterator>::value_type
      orbit_type;
  typedef typename orbit_type::Element Element;

  if (opt.orbit_print_mode == ORBIT_PRINT_MODE::PROTO) {
    OrbitPrinter<Element, ORBIT_PRINT_MODE::PROTO> orbit_printer(opt);
    write_clust(begin, end, json, orbit_printer);
  } else if (opt.orbit_print_mode == ORBIT_PRINT_MODE::FULL) {
    OrbitPrinter<Element, ORBIT_PRINT_MODE::FULL> orbit_printer(opt);
    write_clust(begin, end, json, orbit_printer);
  }
  return json;
}

/// \brief Write Orbit<SymCompareType> to JSON, including 'bspecs'
///
/// Format:
/// \code
/// {
///   "orbits": [
///     {
///       "prototype" : {
///         "min_length" : number,
///         "max_length" : number,
///         "sites" : (JSON array of UnitCellCoord)
///       }
///     },
///     ... for each orbit ...
///   ],
///   "bspecs" : (JSON object, the contents of bspecs.json)
///   "prim" : (JSON object, the contents of prim.json)
/// }
/// \endcode
///
template <typename ClusterOrbitIterator, typename Printer>
jsonParser &write_clust(ClusterOrbitIterator begin, ClusterOrbitIterator end,
                        jsonParser &json, Printer printer,
                        const jsonParser &bspecs) {
  write_clust(begin, end, json, printer);
  json["bspecs"] = bspecs;
  return json;
}

}  // namespace CASM

#endif
