#include "casm/clusterography/IntegralCluster_impl.hh"
#include "casm/clusterography/io/OrbitPrinter_impl.hh"
#include "casm/global/enum/json_io.hh"
#include "casm/global/enum/stream_io.hh"

namespace CASM {

const std::string traits<ORBIT_PRINT_MODE>::name = "orbit_print_mode";

const std::multimap<ORBIT_PRINT_MODE, std::vector<std::string> >
    traits<ORBIT_PRINT_MODE>::strval = {
        {ORBIT_PRINT_MODE::PROTO, {"PROTO", "Proto", "proto"}},
        {ORBIT_PRINT_MODE::FULL, {"FULL", "Full", "full"}}};

ENUM_IO_DEF(ORBIT_PRINT_MODE)
ENUM_JSON_IO_DEF(ORBIT_PRINT_MODE)

jsonParser &to_json(const OrbitPrinterOptions &opt, jsonParser &json) {
  json.put_obj();
  json["indent_space"] = opt.indent_space;
  // just keep default delim
  json["prec"] = opt.prec;
  json[traits<COORD_TYPE>::name] = opt.coord_type;
  json[traits<ORBIT_PRINT_MODE>::name] = opt.orbit_print_mode;
  json["print_coordinates"] = opt.print_coordinates;
  json["print_equivalence_map"] = opt.print_equivalence_map;
  json["print_invariant_group"] = opt.print_invariant_group;
  json["sym_info_opt"] = opt.sym_info_opt;
  return json;
}

/// \brief Read from JSON
void from_json(OrbitPrinterOptions &opt, const jsonParser &json) {
  json.get_if(opt.indent_space, "indent_space");
  // just keep default delim
  json.get_if(opt.prec, "prec");
  json.get_if(opt.coord_type, traits<COORD_TYPE>::name);
  json.get_if(opt.orbit_print_mode, traits<ORBIT_PRINT_MODE>::name);
  json.get_if(opt.print_coordinates, "print_coordinates");
  json.get_if(opt.print_equivalence_map, "print_equivalence_map");
  json.get_if(opt.print_invariant_group, "print_invariant_group");
  json.get_if(opt.sym_info_opt, "sym_info_opt");
}

OrbitPrinterOptions jsonConstructor<OrbitPrinterOptions>::from_json(
    const jsonParser &json) {
  OrbitPrinterOptions res;
  CASM::from_json(res, json);
  return res;
}

PrinterBase::PrinterBase(const OrbitPrinterOptions &_opt) : opt(_opt) {}

void PrinterBase::coord_type(Log &out) {
  out << out.indent_str() << "COORD_MODE = " << opt.coord_type << std::endl
      << std::endl;
}

const std::string Printer<IntegralCluster>::element_name = "Clusters";

void Printer<IntegralCluster>::print(const IntegralCluster &clust,
                                     Log &out) const {
  if (!out.print()) {
    return;
  }

  COORD_TYPE _mode = this->opt.coord_type;
  if (_mode == COORD_DEFAULT) {
    _mode = xtal::COORD_MODE::CHECK();
  }
  xtal::COORD_MODE printer_mode(_mode);
  if (_mode != INTEGRAL) {
    // calculate nice widths
    int prec = this->opt.prec;
    int width = prec;
    Eigen::Vector3d vec;
    out.ostream().precision(prec);
    out.ostream().flags(std::ios::showpoint | std::ios::fixed |
                        std::ios::right);
    for (const auto &coord : clust) {
      if (_mode == CART)
        vec = coord.coordinate(clust.prim()).cart();
      else if (_mode == FRAC)
        vec = coord.coordinate(clust.prim()).frac();
      width = print_matrix_width(out, vec.transpose(), width);
    }

    // calculate nice widths
    Eigen::IOFormat format(prec, width + 1);
    for (const auto &coord : clust) {
      out << out.indent_str();
      coord.site(clust.prim()).print(out, format);
      if (this->opt.delim) out << this->opt.delim;
      out << std::flush;
    }
  } else {
    // calculate nice widths
    int prec = 1;
    int width = prec;
    out.ostream().flags(std::ios::showpoint | std::ios::fixed |
                        std::ios::right);
    for (const auto &coord : clust) {
      width = print_matrix_width(out, coord.unitcell().transpose(), width);
    }

    // print
    Eigen::IOFormat format(prec, width);
    for (const auto &coord : clust) {
      out << out.indent_str() << coord << " ";
      xtal::Site::print_occupant_dof(coord.site(clust.prim()).occupant_dof(),
                                     out);

      out << std::flush;
      if (this->opt.delim) out << this->opt.delim;
      out << std::flush;
    }
  }
}

// // explicit template instantiations
//
// #define PRINT_CLUST_INST(ITERATOR, INSERTER, PRINTER) \
//   template void print_clust<ITERATOR, PRINTER>(ITERATOR begin, ITERATOR end,
//   \
//                                                Log & out, PRINTER printer); \
//   template jsonParser &write_clust<ITERATOR>( \
//       ITERATOR begin, ITERATOR end, jsonParser & json, PRINTER printer); \
//   template jsonParser &write_clust<ITERATOR>( \
//       ITERATOR begin, ITERATOR end, jsonParser & json, PRINTER printer, \
//       const jsonParser &bspecs);
//
// #define ORBIT_CONTAINER_INST(ITERATOR, INSERTER, ORBIT) \
//   PRINT_CLUST_INST(ITERATOR, INSERTER, ProtoSitesPrinter) \
//   PRINT_CLUST_INST(ITERATOR, INSERTER, FullSitesPrinter) \
//   PRINT_CLUST_INST(ITERATOR, INSERTER, ProtoFuncsPrinter) \
//   template void print_clust<ITERATOR>(ITERATOR begin, ITERATOR end, Log &
//   out, \
//                                       const OrbitPrinterOptions &opt); \
//   template INSERTER read_clust<INSERTER, typename ORBIT::SymCompareType>( \
//       INSERTER result, const jsonParser &json, const Structure &prim, \
//       const SymGroup &generating_grp, \ const typename ORBIT::SymCompareType
//       &sym_compare);
//
// #define _VECTOR_IT(ORBIT) std::vector<ORBIT>::iterator
// #define _VECTOR_INSERTER(ORBIT) std::back_insert_iterator<std::vector<ORBIT>
// >
//
// #define _SET_IT(ORBIT) std::set<ORBIT>::iterator
// #define _SET_INSERTER(ORBIT) std::insert_iterator<std::set<ORBIT> >
//
// #define ORBIT_VECTOR_INST(ORBIT) \
//   ORBIT_CONTAINER_INST(_VECTOR_IT(ORBIT), _VECTOR_INSERTER(ORBIT), ORBIT)
// #define ORBIT_SET_INST(ORBIT) \
//   ORBIT_CONTAINER_INST(_SET_IT(ORBIT), _SET_INSERTER(ORBIT), ORBIT)
//
// ORBIT_VECTOR_INST(LocalIntegralClusterOrbit)
// ORBIT_VECTOR_INST(PrimPeriodicIntegralClusterOrbit)
// ORBIT_VECTOR_INST(ScelPeriodicIntegralClusterOrbit)
//
// ORBIT_SET_INST(LocalIntegralClusterOrbit)
// ORBIT_SET_INST(PrimPeriodicIntegralClusterOrbit)
// ORBIT_SET_INST(ScelPeriodicIntegralClusterOrbit)

}  // namespace CASM
