#ifndef CASM_clex_io_ProtoFuncsPrinter_impl
#define CASM_clex_io_ProtoFuncsPrinter_impl

#include "casm/casm_io/json/jsonParser.hh"
#include "casm/clex/ClexBasis.hh"
#include "casm/clex/io/ProtoFuncsPrinter.hh"
#include "casm/clusterography/io/OrbitPrinter_impl.hh"

namespace CASM {

template <typename OrbitType>
void ProtoFuncsPrinter::operator()(const OrbitType &orbit, Log &out,
                                   Index orbit_index, Index Norbits) const {
  out << out.indent_str() << "Prototype"
      << " of " << orbit.size() << " Equivalent " << element_name
      << " in Orbit " << orbit_index << std::endl;

  // out.flags(std::ios::showpoint | std::ios::fixed | std::ios::left);
  // out.precision(5);

  COORD_MODE printer_mode(opt.coord_type);

  auto const &clust = orbit.prototype();
  Index np = 0;
  this->increase_indent(out);
  for (const auto &coord : clust) {
    out.indent();
    if (opt.coord_type == INTEGRAL) {
      out << coord;
      out << " ";
      Site::print_occupant_dof(coord.site(*prim_ptr).occupant_dof(), out);
      out << std::flush;
    } else {
      coord.site(*prim_ptr).print(out);
    }
    out << "  basis_index: " << coord.sublattice() << "  clust_index: " << np++
        << " ";
    if (opt.delim) out << opt.delim;
    out << std::flush;
  }
  this->decrease_indent(out);

  Index func_index = 0;
  for (Index i = 0; i < orbit_index; i++)
    func_index += clex_basis.clust_basis(i, 0).size();

  // From clust:
  out.indent() << "Basis Functions:\n";
  BasisSet tbasis(clex_basis.clust_basis(orbit_index, 0));
  tbasis.accept(OccFuncLabeler("\\phi_%b_%f(s_%n)"));
  this->increase_indent(out);
  for (Index i = 0; i < tbasis.size(); i++) {
    out.indent() << "\\Phi_" << func_index + i << " = "
                 << tbasis[i]->tex_formula() << std::endl;
  }
  this->decrease_indent(out);
  out << "\n\n" << std::flush;
}

template <typename OrbitType>
jsonParser &ProtoFuncsPrinter::to_json(const OrbitType &orbit, jsonParser &json,
                                       Index orbit_index, Index Norbits) const {
  json.put_obj();
  json["prototype"] = orbit.prototype();
  json["linear_orbit_index"] = orbit_index;
  json["mult"] = orbit.size();

  jsonParser &orbitf = json["cluster_functions"];
  orbitf = jsonParser::array();

  // basis function info
  Index func_index = 0;
  for (Index i = 0; i < orbit_index; i++)
    func_index += clex_basis.clust_basis(i, 0).size();

  BasisSet tbasis(clex_basis.clust_basis(orbit_index, 0));
  tbasis.accept(OccFuncLabeler("\\phi_%b_%f(s_%n)"));
  for (auto const &labeler : labelers) {
    tbasis.accept(labeler);
  }

  for (Index nf = 0; nf < tbasis.size(); ++nf) {
    orbitf.push_back(json_pair("\\Phi_" + std::to_string(func_index + nf),
                               tbasis[nf]->tex_formula()));
  }

  return json;
}

}  // namespace CASM

#endif
