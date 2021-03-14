#ifndef CASM_clex_io_ProtoFuncsPrinter_impl
#define CASM_clex_io_ProtoFuncsPrinter_impl

#include "casm/casm_io/json/jsonParser.hh"
#include "casm/clex/ClexBasis.hh"
#include "casm/clex/io/ProtoFuncsPrinter.hh"
#include "casm/crystallography/Site.hh"

namespace CASM {

template <typename OrbitType>
void ProtoFuncsPrinter::operator()(const OrbitType &orbit, Log &out,
                                   Index orbit_index, Index Norbits) const {
  if (m_align) {
    print_tex_tabular_cluster_sites(out, orbit.prototype(), *prim_ptr,
                                    opt.coord_type);

    out << "Prototype basis functions:\n\n";
    Index func_index = 0;
    for (Index i = 0; i < orbit_index; i++)
      func_index += clex_basis.clust_basis(i, 0).size();

    BasisSet tbasis(clex_basis.clust_basis(orbit_index, 0));
    tbasis.accept(OccFuncLabeler("\\phi_{%b,%f}(s_{%n})"));
    out << "\\begin{align*}\n";
    for (Index i = 0; i < tbasis.size(); i++) {
      out << "\\Phi_{" << func_index + i << "} ={}& ";

      // wrap formula across multiple lines if its long
      std::string formula = tbasis[i]->tex_formula();
      int maxl = opt.max_line_width;
      int j = 0;
      int pos = 0;
      while (formula.size() > maxl) {
        if (formula[j] == '-' || formula[j] == '+') {
          if (j > maxl) {
            if (pos == 0) {
              pos = j;
            }
            out << formula.substr(0, pos) << "\\\\\n";
            out << " & ";
            formula = formula.substr(pos, formula.size());
            j = 0;
            pos = 0;
          } else {
            pos = j;
          }
        }
        ++j;
      }
      out << formula;

      if (i + 1 != tbasis.size()) {
        out << "\\\\\n";
      } else {
        out << "\n";
      }
    }
    out << "\\end{align*}\n\n";
  } else {
    out << out.indent_str() << "Prototype"
        << " of " << orbit.size() << " Equivalent " << element_name
        << " in Orbit " << orbit_index << std::endl;

    // out.flags(std::ios::showpoint | std::ios::fixed | std::ios::left);
    // out.precision(5);

    xtal::COORD_MODE printer_mode(opt.coord_type);

    auto const &clust = orbit.prototype();
    this->increase_indent(out);
    if (this->opt.print_coordinates) {
      print_coordinates(*this, orbit.prototype(), out);
    }
    if (this->opt.print_invariant_group) {
      this->print_invariant_group(orbit, orbit.prototype(), out);
    }
    this->decrease_indent(out);

    Index func_index = 0;
    for (Index i = 0; i < orbit_index; i++)
      func_index += clex_basis.clust_basis(i, 0).size();

    // From clust:
    out.indent() << "Prototype basis functions:\n";
    BasisSet tbasis(clex_basis.clust_basis(orbit_index, 0));
    tbasis.accept(OccFuncLabeler("\\phi_{%b,%f}(s_{%n})"));
    this->increase_indent(out);
    for (Index i = 0; i < tbasis.size(); i++) {
      out.indent() << "\\Phi_{" << func_index + i
                   << "} = " << tbasis[i]->tex_formula() << std::endl;
    }
    this->decrease_indent(out);
    out << std::flush;
  }
}

template <typename OrbitType>
jsonParser &ProtoFuncsPrinter::to_json(const OrbitType &orbit, jsonParser &json,
                                       Index orbit_index, Index Norbits) const {
  json.put_obj();
  json["prototype"] = orbit.prototype();
  json["linear_orbit_index"] = orbit_index;
  json["mult"] = orbit.size();
  if (this->opt.print_invariant_group) {
    this->print_invariant_group(orbit, orbit.prototype(), json["prototype"]);
  }

  jsonParser &orbitf = json["cluster_functions"];
  orbitf = jsonParser::array();

  // basis function info
  Index func_index = 0;
  for (Index i = 0; i < orbit_index; i++)
    func_index += clex_basis.clust_basis(i, 0).size();

  BasisSet tbasis(clex_basis.clust_basis(orbit_index, 0));
  tbasis.accept(OccFuncLabeler("\\phi_{%b,%f}(s_{%n})"));
  for (auto const &labeler : labelers) {
    tbasis.accept(labeler);
  }

  for (Index nf = 0; nf < tbasis.size(); ++nf) {
    orbitf.push_back(
        json_pair("\\Phi_{" + std::to_string(func_index + nf) + "}",
                  tbasis[nf]->tex_formula()));
  }

  return json;
}

}  // namespace CASM

#endif
