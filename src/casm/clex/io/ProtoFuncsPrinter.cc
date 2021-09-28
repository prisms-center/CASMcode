#include "casm/basis_set/DoFTraits.hh"
#include "casm/clex/ClexBasis_impl.hh"
#include "casm/clex/io/ProtoFuncsPrinter_impl.hh"
#include "casm/clusterography/ClusterOrbits_impl.hh"
#include "casm/clusterography/ClusterSpecs_impl.hh"

namespace CASM {

namespace {
double pretty(double value, double tol) {
  double pretty_value = value;
  if (std::abs(value) < tol) {
    return 0.;
  }
  if (std::abs(std::round(value) - value) < tol) {
    pretty_value = std::round(value);
  }
  return pretty_value;
}
}  // namespace

ProtoFuncsPrinter::ProtoFuncsPrinter(ClexBasis const &_clex_basis,
                                     PrimType_ptr _prim_ptr, bool align,
                                     OrbitPrinterOptions const &_opt)
    : SitesPrinter(_opt),
      clex_basis(_clex_basis),
      prim_ptr(_prim_ptr),
      m_align(align) {
  for (auto const &dofset : clex_basis.site_bases()) {
    for (BasisSet const &bset : dofset.second) {
      if (dofset.first != "occ" && bset.size() &&
          bset[0]->type_name() != "Variable") {
        labelers.push_back(SubExpressionLabeler(
            bset.name(), "\\phi^{(" + dofset.first + ")}_{%n,%l}"));
      }
    }
  }
}

void print_site_basis_funcs(std::shared_ptr<const Structure> prim_ptr,
                            ClexBasis const &clex_basis, Log &out,
                            Index indent_space, COORD_TYPE mode) {
  std::string indent(indent_space, ' ');

  std::ostream nullstream(0);
  xtal::COORD_MODE printer_mode(mode);
  std::vector<Orbit<PrimPeriodicSymCompare<IntegralCluster>>> asym_unit;
  make_prim_periodic_asymmetric_unit(prim_ptr,
                                     CASM_TMP::ConstantFunctor<bool>(true), TOL,
                                     std::back_inserter(asym_unit), nullstream);

  const Structure &prim = *prim_ptr;

  for (auto const &dofset : clex_basis.site_bases()) {
    out << "Site basis functions for DoF \"" << dofset.first << "\":\n\n";
    for (Index no = 0; no < asym_unit.size(); no++) {
      out << "Asymmetric unit " << no + 1 << ":\n";
      for (Index ne = 0; ne < asym_unit[no].size(); ne++) {
        Index b = asym_unit[no][ne][0].sublattice();
        out << indent << "Basis site " << b << ":  ";
        if (printer_mode.check() == INTEGRAL) {
          out << asym_unit[no][ne][0] << ' ';
          Site::print_occupant_dof(
              asym_unit[no][ne][0].site(prim).occupant_dof(), out);
          out << std::endl;
        } else {
          asym_unit[no][ne][0].site(prim).print(out);
        }

        if (dofset.second[b].size() == 0)
          out << indent << "[No site basis functions]\n\n";
        if (dofset.first == "occ") {
          Index max_length = 0;
          std::vector<std::vector<std::string>> funcstrs;
          for (Index f = 0; f < dofset.second[b].size(); f++) {
            funcstrs.push_back(std::vector<std::string>());
            BasisSet tbasis(dofset.second[b]);

            int s;
            /* std::vector<DoF::RemoteHandle> remote(1, DoF::RemoteHandle("occ",
             * "s", prim.basis()[b].occupant_dof().ID())); */
            std::vector<DoF::RemoteHandle> remote(
                1, DoF::RemoteHandle("occ", "s", b));
            remote[0] = s;
            tbasis.register_remotes(remote);

            for (s = 0; s < prim.basis()[b].occupant_dof().size(); s++) {
              std::stringstream ss;
              ss << "\\phi_{" << b << "," << f << "}["
                 << prim.basis()[b].occupant_dof()[s].name()
                 << "] = " << pretty(tbasis[f]->remote_eval(), 1e-10);
              if (s + 1 != prim.basis()[b].occupant_dof().size()) ss << ",";
              funcstrs[f].push_back(ss.str());
              if (ss.str().size() > max_length) {
                max_length = ss.str().size();
              }
            }
          }
          for (auto const &sublat_funcstrs : funcstrs) {
            out << indent << indent;
            for (auto const &funcstr : sublat_funcstrs) {
              out << std::setw(max_length + 3) << std::left << funcstr;
            }
            out << std::endl;
          }

        } else {
          // std::string formula;
          // bool relabel=false;
          for (Index f = 0; f < dofset.second[b].size(); f++) {
            if (dofset.second[b][f]) {
              out << indent << indent;
              if (dofset.second[b][f]->type_name() != "Variable") {
                out << "\\phi^{" << (dofset.first) << "}_{" << b << "," << f
                    << "} = ";
              }
              out << dofset.second[b][f]->tex_formula() << "\n";
            }
          }
        }
      }
    }

    out << std::endl << std::endl;
  }
}

/// Print prototype cluster sites as a tex tabular
void print_tex_tabular_cluster_sites(Log &out, IntegralCluster const &cluster,
                                     xtal::BasicStructure const &prim,
                                     COORD_TYPE mode) {
  if (!cluster.size()) {
    return;
  }

  if (mode == FRAC) {
    out << "Prototype cluster sites (Fractional coordinates):\n\n";
  } else if (mode == CART) {
    out << "Prototype cluster sites (Cartesian coordinates):\n\n";
  } else if (mode == INTEGRAL) {
    out << "Prototype cluster sites (Integral coordinates):\n\n";
  }

  // for sites in asym_unit:
  // 'tabular: site index & basis index  & u & v & w & occupants \\'
  out << "\\vspace{2mm}\n"
      << "\\begin{tabular}{rrrrrr}\n"
      << "\\hline\n";
  if (mode == FRAC) {
    out << "site index & basis index  & u & v & w & occupants \\\\\n";
  } else if (mode == CART) {
    out << "site index & basis index  & x & y & z & occupants \\\\\n";
  } else if (mode == INTEGRAL) {
    out << "site index & basis index  & i & j & k & occupants \\\\\n";
  }
  out << "\\hline\n";
  //       `0 & 0 & 0.0000000& 0.0000000& 0.0000000&  A B C \\`
  //       `1 & 2 & 0.0000000& 0.0000000& 0.0000000&  A B C \\`
  for (Index site_index = 0; site_index < cluster.size(); site_index++) {
    Index b = cluster[site_index].sublattice();
    Site const &site = cluster[site_index].site(prim);

    std::stringstream occupants_ss;
    xtal::Site::print_occupant_dof(site.occupant_dof(), occupants_ss);
    std::string occupants = occupants_ss.str();

    if (mode == FRAC) {
      Eigen::Vector3d frac = site.frac();
      out << site_index << " & " << b << " & " << frac(0) << " & " << frac(1)
          << " & " << frac(2) << " & " << occupants << " \\\\\n";
    } else if (mode == CART) {
      Eigen::Vector3d cart = site.cart();
      out << site_index << " & " << b << " & " << cart(0) << " & " << cart(1)
          << " & " << cart(2) << " & " << occupants << " \\\\\n";
    } else if (mode == INTEGRAL) {
      UnitCell ijk = cluster[site_index].unitcell();
      out << site_index << " & " << b << " & " << ijk(0) << " & " << ijk(1)
          << " & " << ijk(2) << " & " << occupants << " \\\\\n";
    }
  }
  out << "\\hline\n"
      << "\\end{tabular}\n"
      << "\\vspace{2mm}\n\n";
}

namespace {

void print_tex_tabular_asym_unit(
    Log &out, Orbit<PrimPeriodicSymCompare<IntegralCluster>> const &orbit,
    Structure const &prim, COORD_TYPE mode) {
  if (mode == FRAC) {
    out << "Sites (Fractional coordinates):\n\n";
  } else if (mode == CART) {
    out << "Sites (Cartesian coordinates):\n\n";
  } else if (mode == INTEGRAL) {
    out << "Sites (Integral coordinates):\n\n";
  }

  // for sites in asym_unit:
  // 'tabular: basis index  & u & v & w & occupants \\'
  out << "\\vspace{2mm}\n"
      << "\\begin{tabular}{rrrrr}\n"
      << "\\hline\n";
  if (mode == FRAC) {
    out << "basis index  & u & v & w & occupants \\\\\n";
  } else if (mode == CART) {
    out << "basis index  & x & y & z & occupants \\\\\n";
  } else if (mode == INTEGRAL) {
    out << "basis index  & i & j & k & occupants \\\\\n";
  }
  out << "\\hline\n";
  //       `0 & 0.0000000& 0.0000000& 0.0000000&  A B C \\`
  //       `2 & 0.0000000& 0.0000000& 0.0000000&  A B C \\`
  for (Index ne = 0; ne < orbit.size(); ne++) {
    Index b = orbit[ne][0].sublattice();
    Site const &site = orbit[ne][0].site(prim);

    std::stringstream occupants_ss;
    xtal::Site::print_occupant_dof(site.occupant_dof(), occupants_ss);
    std::string occupants = occupants_ss.str();

    if (mode == FRAC) {
      Eigen::Vector3d frac = site.frac();
      out << b << " & " << frac(0) << " & " << frac(1) << " & " << frac(2)
          << " & " << occupants << " \\\\\n";
    } else if (mode == CART) {
      Eigen::Vector3d cart = site.cart();
      out << b << " & " << cart(0) << " & " << cart(1) << " & " << cart(2)
          << " & " << occupants << " \\\\\n";
    } else if (mode == INTEGRAL) {
      UnitCell ijk = orbit[ne][0].unitcell();
      out << b << " & " << ijk(0) << " & " << ijk(1) << " & " << ijk(2) << " & "
          << occupants << " \\\\\n";
    }
  }
  out << "\\hline\n"
      << "\\end{tabular}\n"
      << "\\vspace{2mm}\n\n";
}

void print_tex_tabular_basis(Log &out, xtal::DoFSet const &dofset,
                             BasisSet const &basis_set) {
  Index dim = dofset.dim();
  Index standard_dim = dofset.basis().rows();

  // DoFSet basis:

  // \vspace{2mm}
  // \begin{tabular}{rccrrrc}
  //              & &   &        dx &       dy &        dz &   \\
  //           da &=& [ &  0.707107 & 0.707107 &       0.0 & ] \\
  //           db &=& [ & -0.707107 & 0.707107 &       0.0 & ]
  // \end{tabular}
  // \vspace{2mm}

  out << "\\vspace{2mm}\n"
      << "\\begin{tabular}{rcc" << std::string(standard_dim, 'r') << "c}\n";

  // header line
  out << " axis name & & & ";
  for (Index i = 0; i < standard_dim; ++i) {
    out << "$" << dofset.traits().standard_var_names()[i] << "$ & ";
  }
  out << " \\\\\n";

  for (Index i = 0; i < dim; i++) {
    out << "$" << dofset.component_names()[i] << "$ &=& [ & ";
    // basis vectors
    for (int s = 0; s < standard_dim; s++) {
      out << pretty(dofset.basis()(s, i), 1e-10) << " & ";
    }
    out << " ] ";
    if (i + 1 != dim) {
      out << "\\\\\n";
    } else {
      out << "\n";
    }
  }

  out << "\\end{tabular}\n"
      << "\\vspace{2mm}\n\n";
}

void print_tex_tabular_occ_basis_funcs(Log &out, Site const &site,
                                       BasisSet const &basis_set, Index b) {
  Index occ_dof_size = site.occupant_dof().size();

  // \vspace{2mm}
  // \begin{tabular}{rccrrrc}
  //              & &   &         A &       B &         C &   \\
  // \phi_{2}_{0} &=& [ &  -1.22474 &       0 &   1.22474 & ] \\
  // \phi_{2}_{1} &=& [ & -0.707107 & 1.41421 & -0.707107 & ]
  // \end{tabular}
  // \vspace{2mm}

  out << "\\vspace{2mm}\n"
      << "\\begin{tabular}{rcc" << std::string(occ_dof_size, 'r') << "c}\n";

  // header line
  out << " & & & ";
  for (Index i = 0; i < occ_dof_size; ++i) {
    out << site.occupant_dof()[i].name() << " & ";
  }
  out << " \\\\\n";

  // function (f) lines:
  for (Index f = 0; f < basis_set.size(); f++) {
    out << "$\\phi_{" << b << "," << f << "}$ &=& [ & ";

    // this gets basis function values:
    BasisSet tbasis(basis_set);
    int s;
    std::vector<DoF::RemoteHandle> remote(1, DoF::RemoteHandle("occ", "s", b));
    remote[0] = s;
    tbasis.register_remotes(remote);
    for (s = 0; s < occ_dof_size; s++) {
      out << pretty(tbasis[f]->remote_eval(), 1e-10) << " & ";
    }

    out << " ] ";
    if (f + 1 != basis_set.size()) {
      out << "\\\\\n";
    } else {
      out << "\n";
    }
  }
  out << "\\end{tabular}\n"
      << "\\vspace{2mm}\n\n";
}

void print_site_basis_set(Log &out, std::string dof_key,
                          BasisSet const &basis_set, Index sublattice_index) {
  Index dim = basis_set.size();

  Index n_functions = 0;
  for (Index f = 0; f < dim; f++) {
    if (basis_set[f]->type_name() != "Variable") {
      n_functions++;
    }
  }
  if (n_functions) {
    // continuous basis functions:
    out << "\\begin{align*}\n";
    // function (f) lines: (basis vectors)
    for (Index f = 0; f < dim; f++) {
      if (basis_set[f]) {
        if (basis_set[f]->type_name() != "Variable") {
          out << "\\phi^{" << dof_key << "}_{" << sublattice_index << "," << f
              << "} = ";
        }
        out << basis_set[f]->tex_formula() << "\n";
      }
      if (f + 1 != dim) {
        out << "\\\\\n";
      }
    }
    out << "\\end{align*}\n\n";
  }
}

void print_global_basis_set(Log &out, std::string dof_key,
                            BasisSet const &basis_set) {
  Index dim = basis_set.size();

  Index n_functions = 0;
  for (Index f = 0; f < dim; f++) {
    if (basis_set[f]->type_name() != "Variable") {
      n_functions++;
    }
  }
  if (n_functions) {
    // continuous basis functions:
    out << "\\begin{align}\n";
    // function (f) lines: (basis vectors)
    for (Index f = 0; f < dim; f++) {
      if (basis_set[f]) {
        if (basis_set[f]->type_name() != "Variable") {
          out << "\\phi^{" << dof_key << "}_{" << f << "} = ";
        }
        out << basis_set[f]->tex_formula() << "\n";
      }
      if (f + 1 != dim) {
        out << "\\\\\n";
      }
    }
    out << "\\end{align}\n\n";
  }
}

}  // namespace

/// Print occ site basis functions, and continuous site and global bases
void print_aligned_site_basis_funcs(std::shared_ptr<const Structure> prim_ptr,
                                    ClexBasis const &clex_basis, Log &out,
                                    Index indent_space, COORD_TYPE mode) {
  std::string indent(indent_space, ' ');

  std::ostream nullstream(0);
  xtal::COORD_MODE printer_mode(mode);
  std::vector<Orbit<PrimPeriodicSymCompare<IntegralCluster>>> asym_unit;
  make_prim_periodic_asymmetric_unit(prim_ptr,
                                     CASM_TMP::ConstantFunctor<bool>(true), TOL,
                                     std::back_inserter(asym_unit), nullstream);

  const Structure &prim = *prim_ptr;
  for (auto const &dofset : clex_basis.site_bases()) {
    out << "\\pagebreak\n";
    if (dofset.first == "occ") {
      out << "Site basis functions for DoF ``" << dofset.first << "\":\n\n";
    } else {
      out << "Site bases for DoF ``" << dofset.first << "\":\n\n";
    }
    out << "\\begin{itemize}\n";

    for (Index no = 0; no < asym_unit.size(); no++) {
      out << "\\item Asymmetric unit " << no + 1 << ":\n"
          << "\n";
      print_tex_tabular_asym_unit(out, asym_unit[no], prim, mode);

      // for basis site in asym_unit:
      for (Index ne = 0; ne < asym_unit[no].size(); ne++) {
        Index b = asym_unit[no][ne][0].sublattice();
        Site const &site = asym_unit[no][ne][0].site(prim);

        if (dofset.second[b].size() == 0) {
          out << "Basis site " << b << ": [No site basis functions]\n\n";
        }
        if (dofset.first == "occ") {
          // occ site DoF basis:
          out << "Basis site " << b << ":\n\n";
          print_tex_tabular_occ_basis_funcs(out, site, dofset.second[b], b);
        } else {
          // continuous site DoF basis:

          CASM::xtal::SiteDoFSet const &current_dofset =
              site.dofs().find(dofset.first)->second;

          out << "Basis site " << b << ":\n\n";
          print_tex_tabular_basis(out, current_dofset, dofset.second[b]);
          // print_site_basis_set(out, dofset.first, dofset.second[b], b);

        }  // end continuous site DoF basis:
      }
    }

    out << "\\end{itemize}\n\n";
  }  // end local bases

  for (auto const &dofset : clex_basis.global_bases()) {
    out << "\\pagebreak\n";
    out << "Global basis for DoF ``" << dofset.first << "\":\n\n";

    // continuous global DoF basis:
    for (Index i = 0; i < dofset.second.size(); ++i) {
      if (dofset.second.size() > 1) {  // this shouldn't happen
        out << "Basis " << i << ":\n\n";
      }

      CASM::xtal::DoFSet const &current_dofset =
          prim.structure().global_dofs().find(dofset.first)->second;

      print_tex_tabular_basis(out, current_dofset, dofset.second[i]);
      // print_global_basis_set(out, dofset.first, dofset.second[i]);
    }
  }  // end global bases
}

void write_site_basis_funcs(std::shared_ptr<const Structure> prim_ptr,
                            ClexBasis const &clex_basis, jsonParser &json) {
  //   "site_functions":[
  //     {
  //       "asym_unit": X,
  //       "sublat": 2,
  //       "basis": {
  //         "phi_b_0": {"Va":0.0, "O":1.0},
  //         "phi_b_1": {"Va":0.0, "O":1.0}
  //       }
  //     },
  //     ...
  //   ],

  jsonParser &sitef = json["site_functions"];
  sitef = jsonParser::array(prim_ptr->basis().size(), jsonParser::object());

  std::ostream nullstream(0);
  std::vector<Orbit<PrimPeriodicSymCompare<IntegralCluster>>> asym_unit;
  make_prim_periodic_asymmetric_unit(prim_ptr,
                                     CASM_TMP::ConstantFunctor<bool>(true), TOL,
                                     std::back_inserter(asym_unit), nullstream);
  const Structure &prim = *prim_ptr;

  for (auto const &dofset : clex_basis.site_bases()) {
    for (Index no = 0; no < asym_unit.size(); no++) {
      for (Index ne = 0; ne < asym_unit[no].size(); ne++) {
        Index b = asym_unit[no][ne][0].sublattice();
        sitef[b]["sublat"] = b;
        sitef[b]["asym_unit"] = no;

        if (dofset.second[b].size() == 0) {
          sitef[b][dofset.first]["basis"].put_null();
          continue;
        }

        for (Index f = 0; f < dofset.second[b].size(); f++) {
          std::stringstream fname;

          if (dofset.first == "occ") {
            fname << "\\phi_{" << b << "," << f << "}";
            BasisSet tbasis(dofset.second[b]);
            int s;
            /* std::vector<DoF::RemoteHandle> remote(1, DoF::RemoteHandle("occ",
             * "s", prim.basis()[b].occupant_dof().ID())); */
            std::vector<DoF::RemoteHandle> remote(
                1, DoF::RemoteHandle("occ", "s", b));
            remote[0] = s;
            tbasis.register_remotes(remote);

            for (s = 0; s < prim.basis()[b].occupant_dof().size(); s++) {
              sitef[b][dofset.first]["basis"][fname.str()]
                   [prim.basis()[b].occupant_dof()[s].name()] =
                       pretty(tbasis[f]->remote_eval(), 1e-10);
            }
          } else {
            if (dofset.second[b][f]) {
              if (dofset.second[b][f]->type_name() != "Variable") {
                fname << "\\phi^{" << (dofset.first) << "}_{" << b << "," << f
                      << "}";
              } else
                fname << "var_{" << b << "}_{" << f << "}";
            }
            sitef[b][dofset.first]["basis"][fname.str()] =
                dofset.second[b][f]->tex_formula();
          }
        }
      }
    }
  }
}

namespace write_clex_basis_impl {

struct WriteClexBasisImpl {
  /// Required input & resulting `json` to be populated
  WriteClexBasisImpl(std::shared_ptr<Structure const> const &_shared_prim,
                     ClexBasisSpecs const &_basis_set_specs, jsonParser &_json)
      : shared_prim(_shared_prim),
        basis_set_specs(_basis_set_specs),
        json(_json) {}

  std::shared_ptr<Structure const> const &shared_prim;
  ClexBasisSpecs const &basis_set_specs;
  jsonParser &json;

  template <typename OrbitVecType>
  void operator()(OrbitVecType const &orbits) const {
    ParsingDictionary<DoFType::Traits> const *dof_dict =
        &DoFType::traits_dict();
    ClexBasis clex_basis{shared_prim, basis_set_specs, dof_dict};
    clex_basis.generate(orbits.begin(), orbits.end());
    write_clex_basis(clex_basis, orbits, json);
  }
};

}  // namespace write_clex_basis_impl

/// Write basis.json format JSON
///
/// - Does not write "bspecs"
jsonParser write_clex_basis(std::shared_ptr<const Structure> const &shared_prim,
                            ClexBasisSpecs const &clex_basis_specs,
                            jsonParser &json) {
  using namespace write_clex_basis_impl;
  WriteClexBasisImpl f{shared_prim, clex_basis_specs, json};
  for_all_orbits(*clex_basis_specs.cluster_specs, log(), f);
  return json;
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
