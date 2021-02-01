#include "casm/clex/io/ProtoFuncsPrinter_impl.hh"
#include "casm/clusterography/ClusterOrbits_impl.hh"

namespace CASM {

ProtoFuncsPrinter::ProtoFuncsPrinter(ClexBasis const &_clex_basis,
                                     PrimType_ptr _prim_ptr,
                                     OrbitPrinterOptions const &_opt)
    : SitesPrinter(_opt), clex_basis(_clex_basis), prim_ptr(_prim_ptr) {
  for (auto const &dofset : clex_basis.site_bases()) {
    for (BasisSet const &bset : dofset.second) {
      if (dofset.first != "occ" && bset.size() &&
          bset[0]->type_name() != "Variable") {
        labelers.push_back(SubExpressionLabeler(
            bset.name(), "\\phi^{(" + dofset.first + ")}_%n_%l"));
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
  std::vector<Orbit<PrimPeriodicSymCompare<IntegralCluster> > > asym_unit;
  make_prim_periodic_asymmetric_unit(prim_ptr,
                                     CASM_TMP::ConstantFunctor<bool>(true), TOL,
                                     std::back_inserter(asym_unit), nullstream);

  const Structure &prim = *prim_ptr;

  for (auto const &dofset : clex_basis.site_bases()) {
    out << indent << indent << "Site basis functions for DoF \"" << dofset.first
        << "\":\n";
    for (Index no = 0; no < asym_unit.size(); no++) {
      out << indent << indent << "Asymmetric unit " << no + 1 << ":\n";
      for (Index ne = 0; ne < asym_unit[no].size(); ne++) {
        Index b = asym_unit[no][ne][0].sublattice();
        out << indent << indent << "  Basis site " << b << ":\n"
            << "  ";
        if (printer_mode.check() == INTEGRAL) {
          out << indent << indent << asym_unit[no][ne][0] << ' ';
          Site::print_occupant_dof(
              asym_unit[no][ne][0].site(prim).occupant_dof(), out);
          out << std::flush;
        } else
          asym_unit[no][ne][0].site(prim).print(out);

        out << "\n";
        if (dofset.second[b].size() == 0)
          out << "        [No site basis functions]\n\n";
        if (dofset.first == "occ") {
          for (Index f = 0; f < dofset.second[b].size(); f++) {
            BasisSet tbasis(dofset.second[b]);

            int s;
            /* std::vector<DoF::RemoteHandle> remote(1, DoF::RemoteHandle("occ",
             * "s", prim.basis()[b].occupant_dof().ID())); */
            std::vector<DoF::RemoteHandle> remote(
                1, DoF::RemoteHandle("occ", "s", b));
            remote[0] = s;
            tbasis.register_remotes(remote);

            for (s = 0; s < prim.basis()[b].occupant_dof().size(); s++) {
              if (s == 0) out << "    ";
              out << "    \\phi_" << b << '_' << f << '['
                  << prim.basis()[b].occupant_dof()[s].name()
                  << "] = " << tbasis[f]->remote_eval();
              if (s + 1 == prim.basis()[b].occupant_dof().size())
                out << "\n";
              else
                out << ",   ";
            }
          }
        } else {
          // std::string formula;
          // bool relabel=false;
          for (Index f = 0; f < dofset.second[b].size(); f++) {
            if (dofset.second[b][f]) {
              out << "        ";
              if (dofset.second[b][f]->type_name() != "Variable") {
                out << "\\phi^{" << (dofset.first) << "}_" << b << '_' << f
                    << " = ";
              }
              out << dofset.second[b][f]->tex_formula() << "\n";
            }
          }
        }
      }
    }
  }
  out << "\n\n";
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
  std::vector<Orbit<PrimPeriodicSymCompare<IntegralCluster> > > asym_unit;
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
            fname << "\\phi_" << b << '_' << f;
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
                       tbasis[f]->remote_eval();
            }
          } else {
            if (dofset.second[b][f]) {
              if (dofset.second[b][f]->type_name() != "Variable") {
                fname << "\\phi^{" << (dofset.first) << "}_" << b << '_' << f;
              } else
                fname << "var_" << b << '_' << f;
            }
            sitef[b][dofset.first]["basis"][fname.str()] =
                dofset.second[b][f]->tex_formula();
          }
        }
      }
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
