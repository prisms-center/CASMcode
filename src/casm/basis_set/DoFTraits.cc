#include "casm/basis_set/DoFTraits.hh"
#include "casm/basis_set/OccupationDoFTraits.hh"
#include "casm/basis_set/StrainDoFTraits.hh"
#include "casm/basis_set/DisplacementDoFTraits.hh"
#include "casm/basis_set/MagSpinDoFTraits.hh"
#include "casm/basis_set/FunctionVisitor.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/symmetry/Orbit_impl.hh"
#include "casm/clusterography/IntegralCluster.hh"
#include "casm/clusterography/ClusterSymCompare_impl.hh"
#include "casm/clusterography/ClusterOrbits_impl.hh"
#include "casm/clex/ClexBasis.hh"
#include "casm/clex/NeighborList.hh"

namespace CASM {

  namespace DoFType {

    //************************************************************
    std::string Traits::clexulator_point_prepare_string(Structure const &_prim,
                                                        std::map<UnitCellCoord, std::set<UnitCellCoord> > const &_nhood,
                                                        PrimNeighborList &_nlist,
                                                        std::vector<BasisSet> const &site_bases,
                                                        std::string const &indent) const {

      std::stringstream ss;
      ss << indent << "switch(nlist_ind) {\n";
      for(auto const &nbor : _nhood) {
        std::stringstream ssvar, ssfunc;
        ss << indent << "case " << _nlist.neighbor_index(nbor.first) << ":\n";
        //Index n = nbor.first;
        std::cout << "neighborhood of nbor.first: " << nbor.first << ": \n";

        //Put neighborhood in a sensible order:
        std::map<Index, std::set<Index> > sublat_nhood;
        for(auto const &ucc : nbor.second) {
          sublat_nhood[ucc.sublat()].insert(_nlist.neighbor_index(ucc));
          std::cout << "ucc : " << ucc << "; n: " << _nlist.neighbor_index(ucc)  << "\n";
        }

        for(auto const &sublat : sublat_nhood) {
          Index b = sublat.first;
          for(Index n : sublat.second) {

            if(!_prim.basis()[b].has_dof(type_name()))
              continue;

            for(Index a = 0; _prim.basis()[b].dof(type_name()).size(); ++a) {
              ssvar << indent << "    m_params.write(m_" << type_name() << "_var_param_key, " << n
                    << ", " << a << ", eval_" << type_name() << "_var_" << b << "_" << a << "(" << n << "));\n";
            }


            if(requires_site_basis()) {
              for(Index f = 0; f < site_bases[b].size(); f++) {
                ssfunc << indent << "    m_params.write(m_" << type_name() << "_func_param_key, " << n
                       << ", " << f << ", eval_" << type_name() << "_func_" << b << "_" << f << "(" << n << "));\n";
              }
            }
          }
        }

        ss <<
           indent << "  if(m_params.eval_mode(m_" << type_name() << "_var_param_key) == ParamPack::DEFAULT) {\n" <<
           ssvar.str() <<
           indent << "  }\n";

        if(requires_site_basis()) {
          ss <<
             indent << "  if(m_params.eval_mode(m_" << type_name() << "_func_param_key) == ParamPack::DEFAULT) {\n" <<
             ssfunc.str() <<
             indent << "  }\n";
        }
        ss << indent << "break;\n";
      }
      ss << indent << "}\n";
      return ss.str();
    }

    //************************************************************

    std::string Traits::clexulator_global_prepare_string(Structure const &_prim,
                                                         std::map<UnitCellCoord, std::set<UnitCellCoord> > const &_nhood,
                                                         PrimNeighborList &_nlist,
                                                         std::vector<BasisSet> const &site_bases,
                                                         std::string const &indent) const {
      std::stringstream ss;

      std::map<Index, std::set<Index> > tot_nhood;
      for(auto const &nbor : _nhood)
        for(auto const &ucc : nbor.second)
          tot_nhood[ucc.sublat()].insert(_nlist.neighbor_index(ucc));

      std::stringstream ssvar, ssfunc;

      for(auto const &nbor : tot_nhood) {
        Index b = nbor.first;
        for(Index n : nbor.second) {

          if(!_prim.basis()[b].has_dof(type_name()))
            continue;

          for(Index a = 0; _prim.basis()[b].dof(type_name()).size(); ++a) {
            ssvar << indent << "    m_params.write(m_" << type_name() << "_var_param_key, " << n
                  << ", " << a << ", eval_" << type_name() << "_var_" << b << "_" << a << "(" << n << "));\n";
          }


          if(requires_site_basis()) {
            for(Index f = 0; f < site_bases[b].size(); f++) {
              ssfunc << indent << "    m_params.write(m_" << type_name() << "_func_param_key, " << n
                     << ", " << f << ", eval_" << type_name() << "_func_" << b << "_" << f << "(" << n << "));\n";
            }
          }
        }
      }
      ss <<
         indent << "if(m_params.eval_mode(m_" << type_name() << "_var_param_key) == ParamPack::DEFAULT) {\n" <<
         ssvar.str() <<
         indent << "}\n";
      if(requires_site_basis()) {
        ss <<
           indent << "if(m_params.eval_mode(m_" << type_name() << "_func_param_key) == ParamPack::DEFAULT) {\n" <<
           ssfunc.str() <<
           indent << "}\n";
      }

      return ss.str();
    }
    //************************************************************

    std::string Traits::clexulator_member_declarations_string(Structure const &_prim,
                                                              std::vector<BasisSet> const &_site_bases,
                                                              std::string const &indent) const {
      std::stringstream stream;
      std::vector<Orbit<IntegralCluster, PrimPeriodicSymCompare<IntegralCluster> > > asym_unit;
      std::ostream nullstream(0);
      make_prim_periodic_asymmetric_unit(_prim,
                                         CASM_TMP::ConstantFunctor<bool>(true),
                                         TOL,
                                         std::back_inserter(asym_unit),
                                         nullstream);
      /*
      for(Index no = 0; no < asym_unit.size(); no++) {
        Index nb = asym_unit[no][0][0].sublat();
        if(_site_bases[nb].size() == 0)
          continue;
        stream <<
               indent << "// Occupation Function tables for basis sites in asymmetric unit " << no << ":\n";
        for(Index ne = 0; ne < asym_unit[no].size(); ne++) {
          nb = asym_unit[no][ne][0].sublat();
          stream <<
                 indent << "//   - basis site " << nb << ":\n";
          for(Index f = 0; f < _site_bases[nb].size(); f++) {
            stream <<
                   indent << "double " << "m_" << type_name() << "_func_" << nb << '_' << f << '[' << _prim.basis()[nb].site_occupant().size() << "];\n";
          }
          stream << '\n';
        }

        }*/
      return stream.str();
    }

    //************************************************************

    std::string Traits::clexulator_private_method_definitions_string(Structure const &_prim,
                                                                     std::vector<BasisSet> const &_site_bases,
                                                                     const std::string &indent) const {
      return "";
    }
    //************************************************************

    std::string Traits::clexulator_private_method_declarations_string(Structure const &_prim,
                                                                      std::vector<BasisSet> const &_site_bases,
                                                                      const std::string &indent) const {
      std::stringstream stream;
      std::vector<Orbit<IntegralCluster, PrimPeriodicSymCompare<IntegralCluster> > > asym_unit;
      std::ostream nullstream(0);
      make_prim_periodic_asymmetric_unit(_prim,
                                         CASM_TMP::ConstantFunctor<bool>(true),
                                         TOL,
                                         std::back_inserter(asym_unit),
                                         nullstream);


      for(Index no = 0; no < asym_unit.size(); no++) {
        Index nb = asym_unit[no][0][0].sublat();
        if(_site_bases[nb].size() == 0)
          continue;


        for(Index ne = 0; ne < asym_unit[no].size(); ne++) {
          nb = asym_unit[no][ne][0].sublat();
          if(!_prim.basis()[nb].has_dof(type_name()))
            continue;
          stream <<
                 indent << "// " << type_name() << " evaluators and accessors for basis site " << nb << ":\n";
          for(Index a = 0; _prim.basis()[nb].dof(type_name()).size(); ++a) {
            stream <<
                   indent << "double const &eval_" << type_name() << "_var_" << nb << '_' << a << "(const int &nlist_ind) const {\n" <<
                   indent << "  return m_local_dof_ptrs[m_" << type_name() << "_var_param_key.index()].site_value(_l(nlist_ind))[" << a << "];\n" <<
                   indent << "}\n\n" <<

                   indent << "double const &" << type_name() << "_var_" << nb << '_' << a << "(const int &nlist_ind) const {\n" <<
                   indent << "  return " << "m_params.read(m_" << type_name() << "_var_param_key, nlist_ind, " << a << ");\n" <<
                   indent << "}\n";
          }

          if(requires_site_basis()) {
            for(Index f = 0; f < _site_bases[nb].size(); f++) {
              stream <<
                     indent << "double const &eval_" << type_name() << "_func_" << nb << '_' << f << "(const int &nlist_ind) const;// {\n" <<
                     //indent << "  return " << "m_" << type_name() << "_func_" << nb << '_' << f << "[_configdof().occ(_l(nlist_ind))];\n" <<
                     //indent << "}\n\n" <<

                     indent << "double const &" << type_name() << "_func_" << nb << '_' << f << "(const int &nlist_ind) const {\n" <<
                     indent << "  return " << "m_params.read(m_" << type_name() << "_func_param_key, nlist_ind, " << f << ");\n" <<
                     indent << "}\n";
            }
            stream << '\n';
          }
        }

      }
      return stream.str();
    }

    //************************************************************

    std::vector<std::tuple<std::string, Index, Index> > Traits::param_pack_allocation(Structure const &_prim,
        std::vector<BasisSet> const &_bases) const {
      std::vector<std::tuple<std::string, Index, Index> > result;
      Index NB = 0, NV = 0;
      bool basis_allocation = false;
      for(BasisSet const &basis : _bases) {
        NB = max(basis.size(), NB);
        for(Index f = 0; f < basis.size() && !basis_allocation; ++f) {
          if(basis[f] && basis[f]->type_name() != "Variable")
            basis_allocation = true;
        }
      }


      for(Site const &site : _prim.basis())
        NV = max(NV, site.dof(type_name()).size());

      //for(Index i = 0; i < NB; i++)
      result.push_back(std::make_tuple(std::string(type_name() + "_var"), Index(-1), Index(NV)));

      if(basis_allocation)
        result.push_back(std::make_tuple(std::string(type_name() + "_func"), Index(-1), Index(NB)));

      return result;

    }

    //************************************************************

    std::string Traits::clexulator_constructor_string(Structure const &_prim,
                                                      std::vector<BasisSet> const &_site_bases,
                                                      const std::string &indent) const {
      std::stringstream stream;
      stream.flags(std::ios::showpoint | std::ios::fixed | std::ios::right);
      stream.precision(10);

      std::vector<Orbit<IntegralCluster, PrimPeriodicSymCompare<IntegralCluster> > > asym_unit;
      std::ostream nullstream(0);
      make_prim_periodic_asymmetric_unit(_prim,
                                         CASM_TMP::ConstantFunctor<bool>(true),
                                         TOL,
                                         std::back_inserter(asym_unit),
                                         nullstream);
      /*
      for(const auto &asym : asym_unit) {
        for(const auto &equiv : asym) {
          Index nb = equiv[0].sublat();
          for(Index f = 0; f < _site_bases[nb].size(); f++) {

            for(Index s = 0; s < _prim.basis()[nb].site_occupant().size(); s++) {
              OccFuncEvaluator t_eval(s);
              _site_bases[nb][f]->accept(t_eval);

              if(s == 0)
                stream << indent;
              stream << "m_" << type_name() << "_func_" << nb << '_' << f << '[' << s << "] = "
                     << t_eval.value();
              if(s + 1 == _prim.basis()[nb].site_occupant().size())
                stream << ";\n\n";
              else
                stream << ", ";
            }
          }
        }
        }*/
      return stream.str();
    }


    std::vector<std::unique_ptr<FunctionVisitor> > Traits::site_function_visitors() const {
      std::vector<std::unique_ptr<FunctionVisitor> > result;
      result.push_back(std::unique_ptr<FunctionVisitor>(new VariableLabeler("%t_%s(%n)")));
      return result;
    }

    std::vector<std::unique_ptr<FunctionVisitor> > Traits::clust_function_visitors() const {
      std::vector<std::unique_ptr<FunctionVisitor> > result;
      result.push_back(std::unique_ptr<FunctionVisitor>(new VariableLabeler("%t_%s(%n)")));
      return result;
    }
  }

  template<>
  DoFType::TraitsDictionary make_parsing_dictionary<DoF::BasicTraits>() {
    //std::cout << "Making Parsing dictionary... \n";
    DoF::register_traits(DoFType::occupation());
    DoF::register_traits(DoFType::displacement());
    DoF::register_traits(DoFType::magspin());
    DoF::register_traits(DoFType::EAstrain());
    DoF::register_traits(DoFType::Hstrain());
    DoF::register_traits(DoFType::GLstrain());
    DoFType::TraitsDictionary dict;

    dict.insert(
      DoFType::occupation(),
      DoFType::displacement(),
      DoFType::magspin(),
      DoFType::EAstrain(),
      DoFType::Hstrain(),
      DoFType::GLstrain());
    return dict;
  }
}
