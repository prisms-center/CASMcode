#include "casm/basis_set/DoFTraits.hh"
#include "casm/basis_set/OccupationDoFTraits.hh"
#include "casm/basis_set/StrainDoFTraits.hh"
#include "casm/basis_set/DisplacementDoFTraits.hh"
#include "casm/basis_set/MagSpinDoFTraits.hh"
#include "casm/basis_set/FunctionVisitor.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/crystallography/SimpleStructure.hh"
#include "casm/symmetry/Orbit_impl.hh"
#include "casm/clusterography/IntegralCluster.hh"
#include "casm/clusterography/ClusterSymCompare_impl.hh"
#include "casm/clusterography/ClusterOrbits_impl.hh"
#include "casm/clex/ClexBasis.hh"
#include "casm/clex/ConfigDoF.hh"
#include "casm/clex/NeighborList.hh"
#include "casm/misc/ParsingDictionary.hh"
#include <memory>


namespace CASM {

  template<>
  DoFType::TraitsDictionary make_parsing_dictionary<DoFType::Traits>() {

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

  namespace DoFType {
    namespace Local {
      static TraitsDictionary &_traits_dict() {
        static TraitsDictionary static_dict = make_parsing_dictionary<Traits>();
        return static_dict;

      }
    }

    void register_traits(Traits const &_traits) {
      Local::_traits_dict().insert(_traits);
    }

    Traits const &traits(std::string const &dof_key) {
      return Local::_traits_dict().lookup(dof_key);
    }


    /// \brief Retrieve the standard values for a DoF from dictionary of properties from properties.calc.json
    ///  Returns matrix with standard values, and names of properties that were used to construct the matrix
    std::pair<Eigen::MatrixXd, std::set<std::string> > Traits::find_values(std::map<std::string, Eigen::MatrixXd> const &values) const {
      std::pair<Eigen::MatrixXd, std::set<std::string> > result;
      for(auto const &val : values) {
        if(AnisoValTraits::name_suffix(val.first) == this->name()) {
          result.first = val.second;
          result.second.insert(val.first);
          return result;
        }
      }
      throw std::runtime_error("Could not identify DoF values for DoF '" + (this->name()) + "' from provided list of tabulated structure properties.");
      return result;
    }

    void Traits::to_json(DoFSet const &_out, jsonParser &_json) const {
      bool simple = false;
      if(_out.dim() == val_traits().standard_var_names().size() && _out.basis().isIdentity()) {
        simple = true;
        for(Index i = 0; i < _out.dim(); ++i) {
          if(_out[i].var_name() != val_traits().standard_var_names()[i]) {
            simple = false;
            break;
          }
        }
      }

      if(simple)
        _json.put_obj();
      else {
        _json["axes"] = _out.basis().transpose();
        _json["axis_names"].put_array();
        for(Index i = 0; i < _out.dim(); ++i) {
          _json["axis_names"].push_back(_out[i].var_name());
        }
      }
    }

    //************************************************************
    void Traits::apply_dof(ConfigDoF const &_dof, BasicStructure<Site> const &_reference, SimpleStructure &_struc) const {
      return;
    }

    //************************************************************
    jsonParser Traits::dof_to_json(ConfigDoF const &_dof, BasicStructure<Site> const &_reference) const {
      jsonParser result;
      result.put_obj();

      if(val_traits().global())
        result["value"] = _dof.global_dof(name()).standard_values();
      else
        result["value"] = _dof.local_dof(name()).standard_values().transpose();
      return result;

    }

    //************************************************************
    std::string Traits::clexulator_point_prepare_string(Structure const &_prim,
                                                        std::map<UnitCellCoord, std::set<UnitCellCoord> > const &_nhood,
                                                        PrimNeighborList &_nlist,
                                                        std::vector<BasisSet> const &site_bases,
                                                        std::string const &indent) const {

      std::stringstream ss;
      if(val_traits().global()) {
        ss <<
           indent << "  if(m_params.eval_mode(m_" << name() << "_var_param_key) != ParamPack::READ) {\n";
        for(Index a = 0; a < _prim.global_dof(name()).size(); ++a) {
          ss << indent << "    ParamPack::Val<Scalar>::set(m_params, m_" << name() << "_var_param_key, " << a
             << ", eval_" << name() << "_var(" << a << "));\n";
        }
        ss << indent << "  }\n";

        if(requires_site_basis()) {
          ss <<
             indent << "  if(m_params.eval_mode(m_" << site_basis_name() << "_param_key) != ParamPack::READ) {\n";
          for(Index f = 0; f < site_bases[0].size(); f++) {
            ss <<
               indent << "    ParamPack::Val<Scalar>::set(m_params, m_" << site_basis_name() << "_param_key, " << f
               << ", eval_" << site_basis_name() << "_" << f << "<Scalar>());\n";
          }
          ss << indent << "  }\n";
        }
      }
      else {
        ss << indent << "switch(nlist_ind) {\n";
        for(auto const &nbor : _nhood) {
          std::stringstream ssvar, ssfunc;
          ss << indent << "case " << _nlist.neighbor_index(nbor.first) << ":\n";
          //Index n = nbor.first;
          //std::cout << "neighborhood of nbor.first: " << nbor.first << ": \n";

          //Put neighborhood in a sensible order:
          std::map<Index, std::set<Index> > sublat_nhood;
          for(auto const &ucc : nbor.second) {
            sublat_nhood[ucc.sublattice()].insert(_nlist.neighbor_index(ucc));
            //std::cout << "ucc : " << ucc << "; n: " << _nlist.neighbor_index(ucc)  << "\n";
          }

          for(auto const &sublat : sublat_nhood) {
            Index b = sublat.first;
            for(Index n : sublat.second) {

              if(!_prim.basis()[b].has_dof(name()))
                continue;

              for(Index a = 0; a < _prim.basis()[b].dof(name()).size(); ++a) {
                ssvar << indent << "    ParamPack::Val<Scalar>::set(m_params, m_" << name() << "_var_param_key, " << a << ", " << n
                      << ", eval_" << name() << "_var_" << b << "_" << a << "(" << n << "));\n";
              }


              if(requires_site_basis()) {
                for(Index f = 0; f < site_bases[b].size(); f++) {
                  ssfunc << indent << "    ParamPack::Val<Scalar>::set(m_params, m_" << site_basis_name() << "_param_key, " << f << ", " << n
                         << ", eval_" << site_basis_name() << '_' << b << '_' << f << "<Scalar>(" << n << "));\n";
                }
              }
            }
          }

          ss <<
             indent << "  if(m_params.eval_mode(m_" << name() << "_var_param_key) != ParamPack::READ) {\n" <<
             ssvar.str() <<
             indent << "  }\n";

          if(requires_site_basis()) {
            ss <<
               indent << "  if(m_params.eval_mode(m_" << site_basis_name() << "_param_key) != ParamPack::READ) {\n" <<
               ssfunc.str() <<
               indent << "  }\n";
          }
          ss << indent << "  break;\n";
        }
        ss << indent << "}\n";
      }
      return ss.str();
    }

    //************************************************************

    std::string Traits::clexulator_global_prepare_string(Structure const &_prim,
                                                         std::map<UnitCellCoord, std::set<UnitCellCoord> > const &_nhood,
                                                         PrimNeighborList &_nlist,
                                                         std::vector<BasisSet> const &site_bases,
                                                         std::string const &indent) const {
      std::stringstream ss;

      if(val_traits().global()) {
        ss <<
           indent << "  if(m_params.eval_mode(m_" << name() << "_var_param_key) != ParamPack::READ) {\n";
        for(Index a = 0; a < _prim.global_dof(name()).size(); ++a) {
          ss << indent << "    ParamPack::Val<Scalar>::set(m_params, m_" << name() << "_var_param_key, " << a
             << ", eval_" << name() << "_var(" << a << "));\n";
        }
        ss << indent << "  }\n";

        if(requires_site_basis()) {
          ss <<
             indent << "  if(m_params.eval_mode(m_" << site_basis_name() << "_param_key) != ParamPack::READ) {\n";
          for(Index f = 0; f < site_bases[0].size(); f++) {
            ss << indent << "    ParamPack::Val<Scalar>::set(m_params, m_" << site_basis_name() << "_param_key, " << f
               << ", eval_" << site_basis_name() << "_" << f << "<Scalar>());\n";
          }
          ss << indent << "  }\n";
        }
      }
      else {

        std::map<Index, std::set<Index> > tot_nhood;
        for(auto const &nbor : _nhood)
          for(auto const &ucc : nbor.second)
            tot_nhood[ucc.sublattice()].insert(_nlist.neighbor_index(ucc));

        std::stringstream ssvar, ssfunc;

        for(auto const &nbor : tot_nhood) {
          Index b = nbor.first;
          for(Index n : nbor.second) {

            if(!_prim.basis()[b].has_dof(name()))
              continue;

            for(Index a = 0; a < _prim.basis()[b].dof(name()).size(); ++a) {
              ssvar << indent << "    ParamPack::Val<Scalar>::set(m_params, m_" << name() << "_var_param_key, " << a << ", " << n
                    << ", eval_" << name() << "_var_" << b << "_" << a << "(" << n << "));\n";
            }

            if(requires_site_basis()) {
              for(Index f = 0; f < site_bases[b].size(); f++) {
                ssfunc << indent << "    ParamPack::Val<Scalar>::set(m_params, m_" << site_basis_name() << "_param_key, " << f << ", " << n
                       << ", eval_" << site_basis_name() << "_" << b << "_" << f << "<Scalar>(" << n << "));\n";
              }
            }
          }
        }
        ss <<
           indent << "  if(m_params.eval_mode(m_" << name() << "_var_param_key) != ParamPack::READ) {\n" <<
           ssvar.str() <<
           indent << "  }\n";
        if(requires_site_basis()) {
          ss <<
             indent << "  if(m_params.eval_mode(m_" << site_basis_name() << "_param_key) != ParamPack::READ) {\n" <<
             ssfunc.str() <<
             indent << "  }\n";
        }
      }
      return ss.str();
    }
    //************************************************************

    std::string Traits::clexulator_member_declarations_string(Structure const &_prim,
                                                              std::vector<BasisSet> const &_site_bases,
                                                              std::string const &indent) const {
      std::stringstream stream;
      std::vector<Orbit<PrimPeriodicSymCompare<IntegralCluster> > > asym_unit;
      std::ostream nullstream(0);

      //TODO: It is not ideal to make a copy of the prim, but everything generated in this function should
      //go out of scope. Solving this will involve rethinking which parts of the prim are needed for the
      //function call, and will affect the implementation of the SymCompare classes
      auto _prim_ptr = std::make_shared<const Structure>(_prim);
      make_prim_periodic_asymmetric_unit(_prim_ptr,
                                         CASM_TMP::ConstantFunctor<bool>(true),
                                         TOL,
                                         std::back_inserter(asym_unit),
                                         nullstream);
      /*
      for(Index no = 0; no < asym_unit.size(); no++) {
        Index nb = asym_unit[no][0][0].sublattice();
        if(_site_bases[nb].size() == 0)
          continue;
        stream <<
               indent << "// Occupation Function tables for basis sites in asymmetric unit " << no << ":\n";
        for(Index ne = 0; ne < asym_unit[no].size(); ne++) {
          nb = asym_unit[no][ne][0].sublattice();
          stream <<
                 indent << "//   - basis site " << nb << ":\n";
          for(Index f = 0; f < _site_bases[nb].size(); f++) {
            stream <<
            indent << "double " << "m_" << site_basis_name() << "_" << nb << '_' << f << '[' << _prim.basis()[nb].occupant_dof().size() << "];\n";
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

    std::string Traits::clexulator_public_method_declarations_string(Structure const &_prim,
                                                                     std::vector<BasisSet> const &site_bases,
                                                                     std::string const &indent) const {
      return "";
    }

    //************************************************************

    std::string Traits::clexulator_public_method_definitions_string(Structure const &_prim,
                                                                    std::vector<BasisSet> const &site_bases,
                                                                    std::string const &indent) const  {
      return "";
    }

    //************************************************************

    std::string Traits::clexulator_private_method_declarations_string(Structure const &_prim,
                                                                      std::vector<BasisSet> const &_site_bases,
                                                                      const std::string &indent) const {
      //std::cout << "PRIVATE METHOD DECLARATIONS FOR DOF " << name() << "\n";
      std::stringstream stream;
      if(val_traits().global()) {
        //std::cout << "**GLOBAL PRIVATE METHOD DECLARATIONS FOR DOF " << name() << "\n";
        stream <<
               indent << "double eval_" << name() << "_var(const int &ind) const {\n" <<
               indent << "  return m_global_dof_ptrs[m_" << name() << "_var_param_key.index()]->values()[ind];\n" <<
               indent << "}\n\n";

        stream <<
               indent << "template<typename Scalar>\n" <<
               indent << "Scalar const &" << name() << "_var(const int &ind) const {\n" <<
               indent << "  return " << "ParamPack::Val<Scalar>::get(m_params, m_" << name() << "_var_param_key, ind);\n" <<
               indent << "}\n";

        if(requires_site_basis()) {
          auto visitors = site_function_visitors("nlist_ind");
          BasisSet site_basis = _site_bases[0];
          for(auto const &vis : visitors)
            site_basis.accept(*vis);

          for(Index f = 0; f < site_basis.size(); f++) {
            stream <<
                   indent << "template<typename Scalar>\n" <<
                   indent << "Scalar eval_" << site_basis_name() << '_' << f << "() const {\n" <<
                   indent << "  return " << site_basis[f]->formula() << ";\n" <<
                   indent << "}\n\n";
          }
        }
        return stream.str();
      }

      std::vector<Orbit<PrimPeriodicSymCompare<IntegralCluster> > > asym_unit;
      std::ostream nullstream(0);

      //TODO: It is not ideal to make a copy of the prim, but everything generated in this function should
      //go out of scope. Solving this will involve rethinking which parts of the prim are needed for the
      //function call, and will affect the implementation of the SymCompare classes
      auto _prim_ptr = std::make_shared<const Structure>(_prim);
      make_prim_periodic_asymmetric_unit(_prim_ptr,
                                         CASM_TMP::ConstantFunctor<bool>(true),
                                         TOL,
                                         std::back_inserter(asym_unit),
                                         nullstream);


      Index max_nf = 0;
      Index max_na = 0;
      for(Index no = 0; no < asym_unit.size(); no++) {
        Index nb = asym_unit[no][0][0].sublattice();
        if(_site_bases[nb].size() == 0)
          continue;


        for(Index ne = 0; ne < asym_unit[no].size(); ne++) {
          nb = asym_unit[no][ne][0].sublattice();
          if(!_prim.basis()[nb].has_dof(name()))
            continue;
          stream <<
                 indent << "// " << name() << " evaluators and accessors for basis site " << nb << ":\n";
          max_na = max(max_na, _prim.basis()[nb].dof(name()).size());
          for(Index a = 0; a < _prim.basis()[nb].dof(name()).size(); ++a) {

            stream <<
                   indent << "double eval_" << name() << "_var_" << nb << '_' << a << "(const int &nlist_ind) const {\n" <<
                   indent << "  return m_local_dof_ptrs[m_" << name() << "_var_param_key.index()]->site_value(_l(nlist_ind))[" << a << "];\n" <<
                   indent << "}\n\n";

          }

          if(requires_site_basis()) {

            max_nf = max(max_nf, _site_bases[nb].size());
            auto visitors = site_function_visitors("nlist_ind");
            BasisSet site_basis = _site_bases[nb];
            for(auto const &vis : visitors)
              site_basis.accept(*vis);

            for(Index f = 0; f < site_basis.size(); f++) {
              stream <<
                     indent << "template<typename Scalar>\n" <<
                     indent << "Scalar eval_" << site_basis_name() << "_" << nb << '_' << f << "(const int &nlist_ind) const {\n" <<
                     indent << "  return " << site_basis[f]->formula() << ";\n" <<
                     indent << "}\n\n";
            }
            stream << '\n';
          }
        }

      }
      for(Index a = 0; a < max_na; ++a) {
        stream <<
               indent << "template<typename Scalar>\n" <<
               indent << "Scalar const &" << name() << "_var_" << a << "(const int &nlist_ind) const {\n" <<
               indent << "  return " << "ParamPack::Val<Scalar>::get(m_params, m_" << name() << "_var_param_key, " << a << ", nlist_ind);\n" <<
               indent << "}\n";
      }
      for(Index f = 0; f < max_nf; ++f) {
        stream <<
               indent << "template<typename Scalar>\n" <<
               indent << "Scalar const &" << site_basis_name() << "_" << f << "(const int &nlist_ind) const {\n" <<
               indent << "  return " << "ParamPack::Val<Scalar>::get(m_params, m_" << site_basis_name() << "_param_key, " << f << ", nlist_ind);\n" <<
               indent << "}\n";

      }
      return stream.str();
    }

    //************************************************************

    std::vector<ParamAllocation > Traits::param_pack_allocation(Structure const &_prim,
                                                                std::vector<BasisSet> const &_bases) const {

      std::vector<ParamAllocation> result;

      if(val_traits().global() && _bases.size()) {
        result.push_back(ParamAllocation(std::string(name() + "_var"), _bases[0].size(), Index(1), true));
        return result;
      }

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
        NV = max(NV, site.dof(name()).size());

      //for(Index i = 0; i < NB; i++)
      result.push_back(ParamAllocation(std::string(name() + "_var"), Index(NV), Index(-1), true));

      if(basis_allocation)
        result.push_back(ParamAllocation(site_basis_name(), Index(NB), Index(-1), false));

      return result;

    }

    //************************************************************

    std::string Traits::clexulator_constructor_string(Structure const &_prim,
                                                      std::vector<BasisSet> const &_site_bases,
                                                      const std::string &indent) const {
      std::stringstream stream;
      stream.flags(std::ios::showpoint | std::ios::fixed | std::ios::right);
      stream.precision(10);

      std::vector<Orbit<PrimPeriodicSymCompare<IntegralCluster> > > asym_unit;
      std::ostream nullstream(0);

      //TODO: It is not ideal to make a copy of the prim, but everything generated in this function should
      //go out of scope. Solving this will involve rethinking which parts of the prim are needed for the
      //function call, and will affect the implementation of the SymCompare classes
      auto _prim_ptr = std::make_shared<const Structure>(_prim);
      make_prim_periodic_asymmetric_unit(_prim_ptr,
                                         CASM_TMP::ConstantFunctor<bool>(true),
                                         TOL,
                                         std::back_inserter(asym_unit),
                                         nullstream);
      return stream.str();
    }


    std::vector<std::unique_ptr<FunctionVisitor> > Traits::site_function_visitors(std::string const &nlist_specifier) const {
      std::vector<std::unique_ptr<FunctionVisitor> > result;
      result.push_back(std::unique_ptr<FunctionVisitor>(new VariableLabeler(name(), "%p_var_%f<Scalar>(" + nlist_specifier + ")")));
      return result;
    }

    std::vector<std::unique_ptr<FunctionVisitor> > Traits::clust_function_visitors() const {
      std::vector<std::unique_ptr<FunctionVisitor> > result;
      if(val_traits().global()) {
        result.push_back(std::unique_ptr<FunctionVisitor>(new VariableLabeler(name(), "%p_var<Scalar>(%f)")));
      }
      else {
        if(requires_site_basis())
          result.push_back(std::unique_ptr<FunctionVisitor>(new SubExpressionLabeler(site_basis_name(), site_basis_name() + "_%l<Scalar>(%n)")));
        else
          result.push_back(std::unique_ptr<FunctionVisitor>(new VariableLabeler(name(), "%p_var_%f<Scalar>(%n)")));
      }
      return result;
    }


    std::string Traits::site_basis_description(BasisSet site_bset, Site site) const {
      return std::string();
    }

  }

}
