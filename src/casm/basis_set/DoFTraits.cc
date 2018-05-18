#include "casm/basis_set/DoFTraits.hh"
#include "casm/basis_set/FunctionVisitor.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/symmetry/Orbit_impl.hh"
#include "casm/clusterography/IntegralCluster.hh"
#include "casm/clusterography/ClusterSymCompare_impl.hh"
#include "casm/clusterography/ClusterOrbits_impl.hh"
#include "casm/clex/ClexBasis.hh"

namespace CASM {
  namespace DoF_impl {
    namespace OccupationDoFTraits_impl {
      typedef std::map<std::string, double> SiteProb;
      void prob_vec_from_json(std::vector<SiteProb> &prob_vec, jsonParser const &_bspecs) {
        auto it = _bspecs["basis_functions"].find("site_basis_functions");
        auto end_it = it;
        ++end_it;

        if(it->is_array()) {
          end_it = it->cend();
          it = it->cbegin();
        }

        bool sublat_spec = true;
        Index num_spec = 0;
        for(; it != end_it; ++it, num_spec++) {
          SiteProb tprob;

          auto it2 = (*it)["composition"].cbegin(), end_it2 = (*it)["composition"].cend();
          for(; it2 != end_it2; ++it2) {
            tprob[it2.name()] = it2->template get<double>();
          }

          if(!(it->contains("sublat_indices")) || !sublat_spec) {
            //we're using this block to check for errors *and* set 'sublat_spec'
            if(num_spec > 0) {
              throw std::runtime_error(std::string("Parse error: If multiple 'site_basis_functions' specifications are provided, 'sublat_indices' must be specified for each.\n")
                                       + "   Example: \"site_basis_functions\" : [\n"
                                       + "                {\n"
                                       + "                    \"sublat_indices\" : [0],\n"
                                       + "                    \"composition\" : [ \"SpeciesA\" : 0.2, \"SpeciesB\" : 0.8]\n"
                                       + "                },\n"
                                       + "                {\n"
                                       + "                    \"sublat_indices\" : [1,2],\n"
                                       + "                    \"composition\" : [ \"SpeciesA\" : 0.7, \"SpeciesB\" : 0.3]\n"
                                       + "                }\n"
                                       + "              ]\n");
            }
            else if(num_spec == 0)
              sublat_spec = false;
          }

          if(!sublat_spec) {
            for(auto &_vec : prob_vec)
              _vec = tprob;
          }
          else {
            it2 = (*it)["sublat_indices"].cbegin();
            end_it2 = (*it)["sublat_indices"].cend();
            for(; it2 != end_it2; ++it2) {
              Index b_ind = it2->template get<long>();
              if(!prob_vec[b_ind].empty())
                throw std::runtime_error("Duplicate sublat_indices specified in BSPECS.JSON\n");

              prob_vec[b_ind] = tprob;
            }
          }
        }
      }
    }


    /// \brief Construct the site basis (if DOF_MODE is LOCAL) for a DoF, given its site
    std::vector<BasisSet> OccupationDoFTraits::construct_site_bases(Structure const &_prim,
                                                                    std::vector<Orbit<IntegralCluster, PrimPeriodicSymCompare<IntegralCluster> > > &_asym_unit,
                                                                    jsonParser const &_bspecs) const {
      std::vector<BasisSet> result(_prim.basis().size());
      std::cout << "Inside construct_site_bases!!!\n";


      //std::cout << "Using " << func_type << " site basis functions." << std::endl << std::endl;

      for(Index i = 0; i < _asym_unit.size(); i++) {
        Site const &_site = _asym_unit[i].prototype()[0].sublat_site();

        std::cout << "SITE OCCUPANT ID FOR SITE " << i << " aka " << _asym_unit[i].prototype()[0].sublat() << ": " << _site.site_occupant().ID() << std::endl;

        if(_site.site_occupant().size() < 2)
          continue;

        Index b_ind = _asym_unit[i].prototype()[0].sublat();


        std::vector<double> tprob;
        if(_bspecs["basis_functions"]["site_basis_functions"].is_string()) {
          std::string func_type = _bspecs["basis_functions"]["site_basis_functions"].template get<std::string>();

          switch(std::tolower(func_type[0])) {
          case 'c': { //chebychev
            tprob = std::vector<double>(_site.site_occupant().size(),
                                        1. / double(_site.site_occupant().size()));
            break;
          }
          case 'o': { //occupation
            tprob = std::vector<double>(_site.site_occupant().size(),
                                        0.);
            tprob[0] = 1.;
            break;
          }
          default: {
            throw std::runtime_error(std::string("Parsing BSPECS.json, the specified 'site_basis_function' option -- \"") + func_type + "\" -- does not exist.\n"
                                     + "valid options are 'chebychev' or 'occupation'.\n");
            break;
          }
          }
        }
        else {
          std::vector<OccupationDoFTraits_impl::SiteProb> prob_vec;
          OccupationDoFTraits_impl::prob_vec_from_json(prob_vec, _bspecs);

          double tsum(0);
          for(Index ns = 0; ns < _site.site_occupant().size(); ns++) {
            if(prob_vec[b_ind].find(_site.site_occupant()[ns].name()) == prob_vec[b_ind].end())
              throw std::runtime_error("In BSPECS.JSON, basis site " + std::to_string(b_ind) + " must have a composition specified for species " + _site.site_occupant()[ns].name() + "\n");

            tprob.push_back(prob_vec[b_ind][_site.site_occupant()[ns].name()]);
            tsum += tprob[ns];
          }
          for(Index j = 0; j < tprob.size(); j++)
            tprob[j] /= tsum;
        }


        result[b_ind].construct_orthonormal_discrete_functions(_site.site_occupant(),
                                                               tprob,
                                                               _asym_unit[i].prototype()[0].sublat(),
                                                               SymGroup(_asym_unit[i].equivalence_map(0).first,
                                                                        _asym_unit[i].equivalence_map(0).second));
        for(Index ne = 1; ne < _asym_unit[i].size(); ne++) {
          result[_asym_unit[i][ne][0].sublat()].apply_sym(_asym_unit[i].equivalence_map()[ne][0]);
          result[_asym_unit[i][ne][0].sublat()].accept(OccFuncBasisIndexer(_asym_unit[i][ne][0].sublat()));
        }
      }



      return result;
    }

    //************************************************************
    std::string OccupationDoFTraits::clexulator_point_prepare_string(Structure const &_prim,
                                                                     std::vector<std::pair<Index, UnitCellCoord> > &_nhood,
                                                                     std::vector<BasisSet> const &site_bases,
                                                                     std::string const &indent) const {
      std::stringstream ss;
      ss
          << indent << "if(eval_mode(occ_func_paramkey)==DEFAULT){\n";
      for(auto const &nbor : _nhood) {
        Index n = nbor.first;
        Index b = nbor.second.sublat();
        for(Index f = 0; f < site_bases[b].size(); f++) {
          ss << indent << "  m_occ_func_vals_f" << f << "[" << n  << "] = m_occ_func_" << b << "_" << f << "[m_config_ptr->occ(" << n << ")];\n";
        }
      }
      ss << "}\n";
      return ss.str();
    }

    //************************************************************

    std::string OccupationDoFTraits::clexulator_global_prepare_string(Structure const &_prim,
                                                                      std::vector<std::pair<Index, UnitCellCoord> > &_nhood,
                                                                      std::vector<BasisSet> const &site_bases,
                                                                      std::string const &indent) const {
      std::stringstream ss;
      ss
          << indent << "if(eval_mode(occ_func_paramkey)==DEFAULT){\n";
      for(auto const &nbor : _nhood) {
        Index n = nbor.first;
        Index b = nbor.second.sublat();
        for(Index f = 0; f < site_bases[n].size(); f++) {
          ss << indent << "  m_occ_func_vals_f" << f << "[" << n  << "] = m_occ_func_" << b << "_" << f << "[m_config_ptr->occ(" << n << ")];\n";
        }
      }
      ss << "}\n";
      return ss.str();
    }
    //************************************************************

    std::string OccupationDoFTraits::clexulator_member_declarations_string(Structure const &_prim,
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
                   indent << "double " << "m_occ_func_" << nb << '_' << f << '[' << _prim.basis()[nb].site_occupant().size() << "];\n";
          }
          stream << '\n';
        }

      }
      return stream.str();
    }

    //************************************************************

    std::string OccupationDoFTraits::clexulator_private_method_declarations_string(Structure const &_prim,
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

          stream <<
                 indent << "// Occupation Function accessors for basis site " << nb << ":\n";
          for(Index f = 0; f < _site_bases[nb].size(); f++) {
            stream <<
                   indent << "const double &occ_func_" << nb << '_' << f
                   << "(const int &nlist_ind)const{return " << "m_occ_func_" << nb << '_' << f << "[*(m_occ_ptr+*(m_nlist_ptr+nlist_ind))];}\n";
          }
          stream << '\n';
        }

      }
      return stream.str();
    }

    //************************************************************

    std::string OccupationDoFTraits::clexulator_constructor_string(Structure const &_prim,
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

      for(Index no = 0; no < asym_unit.size(); no++) {
        for(Index ne = 0; ne < asym_unit[no].size(); ne++) {
          Index nb = asym_unit[no][ne][0].sublat();
          for(Index f = 0; f < _site_bases[nb].size(); f++) {
            for(Index s = 0; s < _prim.basis()[nb].site_occupant().size(); s++) {
              OccFuncEvaluator t_eval(s);
              _site_bases[nb][f]->accept(t_eval);
              if(s == 0)
                stream << indent;
              stream << "m_occ_func_" << nb << '_' << f << '[' << s << "] = "
                     << t_eval.value();
              if(s + 1 == _prim.basis()[nb].site_occupant().size())
                stream << ";\n\n";
              else
                stream << ", ";
            }
          }
        }
      }
      return stream.str();
    }

    std::string OccupationDoFTraits::site_basis_description(BasisSet site_bset, Site site) const {
      std::stringstream ss;
      if(site_bset.size() == 0)
        ss << "        [No site basis functions]\n\n";
      std::vector<DoF::RemoteHandle> handles(1, site.site_occupant().handle());
      int s;
      handles[0] = s;
      site_bset.register_remotes(handles);
      for(Index f = 0; f < site_bset.size(); f++) {
        for(s = 0; s < site.site_occupant().size(); s++) {
          if(s == 0)
            ss << "    ";
          ss << "    \\phi_" << site.basis_ind() << '_' << f << '[' << site.site_occupant()[s].name() << "] = "
             << site_bset[f]->remote_eval();
          if(s + 1 == site.site_occupant().size())
            ss << "\n";
          else
            ss << ",   ";
        }
      }
      return ss.str();
    }

  }
}
