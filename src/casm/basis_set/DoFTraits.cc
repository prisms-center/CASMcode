#include "casm/basis_set/DoFTraits.hh"
#include "casm/basis_set/FunctionVisitor.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/clex/ClexBasis.hh"
#include "casm/clusterography/IntegralCluster.hh"
#include "casm/clusterography/ClusterOrbits.hh"

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
      std::vector<BasisSet> result(_prim.basis.size());

      if(_bspecs["basis_functions"]["site_basis_functions"].is_string()) {
        std::string func_type = _bspecs["basis_functions"]["site_basis_functions"].template get<std::string>();

        //std::cout << "Using " << func_type << " site basis functions." << std::endl << std::endl;
        switch(std::tolower(func_type[0])) {
        case 'c': { //chebychev
          for(Index i = 0; i < _asym_unit.size(); i++) {
            Site const &_site = _prim.basis[_asym_unit[i].prototype()[0].sublat()];
            std::vector<double> tprob(_site.site_occupant().size(),
                                      1.0 / double(_site.site_occupant().size()));
            if(tprob.empty())
              continue;
            BasisSet tbasis;
            tbasis.construct_orthonormal_discrete_functions(_site.site_occupant(),
                                                            tprob,
                                                            _asym_unit[i].prototype()[0].sublat(),
                                                            SymGroup(_asym_unit[i].equivalence_map(0).first,
                                                                     _asym_unit[i].equivalence_map(0).second));
            for(Index ne = 0; ne < _asym_unit[i].size(); ne++) {
              result[_asym_unit[i][ne][0].sublat()] = tbasis;
              result[_asym_unit[i][ne][0].sublat()].apply_sym(_asym_unit[i].equivalence_map()[ne][0]);
              result[_asym_unit[i][ne][0].sublat()].accept(OccFuncBasisIndexer(_asym_unit[i][ne][0].sublat()));
            }
          }
          break;
        }
        case 'o': { //occupation
          for(Index i = 0; i < _asym_unit.size(); i++) {
            Site const &_site = _prim.basis[_asym_unit[i].prototype()[0].sublat()];
            std::vector<double> tprob(_site.site_occupant().size(), 0.0);
            if(tprob.empty())
              continue;

            tprob[0] = 1.0;

            BasisSet tbasis;
            tbasis.construct_orthonormal_discrete_functions(_site.site_occupant(),
                                                            tprob,
                                                            _asym_unit[i].prototype()[0].sublat(),
                                                            SymGroup(_asym_unit[i].equivalence_map(0).first,
                                                                     _asym_unit[i].equivalence_map(0).second));
            for(Index ne = 0; ne < _asym_unit[i].size(); ne++) {
              result[_asym_unit[i][ne][0].sublat()] = tbasis;
              result[_asym_unit[i][ne][0].sublat()].apply_sym(_asym_unit[i].equivalence_map()[ne][0]);
              result[_asym_unit[i][ne][0].sublat()].accept(OccFuncBasisIndexer(_asym_unit[i][ne][0].sublat()));
            }
          }
          break;
        }
        default: {
          throw std::runtime_error(std::string("Parsing BSPECS.json, the specified 'site_basis_function' option -- \"") + func_type + "\" -- does not exist.\n"
                                   + "valid options are 'chebychev' or 'occupation'.\n");
          break;
        }
        }
      }
      else { // composition-optimized functions

        std::vector<OccupationDoFTraits_impl::SiteProb> prob_vec;
        OccupationDoFTraits_impl::prob_vec_from_json(prob_vec, _bspecs);

        for(Index i = 0; i < _asym_unit.size(); i++) {
          Site const &_site = _prim.basis[_asym_unit[i].prototype()[0].sublat()];
          if(_site.site_occupant().size() < 2)
            continue;
          std::vector<double> tprob(_site.site_occupant().size(),
                                    0.0);
          if(tprob.size() == 0)
            continue;
          Index b_ind = _asym_unit[i].prototype()[0].sublat();
          double tsum(0);
          for(Index ns = 0; ns < _site.site_occupant().size(); ns++) {
            if(prob_vec[b_ind].find(_site.site_occupant()[ns].name()) == prob_vec[b_ind].end())
              throw std::runtime_error("In BSPECS.JSON, basis site " + std::to_string(b_ind) + " must have a composition specified for species " + _site.site_occupant()[ns].name() + "\n");

            tprob[ns] = prob_vec[b_ind][_site.site_occupant()[ns].name()];
            tsum += tprob[ns];
          }
          for(Index j = 0; j < tprob.size(); j++)
            tprob[j] /= tsum;

          BasisSet tbasis;
          tbasis.construct_orthonormal_discrete_functions(_site.site_occupant(),
                                                          tprob,
                                                          b_ind,
                                                          SymGroup(_asym_unit[i].equivalence_map(0).first,
                                                                   _asym_unit[i].equivalence_map(0).second));

          for(Index ne = 0; ne < _asym_unit[i].size(); ne++) {
            result[_asym_unit[i][ne][0].sublat()] = tbasis;
            result[_asym_unit[i][ne][0].sublat()].apply_sym(_asym_unit[i].equivalence_map()[ne][0]);
            result[_asym_unit[i][ne][0].sublat()].accept(OccFuncBasisIndexer(_asym_unit[i][ne][0].sublat()));
          }
        }
        //std::cout << "Using concentration-optimized site basis functions." << std::endl << std::endl;

      }
      return result;
    }

    //************************************************************

    void print_clexulator_member_definitions(std::ostream &stream, const SiteOrbitree &tree, const std::string &indent)const {
      const SiteOrbitBranch &asym_unit(tree.asym_unit());
      for(Index no = 0; no < asym_unit.size(); no++) {
        if(asym_unit[no].size() == 0 || asym_unit[no][0].clust_basis.size() == 0)
          continue;

        stream <<
               indent << "// Occupation Function tables for basis sites in asymmetric unit " << no << ":\n";
        for(Index ne = 0; ne < asym_unit[no].size(); ne++) {
          Index b = asym_unit[no][ne][0].basis_ind();
          stream << indent << "//   - basis site " << b << ":\n";
          for(Index f = 0; f < asym_unit[no][ne].clust_basis.size(); f++) {
            stream <<
                   indent << "double " << "m_occ_func_" << b << '_' << f << '[' << asym_unit[no][ne][0].site_occupant().size() << "];\n";
          }
          stream << '\n';
        }

      }

    }
    //************************************************************
    void print_clexulator_private_method_definitions(std::ostream &stream, const SiteOrbitree &tree, const std::string &indent) const {
      const SiteOrbitBranch &asym_unit(tree.asym_unit());
      for(Index no = 0; no < asym_unit.size(); no++) {
        if(asym_unit[no].size() == 0 || asym_unit[no][0].clust_basis.size() == 0)
          continue;

        for(Index ne = 0; ne < asym_unit[no].size(); ne++) {
          Index b = asym_unit[no][ne][0].basis_ind();
          stream <<
                 indent << "// Occupation Function accessors for basis site " << b << ":\n";
          for(Index f = 0; f < asym_unit[no][ne].clust_basis.size(); f++) {
            stream <<
                   indent << "const double &occ_func_" << b << '_' << f << "(const int &nlist_ind)const{return " << "m_occ_func_" << b << '_' << f << "[*(m_occ_ptr+*(m_nlist_ptr+nlist_ind))];}\n";
          }
          stream << '\n';
        }

      }
    }

    //************************************************************

    void print_to_clexulator_constructor(std::ostream &stream, const SiteOrbitree &tree, const std::string &indent) const {
      stream.flags(std::ios::showpoint | std::ios::fixed | std::ios::right);
      stream.precision(10);

      const SiteOrbitBranch &asym_unit(tree.asym_unit());
      for(Index no = 0; no < asym_unit.size(); no++) {
        for(Index ne = 0; ne < asym_unit[no].size(); ne++) {
          Index b = asym_unit[no][ne][0].basis_ind();
          for(Index f = 0; f < asym_unit[no][ne].clust_basis.size(); f++) {
            for(Index s = 0; s < asym_unit[no][ne][0].site_occupant().size(); s++) {
              if(s == 0)
                stream << indent;
              stream << "m_occ_func_" << b << '_' << f << '[' << s << "] = "
                     << asym_unit[no][ne].clust_basis[f]->eval(Array<Index>(1, asym_unit[no][ne][0].site_occupant().ID()), Array<Index>(1, s));
              if(s + 1 == asym_unit[no][ne][0].site_occupant().size())
                stream << ";\n\n";
              else
                stream << ", ";
            }
          }
        }
      }
    }



  }
}
