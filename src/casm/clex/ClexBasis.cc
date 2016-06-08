namespace CASM {


  //********************************************************************

  Index print_clust_basis(std::ostream &stream,
                          BasisSet _clust_basis,
                          IntegralCluster const &_prototype,
                          Index func_ind,
                          int space,
                          char delim) const {
    for(Index np = 0; np < _prototype.size(); np++) {

      stream << std::string(space, ' ');

      stream.setf(std::ios::showpoint, std::ios_base::fixed);
      stream.precision(5);
      stream.width(9);
      _prototype.coordinate(np).print(stream);
      stream << "  basis_index: " << _prototype[np].sublat() << "  clust_index: " << np << " ";
      if(delim)
        stream << delim;
    }
    stream << "\n"
           << "            Basis Functions:\n";

    _clust_basis.set_dof_IDs(sequence<Index>(0, clust.size() - 1));
    _clust_basis.accept(OccFuncLabeler("\\phi_%b_%f(s_%n)"));
    Index i;
    for(i = 0; i < _clust_basis.size(); i++) {
      stream << "              \\Phi_" << func_ind + i << " = " << _clust_basis[i]->tex_formula() << std::endl;
    }
    return _clust_basis.size();
  }


  //************************************************************
  /*
  void print_eci_in(ClexBasis const &_clex_basis,
                    std::ostream &out,
                    std::vector<Orbit<IntegralCluster> > const &_tree) const {
    if(_tree.index.size() != _tree.size())
      get_index();
    if(_tree.subcluster.size() != _tree.size())
      _tree.get_hierarchy();

    out << std::left
        << std::setw(8) << "label"
        << std::setw(8) << "weight"
        << std::setw(8) << "mult"
        << std::setw(8) << "size"
        << std::setw(12) << "length"
        << std::setw(8) << "hierarchy" << std::endl;


    int functioncount = 0;
    for(Index i = 0; i < _tree.size(); i++) {
      for(Index j = 0; j < _tree.size(i); j++) {

        for(Index k = 0; k < m_bset_tree[_tree.orbit(i, j).index()].size(); k++, functioncount++) {

          out << std::left
              << std::setw(8) << functioncount
              << std::setw(8) << 0
              << std::setw(8) << _tree.orbit(i, j).size()
              << std::setw(8) << _tree.prototype(i, j).size()
              << std::setw(12) << _tree.prototype(i, j).max_length();

          // print hierarchy
          out << std::left << std::setw(8) << 0;
          for(Index l = 0; l < _tree.subcluster[ _tree.index[i][j]].size(); l++) {
            out << std::left
                << std::setw(8) << _tree.subcluster[ _tree.index[i][j] ][l];
          }
          out << '\n' << std::flush;
        }

      }
    }

    //std::cout << "finish print_eci_in" << std::endl;

  }
  */
  //********************************************************************

  void ClexBasis::_populate_site_bases(Structure const &_prim) {

    if(bspecs()["basis_functions"]["site_basis_functions"].is_string()) {
      std::string func_type = bspecs()["basis_functions"]["site_basis_functions"].template get<std::string>();

      //std::cout << "Using " << func_type << " site basis functions." << std::endl << std::endl;
      switch(std::tolower(func_type[0])) {
      case 'c': { //chebychev
        for(Index i = 0; i < _asym_unit().size(); i++) {
          std::vector<double> tprob(m_asym_unit.prototype(i)[0].site_occupant().size(),
                                    1.0 / double(_asym_unit().prototype(i)[0].site_occupant().size()));
          m_asym_unit.prototype(i).clust_basis.construct_orthonormal_discrete_functions(_asym_unit().prototype(i)[0].site_occupant(),
                                                                                        tprob,
                                                                                        _asym_unit().prototype(i)[0].sublat(),
                                                                                        asym_unit().prototype(i).clust_group());
          for(Index ne = 0; ne < _asym_unit()[i].size(); ne++)
            m_asym_unit[i][ne].clust_basis.construct_orthonormal_discrete_functions(_asym_unit()[i][ne][0].site_occupant(),
                                                                                    tprob,
                                                                                    _asym_unit()[i][ne][0].sublat(),
                                                                                    asym_unit().prototype(i).clust_group());
        }
        break;
      }
      case 'o': { //occupation
        for(Index i = 0; i < _asym_unit().size(); i++) {
          std::vector<double> tprob(_asym_unit().prototype(i)[0].site_occupant().size(),
                                    0.0);
          if(tprob.size()) {
            tprob[0] = 1.0;
            m_asym_unit.prototype(i).clust_basis.construct_orthonormal_discrete_functions(_asym_unit().prototype(i)[0].site_occupant(),
                                                                                          tprob,
                                                                                          _asym_unit().prototype(i)[0].sublat(),
                                                                                          asym_unit().prototype(i).clust_group());
            for(Index ne = 0; ne < _asym_unit()[i].size(); ne++)
              m_asym_unit[i][ne].clust_basis.construct_orthonormal_discrete_functions(_asym_unit()[i][ne][0].site_occupant(),
                                                                                      tprob,
                                                                                      _asym_unit()[i][ne][0].sublat(),
                                                                                      asym_unit().prototype(i).clust_group());
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
      typedef std::map<std::string, double> SiteProb;
      std::vector<SiteProb> prob_vec(m_b2asym.size());

      auto it = bspecs()["basis_functions"].find("site_basis_functions");
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

      for(Index i = 0; i < _asym_unit().size(); i++) {
        if(_asym_unit().prototype(i)[0].site_occupant().size() < 2)
          continue;
        std::vector<double> tprob(_asym_unit().prototype(i)[0].site_occupant().size(),
                                  0.0);
        if(tprob.size() == 0)
          continue;
        Index b_ind = _asym_unit().prototype(i)[0].sublat();
        double tsum(0);
        for(Index ns = 0; ns < _asym_unit().prototype(i)[0].site_occupant().size(); ns++) {
          if(prob_vec[b_ind].find(_asym_unit().prototype(i)[0].site_occupant()[ns].name) == prob_vec[b_ind].end())
            throw std::runtime_error("In BSPECS.JSON, basis site " + std::to_string(b_ind) + " must have a composition specified for species " + _asym_unit().prototype(i)[0].site_occupant()[ns].name + "\n");

          tprob[ns] = prob_vec[b_ind][_asym_unit().prototype(i)[0].site_occupant()[ns].name];
          tsum += tprob[ns];
        }
        for(Index j = 0; j < tprob.size(); j++)
          tprob[j] /= tsum;
        m_asym_unit.prototype(i).clust_basis.construct_orthonormal_discrete_functions(_asym_unit().prototype(i)[0].site_occupant(),
                                                                                      tprob,
                                                                                      _asym_unit().prototype(i)[0].sublat(),
                                                                                      asym_unit().prototype(i).clust_group());
        for(Index ne = 0; ne < _asym_unit()[i].size(); ne++)
          m_asym_unit[i][ne].clust_basis.construct_orthonormal_discrete_functions(_asym_unit()[i][ne][0].site_occupant(),
                                                                                  tprob,
                                                                                  _asym_unit()[i][ne][0].sublat(),
                                                                                  asym_unit().prototype(i).clust_group());
      }

      //std::cout << "Using concentration-optimized site basis functions." << std::endl << std::endl;
    }
  }

  namespace ClexBasis_impl {

    BasisSet construct_clust_dof_basis(IntegralCluster const &_clust, std::vector<BasisSet const *> const &site_dof_sets) {
      BasisSet result;
      result.set_dof_IDs(_clust.nlist_inds());
      std::vector<SymGroupRep const *> subspace_reps;
      for(BasisSet const *site_bset_ptr : site_dof_sets) {
        if(site_bset_ptr) {
          result.append(*site_bset_ptr);
          subspace_reps.push_back(SymGroupRep::RemoteHandle(_clust.clust_group(),
                                                            site_bset_ptr->basis_symrep_ID()).rep_ptr());
        }
        else {
          subspace_reps.push_back(SymGroupRep::RemoteHandle(_clust.clust_group(),
                                                            SymGroupRepID::identity(0)).rep_ptr());
        }
      }
      result.set_basis_symrep_ID(permuted_direct_sum_rep(*(_clust.permute_rep().rep_ptr()),
                                                         subspace_reps).add_copy_to_master());
      return result;
    }

  }

}
