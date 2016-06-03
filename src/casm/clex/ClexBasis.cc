namespace CASM {

  void ClexBasis::generate(std::vector<Orbit<IntegralCluster> > const &_orbits,
                           std::vector<DoFType> const &_dof_keys,
                           Index _max_poly_order /*= -1*/) {
    std::vector<DoFType> global_keys;
    std::vector<DoFType> local_keys;
    //separate local_args from global_args
    for(DoFType const &key : dof_keys) {
      if(m_global_bases.find(key) != m_global_bases.end()) {
        global_args.push_back(key);
      }
      else if(m_local_bases.find(key) != m_local_bases.end()) {
        local_keys.push_back
      }
      else {
        throw std::runtime_error(std::string("Attempting to build Clex basis set, but missing degree of freedom \"") + key + "\n");
      }
    }
    m_bset_tree.resize(_orbits.size());

    for(Index i = 0; i < _orbits.size(); i++) {
      m_bset_tree[i].reserve(_orbits[i].size());
      m_bset_tree[i].push_back(generate_clust_basis(_orbits[i].prototype(), local_keys, global_keys));
      for(Index j = 1; j < _orbits[i].size(); j++) {
        m_bset_tree[i].push_back((*(_orbits[i].equivalence_map(j).first))*m_bset_tree[i][0]);
      }
    }
  }

  //*************************************************
  // @param local_args[i][j] is BasisSet for i'th DoFspace at j'th site of cluster
  //                                site  0            site 1            site 2
  // DoFType 0: displacement       {x0,y0,z0}       {x1, y1, z1}       {x2,y2,z2}
  // DoFType 1: configuration       {pA0,pB0}            {}            {pA2, pB2}
  //
  // Step 1:  Get the kroenecker product of cluster permutation with DoF symrep
  //
  //  permutation  |    kronecker prod  |    DoF Symrep (e.g., x--y displacement)
  //   [ 0  1 ]              v                     [cos -sin]
  //   [ 1  0 ]             XkX                    [sin  cos]
  //
  //       [x0]      [  0    0   cos -sin ]   [x0]
  //       [y0]      [  0    0   sin  cos ]   [y0]
  //   S * [x1]  =   [ cos -sin   0    0  ]   [x1]
  //       [x2]      [ sin  cos   0    0  ]   [x2]
  //
  // ----------------------------------------------------------------------
  //
  // Step 2: mix-in @param global_args to get all_argsets
  //
  // GLOBAL ARGS
  //
  // strain             { e1, e2, e3, e4, e5, e6}
  // composition        { comp_a, comp_b }
  //
  // arg_subsets =   [ {x0,y0,z0,x1,y1,z1,x2,y2,z2},
  //                   {pA0,pB0,pA2,pB2},
  //                   {e1,e2,e3,e4,e5,e6},
  //                   {comp_a,comp_b}]
  //
  BasisSet ClexBasis::_construct_clust_basis(IntegralCluster const &prototype,
                                             std::vector<DoFType> const &local_keys,
                                             std::vector<DoFType> const &global_keys,
                                             Index max_poly_order) const {
    //std::cout<<"In IntegralCluster::generate_clust_basis, the size of this cluster is:"<<size()<<std::endl;
    //std::cout<<"valid_index evaluates to:"<<valid_index(max_poly_order)<<std::endl;

    if(!valid_index(max_poly_order))
      max_poly_order = prototype.size();
    //std::cout<<"Max_poly_order "<<max_poly_order<<std::endl;

    std::vector<BasisSet const *> arg_subsets;
    for(DoFType const &key : global_keys) {
      arg_subsets.push_back(&(m_global_bases[key]));
    }

    std::vector<BasisSet> all_local;
    all_local.reserve(local_keys.size());

    //Loop over dof's
    for(DoFType const &key : local_keys) {
      // Make copies of local arguments to ensure that they are distinguishable by their DoF_IDs
      // i.e., make copies in 'tlocal' and reset the DoF_IDs to {0,1,2,etc...}
      std::vector<BasisSet> const &arg_vec(m_local_bases[key]);
      std::vector<BasisSet> tlocal;
      tlocal.reserve(prototype.size());
      std::vector<BasisSet const *> site_args(size(), nullptr);
      //Loop over sites
      for(Index i = 0; i < prototype.size(); i++) {
        if(arg_vec[prototype[i].sublat()].size()) {
          tlocal.push_back(arg_vec[prototype[i].sublat()]);
          tlocal.back().set_dof_IDs(std::vector<Index>(1, i));
          site_args[i] = &tlocal.back();
        }
      }
      all_local.push_back(ClexBasis_impl::construct_clust_dof_basis(prototype, site_args));
      if(all_local.back().size())
        arg_subsets.push_back(&(all_local.back()));
    }

    return m_generating_function(prototype, arg_subsets, max_poly_order, 1);
  }

  //********************************************************************

  Index print_clust_basis(std::ostream &stream,
                          ClexBasis const &_basis_set,
                          IntegralCluster const &_prototype,
                          Index orbit_ind,
                          Index func_ind,
                          int space,
                          char delim,
                          COORD_TYPE mode) const {
    if(mode == COORD_DEFAULT)
      mode = COORD_MODE::CHECK();
    COORD_MODE C(mode);
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
    BasisSet tbasis(_basis_set.clust_basis(orbit_ind, 0));
    tbasis.set_dof_IDs(sequence<Index>(0, clust.size() - 1));
    tbasis.accept(OccFuncLabeler("\\phi_%b_%f(s_%n)"));
    Index i;
    for(i = 0; i < tbasis.size(); i++) {
      stream << "              \\Phi_" << func_ind + i << " = " << tbasis[i]->tex_formula() << std::endl;
    }
    return tbasis.size();
  }

  //********************************************************************

  // Divide by multiplicity. Same result as evaluating correlations via orbitree.
  std::vector<std::string> orbit_function_cpp_strings(ClexBasis::BSetOrbit _bset_orbit, // used as temporary
                                                      Orbit<IntegralCluster> const &_clust_orbit,
                                                      PrimNeighborList const &_nlist,
                                                      std::vector<FunctionVisitor *> const &labelers) {
    std::string prefix, suffix;
    std::vector<std::string> formulae(_bset_orbit[0].size(), std::string());
    if(_clust_orbit.size() > 1) {
      prefix = "(";
      suffix = ")/" + std::to_string(_clust_orbit.size()) + ".";
    }

    for(Index ne = 0; ne < _bset_orbit.size(); ne++) {
      _bset_orbit[ne].set_dof_IDs(_nlist.indices(_clust_orbit[ne].elements()));
      for(Index nl = 0; nl < labelers.size(); nl++)
        _bset_orbit[ne].accept(*labelers[nl]);
    }

    for(Index nf = 0; nf < _bset_orbit[0].size(); nf++) {
      for(Index ne = 0; ne < _bset_orbit.size(); ne++) {
        if(!_bset_orbit[ne][nf] || (_bset_orbit[ne][nf]->is_zero()))
          continue;

        if(formulae[nf].empty())
          formulae[nf] += prefix;
        else if((_bset_orbit[ne][nf]->formula())[0] != '-' && (_bset_orbit[ne][nf]->formula())[0] != '+')
          formulae[nf] += " + ";

        formulae[nf] += _bset_orbit[ne][nf]->formula();
      }

      if(!formulae[nf].empty())
        formulae[nf] += suffix;
    }
    return formulae;
  }

  //********************************************************************

  /// nlist_index is the index of the basis site in the neighbor list
  std::vector<std::string>  flower_function_cpp_strings(ClexBasis::BSetOrbit _bset_orbit, // used as temporary
                                                        Orbit<IntegralCluster> const &_clust_orbit,
                                                        PrimNeighborList const &_nlist,
                                                        std::vector<FunctionVisitor *> const &labelers,
                                                        Index sublat_index) {


    std::vector<std::string> formulae(_bset_orbit[0].size(), std::string());
    std::string prefix, suffix;

    //normalize by multiplicity (by convention)
    if(_clust_orbit.size() > 1) {
      prefix = "(";
      suffix = ")/" + std::to_string(_clust_orbit.size()) + ".";
    }

    // loop over equivalent clusters
    for(Index ne = 0; ne < _clust_orbit.size(); ne++) {

      // for each site in the cluster that belongs to the target sublattice,
      // translate the cluster so that that site is within the (0,0,0) cell
      // Use std::set to ensure that each translation is only attempted once
      std::set<UnitCell> attempted_trans;

      // loop over cluster sites
      auto it(_clust_orbit[ne].cbegin()), end_it(_clust_orbit[ne].cbegin());
      for(; it != end_it; ++it) {

        // Continue if the cluster site doesn't belong to the target sublattice, or if we have allready attempted the translation
        if(sublat_index != it -> sublat() || !(attempted_trans.insert(it->unitcell()).second))
          continue;

        _bset_orbit[ne].set_dof_IDs(_nlist.indices((_clust_orbit[ne] - it->unitcell()).elements()));

        for(Index nl = 0; nl < labelers.size(); nl++)
          _bset_orbit[ne].accept(*labelers[nl]);

        for(Index nf = 0; nf < _bset_orbit[ne].size(); nf++) {
          if(!_bset_orbit[ne][nf] || (_bset_orbit[ne][nf]->is_zero()))
            continue;

          if(formulae[nf].empty())
            formulae[nf] += prefix;
          else if((_bset_orbit[ne][nf]->formula())[0] != '-' && (_bset_orbit[ne][nf]->formula())[0] != '+')
            formulae[nf] += " + ";

          formulae[nf] += _bset_orbit[ne][nf]->formula();

        }
      }
    }

    // append suffix to all formulae
    for(Index nf = 0; nf < formulae.size(); nf++) {
      if(!formulae[nf].empty())
        formulae[nf] += suffix;
    }

    return formulae;
  }

  //********************************************************************
  //TODO:
  /// b_index is the basis site index, f_index is the index of the configurational site basis function in Site::occupant_basis
  /// nlist_index is the index of the basis site in the neighbor list
  std::vector<std::string>  delta_occfunc_flower_function_cpp_strings(ClexBasis const &_clex_basis,
                                                                      BasisSet site_basis, // passed by value because we use it as a temporary
                                                                      const std::vector<FunctionVisitor *> &labelers,
                                                                      Index nlist_index,
                                                                      Index b_index,
                                                                      Index f_index) {

    std::vector<BasisSet> const &bset_orbit(m_bset_tree[_orbit.index()]);
    std::vector<std::string> formulae(bset_orbit[0].size(), std::string());
    std::string suffix;
    Index ib;
    //normalize by multiplicity (by convention)
    if(size()*_orbit.prototype.size() > 1) {
      formulae.resize(bset_orbit[0].size(), std::string("("));
      suffix = ")/" + std::to_string(size()) + ".";
    }

    for(Index ne = 0; ne < _orbit.size(); ne++) {
      //std::cout << " **** for ne = " << ne << ":\n";
      for(Index nt = 0; nt < _orbit[ne].trans_nlists().size(); nt++) {
        ib = _orbit[ne].trans_nlist(nt).find(nlist_index);
        if(ib == _orbit[ne].size())
          continue;
        //std::cout << " **** for translist: " << _orbit[ne].trans_nlist(nt) << ":\n";
        _orbit[ne].set_nlist_inds(_orbit[ne].trans_nlist(nt));

        site_basis.set_dof_IDs(std::vector<Index>(1, _orbit[ne][ib].nlist_ind()));

        BasisSet quotient_basis = bset_orbit[ne].poly_quotient_set(site_basis[f_index]);
        for(Index nl = 0; nl < labelers.size(); nl++)
          quotient_basis.accept(*labelers[nl]);

        for(Index nf = 0; nf < quotient_basis.size(); nf++) {

          if((quotient_basis[nf]->formula()) == "0")
            continue;

          if(formulae[nf].size() > 1)
            formulae[nf] += " + ";
          formulae[nf] += "(" + quotient_basis[nf]->formula() + ")";
        }
      }
    }

    // Make sure that formulae that evaluate to zero have an empty sting.
    for(Index nf = 0; nf < bset_orbit[0].size(); nf++) {
      if(formulae[nf].size() <= 1)
        formulae[nf].clear();
      else
        formulae[nf] += suffix;
    }

    return formulae;
  }

  //***********************************************

  void print_proto_clust_funcs(ClexBasis const &_clex_basis,
                               std::ostream &out,
                               std::vector<Orbit<IntegralCluster> > const &_tree) const {
    //Prints out all prototype clusters (CLUST file)

    if(_tree.index.size() != _tree.size())
      _tree.get_index();

    out << "COORD_MODE = " << COORD_MODE::NAME() << std::endl << std::endl;

    out.flags(std::ios::showpoint | std::ios::fixed | std::ios::left);
    out.precision(5);
    for(Index no = 0; no < _tree.asym_unit().size(); no++) {
      out << "Asymmetric unit " << no + 1 << ":\n";
      for(Index ne = 0; ne < _tree.asym_unit()[no].size(); ne++) {
        Index b = _tree.asym_unit()[no][ne][0].sublat();
        out << "  Basis site " << b << ":\n"
            << "  ";
        _tree.asym_unit()[no][ne][0].print(out);
        out << "\n";
        if(_tree.asym_unit()[no][ne].clust_basis.size() == 0)
          out << "        [No site basis functions]\n\n";
        for(Index f = 0; f < _tree.asym_unit()[no][ne].clust_basis.size(); f++) {
          for(Index s = 0; s < _tree.asym_unit()[no][ne][0].site_occupant().size(); s++) {
            if(s == 0)
              out << "    ";
            out << "    \\phi_" << b << '_' << f << '[' << asym_unit()[no][ne][0].site_occupant()[s].name << "] = "
                << asym_unit()[no][ne].clust_basis[f]->eval(std::vector<Index>(1, asym_unit()[no][ne][0].site_occupant().ID()), std::vector<Index>(1, s));
            if(s + 1 == asym_unit()[no][ne][0].site_occupant().size())
              out << "\n";
            else
              out << ",   ";
          }
        }
      }
    }
    out << "\n\n";
    Index nf = 0;
    for(Index i = 0; i < size(); i++) {
      if(size(i) != 0) out << "** Branch " << i << " ** \n" << std::flush;
      for(Index j = 0; j < size(i); j++) { //Loops over all i sized Orbits
        out << "      ** " << index[i][j] << " of " << Norbits << " Orbits **"
            << "  Orbit: " << i << " " << j
            << "  Points: " << orbit(i, j).prototype.size()
            << "  Mult: " << orbit(i, j).size()
            << "  MinLength: " << orbit(i, j).prototype.min_length()
            << "  MaxLength: " << orbit(i, j).prototype.max_length()
            << '\n' << std::flush;

        out << "            " << "Prototype" << " of " << orbit(i, j).size() << " Equivalent Clusters in Orbit " << index[i][j] << '\n' << std::flush;
        prototype(i, j).print_clust_basis(out, nf, 8, '\n');
        nf += prototype(i, j).clust_basis.size();
        out << "\n\n" << std::flush;
      }
      if(size(i) != 0) out << '\n' << std::flush;
    }
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
