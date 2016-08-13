#include "casm/clusterography/SiteCluster.hh"
#include "casm/basis_set/FunctionVisitor.hh"

namespace CASM {

  SiteCluster::SiteCluster(const Lattice &init_home) : GenericCluster<Site>(init_home) {
  }

  //*************************************************

  void SiteCluster::push_back(const Site &new_site) {
    GenericCluster<Site>::push_back(new_site);
    back().set_nlist_ind(size() - 1);
    return;
  }

  //*************************************************
  ReturnArray<Index> SiteCluster::nlist_inds() const {
    Array<Index> inds(size(), 0);
    for(Index i = 0; i < size(); i++)
      inds[i] = at(i).nlist_ind();

    return inds;
  }
  //*************************************************

  void SiteCluster::set_nlist_inds(const Array<Index> &new_indices) {
    Array<Index> old_indices(size(), 0);

    for(Index i = 0; i < size(); i++) {
      old_indices[i] = at(i).nlist_ind();
      at(i).set_nlist_ind(new_indices[i]);
    }
    //std::cout << "Updating indices " << old_indices << " to " << new_indices << "\n";
    clust_basis.set_dof_IDs(new_indices);
    //ccd_basis.update_dof_IDs(old_indices, new_indices);

    return;
  }

  //*************************************************

  const Array<Array<Index> > &SiteCluster::trans_nlists() const {
    return m_trans_nlist_inds;
  }

  //*************************************************

  const Array<Index> &SiteCluster::trans_nlist(Index i) const {
    return m_trans_nlist_inds[i];
  }

  //*************************************************

  void SiteCluster::add_trans_nlist(const Array<Index> &new_nlist) {
    if(new_nlist.size() != size()) {
      std::cerr << "CRITICAL ERROR: Size mismatch in SiteCluster::add_trans_nlist().\n"
                << "                Exiting...\n";
    }
    m_trans_nlist_inds.push_back(new_nlist);
  }


  //*************************************************
  // @param local_args[i][j] is BasisSet for i'th DoFspace at j'th site of cluster
  //                            site  0            site 1            site 2
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
  void SiteCluster::generate_clust_basis(multivector<BasisSet const *>::X<2> const &local_args,
                                         std::vector<BasisSet const *> const &global_args,
                                         Index max_poly_order) {
    //std::cout<<"In SiteCluster::generate_clust_basis, the size of this cluster is:"<<size()<<std::endl;
    //std::cout<<"valid_index evaluates to:"<<valid_index(max_poly_order)<<std::endl;

    clust_basis.set_dof_IDs(nlist_inds());
    if(!valid_index(max_poly_order))
      max_poly_order = size();
    //std::cout<<"Max_poly_order "<<max_poly_order<<std::endl;

    Array<BasisSet const *> arg_subsets;
    for(BasisSet const *global_ptr : global_args) {
      if(global_ptr)
        arg_subsets.push_back(global_ptr);
    }

    Array<BasisSet> all_local;
    all_local.reserve(local_args.size());

    //Loop over dof's
    for(Index d = 0; d < local_args.size(); d++) {
      assert(local_args[d].size() == size() && "In SiteCluster::generate_clust_basis(), local_args must have same size as cluster.");
      // Make copies of local arguments to ensure that they are distinguishable by their DoF_IDs
      // i.e., make copies in 'tlocal' and reset the DoF_IDs to {0,1,2,etc...}
      std::vector<BasisSet> tlocal;
      tlocal.reserve(local_args[d].size());
      std::vector<BasisSet const *> site_args(size(), nullptr);
      //Loop over sites
      for(Index i = 0; i < local_args[d].size(); i++) {
        if(local_args[d][i]) {
          tlocal.push_back(*local_args[d][i]);
          tlocal.back().set_dof_IDs(Array<Index>(1, at(i).nlist_ind()));
          site_args[i] = &tlocal.back();
        }
      }
      all_local.push_back(SiteCluster_impl::construct_clust_dof_basis(*this, site_args));
      if(all_local.back().size())
        arg_subsets.push_back(&(all_local.back()));
    }

    //std::cerr << "WARNING: THIS VERSION OF CASM CANNOT PRODUCE CLUSTER FUNCTIONS!! YOU WILL HAVE NO CORRELATIONS\n";
    // BasisSet::construct_invariant_cluster_polynomials() does the heavy lifting
    //clust_basis.construct_invariant_cluster_polynomials(site_args, global_args, clust_group, permute_group, max_poly_order);
    clust_basis.construct_invariant_polynomials(arg_subsets, clust_group(), max_poly_order, 1);
  }

  //*************************************************
  /*
   * Using the permuation group of the cluster, we
   * count over the allowed occupants of each site
   * and shuffle them with every permutation. We then
   * push back every unique decoration onto decor_map
   * which we can later use to actually make new clusters
   */
  //*************************************************

  //Take decor_map out of the class definiton and have this function return
  //an equivalent object
  ReturnArray<Array<int> > SiteCluster::get_decor_map() const {
    Array<Array<int> > decor_map;
    //I start counting at 1 because I don't want to decorate with the "background" occupant
    //just everything else
    Array<int> init(this->size(), 1);
    Array<int> endocc;

    for(Index i = 0; i < this->size(); i++) {
      endocc.push_back(at(i).site_occupant().size() - 1);
    }

    Counter<Array<int> > ornaments(init, endocc, Array<int>(this->size(), 1));

    if(decor_map.size() != 0) {
      std::cerr << "WARNING in SiteCluster::get_decor_map" << std::endl;
      std::cerr << "You requested populating your decor_map, but it isn't empty.";
      std::cerr << " your current map is about to be obliterated." << std::endl;
      decor_map.clear();
    }

    do {
      bool new_decoration = true;
      for(Index i = 0; i < permute_rep().size(); i++) {
        if(decor_map.contains((permute_rep()[i]->get_permutation())->permute(ornaments()))) {
          new_decoration = false;
          break;
        }
      }

      //If you store permuations of the unique decorations seperately,
      //then you can just compare your new arrangemets(ornaments) to the
      //permuted group instead of the unique group.
      if(new_decoration) {
        decor_map.push_back(ornaments);
      }
    }
    while(++ornaments);

    return decor_map;
  }

  //*************************************************
  /*
   * Using the permuation group of the cluster, we
   * count over the allowed occupants of each site
   * and shuffle them with every permutation. We then
   * push back every unique decoration (and equivalents) onto decor_map
   * which we can use later
   */
  //*************************************************

  ReturnArray< Array<int> > SiteCluster::get_full_decor_map() const {
    //std::cout << "begin SiteCluster::get_full_decor_map()" << std::endl;

    Array< Array<int> > decor_map; //[uniq][site]

    //I start counting at 0 because I want to include all decorations (including background)

    Array<int> init(this->size(), 0);
    Array<int> endocc;

    for(Index i = 0; i < this->size(); i++) {
      endocc.push_back(at(i).site_occupant().size() - 1);
    }

    Counter<Array<int> > ornaments(init, endocc, Array<int>(this->size(), 1));

    do {
      bool new_decoration = true;
      for(Index i = 0; i < permute_rep().size(); i++) {
        // check if a symmetry equivalent to this arrangement has already been enumerated
        if(decor_map.contains((permute_rep()[i]->get_permutation())->permute(ornaments()))) {
          new_decoration = false;
          break;
        }
      }

      //If you store permuations of the unique decorations seperately,
      //then you can just compare your new arrangemets(ornaments) to the
      //permuted group instead of the unique group.
      if(new_decoration) {
        decor_map.push_back(ornaments);
      }
    }
    while(++ornaments);

    //std::cout << "finish SiteCluster::get_full_decor_map()" << std::endl;
    return decor_map;
  }



  //*************************************************
  /*
   * Given a map of how to decorate a cluster, make
   * copies of said cluster and change occ_ind at each
   * site accordingly.
   * Maybe this should return an OrbitBranch? See
   * Structure::bedazzle
   */
  //*************************************************

  ReturnArray<SiteCluster> SiteCluster::get_decorations(const Array<Array<int> > &dmap) const {
    if(dmap[0].size() != this->size()) {
      std::cerr << "ERROR in SiteCluster::get_decorations" << std::endl;
      std::cerr << "The provided map is for a cluster of the wrong size!" << std::endl;
      exit(70);
    }

    Array<SiteCluster> decorations;
    for(Index i = 0; i < dmap.size(); i++) {
      SiteCluster decorclust(*this);

      decorclust.decorate(dmap[i]);

      decorations.push_back(decorclust);
    }
    return decorations;
  }

  //*************************************************

  SiteCluster &SiteCluster::permute(const Array<Index> &iperm) {
    GenericCluster<Site>::permute(iperm);

    return *this;
  }

  //*************************************************

  SiteCluster &SiteCluster::permute(const Permutation &perm) {
    return permute(perm.perm_array());
  }
  //*************************************************

  SiteCluster &SiteCluster::apply_sym_no_trans(const SymOp &op) {
    GenericCluster<Site>::apply_sym_no_trans(op);
    clust_basis.apply_sym(op);
    return *this;
  }

  //*************************************************
  /*
   * Uses basis_ind() of each site in the cluster to extract
   * values from a given bitstring. The extracted values
   * are then returned in a new array that has the same
   * length as the cluster.
   */
  //*************************************************

  ReturnArray<int> SiteCluster::get_occ_array(const Array<int> &bitstring) const {
    Array<int> occ_array;

    for(Index i = 0; i < size(); i++) {
      occ_array.push_back(bitstring[at(i).basis_ind()]);
    }

    return occ_array;
  }

  //\John G 010413

  //*************************************************

  SiteCluster &SiteCluster::apply_sym(const SymOp &op) {
    GenericCluster<Site>::apply_sym(op);
    clust_basis.apply_sym(op);
    return *this;
  }

  //********************************************
  void SiteCluster::print_clust_basis(std::ostream &stream, Index begin_ind, int space, char delim, COORD_TYPE mode) const {
    if(mode == COORD_DEFAULT)
      mode = COORD_MODE::CHECK();
    COORD_MODE C(mode);
    for(Index np = 0; np < size(); np++) {

      stream << std::string(space, ' ');

      stream.setf(std::ios::showpoint, std::ios_base::fixed);
      stream.precision(5);
      stream.width(9);
      at(np).print(stream);
      stream << "  basis_index: " << at(np).basis_ind() << "  clust_index: " << at(np).nlist_ind() << " ";
      if(delim)
        stream << delim;
    }
    stream << "\n"
           << "            Basis Functions:\n";
    BasisSet tbasis(clust_basis);
    tbasis.accept(OccFuncLabeler("\\phi_%b_%f(s_%n)"));
    for(Index i = 0; i < tbasis.size(); i++) {
      stream << "              \\Phi_" << begin_ind + i << " = " << tbasis[i]->tex_formula() << std::endl;
    }
  }

  //*************************************************

  SiteCluster operator*(const SymOp &LHS, const SiteCluster &RHS) {
    return SiteCluster(RHS).apply_sym(LHS);
  }


  jsonParser &SiteCluster::to_json(jsonParser &json) const {

    // class SiteCluster : public GenericCluster<Site>
    GenericCluster<Site>::to_json(json);

    //BasisSet ccd_basis;
    //if(ccd_basis.size() > 0)
    //json["ccd_basis"] = ccd_basis;

    // BasisSet clust_basis;
    if(clust_basis.size() > 0)
      json["clust_basis"] = clust_basis;

    // Array<Tensor<double> >occupation_basis_tensors;
    //if(occupation_basis_tensors.size() > 0)
    //json["occupation_basis_tensors"] = occupation_basis_tensors;

    // Array<double> eci_coeffs;
    //if(eci_coeffs.size() > 0)
    //json["eci_coeffs"] = eci_coeffs;

    return json;
  }

  void SiteCluster::from_json(const jsonParser &json) {
    try {

      // class SiteCluster : public GenericCluster<Site>
      GenericCluster<Site>::from_json(json);

      // No reading BasisSet yet.
      //
      // //BasisSet ccd_basis;
      // if(json.contains("ccd_basis")) {
      //   from_json(ccd_basis, json["ccd_basis"]);
      // }
      //
      // // BasisSet clust_basis;
      // if(json.contains("clust_basis")) {
      //   from_json(clust_basis, json["clust_basis"]);
      // }


      // Array<Tensor<double> >occupation_basis_tensors;
      //if(json.contains("occupation_basis_tensors")) {
      //CASM::from_json(occupation_basis_tensors, json["occupation_basis_tensors"]);
      //}

      // Array<double> eci_coeffs;
      //if(json.contains("eci_coeffs")) {
      //CASM::from_json(eci_coeffs, json["eci_coeffs"]);
      //}
    }
    catch(...) {
      /// re-throw exceptions
      throw;
    }
  }

  jsonParser &to_json(const SiteCluster &clust, jsonParser &json) {
    return clust.to_json(json);
  }
  void from_json(SiteCluster &clust, const jsonParser &json) {
    try {
      clust.from_json(json);
    }
    catch(...) {
      /// re-throw exceptions
      throw;
    }
  }

  namespace SiteCluster_impl {

    BasisSet construct_clust_dof_basis(SiteCluster const &_clust, std::vector<BasisSet const *> const &site_dof_sets) {
      BasisSet result;
      result.set_dof_IDs(_clust.nlist_inds());
      Array<SymGroupRep const *> subspace_reps;
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

