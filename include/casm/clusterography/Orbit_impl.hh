#include <climits>

#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/Site.hh"
#include "casm/symmetry/SymMatrixXd.hh"
#include "casm/clusterography/SiteCluster.hh"

namespace CASM {
  //**************************************************************************************

  template<typename ClustType>
  void GenericOrbit<ClustType>::set_lattice(const Lattice &new_home, COORD_TYPE mode) {
    prototype.set_lattice(new_home, mode);

    for(Index i = 0; i < size(); i++)
      at(i).set_lattice(new_home, mode);

    for(Index i = 0; i < equivalence_map.size(); i++)
      for(Index j = 0; j < equivalence_map[i].size(); j++)
        equivalence_map[i][j].set_lattice(new_home, mode);

    return;

  }


  //****************************************************************************************************
  //void get_equivalent(const Array<SymOp> &sym_group);
  //Take prototype cluster, find all of the equivalent clusters (fills orbit array)
  //Also, fill equivalence_map

  //Apply each SymOp to prototype, check to see if result is already in the orbit array
  // if(contains(prototype.apply_sym(sym_group.get(i))) do:
  //If it is not in the orbit array, add it, and record the SymOp in equivalence_map
  //If it is already in the orbit array, record the SymOp that mapped the prototype onto that cluster

  template<typename ClustType>
  void GenericOrbit<ClustType>::get_equivalent(const Array<SymOp> &sym_group) {
    if(sym_group.size() == 0) {
      std::cerr << "WARNING: In Orbit::get_equivalent, sym_group must at least have one element (identity).\n";
      assert(0);
    }

    Index i, map_ind;
    ClustType t_cluster(prototype.get_home());
    Coordinate map_shift(Vector3<double>(0, 0, 0), prototype.get_home());
    Coordinate within_shift(map_shift);
    clear();
    equivalence_map.clear();

    prototype.within();
    prototype.prepare_prototype();
    // first put prototype into the orbit and preserve open space for sym_op(s) that map prototype onto itself
    // doing resize, basically erases all work history in the past
    this->push_back(prototype);
    equivalence_map.resize(1);

    // go through all symmetry operation in sym_group
    for(i = 0; i < sym_group.size(); i++) {
      t_cluster = sym_group[i] * prototype; // apply current sym_op to prototype, using the overloaded *

      map_ind = find(t_cluster, map_shift);

      //Case 1: Cluster already exists in orbit; update equivalence map to record
      //        symmetry operation that maps prototype onto equivalent cluster
      if(map_ind < size()) {
        equivalence_map[map_ind].push_back(SymOp(map_shift)*sym_group[i]);
      }

      //Case 2: Identified new equivalent cluster; add it to orbit and update equivalence map
      //        to record symmetry operation that maps prototype onto equivalent cluster
      else {
        t_cluster.within(0, within_shift); 	//make sure new equivalent cluster has first site within primitive cell
        this->push_back(t_cluster);   // add it
        SymOp tsymm(SymOp(map_shift + within_shift)*sym_group[i]);
        equivalence_map.push_back(Array<SymOp>(1, tsymm));
      }

    }
    return;
  }

  //****************************************************************************************************

  template<typename ClustType>
  GenericOrbit<ClustType> &GenericOrbit<ClustType>::apply_sym(const SymOp &op) {
    prototype.apply_sym(op);
    for(Index i = 0; i < size(); i++) {
      at(i).apply_sym(op);
      if(equivalence_map.size() < size()) continue;
      for(Index j = 0; j < equivalence_map[i].size(); j++) {
        equivalence_map[i][j].apply_sym(op);
      }
    }
    return *this;
  }

  //****************************************************************************************************
  //bool contains(const Cluster &test_clust);
  //go through all cluster in orbit array and see if test_clust is among them, accounting for translation
  //HINT:  only translate if the periodicity flag is on

  template<typename ClustType>
  bool GenericOrbit<ClustType>::contains(const ClustType  &test_clust) const {
    ClustType  tclust(test_clust);
    for(Index i = 0; i < size(); i++) {
      if(tclust.map_onto(at(i))) {
        return true;
      }
    }

    return false;
  }

  //****************************************************************************************************
  //bool contains(const Site &test_site);
  //go through all cluster in orbit array and see if test_site is contained in any of them, accounting for translation

  template<typename ClustType>
  bool GenericOrbit<ClustType>::contains(const typename ClustType::WhichCoordType &test_coord) const {

    for(Index i = 0; i < size(); i++) {
      if(at(i).contains_periodic(test_coord))
        return true;
    }

    return false;
  }

  //****************************************************************************************************
  //int find(const ClustType  &test_clust);
  //go through all cluster in orbit array and return the index of test_clust, accounting for translation
  //HINT:  only translate if the periodicity flag is on

  template<typename ClustType>
  Index GenericOrbit<ClustType>::find(const ClustType  &test_clust) const {
    Index i;
    ClustType  tclust(test_clust);
    for(i = 0; i < size(); i++) {
      if(tclust.map_onto(at(i)))
        break;
    }
    return i;
  }

  //****************************************************************************************************
  //int find(const ClustType  &test_clust);
  //go through all cluster in orbit array and return the index of test_clust, accounting for translation
  //HINT:  only translate if the periodicity flag is on

  template<typename ClustType>
  Index GenericOrbit<ClustType>::find(const ClustType  &test_clust, Coordinate &trans) const {
    Index i;
    trans(FRAC) = Vector3<double>(0, 0, 0);
    ClustType  tclust(test_clust);
    for(i = 0; i < size(); i++) {
      if(tclust.map_onto(at(i), trans)) {
        break;
      }
    }
    return i;
  }


  //********************************************************************

  template<typename ClustType>
  void GenericOrbit<ClustType>::collect_basis_info(const BasicStructure<Site> &struc, const Coordinate &shift) {
    collect_basis_info(struc.basis, shift);
  }

  //********************************************************************

  template<typename ClustType>
  void GenericOrbit<ClustType>::collect_basis_info(const BasicStructure<Site> &struc) {
    collect_basis_info(struc.basis);
  }

  //********************************************************************

  template<typename ClustType>
  void GenericOrbit<ClustType>::collect_basis_info(const Array<typename ClustType::WhichCoordType> &basis, const Coordinate &shift) {
    prototype.collect_basis_info(basis, shift);
    for(Index ne = 0; ne < size(); ne++)
      at(ne).collect_basis_info(basis, shift);
  }

  //********************************************************************

  template<typename ClustType>
  void GenericOrbit<ClustType>::collect_basis_info(const Array<typename ClustType::WhichCoordType> &basis) {
    prototype.collect_basis_info(basis);
    for(Index ne = 0; ne < size(); ne++)
      at(ne).collect_basis_info(basis);
  }

  //********************************************************************

  template<typename ClustType>
  void GenericOrbit<ClustType>::get_cluster_symmetry() {
    prototype.clust_group.clear();

    for(Index i = 0; i < equivalence_map[0].size(); i++) {
      prototype.clust_group.push_back(equivalence_map[0][i]);
    }
    prototype.get_permute_group();

    for(Index i = 0; i < size(); i++) {
      at(i).clust_group = prototype.clust_group;
      if(i > 0)
        at(i).clust_group.apply_sym(equivalence_map[i][0]);
      at(i).get_permute_group();//JCT - I think this is correct; the relation between equivalence_map and permutations is ill-defined
    }
  }

  //********************************************************************

  /// get permutation representation of every operation in equivalence_map to describe how operations permute site order
  template<typename ClustType>
  SymGroupRep const *GenericOrbit<ClustType>::get_full_permutation_representation() {

    if(!equivalence_map.size())
      return NULL;

    if(permute_rep_ID != -1) return (equivalence_map[0][0].master_group().representation(permute_rep_ID));

    SymGroupRep *tRep(new SymGroupRep((equivalence_map[0][0].master_group())));

    tRep->resize(equivalence_map[0][0].master_group().size(), NULL);

    Array<Index> iperm;
    //This is a hack for doing cubic-tetragonal
    if(equivalence_map.size() > 1) {
      //std::cout << "*********** equivalence_map.size() is " << equivalence_map.size() << '\n';
      equivalence_map[1].swap_elem(equivalence_map[1].size() - 1, 0);
      //std::cout << "I think this matrix is inversion:\n" << equivalence_map[1][0].get_matrix(CART) << "\n";
      at(1) = equivalence_map[1][0] * prototype;
    }
    for(Index i = 0; i < equivalence_map.size(); i++) {
      for(Index j = 0; j < equivalence_map[i].size(); j++) {
        //std::cout << "i= " << i << " and j= " << j << std::endl;
        iperm.clear();
        (equivalence_map[i][j]*prototype).find(at(i), iperm);
        //std::cout << "After mapping, iperm is " << iperm << std::endl;
        tRep->at(equivalence_map[i][j].index()) = new SymPermutation(iperm);
        //std::cout << "inserted new permutations in tRep" << std::endl;
      }
    }

    permute_rep_ID = equivalence_map[0][0].master_group().add_representation(tRep);
    return tRep;
  }

  //********************************************************************

  /// get permutation representation of every operation in equivalence_map to describe how operations permute site order
  template<typename ClustType>
  SymGroupRep const *GenericOrbit<ClustType>::get_full_coord_representation() {

    if(!equivalence_map.size()) return NULL;

    if(coord_rep_ID != -1) return equivalence_map[0][0].master_group().representation(coord_rep_ID);

    if(permute_rep_ID == -1) get_full_permutation_representation();

    if(equivalence_map[0][0].master_group().get_coord_rep_ID() == -1) equivalence_map[0][0].master_group().add_coord_rep();

    coord_rep_ID = equivalence_map[0][0].master_group().add_kronecker_rep(permute_rep_ID, equivalence_map[0][0].master_group().get_coord_rep_ID());

    return equivalence_map[0][0].master_group().representation(coord_rep_ID);

  }

  //********************************************************************
  /**
   * Reads the orbit in the specified mode.
   * If read_tensors is set to true, then the Cluster read will read
   * in the tensor basis.
   *
   * @param stream Input file stream
   * @param mode  Cartesian or fractional mode
   * @param read_tensors Flag determining whether tensors should be read
   *        in or not
   */
  //********************************************************************

  template<typename ClustType>
  void GenericOrbit<ClustType>::read(std::istream &stream, COORD_TYPE mode, const SymGroup &sym_group, bool read_tensors) {

    int num_clust;
    char ch;

    ch = stream.peek();
    while((ch != 'C') && (ch != 'c') && !stream.eof()) {
      stream.ignore(1000, '\n');
      ch = stream.peek();
    }

    stream.ignore(1000, ':');
    ch = stream.peek();
#ifdef DEBUG
    std::cout << "ch in cluster is " << ch << "\n";
#endif //DEBUG
    while(ch < '0' || ch > '9') {
      stream.ignore(1, '\n');
      ch = stream.peek();
    }
    stream >> num_clust;

#ifdef DEBUG
    std::cout << "num_clust is " << num_clust << "\n";
#endif //DEBUG

    for(int i = 0; i < num_clust; i++) {
      if(i == 0) {
        prototype.clear();
        prototype.read(stream, mode, read_tensors);
        (*this).push_back(prototype);
      }
      else {
        SiteCluster tclust(prototype.get_home());
        (*this).push_back(tclust);
        (*this).back().read(stream, mode, read_tensors);
      }
    }
    //    }

    //if ( ((*this).back().max_length() > TOL) ||
    //	 (((*this).back().max_length() < TOL) && ((*this).back().size() == 1)) ){
    //This works only when reading in orbitree
    if((*this).back().size() > 0) {
#ifdef DEBUG
      std::cout << "Inside this loop the cluster is \n";
      std::cout << (*this).back() << "\n";
#endif //DEBUG
      get_equivalent(sym_group);
      //}
    }

    for(Index t = 0; t < size(); t++) {
      at(t).within();
    }

    //get_cluster_symmetry();

    return;
  }

  //********************************************************************

  // Divide by multiplicity. Same result as evaluating correlations via orbitree.
  template<typename ClustType>
  ReturnArray<std::string> GenericOrbit<ClustType>::orbit_function_cpp_strings(const Array<FunctionVisitor *> &labelers) {
    std::string suffix("");
    Array<std::string> formulae(prototype.clust_basis.size(), std::string());
    if(size() > 1) {
      formulae.resize(prototype.clust_basis.size(), std::string("("));
      suffix = ")/" + std::to_string(size()) + ".0";
    }

    for(Index ne = 0; ne < size(); ne++) {
      for(Index nl = 0; nl < labelers.size(); nl++)
        at(ne).clust_basis.accept(*labelers[nl]);
    }

    for(Index nf = 0; nf < prototype.clust_basis.size(); nf++) {
      for(Index ne = 0; ne < size(); ne++) {
        if(!at(ne).clust_basis[nf] || (at(ne).clust_basis[nf]->formula()) == "0")
          continue;
        if(formulae[nf].size() > 1)
          formulae[nf] += " + ";
        formulae[nf] += "(" + at(ne).clust_basis[nf]->formula() + ")";
      }

      if(formulae[nf].size() <= 1)
        formulae[nf].clear();
      else
        formulae[nf] += suffix;
    }
    return formulae;
  }

  //********************************************************************

  /// b_index is the basis site index, f_index is the index of the configurational site basis function in Site::occupant_basis
  template<typename ClustType>
  ReturnArray<std::string>  GenericOrbit<ClustType>::flower_function_cpp_strings(const Array<FunctionVisitor *> &labelers, Index b_index) {
    Array<std::string> formulae(prototype.clust_basis.size(), std::string());
    std::string suffix;
    Index ib;

    //normalize by multiplicity (by convention)
    if(size()*prototype.size() > 1) {
      formulae.resize(prototype.clust_basis.size(), std::string("("));
      suffix = ")/" + std::to_string(size()) + ".0";
    }

    for(Index ne = 0; ne < size(); ne++) {
      //std::cout << "# of translists: " << at(ne).trans_nlists().size() << "\n";
      for(Index nt = 0; nt < at(ne).trans_nlists().size(); nt++) {
        ib = at(ne).trans_nlist(nt).find(b_index);

        //std::cout << "ib is " << ib << " of " << at(ne).size() << "\n";
        if(ib == at(ne).size())
          continue;

        at(ne).set_nlist_inds(at(ne).trans_nlist(nt));

        for(Index nl = 0; nl < labelers.size(); nl++)
          at(ne).clust_basis.accept(*labelers[nl]);

        for(Index nf = 0; nf < at(ne).clust_basis.size(); nf++) {
          if(!at(ne).clust_basis[nf] || (at(ne).clust_basis[nf]->formula()) == "0")
            continue;

          if(formulae[nf].size() > 1)
            formulae[nf] += " + ";
          formulae[nf] += "(" + at(ne).clust_basis[nf]->formula() + ")";
        }
      }
    }

    // Make sure that formulae that evaluate to zero have an empty string.
    for(Index nf = 0; nf < prototype.clust_basis.size(); nf++) {
      if(formulae[nf].size() <= 1)
        formulae[nf].clear();
      else
        formulae[nf] += suffix;
    }
    return formulae;
  }

  //********************************************************************

  /// b_index is the basis site index, f_index is the index of the configurational site basis function in Site::occupant_basis
  template<typename ClustType>
  ReturnArray<std::string>  GenericOrbit<ClustType>::delta_occfunc_flower_function_cpp_strings(const Array<FunctionVisitor *> &labelers, Index b_index, Index f_index) {
    Array<std::string> formulae(prototype.clust_basis.size(), std::string());
    std::string suffix;
    Index ib;

    //normalize by multiplicity (by convention)
    if(size()*prototype.size() > 1) {
      formulae.resize(prototype.clust_basis.size(), std::string("("));
      suffix = ")/" + std::to_string(size()) + ".0";
    }


    for(Index ne = 0; ne < size(); ne++) {
      //std::cout << " **** for ne = " << ne << ":\n";
      for(Index nt = 0; nt < at(ne).trans_nlists().size(); nt++) {
        ib = at(ne).trans_nlist(nt).find(b_index);
        if(ib == at(ne).size())
          continue;
        //std::cout << " **** for translist: " << at(ne).trans_nlist(nt) << ":\n";
        at(ne).set_nlist_inds(at(ne).trans_nlist(nt));
        BasisSet quotient_basis = at(ne).clust_basis.poly_quotient_set(at(ne)[ib].occupant_basis()[f_index]);
        for(Index nf = 0; nf < quotient_basis.size(); nf++) {
          if(!quotient_basis[nf])
            continue;
          for(Index nl = 0; nl < labelers.size(); nl++)
            quotient_basis[nf]->accept(*labelers[nl]);

          if((quotient_basis[nf]->formula()) == "0")
            continue;

          if(formulae[nf].size() > 1)
            formulae[nf] += " + ";
          formulae[nf] += "(" + quotient_basis[nf]->formula() + ")";
        }
      }
    }

    // Make sure that formulae that evaluate to zero have an empty sting.
    for(Index nf = 0; nf < prototype.clust_basis.size(); nf++) {
      if(formulae[nf].size() <= 1)
        formulae[nf].clear();
      else
        formulae[nf] += suffix;
    }

    return formulae;
  }

  //********************************************************************
  /**
   * Calculates the ECI's of the prototype cluster, then calculates
   * the ECI of the clusters in the orbit.  Note that clusters in the
   * same orbit should have the same tensor coefficients.
   *
   * @param[in] rank The rank of the tensors.
   */
  //********************************************************************

  template<typename ClustType>
  void GenericOrbit<ClustType>::calc_eci(Index rank) {

    prototype.eci.redefine(Array<Index>(rank, 3));

    //loop over tensor_basis of each prototype
    for(Index i = 0; i < prototype.tensor_basis.size(); i++) {
      if(std::isnan(prototype.tensor_basis.eci(i))) {
        std::cerr << "WARNING: One of your ECI's have not been entered in the input file. \n"
                  << "It has been set to 0. \n";
        prototype.tensor_basis.eci(i) = 0;
      }
      prototype.eci += (prototype.tensor_basis.eci(i) * prototype.tensor_basis[i]);
    }

    for(Index j = 0; j < size(); j++) {
      at(j).eci.redefine(Array<Index>(rank, 3));

      for(Index k = 0; k < at(j).tensor_basis.size(); k++) {
        if(std::isnan(prototype.tensor_basis.eci(k))) {
          std::cerr << "WARNING: One of your ECI's have not been entered in the input file. \n"
                    << "It has been set to 0. \n";
          prototype.tensor_basis.eci(k) = 0;
        }
        at(j).eci += (prototype.tensor_basis.eci(k) * at(j).tensor_basis[k]);
        at(j).tensor_basis.eci(k) = prototype.tensor_basis.eci(k);
      }

    }
  }

  //********************************************************************
  /**
   * Gets the tensor basis for each prototype and maps it onto the
   * equivalent clusters.
   *
   * Gets the tensor basis for the prototype cluster before using
   * equivalence_map to symmetrically transform the tensor basis of the
   * clusters in the orbit.
   *
   * @param[in] rank The rank of the tensors in the tensor basis.
   */
  //********************************************************************

  template<typename ClustType>
  void GenericOrbit<ClustType>::get_tensor_basis(Index rank) {

    prototype.get_tensor_basis(rank);

    // looping over clusters in orbit
    for(Index i = 0; i < size(); i++) {
      at(i).tensor_basis = prototype.tensor_basis;
      at(i).tensor_basis.apply_sym(equivalence_map[i][0]);
    }
  }

  //********************************************************************

  template<typename ClustType>
  jsonParser &GenericOrbit<ClustType>::to_json(jsonParser &json) const {

    json.put_obj();

    // template<typename ClustType>
    // class GenericOrbit : public Array<ClustType >
    json["clusters"].put_array(size());
    for(Index i = 0; i < size(); i++) {
      json["clusters"][i] = at(i);
    }

    // mutable int index;
    json["index"] = index;

    // int permute_rep_ID, coord_rep_ID;
    // json["permute_rep_ID"] = permute_rep_ID;
    // json["coord_rep_ID"] = coord_rep_ID;

    // Array<GenericOrbit *> sub_cluster;
    // This doesn't seem used right now?

    // ClustType  prototype;
    json["prototype"] = prototype;

    // Array< Array<SymOp> > equivalence_map;
    //json["equivalence_map"] = equivalence_map;

    return json;
  }

  //********************************************************************

  /// Assumes the prototype lattice is already set
  template<typename ClustType>
  void GenericOrbit<ClustType>::from_json(const jsonParser &json) {
    try {

      // template<typename ClustType>
      // class GenericOrbit : public Array<ClustType >
      clear();
      ClustType clust(prototype);
      // // GOING TO GENERATE THE EQUIVALENT CLUSTERS INSTEAD OF READING
      // // THEM IN
      // std::cout<<"Number of clusters:"<<json["clusters"].size()<<std::endl;
      // this->resize(json["clusters"].size(), clust);
      // for(int i = 0; i < json["clusters"].size(); i++) {
      //   std::cout<<"Working on:"<<i;
      //   CASM::from_json(at(i), json["clusters"][i]);
      // }

      // mutable int index;
      CASM::from_json(index, json["index"]);

      // int permute_rep_ID, coord_rep_ID;
      // CASM::from_json(permute_rep_ID, json["permute_rep_ID"]);
      // CASM::from_json(coord_rep_ID, json["coord_rep_ID"]);

      // Array<GenericOrbit *> sub_cluster;
      // This doesn't seem used right now?

      // ClustType  prototype;
      //std::cout<<"Reading in the prototype cluster"<<std::endl;
      CASM::from_json(prototype, json["prototype"]);
      // std::cout<<"Generating the equivalent clusters"<<std::endl;

      // // Array< Array<SymOp> > equivalence_map;
      // SymOp op(prototype.get_home());
      // Array<SymOp> equiv;
      // equivalence_map.clear();
      // std::cout<<"Reading in the equivalence_map"<<std::endl;
      // std::cout<<"i<"<<json["equivalence_map"].size()<<std::endl;
      // for(int i = 0; i < json["equivalence_map"].size(); i++) {
      //   equiv.clear();
      //   std::cout<<"   i:"<<i<<std::endl;
      //   std::cout<<"j<"<<json["equivalence_map"][i].size()<<std::endl;
      //   for(int j = 0; j < json["equivalence_map"][i].size(); j++) {
      //     std::cout<<"       j:"<<j;
      //     CASM::from_json(op, json["equivalence_map"][i][j]);
      //     equiv.push_back(op);
      //   }
      //   equivalence_map.push_back(equiv);
      // }

      //std::cout<<"Done with reading in the Orbit"<<std::endl;
    }
    catch(...) {
      /// re-throw exceptions
      throw;
    }
  }


}


