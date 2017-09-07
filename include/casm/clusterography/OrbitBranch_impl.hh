namespace CASM {
  //********************************************************************

  template<typename ClustType>
  void GenericOrbitBranch<ClustType>::clear() {
    Array<GenericOrbit<ClustType> > :: clear();
    pivot.clear();
    index.clear();
    m_num_sites = 0;
  }

  //********************************************************************

  template<typename ClustType>
  const GenericOrbit<ClustType> &GenericOrbitBranch<ClustType>::orbit(Index no) const {
    return at(no);
  }

  //********************************************************************

  template<typename ClustType>
  GenericOrbit<ClustType> &GenericOrbitBranch<ClustType>::orbit(Index no) {
    return at(no);
  }

  //********************************************************************

  template<typename ClustType>
  const ClustType &GenericOrbitBranch<ClustType>::prototype(Index no) const {
    return at(no).prototype;
  }

  //********************************************************************

  template<typename ClustType>
  ClustType &GenericOrbitBranch<ClustType>::prototype(Index no) {
    return at(no).prototype;
  }

  //********************************************************************

  template<typename ClustType>
  const ClustType &GenericOrbitBranch<ClustType>::equiv(Index no, Index ne) const {
    return at(no).at(ne);
  }

  //********************************************************************

  template<typename ClustType>
  ClustType &GenericOrbitBranch<ClustType>::equiv(Index no, Index ne) {
    return at(no).at(ne);
  }

  //********************************************************************

  template<typename ClustType>
  Index GenericOrbitBranch<ClustType>::size(Index no) const {
    return at(no).size();
  }

  //********************************************************************

  template<typename ClustType>
  void GenericOrbitBranch<ClustType>::set_lattice(const Lattice &new_lat, COORD_TYPE mode) {
    for(Index no = 0; no < size(); no++)
      orbit(no).set_lattice(new_lat, mode);
    pivot.set_lattice(new_lat, mode);
    return;
  }

  //*******************************************************************************************

  template<typename ClustType>
  void GenericOrbitBranch<ClustType>::set_pivot(const ClustType &new_pivot) {
    if(&(pivot.home()) != &(new_pivot.home())) {
      std::cerr << "WARNING!!! In GenericOrbitBranch::set_pivot(), the new 'pivot' clust has a different home lattice than the old pivot cluster.\n"
                << "           This may result in unexpected behavior!\n";
      assert(0);
      set_lattice(new_pivot.home(), CART); // <-- this makes it slightly safer
      pivot = new_pivot;
    }

    pivot = new_pivot;

  }
  //********************************************************************

  template<typename ClustType>
  void GenericOrbitBranch<ClustType>::push_back(const GenericOrbit<ClustType> &new_orbit) {
    if(num_sites() < 0)
      m_num_sites = new_orbit.prototype.size();

    else if(num_sites() != new_orbit.prototype.size()) {
      std::cerr << "WARNING:  Trying to add " << new_orbit.prototype.size()
                << "-site Orbit to an OrbitBranch intended for " << num_sites()
                << "-site Orbits.  Continuing...\n";
      return;
    }

    Array< GenericOrbit<ClustType> >::push_back(new_orbit);
    return;
  }

  //********************************************************************

  template<typename ClustType>
  GenericOrbitBranch<ClustType> &GenericOrbitBranch<ClustType>::apply_sym(const SymOp &op) {
    for(Index no = 0; no < size(); no++)
      at(no).apply_sym(op);
    pivot.apply_sym(op);
    return *this;
  }

  //********************************************************************
  template<typename ClustType>
  void GenericOrbitBranch<ClustType>::print(std::ostream &stream, COORD_TYPE mode) {
    /*Mainly will be used for print out the asymmetric
     * unit, however will print any orbitbranch.
     * First, loop over the outer branch, the size
     * of which is the number of atoms in the asymmetric unit.
     * Then, loop over each point in this branch to print out
     * the symmetrically equivalent points
     */
    stream << "-----------------------------------------------\n";
    for(Index i = 0; i < size(); i++) {
      stream << "Asymmetric Point # " << i + 1 << " of " << size() << ":\n";
      /*Print the unique point (first item in the branch we're currently in)
       * at(i).at(0) points to this, which should be a cluster... How do we just print the cluster?
       * GenericOrbitBranch<Site>, Array<SiteOrbit<Site>> --> Array< Array<GenericCluster<Site>> --> Array<Array<Array<Site>>>
       */
      at(i).at(0).at(0).print_occ(stream); //Since we only have one item in the cluster, can print the 1st thing in it (the second at(0) ) to JUST print the site info, not the cluster info

      stream << "\n\nEquivalent points are: \n";

      for(Index j = 1; j < at(i).size(); j++) { //Loop over all the equivalent points and print them
        stream << " "; //Indent all of the equivalent points
        at(i).at(j).at(0).print_occ(stream);
        stream << "\n";
      }

      stream << "\n";
      stream << "-----------------------------------------------\n";

    }
  }

  //********************************************************************

  template<typename ClustType>
  void GenericOrbitBranch<ClustType>::sort() {
    for(Index i = 0; i < size(); i++) { //Loops over the inner array, which is the subset of all the same-sized clusters
      for(Index j = i + 1; j < size(); j++) { //Takes the next-in-line of the  inner array and compares it
        //	  if(orbit(i).prototype.max_length > orbit(j).prototype.max_length){
        if(std::abs(orbit(i).max_length() - orbit(j).max_length()) < TOL && orbit(i).min_length() > orbit(j).min_length()) {
          GenericOrbit<ClustType> torbit(orbit(j));
          orbit(j) = orbit(i);
          orbit(i) = torbit;
        }
        else if(orbit(i).max_length() > orbit(j).max_length() + TOL) {
          GenericOrbit<ClustType> torbit(orbit(j));
          orbit(j) = orbit(i);
          orbit(i) = torbit;
        }
      }
    }
  }
  //*******************************

  template<typename ClustType>
  Index GenericOrbitBranch<ClustType>::find(const ClustType &test_clust, double tol) const {
    for(Index i = 0; i < size(); i++) {
      if(orbit(i).contains(test_clust, tol))
        return i;
    }
    return size();
  }

  //*******************************

  template<typename ClustType>
  bool GenericOrbitBranch<ClustType>::contains(const ClustType &test_clust, double tol) const {
    if(size() == 0) return false;

    if(find(test_clust, tol) < size()) {
      return true;
    }
    return false;

  }

  //***********************************************************

  template<typename ClustType>
  void GenericOrbitBranch<ClustType>::generate_asymmetric_unit(const Array<typename ClustType::WhichCoordType > &basis, const SymGroup &factor_group, double tol) {
    m_num_sites = 1;
    if(size() != 0) { //Added by Ivy 11/19/13

      std::cerr << "WARNING in GenericOrbitBranch::generate_asymmetric_unit()\n"
                << "This OrbitBranch is not empty\n"
                << "and is about to be overwritten!\n";

      clear();

    }

    if(basis.size() == 0)
      return;
    set_lattice(basis[0].home(), FRAC);
    ClustType tclust(basis[0].home());
    Index nb, no;


    //Loop over basis sites
    for(nb = 0; nb < basis.size(); nb++) {

      tclust.clear();
      //Make point cluster for each site
      tclust.push_back(basis[nb]);
      tclust.within();

      //Check if that point is among the asymmetric points already found
      for(no = 0; no < size(); no++) {
        if(at(no).contains(tclust, tol)) {
          break;
        }
      }


      //If it is not found, no is equal to asym_unit.size()
      //add this new point cluster as a point orbit to asym_unit
      //and get symmetrically equivalent points to fill orbit
      if(no >= size()) {
        push_back(GenericOrbit<ClustType >(tclust));
        back().get_equivalent(factor_group, tol);
      }
    }
    for(no = 0; no < size(); no++) {
      at(no).collect_basis_info(basis);
    }

    //std::cout << "finish generate_asymmeric_unit()" << std::endl;
    return;
  }


  //*******************************
  // Extracts the orbits that include cluster 'pivot', such that each orbit is a 'flower' of clusters
  // with 'pivot' at the center

  template<typename ClustType>
  bool GenericOrbitBranch<ClustType>::extract_orbits_including(const ClustType &pivot, GenericOrbitBranch<ClustType> &flowerbranch, double tol) const {
    Index no, ne;
    ClustType tclust(pivot.home());
    bool found_any = false;

    for(no = 0; no < size(); no++) {
      GenericOrbit<ClustType> torbit(tclust); //TODO check constructor
      torbit.set_index(orbit(no).get_index());
      for(ne = 0; ne < size(no); ne++) {
        //std::cout << "(no,ne) = " << no << ", " << ne << '\n';
        tclust = equiv(no, ne);
        int tmaps(0);
        //Repeatedly try to map the clusters onto the subcluster until it fails. In the process,
        //map_onto_sublcuster applies the transformation to the clusters, yielding new translational
        //clusters that we want.
        while(tclust.map_onto_subcluster(pivot, tmaps)) {
          tmaps++;
          //turn off periodicity, so Orbit::contains preserves locality of pivot
          PERIODICITY_MODE temp_mode(LOCAL);
          if(!torbit.contains(tclust, tol)) {
            if(torbit.size() == 0) {
              torbit.prototype = tclust;
            }
            torbit.push_back(tclust);

            //We use get_equivalent, because it is efficient, and because it produces equivalence_map
            //flowerbranch.back().get_equivalent(pivot.clust_group);
            found_any = true;
          }
        }
      }
      if(torbit.size() != 0) {
        flowerbranch.push_back(torbit);
      }
    }
    return found_any;
  }


  //********************************************************************

  template<typename ClustType>
  jsonParser &GenericOrbitBranch<ClustType>::to_json(jsonParser &json) const {
    //    std::cout<<"In GenericOrbitBranch<ClustType>::to_json"<<std::endl;
    json.put_obj();

    // template<typename ClustType>
    // class GenericOrbitBranch : public Array< GenericOrbit<ClustType> >
    json["orbits"].put_array(size());
    for(Index i = 0; i < size(); i++) {
      json["orbits"][i] = at(i);
    }

    // int m_num_sites;
    json["m_num_sites"] = m_num_sites;

    // ClustType pivot;
    json["pivot"] = pivot;

    // Array<int> index;
    json["index"] = index;

    return json;
  }

  //********************************************************************

  /// Assumes the pivot lattice is already set
  template<typename ClustType>
  void GenericOrbitBranch<ClustType>::from_json(const jsonParser &json) {
    try {

      // template<typename ClustType>
      // class GenericOrbitBranch : public Array< GenericOrbit<ClustType> >
      //std::cout<<"Number of Orbits:"<<json["orbits"].size()<<std::endl;
      GenericOrbit<ClustType> orbit(pivot);
      this->resize(json["orbits"].size(), orbit);
      for(int i = 0; i < json["orbits"].size(); i++) {
        //std::cout<<"Working on orbit:"<<i<<std::endl;
        CASM::from_json(at(i), json["orbits"][i]);
      }

      //std::cout<<"Reading in m_num_sites"<<std::endl;
      // int m_num_sites;
      CASM::from_json(m_num_sites, json["m_num_sites"]);

      //std::cout<<"Reading in pivot"<<std::endl;
      // ClustType pivot;
      CASM::from_json(pivot, json["pivot"]);

      //std::cout<<"Reading in index"<<std::endl;
      // Array<int> index;
      CASM::from_json(index, json["index"]);


    }
    catch(...) {
      /// re-throw exceptions
      throw;
    }
  }

}

