#include "casm/crystallography/Lattice.hh"
#include "casm/symmetry/SymPermutation.hh"

namespace CASM {


  //****************************************************
  template <typename CoordType>
  GenericCluster<CoordType>::GenericCluster(const Lattice &init_home) :
    home(&init_home), min_length_val(0.0), max_length_val(0.0),
    DoF_rep(-1), clust_group(LOCAL), permute_group(SymGroupRep::NO_HOME, 0) {

    //nothing else to do for now
  };

  //****************************************************

  template <typename CoordType>
  void GenericCluster<CoordType>::set_lattice(const Lattice &new_home, COORD_TYPE mode) {
    if(home == &new_home)
      return;
    home = &new_home;
    for(Index i = 0; i < size(); i++)
      at(i).set_lattice(new_home, mode);

    return;
  }

  //****************************************************
  template <typename CoordType>
  void GenericCluster<CoordType>::update_data_members(const BasicStructure<CoordType> &ref_struc) {
    Index basis_ind;
    for(Index i = 0; i < (*this).size(); i++) {
      basis_ind = ref_struc.find(at(i));
      if(basis_ind >= ref_struc.basis.size()) {
        std::cerr << "ERROR in GenericCluster<CoordType>::set_occupant_basis. Could not "
                  << "find a basis site in ref_struc that corresponds to site " << i
                  << "QUITTING" << std::endl;
        std::cerr << "DEBUGGING info: " << std::endl;
        print(std::cerr);
        //std::cerr<<"The structure is:"<<std::endl;
        //ref_struc.print(std::cerr);
        exit(666);
      }
      (*this).at(i).update_data_members(ref_struc.basis[basis_ind]);
    }
    return;
  }

  //********************************************

  template <typename CoordType>
  GenericCluster<CoordType> &GenericCluster<CoordType>::permute(const Array<Index> &iperm) {
    if(iperm.size() != size()) {
      std::cerr << "WARNING: Attempted to permute points of cluster of size  " << size() << " using a permutation of " << iperm.size() << " indicies.\n";
      return *this;
    }
    Array<CoordType>::permute(iperm);

    //Following lines will have to be updated as we add cluster basis, etc.
    for(Index i = 0; i < tensor_basis.size(); i++)
      tensor_basis[i].dim_permute(iperm);
    return *this;
  }

  //********************************************

  template <typename CoordType>
  GenericCluster<CoordType> &GenericCluster<CoordType>::permute(const Permutation &perm) {
    return permute(perm.perm_array());
  }

  //********************************************
  template <typename CoordType>
  GenericCluster<CoordType> &GenericCluster<CoordType>::apply_sym(const SymOp &op) {
    for(Index i = 0; i < size(); i++)
      at(i).apply_sym(op);

    //apply symmetry to eci, which is a tensor.  This will need to be revised to handle non-3d tensors
    tensor_basis.apply_sym(op);
    eci.apply_sym(op);
    clust_group.apply_sym(op);
    return *this;
  }

  //********************************************

  template <typename CoordType>
  GenericCluster<CoordType> &GenericCluster<CoordType>::apply_sym_no_trans(const SymOp &op) {
    for(Index i = 0; i < size(); i++)
      at(i).apply_sym_no_trans(op);

    tensor_basis.apply_sym(op);
    eci.apply_sym(op);
    clust_group.apply_sym(op);
    return *this;
  }

  //********************************************
  // Assuming that we have the factor group (or some similar group),
  //what are the operations that leave the cluster unchanged
  // Find them, alter the translations, and store them in clust_group

  template <typename CoordType>
  void GenericCluster<CoordType>::get_clust_group(const SymGroup &super_group) {
    Coordinate trans(*home);
    clust_group.clear();
    permute_group.clear();
    Array<Index> iperm;
    for(Index ng = 0; ng < super_group.size(); ng++) {
      GenericCluster<CoordType> tclust(*this);
      tclust.apply_sym(super_group[ng]);
      if(tclust.map_onto(*this, trans)) {
        clust_group.push_back(SymOp::translation(trans.cart())*super_group[ng]);
        find(tclust, iperm);

        permute_group.push_back_copy(SymPermutation(iperm));
      }
    }

    return;
  }

  //********************************************

  template <typename CoordType>
  void GenericCluster<CoordType>::get_permute_group() {
    permute_group.clear();
    Array<Index> iperm(size());
    for(Index i = 0; i < clust_group.size(); i++) {
      for(Index j = 0; j < size(); j++)
        iperm[j] = find(clust_group[i] * at(j));
      permute_group.push_back_copy(SymPermutation(iperm));
    }

  }

  //********************************************

  template <typename CoordType>
  void GenericCluster<CoordType>::get_tensor_basis(Index nrank) {
    if(!clust_group.size()) {
      std::cerr << "WARNING: Attempting to find tensor basis of cluster \n";
      print(std::cerr);
      std::cerr << "but cluster point group symmetry has not been specified.  Continuing... \n";
      return;
    }

    tensor_basis.generate_basis(nrank, clust_group, permute_group);

    return;
  }

  //********************************************

  template <typename CoordType>
  void GenericCluster<CoordType>::within(Index pivot_ind) {

    if(!size() || PERIODICITY_MODE::IS_LOCAL()) return;

    if(pivot_ind >= size()) {
      std::cerr << "WARNING: Attempted to map cluster within unit cell with invalid pivot." << std::endl;
      return;
    }

    Coordinate tcoord(*home);
    at(pivot_ind).within(tcoord);
    for(Index i = 0; i < size(); i++) {
      if(i != pivot_ind)
        at(i) += tcoord;
    }

    return;
  }

  //********************************************

  template <typename CoordType>
  void GenericCluster<CoordType>::within(Index pivot_ind, Coordinate &trans) {

    if(!size() || PERIODICITY_MODE::IS_LOCAL()) return;

    if(pivot_ind >= size()) {
      std::cerr << "WARNING: Attempted to map cluster within unit cell with invalid pivot." << std::endl;
      return;
    }

    trans.set_lattice(*home, CART);
    at(pivot_ind).within(trans);
    for(Index i = 0; i < size(); i++) {
      if(i != pivot_ind)
        at(i) += trans;
    }

    return;
  }

  //********************************************

  template <typename CoordType>
  void GenericCluster<CoordType>::all_within() {
    // should all_within check the periodicity mode?
    // It will probably be used in situations where periodicity_mode is local
    // but where normal functionality is desired.
    // e.g., mapping local clusters onto a supercell configuration
    if(PERIODICITY_MODE::IS_LOCAL()) return;
    for(Index i = 0; i < size(); i++)
      at(i).within();
    return;
  }

  //********************************************
  // Returns true if a vector can be found connecting one site of the cluster
  // to a periodic image of the cluster that is shorter than at least one of the
  // vectors connecting two sites of the cluster, assuming the periodic boundaries
  // defined by Lattice given by 'cell'

  template <typename CoordType>
  bool GenericCluster<CoordType>::image_check(const Lattice &cell, int max_nV) const {
    Coordinate tcoord(*home);
    Index i, j;
    int tnV;
    for(i = 0; i < size(); i++) {
      for(j = i + 1; j < size(); j++) {
        tcoord = at(j) - at(i);
        tnV = tcoord.voronoi_number(cell);
        if(tnV > max_nV || tnV < 0)
          return true;
      }
    }
    return false;
  }

  //********************************************
  //Returns the indices of the supercells that are large enough, so
  //that the cluster does not see an image of itself
  template <typename CoordType>
  ReturnArray<int> GenericCluster<CoordType>::SuperScreener(Array<Lattice> &s_cells, const Array<SymOp> &symoplist) {
    Array<int> keepNums;
    for(Index i = 0; i < s_cells.size(); i++) {
      if(!image_check(s_cells[i])) {
        keepNums.push_back(i);
      }
    }

    for(Index i = 0; i < s_cells.size(); i++) {
      for(Index j = 0; j < symoplist.size(); j++) {
        if(!is_integer(s_cells[i].inv_lat_column_mat()*symoplist[j].matrix()*s_cells[i].lat_column_mat(), TOL)) {
          if(keepNums.contains(i)) {
            keepNums.remove(keepNums.find(i));
          }
        }
      }
    }

    return keepNums;
  }

  //********************************************

  template <typename CoordType>
  Coordinate GenericCluster<CoordType>::geometric_center() const {
    if(!size())
      return Coordinate(*home);

    Coordinate tcoord(at(0));
    for(Index i = 1; i < size(); i++) {
      tcoord += at(i);
    }
    tcoord.cart() /= size();
    return tcoord;
  }

  //********************************************

  template <typename CoordType>
  bool GenericCluster<CoordType>::contains(const GenericCluster<CoordType> &test_cluster) const {
    Index i, j;
    for(i = 0; i < test_cluster.size(); i++) {
      for(j = 0; j < size(); j++) {
        if(test_cluster[i] == at(j))
          break;
      }
      if(j == size())
        return false;
    }
    return true;
  }

  //********************************************

  template <typename CoordType>
  bool GenericCluster<CoordType>::contains_periodic(const CoordType &test_site) const {
    for(Index i = 0; i < size(); i++) {
      if(at(i).compare(test_site))
        return true;
    }
    return false;
  }
  //********************************************
  //returns true if test_cluster is a subcluster of cluster (*this)
  //Array 'index' holds indices of (*this) that yield test_cluster
  template <typename CoordType>
  bool GenericCluster<CoordType>::find(const GenericCluster<CoordType> &test_cluster, Array<Index> &index) const {
    Index i, j;
    index.clear();
    index.reserve(test_cluster.size());
    for(i = 0; i < test_cluster.size(); i++) {
      for(j = 0; j < size(); j++) {
        if(test_cluster[i] == at(j) && !index.contains(j))
          break;
      }
      if(j == size()) {
        index.clear();
        return false;
      }
      index.push_back(j);
    }
    return true;
  }

  //********************************************
  //This routine needs to be edited so that when sites are reordered,
  //the cluster basis functions are updated to reflect the change in order

  template <typename CoordType>
  bool GenericCluster<CoordType>::map_onto_subcluster(const GenericCluster<CoordType> &pivot) {
    Index i, j, tsize;

    if(!(get_home() == pivot.get_home())) {
      std::cerr << "WARNING in Cluster::map_onto_subcluser!!\n"
                << "You are trying to map a pivot onto a cluster with\n"
                << "a different lattice! \n";
    }


    if(PERIODICITY_MODE::IS_LOCAL() && size() > 1) tsize = 1;
    else tsize = size();

    Array<Index> tarray, tperm(size(), 0);

    for(i = 0; i < tsize; i++) {
      Coordinate ttrans(at(i) - pivot[0]);
      if(ttrans.is_lattice_shift() && ((at(i) - ttrans) == pivot[0])) {
        (*this) -= ttrans;
      }
      else
        continue; //If not a lattice translation, move onto the next site



      //Following checks true/false condition
      //if condition is true, then sites are reordered so that first N sites of current cluster
      //match the N sites of 'pivot', in order

      if(find(pivot, tarray)) {
        Index t_ind, n_switch(0);
        for(j = 0; j < tperm.size(); j++) {
          t_ind = tarray.find(j);
          if(t_ind < tarray.size()) {
            tperm[t_ind] = j;

          }
          else {
            tperm[tarray.size() + n_switch] = j;
            n_switch++;
          }

        }
        permute(tperm);
        return true;
      }
    }
    return false;
  }

  //********************************************
  //This routine needs to be edited so that when sites are reordered,
  //the cluster basis functions are updated to reflect the change in order
  //num_maps should be from 0 to size of the cluster
  //returns true if successful maps == num_maps
  template <typename CoordType>
  bool GenericCluster<CoordType>::map_onto_subcluster(const GenericCluster<CoordType> &pivot, int num_maps) {
    Index i, j, tsize;

    if(!(get_home() == pivot.get_home())) {
      std::cerr << "WARNING in Cluster::map_onto_subcluser!!\n"
                << "You are trying to map a pivot onto a cluster with\n"
                << "a different lattice! \n";
    }


    if(PERIODICITY_MODE::IS_LOCAL() && size() > 1) tsize = 1;
    else tsize = size();

    Array<Index> tarray, tperm(size(), 0);
    int tmaps(0);
    for(i = 0; i < tsize; i++) {
      Coordinate ttrans(at(i) - pivot[0]);
      if(ttrans.is_lattice_shift() && ((at(i) - ttrans) == pivot[0])) {
        (*this) -= ttrans;
      }
      else continue; //If not a lattice translation, move onto the next site



      //Following checks true/false condition
      //if condition is true, then sites are reordered so that first N sites of current cluster
      //match the N sites of 'pivot', in order

      if(!find(pivot, tarray))continue;
      if(tmaps != num_maps) {
        tmaps++;
        continue;
      }
      Index t_ind, n_switch(0);
      for(j = 0; j < tperm.size(); j++) {
        t_ind = tarray.find(j);
        if(t_ind < tarray.size()) {
          tperm[t_ind] = j;

        }
        else {
          tperm[tarray.size() + n_switch] = j;
          n_switch++;
        }

      }
      permute(tperm);
      return true;

    }
    return false;
  }

  //********************************************************************

  template <typename CoordType>
  void GenericCluster<CoordType>::collect_basis_info(const Array<CoordType> &basis, const Coordinate &shift) {
    for(Index ns = 0; ns < size(); ns++) {
      at(ns).set_basis_ind(-1);
      for(Index nb = 0; nb < basis.size(); nb++) {
        if(at(ns).compare(basis[nb], shift)) {
          at(ns).set_basis_ind(nb);
          break;
        }
      }
      if(!valid_index(at(ns).basis_ind())) {
        std::cerr << "WARNING in GenericCluster::collect_basis_info" << std::endl;
        std::cerr << "Could not get basis info for " << at(ns) << std::endl;
        return;
      }
    }
  }

  //********************************************************************

  template <typename CoordType>
  void GenericCluster<CoordType>::collect_basis_info(const Array<CoordType> &basis) {
    for(Index ns = 0; ns < size(); ns++) {
      at(ns).set_basis_ind(-1);
      for(Index nb = 0; nb < basis.size(); nb++) {
        if(at(ns).compare(basis[nb])) {
          at(ns).set_basis_ind(nb);
          break;
        }
      }
      if(!valid_index(at(ns).basis_ind())) {
        std::cerr << "WARNING in GenericCluster::collect_basis_info" << std::endl;
        std::cerr << "Could not get basis info for " << at(ns) << std::endl;
        return;
      }
    }
  }

  //********************************************
  /**
   * Reads the cluster in the specified mode
   * and if read_tensors is set to true, it will
   * read the corresponding tensor basis and the
   * effective cluster interaction tensor.
   *
   * @param stream Input file stream
   * @param mode Cartesian or fractional mode
   * @param read_tensors Flag determining whether
   *    tensors should be read in or not
   */
  //********************************************
  template <typename CoordType>
  void GenericCluster<CoordType>::read(std::istream &stream, COORD_TYPE mode, bool read_tensors) {

    Index np;
    CoordType t_coord(*home);
    std::string tstring;
    char tchar[256];
    char ch;

    ch = stream.peek();

#ifdef DEBUG
    std::cout << "INSIDE CLUSTER READ \n";
    std::cout << "char(35) is " << char(35) << "\n";
#endif //DEBUG

    while((ch != char(35))) {
#ifdef DEBUG
      std::cout << "looking for # ch is " << ch << "\n";
#endif //DEBUG
      stream.ignore(1000, '\n');
      ch = stream.peek();
    }

    stream.getline(tchar, 1000, ':');
    stream >> np;

#ifdef DEBUG
    std::cout << "np is " << np << "\n";
#endif //DEBUG

    ch = stream.peek();
    while((ch != 'M') && (ch != 'm')) {
      stream.ignore(1000, '\n');
      ch = stream.peek();
    }

#ifdef DEBUG
    std::cout << "GOING TO TRY TO READ MAX MIN LENGTHS\n";
#endif //DEBUG

    stream.getline(tchar, 1000, ':'); //Changed from = to :
    stream >> max_length_val;
    stream.getline(tchar, 1000, ':');
    stream >> min_length_val;

#ifdef DEBUG
    std::cout << "max_length_val is " << max_length_val << "\n";
    std::cout << "min_length_val is " << min_length_val << "\n";
#endif //DEBUG

    // if it's the empty cluster, there isn't a tensor basis
    // to be read
    if(np == 0) {
      //std::cout << "Setting read_tensors to be false \n";
      //read_tensors = false;
    }

    COORD_MODE input_mode(mode);

    for(Index i = 0; i < np; i ++) {
      t_coord.read(stream);
      push_back(t_coord);
    }

    if(read_tensors == true) {
      //starting to read tensor basis
      tensor_basis.clear();
      ch = stream.peek();
#ifdef DEBUG
      std::cout << "ch in cluster::read in read_tensors == true is " << ch << "\n";
#endif //DEBUG
      while((ch != 'T') && (ch != 't')) {
        //stream.ignore(1000, '\n');
        stream.ignore(1000, '\n');
        ch = stream.peek();
        //std::cout << "ch = " << ch << "\n";
      }

      if(!tensor_basis.read(stream)) return;

      ch = stream.peek();
      while((ch != 'F') && (ch != 'f')) {
        stream.ignore(1000, '\n');
        ch = stream.peek();
      }

      stream.ignore(1000, '\n');

      ch = stream.peek();
      if((eci == 0) && (ch != 'N') && (ch != 'n')) {
        eci.redefine(Array<Index>(2, 3));
        eci.read(stream);
      }
    } //end of if read_tensors

    return;


  };

  //********************************************

  template <typename CoordType>
  void GenericCluster<CoordType>::push_back(const CoordType &new_coord) {

    Array<CoordType>::push_back(new_coord);

    if(!home)
      home = &back().home();

    else if(&(back().home()) != home)
      back().set_lattice(*home, CART);

    return;
  }

  //********************************************
  /// Merge adds the unique points of RHS onto the current cluster

  template <typename CoordType>
  void GenericCluster<CoordType>::merge(const GenericCluster<CoordType> &RHS) {
    for(Index i = 0; i < RHS.size(); i++) {
      if(!contains(RHS[i]))
        push_back(RHS[i]);
    }
    return;
  }


  //********************************************

  template <typename CoordType>
  GenericCluster<CoordType> &GenericCluster<CoordType>::operator+=(const Coordinate &RHS) {
    for(Index i = 0; i < size(); i++)
      at(i) += RHS;
    return *this;
  }

  //********************************************

  template <typename CoordType>
  GenericCluster<CoordType> &GenericCluster<CoordType>::operator-=(const Coordinate &RHS) {
    for(Index i = 0; i < size(); i++)
      at(i) -= RHS;
    return *this;
  }

  //********************************************

  template <typename CoordType>
  bool GenericCluster<CoordType>::operator==(const GenericCluster<CoordType> &RHS) const {
    if(size() != RHS.size()) return false;

    Array<bool> check_ind(size(), false);
    Index i, j;
    for(i = 0; i < size(); i++) {
      for(j = 0; j < RHS.size(); j++) {
        if((at(i) == RHS[j]) && (!check_ind[j])) {
          check_ind[j] = true;
          break;
        }
      }
      if(j == RHS.size())
        return false;
    }
    return true;
  }

  //********************************************

  template <typename CoordType>
  bool GenericCluster<CoordType>::is_equivalent(const GenericCluster<CoordType> &test_clust) const {
    if(size() != test_clust.size()) return false;
    if(PERIODICITY_MODE::IS_LOCAL()) return (*this) == test_clust;

    GenericCluster<CoordType> tclust(*this);
    for(Index i = 0; i < tclust.size(); i++) {
      tclust.within(i);
      if(tclust == test_clust) {
        return true;
      }
    }
    return false;
  }

  //********************************************

  template <typename CoordType>
  bool GenericCluster<CoordType>::is_equivalent(const GenericCluster<CoordType> &test_clust, Coordinate &trans) const {
    if(size() != test_clust.size()) return false;
    if(PERIODICITY_MODE::IS_LOCAL()) return (*this) == test_clust;

    GenericCluster<CoordType> tclust(*this);
    for(Index i = 0; i < tclust.size(); i++) {
      tclust.within(i, trans);
      if(tclust == test_clust) {
        return true;
      }
    }
    return false;
  }

  //********************************************
  //Checks for lattice translations that map the cluster onto test_clust
  //map_onto(...) is non-const, and if it returns true, the cluster
  //has been translated onto test_clust, although the points may be ordered differently.
  //If it returns false, the cluster is unchanged.

  template <typename CoordType>
  bool GenericCluster<CoordType>::map_onto(const GenericCluster<CoordType> &test_clust) {
    if(size() != test_clust.size()) return false;
    if(size() == 0) return true;

    //If non-periodic, translations need not be considered
    if(PERIODICITY_MODE::IS_LOCAL()) return (*this) == test_clust;

    //tshift keeps track of translations
    Coordinate tshift(*home), trans(*home);
    for(Index i = 0; i < size(); i++) {
      tshift = test_clust[0] - at(i);
      if(tshift.is_lattice_shift()) {
        (*this) += tshift;
        trans += tshift;
        if((*this) == test_clust)
          return true;
      }
    }
    (*this) -= trans;
    return false;
  }

  //********************************************

  template <typename CoordType>
  bool GenericCluster<CoordType>::map_onto(const GenericCluster<CoordType> &test_clust, Coordinate &trans) {
    if(size() != test_clust.size()) return false;
    trans.frac() = Eigen::Vector3d::Zero();
    if(size() == 0) return true;

    if(PERIODICITY_MODE::IS_LOCAL()) {
      return (*this) == test_clust;
    }

    Coordinate tshift(*home);
    for(Index i = 0; i < size(); i++) {
      tshift = test_clust[0] - at(i);
      if(tshift.is_lattice_shift()) {
        (*this) += tshift;
        trans += tshift;
        if((*this) == test_clust)
          return true;
      }
    }
    (*this) -= trans;
    trans.frac() = Eigen::Vector3d::Zero();
    return false;
  }

  //********************************************

  template <typename CoordType>
  std::complex<double> GenericCluster<CoordType>::get_phase(const Coordinate &k, int i, int j) {
    double targ = k.cart().dot(s2s_vec[i][j].cart());
    return std::complex<double>(cos(targ), sin(targ));
  }

  //********************************************
  //fill s2s_vec with Coordinates that describe vectors
  //pointing from site i to site j (recorded in s2s_vec[i][j])

  template <typename CoordType>
  void GenericCluster<CoordType>::get_s2s_vec() {
    Index i, j;
    s2s_vec.clear();
    s2s_vec.resize(size());
    s2s_norm_vec.resize(size());
    for(i = 0; i < size(); i++) {
      s2s_vec[i].reserve(size());
      s2s_norm_vec[i].reserve(size());
      for(j = 0; j < size(); j++) {
        s2s_vec[i].push_back(at(j) - at(i));
        s2s_norm_vec[i].push_back(s2s_vec[i].back());
        s2s_norm_vec[i].back().cart() /= s2s_norm_vec[i].back().const_cart().norm();
      }
    }
    return;
  }


  //********************************************

  template <typename CoordType>
  void GenericCluster<CoordType>::calc_properties() {
    double tlength;

    //Point clusters don't have a max_length - (is this necessary?)
    if(size() <= 1)
      max_length_val = min_length_val = 0;

    else if(size() > 1) {
      //Establish max and min as distance from first set
      max_length_val = min_length_val = at(0).dist(at(1));

      for(Index i = 0; i < size(); i++) {
        for(Index j = i + 1; j < size(); j++) {
          tlength = at(i).dist(at(j));
          if(tlength < min_length_val)
            min_length_val = tlength;
          if(tlength > max_length_val)
            max_length_val = tlength;
        }
      }
    }
    return;
  }



  //********************************************
  /**
   * \brief   calculate min/max lengths relative the phenom_clust & this cluster
   *
   *   - sets max_length based on lengths between sites in the cluster or lengths between sites in the cluster and the phenom_clust
   *   - only set min_length based on lengths between sites in the cluster
   */
  //************************************************************

  template <typename CoordType>
  void GenericCluster<CoordType>::calc_properties(GenericCluster<CoordType> phenom_clust) {
    double dist;

    // calculate min/max lengths relative the phenom_clust & this cluster
    max_length_val = 0;
    min_length_val = 1e20;

    // first check distances to phenom_clust
    for(Index i = 0; i < size(); i++) {
      for(Index k = 0; k < phenom_clust.size(); k++) {
        dist = phenom_clust[k].dist(at(i));
        if(dist > max_length_val)
          max_length_val = dist;
        if(dist < min_length_val)		// only set min_length for length between sites in the cluster
          min_length_val = dist;

      }
    }

    // then to sites in this cluster
    for(Index i = 0; i < size(); i++) {
      for(Index k = i + 1; k < size(); k++) {
        dist = at(k).dist(at(i));
        if(dist > max_length_val)
          max_length_val = dist;
        if(dist < min_length_val)
          min_length_val = dist;

      }
    }

    //std::cout << "  max_length: " << max_length_val << "  min_length: " << min_length_val << endl;

    return;
  }
  //********************************************

  template <typename CoordType>
  void GenericCluster<CoordType>::update() {
    for(Index i = 0; i < size(); i++)
      at(i).update();
    calc_properties();
    return;
  }

  //********************************************

  template <typename CoordType>
  void GenericCluster<CoordType>::print(std::ostream &stream, char delim, COORD_TYPE mode) const {
    if(mode == COORD_DEFAULT)
      mode = COORD_MODE::CHECK();
    COORD_MODE C(mode);

    stream << "#Points: " << size() << std::endl
           << "MaxLength: " << max_length_val << "  MinLength: " << min_length_val << std::endl;
    for(Index np = 0; np < size(); np++) {
      stream.setf(std::ios::showpoint, std::ios_base::fixed);
      stream.precision(5);
      stream.width(9);
      at(np).print(stream);//Changed by Ivy from at(np).print(stream,mode) -- the "mode" should actually be the SD_flag input 11/04/12
      if(delim)
        stream << delim;
    };
  };

  //********************************************

  template <typename CoordType>
  void GenericCluster<CoordType>::print_shifted(std::ostream &stream, const Coordinate &shift, char delim, COORD_TYPE mode) const {
    if(mode == COORD_DEFAULT)
      mode = COORD_MODE::CHECK();
    COORD_MODE C(mode);

    stream << "#Points: " << size() << std::endl
           << "MaxLength: " << max_length_val << "  MinLength: " << min_length_val << std::endl;
    for(Index np = 0; np < size(); np++) {
      stream.setf(std::ios::showpoint, std::ios_base::fixed);
      stream.precision(5);
      stream.width(9);
      (at(np) + shift).print(stream); //Changed by Ivy from at(np).print(stream,mode) -- the "mode" should actually be the SD_flag input 11/04/12
      if(delim)
        stream << delim;
    };

  };

  //********************************************

  template <typename CoordType>
  void GenericCluster<CoordType>::print_sites(std::ostream &stream, int space, char delim, COORD_TYPE mode) const {
    if(mode == COORD_DEFAULT)
      mode = COORD_MODE::CHECK();
    COORD_MODE C(mode);
    for(Index np = 0; np < size(); np++) {
      for(int i = 0; i < space; i++) {
        stream << ' ';
      }
      stream.setf(std::ios::showpoint, std::ios_base::fixed);
      stream.precision(5);
      stream.width(9);
      at(np).print(stream);//Changed by Ivy from at(np).print(stream,mode) -- the "mode" should actually be the SD_flag input 11/04/12
      if(delim)
        stream << delim;
    };
  };

  //********************************************

  template <typename CoordType>
  void GenericCluster<CoordType>::print_basis_info(std::ostream &stream, int space, char delim, COORD_TYPE mode) const {
    if(mode == COORD_DEFAULT)
      mode = COORD_MODE::CHECK();
    COORD_MODE C(mode);
    for(Index np = 0; np < size(); np++) {
      for(int i = 0; i < space; i++) {
        stream << ' ';
      }
      stream.setf(std::ios::showpoint, std::ios_base::fixed);
      stream.precision(5);
      stream.width(9);
      at(np).print(stream);
      stream << "  " << at(np).basis_ind() << " ";
      if(delim)
        stream << delim;
    };
  };

  //********************************************

  template <typename CoordType>
  void GenericCluster<CoordType>::print_decorated_sites(std::ostream &stream, int space, char delim, COORD_TYPE mode) const {
    if(mode == COORD_DEFAULT)
      mode = COORD_MODE::CHECK();
    COORD_MODE C(mode);
    for(Index np = 0; np < size(); np++) {
      for(int i = 0; i < space; i++) {
        stream << ' ';
      }
      stream.setf(std::ios::showpoint, std::ios_base::fixed);
      stream.precision(5);
      stream.width(9);
      at(np).print_occ(stream);//Changed by Ivy from at(np).print(stream,mode) -- the "mode" should actually be the SD_flag input 11/04/12
      if(delim)
        stream << delim;
    };
  };


  // //****************************************************
  // template <typename CoordType>
  // GenericOrbitree< GenericCluster<CoordType> > GenericCluster<CoordType>::enumerate_subclusters(const Structure &prim, bool verbose) const {
  //   if(prim.factor_group().size() == 0) {
  //     std::cerr << "WARNING: In Orbitree::generate_orbitree, prim's factor_group is empty. It  must at least have one element (identity).\n";
  //     assert(0);
  //   }

  //   //Ensure that the lattices are the same:
  //   if ( !(*home==prim.lattice) ){
  //     std::cerr<<"WARNING in GenericCluster<CoordType>::enumerate_subclusters, the lattice in prim and the lattice"
  //              <<" that was used to construct this cluster are not the same"<<std::endl;
  //     assert(0);
  //   }

  //   Index i,j;
  //   std::string clean(80, ' ');
  //   Array<int> master_choose(size(),0);
  //   GenericOrbitree< GenericCluster<CoordType> > sub_orbitree(prim.lattice);
  //   //Setup Orbitree here?
  //   sub_orbitree.lattice = prim.lattice;
  //   sub_orbitree.resize( size()+1, GenericOrbitBranch< GenericCluster<CoordType> >( GenericCluster<CoordType>(prim.lattice) ) );
  //   // Add in the empty cluster
  //   at(0).push_back(GenericOrbit< GenericCluster<CoordType>  >(GenericCluster<CoordType>(prim.lattice)));
  //   at(0).back().get_equivalent(prim.factor_group());
  //   at(0).back().get_cluster_symmetry();

  //   for(i=1;i<=size();i++){
  //     Array<int> choose = master_choose;
  //     for(j=0;j<i;j++)
  //       choose[choose.size()-i-1] = 1;
  //     GenericCluster<CoordType> test_clust(prim.lattice);
  //     do{
  //       test_clust.clear();
  //       for(j=0;j<choose.size();j++){
  //         if(choose[j]==1)
  //           test_clust.push_back(at(j));
  //       }
  //       test_clust.within();
  //       test_clust.calc_properties();

  //       if(!sub_orbitree.contains(test_clust)){
  //         sub_orbitree[i].push_back( GenericOrbit< GenericCluster<CoordType> >(test_clust)  );
  //         sub_orbitree[i].back().get_equivalent(prim.factor_group());
  //         sub_orbitree[i].back().get_cluster_symmetry();
  //       }
  //     }while(choose.next_permute());
  //   }
  // }


  //********************************************
  /**
   * Allows for printing using the << operator.
   */
  //********************************************
  template <typename CoordType>
  std::ostream &operator<< (std::ostream &stream, const GenericCluster<CoordType> &cluster) {
    cluster.print(stream, '\n'); //Ivy added newline delimiter 07/01/13
    return stream;
  };


  //****************************************************
  template <typename CoordType>
  GenericCluster<CoordType> operator*(const SymOp &LHS, const GenericCluster<CoordType> &RHS) {
    return GenericCluster<CoordType>(RHS).apply_sym(LHS);
  }

  //****************************************************
  template <typename CoordType>
  GenericCluster<CoordType> operator+(const GenericCluster<CoordType> &LHS, const Coordinate &RHS) {
    return GenericCluster<CoordType>(LHS) += RHS;
  }


  //****************************************************
  template <typename CoordType>
  GenericCluster<CoordType> operator-(const GenericCluster<CoordType> &LHS, const Coordinate &RHS) {
    return GenericCluster<CoordType>(LHS) -= RHS;
  }


};

