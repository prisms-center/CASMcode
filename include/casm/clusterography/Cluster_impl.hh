#include "casm/crystallography/Lattice.hh"
#include "casm/symmetry/SymPermutation.hh"

namespace CASM {


  //****************************************************
  template <typename CoordType>
  GenericCluster<CoordType>::GenericCluster(const Lattice &init_home) :
    m_lat_ptr(&init_home),
    m_min_length(0.0),
    m_max_length(0.0),
    m_clust_group(LOCAL) {

    //nothing else to do for now
  }

  //****************************************************

  template <typename CoordType>
  void GenericCluster<CoordType>::set_lattice(const Lattice &new_home, COORD_TYPE mode) {
    if(m_lat_ptr == &new_home)
      return;
    m_lat_ptr = &new_home;
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

    m_clust_group.apply_sym(op);
    return *this;
  }

  //********************************************

  template <typename CoordType>
  GenericCluster<CoordType> &GenericCluster<CoordType>::apply_sym_no_trans(const SymOp &op) {
    for(Index i = 0; i < size(); i++)
      at(i).apply_sym_no_trans(op);

    m_clust_group.apply_sym(op);
    return *this;
  }


  //*******************************************************************************************
  // Assuming that we have the factor group (or some similar group),
  // what are the operations that leave the cluster unchanged
  // Find them, alter the translations, and store them in clust_group

  template <typename _CoordType>
  void GenericCluster<_CoordType>::generate_clust_group(const SymGroup &super_group, std::vector<Permutation> *perm_array_ptr, double tol) {
    Coordinate trans(home());
    _clust_group().clear();
    if(perm_array_ptr)
      perm_array_ptr->clear();
    Array<Index> iperm;
    for(Index ng = 0; ng < super_group.size(); ng++) {
      GenericCluster<_CoordType> tclust(*this);
      tclust.apply_sym(super_group[ng]);
      if(tclust.map_onto(*this, trans)) {
        _clust_group().push_back(SymOp::translation(trans.cart())*super_group[ng]);
        find(tclust, iperm, tol);

        if(perm_array_ptr)
          perm_array_ptr->push_back(Permutation(iperm));
      }
    }

    return;
  }

  //*******************************************************************************************

  template <typename _CoordType>
  std::vector<Permutation> GenericCluster<_CoordType>::clust_group_permutations(double tol) const {
    std::vector<Permutation> tperms;
    tperms.reserve(clust_group().size());
    Array<Index> iperm(size());
    for(Index i = 0; i < clust_group().size(); i++) {
      for(Index j = 0; j < size(); j++)
        iperm[j] = find(clust_group()[i] * at(j), tol);
      tperms.push_back(Permutation(iperm));
    }
    return tperms;
  }

  //********************************************

  template <typename CoordType>
  void GenericCluster<CoordType>::within(Index pivot_ind) {

    if(!size() || PERIODICITY_MODE::IS_LOCAL()) return;

    if(pivot_ind >= size()) {
      std::cerr << "WARNING: Attempted to map cluster within unit cell with invalid pivot." << std::endl;
      return;
    }

    Coordinate tcoord(*m_lat_ptr);
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

    trans.set_lattice(*m_lat_ptr, CART);
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
    Coordinate tcoord(*m_lat_ptr);
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

  template <typename CoordType>
  Coordinate GenericCluster<CoordType>::geometric_center() const {
    if(!size())
      return Coordinate(*m_lat_ptr);

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
  bool GenericCluster<CoordType>::contains_periodic(const CoordType &test_site, double tol) const {
    for(Index i = 0; i < size(); i++) {
      if(at(i).compare(test_site, tol))
        return true;
    }
    return false;
  }
  //********************************************
  //returns true if test_cluster is a subcluster of cluster (*this)
  //Array 'index' holds indices of (*this) that yield test_cluster
  template <typename CoordType>
  Index GenericCluster<CoordType>::find(const CoordType &test_elem, double tol) const {
    Index j;

    for(j = 0; j < size(); j++) {
      if(test_elem.almost_equal(at(j), tol))
        return j;
    }

    return j;
  }
  //********************************************
  //returns true if test_cluster is a subcluster of cluster (*this)
  //Array 'index' holds indices of (*this) that yield test_cluster
  template <typename CoordType>
  bool GenericCluster<CoordType>::find(const GenericCluster<CoordType> &test_cluster, Array<Index> &index, double tol) const {
    Index i, j;
    index.clear();
    index.reserve(test_cluster.size());
    for(i = 0; i < test_cluster.size(); i++) {
      for(j = 0; j < size(); j++) {
        if(test_cluster[i].almost_equal(at(j), tol) && !index.contains(j))
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
  bool GenericCluster<CoordType>::map_onto_subcluster(const GenericCluster<CoordType> &pivot, double tol) {
    Index i, j, tsize;

    if(!(home() == pivot.home())) {
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

      if(find(pivot, tarray, tol)) {
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
  bool GenericCluster<CoordType>::map_onto_subcluster(const GenericCluster<CoordType> &pivot, int num_maps, double tol) {
    Index i, j, tsize;

    if(!(home() == pivot.home())) {
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

      if(!find(pivot, tarray, tol))continue;
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
   *
   * @param stream Input file stream
   * @param mode Cartesian or fractional mode
   */
  //********************************************
  template <typename CoordType>
  void GenericCluster<CoordType>::read(std::istream &stream, COORD_TYPE mode) {

    Index np;
    CoordType t_coord(*m_lat_ptr);
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
    stream >> m_max_length;
    stream.getline(tchar, 1000, ':');
    stream >> m_min_length;

#ifdef DEBUG
    std::cout << "m_max_length is " << m_max_length << "\n";
    std::cout << "m_min_length is " << m_min_length << "\n";
#endif //DEBUG


    COORD_MODE input_mode(mode);

    for(Index i = 0; i < np; i ++) {
      t_coord.read(stream);
      push_back(t_coord);
    }

    return;


  }

  //********************************************

  template <typename CoordType>
  void GenericCluster<CoordType>::push_back(const CoordType &new_coord) {

    Array<CoordType>::push_back(new_coord);

    if(!m_lat_ptr)
      m_lat_ptr = &back().home();

    else if(&(back().home()) != m_lat_ptr)
      back().set_lattice(*m_lat_ptr, CART);

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
  bool GenericCluster<CoordType>::map_onto(const GenericCluster<CoordType> &test_clust, double tol) {
    if(size() != test_clust.size()) return false;
    if(size() == 0) return true;

    //If non-periodic, translations need not be considered
    if(PERIODICITY_MODE::IS_LOCAL()) return (*this) == test_clust;

    //tshift keeps track of translations
    Coordinate tshift(*m_lat_ptr), trans(*m_lat_ptr);
    for(Index i = 0; i < size(); i++) {
      tshift = test_clust[0] - at(i);
      if(tshift.is_lattice_shift(tol)) {
        (*this) += tshift;
        trans += tshift;
        if(almost_equal((*this), test_clust, tol))
          return true;
      }
    }
    (*this) -= trans;
    return false;
  }

  //********************************************

  template <typename CoordType>
  bool GenericCluster<CoordType>::map_onto(const GenericCluster<CoordType> &test_clust, Coordinate &trans, double tol) {
    if(size() != test_clust.size()) return false;
    trans.frac() = Eigen::Vector3d::Zero();
    if(size() == 0) return true;

    if(PERIODICITY_MODE::IS_LOCAL()) {
      return (*this) == test_clust;
    }

    Coordinate tshift(*m_lat_ptr);
    for(Index i = 0; i < size(); i++) {
      tshift = test_clust[0] - at(i);
      if(tshift.is_lattice_shift(tol)) {
        (*this) += tshift;
        trans += tshift;
        if(almost_equal((*this), test_clust, tol))
          return true;
      }
    }
    (*this) -= trans;
    trans.frac() = Eigen::Vector3d::Zero();
    return false;
  }

  //********************************************

  template <typename CoordType>
  void GenericCluster<CoordType>::calc_properties() {
    double tlength;

    //Point clusters don't have a max_length - (is this necessary?)
    if(size() <= 1)
      m_max_length = m_min_length = 0;

    else if(size() > 1) {
      //Establish max and min as distance from first set
      m_max_length = m_min_length = at(0).dist(at(1));

      for(Index i = 0; i < size(); i++) {
        for(Index j = i + 1; j < size(); j++) {
          tlength = at(i).dist(at(j));
          if(tlength < m_min_length)
            m_min_length = tlength;
          if(tlength > m_max_length)
            m_max_length = tlength;
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
    m_max_length = 0;
    m_min_length = 1e20;

    // first check distances to phenom_clust
    for(Index i = 0; i < size(); i++) {
      for(Index k = 0; k < phenom_clust.size(); k++) {
        dist = phenom_clust[k].dist(at(i));
        if(dist > m_max_length)
          m_max_length = dist;
        if(dist < m_min_length)		// only set min_length for length between sites in the cluster
          m_min_length = dist;

      }
    }

    // then to sites in this cluster
    for(Index i = 0; i < size(); i++) {
      for(Index k = i + 1; k < size(); k++) {
        dist = at(k).dist(at(i));
        if(dist > m_max_length)
          m_max_length = dist;
        if(dist < m_min_length)
          m_min_length = dist;

      }
    }

    //std::cout << "  max_length: " << m_max_length << "  min_length: " << m_min_length << endl;

    return;
  }

  //********************************************

  template <typename CoordType>
  void GenericCluster<CoordType>::print(std::ostream &stream, char delim, COORD_TYPE mode) const {
    if(mode == COORD_DEFAULT)
      mode = COORD_MODE::CHECK();
    COORD_MODE C(mode);

    stream << "#Points: " << size() << std::endl
           << "MaxLength: " << m_max_length << "  MinLength: " << m_min_length << std::endl;
    for(Index np = 0; np < size(); np++) {
      stream.setf(std::ios::showpoint, std::ios_base::fixed);
      stream.precision(5);
      stream.width(9);
      at(np).print(stream);//Changed by Ivy from at(np).print(stream,mode) -- the "mode" should actually be the SD_flag input 11/04/12
      if(delim)
        stream << delim;
    }
  }

  //********************************************

  template <typename CoordType>
  void GenericCluster<CoordType>::print_shifted(std::ostream &stream, const Coordinate &shift, char delim, COORD_TYPE mode) const {
    if(mode == COORD_DEFAULT)
      mode = COORD_MODE::CHECK();
    COORD_MODE C(mode);

    stream << "#Points: " << size() << std::endl
           << "MaxLength: " << m_max_length << "  MinLength: " << m_min_length << std::endl;
    for(Index np = 0; np < size(); np++) {
      stream.setf(std::ios::showpoint, std::ios_base::fixed);
      stream.precision(5);
      stream.width(9);
      (at(np) + shift).print(stream); //Changed by Ivy from at(np).print(stream,mode) -- the "mode" should actually be the SD_flag input 11/04/12
      if(delim)
        stream << delim;
    }

  }

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
    }
  }

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
    }
  }

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
    }
  }

  //********************************************
  /**
   * Allows for printing using the << operator.
   */
  //********************************************
  template <typename CoordType>
  std::ostream &operator<< (std::ostream &stream, const GenericCluster<CoordType> &cluster) {
    cluster.print(stream, '\n'); //Ivy added newline delimiter 07/01/13
    return stream;
  }


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

  //****************************************************
  template <typename CoordType>
  bool almost_equal(const GenericCluster<CoordType> &LHS, const GenericCluster<CoordType> &RHS, double tol) {
    if(LHS.size() != RHS.size()) return false;

    Array<bool> check_ind(LHS.size(), false);
    Index i, j;
    for(i = 0; i < LHS.size(); i++) {
      for(j = 0; j < RHS.size(); j++) {
        if((LHS[i].almost_equal(RHS[j], tol)) && (!check_ind[j])) {
          check_ind[j] = true;
          break;
        }
      }
      if(j == RHS.size())
        return false;
    }
    return true;


  }



}
