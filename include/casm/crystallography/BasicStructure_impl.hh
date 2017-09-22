#include <sstream>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/Coordinate.hh"
#include "casm/crystallography/UnitCellCoord.hh"
#include "casm/crystallography/PrimGrid.hh"
#include "casm/symmetry/SymPermutation.hh"
#include "casm/symmetry/SymBasisPermute.hh"
#include "casm/symmetry/SymGroupRep.hh"
#include "casm/crystallography/Niggli.hh"
#include "casm/casm_io/json_io/container.hh"

namespace CASM {
  template<typename CoordType>
  BasicStructure<CoordType>::BasicStructure(const fs::path &filepath) : m_lattice() {
    if(!fs::exists(filepath)) {
      std::cout << "Error in BasicStructure<CoordType>::BasicStructure<CoordType>(const fs::path &filepath)." << std::endl;
      std::cout << "  File does not exist at: " << filepath << std::endl;
      exit(1);
    }
    fs::ifstream infile(filepath);

    read(infile);
  }

  //***********************************************************

  template<typename CoordType>
  BasicStructure<CoordType>::BasicStructure(const BasicStructure &RHS) :
    m_lattice(RHS.lattice()), title(RHS.title), basis(RHS.basis) {
    for(Index i = 0; i < basis.size(); i++) {
      basis[i].set_lattice(lattice(), CART);
    }
  }

  //***********************************************************

  template<typename CoordType>
  BasicStructure<CoordType>::~BasicStructure() {}

  //***********************************************************

  template<typename CoordType>
  BasicStructure<CoordType> &BasicStructure<CoordType>::operator=(const BasicStructure<CoordType> &RHS) {
    m_lattice = RHS.lattice();
    title = RHS.title;
    basis = RHS.basis;
    for(Index i = 0; i < basis.size(); i++) {
      basis[i].set_lattice(lattice(), CART);
    }

    return *this;
  }


  //***********************************************************

  template<typename CoordType>
  void BasicStructure<CoordType>::copy_attributes_from(const BasicStructure<CoordType> &RHS) {

  }

  //***********************************************************
  /*
  template<typename CoordType>
  BasicStructure<CoordType> &BasicStructure<CoordType>::apply_sym(const SymOp &op) {
    for(Index i = 0; i < basis.size(); i++) {
      basis[i].apply_sym(op);
    }
    return *this;
  }
  */
  //***********************************************************

  template<typename CoordType>
  void BasicStructure<CoordType>::reset() {


    for(Index nb = 0; nb < basis.size(); nb++) {
      basis[nb].set_basis_ind(nb);
    }
    within();
    //set_site_internals();

  }

  //***********************************************************

  template<typename CoordType>
  void BasicStructure<CoordType>::update() {
    set_site_internals();
  }

  //***********************************************************

  template<typename CoordType>
  void BasicStructure<CoordType>::within() {
    for(Index i = 0; i < basis.size(); i++) {
      basis[i].within();
    }
    return;
  }

  //***********************************************************

  template<typename CoordType>
  void BasicStructure<CoordType>::generate_factor_group_slow(SymGroup &factor_group) const {
    Array<CoordType> trans_basis;
    Index pg, b0, b1, b2;
    Coordinate t_tau(lattice());
    Index num_suc_maps;

    SymGroup point_group;
    lattice().generate_point_group(point_group);

    if(factor_group.size() != 0) {
      std::cerr << "WARNING in BasicStructure<CoordType>::generate_factor_group_slow" << std::endl;
      std::cerr << "The factor group passed isn't empty and it's about to be rewritten!" << std::endl;
      factor_group.clear();
    }
    factor_group.set_lattice(lattice());
    //Loop over all point group ops of the lattice
    for(pg = 0; pg < point_group.size(); pg++) {
      trans_basis.clear();
      //First, generate the symmetrically transformed basis sites
      //Loop over all sites in basis
      for(b0 = 0; b0 < basis.size(); b0++) {
        trans_basis.push_back(point_group[pg]*basis[b0]);
      }

      //Using the symmetrically transformed basis, find all possible translations
      //that MIGHT map the symmetrically transformed basis onto the original basis
      for(b0 = 0; b0 < trans_basis.size(); b0++) {

        if(!basis[0].compare_type(trans_basis[b0]))
          continue;

        t_tau = basis[0] - trans_basis[b0];

        t_tau.within();
        num_suc_maps = 0; //Keeps track of number of old->new basis site mappings that are found

        double tdist = 0.0;
        double max_error = 0.0;
        std::vector<Index> mappings(basis.size(), basis.size());
        for(b1 = 0; b1 < basis.size(); b1++) { //Loop over original basis sites
          for(b2 = 0; b2 < trans_basis.size(); b2++) { //Loop over symmetrically transformed basis sites

            //see if translation successfully maps the two sites uniquely
            if(basis[b1].compare(trans_basis[b2], t_tau)) {// maybe nuke this with Hungarian
              tdist = basis[b1].min_dist(Coordinate(trans_basis[b2]) + t_tau);
              if(tdist > max_error) {
                max_error = tdist;
              }
              mappings[b1] = b2;
              num_suc_maps++;
              break;
            }
          }

          //break out of outer loop if inner loop finds no successful map
          if(b2 == trans_basis.size()) {
            break;
          }
        }
        std::set<Index> unique_mappings;
        for(auto &e : mappings)
          unique_mappings.insert(e);
        if(num_suc_maps == basis.size() && unique_mappings.size() == mappings.size()) {
          //If all atoms in the basis are mapped successfully, try to add the corresponding
          //symmetry operation to the factor_group
          Coordinate center_of_mass(lattice());
          for(Index b = 0; b < basis.size(); b++) {
            //for each basis site loop through all trans_basis to find the closest one
            double smallest = 1000000;
            Coordinate tshift(lattice());
            double dist = trans_basis[mappings[b]].min_dist(basis[b] - t_tau, tshift);
            //in tshift is stored trans_basis - basis
            tshift.cart() *= (1.0 / basis.size());
            center_of_mass += tshift;
          }

          /*
          point_group[pg] operates on all of basis and save = trans_basis
          //t_shift is the vector from basis to op_basis magnitude of tshift=min_dist
          average all t_shifts and add/subtract to t_tau

          if t_shift is b-> ob then subtract
          if t_shift is ob -> b then add
          matrix * basis + tau = operbasis
          */
          t_tau -= center_of_mass;/**/

          SymOp tSym(SymOp::translation(t_tau.cart())*point_group[pg]);
          tSym.set_map_error(max_error);

          if(!factor_group.contains(tSym)) {
            factor_group.push_back(tSym);
          }
        }
      }
    } //End loop over point_group operations
    factor_group.enforce_group(lattice().tol());
    factor_group.sort();
    factor_group.max_error();

    return;
  }

  //************************************************************

  template<typename CoordType>
  void BasicStructure<CoordType>::generate_factor_group(SymGroup &factor_group) const {
    BasicStructure<CoordType> tprim;
    factor_group.clear();
    factor_group.set_lattice(lattice());
    // CASE 1: Structure is primitive
    if(is_primitive(tprim)) {
      generate_factor_group_slow(factor_group);
      return;
    }


    // CASE 2: Structure is not primitive

    PrimGrid prim_grid(tprim.lattice(), lattice());
    SymGroup prim_fg;
    tprim.generate_factor_group_slow(prim_fg);

    SymGroup point_group;
    lattice().generate_point_group(point_group);
    point_group.enforce_group(lattice().tol());

    for(Index i = 0; i < prim_fg.size(); i++) {
      if(point_group.find_no_trans(prim_fg[i]) == point_group.size()) {
        continue;
      }
      else {
        for(Index j = 0; j < prim_grid.size(); j++) {
          factor_group.push_back(SymOp::translation(prim_grid.coord(j, SCEL).cart())*prim_fg[i]);
          // set lattice, in case SymOp::operator* ever changes
        }
      }
    }

    return;
  }

  //************************************************************

  template<typename CoordType>
  void BasicStructure<CoordType>::fg_converge(double small_tol, double large_tol, double increment) {
    SymGroup factor_group;
    fg_converge(factor_group, small_tol, large_tol, increment);
  }

  //************************************************************

  template<typename CoordType>
  void BasicStructure<CoordType>::fg_converge(SymGroup &factor_group, double small_tol, double large_tol, double increment) {

    Array<double> tols;
    Array<bool> is_group;
    Array<int> num_ops, num_enforced_ops;
    Array<std::string> name;

    double orig_tol = lattice().tol();
    for(double i = small_tol; i < large_tol; i += increment) {
      tols.push_back(i);
      m_lattice.set_tol(i);

      factor_group.clear();
      generate_factor_group(factor_group);
      factor_group.get_multi_table();
      num_ops.push_back(factor_group.size());
      is_group.push_back(factor_group.is_group(i));
      factor_group.enforce_group(i);
      num_enforced_ops.push_back(factor_group.size());
      factor_group.character_table();
      name.push_back(factor_group.get_name());
    }
    m_lattice.set_tol(orig_tol);

    for(Index i = 0; i < tols.size(); i++) {
      std::cout << tols[i] << "\t" << num_ops[i] << "\t" << is_group[i] << "\t" << num_enforced_ops[i] << "\t name: " << name[i] << "\n";
    }

    return;
  }

  //***********************************************************
  /**
   * It is NOT wise to use this function unless you have already
   * initialized a superstructure with lattice vectors.
   *
   * It is more wise to use the two methods that call this method:
   * Either the overloaded * operator which does:
   *  SCEL_Lattice * Prim_Structrue = New_Superstructure
   *       --- or ---
   *  New_Superstructure=Prim_BasicStructure<CoordType>.create_superstruc(SCEL_Lattice);
   *
   *  Both of these will return NEW superstructures.
   */
  //***********************************************************

  template<typename CoordType>
  void BasicStructure<CoordType>::fill_supercell(const BasicStructure<CoordType> &prim) {
    Index i, j;

    copy_attributes_from(prim);

    PrimGrid prim_grid(prim.lattice(), lattice());

    basis.clear();

    //loop over basis sites of prim
    for(j = 0; j < prim.basis.size(); j++) {

      //loop over prim_grid points
      for(i = 0; i < prim_grid.size(); i++) {

        //push back translated basis site of prim onto superstructure basis
        basis.push_back(prim.basis[j] + prim_grid.coord(i, PRIM));

        //reset lattice for most recent superstructure CoordType
        //set_lattice() converts fractional coordinates to be compatible with new lattice
        basis.back().set_lattice(lattice(), CART);

        basis.back().within();
      }
    }

    set_site_internals();

    return;
  }

  //***********************************************************
  /**
   * Operates on the primitive structure and takes as an argument
   * the supercell lattice.  It then returns a new superstructure.
   *
   * This is similar to the Lattice*Primitive routine which returns a
   * new superstructure.  Unlike the fill_supercell routine which takes
   * the primitive structure, this WILL fill the sites.
   */
  //***********************************************************

  template<typename CoordType>
  BasicStructure<CoordType> BasicStructure<CoordType>::create_superstruc(const Lattice &scel_lat) const {
    BasicStructure<CoordType> tsuper(scel_lat);
    tsuper.fill_supercell(*this);
    return tsuper;
  }


  //***********************************************************
  /**
   * Determines if structure is primitive description of the crystal
   */
  //***********************************************************

  template<typename CoordType>
  bool BasicStructure<CoordType>::is_primitive() const {
    Coordinate tshift(lattice());//, bshift(lattice);
    Index b1, b2, b3, num_suc_maps;

    for(b1 = 1; b1 < basis.size(); b1++) {
      if(!basis[0].compare_type(basis[b1])) {
        continue;
      }

      tshift = basis[0] - basis[b1];
      num_suc_maps = 0;
      for(b2 = 0; b2 < basis.size(); b2++) {
        for(b3 = 0; b3 < basis.size(); b3++) {
          //if(basis[b3].compare_type(basis[b2], bshift) && tshift.min_dist(bshift) < lattice().tol()) {
          if(basis[b3].compare(basis[b2], tshift)) {
            num_suc_maps++;
            break;
          }
        }
        if(b3 == basis.size()) {
          break;
        }
      }

      if(num_suc_maps == basis.size()) {
        return false;
      }

    }

    return true;
  }


  //***********************************************************
  /**
   * Determines if structure is primitive description of the crystal
   * If not, finds primitive cell and copies to new_prim
   */
  //***********************************************************

  template<typename CoordType>
  bool BasicStructure<CoordType>::is_primitive(BasicStructure<CoordType> &new_prim) const {

    Coordinate tshift(lattice());//, bshift(lattice);
    Eigen::Vector3d prim_vec0(lattice()[0]), prim_vec1(lattice()[1]), prim_vec2(lattice()[2]);
    Array<Eigen::Vector3d > shift;
    Index b1, b2, b3, sh, sh1, sh2;
    Index num_suc_maps;
    double tvol, min_vol;
    bool prim_flag = true;
    double prim_vol_tol = std::abs(0.5 * lattice().vol() / double(basis.size())); //sets a hard lower bound for the minimum value of the volume of the primitive cell


    for(b1 = 1; b1 < basis.size(); b1++) {

      tshift = basis[0] - basis[b1];

      if(almost_zero(tshift.min_dist(Coordinate::origin(lattice()))))
        continue;

      num_suc_maps = 0;
      for(b2 = 0; b2 < basis.size(); b2++) {
        for(b3 = 0; b3 < basis.size(); b3++) {
          if(basis[b3].compare(basis[b2], tshift)) {
            num_suc_maps++;
            break;
          }
        }
        if(b3 == basis.size()) {
          break;
        }
      }

      if(num_suc_maps == basis.size()) {
        prim_flag = false;
        shift.push_back(tshift.cart());
      }
    }


    if(prim_flag) {
      new_prim = (*this);
      return true;
    }


    shift.push_back(lattice()[0]);
    shift.push_back(lattice()[1]);
    shift.push_back(lattice()[2]);

    //We want to minimize the volume of the primitivized cell, but to make it not a weird shape
    //that leads to noise we also minimize the dot products like reduced cell would
    min_vol = std::abs(lattice().vol());
    for(sh = 0; sh < shift.size(); sh++) {
      for(sh1 = sh + 1; sh1 < shift.size(); sh1++) {
        for(sh2 = sh1 + 1; sh2 < shift.size(); sh2++) {
          tvol = std::abs(triple_prod(shift[sh], shift[sh1], shift[sh2]));
          if(tvol < min_vol && tvol > prim_vol_tol) {
            min_vol = tvol;
            prim_vec0 = shift[sh];
            prim_vec1 = shift[sh1];
            prim_vec2 = shift[sh2];
          }
        }
      }

    }


    Lattice new_lat(prim_vec0, prim_vec1, prim_vec2);
    Lattice reduced_new_lat = niggli(new_lat, lattice().tol());
    //The lattice so far is OK, but it's noisy enough to matter for large
    //superstructures. We eliminate the noise by reconstructing it now via
    //rounded to integer transformation matrix.

    Eigen::Matrix3d transmat, invtransmat, reduced_new_lat_mat;
    SymGroup pgroup;
    reduced_new_lat.generate_point_group(pgroup);

    //Do not check to see if it returned true, it very well may not!
    lattice().is_supercell_of(reduced_new_lat, pgroup, transmat);
    //Round transformation elements to integers
    for(int i = 0; i < 3; i++) {
      for(int j = 0; j < 3; j++) {
        transmat(i, j) = floor(transmat(i, j) + 0.5);
      }
    }
    invtransmat = transmat.inverse();
    reduced_new_lat_mat = lattice().lat_column_mat();
    //When constructing this, why are we using *this as the primitive cell? Seems like I should only specify the vectors
    Lattice reconstructed_reduced_new_lat(reduced_new_lat_mat * invtransmat);
    reconstructed_reduced_new_lat.make_right_handed();
    //Lattice reconstructed_reduced_new_lat(reduced_new_lat_mat*invtransmat,lattice);

    new_prim.set_lattice(reconstructed_reduced_new_lat, CART);
    CoordType tsite(new_prim.lattice());
    for(Index nb = 0; nb < basis.size(); nb++) {
      tsite = basis[nb];
      tsite.set_lattice(new_prim.lattice(), CART);
      if(new_prim.find(tsite) == new_prim.basis.size()) {
        tsite.within();
        new_prim.basis.push_back(tsite);
      }
    }

    return false;
  }

  //***********************************************************

  template<typename CoordType>
  void BasicStructure<CoordType>::set_site_internals() {
    Index nb;


    for(nb = 0; nb < basis.size(); nb++) {
      basis[nb].set_basis_ind(nb);
    }

  }

  //***********************************************************

  template<typename CoordType> template<typename CoordType2>
  Index BasicStructure<CoordType>::find(const CoordType2 &test_site) const {
    for(Index i = 0; i < basis.size(); i++) {
      if(basis[i].compare(test_site)) {
        return i;
      }
    }
    return basis.size();
  }

  //***********************************************************

  template<typename CoordType> template<typename CoordType2>
  Index BasicStructure<CoordType>::find(const CoordType2 &test_site, const Coordinate &shift) const {
    for(Index i = 0; i < basis.size(); i++) {
      if(basis[i].compare(test_site, shift)) {
        return i;
      }
    }
    return basis.size();
  }

  //***********************************************************

  template<typename CoordType>
  void BasicStructure<CoordType>::set_lattice(const Lattice &new_lat, COORD_TYPE mode) {

    m_lattice = new_lat;

    for(Index nb = 0; nb < basis.size(); nb++) {
      basis[nb].set_lattice(lattice(), mode);
    }
  }

  //\Liz D 032514
  //***********************************************************
  /**
   * Allows for the basis elements of a basic structure to be
   * manually set, e.g. as in jsonParser.cc.
   */
  //***********************************************************


  template<typename CoordType>
  void BasicStructure<CoordType>::set_basis(Array<CoordType> basis_in) {
    basis = basis_in;
    set_site_internals();
  }




  //***********************************************************
  /**
   * Given a symmetry group, the basis of the structure will have
   * each operation applied to it. The resulting set of basis
   * from performing these operations will be averaged out,
   * yielding a new average basis that will replace the current one.
   */
  //***********************************************************

  template<typename CoordType>
  void BasicStructure<CoordType>::symmetrize(const SymGroup &relaxed_factors) {
    //First make a copy of your current basis
    //This copy will eventually become the new average basis.
    Array<CoordType> avg_basis = basis;

    //Loop through given symmetry group an fill a temporary "operated basis"
    Array<CoordType> operbasis;

    //assume identity comes first, so we skip it
    for(Index rf = 0; rf < relaxed_factors.size(); rf++) {
      operbasis.clear();
      for(Index b = 0; b < basis.size(); b++) {
        operbasis.push_back(relaxed_factors[rf]*basis[b]);
      }
      //Now that you have a transformed basis, find the closest mapping of atoms
      //Then average the distance and add it to the average basis
      for(Index b = 0; b < basis.size(); b++) {
        double smallest = 1000000;
        Coordinate bshift(lattice()), tshift(lattice());
        for(Index ob = 0; ob < operbasis.size(); ob++) {
          double dist = operbasis[ob].min_dist(basis[b], tshift);
          if(dist < smallest) {
            bshift = tshift;
            smallest = dist;
          }
        }
        bshift.cart() *= (1.0 / relaxed_factors.size());
        avg_basis[b] += bshift;
      }

    }

    return;
  }

  //***********************************************************
  /**
   * Same as the other symmetrize routine, except this one assumes
   * that the symmetry group you mean to use is the factor group
   * of your structure within a certain tolerance.
   *
   * Notice that the tolerance is also used on your point group!!
   */
  //***********************************************************


  template<typename CoordType>
  void BasicStructure<CoordType>::symmetrize(const double &tolerance) {
    SymGroup factor_group;
    double orig_tol = lattice().tol();
    m_lattice.set_tol(tolerance);
    generate_factor_group(factor_group);
    symmetrize(factor_group);
    m_lattice.set_tol(orig_tol);
    return;
  }


  //************************************************************
  /// Counts sites that allow vacancies
  template<typename CoordType>
  Index BasicStructure<CoordType>::max_possible_vacancies()const {
    Index result(0);
    for(Index i = 0; i < basis.size(); i++) {
      if(basis[i].contains("Va"))
        ++result;
    }
    return result;
  }

  //************************************************************
  //read a POSCAR like file and collect all the structure variables
  //modified to read PRIM file and determine which basis to use
  //Changed by Ivy to read new VASP POSCAR format

  template<typename CoordType>
  void BasicStructure<CoordType>::read(std::istream &stream) {
    int i, t_int;
    char ch;
    Array<double> num_elem;
    Array<std::string> elem_array;
    bool read_elem = false;
    std::string tstr;
    std::stringstream tstrstream;

    CoordType tsite(lattice());

    SD_flag = false;
    getline(stream, title);
    if(title.back() == '\r')
      throw std::runtime_error(std::string("Structure file is formatted for DOS. Please convert to Unix format. (This can be done with the dos2unix command.)"));

    m_lattice.read(stream);

    stream.ignore(100, '\n');

    //Search for Element Names
    ch = stream.peek();
    while(ch != '\n' && !stream.eof()) {
      if(isalpha(ch)) {
        read_elem = true;
        stream >> tstr;
        elem_array.push_back(tstr);
        ch = stream.peek();
      }
      else if(ch == ' ' || ch == '\t') {
        stream.ignore();
        ch = stream.peek();
      }
      else if(ch >= '0' && ch <= '9') {
        break;
      }
      else {
        throw std::runtime_error(std::string("Error attempting to read Structure. Error reading atom names."));
      }
    }

    if(read_elem == true) {
      stream.ignore(10, '\n');
      ch = stream.peek();
    }

    //Figure out how many species
    int num_sites = 0;
    while(ch != '\n' && !stream.eof()) {
      if(ch >= '0' && ch <= '9') {
        stream >> t_int;
        num_elem.push_back(t_int);
        num_sites += t_int;
        ch = stream.peek();
      }
      else if(ch == ' ' || ch == '\t') {
        stream.ignore();
        ch = stream.peek();
      }
      else {
        throw std::runtime_error(std::string("Error in line 6 of structure input file. Line 6 of structure input file should contain the number of sites."));
      }
    }
    stream.get(ch);

    // fractional coordinates or cartesian
    COORD_MODE input_mode(FRAC);

    stream.get(ch);
    while(ch == ' ' || ch == '\t') {
      stream.get(ch);
    }

    if(ch == 'S' || ch == 's') {
      SD_flag = true;
      stream.ignore(1000, '\n');
      while(ch == ' ' || ch == '\t') {
        stream.get(ch);
      }
      stream.get(ch);
    }

    if(ch == 'D' || ch == 'd') {
      input_mode.set(FRAC);
    }
    else if(ch == 'C' || ch == 'c') {
      input_mode.set(CART);
    }
    else if(!SD_flag) {
      throw std::runtime_error(std::string("Error in line 7 of structure input file. Line 7 of structure input file should specify Direct, Cartesian, or Selective Dynamics."));
    }
    else if(SD_flag) {
      throw std::runtime_error(std::string("Error in line 8 of structure input file. Line 8 of structure input file should specify Direct or Cartesian when Selective Dynamics is on."));
    }

    stream.ignore(1000, '\n');
    //Clear basis if it is not empty
    if(basis.size() != 0) {
      std::cerr << "The structure is going to be overwritten." << std::endl;
      basis.clear();
    }

    if(read_elem) {
      int j = -1;
      int sum_elem = 0;
      basis.reserve(num_sites);
      for(i = 0; i < num_sites; i++) {
        if(i == sum_elem) {
          j++;
          sum_elem += num_elem[j];
        }

        tsite.read(stream, elem_array[j], SD_flag);
        basis.push_back(tsite);
      }
    }
    else {
      //read the site info
      basis.reserve(num_sites);
      for(i = 0; i < num_sites; i++) {
        tsite.read(stream, SD_flag);
        if((stream.rdstate() & std::ifstream::failbit) != 0) {
          std::cerr << "Error reading site " << i + 1 << " from structure input file." << std::endl;
          exit(1);
        }
        basis.push_back(tsite);
      }
    }

    // Check whether there are additional sites listed in the input file
    std::string s;
    getline(stream, s);
    std::istringstream tmp_stream(s);
    Eigen::Vector3d coord;
    tmp_stream >> coord;
    if(tmp_stream.good()) {
      throw std::runtime_error(std::string("ERROR: too many sites listed in structure input file."));
    }

    update();
    return;

  }

  //***********************************************************

  template<typename CoordType>
  void BasicStructure<CoordType>::print_xyz(std::ostream &stream) const {
    stream << basis.size() << '\n';
    stream << title << '\n';
    stream.precision(7);
    stream.width(11);
    stream.flags(std::ios::showpoint | std::ios::fixed | std::ios::right);

    for(Index i = 0; i < basis.size(); i++) {
      stream << std::setw(2) << basis[i].occ_name() << " ";
      stream << std::setw(12) << basis[i].cart() << '\n';
    }

  }

  //***********************************************************
  /*
    template<typename CoordType>
    void BasicStructure<CoordType>::print_cif(std::ostream &stream) const {
      const char quote = '\'';
      const char indent[] = "   ";

      //double amag, bmag, cmag;
      //double alpha, beta, gamma;

      // Copying format based on VESTA .cif output.

      // Heading text.

      stream << '#';
      for(int i = 0; i < 70; i++) {
        stream << '=';
      }
      stream << "\n\n";
      stream << "# CRYSTAL DATA\n\n";
      stream << '#';
      for(int i = 0; i < 70; i++) {
        stream << '-';
      }
      stream << "\n\n";
      stream << "data_CASM\n\n\n";

      stream.precision(5);
      stream.width(11);
      stream.flags(std::ios::showpoint | std::ios::fixed | std::ios::left);

      stream << std::setw(40) << "_pd_phase_name" << quote << title << quote << '\n';
      stream << std::setw(40) << "_cell_length_a" << lattice().lengths[0] << '\n';
      stream << std::setw(40) << "_cell_length_b" << lattice().lengths[1] << '\n';
      stream << std::setw(40) << "_cell_length_c" << lattice().lengths[2] << '\n';
      stream << std::setw(40) << "_cell_angle_alpha" << lattice().angles[0] << '\n';
      stream << std::setw(40) << "_cell_angle_beta" << lattice().angles[1] << '\n';
      stream << std::setw(40) << "_cell_angle_gamma" << lattice().angles[2] << '\n';
      stream << std::setw(40) << "_symmetry_space_group_name_H-M" << quote << "TBD" << quote << '\n';
      stream << std::setw(40) << "_symmetry_Int_Tables_number" << "TBD" << "\n\n";

      stream << "loop_\n";
      stream << "_symmetry_equiv_pos_as_xyz\n";

      // Equivalent atom positions here. Form: 'x, y, z', '-x, -y, -z', 'x+1/2, y+1/2, z', etc.
      // Use stream << indent << etc.

      stream << '\n';
      stream << "loop_\n";
      stream << indent << "_atom_site_label" << '\n';
      stream << indent << "_atom_site_occupancy" << '\n';
      stream << indent << "_atom_site_fract_x" << '\n';
      stream << indent << "_atom_site_fract_y" << '\n';
      stream << indent << "_atom_site_fract_z" << '\n';
      stream << indent << "_atom_site_adp_type" << '\n';
      stream << indent << "_atom_site_B_iso_or_equiv" << '\n';
      stream << indent << "_atom_site_type_symbol" << '\n';

      // Use stream << indent << etc.
    }
  */
  //***********************************************************

  template<typename CoordType>
  BasicStructure<CoordType> &BasicStructure<CoordType>::operator+=(const Coordinate &shift) {

    for(Index i = 0; i < basis.size(); i++) {
      basis[i] += shift;
    }

    //factor_group += shift;
    //asym_unit += shift;
    return (*this);
  }


  //***********************************************************

  template<typename CoordType>
  BasicStructure<CoordType> &BasicStructure<CoordType>::operator-=(const Coordinate &shift) {

    for(Index i = 0; i < basis.size(); i++) {
      basis[i] -= shift;
    }
    //factor_group -= shift;
    //asym_unit -= shift;
    return (*this);
  }

  //***********************************************************

  template<typename CoordType>
  BasicStructure<CoordType> operator*(const SymOp &LHS, const BasicStructure<CoordType> &RHS) { //AAB

    return BasicStructure<CoordType>(RHS).apply_sym(LHS);
  }

  //***********************************************************

  template<typename CoordType>
  BasicStructure<CoordType> operator*(const Lattice &LHS, const BasicStructure<CoordType> &RHS) {
    BasicStructure<CoordType> tsuper(LHS);
    tsuper.fill_supercell(RHS);
    return tsuper;
  }

  //****************************************************

  template<typename CoordType>
  jsonParser &BasicStructure<CoordType>::to_json(jsonParser &json) const {
    json.put_obj();

    // std::string title;
    json["title"] = title;

    // Lattice lattice;
    json["lattice"] = lattice();

    // Array<CoordType> basis;
    json["basis"] = basis;

    return json;
  }

  //****************************************************

  // Assumes constructor CoordType::CoordType(Lattice) exists
  template<typename CoordType>
  void BasicStructure<CoordType>::from_json(const jsonParser &json) {
    try {

      // std::string title;
      CASM::from_json(title, json["title"]);

      // Lattice lattice;
      CASM::from_json(m_lattice, json["lattice"], lattice().tol());

      // Array<CoordType> basis;
      basis.clear();
      CoordType coordtype(lattice());
      for(int i = 0; i < json["basis"].size(); i++) {
        CASM::from_json(coordtype, json["basis"][i]);
        basis.push_back(coordtype);
      }

    }
    catch(...) {
      /// re-throw exceptions
      throw;
    }

  }

  //****************************************************

  template<typename CoordType>
  jsonParser &to_json(const BasicStructure<CoordType> &basic, jsonParser &json) {
    return basic.to_json(json);
  }

  //****************************************************

  // Assumes constructor CoordType::CoordType(Lattice) exists
  template<typename CoordType>
  void from_json(BasicStructure<CoordType> &basic, const jsonParser &json) {
    basic.from_json(json);
  }

}

