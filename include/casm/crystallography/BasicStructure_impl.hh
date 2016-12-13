#include <sstream>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>

#include "casm/crystallography/Coordinate.hh"
#include "casm/crystallography/UnitCellCoord.hh"
#include "casm/crystallography/PrimGrid.hh"
#include "casm/symmetry/SymPermutation.hh"
#include "casm/symmetry/SymBasisPermute.hh"
#include "casm/symmetry/SymGroupRep.hh"
#include "casm/crystallography/Niggli.hh"

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
  void BasicStructure<CoordType>::generate_factor_group_slow(SymGroup &factor_group, double map_tol) const {
    //std::cout << "SLOW GENERATION OF FACTOR GROUP " << &factor_group << "\n";
    //std::cout << "begin generate_factor_group_slow() " << this << std::endl;

    Array<CoordType> trans_basis;
    Index pg, b0, b1, b2;
    Coordinate t_tau(lattice());
    Index num_suc_maps;

    SymGroup point_group;
    //reset();
    lattice().generate_point_group(point_group, map_tol);

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
        for(b1 = 0; b1 < basis.size(); b1++) { //Loop over original basis sites
          for(b2 = 0; b2 < trans_basis.size(); b2++) { //Loop over symmetrically transformed basis sites

            //see if translation successfully maps the two sites
            if(basis[b1].compare(trans_basis[b2], t_tau, map_tol)) {
              tdist = basis[b1].min_dist(Coordinate(trans_basis[b2]) + t_tau);
              if(tdist > max_error) {
                max_error = tdist;
              }
              num_suc_maps++;
              break;
            }
          }

          //break out of outer loop if inner loop finds no successful map
          if(b2 == trans_basis.size()) {
            break;
          }
        }

        //If all atoms in the basis are mapped successfully, try to add the corresponding
        //symmetry operation to the factor_group
        if(num_suc_maps == basis.size()) {
          SymOp tSym(SymOp::translation(t_tau.cart())*point_group[pg]);
          tSym.set_map_error(max_error);

          if(!factor_group.contains(tSym)) {
            factor_group.push_back(tSym);
          }
        }
      }
    } //End loop over point_group operations
    factor_group.enforce_group(map_tol);
    factor_group.sort();
    factor_group.max_error();

    //std::cout << "finish generate_factor_group_slow() " << this << std::endl;
    return;
  }

  //************************************************************

  template<typename CoordType>
  void BasicStructure<CoordType>::generate_factor_group(SymGroup &factor_group, double map_tol) const {
    //std::cout << "begin generate_factor_group() " << this << std::endl;
    BasicStructure<CoordType> tprim;
    factor_group.clear();
    factor_group.set_lattice(lattice());
    // CASE 1: Structure is primitive
    if(is_primitive(tprim, map_tol)) {
      generate_factor_group_slow(factor_group, map_tol);
      return;
    }


    // CASE 2: Structure is not primitive

    PrimGrid prim_grid(tprim.lattice(), lattice());
    //std::cout << "FAST GENERATION OF FACTOR GROUP " << &factor_group << " FOR STRUCTURE OF VOLUME " << prim_grid.size() << "\n";

    SymGroup prim_fg;
    tprim.generate_factor_group_slow(prim_fg, map_tol);

    SymGroup point_group;
    lattice().generate_point_group(point_group, map_tol);
    point_group.enforce_group(map_tol);

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

    for(double i = small_tol; i < large_tol; i += increment) {
      tols.push_back(i);
      factor_group.clear();
      generate_factor_group(factor_group, i);
      factor_group.get_multi_table();
      num_ops.push_back(factor_group.size());
      is_group.push_back(factor_group.is_group(i));
      factor_group.enforce_group(i);
      num_enforced_ops.push_back(factor_group.size());
      factor_group.character_table();
      name.push_back(factor_group.get_name());
    }

    for(Index i = 0; i < tols.size(); i++) {
      std::cout << tols[i] << "\t" << num_ops[i] << "\t" << is_group[i] << "\t" << num_enforced_ops[i] << "\t name: " << name[i] << "\n";
    }

    return;
  }

  //***********************************************************

  //This function gets the permutation representation of the
  // factor group operations of the structure. It first applies
  // the factor group operation to the structure, and then tries
  // to map the new position of the basis atom to the various positions
  // before symmetry was applied. It only checks the positions after
  // it brings the basis within the crystal.

  template<typename CoordType>
  SymGroupRepID BasicStructure<CoordType>::generate_basis_permutation_representation(const MasterSymGroup &factor_group, bool verbose) const {

    if(factor_group.size() <= 0 || !basis.size()) {
      std::cerr << "ERROR in BasicStructure::generate_basis_permutation_representation" << std::endl;
      std::cerr << "You have NOT generated the factor group, or something is very wrong with your structure. I'm quitting!" << std::endl;;
      exit(1);
    }

    SymGroupRep basis_permute_group(factor_group);
    SymGroupRepID rep_id;

    std::string clr(100, ' ');

    for(Index ng = 0; ng < factor_group.size(); ng++) {
      if(verbose) {
        if(ng % 100 == 0)
          std::cout << '\r' << clr.c_str() << '\r' << "Find permute rep for symOp " << ng << "/" << factor_group.size() << std::flush;
      }

      basis_permute_group.set_rep(ng, SymBasisPermute(factor_group[ng], *this, TOL));
    }
    // Adds the representation into the master sym group of this structure and returns the rep id
    rep_id = factor_group.add_representation(basis_permute_group);

    //std::cerr << "Added basis permutation rep id " << rep_id << '\n';

    if(verbose) std::cout << '\r' << clr.c_str() << '\r' << std::flush;
    return rep_id;
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
  void BasicStructure<CoordType>::fill_supercell(const BasicStructure<CoordType> &prim, double map_tol) {
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
  BasicStructure<CoordType> BasicStructure<CoordType>::create_superstruc(const Lattice &scel_lat, double map_tol) const {
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
  bool BasicStructure<CoordType>::is_primitive(double prim_tol) const {
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
          //if(basis[b3].compare_type(basis[b2], bshift) && tshift.min_dist(bshift) < prim_tol) {
          if(basis[b3].compare(basis[b2], tshift, prim_tol)) {
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
  bool BasicStructure<CoordType>::is_primitive(BasicStructure<CoordType> &new_prim, double prim_tol) const {
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
          if(basis[b3].compare(basis[b2], tshift, prim_tol)) {
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
    //that leads to noise we also minimize the dot products like get_reduced cell would
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
    Lattice reduced_new_lat = niggli(new_lat, prim_tol);

    //The lattice so far is OK, but it's noisy enough to matter for large
    //superstructures. We eliminate the noise by reconstructing it now via
    //rounded to integer transformation matrix.

    Eigen::Matrix3d transmat, invtransmat, reduced_new_lat_mat;
    SymGroup pgroup;
    reduced_new_lat.generate_point_group(pgroup, prim_tol);

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
      if(new_prim.find(tsite, prim_tol) == new_prim.basis.size()) {
        tsite.within();
        new_prim.basis.push_back(tsite);
      }
    }

    //std::cout<<"%%%%%%%%%%%%%%%%"<<std::endl;
    //std::cout<<"new_lat"<<std::endl;
    //new_lat.print(std::cout);
    //std::cout<<"reduced_new_lat"<<std::endl;
    //reduced_new_lat.print(std::cout);
    //std::cout<<"reconstructed_reduced_new_lat"<<std::endl;
    //reconstructed_reduced_new_lat.print(std::cout);
    //std::cout<<"transmat (rounded)"<<std::endl;
    //std::cout<<transmat<<std::endl;
    //std::cout<<"%%%%%%%%%%%%%%%%"<<std::endl;
    return false;
  }

  //***********************************************************

  template<typename CoordType>
  void BasicStructure<CoordType>::set_site_internals() {
    //std::cout << "begin set_site_internals() " << this << std::endl;
    Index nb;


    for(nb = 0; nb < basis.size(); nb++) {
      basis[nb].set_basis_ind(nb);
    }

  }

  //*********************************************************

  template<typename CoordType>
  void BasicStructure<CoordType>::map_superstruc_to_prim(BasicStructure<CoordType> &prim, const SymGroup &point_group) {

    int prim_to_scel = -1;
    CoordType shifted_site(prim.lattice());

    //Check that (*this) is actually a supercell of the prim
    if(!lattice().is_supercell_of(prim.lattice(), point_group)) {
      std::cout << "*******************************************\n"
                << "ERROR in BasicStructure<CoordType>::map_superstruc_to_prim:\n"
                << "The structure \n";
      //print(std::cout);
      std::cout << jsonParser(*this) << std::endl;
      std::cout << "is not a supercell of the given prim!\n";
      //prim.print(std::cout);
      std::cout << jsonParser(prim) << std::endl;
      std::cout << "*******************************************\n";
      exit(1);
    }

    //Get prim grid of supercell to get the lattice translations
    //necessary to stamp the prim in the superstructure
    PrimGrid prim_grid(prim.lattice(), lattice());

    // Translate each of the prim atoms by prim_grid translation
    // vectors, and map that translated atom in the supercell.
    for(Index pg = 0; pg < prim_grid.size(); pg++) {
      for(Index pb = 0; pb < prim.basis.size(); pb++) {
        shifted_site = prim.basis[pb];
        //shifted_site lattice is PRIM, so get prim_grid.coord in PRIM mode
        shifted_site += prim_grid.coord(pg, PRIM);
        shifted_site.set_lattice(lattice(), CART);
        shifted_site.within();

        // invalidate asym_ind and basis_ind because when we use
        // BasicStructure<CoordType>::find, we don't want a comparison using the
        // basis_ind and asym_ind; we want a comparison using the
        // cartesian and Specie type.

        shifted_site.set_basis_ind(-1);
        prim_to_scel = find(shifted_site);

        if(prim_to_scel == basis.size()) {
          std::cout << "*******************************************\n"
                    << "ERROR in BasicStructure<CoordType>::map_superstruc_to_prim:\n"
                    << "Cannot associate site \n"
                    << shifted_site << "\n"
                    << "with a site in the supercell basis. \n"
                    << "*******************************************\n";
          std::cout << "The basis_ind and asym_ind are "
                    << shifted_site.basis_ind() << "\t "
                    << shifted_site.asym_ind() << "\n";
          exit(2);
        }

        // Set ind_to_prim of the basis site
        basis[prim_to_scel].ind_to_prim = pb;
      }
    }
  }

  //***********************************************************

  template<typename CoordType> template<typename CoordType2>
  Index BasicStructure<CoordType>::find(const CoordType2 &test_site, double tol) const {
    for(Index i = 0; i < basis.size(); i++) {
      if(basis[i].compare(test_site, tol)) {
        return i;
      }
    }
    return basis.size();
  }

  //***********************************************************

  template<typename CoordType> template<typename CoordType2>
  Index BasicStructure<CoordType>::find(const CoordType2 &test_site, const Coordinate &shift, double tol) const {
    for(Index i = 0; i < basis.size(); i++) {
      if(basis[i].compare(test_site, shift, tol)) {
        return i;
      }
    }
    return basis.size();
  }

  //John G 070713
  //***********************************************************
  /**
   * Using the lattice of (*this), this function will return
   * a UnitCellCoord that corresponds to a site passed to it
   * within a given tolerance. This function is useful for making
   * a nearest neighbor table from sites that land outside of the
   * primitive cell.
   */
  //***********************************************************

  template<typename CoordType> template<typename CoordType2>
  UnitCellCoord BasicStructure<CoordType>::get_unit_cell_coord(const CoordType2 &bsite, double tol) const {

    CoordType2 tsite = bsite;

    tsite.set_lattice(lattice(), CART);

    Index b;

    b = find(tsite, tol);

    if(b == basis.size()) {
      std::cerr << "ERROR in BasicStructure::get_unit_cell_coord" << std::endl
                << "Could not find a matching basis site." << std::endl
                << "  Looking for: FRAC: " << tsite.const_frac() << "\n"
                << "               CART: " << tsite.const_cart() << "\n";
      exit(1);
    }

    return UnitCellCoord(b, lround(tsite.const_frac() - basis[b].const_frac()));
  };


  //*******************************************************************************************

  template<typename CoordType>
  CoordType BasicStructure<CoordType>::get_site(const UnitCellCoord &ucc) const {
    if(ucc[0] < 0 || ucc[0] >= basis.size()) {
      std::cerr << "CRITICAL ERROR: In BasicStructure<CoordType>::get_site(), UnitCellCoord " << ucc << " is out of bounds!\n"
                << "                Cannot index basis, which contains " << basis.size() << " objects.\n";
      assert(0);
      exit(1);
    }
    Coordinate trans(Eigen::Vector3d(ucc[1], ucc[2], ucc[3]), lattice(), FRAC);
    return basis[ucc[0]] + trans;
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


  //\John G 121212
  //***********************************************************
  /**
   * Goes to a specified site of the basis and makes a flower tree
   * of pairs. It then stores the length and multiplicity of
   * the pairs in a double array, giving you a strict
   * nearest neighbor table. This version also fills up a SiteOrbitree
   * in case you want to keep it.
   * Blatantly copied from Anna's routine in old new CASM
   */
  //***********************************************************
  /*
    template<typename CoordType>
    Array<Array<Array<double> > > BasicStructure<CoordType>::get_NN_table(const double &maxr, GenericOrbitree<GenericCluster<CoordType> > &bouquet) {
      if(!bouquet.size()) {
        std::cerr << "WARNING in BasicStructure<CoordType>::get_NN_table" << std::endl;
        std::cerr << "The provided GenericOrbitree<Cluster<CoordType> > is about to be rewritten!" << std::endl;
      }

      Array<Array<Array<double> > > NN;
      GenericOrbitree<GenericCluster<CoordType> > normtree(lattice());
      GenericOrbitree<GenericCluster<CoordType> > tbouquet(lattice());
      bouquet = tbouquet;
      normtree.min_num_components = 1;
      normtree.max_num_sites = 2;
      normtree.max_length.push_back(0.0);
      normtree.max_length.push_back(0.0);
      normtree.max_length.push_back(maxr);

      normtree.generate_orbitree(*this);
      normtree.print_full_clust(std::cout);
      generate_basis_bouquet(normtree, bouquet, 2);

      Array<Array<double> > oneNN;
      oneNN.resize(2);
      for(Index i = 0; i < bouquet.size(); i++) {
        NN.push_back(oneNN);
        for(Index j = 0; j < bouquet[i].size(); j++) {
          NN[i][0].push_back(bouquet[i][j].size());
          NN[i][1].push_back(bouquet[i][j].max_length());
        }
      }
      return NN;
    }
  */
  //***********************************************************
  /**
   * Goes to a specified site of the basis and makes a flower tree
   * of pairs. It then stores the length and multiplicity of
   * the pairs in a double array, giving you a strict
   * nearest neighbor table. The bouquet used for this
   * falls into the void.
   */
  //***********************************************************
  /*
    template<typename CoordType>
    Array<Array<Array<double> > > BasicStructure<CoordType>::get_NN_table(const double &maxr) {
      GenericOrbitree<GenericCluster<CoordType> > bouquet(lattice());
      return get_NN_table(maxr, bouquet);
    }
  */

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
        operbasis.back().print(std::cout);
        std::cout << std::endl;
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
    generate_factor_group(factor_group, tolerance);
    symmetrize(factor_group);
    return;
  }

  //***********************************************************
  /**
   *  Call this on a structure to get new_surface_struc: the structure with a
   *  layer of vacuum added parallel to the ab plane.
   *  vacuum_thickness: thickness of vacuum layer (Angstroms)
   *  shift:  shift vector from layer to layer, assumes FRAC unless specified.
   *  The shift vector should only have values relative to a and b vectors (eg x, y, 0).
   *  Default shift is zero.
   */
  //***********************************************************

  template<typename CoordType>
  void BasicStructure<CoordType>::add_vacuum_shift(BasicStructure<CoordType> &new_surface_struc, double vacuum_thickness, Eigen::Vector3d shift, COORD_TYPE mode) const {

    Coordinate cshift(shift, lattice(), mode);    //John G 121030
    if(!almost_zero(cshift.frac(2))) {
      std::cout << cshift.const_frac() << std::endl;
      std::cout << "WARNING: You're shifting in the c direction! This will mess with your vacuum and/or structure!!" << std::endl;
      std::cout << "See BasicStructure<CoordType>::add_vacuum_shift" << std::endl;
    }

    Eigen::Vector3d vacuum_vec;                 //unit vector perpendicular to ab plane
    vacuum_vec = lattice()[0].cross(lattice()[1]);
    vacuum_vec.normalize();
    Lattice new_lattice(lattice()[0],
                        lattice()[1],
                        lattice()[2] + vacuum_thickness * vacuum_vec + cshift.cart()); //Add vacuum and shift to c vector

    new_surface_struc = *this;
    new_surface_struc.set_lattice(new_lattice, CART);
    new_surface_struc.initialize();
    return;
  }

  //***********************************************************
  template<typename CoordType>
  void BasicStructure<CoordType>::add_vacuum_shift(BasicStructure<CoordType> &new_surface_struc, double vacuum_thickness, Coordinate shift) const {
    if(&(shift.home()) != &lattice()) {
      std::cout << "WARNING: The lattice from your shift coordinate does not match the lattice of your structure!" << std::endl;
      std::cout << "See BasicStructure<CoordType>::add_vacuum_shift" << std::endl << std::endl;
    }

    add_vacuum_shift(new_surface_struc, vacuum_thickness, shift.cart(), CART);
    return;
  }

  //***********************************************************
  template<typename CoordType>
  void BasicStructure<CoordType>::add_vacuum(BasicStructure<CoordType> &new_surface_struc, double vacuum_thickness) const {
    Eigen::Vector3d shift(0, 0, 0);

    add_vacuum_shift(new_surface_struc, vacuum_thickness, shift, FRAC);

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
      CASM::from_json(m_lattice, json["lattice"]);

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

  // Assumes constructor CoordType::CoordType(Lattice) exists
  template<typename CoordType>
  void from_json(BasicStructure<CoordType> &basic, const jsonParser &json) {
    basic.from_json(json);
  }


}

