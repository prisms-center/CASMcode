#include "casm/crystallography/Structure.hh"

#include <sstream>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include "casm/misc/algorithm.hh"
#include "casm/crystallography/PrimGrid.hh"
#include "casm/basis_set/DoF.hh"
#include "casm/symmetry/SymGroupRep.hh"
#include "casm/symmetry/SymBasisPermute.hh"
#include "casm/casm_io/Log.hh"


namespace CASM {
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  Structure::Structure(const fs::path &filepath) : BasicStructure<Site>() {
    if(!fs::exists(filepath)) {
      default_err_log() << "Error in Structure::Structure(const fs::path &filepath)." << std::endl;
      default_err_log() << "  File does not exist at: " << filepath << std::endl;
      exit(1);
    }
    fs::ifstream infile(filepath);

    read(infile);
  }

  //***********************************************************

  Structure::Structure(const Structure &RHS) :
    BasicStructure<Site>(RHS) {

    copy_attributes_from(RHS);

  };

  //***********************************************************

  Structure::~Structure() {}

  //***********************************************************

  Structure &Structure::operator=(const Structure &RHS) {
    BasicStructure<Site>::operator=(RHS);

    //Following gets done by base class
    //lattice = RHS.lattice;
    //basis = RHS.basis;
    //title = RHS.title;

    /*
    for(Index i = 0; i < basis.size(); i++) {
      basis[i].set_lattice(lattice);
    }
    */

    copy_attributes_from(RHS);

    return *this;
  }


  //***********************************************************

  void Structure::copy_attributes_from(const Structure &RHS) {

    BasicStructure<Site>::copy_attributes_from(RHS);

    SD_flag = RHS.SD_flag;

    basis_perm_rep_ID = RHS.basis_perm_rep_ID; //this *should* work

    m_factor_group = RHS.m_factor_group;
    m_factor_group.set_lattice(lattice());
  }

  //***********************************************************


  void Structure::generate_factor_group_slow() const {
    m_factor_group.clear();
    m_factor_group.set_lattice(lattice());
    BasicStructure<Site>::generate_factor_group_slow(m_factor_group);
    return;
  }

  //************************************************************
  void Structure::generate_factor_group() const {
    m_factor_group.clear();
    m_factor_group.set_lattice(lattice());
    BasicStructure<Site>::generate_factor_group(m_factor_group);
    return;
  }

  //************************************************************
  const MasterSymGroup &Structure::factor_group() const {
    if(!m_factor_group.size())
      generate_factor_group();
    return m_factor_group;
  }

  //************************************************************
  const SymGroup &Structure::point_group() const {
    return factor_group().point_group();
  }

  //***********************************************************

  SymGroupRep const *Structure::basis_permutation_symrep() const {
    return &(factor_group().representation(basis_permutation_symrep_ID()));
  }

  //***********************************************************

  SymGroupRepID Structure::basis_permutation_symrep_ID() const {
    if(basis_perm_rep_ID.empty())
      generate_basis_permutation_representation();

    return basis_perm_rep_ID;
  }

  //************************************************************

  /// Returns an Array of each *possible* Specie in this Structure
  std::vector<AtomSpecies> Structure::struc_species() const {

    std::vector<Molecule> tstruc_molecule = struc_molecule();
    std::vector<AtomSpecies> tstruc_species;

    Index i, j;

    //For each molecule type
    for(i = 0; i < tstruc_molecule.size(); i++) {
      // For each atomposition in the molecule
      for(j = 0; j < tstruc_molecule[i].size(); j++) {
        if(!contains(tstruc_species, tstruc_molecule[i].atom(j).species())) {
          tstruc_species.push_back(tstruc_molecule[i].atom(j).species());
        }
      }
    }

    return tstruc_species;
  }

  //************************************************************

  /// Returns an Array of each *possible* Molecule in this Structure
  std::vector<Molecule> Structure::struc_molecule() const {

    std::vector<Molecule> tstruc_molecule;
    Index i, j;

    //loop over all Sites in basis
    for(i = 0; i < basis.size(); i++) {
      //loop over all Molecules in Site
      for(j = 0; j < basis[i].site_occupant().size(); j++) {
        //Collect unique Molecules
        if(!contains(tstruc_molecule, basis[i].site_occupant()[j])) {
          tstruc_molecule.push_back(basis[i].site_occupant()[j]);
        }
      }
    }//end loop over all Sites

    return tstruc_molecule;
  }

  /// Returns an Array of each *possible* AtomSpecie in this Structure
  std::vector<std::string> Structure::struc_species_name() const {

    // get AtomSpecie allowed in struc
    std::vector<AtomSpecies> struc_spec = struc_species();

    // store AtomSpecie names in vector
    std::vector<std::string> struc_spec_name;
    for(int i = 0; i < struc_spec.size(); i++) {
      struc_spec_name.push_back(struc_spec[i].name());
    }

    return struc_spec_name;
  }

  /// Returns an Array of each *possible* Molecule in this Structure
  std::vector<std::string> Structure::struc_molecule_name() const {

    // get Molecule allowed in struc
    std::vector<Molecule> struc_mol = struc_molecule();

    // store Molecule names in vector
    std::vector<std::string> struc_mol_name;
    for(int i = 0; i < struc_mol.size(); i++) {
      struc_mol_name.push_back(struc_mol[i].name());
    }

    return struc_mol_name;
  }

  //************************************************************

  /// Returns a list of how many of each species exist in this Structure
  ///   The Specie types are ordered according to struc_species()
  Eigen::VectorXi Structure::num_each_species() const {

    std::vector<AtomSpecies> tstruc_species = struc_species();
    Eigen::VectorXi tnum_each_species = Eigen::VectorXi::Zero(tstruc_species.size());

    Index i, j;
    // For each site
    for(i = 0; i < basis.size(); i++) {
      // For each atomposition in the molecule on the site
      for(j = 0; j < basis[i].occ().size(); j++) {
        // Count the present species
        tnum_each_species(find_index(tstruc_species, basis[i].occ().atom(j).species()))++;
      }
    }

    return tnum_each_species;
  }

  //************************************************************

  /// Returns a list of how many of each molecule exist in this Structure
  ///   The molecule types are ordered according to struc_molecule()
  Eigen::VectorXi Structure::num_each_molecule() const {

    std::vector<Molecule> tstruc_molecule = struc_molecule();
    Eigen::VectorXi tnum_each_molecule = Eigen::VectorXi(tstruc_molecule.size());

    Index i;
    // For each site
    for(i = 0; i < basis.size(); i++) {
      // Count the molecule
      tnum_each_molecule(find_index(tstruc_molecule, basis[i].occ()))++;
    }

    return tnum_each_molecule;
  }




  //************************************************************
  void Structure::fg_converge(double small_tol, double large_tol, double increment) {
    BasicStructure<Site>::fg_converge(m_factor_group, small_tol, large_tol, increment);
    return;
  }

  //************************************************************
  void Structure::fg_converge(double large_tol) {
    BasicStructure<Site>::fg_converge(m_factor_group, lattice().tol(), large_tol, (large_tol - lattice().tol()) / 10.0);
    return;
  }

  //************************************************************
  /*void Structure::print_factor_group(std::ostream &stream) const {
    stream << "Factor Group of " << title << ", containing "
           << factor_group().size() << " symmetry operations:\n";

    for(Index i = 0; i < m_factor_group.size(); i++) {
      m_factor_group[i].print(stream);
    }

    return;
  }
  */
  //***********************************************************
  /**
   * It is NOT wise to use this function unless you have already
   * initialized a superstructure with lattice vectors.
   *
   * It is more wise to use the two methods that call this method:
   * Either the overloaded * operator which does:
   *  SCEL_Lattice * Prim_Structrue = New_Superstructure
   *       --- or ---
   *  New_Superstructure=Prim_Structure.create_superstruc(SCEL_Lattice);
   *
   *  Both of these will return NEW superstructures.
   */
  //***********************************************************

  void Structure::fill_supercell(const Structure &prim) {
    Index i, j;

    SymGroup latvec_pg;
    m_lattice.generate_point_group(latvec_pg);

    SD_flag = prim.SD_flag;
    PrimGrid prim_grid(prim.lattice(), lattice());

    basis.clear();
    Coordinate tcoord(lattice());

    //loop over basis sites of prim
    for(j = 0; j < prim.basis.size(); j++) {

      //loop over prim_grid points
      for(i = 0; i < prim_grid.size(); i++) {

        //push back translated basis site of prim onto superstructure basis
        basis.push_back(prim.basis[j] + prim_grid.coord(i, PRIM));

        //reset lattice for most recent superstructure Site
        //set_lattice() converts fractional coordinates to be compatible with new lattice
        basis.back().set_lattice(lattice(), CART);

        basis.back().within();
        for(Index k = 0; k < basis.size() - 1; k++) {
          if(basis[k].compare(basis.back())) {
            basis.pop_back();
            break;
          }
        }
      }
    }
    //trans_and_expand primitive factor_group
    for(i = 0; i < prim.factor_group().size(); i++) {
      if(latvec_pg.find_no_trans(prim.factor_group()[i]) == latvec_pg.size()) {
        continue;
      }
      else {
        for(Index j = 0; j < prim_grid.size(); j++) {
          Coordinate t_tau(prim.factor_group()[i].tau() + prim_grid.coord(j, SCEL).const_cart(), lattice(), CART);
          t_tau.within();
          m_factor_group.push_back(SymOp(prim.factor_group()[i].matrix(),
                                         t_tau.cart()));
        }
      }
    }
    if(m_factor_group.size() > 200) {// how big is too big? this is at least big enough for FCC conventional cell
#ifndef NDEBUG
      default_err_log() << "WARNING: You have a very large factor group of a non-primitive structure. Certain symmetry features will be unavailable.\n";
#endif
      m_factor_group.invalidate_multi_tables();
    }
    update();

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

  Structure Structure::create_superstruc(const Lattice &scel_lat) const {
    Structure tsuper(scel_lat);
    tsuper.fill_supercell(*this);
    return tsuper;
  }

  //***********************************************************
  /*shorttag allows the user to decide whether or not to print the
   * entire symmetry matrix.  If shorttag=0, then only the name of the
   * operation and the eigenvector will be printed.  If it equals anything else,
   * then the symmetry operation matrix will be printed out.
   */
  /*
    void Structure::print_site_symmetry(std::ostream &stream, COORD_TYPE mode, int shorttag = 1) {
      GenericOrbitBranch<SiteCluster> asym_unit(lattice());
      asym_unit.generate_asymmetric_unit(basis, factor_group());

      stream <<  " Printing symmetry operations that leave each site unchanged:\n \n";
      for(Index i = 0; i < asym_unit.size(); i++) {
        for(Index j = 0; j < asym_unit[i].size(); j++) {
          stream << "Site: ";
          asym_unit.at(i).at(j).at(0).print(stream); //Print the site (not cluster) for the asymmetric unit
          stream << " Total Symmetry Operations: " << asym_unit[i][j].clust_group().size() << "\n";

          for(Index nc = 0; nc < asym_unit[i][j].clust_group().size(); nc++) {
            //print_short is a new fxn I added in SymOp to not print out the whole symmetry matrix
            if(shorttag == 0) {
              stream <<  std::setw(4)  << nc + 1 << ": ";
              //asym_unit[i][j].clust_group()[nc].print_short(stream); // I turned this off because the "print_short" function does not appear to exist...
            }


            else {
              stream.flags(std::ios::left);
              stream << "\n" <<  std::setw(3)  << nc + 1 << ": ";
              stream.unsetf(std::ios::left);
              if(mode == CART)
                asym_unit[i][j].clust_group()[nc].print(stream, Eigen::Matrix3d::Identity());
              else
                asym_unit[i][j].clust_group()[nc].print(stream, lattice().inv_lat_column_mat());
            }
          }
          stream << "\n  ---------------------------------------------------------------------   \n\n";
        }
      }
      return;
    }
  */
  //***********************************************************

  void Structure::reset() {
    for(Index nb = 0; nb < basis.size(); nb++) {
      basis[nb].set_basis_ind(nb);
    }
    within();
    m_factor_group.clear();

    /** Should we also invalidate the occupants?
    for(Index i = 0; i < basis.size(); i++) {
      basis[i].site_occupant.set_value(-1);
    }
    */
    return;
  }

  //*********************************************************

  void Structure::map_superstruc_to_prim(Structure &prim) {

    Index prim_to_scel;
    Site shifted_site(prim.lattice());


    //Check that (*this) is actually a supercell of the prim
    if(!lattice().is_supercell_of(prim.lattice(), prim.point_group())) {
      default_err_log() << "*******************************************\n"
                        << "ERROR in Structure::map_superstruc_to_prim:\n"
                        << "The structure \n";
      default_err_log() << jsonParser(*this) << std::endl;
      default_err_log() << "is not a supercell of the given prim!\n";
      default_err_log() << jsonParser(prim) << std::endl;
      default_err_log() << "*******************************************\n";
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
        shifted_site += prim_grid.coord(pg, PRIM);
        shifted_site.set_lattice(lattice(), CART);
        shifted_site.within();

        // invalidate asym_ind and basis_ind because when we use
        // Structure::find, we don't want a comparison using the
        // basis_ind and asym_ind; we want a comparison using the
        // cartesian and Specie type.

        shifted_site.set_basis_ind(-1);
        prim_to_scel = find(shifted_site);

        if(prim_to_scel == basis.size()) {
          default_err_log() << "*******************************************\n"
                            << "ERROR in Structure::map_superstruc_to_prim:\n"
                            << "Cannot associate site \n"
                            << shifted_site << "\n"
                            << "with a site in the supercell basis. \n"
                            << "*******************************************\n";
          default_err_log() << "The basis_ind is "
                            << shifted_site.basis_ind() << "\n ";
          exit(2);
        }

        // Set ind_to_prim of the basis site
        //basis[prim_to_scel].ind_to_prim = pb;
      }
    }
  };

  //***********************************************************
  //
  void Structure::set_lattice(const Lattice &new_lat, COORD_TYPE mode) {
    bool is_equiv(lattice() == new_lat);

    m_lattice = new_lat;

    for(Index nb = 0; nb < basis.size(); nb++) {
      basis[nb].set_lattice(lattice(), mode);
    }

    if(is_equiv)
      m_factor_group.set_lattice(lattice());
    else
      reset();
  }

  //***********************************************************
  /**
   * Loop through basis and rearrange atoms by type. Uses bubble
   * sort algorithm by comparing integer values of the strings
   * assigned to the basis occupants.
   *
   * The basis gets sorted in a sort of alphabetical way, so be
   * mindful of any POTCARs you might have if you run this.
   */
  //***********************************************************

  void Structure::sort_basis() {

    for(Index i = 0; i < basis.size(); i++) {
      for(Index j = 0; j < basis.size() - 1; j++) {

        if(basis[j].occ_name() > basis[j + 1].occ_name()) {
          basis.swap_elem(j, j + 1);
        }
      }
    }

    return;
  }

  //***********************************************************
  /**
   * Given a symmetry group, the basis of the structure will have
   * each operation applied to it. The resulting set of basis
   * from performing these operations will be averaged out,
   * yielding a new average basis that will replace the current one.
   */
  //***********************************************************

  void Structure::symmetrize(const SymGroup &relaxed_factors) {
    //First make a copy of your current basis
    //This copy will eventually become the new average basis.
    reset();
    Array<Site> avg_basis = basis;

    //Loop through given symmetry group an fill a temporary "operated basis"
    Array<Site> operbasis;
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
    basis = avg_basis;
    //generate_factor_group();
    update();
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

  void Structure::symmetrize(const double &tolerance) {
    double orig_tol = lattice().tol();
    m_lattice.set_tol(tolerance);
    generate_factor_group();
    SymGroup g = factor_group();
    symmetrize(g);
    m_lattice.set_tol(orig_tol);
    return;
  }

  // Added by Donghee

  //***********************************************************
  /**
   * Creates N of image POSCAR files in seperate directories by
   * interpolating linearly between structures.
   * For exact interpolation, choose " LOCAL " or "1" in mode
   * For nearest-image interpolation, chosse "PERIODIC" or "0"
   *
   * Stored in Array<Structure> images;
   * images[0] = start structure, images[Nofimage+1] = end structure
   */
  //***********************************************************

  void Structure::intpol(Structure end_struc, int Nofimag, PERIODICITY_TYPE mode, Array<Structure> &images) {

    for(int m = 0; m < Nofimag + 1; m++) {

      Structure tstruc(end_struc);

      // Change title
      std::string header = "POSCAR 0";
      std::stringstream convert; // stringstream used for the conversion;
      convert << m;
      tstruc.title = header + convert.str(); // title of submiges is POSCAR0X


      for(Index i = 0 ; i < basis.size(); i++) {
        Coordinate temp(lattice());

        temp.frac() = (basis[i] - end_struc.basis[i]).const_frac() * m / (Nofimag + 1);

        tstruc.basis[i] = basis[i] + temp;

      }

      images.push_back(tstruc);

    }
    images.push_back(end_struc);

    return ;
  }



  //***********************************************************

  //This function gets the permutation representation of the
  // factor group operations of the structure. It first applies
  // the factor group operation to the structure, and then tries
  // to map the new position of the basis atom to the various positions
  // before symmetry was applied. It only checks the positions after
  // it brings the basis within the crystal.

  SymGroupRepID Structure::generate_basis_permutation_representation(bool verbose) const {

    if(factor_group().size() <= 0 || !basis.size()) {
      default_err_log() << "ERROR in BasicStructure::generate_basis_permutation_representation" << std::endl;
      default_err_log() << "You have NOT generated the factor group, or something is very wrong with your structure. I'm quitting!" << std::endl;;
      exit(1);
    }

    SymGroupRep basis_permute_group(m_factor_group);

    std::string clr(100, ' ');

    for(Index ng = 0; ng < m_factor_group.size(); ng++) {
      if(verbose) {
        if(ng % 100 == 0)
          std::cout << '\r' << clr.c_str() << '\r' << "Find permute rep for symOp " << ng << "/" << m_factor_group.size() << std::flush;
      }

      basis_permute_group.set_rep(ng, SymBasisPermute(m_factor_group[ng], *this, lattice().tol()));
    }
    // Adds the representation into the master sym group of this structure and returns the rep id
    basis_perm_rep_ID = m_factor_group.add_representation(basis_permute_group);

    //default_err_log() << "Added basis permutation rep id " << rep_id << '\n';

    if(verbose) std::cout << '\r' << clr.c_str() << '\r' << std::flush;
    return basis_perm_rep_ID;
  }

  //***********************************************************

  //Sets the occupants in the basis sites to those specified by occ_index
  void Structure::set_occs(Array <int> occ_index) {
    if(occ_index.size() != basis.size()) {
      default_err_log() << "The size of the occ index and basis index do not match!\nEXITING\n";
      exit(1);
    }
    for(Index i = 0; i < basis.size(); i++) {
      basis[i].set_occ_value(occ_index[i]);
    }
  }

  //***********************************************************
  Structure &Structure::operator+=(const Coordinate &shift) {

    for(Index i = 0; i < basis.size(); i++) {
      basis[i] += shift;
    }

    m_factor_group += shift.cart();
    return (*this);
  }


  //***********************************************************
  Structure &Structure::operator-=(const Coordinate &shift) {

    for(Index i = 0; i < basis.size(); i++) {
      basis[i] -= shift;
    }
    m_factor_group -= shift.cart();
    return (*this);
  }

  //****************************************************

  jsonParser &Structure::to_json(jsonParser &json) const {

    // class Structure : public BasicStructure<Site>
    BasicStructure<Site>::to_json(json);

    // mutable MasterSymGroup m_factor_group;
    json["factor_group"] = m_factor_group;

    // bool SD_flag;
    json["SD_flag"] = SD_flag;

    return json;
  }

  //****************************************************

  // Assumes constructor CoordType::CoordType(Lattice) exists
  void Structure::from_json(const jsonParser &json) {
    try {

      // class Structure : public BasicStructure<Site>
      BasicStructure<Site> &basic = *this;
      basic.from_json(json);

      // mutable MasterSymGroup m_factor_group;
      m_factor_group.clear();
      m_factor_group.from_json(json["factor_group"]);

      // bool SD_flag;
      CASM::from_json(SD_flag, json["SD_flag"]);

    }
    catch(...) {
      /// re-throw exceptions
      throw;
    }
    update();
  }

  //****************************************************

  jsonParser &to_json(const Structure &structure, jsonParser &json) {
    return structure.to_json(json);
  }

  void from_json(Structure &structure, const jsonParser &json) {
    structure.from_json(json);
  }

  //***********************************************************
  /*
  Structure operator*(const SymOp &LHS, const Structure &RHS) { //AAB

    return Structure(RHS).apply_sym(LHS);
  }
  */
  //***********************************************************
  Structure operator*(const Lattice &LHS, const Structure &RHS) {
    Structure tsuper(LHS);
    tsuper.fill_supercell(RHS);
    return tsuper;
  }






  /// Helper Functions


  /// Returns 'converter' which converts Site::site_occupant indices to 'mol_list' indices:
  ///   mol_list_index = converter[basis_site][site_occupant_index]
  std::vector< std::vector<Index> > make_index_converter(const Structure &struc, std::vector<Molecule> mol_list) {

    std::vector< std::vector<Index> > converter(struc.basis.size());

    for(Index i = 0; i < struc.basis.size(); i++) {
      converter[i].resize(struc.basis[i].site_occupant().size());

      for(Index j = 0; j < struc.basis[i].site_occupant().size(); j++) {
        converter[i][j] = find_index(mol_list, struc.basis[i].site_occupant()[j]);
      }
    }

    return converter;

  }

  /// Returns 'converter' which converts Site::site_occupant indices to 'mol_name_list' indices:
  ///   mol_name_list_index = converter[basis_site][site_occupant_index]
  std::vector< std::vector<Index> > make_index_converter(const Structure &struc, std::vector<std::string> mol_name_list) {

    std::vector< std::vector<Index> > converter(struc.basis.size());

    for(Index i = 0; i < struc.basis.size(); i++) {
      converter[i].resize(struc.basis[i].site_occupant().size());

      for(Index j = 0; j < struc.basis[i].site_occupant().size(); j++) {
        converter[i][j] = find_index(mol_name_list, struc.basis[i].site_occupant()[j].name());
      }
    }

    return converter;

  }

  /// Returns 'converter_inverse' which converts 'mol_name_list' indices to Site::site_occupant indices:
  ///  site_occupant_index = converter_inverse[basis_site][mol_name_list_index]
  ///
  /// If mol is not allowed on basis_site, return struc.basis[basis_site].site_occupant().size()
  std::vector< std::vector<Index> > make_index_converter_inverse(const Structure &struc, std::vector<std::string> mol_name_list) {

    std::vector< std::vector<Index> > converter_inv(struc.basis.size());

    for(Index i = 0; i < struc.basis.size(); i++) {
      converter_inv[i].resize(mol_name_list.size());

      std::vector<std::string> site_occ_name_list;
      for(Index j = 0; j < struc.basis[i].site_occupant().size(); j++) {
        site_occ_name_list.push_back(struc.basis[i].site_occupant()[j].name());
      }

      for(Index j = 0; j < mol_name_list.size(); j++) {
        converter_inv[i][j] = find_index(site_occ_name_list, mol_name_list[j]);
      }
    }

    return converter_inv;

  }

};


