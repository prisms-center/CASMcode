#include "casm/crystallography/Structure.hh"

#include <sstream>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>

#include "casm/clusterography/SiteCluster.hh"
#include "casm/clusterography/Orbitree.hh"
#include "casm/clusterography/jsonClust.hh"
#include "casm/misc/algorithm.hh"


namespace CASM {
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  Structure::Structure(const fs::path &filepath) : BasicStructure<Site>() {
    if(!fs::exists(filepath)) {
      std::cerr << "Error in Structure::Structure(const fs::path &filepath)." << std::endl;
      std::cerr << "  File does not exist at: " << filepath << std::endl;
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


  void Structure::generate_factor_group_slow(double map_tol) const {
    m_factor_group.clear();
    m_factor_group.set_lattice(lattice());
    BasicStructure<Site>::generate_factor_group_slow(m_factor_group, map_tol);
    return;
  }

  //************************************************************
  void Structure::generate_factor_group(double map_tol) const {
    m_factor_group.clear();
    m_factor_group.set_lattice(lattice());
    //std::cout << "GENERATING STRUCTURE FACTOR GROUP " << &m_factor_group << "\n";
    BasicStructure<Site>::generate_factor_group(m_factor_group, map_tol);
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
  std::vector<Specie> Structure::get_struc_specie() const {

    std::vector<Molecule> struc_molecule = get_struc_molecule();
    std::vector<Specie> struc_specie;

    Index i, j;

    //For each molecule type
    for(i = 0; i < struc_molecule.size(); i++) {
      // For each atomposition in the molecule
      for(j = 0; j < struc_molecule[i].size(); j++) {
        if(!contains(struc_specie, struc_molecule[i][j].specie)) {
          struc_specie.push_back(struc_molecule[i][j].specie);
        }
      }
    }

    return struc_specie;
  }

  //************************************************************

  /// Returns an Array of each *possible* Molecule in this Structure
  std::vector<Molecule> Structure::get_struc_molecule() const {

    std::vector<Molecule> struc_molecule;
    Index i, j;

    //loop over all Sites in basis
    for(i = 0; i < basis.size(); i++) {
      //loop over all Molecules in Site
      for(j = 0; j < basis[i].site_occupant().size(); j++) {
        //Collect unique Molecules
        if(!contains(struc_molecule, basis[i].site_occupant()[j])) {
          struc_molecule.push_back(basis[i].site_occupant()[j]);
        }
      }
    }//end loop over all Sites

    return struc_molecule;
  }

  /// Returns an Array of each *possible* Molecule in this Structure
  std::vector<std::string> Structure::get_struc_molecule_name() const {

    // get Molecule allowed in prim, and how many there are
    std::vector<Molecule> struc_mol = get_struc_molecule();

    // store Molecule names in vector
    std::vector<std::string> struc_mol_name;
    for(int i = 0; i < struc_mol.size(); i++) {
      struc_mol_name.push_back(struc_mol[i].name);
    }

    return struc_mol_name;
  }

  //************************************************************

  /// Returns a list of how many of each specie exist in this Structure
  ///   The Specie types are ordered according to get_struc_specie()
  Eigen::VectorXi Structure::get_num_each_specie() const {

    std::vector<Specie> struc_specie = get_struc_specie();
    Eigen::VectorXi num_each_specie = Eigen::VectorXi::Zero(struc_specie.size());

    Index i, j;
    // For each site
    for(i = 0; i < basis.size(); i++) {
      // For each atomposition in the molecule on the site
      for(j = 0; j < basis[i].occ().size(); j++) {
        // Count the present specie
        num_each_specie(find_index(struc_specie, basis[i].occ()[j].specie))++;
      }
    }

    return num_each_specie;
  }

  //************************************************************

  /// Returns a list of how many of each molecule exist in this Structure
  ///   The molecule types are ordered according to get_struc_molecule()
  Eigen::VectorXi Structure::get_num_each_molecule() const {

    std::vector<Molecule> struc_molecule = get_struc_molecule();
    Eigen::VectorXi num_each_molecule = Eigen::VectorXi(struc_molecule.size());

    Index i;
    // For each site
    for(i = 0; i < basis.size(); i++) {
      // Count the molecule
      num_each_molecule(find_index(struc_molecule, basis[i].occ()))++;
    }

    return num_each_molecule;
  }




  //************************************************************
  void Structure::fg_converge(double small_tol, double large_tol, double increment) {
    BasicStructure<Site>::fg_converge(m_factor_group, small_tol, large_tol, increment);
    return;
  }

  //************************************************************
  void Structure::fg_converge(double large_tol) {
    BasicStructure<Site>::fg_converge(m_factor_group, TOL, large_tol, (large_tol - TOL) / 10.0);
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

  void Structure::fill_supercell(const Structure &prim, double map_tol) {
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
          if(basis[k].compare(basis.back(), map_tol)) {
            basis.pop_back();
            break;
          }
        }
      }
    }
    //std::cout << "WORKING ON FACTOR GROUP " << &m_factor_group << " for structure with volume " << prim_grid.size() << ":\n";
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
      std::cerr << "WARNING: You have a very large factor group of a non-primitive structure. Certain symmetry features will be unavailable.\n";
#endif
      m_factor_group.invalidate_multi_tables();
    }
    //std::cout << "Final size is: " << m_factor_group.size() << "\n";
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

  Structure Structure::create_superstruc(const Lattice &scel_lat, double map_tol) const {
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

  void Structure::print_site_symmetry(std::ostream &stream, COORD_TYPE mode, int shorttag = 1, double tol = TOL) {
    GenericOrbitBranch<SiteCluster> asym_unit(lattice());
    asym_unit.generate_asymmetric_unit(basis, factor_group(), tol);

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

  //***********************************************************

  void Structure::reset() {
    //std::cout << "begin reset() " << this << std::endl;

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
    //std::cout << "finish reset()" << std::endl;
    return;
  }

  //***********************************************************
  /**
   * A flowertree has clusters of every size organized like the
   * branches of a regular orbitree, but every cluster in the
   * entire flowertree has the same pivot.
   */
  //***********************************************************
  void Structure::generate_flowertrees(const SiteOrbitree &in_tree, Array<SiteOrbitree> &out_trees, double tol) {
    if(out_trees.size() != 0) {
      std::cerr << "WARNING in Structure::generate_flowertrees_safe" << std::endl;
      std::cerr << "The provided array of SiteOrbitree wasn't empty! Hope nothing important was there..." << std::endl;
      out_trees.clear();
    }


    //Theres a flowertree for each basis site b
    for(Index b = 0; b < basis.size(); b++) {
      SiteOrbitree ttree(lattice(), tol);   //Weird behavior without this
      ttree.reserve(in_tree.size());
      out_trees.push_back(ttree);
      //We want clusters of every size in the tree, except the 0th term (...right guys?)
      for(Index s = 1; s < in_tree.size(); s++) {
        //Make each site basis into a cluster so we can use extract_orbits_including
        SiteCluster tsiteclust(lattice());
        tsiteclust.within();
        tsiteclust.push_back(basis[b]);
        //Make a branch where all the clusters of a given size s go
        SiteOrbitBranch tbranch(lattice());
        //Push back flower onto tree
        out_trees[b].push_back(tbranch);
        //Get one flower for a given basis site
        in_tree[s].extract_orbits_including(tsiteclust, out_trees[b].back(), tol);  //couldn't get it to work withough directly passing tree.back() (i.e. can't make branch then push_back)
      }
      out_trees[b].get_index();
      out_trees[b].collect_basis_info(*this);
    }

    return;
  }

  //***********************************************************
  // Fills an orbitree such that each OrbitBranch corresponds to a basis site and contains all orbits
  // that include that basis site
  void Structure::generate_basis_bouquet(const SiteOrbitree &in_tree, SiteOrbitree &out_tree, Index num_sites, double tol) {
    Index na, np, nb, ne;

    out_tree.clear(); //Added by Ivy 11/07/12

    //make room in out_tree for flowers
    //out_tree.reserve(basis.size());

    out_tree.lattice = in_tree.lattice;
    //out_tree.set_lattice(in_tree.lattice,FRAC);   //better?

    //get asymmetric unit
    GenericOrbitBranch<SiteCluster> asym_unit(lattice());
    asym_unit.generate_asymmetric_unit(basis, factor_group(), tol);


    for(na = 0; na < asym_unit.size(); na++) {
      asym_unit.prototype(na).clust_group().character_table();
      //std::cout << "Working on asym_unit element " << na + 1 << '\n';
      nb = out_tree.size();

      //Add new branch with asym_unit.prototype(na) as pivot
      out_tree.push_back(SiteOrbitBranch(asym_unit.prototype(na)));

      for(np = 0; np < in_tree.size(); np++) {
        if(num_sites != in_tree[np].num_sites()) {
          continue;
        }
        //std::cout << "Starting to extract petals from OrbitBranch " << np << '\n';
        in_tree[np].extract_orbits_including(asym_unit.prototype(na), out_tree.back(), tol);
      }


      for(ne = 1; ne < asym_unit[na].equivalence_map.size(); ne++) {
        out_tree.push_back(out_tree[nb]);
        out_tree.back().apply_sym(asym_unit[na].equivalence_map[ne][0]);
      }
    }
    //out_tree.get_index(); //Should not get_index because we want to keep the in_tree indices Ivy 01/15/14
    out_tree.collect_basis_info(*this);
    return;
  }

  //***********************************************************
  /**
   * Fills out_tree with OrbitBranches with pivots corresponding
   * to the different asym_unit sites.  The orbits of each branch
   * are "petals" radiating from that branch's pivot.
   *
   * @param in_tree General tree that holds clusters of all sizes
   * @param out_tree Holds all flowers (orbits) radiating from asym_unit sites
   * @param num_sites Size of clusters filling up out_tree
   */
  //***********************************************************
  // Added by Ivy 10/17/12
  void Structure::generate_asym_bouquet(const SiteOrbitree &in_tree, SiteOrbitree &out_tree, Index num_sites, double tol) {
    Index na, np;//, nb;

    out_tree.clear();

    //get asymmetric unit
    GenericOrbitBranch<SiteCluster> asym_unit(lattice());
    asym_unit.generate_asymmetric_unit(basis, factor_group(), tol);

    //make room in out_tree for flowers
    out_tree.reserve(asym_unit.size());

    out_tree.lattice = in_tree.lattice;


    for(na = 0; na < asym_unit.size(); na++) {

      //nb = out_tree.size();

      //Add new branch with asym_unit.prototype(na) as pivot
      out_tree.push_back(SiteOrbitBranch(asym_unit.prototype(na)));

      for(np = 0; np < in_tree.size(); np++) {
        if(num_sites != in_tree[np].num_sites()) {
          continue;
        }

        //prototype of na orbit
        in_tree[np].extract_orbits_including(asym_unit.prototype(na), out_tree.back(), tol);
      }

    }
    //out_tree.get_index(); //Should not get_index because we want to keep the in_tree indices Ivy 01/15/14
    out_tree.collect_basis_info(*this);
    return;


  }

  //John G 230913

  //***********************************************************
  /**
   * A flowertree is NOT a bouquet.
   * Bouquets are constructed from an orbitree (generate_basis_bouquet)
   * by taking every cluster of a given size and storing them
   * in branches, where each branch contains clusters with a
   * common pivot.
   * A flowertree has clusters of every size organized like the
   * branches of a regular orbitree, but every cluster in the
   * entire flowertree has the same pivot.
   *
   * Begin by making n bouquets, where n is the size of the
   * larges cluster in the provided orbitree. Then shuffle
   * all the bouquets together and sort them by basis site.
   * The result is an array of "flowertrees", which have been
   * made by picking out flowers of every bouquet that have
   * the same pivot. You end up with one flowertree per basis site
   * each of which has clusters from size 1 to n.
   */
  //***********************************************************
  void Structure::generate_flowertrees_safe(const SiteOrbitree &in_tree, Array<SiteOrbitree> &out_trees, double tol) { //can we make it const? It would require generate_basis_bouquet to also be const
    if(out_trees.size() != 0) {
      std::cerr << "WARNING in Structure::generate_flowertrees_safe" << std::endl;
      std::cerr << "The provided array of SiteOrbitree wasn't empty! Hope nothing important was there..." << std::endl;
      out_trees.clear();
    }

    //Determine what the largest cluster in the given orbitree is, accounting for the fact that there's a 0 point cluster
    Index max_clust_size = in_tree.size() - 1;
    //Plant all the bouquets, one for each cluster size and store them in
    Array<SiteOrbitree> bouquets;

    for(Index i = 1; i < (max_clust_size + 1); i++) {
      SiteOrbitree tbouquet(lattice(), tol);
      generate_basis_bouquet(in_tree, tbouquet, i, tol);
      bouquets.push_back(tbouquet);
    }

    //Pluck a branch off each bouquet for a given basis site and stuff it into a SiteOribitree (this is a flowertree)
    //SiteOrbitree tbouquet=in_tree;   //This is so I don't have to worry about initializing things I'm scared of
    SiteOrbitree tbouquet(lattice(), tol);

    //How many flowertrees do we need? As many as there are basis sites, i.e. bouquet[i].size();
    for(Index i = 0; i < basis.size(); i++) {
      tbouquet.clear();
      //How many trees, i.e. cluster sizes
      for(Index j = 0; j < bouquets.size(); j++) {
        //Not pushing back empty cluster!!
        tbouquet.push_back(bouquets[j][i]);
      }
      tbouquet.get_index();
      tbouquet.collect_basis_info(*this);
      out_trees.push_back(tbouquet);
    }

    return;
  }

  //\John G 230913


  //*********************************************************

  void Structure::map_superstruc_to_prim(Structure &prim) {

    Index prim_to_scel;
    Site shifted_site(prim.lattice());


    //Check that (*this) is actually a supercell of the prim
    if(!lattice().is_supercell_of(prim.lattice(), prim.point_group())) {
      std::cerr << "*******************************************\n"
                << "ERROR in Structure::map_superstruc_to_prim:\n"
                << "The structure \n";
      std::cerr << jsonParser(*this) << std::endl;
      std::cerr << "is not a supercell of the given prim!\n";
      std::cerr << jsonParser(prim) << std::endl;
      std::cerr << "*******************************************\n";
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
          std::cerr << "*******************************************\n"
                    << "ERROR in Structure::map_superstruc_to_prim:\n"
                    << "Cannot associate site \n"
                    << shifted_site << "\n"
                    << "with a site in the supercell basis. \n"
                    << "*******************************************\n";
          std::cerr << "The basis_ind is "
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
   * Returns a heterostructre. The structure this function is
   * called onto will be stretched to match a and b of the
   * structure is it passed (understruc). The returned strucure
   * will contain the basis of both structures, stacked on top
   * of each other in the c direction.
   *
   * This function only matches the ab planes.
   */
  //***********************************************************

  Structure Structure::stack_on(const Structure &understruc, bool override) const {
    Structure overstruc(*this);

    //Before doing anything check to see that lattices are aligned correctly (paralled ab planes)
    if(!override) {
      double axbangle = angle(understruc.lattice()[0].cross(understruc.lattice()[1]), overstruc.lattice()[0].cross(overstruc.lattice()[1]));

      if(!almost_zero(axbangle)) {
        std::cerr << "ERROR in Structure::stack_on" << std::endl;
        std::cerr << "Lattice mismatch! Your structures are oriented differently in space. ab planes of structures must be parallel before stacking.\
                    Redefine your structure or use Structure::align_with to fix issue." << std::endl;
        std::cerr << "Your vectors were:\n" << understruc.lattice()[0].cross(understruc.lattice()[1]) << "\n" << overstruc.lattice()[0].cross(overstruc.lattice()[1]) << std::endl;
        exit(13);
      }

      //Also check if c axes are on the same side of abplane
      //else if
      //{
      //}
    }

    else {
      std::cerr << "WARNING in Structure::stack_on" << std::endl;
      std::cerr << "You've chosen to ignore if ab planes of your structures are not parallel and/or c axes point the right way.\
                This will almost surely result in a malformed structure. Do not take the output for granted." << std::endl;
    }


    //First stretch a and b of overstruc to match a and b of understruc
    Lattice newoverstruclat(understruc.lattice()[0],
                            understruc.lattice()[1],
                            overstruc.lattice()[2]);

    overstruc.set_lattice(newoverstruclat, FRAC);


    //Create a lattice big enough to fit both structures
    Lattice heterolat(understruc.lattice()[0], //=overstruc.lattice()[0]
                      understruc.lattice()[1], //=overstruc.lattice()[1]
                      understruc.lattice()[2] + overstruc.lattice()[2]);


    //Copy understruc and set its lattice to heterolat, keeping its cartesian coordinates.
    //This leaves the basis intact but adds extra space where the overstruc basis can be placed.
    Structure heterostruc(understruc);

    heterostruc.set_lattice(heterolat, CART);

    //Fill up empty space in heterostruc with overstruc
    for(Index i = 0; i < overstruc.basis.size(); i++) {
      Site tsite(overstruc.basis[i]);
      tsite.cart() = tsite.const_cart() + understruc.lattice()[2];
      tsite.set_lattice(heterostruc.lattice(), CART);
      heterostruc.basis.push_back(tsite);

    }

    for(Index i = 0; i < heterostruc.basis.size(); i++) {
      heterostruc.basis[i].within();
    }

    heterostruc.update();

    return heterostruc;
  }
  //\John G 051112

  //John G 121212
  //***********************************************************
  /**
   * Makes reflection of structure and returns new reflected strucutre. Meant
   * for primitive structures, but can be forced on superstructures. If forced
   * on non primitive structure then the reflected structure will think it's
   * primitive. To avoid this issue reflect  primitive structure and then
   * make a reflected superlattice to fill it with.
   */
  //***********************************************************
  Structure Structure::get_reflection() const {

    Structure reflectstruc(*this);
    Eigen::Vector3d zreflect;
    zreflect << 1, 1, -1;

    // resets the lattice and reflects cartesian coordinates of the basis atoms
    reflectstruc.set_lattice(Lattice(zreflect.asDiagonal() * lattice().lat_column_mat()), FRAC);

    reflectstruc.update();

    return reflectstruc;
  }

  //***********************************************************
  /**
   * Given a minimum distance, find pair clusters within this length
   * and average out the distance between the two atoms and put
   * one atom there, eliminating the other two. Meant to be used
   * for grain boundaries.
   *
   * This function is meant for averaging atoms of the same type!
   */
  //***********************************************************
  void Structure::clump_atoms(double mindist, double tol) {

    // Warning, I don't think we want to do this here, but I'm leaving it in for now - JCT 03/07/14
    update();

    //Define Orbitree just for pairs
    SiteOrbitree siamese(lattice(), tol);
    siamese.max_num_sites = 2;
    siamese.max_length.push_back(0);
    siamese.max_length.push_back(0);
    siamese.max_length.push_back(mindist);
    siamese.min_length = 0.0;
    siamese.min_num_components = 1;
    siamese.generate_orbitree(*this);

    //Set one atom of the cluster pairs to be at the average distance and flag
    //the other one for removal later
    for(Index i = 0; i < siamese[2].size(); i++) {
      for(Index j = 0; j < siamese[2][i].size(); j++) {
        Site avgsite(siamese[2][i][j][1]);
        avgsite.cart() = (siamese[2][i][j][0].const_cart() + siamese[2][i][j][1].const_cart()) * 0.5;
        avgsite.set_lattice(lattice(), CART); //It's dumb that I need to do this
        avgsite.within();

        std::cout << "###############" << std::endl;
        std::cout << siamese[2][i][j][0].basis_ind() << "    " << siamese[2][i][j][0] << std::endl;
        std::cout << siamese[2][i][j][1].basis_ind() << "    " << siamese[2][i][j][1] << std::endl;
        std::cout << "-------" << std::endl;
        std::cout << basis[siamese[2][i][j][0].basis_ind()].basis_ind() << "    " << basis[siamese[2][i][j][0].basis_ind()] << std::endl;
        std::cout << basis[siamese[2][i][j][1].basis_ind()].basis_ind() << "    " << basis[siamese[2][i][j][1].basis_ind()] << std::endl;
        std::cout << "###############" << std::endl;


        basis[siamese[2][i][j][0].basis_ind()].set_basis_ind(-99);
        basis[siamese[2][i][j][1].basis_ind()].set_basis_ind(-99);
        basis.push_back(avgsite);
      }
    }

    std::cout << "BASIS IND:" << std::endl;
    //Remove unwanted atoms.
    int bsize = basis.size();
    for(int i = bsize - 1; i >= 0; i--) {
      std::cout << basis[i].basis_ind() << "   " << basis[i] << std::endl;
      if(basis[i].basis_ind() == Index(-99)) {
        basis.remove(i);
      }
    }

    update();

    return;
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



  //John G 050513

  //***********************************************************
  /**
   * Decorate your structure easily! Given a cluster, this
   * routine will copy your structure and replace sites of the
   * copy with sites from the cluster. Bedazzle!
   * The lattice of your structure and cluster must be related
   * by an integer matrix (supercell). Be mindful with the size
   * of your stamped structure, it has to be big enough to fit
   * the cluster.
   */
  //***********************************************************

  Structure Structure::stamp_with(SiteCluster stamp, bool lat_override, bool im_override) const {
    //The following check may not work properly
    if(!lat_override && !lattice().is_supercell_of((stamp.home()), point_group())) {
      std::cerr << "ERROR in Structure::stamp_with (are you using Structure::bedazzle?)" << std::endl;
      std::cerr << "The lattice of your cluster is not related to the lattice of the structure you want to stamp!" << std::endl;
      exit(60);
    }

    //Some sort of check should happen here to make sure your superstructure
    //is large enough to accomodate the stamp. Otherwise shit will get messed up
    //I think image_check is the way to go, but I'm not sure I'm using it correctly
    //I'm using the default value for nV (1) in image_check. This will allow the cluster
    //to touch the faces of the voronoi cell.

    if(!im_override && stamp.size() > 1 && stamp.image_check(lattice(), 1)) {    //Skip point clusters?
      std::cerr << "ERROR in Structure::stamp_with" << std::endl;
      std::cerr << "Your superstructure isn't big enough to fit your cluster!" << std::endl;
      std::cerr << "Culprit:" << std::endl;
      stamp.print(std::cerr);
      exit(62);
    }

    stamp.set_lattice(lattice(), CART);
    stamp.within();

    Structure stampedstruc = *this;
    Index sanity_count = 0;
    for(Index i = 0; i < stamp.size(); i++) {
      for(Index j = 0; j < basis.size(); j++) {
        if(Coordinate(stamp[i]) == Coordinate(basis[j])) {
          sanity_count++;
          stampedstruc.basis[j] = stamp[i];
        }
      }
    }

    if(sanity_count != stamp.size()) {  //Allow if override?
      std::cerr << "ERROR in Structure::stamp_with (are you using Structure::bedazzle?)" << std::endl;
      std::cerr << stamp.size() - sanity_count << " sites in the cluster (size " << stamp.size() << " did not map onto your structure.\
                If you got this far your lattices passed the test, but it seems your bases are unrelated. My guesses:" << std::endl;
      std::cerr << "You're trying to map a relaxed site onto an ideal one." << std::endl;
      std::cerr << "You have a vacancy in the structure you're stamping that's not in the basis." << std::endl;
      exit(61);
    }

    stampedstruc.reset();
    return stampedstruc;
  }

  //***********************************************************
  /**
   * Same as stamp_with, but takes an array of clusters and returns
   * and array of structures. Basically saves you a for loop.
   * Perhaps this should actually take an OrbitBranch?
   */
  //***********************************************************

  Array<Structure> Structure::bedazzle(Array<SiteCluster> stamps, bool lat_override, bool im_override) const {
    Array<Structure> all_decorations;
    for(Index i = 0; i < stamps.size(); i++) {
      all_decorations.push_back(stamp_with(stamps[i], lat_override, im_override));
    }

    return all_decorations;
  }

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

  Array<Array<Array<double> > > Structure::get_NN_table(const double &maxr, SiteOrbitree &bouquet, double tol) {
    if(!bouquet.size()) {
      std::cerr << "WARNING in Structure::get_NN_table" << std::endl;
      std::cerr << "The provided SiteOrbitree is about to be rewritten!" << std::endl;
    }

    Array<Array<Array<double> > > NN;
    SiteOrbitree normtree(lattice(), tol);
    SiteOrbitree tbouquet(lattice(), tol);
    bouquet = tbouquet;
    normtree.min_num_components = 1;
    normtree.max_num_sites = 2;
    normtree.max_length.push_back(0.0);
    normtree.max_length.push_back(0.0);
    normtree.max_length.push_back(maxr);

    normtree.generate_orbitree(*this);
    normtree.print_full_clust(std::cout);
    generate_basis_bouquet(normtree, bouquet, 2, tol);

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

  //***********************************************************
  /**
   * Goes to a specified site of the basis and makes a flower tree
   * of pairs. It then stores the length and multiplicity of
   * the pairs in a double array, giving you a strict
   * nearest neighbor table. The bouquet used for this
   * falls into the void.
   */
  //***********************************************************

  Array<Array<Array<double> > > Structure::get_NN_table(const double &maxr, double tol) {
    SiteOrbitree bouquet(lattice(), tol);
    return get_NN_table(maxr, bouquet, tol);
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
    generate_factor_group(tolerance);
    symmetrize(factor_group());
    return;
  }

  //\John G 050513

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

  void Structure::add_vacuum_shift(Structure &new_surface_struc, double vacuum_thickness, Eigen::Vector3d shift, COORD_TYPE mode) const {

    Coordinate cshift(shift, lattice(), mode);    //John G 121030
    if(!almost_zero(cshift.frac(2))) {
      std::cerr << cshift.const_frac() << std::endl;
      std::cerr << "WARNING: You're shifting in the c direction! This will mess with your vacuum and/or structure!!" << std::endl;
      std::cerr << "See Structure::add_vacuum_shift" << std::endl;
    }

    Eigen::Vector3d vacuum_vec;                 //unit vector perpendicular to ab plane
    vacuum_vec = lattice()[0].cross(lattice()[1]);
    vacuum_vec.normalize();
    Lattice new_lattice(lattice()[0],
                        lattice()[1],
                        lattice()[2] + vacuum_thickness * vacuum_vec + cshift.const_cart()); //Add vacuum and shift to c vector

    new_surface_struc = *this;
    new_surface_struc.set_lattice(new_lattice, CART);
    return;
  }

  //***********************************************************
  void Structure::add_vacuum_shift(Structure &new_surface_struc, double vacuum_thickness, Coordinate shift) const {

    add_vacuum_shift(new_surface_struc, vacuum_thickness, shift.cart(), CART);
    return;
  }

  //***********************************************************
  void Structure::add_vacuum(Structure &new_surface_struc, double vacuum_thickness) const {
    Eigen::Vector3d shift(0, 0, 0);

    add_vacuum_shift(new_surface_struc, vacuum_thickness, shift, FRAC);

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

  //This function gets the basis_permutation representation of the
  // factor group operations of the structure. It first applies
  // the factor group operation to the structure, and then tries
  // to map the new position of the basis atom to the various positions
  // before symmetry was applied. It only checks the positions after
  // it brings the basis within the crystal.

  // ROUTINE STILL NEEDS TO BE TESTED!

  SymGroupRepID Structure::generate_basis_permutation_representation(bool verbose) const {
    basis_perm_rep_ID = BasicStructure<Site>::generate_basis_permutation_representation(factor_group(), verbose);
    return basis_perm_rep_ID;
  }

  //***********************************************************

  //Sets the occupants in the basis sites to those specified by occ_index
  void Structure::set_occs(Array <int> occ_index) {
    if(occ_index.size() != basis.size()) {
      std::cerr << "The size of the occ index and basis index do not match!\nEXITING\n";
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
  std::vector< std::vector<Index> > get_index_converter(const Structure &struc, std::vector<Molecule> mol_list) {

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
  std::vector< std::vector<Index> > get_index_converter(const Structure &struc, std::vector<std::string> mol_name_list) {

    std::vector< std::vector<Index> > converter(struc.basis.size());

    for(Index i = 0; i < struc.basis.size(); i++) {
      converter[i].resize(struc.basis[i].site_occupant().size());

      for(Index j = 0; j < struc.basis[i].site_occupant().size(); j++) {
        converter[i][j] = find_index(mol_name_list, struc.basis[i].site_occupant()[j].name);
      }
    }

    return converter;

  }

  /// Returns 'converter_inverse' which converts 'mol_name_list' indices to Site::site_occupant indices:
  ///  site_occupant_index = converter_inverse[basis_site][mol_name_list_index]
  ///
  /// If mol is not allowed on basis_site, return struc.basis[basis_site].site_occupant().size()
  std::vector< std::vector<Index> > get_index_converter_inverse(const Structure &struc, std::vector<std::string> mol_name_list) {

    std::vector< std::vector<Index> > converter_inv(struc.basis.size());

    for(Index i = 0; i < struc.basis.size(); i++) {
      converter_inv[i].resize(mol_name_list.size());

      std::vector<std::string> site_occ_name_list;
      for(Index j = 0; j < struc.basis[i].site_occupant().size(); j++) {
        site_occ_name_list.push_back(struc.basis[i].site_occupant()[j].name);
      }

      for(Index j = 0; j < mol_name_list.size(); j++) {
        converter_inv[i][j] = find_index(site_occ_name_list, mol_name_list[j]);
      }
    }

    return converter_inv;

  }

};


