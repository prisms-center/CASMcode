#ifndef STRUCTURE_HH
#define STRUCTURE_HH

#include <iostream>
#include <cstdlib>
#include <cmath>

#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/Site.hh"

namespace CASM {
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  class MasterSymGroup;

  class SiteCluster;

  template<typename ClustType>
  class GenericOrbitree;

  typedef GenericOrbitree<SiteCluster> SiteOrbitree;

  /** \ingroup Structure
   *  @{
   */

  ///\brief Structure specifies the lattice and atomic basis of a crystal
  class Structure : public BasicStructure<Site> {
  protected: //PROTECTED DATA MEMBERS

    /// Group symmetry operations that map the lattice and basis of Structure onto themselves,
    /// assuming that the crystal is periodic
    mutable MasterSymGroup m_factor_group;
    /// This holds the representation id of the permutation representation
    mutable SymGroupRepID basis_perm_rep_ID;
    ///Specifies whether selectice dynamics is on or of for DFT calculations
    bool SD_flag;

  public: //PUBLIC DATA MEMBERS (Public for now)

    /// *** Inherited from BasicStructure<Site>
    //std::string title;     // user-specified name of this structure
    //Lattice lattice;     // Lattice vectors that specifies periodicity of the crystal
    //Array<Site> basis;   // List of basis sites of the crystal.  Each site may host a molecule or a specie



  private: //PRIVATE METHODS

    void main_print(std::ostream &stream, COORD_TYPE mode, bool version5, int option) const;


  public: //PUBLIC METHODS

    //  ****Constructors****
    Structure() : BasicStructure<Site>() {}
    explicit Structure(const Lattice &init_lat) : BasicStructure<Site>(init_lat) {}
    explicit Structure(const BasicStructure<Site> &base) : BasicStructure<Site>(base) {}
    explicit Structure(const fs::path &filepath);

    /// Have to explicitly define the copy constructor so that factor_group
    /// does not depend on the lattice of 'RHS'
    Structure(const Structure &RHS);

    //  ****Inspectors/Accessors****
    //      - non-const versions try to populate data members before access
    const MasterSymGroup &factor_group() const;
    //const MasterSymGroup &factor_group();
    const SymGroup &point_group() const;
    //const SymGroup &point_group();
    SymGroupRep const *basis_permutation_symrep()const;
    SymGroupRepID basis_permutation_symrep_ID()const;

    std::vector<Specie> get_struc_specie() const;
    std::vector<Molecule> get_struc_molecule() const;
    std::vector<std::string> get_struc_molecule_name() const;
    Eigen::VectorXi get_num_each_specie() const;
    Eigen::VectorXi get_num_each_molecule() const;

    // ****Mutators****

    //   - Basic assignment/bookkeeping

    /// Have to explicitly define the assignment operator so that sites in this structure
    /// do not depend on the lattice of 'RHS'
    Structure &operator=(const Structure &RHS);

    /// copy all non-derived members
    void copy_attributes_from(const Structure &RHS);

    /// clears symmetry, site internals, and other attributes
    void reset();

    ///change the lattice and update site coordinates.  Argument 'mode' specifies which mode is preserved
    /// e.g.: struc.set_lattice(new_lat, CART) calculates all Cartesian coordinates,
    ///       invalidates the FRAC coordinates, and changes the lattice
    void set_lattice(const Lattice &lattice, COORD_TYPE mode);

    //   - Symmetry

    /// apply a symmetry operation to the current structure (does not change lattice vectors, because
    /// 'op' is assumed to be a lattice homomorphism)
    //Structure &apply_sym(const SymOp &op);  //Anna do this

    /// determines primitive cell, finds its factor group using generate_factor_group_slow, and then
    /// expands the factor group into the supercell using lattice translations
    void generate_factor_group(double map_tol = TOL) const; // TOL is max distance for site equivalence, in Angstr.

    /// Obtain factor group by testing all operations of the lattice point_group and keep
    void generate_factor_group_slow(double map_tol = TOL) const; // TOL is max distance for site equivalence, in Angstr.

    /// generate factor groups for a range of tol values, prints results to screen (for now)
    void fg_converge(double large_tol);
    void fg_converge(double small_tol, double large_tol, double increment);

    /// Obtain the basis permutation representation of factor_group, returns its rep_id, and sets internal basis_perm_rep_ID
    SymGroupRepID generate_basis_permutation_representation(bool verbose = false) const;

    void symmetrize(const SymGroup &relaxed_factors);
    void symmetrize(const double &tolerace);


    /// Returns true if the structure describes a crystal primitive cell
    /// i.e., no translation smaller than a lattice vector can map the structure onto itself
    //bool is_primitive(double prim_tol = TOL) const;          //Donghee do this

    /// Returns true if the structure describes a crystal primitive cell
    /// and finds the primitive cell and stores it in 'new_prim'
    //bool is_primitive(Structure &new_prim, double prim_tol = TOL) const;   //Donghee do this

    /// fill an empty structure with the basis of its corresponding primitive cell - performs optimized factor_group expansion
    void fill_supercell(const Structure &prim, double map_tol = TOL); //Ivy

    ///  Shortcut routine to create a supercell structure and fill it with sites
    Structure create_superstruc(const Lattice &scel_lat, double map_tol = TOL) const;

    /// Figures out which prim basis each superstructure basis corresponds to
    void map_superstruc_to_prim(Structure &prim); //Added by Ivy 06/29/2013

    /// Setting the current occupants of the structure to those specified by an array of integers
    void set_occs(Array <int> occ_index);

    //    - Cluster related routines

    /// For each basis site, find each cluster in 'in_tree' that contain that basis site,
    /// create an orbit of all equivalent clusters that also contain the basis site
    /// and store each of these orbits in out_tree
    /// num_sites specifies what size of clusters are desired (e.g., pairs, triplets, etc)
    void generate_basis_bouquet(const SiteOrbitree &in_tree, SiteOrbitree &out_tree, Index num_sites, double tol);

    /// For each asymmetric unit site, find each cluster in 'in_tree' that contain that basis site,
    /// create an orbit of all equivalent clusters that also contain the basis site
    /// and store each of these orbits in out_tree
    /// num_sites specifies what size of clusters are desired (e.g., pairs, triplets, etc)
    void generate_asym_bouquet(const SiteOrbitree &in_tree, SiteOrbitree &out_tree, Index num_sites, double tol);

    //John G 230913
    /// Gets clusters of every size radiating from one site and saves them to a flowertree. A garland for each site is constructed.
    void generate_flowertrees_safe(const SiteOrbitree &in_tree, Array<SiteOrbitree> &out_trees, double tol = TOL);
    void generate_flowertrees(const SiteOrbitree &in_tree, Array<SiteOrbitree> &out_trees, double tol = TOL);
    //\John G 230913


    //John G 051112
    Structure stack_on(const Structure &understruc, bool override = 0) const;
    //\John G 051112

    //John G 121212
    /// Return reflection of structure
    Structure get_reflection() const;
    /// If atoms are too close together, average their distance and make them one
    void clump_atoms(double maxdist, double tol); //Only for same atom types
    /// Rearrange basis by grouping atoms by type
    void sort_basis();
    //\John G 121212

    //John G 050513
    Structure stamp_with(SiteCluster stamp, bool lat_override = 0, bool im_override = 0) const;
    Array<Structure> bedazzle(Array<SiteCluster> stamps, bool lat_override = 0, bool im_override = 0) const;
    Array<Array<Array<double> > > get_NN_table(const double &maxr, SiteOrbitree &bouquet, double tol);
    Array<Array<Array<double> > > get_NN_table(const double &maxr, double tol);
    //\John G 050513

    ///Add vacuum and shift c vector. The vacuum is always added parallel to c, and the shift vector should also be parallel to the ab plane (x,y,0)
    void add_vacuum_shift(Structure &new_surface_struc, double vacuum_thickness, Eigen::Vector3d shift, COORD_TYPE mode) const;
    void add_vacuum_shift(Structure &new_surface_struc, double vacuum_thickness, Coordinate shift) const;  //Because Anton thought a coordinate would be better
    ///Adds vacuum layer on top of ab plane
    void add_vacuum(Structure &new_surface_struc, double vacuum_thickness) const;
    ///Translates all atoms in cell
    Structure &operator+=(const Coordinate &shift);
    Structure &operator-=(const Coordinate &shift);

    /// Creates Nofimag POSCAR files by interpolating linearly between structures in star_stru and end_stru
    /// for exact interpolation, choose "LOCAL" or "1", for nearest-image interpolation, choose "PERIODIC" or "0"
    void intpol(Structure end_struc, int Nofimag, PERIODICITY_TYPE mode, Array<Structure> &images);

    // ****Input/Output****

    /// For each symmetrically distinct site, print the symmetry operations that map it onto itself
    void print_site_symmetry(std::ostream &stream, COORD_TYPE mode, int shorttag, double tol);
    //void print_factor_group(std::ostream &stream) const;

    bool read_species(); //Ivy 11/27/12
    void assign_species(Array<std::string> &names, Array<double> &masses, Array<double> &magmoms, Array<double> &Us, Array<double> &Js); //Added by Ivy

    jsonParser &to_json(jsonParser &json) const;
    void from_json(const jsonParser &json);
  };

  //Structure operator*(const SymOp &LHS, const Structure &RHS);

  jsonParser &to_json(const Structure &structure, jsonParser &json);
  void from_json(Structure &structure, const jsonParser &json);

  Structure operator*(const Lattice &LHS, const Structure &RHS);

  //Translation operators -- not yet defined
  Structure operator+(const Coordinate &LHS, const Structure &RHS);
  Structure operator+(const Structure &LHS, const Coordinate &RHS);

  //Not yet sure how these will work
  Structure operator+(const Structure &LHS, const Structure &RHS);
  Structure operator+(const Structure &LHS, const Lattice &RHS);
  Structure operator+(const Lattice &LHS, const Structure &RHS);


  /// Helper Functions

  /// Returns 'converter' which converts site_occupant indices to 'mol_list' indices:
  ///   mol_list_index = converter[basis_site][site_occupant_index]
  std::vector< std::vector<Index> > get_index_converter(const Structure &struc, std::vector<Molecule> mol_list);

  /// Returns 'converter' which converts site_occupant indices to 'mol_name_list' indices:
  ///   mol_name_list_index = converter[basis_site][site_occupant_index]
  std::vector< std::vector<Index> > get_index_converter(const Structure &struc, std::vector<std::string> mol_name_list);

  /// Returns 'converter_inverse' which converts 'mol_name_list' indices to Site::site_occupant indices:
  ///  site_occupant_index = converter_inverse[basis_site][mol_name_list_index]
  std::vector< std::vector<Index> > get_index_converter_inverse(const Structure &struc, std::vector<std::string> mol_name_list);

  /** @} */
};

//#include "casm/clusterography/Orbitree_impl.hh"

#endif
