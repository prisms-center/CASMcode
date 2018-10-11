#ifndef STRUCTURE_HH
#define STRUCTURE_HH

#include <iostream>
#include <cstdlib>
#include <cmath>

#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/Site.hh"
#include "casm/crystallography/Molecule.hh"
#include "casm/symmetry/SymGroup.hh"

namespace CASM {
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  class SiteCluster;



  //  template<typename ClustType>
  //  class GenericOrbitree;

  //  typedef GenericOrbitree<SiteCluster> SiteOrbitree;

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
    mutable SymGroupRepID m_basis_perm_rep_ID;

  private: //PRIVATE METHODS

    void main_print(std::ostream &stream, COORD_TYPE mode, bool version5, int option) const;

    /// Obtain the basis permutation symrep and site dof symreps of factor_group
    /// sets internal m_basis_perm_rep_ID
    void _generate_basis_symreps(bool verbose = false) const;

    /// Obtain global dof symreps of factor_group
    void _generate_global_symreps(bool verbose = false) const;


  public: //PUBLIC METHODS

    //  ****Constructors****
    Structure() : BasicStructure<Site>() {}
    explicit Structure(const Lattice &init_lat) : BasicStructure<Site>(init_lat) {}
    explicit Structure(const BasicStructure<Site> &base) : BasicStructure<Site>(base) {}
    explicit Structure(const fs::path &filepath);

    /// Have to explicitly define the copy constructor so that factor_group
    /// does not depend on the lattice of 'RHS'
    Structure(const Structure &RHS);

    ~Structure();

    //  ****Inspectors/Accessors****
    //      - non-const versions try to populate data members before access
    const MasterSymGroup &factor_group() const;
    //const MasterSymGroup &factor_group();
    const SymGroup &point_group() const;
    //const SymGroup &point_group();
    SymGroupRep const *basis_permutation_symrep()const;
    SymGroupRepID basis_permutation_symrep_ID()const override;

    std::vector<AtomSpecies> struc_species() const;
    std::vector<Molecule> struc_molecule() const;
    std::vector<std::string> struc_species_name() const;
    std::vector<std::string> struc_molecule_name() const;
    Eigen::VectorXi num_each_species() const;
    Eigen::VectorXi num_each_molecule() const;

    // ****Mutators****

    //   - Basic assignment/bookkeeping

    /// Have to explicitly define the assignment operator so that sites in this structure
    /// do not depend on the lattice of 'RHS'
    Structure &operator=(const Structure &RHS);

    /// copy all non-derived members
    void copy_attributes_from(const Structure &RHS);

    /// clears symmetry, site internals, and other attributes
    void reset() override;

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
    void generate_factor_group() const; // TOL is max distance for site equivalence, in Angstr.

    /// Obtain factor group by testing all operations of the lattice point_group and keep
    void generate_factor_group_slow() const; // TOL is max distance for site equivalence, in Angstr.

    /// generate factor groups for a range of tol values, prints results to screen (for now)
    void fg_converge(double large_tol);
    void fg_converge(double small_tol, double large_tol, double increment);

    void symmetrize(const SymGroup &relaxed_factors);
    void symmetrize(const double &tolerace);


    /// fill an empty structure with the basis of its corresponding primitive cell - performs optimized factor_group expansion
    void fill_supercell(const Structure &prim); //Ivy

    ///  Shortcut routine to create a supercell structure and fill it with sites
    Structure create_superstruc(const Lattice &scel_lat) const;

    /// Figures out which prim basis each superstructure basis corresponds to
    void map_superstruc_to_prim(Structure &prim); //Added by Ivy 06/29/2013

    /// Setting the current occupants of the structure to those specified by an array of integers
    void set_occs(Array <int> occ_index);

    /// Rearrange basis by grouping atoms by type
    void sort_basis();
    //\John G 121212

    ///Translates all atoms in cell
    Structure &operator+=(const Coordinate &shift);
    Structure &operator-=(const Coordinate &shift);

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
  std::vector< std::vector<Index> > make_index_converter(const Structure &struc, std::vector<Molecule> mol_list);

  /// Returns 'converter' which converts site_occupant indices to 'mol_name_list' indices:
  ///   mol_name_list_index = converter[basis_site][site_occupant_index]
  std::vector< std::vector<Index> > make_index_converter(const Structure &struc, std::vector<std::string> mol_name_list);

  /// Returns 'converter_inverse' which converts 'mol_name_list' indices to Site::site_occupant indices:
  ///  site_occupant_index = converter_inverse[basis_site][mol_name_list_index]
  std::vector< std::vector<Index> > make_index_converter_inverse(const Structure &struc, std::vector<std::string> mol_name_list);

  /** @} */
};

#endif
