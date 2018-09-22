#ifndef BASICSTRUCTURE_HH
#define BASICSTRUCTURE_HH

#include <iostream>
#include <cstdlib>
#include <cmath>

#include "casm/crystallography/Lattice.hh"

namespace CASM {
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  class SymGroupRepID;
  class Coordinate;
  class UnitCellCoord;
  class SiteCluster;
  class MasterSymGroup;
  template<typename ClustType> class GenericOrbitree;
  typedef GenericOrbitree<SiteCluster> SiteOrbitree;

  /** \defgroup Structure
   *  \ingroup Crystallography
   *  \brief Specifies the lattice and basis of a crystal
   *  @{
   */

  ///\brief BasicStructure specifies the lattice and atomic basis of a crystal
  template<typename CoordType>
  class BasicStructure {
  protected:
    Lattice m_lattice;

    ///Specifies whether selectice dynamics is on or of for DFT calculations
    bool SD_flag;

  public: // PUBLIC DATA MEMBERS -- (long-term, at least lattice should be made private and only updated via Structure::set_lattice)
    /// User-specified name of this Structure
    std::string title;

    /// Lattice vectors that specifies periodicity of the crystal
    Array<CoordType> basis;


  private: // PRIVATE METHODS

    void main_print(std::ostream &stream, COORD_TYPE mode, bool version5, int option) const;

  public: // PUBLIC METHODS

    // ****Constructors****
    BasicStructure(const Lattice &init_lat) : m_lattice(init_lat) {};
    BasicStructure() : m_lattice() {}; //added by Ivy (do we need/want this??)
    BasicStructure(const fs::path &filepath);

    /// Have to explicitly define the copy constructor so that sites in the new structure
    /// do not depend on the lattice of 'RHS'
    BasicStructure(const BasicStructure &RHS);

    virtual ~BasicStructure() {}

    //  ****Inspectors/Accessors****

    /// return basis index of site that matches test_coord, if it is in basis
    /// otherwise, returns basis.size()
    template<typename CoordType2>
    Index find(const CoordType2 &test_site, double tol = TOL) const;

    /// return basis index of site that matches test_site+shift, if it is in basis
    /// otherwise, returns basis.size()
    template<typename CoordType2>
    Index find(const CoordType2 &test_site, const Coordinate &shift, double tol) const;

    const Lattice &lattice() const {
      return m_lattice;
    }

    /// Return the UnitCellCoord corresponding to test_site (i.e., finds the basis index and
    /// the lattice translation)
    template<typename CoordType2>
    UnitCellCoord get_unit_cell_coord(const CoordType2 &test_site, double tol = TOL)const;

    // ****Mutators****

    //   - Basic assignment/bookkeeping

    /// Have to explicitly define the assignment operator so that sites in this structure
    /// do not depend on the lattice of 'RHS'
    virtual BasicStructure &operator=(const BasicStructure &RHS);

    /// Use this is the copy interface for things that derive from BasicStructure<>
    /// It should be overloaded in derived classes so that all important attributes besides lattice, basis, and title
    /// get copied
    void copy_attributes_from(const BasicStructure &RHS);
    //virtual void copy_attributes_from(const BasicStructure &RHS); <---- should be virtual, but will have to use some sort of visitor pattern to make work

    // update does reset() first, and then calls set_site_internals
    void update();

    // clears site_internals and does within()
    virtual void reset();

    /// Associate each site with its basis index by setting its internal flags (asym_ind -> -1)
    void set_site_internals();

    /// Translate all basis sites so that they are inside the unit cell
    void within();

    CoordType get_site(const UnitCellCoord &ucc) const;

    ///change the lattice and update site coordinates.  Argument 'mode' specifies which mode is preserved
    /// e.g.: struc.set_lattice(new_lat, CART) calculates all Cartesian coordinates,
    ///       invalidates the FRAC coordinates, and changes the lattice
    void set_lattice(const Lattice &lattice, COORD_TYPE mode);

    /// Manually set the basis sites
    void set_basis(Array<CoordType> basis_in);

    //  - Symmetry

    ///apply a symmetry operation to the current structure (may change lattice vectors and order of basis atoms)
    //BasicStructure &apply_sym(const SymOp &op);  //Anna do this

    void generate_factor_group(SymGroup &factor_group, double map_tol) const;
    void generate_factor_group_slow(SymGroup &factor_group, double map_tol) const;
    void fg_converge(double small_tol, double large_tol, double increment);
    void fg_converge(SymGroup &factor_group, double small_tol, double large_tol, double increment);

    SymGroupRepID generate_basis_permutation_representation(const MasterSymGroup &factor_group, bool verbose)const;

    void symmetrize(const SymGroup &relaxed_factors);
    void symmetrize(const double &tolerance);


    /// Returns true if the structure describes a crystal primitive cell
    /// i.e., no translation smaller than a lattice vector can map the structure onto itself
    bool is_primitive(double prim_tol = TOL) const;          //Donghee do this

    /// Returns true if the structure describes a crystal primitive cell
    /// and finds the primitive cell and stores it in 'new_prim'
    bool is_primitive(BasicStructure &new_prim, double prim_tol = TOL) const;   //Donghee do this

    /// fill an empty structure with the basis of its corresponding primitive cell
    void fill_supercell(const BasicStructure &prim, double map_tol = TOL); //Ivy

    ///  Shortcut routine to create a supercell structure and fill it with sites
    BasicStructure create_superstruc(const Lattice &scel_lat, double map_tol = TOL) const;

    //John G 230913
    /// Gets clusters of every size radiating from one site and saves them to a flowertree. A garland for each site is constructed.
    void generate_flowertrees_safe(const SiteOrbitree &in_tree, Array<SiteOrbitree> &out_trees);
    void generate_flowertrees(const SiteOrbitree &in_tree, Array<SiteOrbitree> &out_trees);
    //\John G 230913

    /// Figures out which prim basis each superstructure basis corresponds to
    void map_superstruc_to_prim(BasicStructure &prim, const SymGroup &point_group);

    /// If atoms are too close together, average their distance and make them one
    void merge_sites(double maxdist); //Only for same atom types

    //Array<Array<Array<double> > > get_NN_table(const double &maxr, GenericOrbitree<GenericCluster<CoordType> > &bouquet);
    //Array<Array<Array<double> > > get_NN_table(const double &maxr);

    ///Add vacuum and shift c vector. The vacuum is always added parallel to c, and the shift vector should also be parallel to the ab plane (x,y,0)
    void add_vacuum_shift(BasicStructure &new_surface_struc, double vacuum_thickness, Eigen::Vector3d shift, COORD_TYPE mode) const;
    void add_vacuum_shift(BasicStructure &new_surface_struc, double vacuum_thickness, Coordinate shift) const;  //Because Anton thought a coordinate would be better
    ///Adds vacuum layer on top of ab plane
    void add_vacuum(BasicStructure &new_surface_struc, double vacuum_thickness) const;
    ///Translates all atoms in cell
    BasicStructure &operator+=(const Coordinate &shift);
    BasicStructure &operator-=(const Coordinate &shift);

    /// Counts sites that allow vacancies
    Index max_possible_vacancies()const;

    ///Return a reference to the lattice
    const Lattice &get_lattice() const {
      return lattice;
    };

    /// Print intpolated images in seperate directries
    //CASM canonical input/output
    virtual void read(std::istream &stream);  //John do this

    /// Output other formats
    void print_xyz(std::ostream &stream) const;
    void print_cif(std::ostream &stream) const;

    jsonParser &to_json(jsonParser &json) const;

    // Assumes constructor CoordType::CoordType(Lattice) exists
    void from_json(const jsonParser &json);
  };

  template<typename CoordType>
  BasicStructure<CoordType> operator*(const SymOp &LHS, const BasicStructure<CoordType> &RHS);

  template<typename CoordType>
  BasicStructure<CoordType> operator*(const Lattice &LHS, const BasicStructure<CoordType> &RHS);

  //Translation operators -- not yet defined
  template<typename CoordType>
  BasicStructure<CoordType> operator+(const Coordinate &LHS, const BasicStructure<CoordType> &RHS);

  template<typename CoordType>
  BasicStructure<CoordType> operator+(const BasicStructure<CoordType> &LHS, const Coordinate &RHS);

  template<typename CoordType>
  jsonParser &to_json(const BasicStructure<CoordType> &basic, jsonParser &json);


  // Assumes constructor CoordType::CoordType(Lattice) exists
  template<typename CoordType>
  void from_json(BasicStructure<CoordType> &basic, const jsonParser &json);

  /** @} */

};

#include "casm/crystallography/BasicStructure_impl.hh"

#endif
