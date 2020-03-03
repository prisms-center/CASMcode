#ifndef STRUCTURE_HH
#define STRUCTURE_HH

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>

#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/Site.hh"
#include "casm/crystallography/Molecule.hh"
#include "casm/symmetry/SymGroup.hh"

namespace CASM {
  namespace xtal {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    //  template<typename ClustType>
    //  class GenericOrbitree;

    //  typedef GenericOrbitree<SiteCluster> SiteOrbitree;

    /** \ingroup Structure
     *  @{
     */

    ///\brief Structure specifies the lattice and atomic basis of a crystal
    class Structure : public BasicStructure {
    protected:

      /// Group symmetry operations that map the lattice and basis of Structure onto themselves,
      /// assuming that the crystal is periodic
      MasterSymGroup m_factor_group;

      /// This holds the representation id of the permutation representation
      SymGroupRepID m_basis_perm_rep_ID;

      /// Hold the SymRepIDs for the occupant DoF, one for each of the basis sites
      //Has to be mutable because everything is const for some reason
      std::vector<SymGroupRepID> m_occupant_symrepIDs;

      /// Hold the SymRepIDs for the continuous DoF, one for each of the basis sites
      //Has to be mutable because everything is const for some reason
      std::vector<SymGroupRepID> m_site_dof_symrepIDs;

    private:

      //Flushes out every SymGroupRepID for each site (occupant DoF) and gives it a default value of identity
      void _reset_occupant_symrepIDs();
      //Flushes out every SymGroupRepID for each site (continuous DoF) and gives it a default values
      void _reset_site_dof_symrepIDs();

      void main_print(std::ostream &stream, COORD_TYPE mode, bool version5, int option) const;

      /// Obtain the basis permutation symrep and site dof symreps of factor_group
      /// sets internal m_basis_perm_rep_ID
      void _generate_basis_symreps(bool verbose = false);

      /// Obtain global dof symreps of factor_group
      void _generate_global_symreps(bool verbose = false);

      void _fg_converge(SymGroup &factor_group, double small_tol, double large_tol, double increment);

      /// determines primitive cell, finds its factor group using generate_factor_group_slow, and then
      /// expands the factor group into the supercell using lattice translations
      void generate_factor_group(); // TOL is max distance for site equivalence, in Angstr.

      /// copy all non-derived members
      void copy_attributes_from(const Structure &RHS);

      /// clears symmetry, site internals, and other attributes
      /* void reset() override; */


    public:

      //TODO:
      //This might seem weird right now, but the plan is to eventually have Structure become a composition
      //of BasicStructure and its related symmetry information (e.g. factor group, symmetry representation
      //IDs, etc).
      /// Returns constant reference to the structure
      const BasicStructure &structure() const;

      //  ****Constructors****
      //TODO: I tried getting rid of this but it seems to get used in almost every casm command right now
      Structure() : BasicStructure() {}
      explicit Structure(const BasicStructure &base);
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
      SymGroupRepID basis_permutation_symrep_ID()const;
      std::vector<SymGroupRepID> occupant_symrepIDs() const;
      std::vector<SymGroupRepID> site_dof_symrepIDs() const;

      // ****Mutators****

      //   - Basic assignment/bookkeeping

      /// Have to explicitly define the assignment operator so that sites in this structure
      /// do not depend on the lattice of 'RHS'
      Structure &operator=(const Structure &RHS);

      ///change the lattice and update site coordinates.  Argument 'mode' specifies which mode is preserved
      /// e.g.: struc.set_lattice(new_lat, CART) calculates all Cartesian coordinates,
      ///       invalidates the FRAC coordinates, and changes the lattice
      /* void set_lattice(const Lattice &lattice, COORD_TYPE mode); */

      //   - Symmetry

      /// generate factor groups for a range of tol values, prints results to screen (for now)
      void fg_converge(double large_tol);
      void fg_converge(double small_tol, double large_tol, double increment);


      /// fill an empty structure with the basis of its corresponding primitive cell - performs optimized factor_group expansion
      void fill_supercell(const Structure &prim); //Ivy

      ///  Shortcut routine to create a supercell structure and fill it with sites
      Structure create_superstruc(const Lattice &scel_lat) const;

      ///Translates all atoms in cell
      Structure &operator+=(const Coordinate &shift);
      Structure &operator-=(const Coordinate &shift);

    };

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
  }
}

#endif
