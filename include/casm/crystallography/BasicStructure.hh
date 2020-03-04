#ifndef BASICSTRUCTURE_HH
#define BASICSTRUCTURE_HH

#include <iostream>
#include <vector>
#include <string>
#include <cstdlib>
#include <cmath>

#include "casm/crystallography/Adapter.hh"
#include "casm/global/enum.hh"
#include "casm/crystallography/Lattice.hh"
#include "casm/crystallography/Site.hh"
#include "casm/basis_set/DoFSet.hh"

namespace CASM {
  namespace xtal {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    class Coordinate;
    class UnitCellCoord;
    class Molecule;
    //  class DoFSet;

    /** \defgroup Structure
     *  \ingroup Crystallography
     *  \brief Specifies the lattice and basis of a crystal
     *  @{
     */

    ///\brief BasicStructure specifies the lattice and atomic basis of a crystal
    class BasicStructure {
    protected:
      Lattice m_lattice;

      /// User-specified name of this Structure
      std::string m_title;

      /// Lattice vectors that specifies periodicity of the crystal
      std::vector<Site> m_basis;

      /// continuous global degrees of freedom
      std::map <DoFKey, CASM::DoFSet> m_global_dof_map;

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

      virtual ~BasicStructure() {};

      //  ****Inspectors/Accessors****

      const Lattice &lattice() const {
        return m_lattice;
      }

      const std::vector<Site> &basis() const {
        return m_basis;
      }

      std::vector<Site> &set_basis() {
        reset();
        return m_basis;
      }
      const std::string &title() const {
        return m_title;
      }

      CASM::DoFSet const &global_dof(std::string const &dof_type) const;

      std::map<DoFKey, CASM::DoFSet> const &global_dofs() const {
        return m_global_dof_map;
      }

      /// Have to explicitly define the assignment operator so that sites in this structure
      /// do not depend on the lattice of 'RHS'
      virtual BasicStructure &operator=(const BasicStructure &RHS);

      // clears site_internals and does within()
      virtual void reset();

      /// Translate all basis sites so that they are inside the unit cell
      void within();

      //Site site(const UnitCellCoord &ucc) const;

      ///change the lattice and update site coordinates.  Argument 'mode' specifies which mode is preserved
      /// e.g.: struc.set_lattice(new_lat, CART) calculates all Cartesian coordinates,
      ///       invalidates the FRAC coordinates, and changes the lattice
      void set_lattice(const Lattice &lattice, COORD_TYPE mode);

      /// Set the title of the structure
      void set_title(std::string const &_title);

      /// Manually set the global DoFs
      void set_global_dofs(std::map <DoFKey, CASM::DoFSet> const &new_dof_map) {
        m_global_dof_map = new_dof_map;
      }

      /// Manually set the basis sites
      void set_basis(std::vector<Site> const &_basis, COORD_TYPE mode = CART);

      /// Clear the basis atoms
      void clear_basis();

      /// Manually set the basis sites
      void push_back(Site const &_site, COORD_TYPE mode = CART);

      //  - Symmetry

      /// \brief Returns true if structure has attributes affected by time reversal
      bool is_time_reversal_active() const;

      /// fill an empty structure with the basis of its corresponding primitive cell
      void fill_supercell(const BasicStructure &prim); //Ivy

      ///  Shortcut routine to create a supercell structure and fill it with sites
      BasicStructure create_superstruc(const Lattice &scel_lat) const;

      ///Translates all atoms in cell
      BasicStructure &operator+=(const Coordinate &shift);
      BasicStructure &operator-=(const Coordinate &shift);

      /// Counts sites that allow vacancies
      Index max_possible_vacancies()const;

      //CASM canonical input/output
      virtual void read(std::istream &stream);  //John do this

      /// Output other formats
      void print_xyz(std::ostream &stream, bool frac = false) const;
    };

    BasicStructure operator*(const Lattice &LHS, const BasicStructure &RHS);

    //Translation operators -- not yet defined
    BasicStructure operator+(const Coordinate &LHS, const BasicStructure &RHS);

    BasicStructure operator+(const BasicStructure &LHS, const Coordinate &RHS);

    std::vector<UnitCellCoord> symop_site_map(SymOp const &_op, BasicStructure const &_struc);
    template<typename ExternSymOp>
    std::vector<UnitCellCoord> symop_site_map(ExternSymOp const &_op, BasicStructure const &_struc) {
      return symop_site_map(adapter::Adapter<SymOp, ExternSymOp>()(_op), _struc);
    }

    std::vector<UnitCellCoord> symop_site_map(SymOp const &_op, BasicStructure const &_struc, double _tol);
    template<typename ExternSymOp>
    std::vector<UnitCellCoord> symop_site_map(ExternSymOp const &_op, BasicStructure const &_struc, double _tol) {
      return symop_site_map(adapter::Adapter<SymOp, ExternSymOp>()(_op), _struc, _tol);
    }

    /// Returns an Array of each *possible* Molecule in this Structure
    std::vector<Molecule> struc_molecule(BasicStructure const &_struc);

    /// Returns an Array of each *possible* AtomSpecie in this Structure
    std::vector<std::string> struc_species(BasicStructure const &_struc);

    /// Returns an Array of each *possible* Molecule in this Structure
    std::vector<std::string> struc_molecule_name(BasicStructure const &_struc);

    /// Returns an Array of each *possible* Molecule in this Structure
    std::vector<std::vector<std::string> > allowed_molecule_unique_names(BasicStructure const &_struc);

    /// Returns a vector with a list of allowed molecule names at each site
    std::vector<std::vector<std::string> > allowed_molecule_names(BasicStructure const &_struc);

    //************************************************************
    // Assumes constructor Site::Site(Lattice) exists
    //void from_json(BasicStructure &basic, const jsonParser &json);

    std::vector<DoFKey> all_local_dof_types(BasicStructure const &_struc);

    std::vector<DoFKey> continuous_local_dof_types(BasicStructure const &_struc);

    std::vector<DoFKey> global_dof_types(BasicStructure const &_struc);

    std::map<DoFKey, Index> local_dof_dims(BasicStructure const &_struc);

    std::map<DoFKey, Index> global_dof_dims(BasicStructure const &_struc);

    std::map<DoFKey, DoFSetInfo> global_dof_info(BasicStructure const &_struc);

    std::map<DoFKey, std::vector<DoFSetInfo> > local_dof_info(BasicStructure const &_struc);

    Index local_dof_dim(DoFKey const &_name, BasicStructure const &_struc);

    /** @} */
  }
}

#endif
