#ifndef BASICSTRUCTURE_HH
#define BASICSTRUCTURE_HH

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <istream>
#include <string>
#include <vector>

#include "casm/crystallography/Adapter.hh"
#include "casm/crystallography/DoFDecl.hh"
#include "casm/crystallography/DoFSet.hh"
#include "casm/crystallography/Lattice.hh"
#include "casm/crystallography/Site.hh"
#include "casm/global/enum.hh"

namespace CASM {
namespace xtal {

class Coordinate;
class UnitCellCoord;
class Molecule;

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
  std::map<DoFKey, DoFSet> m_global_dof_map;

  /// Names to distinguish occupants with the same chemical name
  std::vector<std::vector<std::string>> m_unique_names;

 private:
  void main_print(std::ostream &stream, COORD_TYPE mode, bool version5,
                  int option) const;

  // TODO: Extract
  // CASM canonical input/output
  void read(std::istream &stream, double tol = TOL);  // John do this

 public:
  // ****Constructors****
  BasicStructure(const Lattice &init_lat) : m_lattice(init_lat){};
  BasicStructure() : m_lattice(){};  // added by Ivy (do we need/want this??)

  /// Have to explicitly define the copy constructor so that sites in the new
  /// structure do not depend on the lattice of 'RHS'
  BasicStructure(const BasicStructure &RHS);

  static BasicStructure from_poscar_stream(std::istream &poscar_stream,
                                           double tol = TOL);

  ~BasicStructure(){};

  //  ****Inspectors/Accessors****

  const Lattice &lattice() const { return m_lattice; }

  const std::vector<Site> &basis() const { return m_basis; }

  std::vector<Site> &set_basis() { return m_basis; }

  const std::string &title() const { return m_title; }

  DoFSet const &global_dof(std::string const &dof_type) const;

  std::map<DoFKey, DoFSet> const &global_dofs() const {
    return m_global_dof_map;
  }

  /// Have to explicitly define the assignment operator so that sites in this
  /// structure do not depend on the lattice of 'RHS'
  BasicStructure &operator=(const BasicStructure &RHS);

  /// Translate all basis sites so that they are inside the unit cell
  void within();

  /// change the lattice and update site coordinates.  Argument 'mode' specifies
  /// which mode is preserved
  /// e.g.: struc.set_lattice(new_lat, CART) calculates all Cartesian
  /// coordinates,
  ///       invalidates the FRAC coordinates, and changes the lattice
  void set_lattice(const Lattice &lattice, COORD_TYPE mode);

  /// Set the title of the structure
  void set_title(std::string const &_title);

  /// Manually set the global DoFs
  void set_global_dofs(std::map<DoFKey, DoFSet> const &new_dof_map) {
    m_global_dof_map = new_dof_map;
  }

  /// Manually set the global DoFs
  void set_global_dofs(std::vector<DoFSet> const &new_dof_vec) {
    m_global_dof_map = make_dofset_map(new_dof_vec);
  }

  /// Manually set the basis sites
  void set_basis(std::vector<Site> const &_basis, COORD_TYPE mode = CART);

  /// Set the names used to distinguish occupants with the same chemical name
  void set_unique_names(
      std::vector<std::vector<std::string>> const &_unique_names);

  /// Get the unique names
  std::vector<std::vector<std::string>> const &unique_names() const;

  /// Manually set the basis sites
  void push_back(Site const &_site, COORD_TYPE mode = CART);

  /// \brief Returns true if structure has DoF or properties affected by time
  /// reversal
  bool is_time_reversal_active() const;

  /// Translates all atoms in cell
  BasicStructure &operator+=(const Coordinate &shift);
  BasicStructure &operator-=(const Coordinate &shift);

  /// Counts sites that allow vacancies
  Index max_possible_vacancies() const;

  // TODO: Extract
  /// Output other formats
  void print_xyz(std::ostream &stream, bool frac = false) const;
};

/* BasicStructure operator*(const Lattice &LHS, const BasicStructure &RHS); */

/// To which site a SymOp transforms each basis site
std::vector<UnitCellCoord> symop_site_map(SymOp const &_op,
                                          BasicStructure const &_struc);

/// To which site a SymOp transforms each basis site
template <typename ExternSymOp>
std::vector<UnitCellCoord> symop_site_map(ExternSymOp const &_op,
                                          BasicStructure const &_struc) {
  return symop_site_map(adapter::Adapter<SymOp, ExternSymOp>()(_op), _struc);
}

/// To which site a SymOp transforms each basis site
std::vector<UnitCellCoord> symop_site_map(SymOp const &_op,
                                          BasicStructure const &_struc,
                                          double _tol);

/// To which site a SymOp transforms each basis site
template <typename ExternSymOp>
std::vector<UnitCellCoord> symop_site_map(ExternSymOp const &_op,
                                          BasicStructure const &_struc,
                                          double _tol) {
  return symop_site_map(adapter::Adapter<SymOp, ExternSymOp>()(_op), _struc,
                        _tol);
}

/// Returns an Array of each *possible* Molecule in this Structure
std::vector<Molecule> struc_molecule(BasicStructure const &_struc);

/// Returns an Array of each *possible* AtomSpecie in this Structure
std::vector<std::string> struc_species(BasicStructure const &_struc);

/// Returns an Array of each *possible* Molecule in this Structure
std::vector<std::string> struc_molecule_name(BasicStructure const &_struc);

/// Returns an Array of each *possible* Molecule in this Structure
std::vector<std::vector<std::string>> allowed_molecule_unique_names(
    BasicStructure const &_struc);

/// Returns a vector with a list of allowed molecule names at each site
std::vector<std::vector<std::string>> allowed_molecule_names(
    BasicStructure const &_struc);

std::vector<DoFKey> all_local_dof_types(BasicStructure const &_struc);

std::vector<DoFKey> continuous_local_dof_types(BasicStructure const &_struc);

std::vector<DoFKey> global_dof_types(BasicStructure const &_struc);

std::vector<DoFKey> all_dof_types(BasicStructure const &_struc);

std::map<DoFKey, Index> local_dof_dims(BasicStructure const &_struc);

std::map<DoFKey, Index> global_dof_dims(BasicStructure const &_struc);

Index local_dof_dim(DoFKey const &_name, BasicStructure const &_struc);

bool has_strain_dof(BasicStructure const &structure);

DoFKey get_strain_dof_key(BasicStructure const &structure);

std::string get_strain_metric(DoFKey strain_dof_key);

/// Returns 'converter' which converts site_occupant indices to 'mol_list'
/// indices:
///   mol_list_index = converter[basis_site][site_occupant_index]
std::vector<std::vector<Index>> make_index_converter(
    const BasicStructure &struc, std::vector<xtal::Molecule> mol_list);

/// Returns 'converter' which converts site_occupant indices to 'mol_name_list'
/// indices:
///   mol_name_list_index = converter[basis_site][site_occupant_index]
std::vector<std::vector<Index>> make_index_converter(
    const BasicStructure &struc, std::vector<std::string> mol_name_list);

/// Returns 'converter_inverse' which converts 'mol_name_list' indices to
/// Site::site_occupant indices:
///  site_occupant_index = converter_inverse[basis_site][mol_name_list_index]
std::vector<std::vector<Index>> make_index_converter_inverse(
    const BasicStructure &struc, std::vector<std::string> mol_name_list);

/** @} */
}  // namespace xtal
}  // namespace CASM

#endif
