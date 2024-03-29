#ifndef STRUCTURE_HH
#define STRUCTURE_HH

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <unordered_map>
#include <vector>

#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/Molecule.hh"
#include "casm/crystallography/Site.hh"
#include "casm/symmetry/SymGroup.hh"

namespace CASM {
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//  template<typename ClustType>
//  class GenericOrbitree;

//  typedef GenericOrbitree<SiteCluster> SiteOrbitree;

/** \ingroup Structure
 *  @{
 */

///\brief Structure specifies the lattice and atomic basis of a crystal
class Structure {
 private:
  // TODO: Collect all the symmetry information and package it into a
  // StructureSymInfo class or the like. Then store alongside the BasicStructure
  // as a shared pointer. This will save you regenerating the factor group each
  // time you make a copy of *this. It'll also make the custom copy constructor
  // and assignment operators unnecessary.

  /// Group symmetry operations that map the lattice and basis of Structure onto
  /// themselves, assuming that the crystal is periodic
  MasterSymGroup m_factor_group;

  /// This holds the representation id of the permutation representation
  SymGroupRepID m_basis_perm_rep_ID;

  /// Hold the SymRepIDs for the occupant DoF, one for each of the basis sites
  std::vector<SymGroupRepID> m_occupant_symrep_IDs;

  /// Hold the SymRepIDs for the continuous DoFs, map of DoFKey : SymGroupRepIDs
  /// for each basis site
  std::vector<std::map<DoFKey, SymGroupRepID>> m_site_dof_symrep_IDs;

  /// Holds SymRepIDs for each of the global DoFs, using the DoF name
  /// as the key to access the SymRepID value
  std::unordered_map<std::string, SymGroupRepID> m_global_dof_symrep_IDs;

  // Flushes out every SymGroupRepID for each site (occupant DoF) and gives it a
  // default value of identity
  void _reset_occupant_symrep_IDs();

  // Flushes out every SymGroupRepID for each site (continuous DoF) and gives it
  // a default values
  void _reset_site_dof_symrep_IDs();

  /// Obtain the basis permutation symrep and site dof symreps of factor_group
  /// sets internal m_basis_perm_rep_ID
  void _generate_basis_symreps();

  /// Obtain global dof symreps of factor_group
  void _generate_global_symreps();

  void _fg_converge(SymGroup &factor_group, double small_tol, double large_tol,
                    double increment);

  /// determines primitive cell, finds its factor group using
  /// generate_factor_group_slow, and then expands the factor group into the
  /// supercell using lattice translations
  void generate_factor_group();  // TOL is max distance for site equivalence, in
                                 // Angstr.

  /// copy all non-derived members
  void copy_attributes_from(const Structure &RHS);

  // Don't ever even think about making this non-const, assignable or modifiable
  // in any way.
  std::shared_ptr<const xtal::BasicStructure> m_structure_ptr;

 public:
  // TODO: Do we even want automatic conversion?
  /// Returns constant reference to the structure
  operator const xtal::BasicStructure &() const;

  const xtal::BasicStructure &structure() const {
    return *this->m_structure_ptr;
  }

  const std::shared_ptr<const xtal::BasicStructure> &shared_structure() const {
    return this->m_structure_ptr;
  }

  const Lattice &lattice() const { return this->structure().lattice(); }

  const std::vector<xtal::Site> &basis() const {
    return this->structure().basis();
  }

  // TODO: I tried getting rid of this but it seems to get used in almost every
  // casm command right now
  Structure();
  explicit Structure(const xtal::BasicStructure &base);
  explicit Structure(const fs::path &filepath);

  /// Have to explicitly define the copy constructor so that factor_group
  /// does not depend on the lattice of 'RHS'
  Structure(const Structure &RHS);

  ~Structure(){};

  const MasterSymGroup &factor_group() const;
  const SymGroup &point_group() const;
  SymGroupRep const *basis_permutation_symrep() const;
  SymGroupRepID basis_permutation_symrep_ID() const;
  std::vector<SymGroupRepID> occupant_symrep_IDs() const;
  std::vector<std::map<DoFKey, SymGroupRepID>> site_dof_symrep_IDs() const;
  SymGroupRepID global_dof_symrep_ID(const std::string dof_name) const;

  /// Have to explicitly define the assignment operator so that sites in this
  /// structure do not depend on the lattice of 'RHS'
  Structure &operator=(const Structure &RHS);
};

//************************************************************************************

struct DoFSetInfo;

std::map<DoFKey, CASM::DoFSetInfo> global_dof_info(Structure const &_struc);

std::map<DoFKey, std::vector<CASM::DoFSetInfo>> local_dof_info(
    Structure const &_struc);

class SupercellSymInfo;
SupercellSymInfo make_supercell_sym_info(Structure const &prim,
                                         Lattice const &super_lattice);

/** @} */
}  // namespace CASM

#endif
