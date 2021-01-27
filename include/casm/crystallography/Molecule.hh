#ifndef MOLECULE_HH
#define MOLECULE_HH

#include <array>
#include <iostream>
#include <string>
#include <vector>

#include "casm/crystallography/Adapter.hh"
#include "casm/crystallography/SpeciesAttribute.hh"
#include "casm/global/definitions.hh"
#include "casm/misc/Comparisons.hh"

namespace CASM {
namespace xtal {
class Molecule;

/** \defgroup Molecule
 *  \ingroup Crystallography
 *  \brief Relates to Molecule
 *  @{
 */

//****************************************************

/// \brief An atomic species associated with a position in space
class AtomPosition {
 public:
  /// \brief Construct with x,y,z position coordinates and atom name
  AtomPosition(double _pos1, double _pos2, double _pos3,
               std::string const &_species)
      : m_species(_species), m_position(_pos1, _pos2, _pos3) {}

  /// \brief Construct with vector position and atom name
  AtomPosition(Eigen::Ref<const Eigen::Vector3d> const &_pos,
               std::string const &_species)
      : m_species(_species), m_position(_pos) {}

  /// Const access of species name
  std::string const &name() const { return m_species; }

  /// \brief Const access of Cartesian position of atom
  Eigen::Vector3d const &cart() const { return m_position; }

  bool time_reversal_active() const {
    for (auto const &_attr : attributes()) {
      if (_attr.second.traits().time_reversal_active()) return true;
    }
    return false;
  }

  std::map<std::string, SpeciesAttribute> const &attributes() const {
    return m_attribute_map;
  }

  void set_attributes(std::map<std::string, SpeciesAttribute> _attr) {
    m_attribute_map = std::move(_attr);
  }

  /// \brief Comparison with tolerance (max allowed distance between LHS and
  /// RHS, in Angstr.)
  bool identical(AtomPosition const &RHS, double _tol) const;

 private:
  /// Atomic species
  std::string m_species;

  /// Cartesian position; origin is centered at site
  Eigen::Vector3d m_position;

  std::map<std::string, SpeciesAttribute> m_attribute_map;
};

bool compare_type(AtomPosition const &A, AtomPosition const &B, double tol);
//****************************************************

/** \defgroup Molecule
 *  \ingroup Crystallography
 *  \brief Relates to Molecule
 *  @{
 */

/// \brief Class representing a Molecule
///
/// - A Molecule is a vector of AtomPosition, with a name
/// - Vacancies are represented as a single atom Molecule, with molecule name ==
/// atom name == "Va"
/// - "make_atom" makes a Molecule with a single atom, with molecule name same
/// as atom name
/// - "make_vacancy" makes a Molecule with a single atom, with molecule name ==
/// atom name == "Va"
///
class Molecule {
 public:
  /// \brief Return an atomic Molecule with specified name
  static Molecule make_atom(std::string const &atom_name);

  /// \brief Return an atomic Molecule with specified name
  static Molecule make_unknown() { return make_atom("UNKNOWN"); }

  /// \brief Return a vacancy Molecule
  static Molecule make_vacancy();

  ///\brief Construct with designated name, a list of atoms, and whether
  ///molecule is chemically divisible
  Molecule(std::string const &_name, std::vector<AtomPosition> _atoms = {},
           bool _divisible = false)
      : m_name(_name), m_atoms(std::move(_atoms)), m_divisible(_divisible) {
    if (m_atoms.empty()) m_atoms.emplace_back(0., 0., 0., m_name);
  }

  /// \brief Number of atoms contained Molecule
  Index size() const { return m_atoms.size(); }

  ///\brief Designated name of Molecule (may be unrelated to constituent
  ///species)
  std::string const &name() const { return m_name; }

  ///\brief Const access of all contained AtomPositions
  std::vector<AtomPosition> const &atoms() const { return m_atoms; }

  ///\brief returns i'th atom position
  AtomPosition const &atom(Index i) const { return m_atoms[i]; }

  ///\brief True if Molecule is atom with no other attributes
  bool is_atomic() const;

  ///\brief True if Molecule represents vacancy
  bool is_vacancy() const;

  ///\brief True if Molecule contains attributes that are affected by time
  ///reversal
  bool time_reversal_active() const {
    for (auto const &_atom : atoms()) {
      if (_atom.time_reversal_active()) return true;
    }
    for (auto const &_attr : attributes()) {
      if (_attr.second.traits().time_reversal_active()) return true;
    }
    return false;
  }

  ///\brief Returns dictionary of all constituent attributes of the Molecule
  /// Does not include attributes associated with individual atoms
  std::map<std::string, SpeciesAttribute> const &attributes() const {
    return m_attribute_map;
  }

  ///\brief Set all constitutent attributes of Molecule
  /// overwrites any existing attributes
  void set_attributes(std::map<std::string, SpeciesAttribute> _attr) {
    m_attribute_map = std::move(_attr);
  }

  ///\brief set all constituent atoms of Molecule
  /// overwrites any existing atoms
  void set_atoms(std::vector<AtomPosition> _atoms) {
    m_atoms = std::move(_atoms);
  }

  /// \brief Check equality of two molecules, within specified tolerance.
  /// Compares atoms, irrespective of order, and attributes (name is not
  /// checked)
  bool identical(Molecule const &RHS, double _tol) const;

  /// \brief Returns true of molecule contains atom of specified name
  bool contains(std::string const &atom_name) const;

  bool is_divisible() const { return m_divisible; }

  bool is_indivisible() const { return !m_divisible; }

 private:
  std::string m_name;
  std::vector<AtomPosition> m_atoms;
  bool m_divisible;

  std::map<std::string, SpeciesAttribute> m_attribute_map;
};

inline bool operator==(Molecule const &A, Molecule const &B) {
  return A.identical(B, TOL);
}

/// \brief A vacancy is any Specie/Molecule with (name == "VA" || name == "va"
/// || name == "Va")
inline bool is_vacancy(const std::string &name) {
  return (name == "VA" || name == "va" || name == "Va");
}

/// \brief Return true if Molecule name matches 'name', including Va checks
inline bool is_molecule_name(const Molecule &mol, std::string name) {
  return mol.name() == name || (mol.is_vacancy() && is_vacancy(name));
}

/** @} */
}  // namespace xtal
}  // namespace CASM

namespace CASM {
namespace sym {
xtal::AtomPosition &apply(const xtal::SymOp &op,
                          xtal::AtomPosition &mutating_atom_pos);
xtal::AtomPosition copy_apply(const xtal::SymOp &op,
                              xtal::AtomPosition atom_pos);

template <typename ExternSymOp>
xtal::AtomPosition copy_apply(const ExternSymOp &op,
                              xtal::AtomPosition atom_pos) {
  return sym::copy_apply(adapter::Adapter<xtal::SymOp, ExternSymOp>()(op),
                         atom_pos);
}

xtal::Molecule &apply(const xtal::SymOp &op, xtal::Molecule &mutating_mol);
xtal::Molecule copy_apply(const xtal::SymOp &op, xtal::Molecule mol);

template <typename ExternSymOp>
xtal::Molecule copy_apply(const ExternSymOp &op, xtal::Molecule mol) {
  return sym::copy_apply(adapter::Adapter<xtal::SymOp, ExternSymOp>()(op), mol);
}
}  // namespace sym
}  // namespace CASM

#endif
