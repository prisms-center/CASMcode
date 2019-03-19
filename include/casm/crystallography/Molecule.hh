#ifndef MOLECULE_HH
#define MOLECULE_HH

#include <iostream>
#include <array>
#include <vector>

#include "casm/misc/Comparisons.hh"
#include "casm/crystallography/SpeciesAttribute.hh"

namespace CASM {

  /** \defgroup Molecule
   *  \ingroup Crystallography
   *  \brief Relates to Molecule
   *  @{
   */

  class SymOp;
  class Molecule;


  //****************************************************

  /// \brief An atomic species associated with a position in space
  class AtomPosition {
  public:

    /// Typedef for selective dynamics array
    typedef std::array<bool, 3> sd_type;


    /// \brief Construct with x,y,z position coordinates and atom name
    AtomPosition(double _pos1,
                 double _pos2,
                 double _pos3,
                 std::string const &_species,
    sd_type const &_sd_flag = sd_type{{false, false, false}}) :
      m_species(_species),
      m_position(_pos1, _pos2, _pos3),
      m_sd_flag(_sd_flag) { }

    /// \brief Construct with vector position and atom name
    AtomPosition(Eigen::Ref<const Eigen::Vector3d> const &_pos,
                 std::string const &_species,
    sd_type const &_sd_flag = sd_type{{false, false, false}}) :
      m_species(_species),
      m_position(_pos),
      m_sd_flag(_sd_flag) { }

    /// Const access of species name
    std::string const &name() const {
      return m_species;
    }

    /// \brief Const access of Cartesian position of atom
    Eigen::Vector3d const &cart() const {
      return m_position;
    }

    /// \brief Const access of selective dynamics flags
    sd_type const &sd_flag() const {
      return m_sd_flag;
    }

    bool time_reversal_active() const {
      for(auto const &_attr : attributes()) {
        if(_attr.second.traits().time_reversal_active())
          return true;
      }
      return false;
    }

    std::map<std::string, SpeciesAttribute> const &attributes() const {
      return m_attribute_map;
    }

    void set_attributes(std::map<std::string, SpeciesAttribute> _attr) {
      m_attribute_map = std::move(_attr);
    }

    /// \brief Comparison with tolerance (max allowed distance between LHS and RHS, in Angstr.)
    bool identical(AtomPosition const &RHS, double _tol)const;

    /// \brief Print AtomPosition after applying affine transformation cart2frac*cart()+trans
    void print(std::ostream &stream,
               Eigen::Ref<const Eigen::Vector3d> const &trans,
               Eigen::Ref<const Eigen::Matrix3d> const &cart2frac,
               int spaces,
               bool print_sd_flags = false) const;

    /// \brief Apply symmetry (translation is ignored)
    AtomPosition &apply_sym(const SymOp &op);

  private:
    /// Atomic species
    std::string m_species;

    /// Cartesian position; origin is centered at site
    Eigen::Vector3d m_position;

    std::map<std::string, SpeciesAttribute> m_attribute_map;

    /// selective dynamics flags
    sd_type m_sd_flag;
  };

  //****************************************************

  /** \defgroup Molecule
   *  \ingroup Crystallography
   *  \brief Relates to Molecule
   *  @{
   */

  /// \brief Class representing a Molecule
  ///
  /// - A Molecule is a vector of AtomPosition, with a name
  /// - Vacancies are represented as a single atom Molecule, with molecule name == atom name == "Va"
  /// - "make_atom" makes a Molecule with a single atom, with molecule name same as atom name
  /// - "make_vacancy" makes a Molecule with a single atom, with molecule name == atom name == "Va"
  ///
  class Molecule {
  public:
    /// \brief Return an atomic Molecule with specified name
    static Molecule make_atom(std::string const &atom_name,
                              AtomPosition::sd_type const &_sd_flags);

    /// \brief Return an atomic Molecule with specified name
    static Molecule make_atom(std::string const &atom_name) {
      return make_atom(atom_name, AtomPosition::sd_type{{false, false, false}});
    }

    /// \brief Return an atomic Molecule with specified name
    static Molecule make_unknown() {
      return make_atom("UNKNOWN");
    }

    /// \brief Return a vacancy Molecule
    static Molecule make_vacancy();

    /*
    Molecule(std::string const &_name,
             std::initializer_list<AtomPosition> const &_atoms = {},
             bool _divisible = false) :
      m_name(_name),
      m_atoms(_atoms),
      m_divisible(_divisible) {}
    */
    Molecule(std::string const &_name,
             std::vector<AtomPosition> _atoms = {},
             bool _divisible = false) :
      m_name(_name),
      m_atoms(std::move(_atoms)),
      m_divisible(_divisible) {}

    Index size() const {
      return m_atoms.size();
    }

    std::string const &name() const {
      return m_name;
    }

    std::vector<AtomPosition> const &atoms() const {
      return m_atoms;
    }

    AtomPosition const &atom(Index i) const {
      return m_atoms[i];
    }

    bool is_vacancy() const;

    bool time_reversal_active() const {
      for(auto const &_atom : atoms()) {
        if(_atom.time_reversal_active())
          return true;
      }
      for(auto const &_attr : attributes()) {
        if(_attr.second.traits().time_reversal_active())
          return true;
      }
      return false;
    }


    std::map<std::string, SpeciesAttribute> const &attributes() const {
      return m_attribute_map;
    }

    void set_attributes(std::map<std::string, SpeciesAttribute> _attr) {
      m_attribute_map = std::move(_attr);
    }

    void set_atoms(std::vector<AtomPosition> _atoms) {
      m_atoms = std::move(_atoms);
    }

    Molecule &apply_sym(SymOp const &op);

    /// \brief Check equality of two molecules, within specified tolerance.
    /// Compares atoms, irrespective of order, and attributes (name is not checked)
    bool identical(Molecule const &RHS, double _tol) const;

    /// \brief Returns true of molecule contains atom of specified name
    bool contains(std::string const &atom_name) const;

    void read(std::istream &stream);

    void print(std::ostream &stream,
               Eigen::Ref<const Eigen::Vector3d> const &trans,
               Eigen::Ref<const Eigen::Matrix3d> const &cart2frac,
               int spaces,
               char delim,
               bool print_sd_flags  = false) const;

    bool is_divisible() const {
      return m_divisible;
    }

    bool is_indivisible() const {
      return !m_divisible;
    }

  private:
    std::string m_name;
    std::vector<AtomPosition> m_atoms;
    bool m_divisible;

    std::map<std::string, SpeciesAttribute> m_attribute_map;
  };

  inline
  bool operator==(Molecule const &A, Molecule const &B) {
    return A.identical(B, TOL);
  }

  /// \brief A vacancy is any Specie/Molecule with (name == "VA" || name == "va" || name == "Va")
  inline bool is_vacancy(const std::string &name) {
    return (name == "VA" || name == "va" || name == "Va");
  }

  /// \brief Return true if Molecule name matches 'name', including Va checks
  inline
  bool is_molecule_name(const Molecule &mol, std::string name) {
    return mol.name() == name || (mol.is_vacancy() && is_vacancy(name));
  }

  /** @} */
}

#endif
