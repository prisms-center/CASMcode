#ifndef MOLECULE_HH
#define MOLECULE_HH

#include <iostream>
#include <array>
#include <vector>

#include "casm/misc/Comparisons.hh"
#include "casm/crystallography/MoleculeAttribute.hh"

namespace CASM {

  /** \defgroup Molecule
   *  \ingroup Crystallography
   *  \brief Relates to Molecule
   *  @{
   */

  class SymOp;
  class Molecule;
  template<typename T> struct jsonConstructor;

  //****************************************************
  ///\brief Lightweight container for intrinsic atom properties.

  /// For now, AtomSpecie only contains the name, but in future other properties
  /// may be needed (mass, atomic number, etc).
  // Additional fields should only be added if absolutely necessary!
  class AtomSpecie : public Comparisons<CRTPBase<AtomSpecie>> {
  public:

    /// \brief Constructor
    AtomSpecie(std::string const &_name) :
      m_name(_name) {}

    /// \brief Return name of species
    std::string const &name() const {
      return m_name;
    }

    /// \brief Equality comparison for two AtomSpecies
    bool operator==(AtomSpecie const &RHS) const {
      return name() == RHS.name();
    }

    bool operator<(AtomSpecie const &RHS) const {
      return name() < RHS.name();
    }

  private:
    std::string m_name;
  };

  void from_json(AtomSpecie &_specie, jsonParser const &json);

  jsonParser &to_json(AtomSpecie const &_specie, jsonParser &json);

  //****************************************************

  /// \brief An atomic specie associated with a position in space
  /// Also tracks selective dynamics flags for atom
  class AtomPosition {
  public:

    /// Typedef for selective dynamics array
    typedef std::array<bool, 3> sd_type;

    /// \brief Construct with x,y,z position coordinates and AtomSpecie
    template<typename AtomSpecieConvertible>
    AtomPosition(double _pos1,
                 double _pos2,
                 double _pos3,
                 AtomSpecieConvertible _specie,
                 sd_type const &_sd_flag = sd_type{false, false, false}) :
      m_position(_pos1, _pos2, _pos3),
      m_specie(_specie),
      m_sd_flag(_sd_flag) { }

    /// \brief Construct with vector position and AtomSpecie
    template<typename AtomSpecieConvertible>
    AtomPosition(Eigen::Ref<const Eigen::Vector3d> const &_pos,
                 AtomSpecieConvertible _specie,
                 sd_type const &_sd_flag = sd_type{false, false, false}) :
      m_position(_pos),
      m_specie(_specie),
      m_sd_flag(_sd_flag) { }

    /// \brief Const access of specie name
    std::string const &name() const {
      return m_specie.name();
    }

    /// \brief Const access of atomic specie
    AtomSpecie const &specie() const {
      return m_specie;
    }

    /// \brief Const access of Cartesian position of atom
    Eigen::Vector3d const &cart() const {
      return m_position;
    }

    /// \brief Const access of selective dynamics flags
    sd_type const &sd_flag() const {
      return m_sd_flag;
    }

    /// \brief Print AtomPosition after applying affine transformation cart2frac*cart()+trans
    void print(std::ostream &stream,
               Eigen::Ref<const Eigen::Vector3d> const &trans,
               Eigen::Ref<const Eigen::Matrix3d> const &cart2frac,
               int spaces,
               bool print_sd_flags = false) const;

    /// \brief Apply symmetry (translation is ignored)
    AtomPosition &apply_sym(const SymOp &op);

  private:
    /// Atomic specie
    AtomSpecie m_specie;

    /// Cartesian position; origin is centered at site
    Eigen::Vector3d m_position;

    /// selective dynamics flags
    sd_type m_sd_flag;
  };

  /// \brief Comparison with tolerance (max allowed distance between LHS and RHS, in Angstr.)
  bool identical(AtomPosition const &LHS, AtomPosition const &RHS, double _tol);

  /// \brief Print AtomPosition to json after applying affine transformation cart2frac*cart()+trans
  jsonParser &to_json(const AtomPosition &apos, jsonParser &json, Eigen::Ref<const Eigen::Matrix3d> const &cart2frac);

  /// \brief Read AtomPosition from json and then apply affine transformation cart2frac*cart()+trans
  void from_json(AtomPosition &apos, const jsonParser &json, Eigen::Ref<const Eigen::Matrix3d> const &frac2cart);

  template<>
  struct jsonConstructor<AtomPosition> {

    /// \brief Read from json [b, i, j, k], using 'unit' for AtomPosition::unit()
    static AtomPosition from_json(const jsonParser &json, Eigen::Matrix3d const &f2c_mat);
  };

  //****************************************************

  /** \defgroup Molecule
   *  \ingroup Crystallography
   *  \brief Relates to Molecule
   *  @{
   */

  class Molecule {
  public:
    /// \brief Return an atomic Molecule with specified name
    static Molecule make_atom(std::string const &atom_name,
                              AtomPosition::sd_type const &_sd_flags);

    /// \brief Return an atomic Molecule with specified name
    static Molecule make_atom(std::string const &atom_name) {
      return make_atom(atom_name, AtomPosition::sd_type{false, false, false});
    }

    /// \brief Return a vacancy Molecule
    static Molecule make_vacancy();

    Molecule(std::string const &_name,
             std::initializer_list<AtomPosition> const &_atoms,
             bool _divisible = false) :
      m_name(_name),
      m_atoms(_atoms),
      m_divisible(_divisible) {}

    Molecule(std::string const &_name,
             std::vector<AtomPosition> const &_atoms,
             bool _divisible = false) :
      m_name(_name),
      m_atoms(_atoms),
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

    jsonParser &to_json(jsonParser &json, Eigen::Matrix3d const &c2f_mat) const;

    void from_json(const jsonParser &json, Eigen::Matrix3d const &f2c_mat);

    bool is_divisible() const {
      return m_divisible;
    }

    bool is_indivisible() const {
      return !m_divisible;
    }

  private:
    std::string m_name;
    std::vector<AtomPosition> m_atoms;
    std::map<std::string, MoleculeAttribute> m_attribute_map;
    bool m_divisible;
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

  jsonParser &to_json(const Molecule &mol, jsonParser &json, Eigen::Ref<const Eigen::Matrix3d> const &c2f_mat);

  void from_json(Molecule &mol, const jsonParser &json, Eigen::Matrix3d const &f2c_mat);

  template<>
  struct jsonConstructor<Molecule> {
    static Molecule from_json(const jsonParser &json, Eigen::Ref<const Eigen::Matrix3d> const &f2c_mat);
  };
  /** @} */
}

#endif
