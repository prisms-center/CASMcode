#ifndef CASM_StrucMapCalculatorInterface
#define CASM_StrucMapCalculatorInterface

#include <iostream>
#include <unordered_set>
#include <vector>

#include "casm/crystallography/Adapter.hh"
#include "casm/crystallography/SimpleStructure.hh"
#include "casm/crystallography/SymType.hh"
#include "casm/global/definitions.hh"
#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/misc/CASM_math.hh"

namespace CASM {
namespace xtal {
class SimpleStructure;

// In this file:
struct MappingNode;
class StrucMapCalculatorInterface;

namespace StrucMapping {

// Denotes (name, number composition) of any species whose number-composition is
// fixed for the parent primitive cell
using FixedSpecies = std::map<std::string, Index>;

// List of species allowed at each site of primitive
using AllowedSpecies = std::vector<std::vector<std::string>>;

}  // namespace StrucMapping

class StrucMapCalculatorInterface {
 public:
  StrucMapCalculatorInterface(
      SimpleStructure _parent,
      SymOpVector const &_factor_group = {SymOp::identity()},
      SimpleStructure::SpeciesMode _species_mode =
          SimpleStructure::SpeciesMode::ATOM,
      StrucMapping::AllowedSpecies allowed_species = {})
      : m_parent(std::move(_parent)), m_species_mode(_species_mode) {
    this->_set_sym_info(_factor_group);

    if (allowed_species.empty()) {
      auto const &p_info(this->struc_info(parent()));
      allowed_species.resize(p_info.size());
      for (Index i = 0; i < p_info.size(); ++i) {
        allowed_species[i].push_back(p_info.names[i]);
      }
    }
    _set_allowed_species(allowed_species);
  }

  SimpleStructure::Info const &struc_info(SimpleStructure const &_struc) const {
    return _struc.info(m_species_mode);
  }

  SimpleStructure::Info &struc_info(SimpleStructure &_struc) const {
    return _struc.info(m_species_mode);
  }

  SimpleStructure const &parent() const { return m_parent; }

  ///\brief Crystallographic tolerance, for now just return CASM::TOL
  double xtal_tol() const { return TOL; }

  ///\brief List of point group operations that map parent onto itself
  ///(neglecting internal translation)
  SymOpVector const &point_group() const { return m_point_group; }

  ///\brief List of internal translations that map parent onto itself
  std::vector<Eigen::Vector3d> const &internal_translations() const {
    return m_internal_translations;
  }

  std::unordered_set<Index> const &va_allowed() const { return m_va_allowed; }

  StrucMapping::FixedSpecies const &fixed_species() const {
    return m_fixed_species;
  }

  /// \brief Return maximum possible number of vacancies in underlying primitive
  /// structure
  Index max_n_va() const { return va_allowed().size(); }

  virtual ~StrucMapCalculatorInterface() {}

  /// \brief construct list of prospective mapping translations
  virtual std::vector<Eigen::Vector3d> translations(
      MappingNode const &_node, SimpleStructure const &child_struc) const = 0;

  /// \brief Creates copy of _child_struc by applying isometry, lattice
  /// transformation, translation, and site permutation of _node
  virtual SimpleStructure resolve_setting(
      MappingNode const &_node, SimpleStructure const &_child_struc) const = 0;

  /// \brief Calculates final mapping score and sets _node.is_valid
  virtual void finalize(MappingNode &_node,
                        SimpleStructure const &child_struc) const = 0;

  virtual bool populate_cost_mat(MappingNode &_node,
                                 SimpleStructure const &child_struc) const = 0;

  /// \brief Make an exact copy of the calculator (including any initialized
  /// members)
  std::unique_ptr<StrucMapCalculatorInterface> clone() const {
    return std::unique_ptr<StrucMapCalculatorInterface>(this->_clone());
  }

  /// \brief Make an exact copy of the calculator (including any initialized
  /// members)
  std::unique_ptr<StrucMapCalculatorInterface> quasi_clone(
      SimpleStructure _parent,
      SymOpVector const &_factor_group = {SymOp::identity()},
      SimpleStructure::SpeciesMode _species_mode =
          SimpleStructure::SpeciesMode::ATOM,
      StrucMapping::AllowedSpecies _allowed_species = {}) const {
    return std::unique_ptr<StrucMapCalculatorInterface>(
        this->_quasi_clone(std::move(_parent), _factor_group, _species_mode,
                           std::move(_allowed_species)));
  }

  template <typename ExternSymOpVector>
  std::unique_ptr<StrucMapCalculatorInterface> quasi_clone(
      SimpleStructure _parent,
      ExternSymOpVector const &_factor_group = {SymOp::identity()},
      SimpleStructure::SpeciesMode _species_mode =
          SimpleStructure::SpeciesMode::ATOM,
      StrucMapping::AllowedSpecies _allowed_species = {}) const {
    return this->quasi_clone(
        _parent,
        adapter::Adapter<SymOpVector, ExternSymOpVector>()(_factor_group),
        _species_mode, _allowed_species);
  }

 protected:
  void _set_allowed_species(StrucMapping::AllowedSpecies allowed_species) {
    m_allowed_species = std::move(allowed_species);

    // Analyze allowed_species:
    for (Index i = 0; i < m_allowed_species.size(); ++i) {
      for (std::string &sp : m_allowed_species[i]) {
        if (sp == "Va" || sp == ("VA") || sp == "va") {
          // Is vacancy
          sp = "Va";
          m_va_allowed.insert(i);
        } else if (m_allowed_species[i].size() > 1) {
          // Not vacancy, but variable species
          m_fixed_species[sp] = 0;
        } else {
          // potentially fixed species
          auto it = m_fixed_species.find(sp);
          if (it == m_fixed_species.end())
            m_fixed_species[sp] = 1;
          else if ((it->second) > 0)
            ++(it->second);
        }
        if (!m_max_n_species.count(sp)) m_max_n_species[sp] = 0;
        ++m_max_n_species[sp];
      }
    }

    for (auto it = m_fixed_species.begin(); it != m_fixed_species.end();) {
      auto curr = it;
      ++it;
      if ((curr->second) == 0) m_fixed_species.erase(curr);
    }
  }

  ///\brief Sets point_group and internal_translations by breaking factor group
  ///into pure translations and rotations/rotoreflections
  /// _factor_group should be sorted in order of decreasing character
  void _set_sym_info(SymOpVector const &_factor_group) {
    m_point_group.clear();
    m_internal_translations.clear();

    // Internal translations are any translation associated with identity
    m_point_group.push_back(SymOp::identity());

    for (SymOp const &op : _factor_group) {
      if (get_matrix(op).isIdentity(TOL) && !get_time_reversal(op)) {
        m_internal_translations.push_back(get_translation(op));
      }
      if (!almost_equal(get_matrix(op), get_matrix(m_point_group.back()),
                        TOL) ||
          get_time_reversal(op) != get_time_reversal(m_point_group.back())) {
        m_point_group.push_back(op);
        m_point_group.back().translation.setZero();
      }
    }
  }

  StrucMapping::AllowedSpecies const &_allowed_species() const {
    return m_allowed_species;
  }

  StrucMapping::FixedSpecies const &_fixed_species() const {
    return m_fixed_species;
  }

  bool _sublat_allows_va(Index b) const { return m_va_allowed.count(b); }

  /// \brief maximum allowed number of each species
  std::map<std::string, Index> const &_max_n_species() const {
    return m_max_n_species;
  }

 private:
  SimpleStructure m_parent;

  // Point group of parent crystal, must contain at least Identity
  SymOpVector m_point_group;

  // internal translations that map parent crystal onto itself, must contain at
  // least (0,0,0)
  std::vector<Eigen::Vector3d> m_internal_translations;

  SimpleStructure::SpeciesMode m_species_mode;

  StrucMapping::AllowedSpecies m_allowed_species;

  StrucMapping::FixedSpecies m_fixed_species;

  /// \brief maximum allowed number of each species
  std::map<std::string, Index> m_max_n_species;

  std::unordered_set<Index> m_va_allowed;

  /// \brief Make an exact copy of the calculator (including any initialized
  /// members)
  virtual StrucMapCalculatorInterface *_clone() const = 0;

  /// \brief Make an exact copy of the calculator (including any initialized
  /// members)
  virtual StrucMapCalculatorInterface *_quasi_clone(
      SimpleStructure _parent,
      SymOpVector const &_factor_group = {SymOp::identity()},
      SimpleStructure::SpeciesMode _species_mode =
          SimpleStructure::SpeciesMode::ATOM,
      StrucMapping::AllowedSpecies _allowed_species = {}) const = 0;
};
}  // namespace xtal
}  // namespace CASM

#endif
