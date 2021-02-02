#ifndef CASM_SimpleStrucMapCalculator
#define CASM_SimpleStrucMapCalculator

#include "casm/crystallography/Adapter.hh"
#include "casm/crystallography/StrucMapCalculatorInterface.hh"

namespace CASM {
namespace xtal {
class SimpleStructure;
struct SymOp;
typedef std::vector<SymOp> SymOpVector;

// In this file:
struct MappingNode;
class StrucMapCalculatorInterface;
class SimpleStrucMapCalculator;

class SimpleStrucMapCalculator : public StrucMapCalculatorInterface {
 public:
  SimpleStrucMapCalculator(
      SimpleStructure _parent,
      SymOpVector const &_factor_group = {SymOp::identity()},
      SimpleStructure::SpeciesMode species_mode =
          SimpleStructure::SpeciesMode::ATOM,
      StrucMapping::AllowedSpecies allowed_species = {})
      : StrucMapCalculatorInterface(std::move(_parent), _factor_group,
                                    species_mode, std::move(allowed_species)) {}

  template <typename ExternSymOpVector>
  SimpleStrucMapCalculator(
      SimpleStructure _parent,
      ExternSymOpVector const &_factor_group = {SymOp::identity()},
      SimpleStructure::SpeciesMode species_mode =
          SimpleStructure::SpeciesMode::ATOM,
      StrucMapping::AllowedSpecies allowed_species = {})
      : SimpleStrucMapCalculator(
            _parent,
            adapter::Adapter<SymOpVector, ExternSymOpVector>()(_factor_group),
            species_mode, allowed_species) {}

  virtual ~SimpleStrucMapCalculator() {}

  std::vector<Eigen::Vector3d> translations(
      MappingNode const &_node,
      SimpleStructure const &child_struc) const override;

  /// \brief Creates copy of _child_struc by applying isometry, lattice
  /// transformation, translation, and site permutation of _node Result has all
  /// sites within the unit cell. After setting resolution, the lattice and
  /// sites of _child_struc match the setting of the parent structure onto which
  /// it has been mapped (as defined by '_node')
  virtual SimpleStructure resolve_setting(
      MappingNode const &_node,
      SimpleStructure const &_child_struc) const override;

  void finalize(MappingNode &_node, SimpleStructure const &child_struc,
                bool const &symmetrize_atomic_cost = false) const override;

  bool populate_cost_mat(MappingNode &_node,
                         SimpleStructure const &child_struc) const override;

  void populate_displacement(MappingNode &_node,
                             SimpleStructure const &child_struc) const;

 private:
  /// \brief Make an exact copy of the calculator (including any initialized
  /// members)
  virtual StrucMapCalculatorInterface *_clone() const override {
    return new SimpleStrucMapCalculator(*this);
  }

  /// \brief Make an exact copy of the calculator (including any initialized
  /// members)
  virtual StrucMapCalculatorInterface *_quasi_clone(
      SimpleStructure _parent,
      SymOpVector const &_factor_group = {SymOp::identity()},
      SimpleStructure::SpeciesMode _species_mode =
          SimpleStructure::SpeciesMode::ATOM,
      StrucMapping::AllowedSpecies _allowed_species = {}) const override {
    return new SimpleStrucMapCalculator(std::move(_parent), _factor_group,
                                        _species_mode,
                                        std::move(_allowed_species));
  }

  /// \brief Initializes child_struc.mol_info based on child_struc.atom_info and
  /// _node. Default behavior simply copies atom_info to mol_info
  virtual bool _assign_molecules(MappingNode &_node,
                                 SimpleStructure const &child_struc) const;
};
}  // namespace xtal
}  // namespace CASM
#endif
