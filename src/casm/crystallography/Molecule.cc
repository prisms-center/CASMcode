#include "casm/crystallography/Molecule.hh"

#include <type_traits>
#include <vector>

#include "casm/crystallography/SpeciesAttribute.hh"
#include "casm/misc/CASM_Eigen_math.hh"

namespace CASM {
namespace xtal {

bool AtomPosition::identical(AtomPosition const &RHS, double _tol) const {
  return almost_equal(cart(), RHS.cart(), _tol) &&
         compare_type(*this, RHS, _tol);
}

//****************************************************

bool compare_type(AtomPosition const &A, AtomPosition const &B, double tol) {
  // compare number of attributes
  if (A.attributes().size() != B.attributes().size()) return false;
  if (A.name() != B.name()) return false;
  // compare attributes
  auto it_A(A.attributes().cbegin()), end_it(A.attributes().cend());
  for (; it_A != end_it; ++it_A) {
    auto it_B = B.attributes().find(it_A->first);
    if (it_B == B.attributes().cend() ||
        !(it_A->second).identical(it_B->second, tol))
      return false;
  }
  return true;
}

bool Molecule::is_atomic() const {
  if (size() != 1) return false;
  if (!attributes().empty()) return false;
  for (AtomPosition const &atom : atoms()) {
    if (atom.cart().norm() > TOL) return false;
    if (!atom.attributes().empty()) return false;
  }
  return true;
}

bool Molecule::is_vacancy() const {
  return ::CASM::xtal::is_vacancy(m_atoms[0].name());
}

bool Molecule::identical(Molecule const &RHS, double _tol) const {
  // compare number of attributes
  if (m_attribute_map.size() != RHS.m_attribute_map.size()) return false;

  // compare number of atoms
  if (size() != RHS.size()) return false;

  // compare atoms, irrespective of order
  for (Index i = 0; i < RHS.size(); i++) {
    Index j = 0;
    for (j = 0; j < size(); j++) {
      if (atom(i).identical(RHS.atom(j), _tol)) break;
    }
    if (j == size()) return false;
  }

  // compare attributes
  auto it(m_attribute_map.cbegin()), end_it(m_attribute_map.cend());
  for (; it != end_it; ++it) {
    auto it_RHS = RHS.m_attribute_map.find(it->first);
    if (it_RHS == RHS.m_attribute_map.cend() ||
        !(it->second).identical(it_RHS->second, _tol))
      return false;
  }

  return true;
}

bool Molecule::contains(std::string const &_name) const {
  for (Index i = 0; i < size(); i++)
    if (atom(i).name() == _name) return true;
  return false;
}

}  // namespace xtal
}  // namespace CASM

namespace CASM {
namespace sym {
xtal::AtomPosition &apply(const xtal::SymOp &op,
                          xtal::AtomPosition &mutating_atom_pos) {
  xtal::AtomPosition transformed_atom_pos = copy_apply(op, mutating_atom_pos);
  std::swap(mutating_atom_pos, transformed_atom_pos);
  return mutating_atom_pos;
}

xtal::AtomPosition copy_apply(const xtal::SymOp &op,
                              xtal::AtomPosition atom_pos) {
  Eigen::Vector3d transformed_position = get_matrix(op) * atom_pos.cart();
  xtal::AtomPosition transformed_atom_pos(transformed_position,
                                          atom_pos.name());
  std::map<std::string, xtal::SpeciesAttribute> transformed_attribute_map;
  for (const auto &name_attr_pr : atom_pos.attributes()) {
    transformed_attribute_map.emplace(name_attr_pr.first,
                                      sym::copy_apply(op, name_attr_pr.second));
  }
  transformed_atom_pos.set_attributes(transformed_attribute_map);

  return transformed_atom_pos;
}

xtal::Molecule &apply(const xtal::SymOp &op, xtal::Molecule &mutating_mol) {
  std::vector<xtal::AtomPosition> transformed_atoms;
  for (const xtal::AtomPosition &atom_pos : mutating_mol.atoms()) {
    transformed_atoms.emplace_back(copy_apply(op, atom_pos));
  }
  mutating_mol.set_atoms(transformed_atoms);

  std::map<std::string, xtal::SpeciesAttribute> transformed_attribute_map;
  for (const auto &name_attr_pr : mutating_mol.attributes()) {
    transformed_attribute_map.emplace(name_attr_pr.first,
                                      copy_apply(op, name_attr_pr.second));
  }
  mutating_mol.set_attributes(transformed_attribute_map);

  return mutating_mol;
}

xtal::Molecule copy_apply(const xtal::SymOp &op, xtal::Molecule mol) {
  xtal::Molecule transformed_mol = apply(op, mol);
  return transformed_mol;
}
}  // namespace sym
}  // namespace CASM
