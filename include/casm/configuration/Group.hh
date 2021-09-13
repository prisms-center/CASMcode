#ifndef CASM_config_Group
#define CASM_config_Group

#include <memory>

#include "casm/configuration/definitions.hh"
#include "casm/crystallography/SymTypes.hh"

namespace CASM {
namespace config {

/// \brief Holds group elements and multiplication table
struct Group {
  /// \brief Construct a head group
  Group(std::vector<SymOp> const &_element,
        MultiplicationTable const &_multiplication_table);

  /// \brief Construct a sub group
  Group(std::shared_ptr<Group const> const &_head_group,
        std::set<Index> const &_head_group_index);

  /// \brief If this is a subgroup, indicates the head group; if this is a head
  /// group, then this is empty
  std::shared_ptr<Group const> const head_group;

  /// \brief Specifies the group elements
  std::vector<Element> const element;

  /// \brief Specifies the head group index for each element (guaranteed sorted)
  ///
  /// If this is the head group, then:
  ///
  ///     this->head_group_index = [0, 1, 2, ...]
  ///
  /// If this is a sub group, then:
  ///
  ///     this->element[i] == head_group->element[this->head_group_index[i]]
  ///
  std::vector<Index> const head_group_index;

  /// \brief Specifies the multiplication table for the elements
  ///
  /// element[k] == element[i] * element[j],
  /// where k = multiplication_table[i][j]
  MultiplicationTable const multiplication_table;

  /// \brief Specifies the index of the inverse element
  ///
  ///     I == element[i] * element[inverse_index[i]]
  ///       == element[inverse_index[i]] * element[i]
  std::vector<Index> const inverse_index;
};

}  // namespace config
}  // namespace CASM

#endif
