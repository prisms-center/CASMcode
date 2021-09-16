#ifndef CASM_config_Group_impl
#define CASM_config_Group_impl

#include <numeric>

#include "casm/configuration/Group.hh"

namespace CASM {
namespace config {

namespace Group_impl {

std::vector<Index> _identity_indices(Index n) {
  std::vector<Index> result(n);
  std::iota(result.begin(), result.end(), 0);
  return result;
}

std::vector<SymOp> _make_subgroup_elements(
    std::shared_ptr<Group const> const &_head_group,
    std::set<Index> const &_head_group_index) {
  std::vector<SymOp> result;
  for (Index index : _head_group_index) {
    result.push_back(_head_group->element[index]);
  }
  return result;
}

std::vector<SymOp> _make_subgroup_multiplication_table(
    std::shared_ptr<Group const> const &_head_group,
    std::set<Index> const &_head_group_index) {
  MultiplicationTable result(_head_group_index.size());
  MultiplicationTable const &head_group_table =
      _head_group.multiplication_table;
  Index N = head_group_table.size();

  for (Index index : _head_group_index) {
    if (index >= N) {
      throw std::runtime_error(
          "Error in Group constructor: head group index >= head group "
          "multiplication table size");
    }
  }

  Index row = 0;
  for (Index i = 0; i < N; ++i) {
    if (head_group_table[i].size() != N) {
      throw std::runtime_error(
          "Error in Group constructor: head group multiplication table is not "
          "square");
    }
    if (!_head_group_index.count(i)) {
      continue;
    }

    for (Index j = 0; j < N; ++j) {
      if (!_head_group_index.count(j)) {
        continue;
      }
      auto it = _head_group_index.find(head_group_table[i][j]);
      if (it == _head_group_index.end()) {
        throw std::runtime_error(
            "Error in Group constructor: subgroup is not closed according to "
            "the head group multiplication table.");
      }
      Index subgroup_entry = std::distance(_head_group_index.begin(), it);
      result[row].push_back(subgroup_entry);
    }

    ++row;
  }
  for (Index index : _head_group_index) {
    result.push_back(_head_group->element[index]);
  }
  return result;
}

/// \brief Collect indices of the inverse elements in a group using the
/// multiplication table
///
/// Notes:
/// - Requires that identity element corresponds to index 0
std::vector<Index> make_inverse_index(
    MultiplicationTable const &multiplication_table) {
  std::vector<Index> index_inverse;

  Index N = multiplication_table.size();
  for (Index i = 0; i < N; ++i) {
    if (multiplication_table[i].size() != N) {
      throw std::runtime_error(
          "Error in Group constructor: multiplication table is not square");
    }
  }

  for (auto const &row : multiplication_table) {
    auto begin = std::begin(row);
    auto end = std::end(row);
    auto it = std::find(begin, end, 0);
    if (it == end) {
      throw std::runtime_error(
          "Error in make_inverse_index: no inverse element");
    }
    index_inverse.push_back(std::distance(begin, it));
  }
  return index_inverse;
}

}  // namespace Group_impl

/// \brief Construct a head group
///
/// \params _element Group elements, expected to be closed and sorted as desired
///
Group::Group(std::vector<SymOp> const &_element,
             MultiplicationTable const &_multiplication_table)
    : head_group(nullptr),
      element(_element),
      head_group_index(Group_impl::_identity_indices(element.size())),
      multiplication_table(_multiplication_table),
      inverse_index(Group_impl::_make_inverse_index(multiplication_table)) {}

/// \brief Construct a sub group
Group::Group(std::shared_ptr<Group const> const &_head_group,
             std::set<Index> const &_head_group_index)
    : head_group(_head_group),
      element(
          Group_impl::_make_subgroup_elements(_head_group, _head_group_index)),
      head_group_index(_head_group_index),
      multiplication_table(Group_impl::_make_subgroup_multiplication_table(
          _head_group, _head_group_index)),
      inverse_index(Group_impl::_make_inverse_index(multiplication_table)) {}

}  // namespace config
}  // namespace CASM

#endif
