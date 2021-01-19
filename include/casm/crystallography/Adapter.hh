#ifndef CRYSTALLOGRAPHYADAPTER_HH
#define CRYSTALLOGRAPHYADAPTER_HH

#include "casm/crystallography/SymType.hh"

namespace CASM {

namespace adapter {
template <typename ToType, typename FromType>
struct Adapter;

/// Convertes any symmetry type to SymOpType, as defined by the crystallography
/// module. Works for any symmetry type that has the get_matrix,
/// get_translation, and get_time_reversal accessors defined.
///
/// If there is a different symmetry type you would like to adapt for the
/// crystallography module, simply declare and define accessors get_matrix,
/// get_translation, and get_time_reversal
template <typename FromType>
struct Adapter<xtal::SymOp, FromType> {
  xtal::SymOp operator()(const FromType &adaptable) {
    return xtal::SymOp(get_matrix(adaptable), get_translation(adaptable),
                       get_time_reversal(adaptable));
  }
};

/// Converts any container of any symmetry type with begin() and end() defined
/// into SymGroupType, as defined by the crystallography module.
template <typename FromType>
struct Adapter<xtal::SymOpVector, FromType> {
  template <typename FromTypeIt>
  xtal::SymOpVector operator()(FromTypeIt begin, FromTypeIt end) {
    xtal::SymOpVector casted_group;
    Adapter<xtal::SymOp, typename FromTypeIt::value_type> to_symop_type;
    for (auto it = begin; it != end; ++it) {
      casted_group.emplace_back(to_symop_type(*it));
    }
    return casted_group;
  }

  xtal::SymOpVector operator()(const FromType &adaptable) {
    return this->operator()(adaptable.begin(), adaptable.end());
  }
};
}  // namespace adapter
}  // namespace CASM

#endif
