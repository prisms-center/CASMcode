#ifndef SYMMETRYTOCRYSTALLOGRAPHY_HH
#define SYMMETRYTOCRYSTALLOGRAPHY_HH

#include "casm/adapters/AdapterInterface.hh"
#include "casm/crystallography/SymType.hh"

namespace CASM {

  namespace adapter {
    /// Convertes any symmetry type to SymOpType, as defined by the crystallography module.
    /// Works for any symmetry type that has the get_matrix, get_translation, and
    /// get_time_reversal accessors defined.
    ///
    /// If there is a different symmetry type you would like to adapt for the crystallography
    /// module, simply declare and define accessors get_matrix, get_translation, and
    /// get_time_reversal
    template <typename FromType>
    struct Adapter<xtal::SymOpType, FromType> {
      xtal::SymOpType operator()(const FromType &adaptable) {
        return xtal::SymOpType(get_matrix(adaptable), get_translation(adaptable), get_time_reversal(adaptable));
      }
    };

    /// Converts any container of any symmetry type with begin() and end() defined into
    /// SymGroupType, as defined by the crystallography module.
    template <typename FromType>
    struct Adapter<xtal::SymGroupType, FromType> {
      template <typename FromTypeIt>
      xtal::SymGroupType operator()(FromTypeIt begin, FromTypeIt end) {
        xtal::SymGroupType casted_group;
        Adapter<xtal::SymOpType, typename FromTypeIt::value_type> to_symop_type;
        for(auto it = begin; it != end; ++it) {
          casted_group.emplace_back(to_symop_type(*it));
        }
        return casted_group;
      }

      xtal::SymGroupType operator()(const FromType &adaptable) {
        return this->operator()(adaptable.begin(), adaptable.end());
      }
    };
  } // namespace adapter
} // namespace CASM

#endif
