#ifndef CASM_monte2_LocalCalculatorTemplate
#define CASM_monte2_LocalCalculatorTemplate

#include "casm/monte2/definitions.hh"

namespace CASM {
namespace Monte2 {

/// \brief Some Monte Carlo methods require "LocalCalculator" methods that
///     implement the interface demonstrated by this class to calculate local
///     quantities such as the activation barrier for a KMC event.
///
/// Notes:
/// - This class is not required to be used, it is meant to demonstate the
///   interface expected by methods that require a "LocalCalculator" type.
/// - Calculators can be `reset` to point at a particular configuration
/// - Then, the calculator methods return values calculated from the
///   configuration currently being pointed at.
template <template ConfigType>
class LocalCalculatorTemplate {
  /// \brief Reset pointer to configuration currently being calculated
  void set(ConfigType const *configuration);

  /// \brief Pointer to configuration currently being calculated
  ConfigType const *get() const;

  /// \brief Value at particular unit cell and phenomenal cluster
  double value(Index unitcell_index, Index equivalent_index);

  /// \brief Return the coordinates of sites (relative to the origin unit cell)
  ///     where a change in DoF values requires this calculators values to be
  ///     re-calculated.
  std::set<UnitCellCoord> required_update_neighborhood(
      Index equivalent_index) const;
};

}  // namespace Monte2
}  // namespace CASM

#endif
