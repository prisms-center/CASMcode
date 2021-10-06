#ifndef CASM_monte2_CalculatorTemplate
#define CASM_monte2_CalculatorTemplate

#include "casm/monte2/definitions.hh"

namespace CASM {
namespace Monte2 {

/// \brief Monte Carlo methods require "Calculator" methods that implement the
///     interface demonstrated by this class to calculate quantities such as
///     the formation energy and change in the formation energy.
///
/// Notes:
/// - This class is not required to be used, it is meant to demonstate the
///   interface expected by methods that require a "Calculator" type.
/// - Calculators can be `reset` to point at a particular configuration
/// - Then, the calculator methods return values calculated from the
///   configuration currently being pointed at.
/// - Monte Carlo methods may not require all methods (i.e. `local_delta_value`
///   does not need to be implemented for an occupation-only Monte Carlo method)
template <template ConfigType>
class CalculatorTemplate {
  /// \brief Reset pointer to configuration currently being calculated
  void set(ConfigType const *configuration);

  /// \brief Pointer to configuration currently being calculated
  ConfigType const *get() const;

  /// \brief Intensive value (value per primitive cell)
  double intensive_value();

  /// \brief Extensive value (sum over entire configuration)
  double extensive_value();

  /// \brief Change in extensive value due to occupation change
  double occ_delta_value(Index linear_site_index, int new_occ);

  /// \brief Change in extensive value due to local DoF value change
  double local_delta_value(DoFKey const &key, Index linear_site_index,
                           Eigen::VectorXd const &new_value);

  /// \brief Change in extensive value due to global DoF value change
  double global_delta_value(DoFKey const &key,
                            Eigen::VectorXd const &new_value);

  /// \brief Return the coordinates of sites (relative to the origin unit cell)
  ///     where a change in DoF values requires this calculators values to be
  ///     re-calculated.
  std::set<UnitCellCoord> required_update_neighborhood() const;
};

}  // namespace Monte2
}  // namespace CASM

#endif
