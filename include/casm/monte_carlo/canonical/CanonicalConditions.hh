#ifndef CASM_CanonicalConditions_HH
#define CASM_CanonicalConditions_HH

#include <optional>

#include "casm/global/definitions.hh"
#include "casm/global/eigen.hh"
#include "casm/monte_carlo/CorrMatchingPotential.hh"

namespace CASM {
class PrimClex;
}

namespace CASM {
namespace Monte {

struct CorrMatchingParams;
class Settings;

/// Conditions for a Canonical run:
/// Temperature
/// Composition
/// Tolerance (for comparing conditions)
///
class CanonicalConditions {
 public:
  /// \brief Default constructor
  CanonicalConditions() {}

  /// \brief Constructor
  ///
  /// \param _primclex PrimClex
  /// \param _temperature in K
  /// \param _param_comp Parametric composition
  /// \param _tol tolerance for comparing conditions
  ///
  CanonicalConditions(
      const PrimClex &_primclex, double _temperature,
      const Eigen::VectorXd &_param_comp, double _tol,
      bool _include_formation_energy = true,
      std::optional<Eigen::VectorXd> _order_parameter_pot = std::nullopt,
      std::optional<Eigen::VectorXd> _order_parameter_quad_pot_target =
          std::nullopt,
      std::optional<Eigen::VectorXd> _order_parameter_quad_pot_vector =
          std::nullopt,
      std::optional<Eigen::MatrixXd> _order_parameter_quad_pot_matrix =
          std::nullopt,
      std::optional<CorrMatchingParams> _corr_matching_pot = std::nullopt,
      std::optional<RandomAlloyCorrMatchingParams>
          _random_alloy_corr_matching_pot = std::nullopt);

  // ***************************************ACCESSORS**********************************************
  // //

  const PrimClex &primclex() const;

  double temperature() const;

  double beta() const;

  /// \brief parametric composition: comp_x
  Eigen::VectorXd param_composition() const;

  /// \brief parametric composition: dcomp_x(index)
  double param_composition(Index index) const;

  /// \brief mol composition: comp_n
  Eigen::VectorXd mol_composition() const;

  /// \brief mol composition: comp_n(index)
  double mol_composition(Index index) const;

  double tolerance() const;

  bool include_formation_energy() const { return m_include_formation_energy; }

  std::optional<Eigen::VectorXd> const &order_parameter_pot() const {
    return m_order_parameter_pot;
  }

  std::optional<Eigen::VectorXd> const &order_parameter_quad_pot_target()
      const {
    return m_order_parameter_quad_pot_target;
  }

  std::optional<Eigen::VectorXd> const &order_parameter_quad_pot_vector()
      const {
    return m_order_parameter_quad_pot_vector;
  }

  std::optional<Eigen::MatrixXd> const &order_parameter_quad_pot_matrix()
      const {
    return m_order_parameter_quad_pot_matrix;
  }

  std::optional<CorrMatchingParams> const &corr_matching_pot() const {
    return m_corr_matching_pot;
  }

  std::optional<RandomAlloyCorrMatchingParams> const &
  random_alloy_corr_matching_pot() const {
    return m_random_alloy_corr_matching_pot;
  }

  // ***************************************OPERATORS**********************************************
  // //

  /// Add temperature and all chemical potentials to *this
  CanonicalConditions &operator+=(const CanonicalConditions &RHS);

  /// Add temperature and all chemical potentials together and return a new
  /// Condition
  CanonicalConditions operator+(const CanonicalConditions &RHS) const;

  /// Subtract temperature and all chemical potentials to *this
  CanonicalConditions &operator-=(const CanonicalConditions &RHS);

  /// Subtract temperature and all chemical potentials together and return a new
  /// Condition
  CanonicalConditions operator-(const CanonicalConditions &RHS) const;

  /// Compare temperature and all chemical potentials to *this
  bool operator==(const CanonicalConditions &RHS) const;

  /// Compare temperature and all chemical potentials to *this
  bool operator!=(const CanonicalConditions &RHS) const;

  /// Divide ALL parameters and return the greatest number in absolute value
  int operator/(const CanonicalConditions &RHS_inc) const;

 protected:
  const PrimClex *m_primclex;

  /// Temperature
  double m_temperature;

  /// Inverse temperature. Includes Boltzmann term
  double m_beta;

  /// Vector of the param composition
  Eigen::VectorXd m_param_composition;

  /// Tolerance for comparison operators == and !=
  double m_tolerance;

  bool m_include_formation_energy;

  std::optional<Eigen::VectorXd> m_order_parameter_pot;
  std::optional<Eigen::VectorXd> m_order_parameter_quad_pot_target;
  std::optional<Eigen::VectorXd> m_order_parameter_quad_pot_vector;
  std::optional<Eigen::MatrixXd> m_order_parameter_quad_pot_matrix;

  std::optional<CorrMatchingParams> m_corr_matching_pot;
  std::optional<RandomAlloyCorrMatchingParams> m_random_alloy_corr_matching_pot;
};

std::ostream &operator<<(std::ostream &sout, const CanonicalConditions &cond);

}  // namespace Monte
}  // namespace CASM

#endif
