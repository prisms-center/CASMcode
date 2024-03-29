#ifndef CASM_GrandCanonicalConditions_HH
#define CASM_GrandCanonicalConditions_HH

#include <optional>

#include "casm/clex/PrimClex.hh"
#include "casm/global/eigen.hh"
#include "casm/monte_carlo/CorrMatchingPotential.hh"

namespace CASM {
namespace Monte {

class Settings;

/// Conditions for a Grand Canonical run:
/// Temperature
/// Chemical potential
/// Tolerance (for comparing conditions)
///
class GrandCanonicalConditions {
 public:
  /// \brief Default constructor
  GrandCanonicalConditions() {}

  /// \brief Constructor
  ///
  /// \param _primclex PrimClex
  /// \param _temperature in K
  /// \param _param_chem_pot Parametric composition chemical potential
  /// \param _tol tolerance for comparing conditions
  ///
  GrandCanonicalConditions(
      const PrimClex &_primclex, double _temperature,
      const Eigen::VectorXd &_param_chem_pot, double _tol,
      bool _include_formation_energy = true,
      bool _include_param_chem_pot = true,
      std::optional<Eigen::VectorXd> _param_comp_quad_pot_target = std::nullopt,
      std::optional<Eigen::VectorXd> _param_comp_quad_pot_vector = std::nullopt,
      std::optional<Eigen::MatrixXd> _param_comp_quad_pot_matrix = std::nullopt,
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

  /// \brief matrix of exchange chemical potential, M(new, curr) = chem_pot(new)
  /// - chem_pot(curr)
  Eigen::MatrixXd exchange_chem_pot() const;

  /// \brief exchange chemical potential: chem_pot(new) - chem_pot(curr)
  double exchange_chem_pot(Index index_new, Index index_curr) const;

  /// \brief parametric chemical potential: dg/dcomp_x
  Eigen::VectorXd param_chem_pot() const;

  /// \brief parametric chemical potential: dg/dcomp_x(index)
  double param_chem_pot(Index index) const;

  double tolerance() const;

  bool include_formation_energy() const { return m_include_formation_energy; }

  bool include_param_chem_pot() const { return m_include_param_chem_pot; }

  std::optional<Eigen::VectorXd> const &param_comp_quad_pot_target() const {
    return m_param_comp_quad_pot_target;
  }

  std::optional<Eigen::VectorXd> const &param_comp_quad_pot_vector() const {
    return m_param_comp_quad_pot_vector;
  }

  std::optional<Eigen::MatrixXd> const &param_comp_quad_pot_matrix() const {
    return m_param_comp_quad_pot_matrix;
  }

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

  // ***************************************MUTATORS***********************************************
  // //

  /// Set the temperature of the current grand canonical condition.
  void set_temperature(double in_temp);

  /// Set all the parametric chemical potentials of the current grand canonical
  /// condition.
  void set_param_chem_pot(const Eigen::VectorXd &in_chem_pot);

  /// Set a single parametric chemical potential by specifying an index and a
  /// value.
  void set_param_chem_pot(Index ind, double in_chem_pot);

  // ***************************************OPERATORS**********************************************
  // //

  /// Add temperature and all chemical potentials to *this
  GrandCanonicalConditions &operator+=(const GrandCanonicalConditions &RHS);

  /// Add temperature and all chemical potentials together and return a new
  /// Condition
  GrandCanonicalConditions operator+(const GrandCanonicalConditions &RHS) const;

  /// Subtract temperature and all chemical potentials to *this
  GrandCanonicalConditions &operator-=(const GrandCanonicalConditions &RHS);

  /// Subtract temperature and all chemical potentials together and return a new
  /// Condition
  GrandCanonicalConditions operator-(const GrandCanonicalConditions &RHS) const;

  /// Compare temperature and all chemical potentials to *this
  bool operator==(const GrandCanonicalConditions &RHS) const;

  /// Compare temperature and all chemical potentials to *this
  bool operator!=(const GrandCanonicalConditions &RHS) const;

  /// Divide ALL parameters and return the greatest number in absolute value
  int operator/(const GrandCanonicalConditions &RHS_inc) const;

 protected:
  const PrimClex *m_primclex;

  /// Temperature
  double m_temperature;

  /// Inverse temperature. Includes Boltzmann term
  double m_beta;

  /// Vector of the parametric chemical potentials conjugate to the parametric
  /// compositions.
  Eigen::VectorXd m_param_chem_pot;

  /// Matrix(i,j) of chem_pot(i) - chem_pot(j)
  Eigen::MatrixXd m_exchange_chem_pot;

  /// Tolerance for comparison operators == and !=
  double m_tolerance;

  bool m_include_formation_energy;
  bool m_include_param_chem_pot;

  std::optional<Eigen::VectorXd> m_param_comp_quad_pot_target;
  std::optional<Eigen::VectorXd> m_param_comp_quad_pot_vector;
  std::optional<Eigen::MatrixXd> m_param_comp_quad_pot_matrix;

  std::optional<Eigen::VectorXd> m_order_parameter_pot;
  std::optional<Eigen::VectorXd> m_order_parameter_quad_pot_target;
  std::optional<Eigen::VectorXd> m_order_parameter_quad_pot_vector;
  std::optional<Eigen::MatrixXd> m_order_parameter_quad_pot_matrix;

  std::optional<CorrMatchingParams> m_corr_matching_pot;
  std::optional<RandomAlloyCorrMatchingParams> m_random_alloy_corr_matching_pot;
};

std::ostream &operator<<(std::ostream &sout,
                         const GrandCanonicalConditions &cond);

}  // namespace Monte
}  // namespace CASM

#endif
