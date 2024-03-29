#include "casm/monte_carlo/canonical/CanonicalConditions.hh"

#include "casm/casm_io/container/json_io.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/monte_carlo/MonteCarlo.hh"
#include "casm/monte_carlo/MonteSettings.hh"
#include "casm/monte_carlo/canonical/CanonicalIO.hh"
#include "casm/monte_carlo/conditions_functions.hh"

namespace CASM {
namespace Monte {

/// \brief Constructor
///
/// \param _primclex PrimClex
/// \param _temperature in K
/// \param _param_chem_pot Parametric composition chemical potential
/// \param _tol tolerance for comparing conditions
///
CanonicalConditions::CanonicalConditions(
    const PrimClex &_primclex, double _temperature,
    const Eigen::VectorXd &_param_comp, double _tol,
    bool _include_formation_energy,
    std::optional<Eigen::VectorXd> _order_parameter_pot,
    std::optional<Eigen::VectorXd> _order_parameter_quad_pot_target,
    std::optional<Eigen::VectorXd> _order_parameter_quad_pot_vector,
    std::optional<Eigen::MatrixXd> _order_parameter_quad_pot_matrix,
    std::optional<CorrMatchingParams> _corr_matching_pot,
    std::optional<RandomAlloyCorrMatchingParams>
        _random_alloy_corr_matching_pot)
    : m_primclex(&_primclex),
      m_temperature(_temperature),
      m_beta(1.0 / (KB * m_temperature)),
      m_param_composition(_param_comp),
      m_tolerance(_tol),
      m_include_formation_energy(_include_formation_energy),
      m_order_parameter_pot(_order_parameter_pot),
      m_order_parameter_quad_pot_target(_order_parameter_quad_pot_target),
      m_order_parameter_quad_pot_vector(_order_parameter_quad_pot_vector),
      m_order_parameter_quad_pot_matrix(_order_parameter_quad_pot_matrix),
      m_corr_matching_pot(_corr_matching_pot),
      m_random_alloy_corr_matching_pot(_random_alloy_corr_matching_pot) {}

// ***************************************ACCESSORS**********************************************
// //

const PrimClex &CanonicalConditions::primclex() const { return *m_primclex; }

double CanonicalConditions::temperature() const { return m_temperature; }

double CanonicalConditions::beta() const { return m_beta; }

/// \brief parametric composition: comp_x
Eigen::VectorXd CanonicalConditions::param_composition() const {
  return m_param_composition;
}

/// \brief parametric composition: dcomp_x(index)
double CanonicalConditions::param_composition(Index index) const {
  return m_param_composition(index);
}

/// \brief mol composition: comp_n
Eigen::VectorXd CanonicalConditions::mol_composition() const {
  return primclex().composition_axes().mol_composition(m_param_composition);
}

/// \brief mol composition: comp_n(index)
double CanonicalConditions::mol_composition(Index index) const {
  return mol_composition()(index);
}

double CanonicalConditions::tolerance() const { return m_tolerance; }

// ***************************************OPERATORS**********************************************
// //

CanonicalConditions &CanonicalConditions::operator+=(
    const CanonicalConditions &RHS) {
  m_temperature += RHS.m_temperature;
  m_beta = 1.0 / (CASM::KB * m_temperature);
  m_param_composition += RHS.m_param_composition;
  conditions::incr(m_order_parameter_pot, RHS.m_order_parameter_pot);
  conditions::incr(m_order_parameter_quad_pot_target,
                   RHS.m_order_parameter_quad_pot_target);
  conditions::incr(m_order_parameter_quad_pot_vector,
                   RHS.m_order_parameter_quad_pot_vector);
  conditions::incr(m_order_parameter_quad_pot_matrix,
                   RHS.m_order_parameter_quad_pot_matrix);

  Monte::incr(m_corr_matching_pot, RHS.m_corr_matching_pot);
  Monte::incr(m_random_alloy_corr_matching_pot,
              RHS.m_random_alloy_corr_matching_pot);
  return *this;
}

CanonicalConditions CanonicalConditions::operator+(
    const CanonicalConditions &RHS) const {
  return CanonicalConditions(*this) += RHS;
}

/// Subtract temperature and all chemical potentials to *this
CanonicalConditions &CanonicalConditions::operator-=(
    const CanonicalConditions &RHS) {
  m_temperature -= RHS.m_temperature;
  m_beta = 1.0 / (CASM::KB * m_temperature);
  m_param_composition -= RHS.m_param_composition;
  conditions::decr(m_order_parameter_pot, RHS.m_order_parameter_pot);
  conditions::decr(m_order_parameter_quad_pot_target,
                   RHS.m_order_parameter_quad_pot_target);
  conditions::decr(m_order_parameter_quad_pot_vector,
                   RHS.m_order_parameter_quad_pot_vector);
  conditions::decr(m_order_parameter_quad_pot_matrix,
                   RHS.m_order_parameter_quad_pot_matrix);

  Monte::decr(m_corr_matching_pot, RHS.m_corr_matching_pot);
  Monte::decr(m_random_alloy_corr_matching_pot,
              RHS.m_random_alloy_corr_matching_pot);
  return *this;
}

CanonicalConditions CanonicalConditions::operator-(
    const CanonicalConditions &RHS) const {
  return CanonicalConditions(*this) -= RHS;
}

bool CanonicalConditions::operator==(const CanonicalConditions &RHS) const {
  if (!CASM::almost_equal(m_temperature, RHS.m_temperature, m_tolerance)) {
    return false;
  }
  if (!CASM::almost_equal(m_param_composition, RHS.m_param_composition,
                          m_tolerance)) {
    return false;
  }
  if (!conditions::almost_equal(m_order_parameter_pot,
                                RHS.m_order_parameter_pot, m_tolerance)) {
    return false;
  }
  if (!conditions::almost_equal(m_order_parameter_quad_pot_target,
                                RHS.m_order_parameter_quad_pot_target,
                                m_tolerance)) {
    return false;
  }
  if (!conditions::almost_equal(m_order_parameter_quad_pot_vector,
                                RHS.m_order_parameter_quad_pot_vector,
                                m_tolerance)) {
    return false;
  }
  if (!conditions::almost_equal(m_order_parameter_quad_pot_matrix,
                                RHS.m_order_parameter_quad_pot_matrix,
                                m_tolerance)) {
    return false;
  }

  if (!Monte::almost_equal(m_corr_matching_pot, RHS.m_corr_matching_pot,
                           m_tolerance)) {
    return false;
  }
  if (!Monte::almost_equal(m_random_alloy_corr_matching_pot,
                           RHS.m_random_alloy_corr_matching_pot, m_tolerance)) {
    return false;
  }

  return true;
}

bool CanonicalConditions::operator!=(const CanonicalConditions &RHS) const {
  return !(*this == RHS);
}

int CanonicalConditions::operator/(const CanonicalConditions &RHS_inc) const {
  int max_division = 0;

  if (!almost_zero(RHS_inc.temperature())) {
    max_division = round(temperature() / RHS_inc.temperature());
  }

  for (Index i = 0; i < m_param_composition.size(); i++) {
    int temp_division =
        round(m_param_composition(i) / RHS_inc.m_param_composition(i));

    if (temp_division > max_division &&
        !almost_zero(RHS_inc.m_param_composition(i))) {
      max_division = temp_division;
    }
  }

  conditions::find_max_division(max_division, m_order_parameter_pot,
                                RHS_inc.m_order_parameter_pot);

  conditions::find_max_division(max_division, m_order_parameter_quad_pot_target,
                                RHS_inc.m_order_parameter_quad_pot_target);
  conditions::find_max_division(max_division, m_order_parameter_quad_pot_vector,
                                RHS_inc.m_order_parameter_quad_pot_vector);
  conditions::find_max_division(max_division, m_order_parameter_quad_pot_matrix,
                                RHS_inc.m_order_parameter_quad_pot_matrix);

  Monte::find_max_division(max_division, m_corr_matching_pot,
                           RHS_inc.m_corr_matching_pot);
  Monte::find_max_division(max_division, m_random_alloy_corr_matching_pot,
                           RHS_inc.m_random_alloy_corr_matching_pot);

  return max_division;
}

std::ostream &operator<<(std::ostream &sout, const CanonicalConditions &cond) {
  jsonParser json;
  to_json(cond, json);
  sout << json << std::endl;
  return sout;
}

}  // namespace Monte
}  // namespace CASM
