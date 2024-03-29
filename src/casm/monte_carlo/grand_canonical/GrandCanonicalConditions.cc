#include "casm/monte_carlo/grand_canonical/GrandCanonicalConditions.hh"

#include "casm/casm_io/container/json_io.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/monte_carlo/MonteCarlo.hh"
#include "casm/monte_carlo/MonteSettings.hh"
#include "casm/monte_carlo/conditions_functions.hh"
#include "casm/monte_carlo/grand_canonical/GrandCanonicalIO.hh"

namespace CASM {
namespace Monte {

/// \brief Constructor
///
/// \param _primclex PrimClex
/// \param _temperature in K
/// \param _param_chem_pot Parametric composition chemical potential
/// \param _tol tolerance for comparing conditions
///
GrandCanonicalConditions::GrandCanonicalConditions(
    const PrimClex &_primclex, double _temperature,
    const Eigen::VectorXd &_param_chem_pot, double _tol,
    bool _include_formation_energy, bool _include_param_chem_pot,
    std::optional<Eigen::VectorXd> _param_comp_quad_pot_target,
    std::optional<Eigen::VectorXd> _param_comp_quad_pot_vector,
    std::optional<Eigen::MatrixXd> _param_comp_quad_pot_matrix,
    std::optional<Eigen::VectorXd> _order_parameter_pot,
    std::optional<Eigen::VectorXd> _order_parameter_quad_pot_target,
    std::optional<Eigen::VectorXd> _order_parameter_quad_pot_vector,
    std::optional<Eigen::MatrixXd> _order_parameter_quad_pot_matrix,
    std::optional<CorrMatchingParams> _corr_matching_pot,
    std::optional<RandomAlloyCorrMatchingParams>
        _random_alloy_corr_matching_pot)
    : m_primclex(&_primclex),
      m_tolerance(_tol),
      m_include_formation_energy(_include_formation_energy),
      m_include_param_chem_pot(_include_param_chem_pot),
      m_param_comp_quad_pot_target(_param_comp_quad_pot_target),
      m_param_comp_quad_pot_vector(_param_comp_quad_pot_vector),
      m_param_comp_quad_pot_matrix(_param_comp_quad_pot_matrix),
      m_order_parameter_pot(_order_parameter_pot),
      m_order_parameter_quad_pot_target(_order_parameter_quad_pot_target),
      m_order_parameter_quad_pot_vector(_order_parameter_quad_pot_vector),
      m_order_parameter_quad_pot_matrix(_order_parameter_quad_pot_matrix),
      m_corr_matching_pot(_corr_matching_pot),
      m_random_alloy_corr_matching_pot(_random_alloy_corr_matching_pot) {
  // -- set T ----
  set_temperature(_temperature);

  // -- set param_chem_pot ----
  set_param_chem_pot(_param_chem_pot);
}

// ***************************************ACCESSORS**********************************************
// //

const PrimClex &GrandCanonicalConditions::primclex() const {
  return *m_primclex;
}

double GrandCanonicalConditions::temperature() const { return m_temperature; }

double GrandCanonicalConditions::beta() const { return m_beta; }

Eigen::MatrixXd GrandCanonicalConditions::exchange_chem_pot() const {
  return m_exchange_chem_pot;
}

double GrandCanonicalConditions::exchange_chem_pot(Index index_new,
                                                   Index index_curr) const {
  return m_exchange_chem_pot(index_new, index_curr);
}

Eigen::VectorXd GrandCanonicalConditions::param_chem_pot() const {
  return m_param_chem_pot;
}

double GrandCanonicalConditions::param_chem_pot(Index index) const {
  return m_param_chem_pot(index);
}

double GrandCanonicalConditions::tolerance() const { return m_tolerance; }

// ***************************************MUTATORS***********************************************
// //

void GrandCanonicalConditions::set_temperature(double in_temp) {
  m_temperature = in_temp;
  m_beta = 1.0 / (KB * m_temperature);
  return;
}

void GrandCanonicalConditions::set_param_chem_pot(
    const Eigen::VectorXd &in_param_chem_pot) {
  m_param_chem_pot = in_param_chem_pot;

  int Ncomp = primclex().composition_axes().components().size();
  m_exchange_chem_pot = Eigen::MatrixXd(Ncomp, Ncomp);
  for (int index_new = 0; index_new < Ncomp; ++index_new) {
    for (int index_curr = 0; index_curr < Ncomp; ++index_curr) {
      Eigen::VectorXl dn = Eigen::VectorXl::Zero(Ncomp);
      dn(index_new) += 1;
      dn(index_curr) -= 1;
      m_exchange_chem_pot(index_new, index_curr) =
          m_param_chem_pot.transpose() *
          primclex().composition_axes().dparam_dmol() * dn.cast<double>();
    }
  }

  return;
}

void GrandCanonicalConditions::set_param_chem_pot(Index ind,
                                                  double in_param_chem_pot) {
  m_param_chem_pot(ind) = in_param_chem_pot;
  set_param_chem_pot(m_param_chem_pot);
  return;
}

// ***************************************OPERATORS**********************************************
// //

GrandCanonicalConditions &GrandCanonicalConditions::operator+=(
    const GrandCanonicalConditions &RHS) {
  m_temperature += RHS.m_temperature;
  m_beta = 1.0 / (CASM::KB * m_temperature);
  set_param_chem_pot(m_param_chem_pot + RHS.m_param_chem_pot);

  conditions::incr(m_param_comp_quad_pot_target,
                   RHS.m_param_comp_quad_pot_target);
  conditions::incr(m_param_comp_quad_pot_vector,
                   RHS.m_param_comp_quad_pot_vector);
  conditions::incr(m_param_comp_quad_pot_matrix,
                   RHS.m_param_comp_quad_pot_matrix);

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

GrandCanonicalConditions GrandCanonicalConditions::operator+(
    const GrandCanonicalConditions &RHS) const {
  return GrandCanonicalConditions(*this) += RHS;
}

/// Subtract temperature and all chemical potentials to *this
GrandCanonicalConditions &GrandCanonicalConditions::operator-=(
    const GrandCanonicalConditions &RHS) {
  m_temperature -= RHS.m_temperature;
  m_beta = 1.0 / (CASM::KB * m_temperature);
  set_param_chem_pot(m_param_chem_pot - RHS.m_param_chem_pot);

  conditions::decr(m_param_comp_quad_pot_target,
                   RHS.m_param_comp_quad_pot_target);
  conditions::decr(m_param_comp_quad_pot_vector,
                   RHS.m_param_comp_quad_pot_vector);
  conditions::decr(m_param_comp_quad_pot_matrix,
                   RHS.m_param_comp_quad_pot_matrix);

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

GrandCanonicalConditions GrandCanonicalConditions::operator-(
    const GrandCanonicalConditions &RHS) const {
  return GrandCanonicalConditions(*this) -= RHS;
}

bool GrandCanonicalConditions::operator==(
    const GrandCanonicalConditions &RHS) const {
  if (!CASM::almost_equal(m_temperature, RHS.m_temperature, m_tolerance)) {
    return false;
  }
  if (!CASM::almost_equal(m_param_chem_pot, RHS.m_param_chem_pot,
                          m_tolerance)) {
    return false;
  }

  if (!conditions::almost_equal(m_param_comp_quad_pot_target,
                                RHS.m_param_comp_quad_pot_target,
                                m_tolerance)) {
    return false;
  }
  if (!conditions::almost_equal(m_param_comp_quad_pot_vector,
                                RHS.m_param_comp_quad_pot_vector,
                                m_tolerance)) {
    return false;
  }
  if (!conditions::almost_equal(m_param_comp_quad_pot_matrix,
                                RHS.m_param_comp_quad_pot_matrix,
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

bool GrandCanonicalConditions::operator!=(
    const GrandCanonicalConditions &RHS) const {
  return !(*this == RHS);
}

int GrandCanonicalConditions::operator/(
    const GrandCanonicalConditions &RHS_inc) const {
  int max_division = 0;

  if (!almost_zero(RHS_inc.temperature())) {
    max_division = round(temperature() / RHS_inc.temperature());
  }

  for (Index i = 0; i < param_chem_pot().size(); i++) {
    int temp_division = round(param_chem_pot(i) / RHS_inc.param_chem_pot(i));

    if (temp_division > max_division &&
        !almost_zero(RHS_inc.param_chem_pot(i))) {
      max_division = temp_division;
    }
  }

  conditions::find_max_division(max_division, m_param_comp_quad_pot_target,
                                RHS_inc.m_param_comp_quad_pot_target);
  conditions::find_max_division(max_division, m_param_comp_quad_pot_vector,
                                RHS_inc.m_param_comp_quad_pot_vector);
  conditions::find_max_division(max_division, m_param_comp_quad_pot_matrix,
                                RHS_inc.m_param_comp_quad_pot_matrix);

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

std::ostream &operator<<(std::ostream &sout,
                         const GrandCanonicalConditions &cond) {
  jsonParser json;
  to_json(cond, json);
  sout << json << std::endl;
  return sout;
}

}  // namespace Monte
}  // namespace CASM
