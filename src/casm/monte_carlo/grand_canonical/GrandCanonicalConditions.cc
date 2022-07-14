#include "casm/monte_carlo/grand_canonical/GrandCanonicalConditions.hh"

#include "casm/casm_io/container/json_io.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/monte_carlo/MonteCarlo.hh"
#include "casm/monte_carlo/MonteSettings.hh"
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
GrandCanonicalConditions::GrandCanonicalConditions(
    const PrimClex &_primclex, double _temperature,
    const Eigen::VectorXd &_param_chem_pot, double _tol,
    std::optional<Eigen::VectorXd> _param_comp_quad_pot_target,
    std::optional<Eigen::VectorXd> _param_comp_quad_pot_vector,
    std::optional<Eigen::MatrixXd> _param_comp_quad_pot_matrix,
    std::optional<Eigen::VectorXd> _order_parameter_pot,
    std::optional<Eigen::VectorXd> _order_parameter_quad_pot_target,
    std::optional<Eigen::VectorXd> _order_parameter_quad_pot_vector,
    std::optional<Eigen::MatrixXd> _order_parameter_quad_pot_matrix)
    :

      m_primclex(&_primclex),
      m_tolerance(_tol),
      m_param_comp_quad_pot_target(_param_comp_quad_pot_target),
      m_param_comp_quad_pot_vector(_param_comp_quad_pot_vector),
      m_param_comp_quad_pot_matrix(_param_comp_quad_pot_matrix),
      m_order_parameter_pot(_order_parameter_pot),
      m_order_parameter_quad_pot_target(_order_parameter_quad_pot_target),
      m_order_parameter_quad_pot_vector(_order_parameter_quad_pot_vector),
      m_order_parameter_quad_pot_matrix(_order_parameter_quad_pot_matrix) {
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
  set_param_chem_pot(m_param_chem_pot + RHS.m_param_chem_pot);
  incr(m_order_parameter_pot, RHS.m_order_parameter_pot);
  incr(m_order_parameter_quad_pot_target,
       RHS.m_order_parameter_quad_pot_target);
  incr(m_order_parameter_quad_pot_vector,
       RHS.m_order_parameter_quad_pot_vector);
  incr(m_order_parameter_quad_pot_matrix,
       RHS.m_order_parameter_quad_pot_matrix);
  m_beta = 1.0 / (CASM::KB * m_temperature);

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
  set_param_chem_pot(m_param_chem_pot - RHS.m_param_chem_pot);
  decr(m_order_parameter_pot, RHS.m_order_parameter_pot);
  decr(m_order_parameter_quad_pot_target,
       RHS.m_order_parameter_quad_pot_target);
  decr(m_order_parameter_quad_pot_vector,
       RHS.m_order_parameter_quad_pot_vector);
  decr(m_order_parameter_quad_pot_matrix,
       RHS.m_order_parameter_quad_pot_matrix);
  m_beta = 1.0 / (CASM::KB * m_temperature);
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
  if (!almost_equal(m_param_chem_pot, RHS.m_param_chem_pot, m_tolerance)) {
    return false;
  }
  if (!almost_equal(m_order_parameter_pot, RHS.m_order_parameter_pot,
                    m_tolerance)) {
    return false;
  }
  if (!almost_equal(m_order_parameter_quad_pot_target,
                    RHS.m_order_parameter_quad_pot_target, m_tolerance)) {
    return false;
  }
  if (!almost_equal(m_order_parameter_quad_pot_vector,
                    RHS.m_order_parameter_quad_pot_vector, m_tolerance)) {
    return false;
  }
  if (!almost_equal(m_order_parameter_quad_pot_matrix,
                    RHS.m_order_parameter_quad_pot_matrix, m_tolerance)) {
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

  find_max_division(max_division, m_order_parameter_pot,
                    RHS_inc.m_order_parameter_pot);
  find_max_division(max_division, m_order_parameter_quad_pot_target,
                    RHS_inc.m_order_parameter_quad_pot_target);
  find_max_division(max_division, m_order_parameter_quad_pot_vector,
                    RHS_inc.m_order_parameter_quad_pot_vector);
  find_max_division(max_division, m_order_parameter_quad_pot_matrix,
                    RHS_inc.m_order_parameter_quad_pot_matrix);

  return max_division;
}

std::ostream &operator<<(std::ostream &sout,
                         const GrandCanonicalConditions &cond) {
  sout << "T: " << cond.temperature() << "\n";
  for (int i = 0; i < cond.param_chem_pot().size(); i++) {
    jsonParser json;
    sout << "param_chem_pot: " << to_json_array(cond.param_chem_pot(), json)
         << "\n";
  }
  return sout;
}

}  // namespace Monte
}  // namespace CASM
