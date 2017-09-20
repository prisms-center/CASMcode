#include "casm/monte_carlo/grand_canonical/GrandCanonicalConditions.hh"
#include "casm/monte_carlo/MonteCarlo.hh"
#include "casm/monte_carlo/MonteSettings.hh"

#include "casm/clex/PrimClex.hh"

namespace CASM {


  /// \brief Constructor
  ///
  /// \param _primclex PrimClex
  /// \param _temperature in K
  /// \param _param_chem_pot Parametric composition chemical potential
  /// \param _tol tolerance for comparing conditions
  ///
  GrandCanonicalConditions::GrandCanonicalConditions(
    const PrimClex &_primclex,
    double _temperature,
    const Eigen::VectorXd &_param_chem_pot,
    double _tol) :

    m_primclex(&_primclex),
    m_tolerance(_tol) {

    // -- set T ----
    set_temperature(_temperature);


    // -- set param_chem_pot ----
    set_param_chem_pot(_param_chem_pot);

  }

  // ***************************************ACCESSORS********************************************** //

  const PrimClex &GrandCanonicalConditions::primclex() const {
    return *m_primclex;
  }

  double GrandCanonicalConditions::temperature() const {
    return m_temperature;
  }

  double GrandCanonicalConditions::beta() const {
    return m_beta;
  }

  Eigen::MatrixXd GrandCanonicalConditions::exchange_chem_pot() const {
    return m_exchange_chem_pot;
  }

  double GrandCanonicalConditions::exchange_chem_pot(Index index_new, Index index_curr) const {
    return m_exchange_chem_pot(index_new, index_curr);
  }

  Eigen::VectorXd GrandCanonicalConditions::param_chem_pot() const {
    return m_param_chem_pot;
  }

  double GrandCanonicalConditions::param_chem_pot(Index index) const {
    return m_param_chem_pot(index);
  }

  double GrandCanonicalConditions::tolerance() const {
    return m_tolerance;
  }


  // ***************************************MUTATORS*********************************************** //

  void GrandCanonicalConditions::set_temperature(double in_temp) {
    m_temperature = in_temp;
    m_beta = 1.0 / (KB * m_temperature);
    return;
  }

  void GrandCanonicalConditions::set_param_chem_pot(const Eigen::VectorXd &in_param_chem_pot) {
    m_param_chem_pot = in_param_chem_pot;
    m_chem_pot = primclex().composition_axes().dparam_dmol().transpose() * m_param_chem_pot;

    int Ncomp = primclex().composition_axes().components().size();
    m_exchange_chem_pot = Eigen::MatrixXd(Ncomp, Ncomp);
    for(int index_new = 0; index_new < Ncomp; ++index_new) {
      for(int index_curr = 0; index_curr < Ncomp; ++index_curr) {
        Eigen::VectorXl dn = Eigen::VectorXl::Zero(Ncomp);
        dn(index_new) += 1;
        dn(index_curr) -= 1;
        m_exchange_chem_pot(index_new, index_curr) =
          m_param_chem_pot.transpose() * primclex().composition_axes().dparam_dmol() * dn.cast<double>();
      }
    }

    return;
  }

  void GrandCanonicalConditions::set_param_chem_pot(Index ind, double in_param_chem_pot) {
    m_param_chem_pot(ind) = in_param_chem_pot;
    set_param_chem_pot(m_param_chem_pot);
    return;
  }


  // ***************************************OPERATORS********************************************** //

  GrandCanonicalConditions &GrandCanonicalConditions::operator+=(const GrandCanonicalConditions &RHS) {
    m_temperature += RHS.m_temperature;
    for(int i = 0; i < m_param_chem_pot.size(); i++) {
      m_param_chem_pot(i) += RHS.m_param_chem_pot(i);
    }
    set_param_chem_pot(m_param_chem_pot);
    m_beta = 1.0 / (CASM::KB * m_temperature);
    return *this;
  }

  GrandCanonicalConditions GrandCanonicalConditions::operator+(const GrandCanonicalConditions &RHS) const {
    return GrandCanonicalConditions(*this) += RHS;
  }

  ///Subtract temperature and all chemical potentials to *this
  GrandCanonicalConditions &GrandCanonicalConditions::operator-=(const GrandCanonicalConditions &RHS) {
    m_temperature -= RHS.m_temperature;
    for(int i = 0; i < m_param_chem_pot.size(); i++) {
      m_param_chem_pot(i) -= RHS.m_param_chem_pot(i);
    }
    set_param_chem_pot(m_param_chem_pot);
    m_beta = 1.0 / (CASM::KB * m_temperature);
    return *this;
  }

  GrandCanonicalConditions GrandCanonicalConditions::operator-(const GrandCanonicalConditions &RHS) const {
    return GrandCanonicalConditions(*this) -= RHS;
  }

  bool GrandCanonicalConditions::operator==(const GrandCanonicalConditions &RHS) const {
    if(!almost_zero(m_temperature - RHS.m_temperature, m_tolerance)) {
      return false;
    }

    for(int i = 0; i < m_param_chem_pot.size(); i++) {
      if(!almost_zero(m_param_chem_pot(i) - RHS.m_param_chem_pot(i), m_tolerance)) {
        return false;
      }
    }

    return true;
  }

  bool GrandCanonicalConditions::operator!=(const GrandCanonicalConditions &RHS) const {
    return !(*this == RHS);
  }

  int GrandCanonicalConditions::operator/(const GrandCanonicalConditions &RHS_inc) const {
    int max_division = 0;

    if(!almost_zero(RHS_inc.temperature())) {
      max_division = round(temperature() / RHS_inc.temperature());
    }

    for(Index i = 0; i < param_chem_pot().size(); i++) {
      int temp_division = round(param_chem_pot(i) / RHS_inc.param_chem_pot(i));

      if(temp_division > max_division && !almost_zero(RHS_inc.param_chem_pot(i))) {
        max_division = temp_division;
      }
    }

    return max_division;
  }

  std::ostream &operator<<(std::ostream &sout, const GrandCanonicalConditions &cond) {
    sout << "T: " << cond.temperature() << "\n";
    for(int i = 0; i < cond.param_chem_pot().size(); i++) {
      jsonParser json;
      sout << "param_chem_pot: " << to_json_array(cond.param_chem_pot(), json) << "\n";
    }
    return sout;
  }


}


