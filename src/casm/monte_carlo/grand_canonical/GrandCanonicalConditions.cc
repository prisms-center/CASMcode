#include "casm/monte_carlo/grand_canonical/GrandCanonicalConditions.hh"
#include "casm/monte_carlo/MonteCarlo.hh"
#include "casm/monte_carlo/MonteSettings.hh"

#include "casm/clex/PrimClex.hh"

namespace CASM {
  
  /// \brief Construct from MonteSettings
  ///
  /// Use MonteSettings::initial_conditions() as the input to this so that that the 
  /// correct T, mu, and tol settings can be read
  ///
  GrandCanonicalConditions::GrandCanonicalConditions(const PrimClex &primclex, const MonteSettings &settings) : 
    m_tolerance(settings.tolerance()) {
    
    // -- set T ----
    
    set_temperature(settings.temperature());

    
    // -- set mu ----
    
    if(!primclex.has_composition_axes()) {
      throw std::runtime_error(std::string("ERROR constructing GrandCanonicalConditions: No composition axes."));
    }
    m_comp_converter = primclex.composition_axes();
    
    int Nparam = primclex.composition_axes().independent_compositions();
    Eigen::VectorXd param_mu(Nparam);

    int start = (int) 'a';

    for(int i = 0; i < Nparam; i++) {
      std::stringstream ss;
      ss << (char)(start + i);
      param_mu(i) = settings.chemical_potential(ss.str());
    }
    
    
    Eigen::VectorXd atomic_mu = primclex.composition_axes().atomic_mu(param_mu);
    
    set_mu(atomic_mu);

  }
  
  // ***************************************ACCESSORS********************************************** //

  double GrandCanonicalConditions::temperature() const {
    return m_temperature;
  }

  double GrandCanonicalConditions::beta() const {
    return m_beta;
  }

  const Eigen::VectorXd& GrandCanonicalConditions::mu() const {
    return m_mu;
  }

  double GrandCanonicalConditions::mu(Index mu_index) const {
    return m_mu(mu_index);
  }
  
  Eigen::VectorXd GrandCanonicalConditions::param_mu() const {
    return m_comp_converter.param_mu(m_mu);
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

  void GrandCanonicalConditions::set_mu(const Eigen::VectorXd &in_mu) {
    m_mu = in_mu;
    return;
  }

  void GrandCanonicalConditions::set_mu(Index ind, double in_mu) {
    m_mu(ind) = in_mu;
    return;
  }

  void GrandCanonicalConditions::increment_by(const GrandCanonicalConditions &cond_increment) {
    m_temperature += cond_increment.m_temperature;
    for(int i = 0; i < m_mu.size(); i++) {
      m_mu(i) += cond_increment.m_mu(i);
    }
    m_beta = 1.0 / (CASM::KB * m_temperature);
    return;
  }

  void GrandCanonicalConditions::decrement_by(const GrandCanonicalConditions &cond_decrement) {
    m_temperature -= cond_decrement.m_temperature;
    for(int i = 0; i < m_mu.size(); i++) {
      m_mu(i) -= cond_decrement.m_mu(i);
    }

    m_beta = 1.0 / (CASM::KB * m_temperature);
    return;
  }


  // ***************************************OPERATORS********************************************** //

  GrandCanonicalConditions &GrandCanonicalConditions::operator+=(const GrandCanonicalConditions &RHS) {
    increment_by(RHS);
    return *this;
  }

  GrandCanonicalConditions GrandCanonicalConditions::operator+(const GrandCanonicalConditions &RHS) const {
    return GrandCanonicalConditions(*this) += RHS;
  }

  ///Subtract temperature and all chemical potentials to *this
  GrandCanonicalConditions &GrandCanonicalConditions::operator-=(const GrandCanonicalConditions &RHS) {
    decrement_by(RHS);
    return *this;
  }

  GrandCanonicalConditions GrandCanonicalConditions::operator-(const GrandCanonicalConditions &RHS) const {
    return GrandCanonicalConditions(*this) -= RHS;
  }

  bool GrandCanonicalConditions::operator==(const GrandCanonicalConditions &RHS) const {
    if(!almost_zero(m_temperature - RHS.m_temperature, m_tolerance)) {
      return false;
    }

    for(int i = 0; i < m_mu.size(); i++) {
      if(!almost_zero(m_mu(i) - RHS.m_mu(i), m_tolerance)) {
        return false;
      }
    }

    return true;
  }

  bool GrandCanonicalConditions::operator!=(const GrandCanonicalConditions &RHS) const {
    return !(*this == RHS);
  }

  bool GrandCanonicalConditions::operator<(const GrandCanonicalConditions &RHS) const {
    // all parameters <

    //equal is not less than
    if(*this == RHS) {
      return false;
    }

    //check temperature
    if(m_temperature > RHS.m_temperature) {
      return false;
    }

    //check all chemical potentials
    for(int i = 0; i < m_mu.size(); i++) {
      if(m_mu(i) > RHS.m_mu(i)) {
        return false;
      }
    }
    return true;
  };

  bool GrandCanonicalConditions::operator<=(const GrandCanonicalConditions &RHS) const {
    // all parameters <=
    if(*this == RHS || *this < RHS) {
      return true;
    }
    return false;
  }

  bool GrandCanonicalConditions::operator>(const GrandCanonicalConditions &RHS) const {
    // all parameters >

    //equal is not greater
    if(*this == RHS) {
      return false;
    }
    if(m_temperature < RHS.m_temperature) {
      return false;
    }
    for(int i = 0; i < m_mu.size(); i++) {
      if(m_mu(i) < RHS.m_mu(i)) {
        return false;
      }
    }

    return true;
  }

  bool GrandCanonicalConditions::operator>=(const GrandCanonicalConditions &RHS) const {
    // all parameters >=
    if(*this == RHS || *this > RHS) {
      return true;
    }
    return false;
  }

  int GrandCanonicalConditions::operator/(const GrandCanonicalConditions &RHS_inc) const {
    int max_division = 0;

    if(!almost_zero(RHS_inc.temperature())) {
      max_division = round(temperature() / RHS_inc.temperature());
    }

    for(Index i = 0; i < mu().size(); i++) {
      int temp_division = round(mu(i) / RHS_inc.mu(i));

      if(temp_division > max_division && !almost_zero(RHS_inc.mu(i))) {
        max_division = temp_division;
      }
    }

    return max_division;
  }
  
  

}


