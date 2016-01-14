#include "casm/monte_carlo/grand_canonical/GrandCanonicalConditions.hh"
#include "casm/monte_carlo/MonteCarlo.hh"
#include "casm/monte_carlo/MonteSettings.hh"

#include "casm/clex/PrimClex.hh"

namespace CASM {
  
  
  /// \brief Constructor
  ///
  /// \param _temperature in K
  /// \param _param_chem_pot Parametric composition chemical potential
  /// \param _comp_converter CompositionConverter for converting from parametric chem_pot to atomic chem_pot
  /// \param _tol tolerance for comparing conditions
  ///
  GrandCanonicalConditions::GrandCanonicalConditions(double _temperature, 
                                                     const Eigen::VectorXd &_param_chem_pot,
                                                     const CompositionConverter& _comp_converter,
                                                     double _tol) : 
    m_comp_converter(_comp_converter),
    m_tolerance(_tol) {
    
    // -- set T ----
    set_temperature(_temperature);

    
    // -- set chem_pot ----
    set_chem_pot(m_comp_converter.chem_pot(_param_chem_pot));

  }
  
  // ***************************************ACCESSORS********************************************** //

  double GrandCanonicalConditions::temperature() const {
    return m_temperature;
  }

  double GrandCanonicalConditions::beta() const {
    return m_beta;
  }

  const Eigen::VectorXd& GrandCanonicalConditions::chem_pot() const {
    return m_chem_pot;
  }

  double GrandCanonicalConditions::chem_pot(Index index) const {
    return m_chem_pot(index);
  }
  
  Eigen::VectorXd GrandCanonicalConditions::param_chem_pot() const {
    return m_comp_converter.param_chem_pot(m_chem_pot);
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

  void GrandCanonicalConditions::set_chem_pot(const Eigen::VectorXd &in_chem_pot) {
    m_chem_pot = in_chem_pot;
    return;
  }

  void GrandCanonicalConditions::set_chem_pot(Index ind, double in_chem_pot) {
    m_chem_pot(ind) = in_chem_pot;
    return;
  }

  void GrandCanonicalConditions::increment_by(const GrandCanonicalConditions &cond_increment) {
    m_temperature += cond_increment.m_temperature;
    for(int i = 0; i < m_chem_pot.size(); i++) {
      m_chem_pot(i) += cond_increment.m_chem_pot(i);
    }
    m_beta = 1.0 / (CASM::KB * m_temperature);
    return;
  }

  void GrandCanonicalConditions::decrement_by(const GrandCanonicalConditions &cond_decrement) {
    m_temperature -= cond_decrement.m_temperature;
    for(int i = 0; i < m_chem_pot.size(); i++) {
      m_chem_pot(i) -= cond_decrement.m_chem_pot(i);
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

    for(int i = 0; i < m_chem_pot.size(); i++) {
      if(!almost_zero(m_chem_pot(i) - RHS.m_chem_pot(i), m_tolerance)) {
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
    for(int i = 0; i < m_chem_pot.size(); i++) {
      if(m_chem_pot(i) > RHS.m_chem_pot(i)) {
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
    for(int i = 0; i < m_chem_pot.size(); i++) {
      if(m_chem_pot(i) < RHS.m_chem_pot(i)) {
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

    for(Index i = 0; i < chem_pot().size(); i++) {
      int temp_division = round(chem_pot(i) / RHS_inc.chem_pot(i));

      if(temp_division > max_division && !almost_zero(RHS_inc.chem_pot(i))) {
        max_division = temp_division;
      }
    }

    return max_division;
  }
  
  

}


