#ifndef CASM_GrandCanonicalConditions_HH
#define CASM_GrandCanonicalConditions_HH

#include "casm/clex/PrimClex.hh"

namespace CASM {
  
  class MonteSettings;
  
  /// Conditions for a Grand Canonical run:
  /// Temperature
  /// Chemical potential
  /// Tolerance (for comparing conditions)
  /// 
  class GrandCanonicalConditions {
  public:
    
    /// \brief Default constructor
    GrandCanonicalConditions() {}
    
    /// \brief Construct from MonteSettings
    GrandCanonicalConditions(const PrimClex &primclex, const MonteSettings &settings);

    // ***************************************ACCESSORS********************************************** //

    double temperature() const;

    double beta() const;

    /// \brief 'atomic' mu
    const Eigen::VectorXd& mu() const;

    /// \brief 'atomic' mu
    double mu(Index mu_index) const;
    
    /// \brief parametric composition mu
    Eigen::VectorXd param_mu() const;


    double tolerance() const;

    // ***************************************MUTATORS*********************************************** //

    ///Set the temperature of the current grand canonical condition. Sets m_temperature and m_beta.
    void set_temperature(double in_temp);

    ///Set all the chemical potentials of the current grand canonical condition. Sets array m_mu.
    void set_mu(const Eigen::VectorXd &in_mu);

    ///Set a single 'atomic' chemical potential by specifying an index and a value. Sets one value in m_mu.
    void set_mu(Index ind, double in_mu);

    ///Acts as +=
    void increment_by(const GrandCanonicalConditions &cond_increment);

    ///Acts as -=
    void decrement_by(const GrandCanonicalConditions &cond_decrement);

    // ***************************************OPERATORS********************************************** //

    ///Assignment operator
    //GrandCanonicalConditions &operator=(const GrandCanonicalConditions &RHS);

    ///Add temperature and all chemical potentials to *this
    GrandCanonicalConditions &operator+=(const GrandCanonicalConditions &RHS);

    ///Add temperature and all chemical potentials together and return a new Condition
    GrandCanonicalConditions operator+(const GrandCanonicalConditions &RHS) const;

    ///Subtract temperature and all chemical potentials to *this
    GrandCanonicalConditions &operator-=(const GrandCanonicalConditions &RHS);

    ///Subtract temperature and all chemical potentials together and return a new Condition
    GrandCanonicalConditions operator-(const GrandCanonicalConditions &RHS) const;

    ///Compare temperature and all chemical potentials to *this
    bool operator==(const GrandCanonicalConditions &RHS) const;

    ///Compare temperature and all chemical potentials to *this
    bool operator!=(const GrandCanonicalConditions &RHS) const;

    ///Returns true if ALL parameters satisfy the inequality
    bool operator<(const GrandCanonicalConditions &RHS) const;

    ///Returns true if ALL parameters satisfy the inequality
    bool operator<=(const GrandCanonicalConditions &RHS) const;

    ///Returns true if ALL parameters satisfy the inequality
    bool operator>(const GrandCanonicalConditions &RHS) const;

    ///Returns true if ALL parameters satisfy the inequality
    bool operator>=(const GrandCanonicalConditions &RHS) const;

    ///Divide ALL parameters and return the greatest number in absolute value
    int operator/(const GrandCanonicalConditions &RHS_inc) const;

  
  protected:

    
    /// Convert between num atoms per prim cell and parametric composition
    CompositionConverter m_comp_converter;
    
    ///Temperature
    double m_temperature;

    ///Inverse temperature. Includes Boltzmann term
    double m_beta;

    ///Vector of the 'atomic' chemical potentials for each species. Ordered as primclex.get_param_comp().get_components()
    Eigen::VectorXd m_mu;

    ///Tolerance for comparison operators == and !=
    double m_tolerance;
    
    
  };

}

#endif

