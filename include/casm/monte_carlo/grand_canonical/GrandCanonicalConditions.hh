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
    
    /// \brief Constructor
    ///
    /// \param _temperature in K
    /// \param _param_chem_pot Parametric composition chemical potential
    /// \param _comp_converter CompositionConverter for converting from parametric chem_pot to species chem_pot
    /// \param _tol tolerance for comparing conditions
    ///
    GrandCanonicalConditions(double _temperature, 
                             const Eigen::VectorXd &_param_chem_pot,
                             const CompositionConverter& _comp_converter,
                             double _tol);
    
    // ***************************************ACCESSORS********************************************** //

    double temperature() const;

    double beta() const;

    /// \brief chemical potential: dg/dcomp_n
    const Eigen::VectorXd& chem_pot() const;

    /// \brief chemical potential: dg/dcomp_n(index)
    double chem_pot(Index index) const;
    
    /// \brief exchange chemical potential: chem_pot(new) - chem_pot(curr)
    double exchange_chem_pot(Index index_new, Index index_curr) const;

    /// \brief parametric chemical potential: dg/dcomp_x
    Eigen::VectorXd param_chem_pot() const;
    
    /// \brief parametric chemical potential: dg/dcomp_x(index)
    double param_chem_pot(Index index) const;
    

    double tolerance() const;

    // ***************************************MUTATORS*********************************************** //

    ///Set the temperature of the current grand canonical condition.
    void set_temperature(double in_temp);

    ///Set all the chemical potentials of the current grand canonical condition.
    void set_chem_pot(const Eigen::VectorXd &in_chem_pot);

    ///Set a single 'atomic' chemical potential by specifying an index and a value.
    void set_chem_pot(Index ind, double in_chem_pot);
    
    ///Set all the parametric chemical potentials of the current grand canonical condition.
    void set_param_chem_pot(const Eigen::VectorXd &in_chem_pot);

    ///Set a single parametric chemical potential by specifying an index and a value.
    void set_param_chem_pot(Index ind, double in_chem_pot);

    
    // ***************************************OPERATORS********************************************** //

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
    Eigen::VectorXd m_chem_pot;
    
    ///Vector of the parametric chemical potentials conjugate to the parametric compositions.
    Eigen::VectorXd m_param_chem_pot;

    ///Tolerance for comparison operators == and !=
    double m_tolerance;
    
    
  };

}

#endif

