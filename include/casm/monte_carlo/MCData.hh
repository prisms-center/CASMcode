#ifndef MCData_HH
#define MCData_HH

#include "casm/CASM_global_definitions.hh"

namespace CASM {

  /// \brief MCData stores observations of properties
  class MCData {

  public:

    typedef unsigned long int size_type;

    /// \brief Default constructor
    MCData() :
      MCData(1) {}

    /// \brief Constructor with initial buffer size 'count'
    MCData(size_type count) :
      m_observation(Eigen::VectorXd::Zero(count)),
      m_size(0) {}

    /// \brief Forget all the observed values (does not resize reserved space)
    void clear() {
      m_size = 0;
    }

    /// \brief Add an observation
    void push_back(double value) {

      m_observation(m_size) = value;
      ++m_size;

      // re-size as necessary by doubling reserved space
      if(m_size == m_observation.size()) {
        Eigen::VectorXd tmp = Eigen::VectorXd::Zero(m_observation.size() * 2);
        tmp.head(m_observation.size()) = m_observation;
        m_observation = tmp;
      }
    }

    /// \brief Return all observations
    Eigen::VectorBlock<const Eigen::VectorXd> observations() const {
      return m_observation.head(m_size);
    }

    /// \brief Number of observations
    size_type size() const {
      return m_size;
    }


  private:

    /// \brief vector of all observations (includes m_size observations, and the rest is reserved space)
    Eigen::VectorXd m_observation;

    /// \brief The number of observations
    size_type m_size;

  };

  /// \brief Checks if a range of observations have equilibrated
  ///
  class MCDataEquilibration {

  public:

    typedef MCData::size_type size_type;


    /// \brief Default constructor
    MCDataEquilibration() {}

    /// \brief Check if a range of observations have equilibrated
    MCDataEquilibration(const Eigen::VectorXd &observations, double prec);

    bool is_equilibrated() const {
      return m_is_equilibrated;
    }

    size_type equilibration_samples() const {
      return m_equil_samples;
    }


  private:

    bool m_is_equilibrated;
    size_type m_equil_samples;

  };

  /// \brief Checks if a range of observations have converged
  ///
  class MCDataConvergence {

  public:

    typedef MCData::size_type size_type;


    /// \brief Default constructor
    MCDataConvergence() {}

    /// \brief Construct a MCDataConvergence object
    MCDataConvergence(const Eigen::VectorXd &observations, double conf);

    /// \brief Returns true if converged to the requested level
    ///
    /// \param prec Desired absolute precision (<X> +/- prec)
    ///
    /// \returns true if Var(<X>) <= pow(prec/(sqrt(2.0)*inv_erf(1.0-conf)), 2.0)
    ///
    /// \seealso MCDataConvergence(const Eigen::VectorXd& observations, double conf)
    ///
    bool is_converged(double prec) const {
      return m_calculated_prec <= prec;
    }

    /// \brief <X>
    double mean() const {
      return m_mean;
    }

    /// \brief <X*X>
    double squared_norm() const {
      return m_squared_norm;
    }

    /// \brief Calculated precision of <X>
    ///
    /// \returns sqrt(2.0*Var(<X>))*inv_erf(1.0-conf)
    ///
    /// - This is the critical value is used for comparison in is_converged
    /// - Var(<X>), the variance in the mean of X, is calculated as described by MCDataConvergence
    ///
    double calculated_precision() const {
      return m_calculated_prec;
    }

  private:

    /// \brief Try to find rho = pow(2.0, -1.0/i), using min i such that CoVar[i]/CoVar[0] <= 0.5
    std::tuple<bool, double, double> _calc_rho(const Eigen::VectorXd &observations);

    bool m_is_converged;
    double m_mean;
    double m_squared_norm;
    double m_calculated_prec;

  };

}

#endif
