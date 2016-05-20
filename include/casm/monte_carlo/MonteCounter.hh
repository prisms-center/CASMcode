#ifndef CASM_MonteCounter_HH
#define CASM_MonteCounter_HH

#include "casm/CASM_global_definitions.hh"
#include "casm/monte_carlo/MonteDefinitions.hh"

namespace CASM {

  class EquilibriumMonteSettings;
  class MonteCarlo;

  /// \brief Track the number of passes, steps and samples taken in a Monte Carlo calculation.
  ///
  /// It's a convenience for figuring out which pass you're on by only having to think of
  /// incrementing after every step.
  ///
  /// It also has information about the period at which certain routines should be called,
  /// such as how often to take a snapshot or write out the state of the run to a file.
  ///
  class MonteCounter {
  public:

    typedef Index size_type;


    /// \brief Construct a MonteCounter and initialize to all counts to zero
    MonteCounter(const EquilibriumMonteSettings &settings, size_type _steps_per_pass);


    /// \brief Number of complete passes performed
    size_type pass() const;

    /// \brief Number of steps into the current pass
    size_type step() const;


    /// \brief Number of steps per pass
    size_type steps_per_pass() const;


    /// \brief Number of samples taken
    size_type samples() const;


    /// \brief Increments the number of samples taken and resets counter until the next sample should be taken
    void increment_samples();


    /// \brief Set all counter variables back to 0
    void reset();

    /// \brief Set all counter variables for performing a restart
    void set(size_type _pass,
             size_type _step,
             size_type _samples);


    /// \brief Returns true if based on period and current number of steps it is time to take a sample
    bool sample_time() const;


    /// \brief Prefix increment step and updates pass
    MonteCounter &operator++();

    /// \brief Postfix increment step and updates pass
    MonteCounter operator++(int);


    /// \brief Check if requested number of pass, step, or samples has been met
    bool is_complete() const;

    /// \brief Check if minimum number of pass, step, and samples has been met
    bool minimums_met() const;

    /// \brief Check if maximum number of pass, step, and samples has been met
    bool maximums_met() const;

    void debugprint(std::ostream &sout) const;

  private:

    /// One pass is equal to the number of sites with variable degrees of freedom in the supercell
    size_type m_pass = 0;

    /// Number of steps taken during the current pass
    size_type m_step = 0;


    /// Number of data samples taken.
    size_type m_samples = 0;


    /// Number of passes or steps since last sample
    size_type m_since_last_sample = 0;

    // --- If not in 'must converge mode' -----

    /// Requested number of steps
    bool m_is_N_step = false;
    size_type m_N_step = 0;

    /// Requested number of passes
    bool m_is_N_pass = false;
    size_type m_N_pass = 0;

    /// Requested number of samples
    bool m_is_N_sample = false;
    size_type m_N_sample = 0;



    // --- Imposed limits ------------

    /// Maximum allowed number of steps
    bool m_is_max_step = false;
    size_type m_max_step = 0;

    /// Minimum allowed number of steps
    bool m_is_min_step = false;
    size_type m_min_step = 0;

    /// Maximum allowed number of passes
    bool m_is_max_pass = false;
    size_type m_max_pass = 0;

    /// Minimum allowed number of passes
    bool m_is_min_pass = false;
    size_type m_min_pass = 0;

    /// Maximum allowed number of samples to take
    bool m_is_max_sample = false;
    size_type m_max_sample = 0;

    /// Minimum allowed number of samples to take
    bool m_is_min_sample = false;
    size_type m_min_sample = 0;


    // --- Periods at which certain things should happen --------

    /// Whether sampling mode is by Monte::SAMPLE_MODE::PASS or Monte::SAMPLE_MODE::STEP
    Monte::SAMPLE_MODE m_sample_mode = Monte::SAMPLE_MODE::PASS;

    /// How often to sample data. In terms of m_sample_mode. If == 1 sample every pass/step.
    int m_sample_period;

    /// Number of steps per pass is equal to the number of sites with variable degrees of freedom in the supercell
    int m_steps_per_pass;

  };

}

#endif


