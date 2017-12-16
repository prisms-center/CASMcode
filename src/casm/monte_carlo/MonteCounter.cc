#include "casm/monte_carlo/MonteCounter.hh"

#include "casm/monte_carlo/MonteSettings.hh"
#include "casm/monte_carlo/MonteCarlo.hh"

namespace CASM {

  // ---- MonteCounter Definitions ---------------------------------

  /// \brief Construct a MonteCounter and initialize to all counts to zero
  MonteCounter::MonteCounter(const EquilibriumMonteSettings &settings, size_type _steps_per_pass) :
    m_steps_per_pass(_steps_per_pass) {

    // --- If not in 'must converge mode' -----

    if(((int) settings.is_N_step()) + ((int) settings.is_N_pass()) + ((int) settings.is_N_sample()) > 1) {
      throw std::runtime_error(
        std::string("Error in MonteCounter constructor\n") +
        "  Zero or one of 'N_step', 'N_pass', and 'N_sample' should be given");
    }

    /// Requested number of steps
    if(settings.is_N_step()) {
      m_is_N_step = settings.is_N_step();
      m_N_step = settings.N_step();
    }

    /// Requested number of passes
    if(settings.is_N_pass()) {
      m_is_N_pass = settings.is_N_pass();
      m_N_pass = settings.N_pass();
    }

    /// Requested number of samples
    if(settings.is_N_sample()) {
      m_is_N_sample = settings.is_N_sample();
      m_N_sample = settings.N_sample();
    }


    // --- Imposed limits ------------

    /// Maximum allowed number of steps
    if(settings.is_max_step()) {
      m_is_max_step = true;
      m_max_step = settings.max_step();
    }

    /// Minimum allowed number of steps
    if(settings.is_min_step()) {
      m_is_min_step = true;
      m_min_step = settings.min_step();
    }

    /// Maximum allowed number of passes
    if(settings.is_max_pass()) {
      m_is_max_pass = true;
      m_max_pass = settings.max_pass();
    }

    /// Minimum allowed number of passes
    if(settings.is_min_pass()) {
      m_is_min_pass = true;
      m_min_pass = settings.min_pass();
    }

    /// Maximum allowed number of samples
    if(settings.is_max_sample()) {
      m_is_max_sample = true;
      m_max_sample = settings.max_sample();
    }

    /// Minimum allowed number of samples
    if(settings.is_min_sample()) {
      m_is_min_sample = true;
      m_min_sample = settings.min_sample();
    }


    // --- Periods at which certain things should happen --------

    /// Whether sampling mode is by Monte::SAMPLE_MODE::PASS or Monte::SAMPLE_MODE::STEP
    m_sample_mode = settings.sample_by_pass() ? Monte::SAMPLE_MODE::PASS : Monte::SAMPLE_MODE::STEP;

    /// How often to sample data. In terms of m_sample_mode. If == 1 sample every pass/step.
    m_sample_period = settings.sample_period();

  }


  MonteCounter::size_type MonteCounter::pass() const {
    return m_pass;
  }

  MonteCounter::size_type MonteCounter::step() const {
    return m_step;
  }


  /// \brief Number of steps per pass
  MonteCounter::size_type MonteCounter::steps_per_pass() const {
    return m_steps_per_pass;
  }


  MonteCounter::size_type MonteCounter::samples() const {
    return m_samples;
  }


  void MonteCounter::increment_samples() {
    m_samples++;
    m_since_last_sample = 0;
  }


  void MonteCounter::reset() {
    m_pass = 0;
    m_step = 0;

    m_samples = 0;

    m_since_last_sample = 0;

    return;
  }

  /// \brief Set all counter variables for performing a restart
  void MonteCounter::set(size_type _pass, size_type _step, size_type _samples) {

    reset();

    m_pass = _pass;
    m_step = _step;

    m_samples = _samples;

  }


  bool MonteCounter::sample_time() const {
    if(m_since_last_sample == m_sample_period) {
      return true;
    }
    return false;
  }

  /// \brief Prefix increment step and updates pass
  MonteCounter &MonteCounter::operator++() {

    ++m_step;
    if(m_step == m_steps_per_pass) {
      ++m_pass;
      m_step = 0;

      if(m_sample_mode == Monte::SAMPLE_MODE::PASS) {
        ++m_since_last_sample;
      }
    }

    if(m_sample_mode == Monte::SAMPLE_MODE::STEP) {
      ++m_since_last_sample;
    }

    return *this;
  }

  /// \brief Postfix increment step and updates pass
  MonteCounter MonteCounter::operator++(int) {

    MonteCounter init(*this);

    ++(*this);

    return init;
  }


  /// \brief Check if requested number of pass, step, or samples has been met
  bool MonteCounter::is_complete() const {
    if(m_is_N_step && step() >= m_N_step) {
      return true;
    }
    if(m_is_N_pass && pass() >= m_N_pass) {
      return true;
    }
    if(m_is_N_sample && samples() >= m_N_sample) {
      return true;
    }
    return false;
  }

  /// \brief Check if minimum number of pass, step, and samples has been met
  bool MonteCounter::minimums_met() const {

    if(m_is_min_sample && (m_samples < m_min_sample)) {
      return false;
    }

    if(m_is_min_pass && (m_pass < m_min_pass)) {
      return false;
    }

    if(m_is_min_step && (m_pass * m_steps_per_pass + m_step < m_min_step)) {
      return false;
    }

    return true;
  }

  /// \brief Check if maximum number of pass, step, and samples has been met
  bool MonteCounter::maximums_met() const {

    if(m_is_max_sample && (m_samples >= m_max_sample)) {
      return true;
    }

    if(m_is_max_pass && (m_pass >= m_max_pass)) {
      return true;
    }

    if(m_is_max_step && (m_pass * m_steps_per_pass + m_step >= m_max_step)) {
      return true;
    }

    return false;
  }

  void MonteCounter::debugprint(std::ostream &sout) const {
    sout << "pass: " << m_pass << "\n";
    sout << "step: " << m_step << "\n";
    sout << "samples: " << m_samples << "\n";
    sout << "since_last_sample: " << m_since_last_sample << "\n";
    sout << "is_max_pass: " << m_is_max_pass << "\n";
    sout << "m_max_pass: " << m_max_pass << "\n\n";

    sout << "is_min_pass: " << m_is_min_pass << "\n";
    sout << "m_min_pass: " << m_min_pass << "\n\n";

    sout << "is_max_step: " << m_is_max_step << "\n";
    sout << "m_max_step: " << m_max_step << "\n\n";

    sout << "is_min_step: " << m_is_min_step << "\n";
    sout << "m_min_step: " << m_min_step << "\n\n";

    sout << "is_max_sample: " << m_is_max_sample << "\n";
    sout << "m_max_sample: " << m_max_sample << "\n\n";

    sout << "is_min_sample: " << m_is_min_sample << "\n" << std::endl;
    sout << "m_min_sample: " << m_min_sample << "\n\n" << std::endl;

  }


}

