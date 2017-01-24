#include "casm/monte_carlo/MonteCarlo.hh"
#include "casm/clex/Configuration.hh"

namespace CASM {


  bool SamplerNameCompare::operator()(const std::string &A, const std::string &B) const {
    std::string::size_type Apos1 = A.find_first_of("([");
    std::string::size_type Bpos1 = B.find_first_of("([");
    if(A.substr(0, Apos1) == B.substr(0, Bpos1)) {
      std::string::size_type Apos2 = A.find_first_of("])");
      std::string::size_type Bpos2 = B.find_first_of("])");

      std::string Aindex = A.substr(Apos1 + 1, Apos2 - Apos1 - 1);
      std::string Bindex = B.substr(Bpos1 + 1, Bpos2 - Bpos1 - 1);

      for(int i = 0; i < Aindex.size(); i++) {
        if(!std::isdigit(Aindex[i])) {
          return Aindex < Bindex;
        }
      }

      return std::stoi(Aindex) < std::stoi(Bindex);
    }
    return A.substr(0, Apos1) < B.substr(0, Bpos1);
  }

  /// \brief Samples all requested property data, and stores pass and step number sample was taken at
  void MonteCarlo::sample_data(const MonteCounter &counter) {

    // call MonteSamper::sample(*this) for all samplers
    for(auto it = m_sampler.begin(); it != m_sampler.end(); ++it) {
      it->second->sample(*this, counter);
    }
    m_sample_time.push_back(std::make_pair(counter.pass(), counter.step()));

    if(m_write_trajectory) {
      m_trajectory.push_back(configdof());
    }

    m_is_equil_uptodate = false;
    m_is_converged_uptodate = false;
  }

  /// \brief Clear all data from all samplers
  void MonteCarlo::clear_samples() {
    for(auto it = m_sampler.begin(); it != m_sampler.end(); ++it) {
      it->second->clear();
    }
    m_trajectory.clear();
    m_sample_time.clear();

    m_is_equil_uptodate = false;
    m_is_converged_uptodate = false;
    m_next_convergence_check = m_convergence_check_period;
  }

  /// \brief Returns pair(true, equil_samples) if required equilibration has occured for all samplers that must converge
  ///
  /// - equil_samples is the number of samples required for all samplers that must equilibrate to equilibrate
  /// - If must_converge() == false, always returns false
  ///
  std::pair<bool, MonteSampler::size_type> MonteCarlo::is_equilibrated() const {

    if(m_is_equil_uptodate) {
      return m_is_equil;
    }

    m_is_equil_uptodate = true;


    _log().check<Log::verbose>("Equilibration");
    _log() << std::boolalpha
           << std::setw(24) << "quantity"
           << std::setw(20) << "is_equilibrated"
           << std::setw(16) << "at_sample"
           << std::endl;
    // check if all samplers that must converge have equilibrated, and find the maximum
    //   of the number of samples needed to equilibrate
    MonteSampler::size_type max_equil_samples = 0;
    for(auto it = m_sampler.cbegin(); it != m_sampler.cend(); ++it) {
      if(it->second->must_converge()) {

        // is_equilibrated returns std::pair(is_equilibrated?, equil_samples)
        auto equil = it->second->is_equilibrated();

        _log() << std::setw(24) << it->second->name()
               << std::setw(20) << (equil.first ? "true" : "false"); // why isn't boolapha working?
        if(equil.first) {
          _log() << std::setw(16) << equil.second << std::endl;
        }
        else {
          _log() << std::setw(16) << "unknown" << std::endl;
        }

        if(!equil.first) {
          return m_is_equil = std::make_pair(false, MonteSampler::size_type(0));
        }
        if(equil.second > max_equil_samples) {
          max_equil_samples = equil.second;
        }
      }
      else {
        //        _log() << std::setw(24) << it->second->name()
        //               << std::setw(20) << "unknown"
        //               << std::setw(16) << "unknown" << std::endl;
      }
    }

    _log() << "Overall equilibration at sample: " << max_equil_samples << "\n" << std::endl;

    return m_is_equil = std::make_pair(true, max_equil_samples);
  }

  /// \brief Check to see if all the properties required to converge have converged
  ///
  /// - Calls is_equilibrated() to check if required equilibration has occured
  /// - If equilibrated, checks if the required convergence has occured for all sampler that must converge
  /// - Convergence is checked only for the range of observations that were taken after all
  ///   samplers that must converge have equilibrated
  ///
  bool MonteCarlo::is_converged() const {

    // if we've already calculated convergence, return result
    if(m_is_converged_uptodate) {
      return m_is_converged;
    }

    m_is_converged_uptodate = true;

    // set the next time for a convergence check
    _set_check_convergence_time();

    if(!m_sampler.size()) {
      m_is_converged = false;
      return m_is_converged;
    }

    // check if required equilibration has occured, and get the number of samples required to equilibrate
    auto equil = is_equilibrated();
    if(!equil.first) {
      m_is_converged = false;
      return m_is_converged;
    }

    _log().check<Log::verbose>("Convergence");
    _log() << std::boolalpha
           << std::setw(24) << "quantity"
           << std::setw(16) << "mean"
           << std::setw(16) << "req_prec"
           << std::setw(16) << "calc_prec"
           << std::setw(16) << "is_converged"
           << std::endl;

    // check if all samplers that must converge have converged using data in range [max_equil_samples, end)
    m_is_converged = std::all_of(m_sampler.cbegin(),
                                 m_sampler.cend(),
    [ = ](const SamplerMap::value_type & val) {
      if(val.second->must_converge()) {
        bool result = val.second->is_converged(equil.second);
        _log() << std::setw(24) << val.second->name()
               << std::setw(16) << val.second->mean(equil.second)
               << std::setw(16) << val.second->requested_precision()
               << std::setw(16) << val.second->calculated_precision(equil.second)
               << std::setw(16) << (result ? "true" : "false") // why isn't boolapha working?
               << std::endl;
        return result;
      }
      else {
        //        _log() << std::setw(24) << val.second->name()
        //               << std::setw(16) << val.second->mean(equil.second)
        //               << std::setw(16) << "none"
        //               << std::setw(16) << "unknown"
        //               << std::setw(16) << "unknown"
        //               << std::endl;
        return true;
      }
    });

    _log() << "Overall convergence?: " << std::boolalpha << m_is_converged << "\n" << std::endl;

    return m_is_converged;
  }

  /// \brief Returns true if a convergence check is due
  ///
  /// Currently set to every 10 samples
  bool MonteCarlo::check_convergence_time() const {

    if(m_sample_time.size() >= m_next_convergence_check) {
      return true;
    }
    return false;
  }

  /// \brief Set the next time convergence is due to be checked
  ///
  /// Currently set to 10 samples from the time this function is called
  void MonteCarlo::_set_check_convergence_time() const {

    m_next_convergence_check = m_sample_time.size() + m_convergence_check_period;

  }

}