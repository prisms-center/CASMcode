#include "casm/monte_carlo/MonteCarlo.hh"
#include "casm/clex/Configuration.hh"

namespace CASM {

  /// \brief Samples all requested property data, and stores pass and step number sample was taken at
  void MonteCarlo::sample_data(MonteCounter::size_type pass, MonteCounter::size_type step) {
    // call MonteSamper::sample(*this) for all samplers
    for(auto it = m_sampler.begin(); it != m_sampler.end(); ++it) {
      it->second->sample(*this);
    }
    m_sample_time.push_back(std::make_pair(pass, step));

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

    // check if all samplers that must converge have equilibrated, and find the maximum
    //   of the number of samples needed to equilibrate
    MonteSampler::size_type max_equil_samples = 0;
    for(auto it = m_sampler.cbegin(); it != m_sampler.cend(); ++it) {
      if(it->second->must_converge()) {

        // is_equilibrated returns std::pair(is_equilibrated?, equil_samples)
        auto equil = it->second->is_equilibrated();
        if(!equil.first) {
          return m_is_equil = std::make_pair(false, MonteSampler::size_type(0));
        }
        if(equil.second > max_equil_samples) {
          max_equil_samples = equil.second;
        }
      }
    }

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

    // check if all samplers that must converge have converged using data in range [max_equil_samples, end)
    m_is_converged = std::all_of(m_sampler.cbegin(),
                                 m_sampler.cend(),
    [ = ](const SamplerMap::value_type & val) {
      if(val.second->must_converge()) {
        return val.second->is_converged(equil.second);
      }
      else {
        return true;
      }
    });
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

  /// \brief Fill supercell with motif, applying a factor group operation if necessary
  ConfigDoF fill_supercell(Supercell &mc_scel, const Configuration &motif) {

    const Lattice &motif_lat = motif.get_supercell().get_real_super_lattice();
    const Lattice &scel_lat = mc_scel.get_real_super_lattice();
    auto begin = mc_scel.get_primclex().get_prim().factor_group().begin();
    auto end = mc_scel.get_primclex().get_prim().factor_group().end();

    auto res = is_supercell(scel_lat, motif_lat, begin, end, TOL);
    if(res.first == end) {

      std::cerr << "Requested supercell transformation matrix: \n"
                << mc_scel.get_transf_mat() << "\n";
      std::cerr << "Requested motif Configuration: " <<
                motif.name() << "\n";
      std::cerr << "Configuration transformation matrix: \n"
                << motif.get_supercell().get_transf_mat() << "\n";

      throw std::runtime_error(
        "Error in 'fill_supercell(const Supercell &mc_scel, const Configuration& motif)'\n"
        "  The motif cannot be tiled onto the specified supercell."
      );
    }

    ConfigTransform f(mc_scel, *res.first);
    return copy_apply(f, motif).configdof();
  }

}