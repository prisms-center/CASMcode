#ifndef CASM_MonteCarlo_impl
#define CASM_MonteCarlo_impl

#include "casm/monte_carlo/MonteCarlo.hh"
#include "casm/monte_carlo/grand_canonical/GrandCanonicalSettings_impl.hh"
#include "casm/monte_carlo/canonical/CanonicalSettings_impl.hh"
#include "casm/casm_io/Log.hh"
#include "casm/app/ProjectSettings.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/Supercell.hh"
#include "casm/clex/ConfigDoF.hh"
#include "casm/clex/Configuration.hh"

namespace CASM {
  namespace Monte {

    /// \brief Construct with a starting ConfigDoF as specified the given Settings and prepare data samplers
    template<typename MonteTypeSettings>
    MonteCarlo::MonteCarlo(const PrimClex &primclex, const MonteTypeSettings &settings, Log &_log) :
      m_settings(settings),
      m_primclex(primclex),
      m_scel(&primclex, settings.simulation_cell_matrix()),
      m_config(m_scel),
      m_configdof(m_config.configdof()),
      m_write_trajectory(settings.write_trajectory()),
      m_debug(m_settings.debug()),
      m_log(_log) {

      settings.samplers(primclex, std::inserter(m_sampler, m_sampler.begin()));

      m_must_converge = false;
      for(auto it = m_sampler.cbegin(); it != m_sampler.cend(); ++it) {
        if(it->second->must_converge()) {
          m_must_converge = true;
          break;
        }
      }
    }

  }
}

#endif
