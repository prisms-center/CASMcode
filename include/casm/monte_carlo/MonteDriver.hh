#ifndef CASM_MonteDriver_HH
#define CASM_MonteDriver_HH

#include <string>
#include <vector>
#include "casm/global/definitions.hh"
#include "casm/monte_carlo/MonteIO.hh"

namespace CASM {
  class PrimClex;
  class Log;

  namespace Monte {
    class MonteCarloEnum;

    /**
     * MonteDriver consists of a specialized MonteCarlo object and a list of conditions
     * to equilibrate at. The first condition in the list requires a starting configuration
     * (read from the setting), while subsequent conditions are calculated using the
     * final state of the previous condition.
     *
     * The different kinds of drive modes the user can specify are:
     * INCREMENTAL:   Given a delta in condition values, increment the conditions by the delta after each point
     * CUSTOM:        Calculate for a list of condition values
     */

    template<typename RunType>
    class MonteDriver {

    public:
      typedef typename RunType::CondType CondType;
      typedef typename RunType::SettingsType SettingsType;

      /// \brief Constructor via MonteSettings
      MonteDriver(const PrimClex &primclex, const SettingsType &settings, Log &_log, Log &_err_log);

      /// \brief Run everything requested by the MonteSettings
      void run();

    private:

      /// run in debug mode?
      bool debug() const {
        return m_debug;
      }

      ///Return the appropriate std::vector of conditions to visit based from settings. Use for construction.
      std::vector<CondType> make_conditions_list(const PrimClex &primclex, const SettingsType &settings);

      ///Converge the MonteCarlo for conditions 'cond_index'
      void single_run(Index cond_index);

      ///Check for existing calculations to find starting conditions
      Index _find_starting_conditions() const;


      /// target for log messages
      Log &m_log;

      /// target for error messages
      Log &m_err_log;

      ///Copy of initial settings given at construction. Will expand to have MonteCarlo states dumped into it.
      SettingsType m_settings;

      /// describes where to write output
      MonteCarloDirectoryStructure m_dir;

      ///Specifies how to build the conditions list from the settings
      const Monte::DRIVE_MODE m_drive_mode;

      ///Specialized Monte Carlo object to use throughout
      RunType m_mc;

      ///List of specialized conditions to visit in sequential order. Does not include initial conditions.
      const std::vector<CondType> m_conditions_list;

      /// run in debug mode?
      bool m_debug;

      /// Enumerated configurations encountered during Monte Carlo calculations
      notstd::cloneable_ptr<MonteCarloEnum> m_enum;
    };


    /// Perform a single monte carlo step, return true if accepted
    template<typename RunType>
    bool monte_carlo_step(RunType &monte_run);

  }
}

#endif
