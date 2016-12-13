#ifndef CASM_MonteCarloEnum
#define CASM_MonteCarloEnum

#include <functional>
#include <utility>

#include "casm/misc/HallOfFame.hh"
#include "casm/misc/cloneable_ptr.hh"
#include "casm/casm_io/DataFormatter.hh"
#include "casm/casm_io/Log.hh"
#include "casm/clex/Configuration.hh"
#include "casm/monte_carlo/MonteDefinitions.hh"

namespace CASM {

  class Configuration;


  class MonteCarloEnumMetric {

  public:

    MonteCarloEnumMetric(const DataFormatter<Configuration> &_formatter) :
      m_formatter(_formatter) {}

    double operator()(const Configuration &config) {
      return m_formatter.evaluate_as_scalar<double>(config);
    }

  private:

    DataFormatter<Configuration> m_formatter;
  };


  class MonteCarloEnumCheck {

  public:

    MonteCarloEnumCheck(const DataFormatter<Configuration> &_formatter) :
      m_formatter(_formatter) {}

    bool operator()(const Configuration &config) {
      return m_formatter.evaluate_as_scalar<bool>(config);
    }

  private:

    DataFormatter<Configuration> m_formatter;
  };


  class MonteCarloEnum {

  public:

    typedef HallOfFame<Configuration, MonteCarloEnumMetric> HallOfFameType;

    template<typename MonteTypeSettings, typename MonteCarloType>
    MonteCarloEnum(PrimClex &primclex, const MonteTypeSettings &settings, Log &log, MonteCarloType &mc);

    /// \brief const Access the PrimClex that *this is based on
    const PrimClex &primclex() const {
      return m_primclex;
    }

    /// \brief Check if enumeration is requested after every acceptance
    bool on_accept() const {
      return m_sample_mode == Monte::ENUM_SAMPLE_MODE::ON_ACCEPT;
    }

    /// \brief Check if enumeration is requested after every sample
    bool on_sample() const {
      return m_sample_mode == Monte::ENUM_SAMPLE_MODE::ON_SAMPLE;
    }

    /// \brief Attempt to insert (canonical) Configuration into enumeration hall of fame
    HallOfFameType::InsertResult insert(const Configuration &config);

    /// \brief Clear hall of fame
    void clear() {
      if(m_halloffame) {
        m_halloffame->clear();
      }
    }

    /// \brief const Access the enumeration hall of fame
    const HallOfFameType &halloffame() const;

    /// \brief Save configurations in the hall of fame to the config list
    void save_configs();

    std::string check_args() const {
      return m_check_args;
    }

    std::string metric_args() const {
      return m_metric_args;
    }

    void print_info() const;

    /// \brief return true if running in debug mode
    bool debug() const {
      return m_debug;
    }

    /// \brief If true, insert configurations in canonical form.
    ///
    /// If m_check_existence == true, this must be true
    bool check_existence() const {
      return m_check_existence;
    }

    /// \brief Map for faster? access of PrimClex's supercells
    bool insert_canonical() const {
      return m_insert_canonical;
    }

    /// \brief Clear hall of fame and reset excluded
    void reset();

  private:

    /// \brief Insert in hall of fame if 'check' passes
    MonteCarloEnum::HallOfFameType::InsertResult _insert(const Configuration &config);

    HallOfFameType &_halloffame();

    /// \brief Access the PrimClex that *this is based on
    PrimClex &_primclex() const {
      return m_primclex;
    }

    Log &_log() const {
      return m_log;
    }


    /// \brief PrimClex for this system
    PrimClex &m_primclex;

    /// \brief Target for messages
    Log &m_log;

    /// when to attempt to insert configurations in the hall of fame
    Monte::ENUM_SAMPLE_MODE m_sample_mode;

    // in debug mode, allow printing or checking extra things
    bool m_debug;

    std::string m_check_args;
    std::string m_metric_args;

    /// \brief Use for enumerating configurations via Monte Carlo
    mutable notstd::cloneable_ptr<HallOfFameType> m_halloffame;

    /// \brief Use for enumerating configurations via Monte Carlo
    notstd::cloneable_ptr<MonteCarloEnumCheck> m_enum_check;

    /// \brief If true, only keep configurations that are not enumerated already
    bool m_check_existence;

    /// \brief If true, insert configurations in canonical form.
    ///
    /// If m_check_existence == true, this must be true
    bool m_insert_canonical;

    /// \brief Map for faster? access of PrimClex's supercells
    mutable std::map<std::string, Supercell *> m_canon_scel;

    /// \brief Used for various purposes
    DataFormatterDictionary<Configuration> m_dict;

    /// \brief holds 'is_new, score' data
    std::map<std::string, std::pair<bool, double> > m_data;


  };
}

#endif
