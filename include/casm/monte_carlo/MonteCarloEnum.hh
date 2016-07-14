#ifndef CASM_MonteCarloEnum
#define CASM_MonteCarloEnum

#include <functional>
#include <utility>

#include "casm/misc/HallOfFame.hh"
#include "casm/misc/cloneable_ptr.hh"

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

  struct MonteCarloEnumCheck {

  public:

    MonteCarloEnumCheck(const DataFormatter<Configuration> &_formatter) :
      m_formatter(_formatter) {}

    bool operator()(const Configuration &config) {
      return m_formatter.evaluate_as_scalar<bool>(config);
    }

  private:

    DataFormatter<Configuration> m_formatter;
  };

  template<typename Base>
  class MonteCarloEnum : public Base {

  public:

    typedef HallOfFame<Configuration, MonteCarloEnumMetric, ConfigDoFOccCompare> HallOfFameType;

    template<typename MonteTypeSettings>
    MonteCarloEnum(PrimClex &primclex, const MonteTypeSettings &settings, Log &_log);

    /// \brief Check if enumeration is requested
    bool enum_configs() const;

    /// \brief Attempt to insert (canonical) Configuration into enumeration hall of fame
    HallOfFameType::InsertResult enum_insert(const Configuration &config);

    /// \brief const Access the enumeration hall of fame
    const HallOfFameType &enum_halloffame() const;

    /// \brief Save configurations in the hall of fame to the config list
    void save_configs() const;

    std::string enum_check_args() const {
      return m_check_args;
    }

    std::string enum_metric_args() const {
      return m_metric_args;
    }

    void print_info() const;

  protected:

    HallOfFameType &_enum_halloffame();


  private:

    std::string m_check_args;
    std::string m_metric_args;

    /// \brief Use for enumerating configurations via Monte Carlo
    mutable notstd::cloneable_ptr<HallOfFameType> m_halloffame;

    /// \brief Use for enumerating configurations via Monte Carlo
    notstd::cloneable_ptr<MonteCarloEnumCheck> m_enum_check;

  };
}

#endif
