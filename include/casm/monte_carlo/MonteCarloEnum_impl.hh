#ifndef CASM_MonteCarloEnum_impl
#define CASM_MonteCarloEnum_impl

#include "casm/monte_carlo/MonteCarloEnum.hh"
#include "casm/clex/ConfigIterator.hh"
#include "casm/clex/ConfigMapping.hh"
#include "casm/clex/PrimClex.hh"

namespace CASM {

  template<typename MonteTypeSettings, typename MonteCarloType>
  MonteCarloEnum::MonteCarloEnum(PrimClex &primclex, const MonteTypeSettings &set, Log &log, MonteCarloType &mc) :
    m_primclex(primclex),
    m_log(log),
    m_sample_mode(set.enumeration_sample_mode()),
    m_debug(set.debug()),
    m_check_args(set.enumeration_check_args()),
    m_metric_args(set.enumeration_metric_args()),
    m_check_existence(set.enumeration_check_existence()),
    m_insert_canonical(set.enumeration_insert_canonical()),
    m_dict(primclex.settings().query_handler<Configuration>().dict()) {

    m_dict.insert(
      ConfigIO::GenericConfigFormatter<double>(
        "potential_energy",
        "potential_energy",
    [&](const Configuration & config) {
      return mc.potential_energy(config);
    },
    [&](const Configuration & config) {
      return true;
    }));

    m_dict.insert(
      ConfigIO::GenericConfigFormatter<bool>(
        "is_new",
        "is_new",
    [&](const Configuration & config) {
      return m_data[config.name()].first;
    },
    [&](const Configuration & config) {
      return m_data.find(config.name()) != m_data.end();
    }
      ));

    m_dict.insert(
      ConfigIO::GenericConfigFormatter<double>(
        "score",
        "score",
    [&](const Configuration & config) {
      return m_data[config.name()].second;
    },
    [&](const Configuration & config) {
      return m_data.find(config.name()) != m_data.end();
    }
      ));

    m_halloffame.unique().reset(
      new HallOfFameType(
        MonteCarloEnumMetric(m_dict.parse(metric_args())),
        std::less<Configuration>(),
        set.enumeration_N_halloffame(),
        set.enumeration_tol()));

    m_enum_check.unique().reset(new MonteCarloEnumCheck(m_dict.parse(set.enumeration_check_args())));

    reset();

    _log().custom("Configuration enumeration");
    _log() << "  check: " << check_args() << "\n";
    _log() << "  metric: " << metric_args() << "\n";
    _log() << "  check_existence: " << std::boolalpha << m_check_existence << "\n";
    _log() << "  insert_canonical: " << std::boolalpha << m_insert_canonical << "\n";
    _log() << "  sample_mode: " << m_sample_mode << "\n";
    _log() << std::endl;

  }

}

#endif

