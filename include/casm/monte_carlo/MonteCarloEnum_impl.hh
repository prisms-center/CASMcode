#ifndef CASM_MonteCarloEnum_impl
#define CASM_MonteCarloEnum_impl

#include "casm/app/ProjectSettings.hh"
#include "casm/app/QueryHandler.hh"
#include "casm/clex/ConfigIO.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/monte_carlo/MonteCarloEnum.hh"

namespace CASM {
namespace Monte {

template <typename MonteTypeSettings, typename MonteCarloType>
MonteCarloEnum::MonteCarloEnum(const PrimClex &primclex,
                               const MonteTypeSettings &set, Log &log,
                               MonteCarloType &mc)
    : m_primclex(primclex),
      m_log(log),
      m_sample_mode(set.enumeration_sample_mode()),
      m_debug(set.debug()),
      m_check_args(set.enumeration_check_args()),
      m_metric_args(set.enumeration_metric_args()),
      m_check_existence(set.enumeration_check_existence()),
      m_insert_canonical(set.enumeration_insert_canonical()),
      m_insert_primitive_only(set.enumeration_insert_primitive_only()),
      m_save_primitive_only(set.enumeration_save_primitive_only()),
      m_order_parameter(
          mc.order_parameter() == nullptr
              ? std::shared_ptr<OrderParameter>()
              : std::make_shared<OrderParameter>(*mc.order_parameter())) {
  // Configuration query dict
  auto config_query_dict =
      primclex.settings().query_handler<Configuration>().dict();
  config_query_dict.insert(GenericDatumFormatter<double, Configuration>(
      "potential_energy", "potential_energy",
      [&](Configuration const &config) { return mc.potential_energy(config); },
      [&](Configuration const &config) { return true; }));

  config_query_dict.insert(
      Generic1DDatumFormatter<Eigen::VectorXd, Configuration>(
          "corr_matching_error", "corr_matching_error",
          [&](Configuration const &config) {
            return make_corr_matching_error(
                correlations(config, mc.clexulator()),
                *mc.conditions().corr_matching_pot());
          },
          [&](Configuration const &config) {
            return mc.conditions().corr_matching_pot().has_value();
          }));

  config_query_dict.insert(
      Generic1DDatumFormatter<Eigen::VectorXd, Configuration>(
          "random_alloy_corr_matching_error",
          "random_alloy_corr_matching_error",
          [&](Configuration const &config) {
            return make_corr_matching_error(
                correlations(config, mc.clexulator()),
                *mc.conditions().random_alloy_corr_matching_pot());
          },
          [&](Configuration const &config) {
            return mc.conditions().random_alloy_corr_matching_pot().has_value();
          }));

  config_query_dict.insert(
      Generic1DDatumFormatter<Eigen::VectorXd, Configuration>(
          "order_parameter", "order_parameter",
          [&](Configuration const &config) {
            return (*m_order_parameter)(config);
          },
          [&](Configuration const &config) {
            return m_order_parameter != nullptr;
          }));

  config_query_dict.insert(GenericDatumFormatter<bool, Configuration>(
      "is_new", "is_new",
      [&](Configuration const &config) {
        return std::get<0>(m_data[config.name()]);
      },
      [&](Configuration const &config) {
        return m_data.find(config.name()) != m_data.end();
      }));

  config_query_dict.insert(GenericDatumFormatter<bool, Configuration>(
      "is_new_primitive", "is_new_primitive",
      [&](Configuration const &config) {
        return std::get<1>(m_data[config.name()]);
      },
      [&](Configuration const &config) {
        return m_data.find(config.name()) != m_data.end();
      }));

  config_query_dict.insert(GenericDatumFormatter<bool, Configuration>(
      "selected", "selected", [](Configuration const &config) { return true; },
      [](Configuration const &config) { return true; }));

  // set hall of fame metric
  m_halloffame.unique().reset(new HallOfFameType(
      MonteCarloEnumMetric(config_query_dict.parse(metric_args())),
      std::less<Configuration>(), set.enumeration_N_halloffame(),
      set.enumeration_tol()));

  // set pre-check function
  m_enum_check.unique().reset(new MonteCarloEnumCheck(
      config_query_dict.parse(set.enumeration_check_args())));

  // std::pair<double, Configuration> query dict
  for (const auto &datum_formatter : config_query_dict) {
    m_dict.insert(
        make_datum_formatter_adapter<PairType, Configuration>(datum_formatter));
  }

  m_dict.insert(GenericDatumFormatter<double, PairType>(
      "score", "score",
      [&](PairType const &score_and_config) { return score_and_config.first; },
      [&](PairType const &score_and_config) { return true; }));

  reset();

  _log().custom("Configuration enumeration");
  _log() << "  check: " << check_args() << "\n";
  _log() << "  metric: " << metric_args() << "\n";
  _log() << "  check_existence: " << std::boolalpha << m_check_existence
         << "\n";
  _log() << "  insert_canonical: " << std::boolalpha << m_insert_canonical
         << "\n";
  _log() << "  insert_primitive_only: " << std::boolalpha
         << m_insert_primitive_only << "\n";
  _log() << "  save_primitive_only: " << std::boolalpha << m_save_primitive_only
         << "\n";
  _log() << "  sample_mode: " << m_sample_mode << "\n";
  _log() << std::endl;
}

}  // namespace Monte
}  // namespace CASM

#endif
