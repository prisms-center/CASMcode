#include "casm/monte2/io/json/Results_json_io.hh"

#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/monte2/checks/CompletionCheck.hh"

namespace CASM {
namespace Monte2 {

// namespace {
//
// jsonParser &to_json(SamplerComponent const &key,
//                     IndividualEquilibrationCheckResult &equilibration_value,
//                     IndividualConvergenceCheckResult &convergence_value,
//                     jsonParser &json) {
//   auto &_json = json[key.sampler_name][std::to_string(key.component_index];
//
//   _json.put_obj();
//   to_json(equilibration_value.is_equilibrated, _json["is_equilibrated"]);
//   to_json(equilibration_value.N_samples_for_equilibration,
//   _json["N_samples_for_equilibration"]); if
//   (convergence_value.N_samples_for_statistics > 0) {
//     to_json(convergence_value.is_converged, _json["is_converged"]);
//     to_json(convergence_value.mean, _json["mean"]);
//     to_json(convergence_value.squared_norm, _json["squared_norm"]);
//     to_json(convergence_value.calculated_precision,
//     _json["calculated_precision"]);
//   }
//   return json;
// }
//
// void from_json(SamplerComponent &key,
//                IndividualEquilibrationCheckResult &equilibration_value,
//                IndividualConvergenceCheckResult &convergence_value,
//                jsonParser const &json) {
//   auto &_json = json[key.sampler_name][std::to_string(key.component_index];
//
//   from_json(equilibration_value.is_equilibrated, _json["is_equilibrated"]);
//   from_json(equilibration_value.N_samples_for_all_to_equilibrate,
//       _json["N_samples_for_all_to_equilibrate"]);
//   if (convergence_value.N_samples_for_statistics > 0) {
//     from_json(convergence_value.all_converged, _json["is_converged"]);
//     from_json(convergence_value.mean, _json["mean"]);
//     from_json(convergence_value.squared_norm, _json["squared_norm");
//     from_json(convergence_value.calculated_precision,
//     _json["calculated_precision"]);
//   }
// }
//
// jsonParser &to_json(CompletionCheckResults const &completion_check,
//                     jsonParser &json) {
//   json.put_obj();
//   to_json(completion_check.is_complete, json["is_complete"]);
//   to_json(completion_check.all_equilibrated, json["all_equilibrated"]);
//   to_json(completion_check.N_samples_for_all_to_equilibrate,
//           json["N_samples_for_all_to_equilibrate"]);
//   to_json(completion_check.all_converged, json["all_converged"]);
//   to_json(completion_check.N_samples_for_statistics,
//           json["N_samples_for_statistics"]);
//
//   auto &_json = json["statistics"];
//   auto it = e.individual_results.begin();
//   auto end = e.individual_results.end();
//   for (; it != end; ++it) {
//     SamplerComponent const &key = it->first;
//     IndividualEquilibrationCheckResult const &equilibration_value =
//     it->second; IndividualConvergenceCheckResult const &convergence_value =
//     completion_check.individual_results[key]; to_json(key,
//     equilibration_value, convergence_value, json["statistics"]);
//   }
//   return json;
// }
//
// void from_json(CompletionCheckResults &completion_check,
//                jsonParser const &json) {
//   from_json(completion_check.is_complete, json["is_complete"]);
//
//   EquilibrationCheckResults const &e =
//       completion_check.equilibration_check_results;
//   ConvergenceCheckResults &c =
//       completion_check.convergence_check_results;
//
//   from_json(e.all_equilibrated, json["all_equilibrated"]);
//   from_json(e.N_samples_for_all_to_equilibrate,
//             json["N_samples_for_all_to_equilibrate"]);
//
//   from_json(c.all_converged, json["all_converged"]);
//   from_json(c.N_samples_for_statistics, json["N_samples_for_statistics"]);
//
//   auto it = json["statistics"].begin();
//   auto end = json["statistics"].end();
//   for (; it != end; ++it) {
//     SamplerComponent key;
//     IndividualEquilibrationCheckResult equilibration_value;
//     IndividualConvergenceCheckResult convergence_value;
//     from_json(key, equilibration_value, convergence_value, *it);
//     e.individual_results.emplace(key, equilibration_value);
//     c.individual_results.emplace(key, convergence_value);
//   }
// }
//
// jsonParser &to_json(SampledData const &data, jsonParser &json) {
//   json.put_obj();
//   for (auto const &ref : data.samplers) {
//     json[ref.first] = ref.second.values()
//   }
//   if (data.count.size()) {
//     json["count"].put_array(data.count.begin(), data.count.end());
//   }
//   if (data.time.size()) {
//     json["time"].put_array(data.time.begin(), data.time.end());
//   }
//   if (data.trajectory.size()) {
//     results["trajectory"] = results.trajectory;
//   }
// }
//
// void parse(
//     InputParser<SampledData> &parser,
//     std::shared_ptr<Structure const> const &shared_prim,
//     std::map<std::string, StateSamplingFunction> const &sampling_functions) {
//   parser.value = notstd::make_unique<SampledData>();
//   SampledData &data = *parser.value;
//   jsonParser &samples = json["samples"];
//   for (auto it = samples.begin; it != samples.end; ++it) {
//     auto function_it = sampling_functions.find(it.name());
//     if (function_it == sampling_functions.end()) {
//       std::stringstream msg;
//       msg << "Sampling function not found: '" << it.name() << "'";
//       throw std::runtime_error(msg.str());
//     }
//     auto result =
//         data.samplers.emplace(it.name(), Sampler(function_it->second));
//     result.first->set_values(it->get<Eigen::MatrixXd>());
//   }
//   parser.optional(data.count, "count");
//   parser.optional(data.time, "time");
//   parser.optional(data.trajectory, "trajectory", shared_prim);
// }
//
// jsonParser &to_json(Results const &results, jsonParser &json) {
//   json.put_obj();
//   to_json(results.initial_configuration, json["initial_configuration"]);
//   to_json(results.conditions, json["conditions"]);
//   to_json(results.final_configuration, json["final_configuration"]);
//   to_json(results.statistics, json["statistics"]);
//   to_json(results.additional_properties, json["additional_properties"]);
//   to_json(results.completion_check_results, json["completion_check"]);
//   if (results.sampled_data) {
//     to_json(*results.sampled_data, json["samples"]);
//   }
//   return json;
// }
//
// void parse(
//     InputParser<Results> &parser,
//     std::shared_ptr<Structure const> const &shared_prim,
//     std::map<std::string, StateSamplingFunction> &state_sampling_functions) {
//   parser.value = notstd::make_unique<Results>();
//   Results &results = *parser.value;
//
//   parser.optional(results.initial_configuration, "initial_configuration",
//                   shared_prim);
//   parser.optional(results.conditions, "conditions", shared_prim);
//   parser.optional(results.final_configuration, "final_configuration",
//                   shared_prim);
//   parser.optional(results.statistics, "statistics");
//   parser.optional(results.additional_properties, "additional_properties");
//   parser.optional(results.completion_check_results, "completion_check");
//
//   {
//     auto subparser =
//         parser.subparse_if<SampledData>("samples", state_sampling_functions);
//     if (subparser->value) {
//       results.sampled_data = std::move(subparser->value);
//     }
//   }
// }
//
// std::string _filename(CountType run_index) {
//   return std::string("results.") + std::to_string(run_index) + ".json";
// }
// }  // namespace
//
// JSONResultsIO::JSONResultsIO(
//     fs::path output_dir, bool write_observations,
//     std::shared_ptr<Structure const> const &shared_prim)
//     : m_output_dir(output_dir),
//       m_write_observations(write_observations),
//       m_shared_prim(shared_prim) {}
//
// void JSONResultsIO::write_results(Results const &results,
//                                   Index run_index) const {
//   jsonParser json;
//   to_json(results, json);
//   fs::path file = m_output_dir / _filename(run_index);
//   results.write(file);
// }
//
// JSONResultsIO::read_results(Results &results, Index run_index) const {
//   fs::path file = m_output_dir / _filename(run_index);
//   jsonParser json{file};
//
//   InputParser<Results> parser{json, m_shared_prim};
//   report_and_throw_if_invalid();
//   results = *parser.value;
//
//   // results.initial_configuration =
//   //     json["initial_configuration"].get<Configuration>();
//   // results.final_configuration =
//   //     json["final_configuration"].get<Configuration>();
//   // results.conditions = json["conditions"].get<VectorValueMap>();
//   //
//   // if (json.contains("observations")) {
//   //   Index n_samples;
//   //   auto it = json["observations"].begin();
//   //   n_samples = it->size();
//   //   for (; it != json["observations"].end(); ++it) {
//   //     auto function_it = sampling_functions.find(it.name());
//   //     if (function_it == sampling_functions.end()) {
//   //       std::stringstream msg;
//   //       msg << "Sampling function not found: '" << it.name() << "'";
//   //       throw std::runtime_error(msg.str());
//   //     }
//   //     auto result = results.sampled_data.samplers->emplace(
//   //         it.name(), Sampler(function_it->second));
//   //     result.first->set_values(it->get<Eigen::MatrixXd>());
//   //   }
//   //
//   //   _from_json(results.data.sample_times, json);
//   // }
//   // return Results(initial_configuration, final_configuration, conditions,
//   // data);
// }
//
// void JSONResultsIO::write_completed_runs_summary(
//     std::vector<Results> const &results_summary) const {}
//
// void JSONResultsIO::read_completed_runs_summary(
//     std::vector<Results> &results_summary) const {}

}  // namespace Monte2
}  // namespace CASM
