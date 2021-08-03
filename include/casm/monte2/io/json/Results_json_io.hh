#ifndef CASM_monte2_Results_json_io
#define CASM_monte2_Results_json_io

#include <map>
#include <vector>

#include "casm/global/definitions.hh"
#include "casm/monte2/io/ResultsIO.hh"

namespace CASM {

class Structure;
template <typename T>
class InputParser;
class jsonParser;

namespace Monte2 {

// struct CompletionCheckResults;
// struct IndividualConvergenceCheckResult;
// struct IndividualEquilibrationCheckResult;
// struct Results;
// struct SampledData;
// struct SamplerComponent;
// struct StateSamplingFunction;
//
// jsonParser &to_json(SamplerComponent const &key,
//                     IndividualEquilibrationCheckResult &equilibration_value,
//                     IndividualConvergenceCheckResult &convergence_value,
//                     jsonParser &json);
// void from_json(SamplerComponent &key,
//                IndividualEquilibrationCheckResult &equilibration_value,
//                IndividualConvergenceCheckResult &convergence_value,
//                jsonParser const &json);
//
// jsonParser &to_json(CompletionCheckResults const &completion_check,
//                     jsonParser &json);
// void from_json(CompletionCheckResults &completion_check,
//                jsonParser const &json);
//
// jsonParser &to_json(SampledData const &data, jsonParser &json);
// void parse(
//     InputParser<SampledData> &parser,
//     std::shared_ptr<Structure const> const &shared_prim,
//     std::map<std::string, StateSamplingFunction> const &sampling_functions);
//
// jsonParser &to_json(Results const &results, jsonParser &json);
// void parse(
//     InputParser<Results> &parser,
//     std::shared_ptr<Structure const> const &shared_prim,
//     std::map<std::string, StateSamplingFunction> &state_sampling_functions);
//
// class JSONResultsIO : public ResultsIO {
//  public:
//   JSONResultsIO(fs::path output_dir, bool write_observations,
//                 bool write_trajectory,
//                 std::shared_ptr<Structure const> const &shared_prim);
//
//   /// \brief Write results of single conditions, including SampledData
//   void write_results(Results const &results, Index run_index) const override;
//
//   /// \brief Read results of single conditions, including SampledData
//   void read_results(Results &results, Index run_index) const override;
//
//   /// \brief Write combined results from all conditions, excluding
//   SampledData void write_results_summary(
//       std::vector<Results> const &results) const override;
//
//   /// \brief Read combined results from all conditions, excluding SampledData
//   void read_results_summary(std::vector<Results> &results) const override;
//
//  private:
//   fs::path m_output_dir;
//
//   bool m_write_observations;
//
//   std::shared_ptr<Structure const> m_shared_prim;
// };

}  // namespace Monte2
}  // namespace CASM

#endif
