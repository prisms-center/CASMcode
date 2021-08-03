#ifndef CASM_monte2_ResultsIO
#define CASM_monte2_ResultsIO

#include <boost/filesystem.hpp>
#include <vector>

#include "casm/global/definitions.hh"

namespace CASM {
namespace Monte2 {

// struct Results;
//
// /// \brief Abstract base class interface for results input / output
// class ResultsIO {
//  public:
//   /// \brief Write results of single conditions, including SampledData
//   virtual void write_results(Results const &results, Index run_index) const =
//   0;
//
//   /// \brief Read results of single conditions, including SampledData
//   virtual void read_results(Results &results, Index run_index) const = 0;
//
//   /// \brief Write combined results from all conditions, excluding
//   SampledData virtual void write_results_summary(
//       std::vector<Results> const &results) const = 0;
//
//   /// \brief Read combined results from all conditions, excluding SampledData
//   virtual void read_results_summary(std::vector<Results> &results) const = 0;
// };

}  // namespace Monte2
}  // namespace CASM

#endif
