#include "casm/monte_carlo/MonteDefinitions.hh"

namespace CASM {

  namespace Monte {

    const std::string traits<ENSEMBLE>::name = "ensemble";

    const std::multimap<ENSEMBLE, std::vector<std::string> > traits<ENSEMBLE>::strval = {
      {ENSEMBLE::GrandCanonical, {"GrandCanonical", "grand_canonical"} }
    };


    const std::string traits<METHOD>::name = "method";

    const std::multimap<METHOD, std::vector<std::string> > traits<METHOD>::strval = {
      {METHOD::Metropolis, {"Metropolis", "metropolis"} },
      {METHOD::LTE1, {"LTE1", "lte1"} }
    };


    const std::string traits<SAMPLE_MODE>::name = "sample_by";

    const std::multimap<SAMPLE_MODE, std::vector<std::string> > traits<SAMPLE_MODE>::strval = {
      {SAMPLE_MODE::STEP, {"Step", "step"} },
      {SAMPLE_MODE::PASS, {"Pass", "pass"} }
    };


    const std::string traits<DRIVE_MODE>::name = "mode";

    const std::multimap<DRIVE_MODE, std::vector<std::string> > traits<DRIVE_MODE>::strval = {
      {DRIVE_MODE::INCREMENTAL, {"Incremental", "incremental"} },
      {DRIVE_MODE::CUSTOM, {"Custom", "custom"} }
    };
  }
}

