#include "casm/monte_carlo/MonteDefinitions.hh"

namespace CASM {

  namespace CASM_TMP {

    const std::string traits<Monte::ENSEMBLE>::name = "ensemble";

    const std::multimap<Monte::ENSEMBLE, std::vector<std::string> > traits<Monte::ENSEMBLE>::strval = {
      {Monte::ENSEMBLE::GrandCanonical, {"GrandCanonical", "grand_canonical"} }
    };


    const std::string traits<Monte::METHOD>::name = "method";

    const std::multimap<Monte::METHOD, std::vector<std::string> > traits<Monte::METHOD>::strval = {
      {Monte::METHOD::Metropolis, {"Metropolis", "metropolis"} },
      {Monte::METHOD::LTE1, {"LTE1", "lte1"} }
    };


    const std::string traits<Monte::SAMPLE_MODE>::name = "sample_by";

    const std::multimap<Monte::SAMPLE_MODE, std::vector<std::string> > traits<Monte::SAMPLE_MODE>::strval = {
      {Monte::SAMPLE_MODE::STEP, {"Step", "step"} },
      {Monte::SAMPLE_MODE::PASS, {"Pass", "pass"} }
    };


    const std::multimap<Monte::DRIVE_MODE, std::vector<std::string> > traits<Monte::DRIVE_MODE>::strval = {
      {Monte::DRIVE_MODE::INCREMENTAL, {"Incremental", "incremental"} },
      {Monte::DRIVE_MODE::CUSTOM, {"Custom", "custom"} }
    };

    const std::string traits<Monte::ENUM_SAMPLE_MODE>::name = "sample_mode";

    const std::multimap<Monte::ENUM_SAMPLE_MODE, std::vector<std::string> > traits<Monte::ENUM_SAMPLE_MODE>::strval = {
      {Monte::ENUM_SAMPLE_MODE::ON_ACCEPT, {"on_accept"} },
      {Monte::ENUM_SAMPLE_MODE::ON_SAMPLE, {"on_sample"} }
    };
  }
}

