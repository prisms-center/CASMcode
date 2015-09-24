#ifndef CASM_MonteDefinitions_HH
#define CASM_MonteDefinitions_HH

#include <string>

namespace CASM {

  namespace Monte {
    
    /// \brief Type of Monte Carlo calculation
    enum class TYPE {
      GrandCanonical
    };
    
    ///How often to sample runs
    enum class SAMPLE_MODE {
      STEP, PASS
    };
    
    ///How to change conditions after a point has equilibrated
    enum class DRIVE_MODE {
      SINGLE, INCREMENTAL
    };
    
    ///What kind of sampling to get for correlations (none, partial, only non-zero, all)
    enum class CORR_SAMPLE_MODE {
      NONE, CUSTOM, NONZERO, ALL
    };
    
    inline Monte::TYPE monte_type(std::string type) {
      if(type == "GrandCanonical" || "grand_canonical") {
        return Monte::TYPE::GrandCanonical;
      }
      else {
        throw std::runtime_error(
          std::string("Error in 'monte_type(std::string type)'\n") +
                      "  type = " + type + " is not allowed.\n" +
                      "  Options are: \n" +
                      "    'GrandCanonical' or 'grand_canonical'");
      }
    }
  }
}
#endif


