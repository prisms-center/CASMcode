#ifndef CASM_MonteDefinitions_HH
#define CASM_MonteDefinitions_HH

#include <string>
#include <stdexcept>
#include <iostream>
#include "casm/external/boost.hh"

namespace CASM {

  namespace Monte {

    /// \brief Monte Carlo ensemble type
    enum class ENSEMBLE {
      GrandCanonical
    };

    inline std::ostream &operator<<(std::ostream &sout, ENSEMBLE ensemble) {
      if(ensemble == ENSEMBLE::GrandCanonical) {
        sout << "GrandCanonical";
      }
      return sout;
    }

    inline Monte::ENSEMBLE monte_ensemble(std::string ensemble) {
      if(ensemble == "GrandCanonical" || "grand_canonical") {
        return Monte::ENSEMBLE::GrandCanonical;
      }
      else {
        throw std::runtime_error(
          std::string("Error in 'monte_ensemble(std::string ensemble)'\n") +
          "  ensemble = " + ensemble + " is not allowed.\n" +
          "  Options are: \n" +
          "    'GrandCanonical' or 'grand_canonical'");
      }
    }


    /// \brief Monte Carlo method type
    enum class METHOD {
      Metropolis, LTE1
    };

    inline std::ostream &operator<<(std::ostream &sout, METHOD method) {
      if(method == METHOD::Metropolis) {
        sout << "Metropolis";
      }
      else if(method == METHOD::LTE1) {
        sout << "LTE1";
      }
      return sout;
    }

    inline Monte::METHOD monte_method(std::string method) {
      if(method == "Metropolis" || method == "metropolis") {
        return Monte::METHOD::Metropolis;
      }
      else if(method == "LTE1" || method == "lte1") {
        return Monte::METHOD::LTE1;
      }
      else {
        throw std::runtime_error(
          std::string("Error in 'monte_method(std::string type)'\n") +
          "  method = " + method + " is not allowed.\n" +
          "  Options are: \n" +
          "    'Metropolis' or 'metropolis'\n" +
          "    'LTE1' or 'lte1'");
      }
    }


    ///How often to sample runs
    enum class SAMPLE_MODE {
      STEP, PASS
    };

    inline std::ostream &operator<<(std::ostream &sout, SAMPLE_MODE mode) {
      if(mode == SAMPLE_MODE::STEP) {
        sout << "step";
      }
      else if(mode == SAMPLE_MODE::PASS) {
        sout << "pass";
      }
      return sout;
    }


    ///How to change conditions after a point has equilibrated
    enum class DRIVE_MODE {
      SINGLE, INCREMENTAL
    };

    inline std::ostream &operator<<(std::ostream &sout, DRIVE_MODE mode) {
      if(mode == DRIVE_MODE::SINGLE) {
        sout << "single";
      }
      else if(mode == DRIVE_MODE::INCREMENTAL) {
        sout << "incremental";
      }
      return sout;
    }


    ///What kind of sampling to get for correlations (none, partial, only non-zero, all)
    enum class CORR_SAMPLE_MODE {
      NONE, CUSTOM, NONZERO, ALL
    };


  }

}
#endif


