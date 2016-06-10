#ifndef CASM_MonteDefinitions_HH
#define CASM_MonteDefinitions_HH

#include "casm/casm_io/EnumIO.hh"

namespace CASM {

  namespace Monte {

    /// \brief Monte Carlo ensemble type
    enum class ENSEMBLE {
      GrandCanonical
    };

    ENUM_IO(CASM::Monte::ENSEMBLE)


    /// \brief Monte Carlo method type
    enum class METHOD {
      Metropolis, LTE1
    };

    ENUM_IO(CASM::Monte::METHOD)


    ///How often to sample runs
    enum class SAMPLE_MODE {
      STEP, PASS
    };

    ENUM_IO(CASM::Monte::SAMPLE_MODE)


    ///How to change conditions
    enum class DRIVE_MODE {
      INCREMENTAL, CUSTOM
    };

    ENUM_IO(CASM::Monte::DRIVE_MODE)


  }

  ENUM_TRAITS(Monte::ENSEMBLE)

  ENUM_TRAITS(Monte::METHOD)

  ENUM_TRAITS(Monte::SAMPLE_MODE)

  ENUM_TRAITS(Monte::DRIVE_MODE)


}
#endif


