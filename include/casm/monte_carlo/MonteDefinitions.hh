#ifndef CASM_MonteDefinitions_HH
#define CASM_MonteDefinitions_HH

#include "casm/CASM_global_definitions.hh"
#include "casm/casm_io/EnumIO.hh"

namespace CASM {
  namespace Monte {

    typedef Index size_type;

    /// \brief Monte Carlo ensemble type
    enum class ENSEMBLE {
      GrandCanonical,
      Canonical
    };

    ENUM_IO_DECL(CASM::Monte::ENSEMBLE)


    /// \brief Monte Carlo method type
    enum class METHOD {
      Metropolis, LTE1
    };

    ENUM_IO_DECL(CASM::Monte::METHOD)


    ///How often to sample runs
    enum class SAMPLE_MODE {
      STEP, PASS
    };

    ENUM_IO_DECL(CASM::Monte::SAMPLE_MODE)


    ///How to change conditions
    enum class DRIVE_MODE {
      INCREMENTAL, CUSTOM
    };

    ENUM_IO_DECL(CASM::Monte::DRIVE_MODE)

    ///How often to sample runs
    enum class ENUM_SAMPLE_MODE {
      ON_SAMPLE, ON_ACCEPT
    };

    ENUM_IO_DECL(CASM::Monte::ENUM_SAMPLE_MODE)

  }
}

namespace CASM {

  ENUM_TRAITS(Monte::ENSEMBLE)

  ENUM_TRAITS(Monte::METHOD)

  ENUM_TRAITS(Monte::SAMPLE_MODE)

  ENUM_TRAITS(Monte::DRIVE_MODE)

  ENUM_TRAITS(Monte::ENUM_SAMPLE_MODE)


}
#endif


