#ifndef CASM_RemovePrimPeriodicDiffTransOrbit
#define CASM_RemovePrimPeriodicDiffTransOrbit

#include <string>
#include "casm/global/definitions.hh"
#include "casm/kinetics/PrimPeriodicDiffTransOrbitTraits.hh"
#include "casm/database/Remove.hh"

namespace CASM {

  class PrimClex;
  class Log;

  namespace Completer {
    class RmOption;
  }

  namespace DB {

    template<>
    class Remove<PrimPeriodicDiffTransOrbit> {
      /// THIS IS LIKELY BUGGY/ NOT IMPLEMENTED
    public:
      Remove(const PrimClex &primclex, std::string report_dir, Log &_file_log);

      /// \brief Erase all enumerated difftransconfigs that have no data; if no data, erase PrimPeriodicDiffTransOrbit
      void erase(const DB::Selection<PrimPeriodicDiffTransOrbit> &selection, bool dry_run);

      /// \brief Erase data and files (permanently), but not difftransconfigs / PrimPeriodicDiffTransOrbit
      void erase_data(const DB::Selection<PrimPeriodicDiffTransOrbit> &selection, bool dry_run);

      /// \brief Removes PrimPeriodicDiffTransorbit, including all difftransconfigurations, data, and files (permanently)
      void erase_all(const DB::Selection<PrimPeriodicDiffTransOrbit> &selection, bool dry_run);

      static std::string desc();
      static int run(const PrimClex &primclex, const Completer::RmOption &opt);

      const PrimClex &primclex() const;

      std::string report_dir() const;

      Log &file_log() const;

    private:
      const PrimClex &m_primclex;
      std::string m_report_dir;
      Log &m_file_log;
    };

  }
}

#endif
