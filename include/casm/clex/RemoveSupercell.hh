#ifndef CASM_RemoveSupercell
#define CASM_RemoveSupercell

#include <string>
#include "casm/CASM_global_definitions.hh"
#include "casm/database/Remove.hh"

namespace CASM {

  class PrimClex;
  class Supercell;
  class Log;

  namespace Completer {
    class RmOption;
  }

  namespace DB {

    template<>
    class Remove<Supercell> {

    public:
      Remove(const PrimClex &primclex, fs::path report_dir, Log &_file_log);

      /// \brief Erase all enumerated configs that have no data; if no data, erase supercell
      void erase(const DB::Selection<Supercell> &selection, bool dry_run);

      /// \brief Erase data and files (permanently), but not configs / supercells
      void erase_data(const DB::Selection<Supercell> &selection, bool dry_run);

      /// \brief Removes supercell, including all configurations, data, and files (permanently)
      void erase_all(const DB::Selection<Supercell> &selection, bool dry_run);

      static const std::string desc;
      static int run(const PrimClex &primclex, const Completer::RmOption &opt);

      const PrimClex &primclex() const;

      fs::path report_dir() const;

      Log &file_log() const;

    private:
      const PrimClex &m_primclex;
      fs::path m_report_dir;
      Log &m_file_log;
    };

  }
}

#endif
