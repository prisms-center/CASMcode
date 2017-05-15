#ifndef CASM_DB_Remove
#define CASM_DB_Remove

#include "casm/database/ConfigData.hh"

namespace CASM {
  namespace Completer {
    class RmOption;
  }
}


// To be specialized for calculable 'ConfigType' classes:
//   StructureMap<ConfigType>::map(fs::path, DatabaseIterator<ConfigType> hint, inserter result)
//   Import<ConfigType>::desc
//   Import<ConfigType>::run
//   Import<ConfigType>::_import_formatter
//   Update<ConfigType>::desc
//   Update<ConfigType>::run
//   Update<ConfigType>::_update_formatter
//   Remove<ConfigType>::desc
//   Remove<ConfigType>::run

namespace CASM {
  namespace DB {

    /// Generic ConfigType-dependent part of Remove
    template<typename _ConfigType>
    class RemoveT : protected ConfigData<_ConfigType> {

    public:

      typedef _ConfigType ConfigType;
      using ConfigData<_ConfigType>::db_config;
      using ConfigData<_ConfigType>::db_props;
      using ConfigData<_ConfigType>::has_existing_data_or_files;
      using ConfigData<_ConfigType>::rm_files;
      using Logging::log;

      RemoveT(const PrimClex &primclex, fs::path report_dir, Log &_file_log) :
        ConfigData<_ConfigType>(primclex, _file_log),
        m_report_dir(report_dir) {}

      /// \brief Erase Configurations that have no data
      void erase(const DB::Selection<ConfigType> &selection, bool dry_run);

      /// \brief Erase data and files (permanently), but not Configuration
      void erase_data(const DB::Selection<ConfigType> &selection, bool dry_run);

      /// \brief Removes Configurations and data and files (permanently)
      ///
      /// - Data are always associated with one 'from' configuration, so the
      ///   selection here indicates 'from' configurations
      /// - The 'to' configurations are updated with the new best mapping properties
      void erase_all(const DB::Selection<ConfigType> &selection, bool dry_run);


    private:

      void _erase_report(const std::vector<std::string> &fail);

      fs::path m_report_dir;
    };

    // To be specialized for ConfigType (default implemenation exists for ConfigTypes)
    template<typename ConfigType>
    class Remove : public RemoveT<ConfigType> {
    public:
      Remove(const PrimClex &primclex, fs::path report_dir, Log &_file_log);
      static const std::string desc;
      static int run(const PrimClex &, const Completer::RmOption &opt);
    };
  }
}

#endif
