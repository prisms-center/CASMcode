#ifndef CASM_DB_Import
#define CASM_DB_Import

#include "casm/database/ConfigData.hh"

namespace CASM {
  namespace Completer {
    class ImportOption;
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

    /// Construct pos_paths from input args --pos && --batch
    template<typename OutputIterator>
    std::pair<OutputIterator, int> construct_pos_paths(
      const PrimClex &primclex,
      const Completer::ImportOption &import_opt,
      OutputIterator result);


    /// Generic ConfigType-dependent part of Import
    template<typename _ConfigType>
    class ImportT : protected ConfigData<_ConfigType> {

    public:

      typedef _ConfigType ConfigType;

      using ConfigData<_ConfigType>::has_existing_data_or_files;
      using ConfigData<_ConfigType>::db_props;
      using ConfigData<_ConfigType>::db_config;
      using ConfigData<_ConfigType>::db_supercell;
      using ConfigData<_ConfigType>::cp_files;
      using ConfigData<_ConfigType>::rm_files;
      using ConfigData<_ConfigType>::primclex;
      using Logging::log;

      /// \brief Constructor
      ImportT(
        const PrimClex &primclex,
        const StructureMap<ConfigType> &mapper,
        bool import_data,
        bool copy_additional_files,
        bool overwrite,
        fs::path report_dir,
        Log &_file_log) :
        ConfigData<_ConfigType>(primclex, _file_log),
        m_structure_mapper(mapper),
        m_import_data(import_data),
        m_copy_additional_files(copy_additional_files),
        m_overwrite(overwrite),
        m_report_dir(report_dir) {}

      template<typename PathIterator>
      void import(PathIterator begin, PathIterator end);

    private:

      /// Allow ConfigType to specialize the report formatting for 'import'
      virtual DataFormatter<ConfigIO::Result> _import_formatter(
        const std::map<std::string, ConfigIO::ImportData> &data_results) const = 0;

      void _import_report(
        std::vector<ConfigIO::Result> &results,
        const std::map<std::string, ConfigIO::ImportData> &data_results);


      const StructureMap<ConfigType> &m_structure_mapper;

      // attempt to import calculation results into database, else just insert
      // configurations w/out data
      bool m_import_data;

      // attempt to copy extra files from the directory where the structure is
      // being imported from to the training_data directory
      bool m_copy_additional_files;

      // Allow overwriting of existing data by 'casm import'
      bool m_overwrite;

      fs::path m_report_dir;
    };


    // To be specialized for ConfigType (no default implemenation exists)
    template<typename ConfigType>
    class Import : public ImportT<ConfigType> {
    public:
      static const std::string desc;
      int run(const PrimClex &, const jsonParser &input, const Completer::ImportOption &opt);

    private:
      // Allow ConfigType to specialize the report formatting for 'import'
      DataFormatter<ConfigIO::Result> _import_formatter(
        const std::map<std::string, ConfigIO::ImportData> &data_results) const override;
    };

  }
}

#endif
