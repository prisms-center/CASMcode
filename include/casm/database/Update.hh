#ifndef CASM_DB_Update
#define CASM_DB_Update

#include "casm/database/ConfigData.hh"

namespace CASM {
  namespace Completer {
    class UpdateOption;
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

    template<typename T> class Selection;

    /// Generic ConfigType-dependent part of Import
    template<typename _ConfigType>
    class UpdateT : protected ConfigData<_ConfigType> {

    public:

      typedef _ConfigType ConfigType;
      using ConfigData<_ConfigType>::has_existing_files;
      using ConfigData<_ConfigType>::no_change;
      using ConfigData<_ConfigType>::db_props;
      using ConfigData<_ConfigType>::db_config;
      using ConfigData<_ConfigType>::db_supercell;
      using ConfigData<_ConfigType>::rm_files;
      using ConfigData<_ConfigType>::primclex;
      using Logging::log;

      /// \brief Constructor
      UpdateT(
        const PrimClex &primclex,
        const StructureMap<ConfigType> &mapper,
        fs::path report_dir) :
        ConfigData<_ConfigType>(primclex, null_log()),
        m_structure_mapper(mapper),
        m_report_dir(report_dir) {}

      /// \brief Re-parse calculations 'from' all selected configurations
      void update(const DB::Selection<ConfigType> &selection, bool force);

    protected:

      // Allow ConfigType to specialize the report formatting for 'update'
      virtual DataFormatter<ConfigIO::Result> _update_formatter() const = 0;

      void _update_report(std::vector<ConfigIO::Result> &results, const DB::Selection<ConfigType> &selection) const;

    private:

      const StructureMap<ConfigType> &m_structure_mapper;

      fs::path m_report_dir;
    };


    // To be specialized for ConfigType (no default implemenation exists)
    template<typename ConfigType>
    class Update : public UpdateT<ConfigType> {
    public:
      static const std::string desc;
      int run(const PrimClex &, const jsonParser &input, const Completer::UpdateOption &opt);

    private:
      // Allow ConfigType to specialize the report formatting for 'update'
      DataFormatter<ConfigIO::Result> _update_formatter() const override;
    };

  }
}

#endif
