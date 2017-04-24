#ifndef CASM_DB_DiffTransConfigImport
#define CASM_DB_DiffTransConfigImport

#include "casm/database/Import.hh"

namespace CASM {
  namespace DB {

    /*
    /// DiffTransConfiguration-specialized Import
    template<>
    class Import<DiffTransConfiguration> : public ImportT<DiffTransConfiguration> {
    public:

      /// \brief Constructor
      Import(
        const PrimClex& primclex,
        const DiffTransConfigMapper& configmapper,
        bool import_data,
        bool copy_additional_files,
        bool overwrite,
        fs::path report_dir = primclex.root_dir() / "import_report");

      static const std::string import_help;
      static int import(PrimClex &primclex, const jsonParser &kwargs, const Completer::ImportOption &import_opt);

      static const std::string update_help;
      static int update(PrimClex &primclex, const jsonParser &kwargs, const Completer::UpdateOption &import_opt);

      static const std::string remove_help;
      static int remove(PrimClex &primclex, const jsonParser &kwargs, const Completer::RemoveOption &import_opt);

    private:

      /// \brief Specialized import method for ConfigType
      import_inserter _import(
        fs::path p,
        DataBaseIterator<DiffTransConfiguration> hint,
        import_inserter result) override;

      /// Allow ConfigType to specialize the report formatting for 'import'
      DataFormatter<Result> _import_formatter(
        const std::map<std::string, ImportData>& data_results) const override;

      // Allow ConfigType to specialize the report formatting for 'update'
      DataFormatter<Result> _update_formatter(
        const std::map<std::string, ImportData>& data_results) const override;

    private:
      const DiffTransConfigMapper& m_configmapper;

    };
    */

  }
}

#endif
