#ifndef CASM_DB_ConfigImport
#define CASM_DB_ConfigImport

#include <vector>
#include <map>
#include <string>
#include <memory>
#include "casm/CASM_global_definitions.hh"
#include "casm/database/Import.hh"

namespace CASM {
  class ConfigMapper;
  class Configuration;
  class PrimClex;
  class jsonParser;
  template<typename T> class BasicStructure;
  class Site;
}

namespace CASM {
  namespace DB {

    template<typename T> class DatabaseIterator;

    /// Configuration-specialized Import
    template<>
    class Import<Configuration> : public ImportT<Configuration> {
    public:

      /// \brief Constructor
      Import(
        const PrimClex &primclex,
        const ConfigMapper &configmapper,
        std::vector<std::string> dof,
        bool primitive_only,
        bool import_data,
        bool copy_additional_files,
        bool overwrite,
        fs::path report_dir);

      static const std::string import_help;
      static int import(const PrimClex &primclex, const jsonParser &kwargs, const Completer::ImportOption &import_opt);

      static const std::string update_help;
      static int update(const PrimClex &primclex, const jsonParser &kwargs, const Completer::UpdateOption &import_opt);

      static const std::string remove_help;
      static int remove(const PrimClex &primclex, const jsonParser &kwargs, const Completer::RemoveOption &import_opt);

      using ImportT<Configuration>::import;
      using ImportT<Configuration>::update;

    private:

      /// Construct ConfigMapper from input args
      static std::pair<ConfigMapper, jsonParser> _make_configmapper(
        const PrimClex &primclex,
        const jsonParser &kwargs);

      /// \brief Path to default calctype training_data directory for config
      fs::path _calc_dir(const std::string configname) const override;

      /// \brief Specialized import method for ConfigType
      import_inserter _import(
        fs::path p,
        DatabaseIterator<Configuration> hint,
        import_inserter result) const override;

      /// Allow ConfigType to specialize the report formatting for 'import'
      DataFormatter<Result> _import_formatter(
        const std::map<std::string, ImportData> &data_results) const override;

      // Allow ConfigType to specialize the report formatting for 'update'
      DataFormatter<Result> _update_formatter(
        const std::map<std::string, ImportData> &data_results) const override;

      /// \brief Read BasicStructure<Site> to be imported
      BasicStructure<Site> _make_structure(const fs::path &p) const;

      /// \brief Import Configuration with only occupation DoF
      bool _occupation_only() const;


      const ConfigMapper &m_configmapper;
      std::vector<std::string> m_dof;
      bool m_primitive_only;
    };

  }
}

#endif
