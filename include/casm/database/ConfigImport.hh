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

    template<>
    class StructureMap<Configuration> : public ConfigData<Configuration> {

    public:

      /// Construct with PrimClex and by moving a ConfigMapper
      StructureMap<Configuration>(
        const PrimClex &_primclex,
        std::unique_ptr<ConfigMapper> mapper,
        bool primitive_only,
        std::vector<std::string> dof);

      /// Construct with PrimClex and settings (see Import / Update desc)
      StructureMap<Configuration>(
        const PrimClex &_primclex,
        const jsonParser &kwargs,
        bool primitive_only,
        std::vector<std::string> dof);

      typedef std::back_insert_iterator<std::vector<ConfigIO::Result> > map_result_inserter;

      /// \brief Specialized mapping method for Configuration
      ///
      /// \param p Path to structure or properties.calc.json file. Not guaranteed to exist or be valid.
      /// \param hint Iterator to 'from' config for 'casm update', or 'end' if unknown as with 'casm import'.
      /// \param result Insert iterator of Result objects to output mapping results
      ///
      /// - Should output one or more mapping results from the structure located at specied path
      /// - >1 result handles case of non-primitive configurations
      /// - responsible for filling in Result data structure
      /// - If 'hint' is not nullptr, use hint as 'from' config, else 'from' == 'to'
      map_result_inserter map(
        fs::path p,
        DatabaseIterator<ConfigType> hint,
        map_result_inserter result) const;

      /// Returns JSON with settings used after combing constructor input and defaults
      const jsonParser &used() const;


    private:

      /// \brief Import Configuration with only occupation DoF
      bool _occupation_only() const;

      /// \brief Read BasicStructure<Site> to be imported
      BasicStructure<Site> _make_structure(const fs::path &p) const;

      std::unique_ptr<ConfigMapper> m_configmapper;
      jsonParser m_used;
      bool m_primitive_only;
      std::vector<std::string> m_dof;
    };

    /// Configuration-specialized Import
    template<>
    class Import<Configuration> : public ImportT<Configuration> {
    public:

      /// \brief Constructor
      Import(
        const PrimClex &primclex,
        const StructureMap<Configuration> &mapper,
        bool import_data,
        bool copy_additional_files,
        bool overwrite,
        fs::path report_dir,
        Log &file_log);

      static const std::string desc;
      static int run(const PrimClex &primclex, const jsonParser &kwargs, const Completer::ImportOption &import_opt);

      using ImportT<Configuration>::import;

    protected:

      /// Allow ConfigType to specialize the report formatting for 'import'
      DataFormatter<ConfigIO::Result> _import_formatter(
        const std::map<std::string, ConfigIO::ImportData> &data_results) const override;

    };

    /// Configuration-specialized Import
    template<>
    class Update<Configuration> : public UpdateT<Configuration> {
    public:

      /// \brief Constructor
      Update(
        const PrimClex &primclex,
        const StructureMap<Configuration> &mapper,
        fs::path report_dir);

      static const std::string desc;
      static int run(const PrimClex &primclex, const jsonParser &kwargs, const Completer::UpdateOption &import_opt);

      using UpdateT<Configuration>::update;

    private:

      // Allow ConfigType to specialize the report formatting for 'update'
      DataFormatter<ConfigIO::Result> _update_formatter(
        const std::map<std::string, ConfigIO::ImportData> &data_results) const override;

    };

    template<>
    class Remove<Configuration> : public RemoveT<Configuration> {
    public:
      Remove(const PrimClex &primclex, fs::path report_dir, Log &_file_log);
      static const std::string desc;
      static int run(const PrimClex &primclex, const Completer::RmOption &import_opt);

    };

  }
}

#endif
