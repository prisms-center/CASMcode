#ifndef CASM_DB_ConfigImport
#define CASM_DB_ConfigImport

#include <vector>
#include <map>
#include <string>
#include <memory>
#include "casm/global/definitions.hh"
#include "casm/database/Update.hh"
#include "casm/database/Import.hh"

namespace CASM {
  namespace xtal {
    template<typename T> class BasicStructure;
    class Site;
    class SimpleStructure;
  }
  using xtal::BasicStructure;
  using xtal::Site;
  using xtal::SimpleStructure;


  class ConfigMapper;
  class Configuration;
  class PrimClex;
  class jsonParser;
}

namespace CASM {
  namespace DB {

    template<>
    class StructureMap<Configuration> {

    public:

      /// Construct with PrimClex and by moving a ConfigMapper
      StructureMap<Configuration>(MappingSettings const &_set,
                                  std::unique_ptr<ConfigMapper> mapper);


      /// Construct with PrimClex and settings (see Import / Update desc)
      StructureMap<Configuration>(MappingSettings const &_set,
                                  const PrimClex &_primclex);



      typedef std::back_insert_iterator<std::vector<ConfigIO::Result> > map_result_inserter;

      /// \brief Specialized mapping method for Configuration
      ///
      /// \param p Path to structure or properties.calc.json file. Not guaranteed to exist or be valid.
      /// \param hint std::unique_ptr<Configuration> to 'from' config for 'casm update', or 'end' if unknown as with 'casm import'.
      /// \param result Insert iterator of Result objects to output mapping results
      ///
      /// - Should output one or more mapping results from the structure located at specied path
      /// - >1 result handles case of non-primitive configurations
      /// - responsible for filling in Result data structure
      /// - If 'hint' is not nullptr, use hint as 'from' config, else 'from' == 'to'
      map_result_inserter map(fs::path p,
                              std::unique_ptr<Configuration> const &hint_config,
                              map_result_inserter result) const;

      /// Returns settings used for import
      const MappingSettings &settings() const {
        return m_set;
      }


    private:

      /// \brief Read SimpleStructure to be imported
      SimpleStructure _make_structure(const fs::path &p) const;

      MappingSettings m_set;
      std::unique_ptr<ConfigMapper> m_configmapper;
    };

    /// Configuration-specialized Import
    template<>
    class Import<Configuration> : public ImportT<Configuration> {
    public:
      using ImportT<ConfigType>::settings;
      using ImportT<ConfigType>::import;

      /// \brief Constructor
      Import(
        const PrimClex &primclex,
        const StructureMap<ConfigType> &mapper,
        ImportSettings const &_set,
        fs::path const &report_dir,
        Log &file_log);

      static const std::string desc;
      static int run(const PrimClex &primclex, const jsonParser &kwargs, const Completer::ImportOption &import_opt);


    protected:

      /// Allow ConfigType to specialize the report formatting for 'import'
      DataFormatter<ConfigIO::Result> _import_formatter(
        const std::map<std::string, ConfigIO::ImportData> &data_results) const override;

    };

    /// Configuration-specialized Import
    template<>
    class Update<Configuration> : public UpdateT<Configuration> {
    public:

      using UpdateT<ConfigType>::update;

      /// \brief Constructor
      Update(
        const PrimClex &primclex,
        const StructureMap<ConfigType> &mapper,
        fs::path const &report_dir);

      static const std::string desc;
      static int run(const PrimClex &primclex, const jsonParser &kwargs, const Completer::UpdateOption &import_opt);


    private:

      // Allow ConfigType to specialize the report formatting for 'update'
      DataFormatter<ConfigIO::Result> _update_formatter() const override;

    };

  }
}

#endif
