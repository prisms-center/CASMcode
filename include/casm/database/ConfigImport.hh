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
    class Site;
    class SimpleStructure;
  }
  using xtal::SimpleStructure;

  namespace ConfigMapping {
    struct Settings;
  }

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


      /// Construct with PrimClex and settings (see Import / Update desc)
      StructureMap<Configuration>(ConfigMapping::Settings const &_set,
                                  const PrimClex &_primclex);



      typedef std::back_insert_iterator<std::vector<ConfigIO::Result> > map_result_inserter;

      /// \brief Specialized mapping method for Configuration
      ///
      /// \param p Path to structure or properties.calc.json file. Not guaranteed to exist or be valid.
      /// \param req_properties, list of names of properties that are required for mapped data to be considered 'complete'
      /// \param hint std::unique_ptr<Configuration> to 'from' config for 'casm update', or 'end' if unknown as with 'casm import'.
      /// \param result Insert iterator of Result objects to output mapping results
      ///
      /// - Should output one or more mapping results from the structure located at specied path
      /// - >1 result handles case of non-primitive configurations
      /// - responsible for filling in Result data structure
      /// - If 'hint' is not nullptr, use hint as 'from' config, else 'from' == 'to'
      // TODO: get rid of req_properties, and have it checked on-the-fly at point where it is needed
      map_result_inserter map(fs::path p,
                              std::vector<std::string> const &req_properties,
                              std::unique_ptr<Configuration> const &hint_config,
                              map_result_inserter result) const;

      /// Returns settings used for mapping
      const ConfigMapping::Settings &settings() const;

    private:

      /// \brief Read SimpleStructure to be imported
      SimpleStructure _make_structure(const fs::path &p) const;

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
        std::string const &report_dir,
        Log &file_log);

      static const std::string desc;
      static int run(const PrimClex &primclex, const jsonParser &kwargs, const Completer::ImportOption &import_opt);


    protected:

      /// Allow ConfigType to specialize the report formatting for 'import'
      DataFormatter<ConfigIO::Result> _import_formatter() const override;

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
        std::string const &report_dir);

      static const std::string desc;
      static int run(const PrimClex &primclex, const jsonParser &kwargs, const Completer::UpdateOption &import_opt);


    private:

      // Allow ConfigType to specialize the report formatting for 'update'
      DataFormatter<ConfigIO::Result> _update_formatter() const override;

    };

  }
}

#endif
