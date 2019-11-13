#ifndef CASM_DB_DiffTransConfigImport
#define CASM_DB_DiffTransConfigImport

#include <vector>
#include <map>
#include <string>
#include <memory>
#include "casm/global/definitions.hh"
#include "casm/database/Update.hh"
#include "casm/database/Import.hh"

namespace CASM {
  namespace Kinetics {
    class DiffTransConfigMapper;
    class DiffTransConfiguration;
  }
  namespace xtal {
    template<typename T> class BasicStructure;
    class Site;
  }
  using xtal::BasicStructure;
  using xtal::Site;;

  class PrimClex;
  class jsonParser;
}

namespace CASM {
  namespace DB {

    template<typename T> class DatabaseIterator;

    template<>
    class StructureMap<Kinetics::DiffTransConfiguration> {

    public:

      using ConfigType = Kinetics::DiffTransConfiguration;

      /// Construct with PrimClex and by moving a ConfigMapper
      StructureMap(MappingSettings const &_set,
                   std::unique_ptr<Kinetics::DiffTransConfigMapper> mapper);

      /// Construct with PrimClex and settings (see Import / Update desc)
      StructureMap(MappingSettings const &_set,
                   const PrimClex &primclex);

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
        std::unique_ptr<ConfigType> const &hint_config,
        map_result_inserter result) const;

      /// Returns JSON with settings used after combing constructor input and defaults
      const MappingSettings &settings() const {
        return m_set;
      }


    private:

      /// \brief Import Configuration with only occupation DoF
      bool _occupation_only() const;

      /// \brief Read SimpleStructure to be imported
      SimpleStructure _make_structure(const fs::path &p) const;

      MappingSettings m_set;

      std::unique_ptr<Kinetics::DiffTransConfigMapper> m_difftransconfigmapper;
    };

    /// DiffTransConfiguration-specialized Import
    template<>
    class Import<Kinetics::DiffTransConfiguration> : public ImportT<Kinetics::DiffTransConfiguration> {
    public:
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

    /// DiffTransConfiguration-specialized Update
    template<>
    class Update<Kinetics::DiffTransConfiguration> : public UpdateT<Kinetics::DiffTransConfiguration> {
    public:

      /// \brief Constructor
      Update(
        const PrimClex &primclex,
        const StructureMap<ConfigType> &mapper,
        fs::path report_dir);

      static const std::string desc;
      static int run(const PrimClex &primclex, const jsonParser &kwargs, const Completer::UpdateOption &import_opt);

      using UpdateT<ConfigType>::update;

    private:

      // Allow ConfigType to specialize the report formatting for 'update'
      DataFormatter<ConfigIO::Result> _update_formatter() const override;

    };


  }
}
#endif
