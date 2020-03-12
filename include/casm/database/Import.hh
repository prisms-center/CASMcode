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


    /// \brief Struct with optional parameters for Config/Data Import
    /// Specifies default parameters for all values, in order to simplify
    /// parsing from JSON
    struct ImportSettings {
      ImportSettings(bool _import = false,
                     bool _copy_files = false,
                     bool _additional_files = false,
                     bool _overwrite = false) :
        import(_import),
        copy_files(_copy_files),
        additional_files(_additional_files),
        overwrite(_overwrite) {}

      void set_default() {
        *this = ImportSettings();
      }

      // attempt to import calculation results into database, else just insert
      // configurations w/out data
      bool import;

      // attempt to copy properties.calc.json file to training directory of best configuration
      // else just import data
      bool copy_files;

      // attempt to copy extra files from the directory where the structure is
      // being imported from to the training_data directory
      bool additional_files;

      // Allow overwriting of existing data by 'casm import'
      bool overwrite;
    };

    jsonParser &to_json(ImportSettings const &_set, jsonParser &_json);

    jsonParser const &from_json(ImportSettings &_set, jsonParser const &_json);

    /// Generic ConfigType-dependent part of Import
    template<typename _ConfigType>
    class ImportT : protected ConfigData {

    public:

      using ConfigType = _ConfigType;

      /// \brief Constructor
      ImportT(
        const PrimClex &primclex,
        const StructureMap<ConfigType> &mapper,
        ImportSettings const  &_set,
        std::string report_dir,
        Log &_file_log);

      template<typename PathIterator>
      void import(PathIterator begin, PathIterator end);

      ImportSettings const &settings() const {
        return m_set;
      }

    protected:
      /// Allow ConfigType to specialize the report formatting for 'import'
      virtual DataFormatter<ConfigIO::Result> _import_formatter() const = 0;

    private:
      void _import_report(
        std::vector<ConfigIO::Result> &results)const;

      // Copies file for each element of results that corresponds to a succesful optimal mapping
      void _copy_files(std::vector<ConfigIO::Result> &results) const;

      const StructureMap<ConfigType> &m_structure_mapper;

      ImportSettings m_set;

      std::string m_report_dir;
    };


    // To be specialized for ConfigType (no default implemenation exists)
    template<typename ConfigType>
    class Import;// : public ImportT<ConfigType>;
    /*{
    public:
      static const std::string desc;
      int run(const PrimClex &, const jsonParser &input, const Completer::ImportOption &opt);

    private:
      // Allow ConfigType to specialize the report formatting for 'import'
      DataFormatter<ConfigIO::Result> _import_formatter() const override;
        };*/

  }
}

#endif
