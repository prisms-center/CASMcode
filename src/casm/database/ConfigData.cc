#include "casm/database/ConfigData_impl.hh"

#include <ctime>
#include "casm/app/DirectoryStructure.hh"
#include "casm/app/ProjectSettings.hh"
#include "casm/database/PropertiesDatabase.hh"
#include "casm/database/DatabaseTypes_impl.hh"

namespace CASM {
  namespace DB {

    jsonParser &to_json(MappingSettings const &_set, jsonParser &_json) {
      _json["lattice_weight"] = _set.lattice_weight;
      _json["ideal"] = _set.ideal;
      _json["strict"] = _set.strict;
      _json["primitive_only"] = _set.primitive_only;
      if(!_set.forced_lattices.empty())
        _json["forced_lattices"] = _set.forced_lattices;
      if(!_set.filter.empty())
        _json["filter"] = _set.filter;
      _json["cost_tol"] = _set.cost_tol;
      _json["min_va_frac"] = _set.min_va_frac;
      _json["max_va_frac"] = _set.max_va_frac;
      _json["max_vol_change"] = _set.max_vol_change;

      return _json;
    }

    jsonParser const &from_json(MappingSettings &_set, jsonParser const &_json) {
      _set.set_default();

      if(_json.contains("lattice_weight"))
        _set.lattice_weight = _json["lattice_weight"].get<double>();

      if(_json.contains("ideal"))
        _set.ideal = _json["ideal"].get<bool>();

      if(_json.contains("strict"))
        _set.strict = _json["strict"].get<bool>();

      if(_json.contains("primitive_only"))
        _set.primitive_only = _json["primitive_only"].get<bool>();

      if(_json.contains("forced_lattices"))
        _set.forced_lattices = _json["forced_lattices"].get<std::vector<std::string> >();

      if(_json.contains("filter"))
        _set.filter = _json["filter"].get<std::string>();

      if(_json.contains("cost_tol"))
        _set.cost_tol = _json["cost_tol"].get<double>();

      if(_json.contains("min_va_frac"))
        _set.min_va_frac = _json["min_va_frac"].get<double>();

      if(_json.contains("max_va_frac"))
        _set.max_va_frac = _json["max_va_frac"].get<double>();

      if(_json.contains("max_vol_change"))
        _set.max_vol_change = _json["max_vol_change"].get<double>();

      return _json;
    }

    /// Create a new report directory to avoid overwriting existing results
    fs::path create_report_dir(fs::path report_dir) {
      Index i = 0;
      while(fs::exists(report_dir.string() + "." + std::to_string(i))) {
        ++i;
      }
      report_dir = report_dir.string() + "." + std::to_string(i);
      fs::create_directory(report_dir);
      return report_dir;
    }

    namespace ConfigIO {

      GenericDatumFormatter<std::string, ConfigIO::Result> path() {
        return GenericDatumFormatter<std::string, Result>("path", "", [&](const Result & res) {
          return res.properties.file_data.path();
        });
      }

      GenericDatumFormatter<std::string, ConfigIO::Result> fail_msg() {
        return GenericDatumFormatter<std::string, Result>("fail_msg", "", [&](const Result & res) {
          return res.fail_msg;
        });
      }

      /// Use 'from_configname' as 'configname'
      GenericDatumFormatter<std::string, ConfigIO::Result> configname() {
        return GenericDatumFormatter<std::string, Result>("configname", "", [&](const Result & res) {
          return res.properties.from;
        });
      }

      GenericDatumFormatter<std::string, ConfigIO::Result> from_configname() {
        return GenericDatumFormatter<std::string, Result>("from_configname", "", [&](const Result & res) {
          return res.properties.from;
        });
      }

      GenericDatumFormatter<std::string, ConfigIO::Result> to_configname() {
        return GenericDatumFormatter<std::string, Result>("to_configname", "", [&](const Result & res) {
          return res.properties.to;
        });
      }

      GenericDatumFormatter<bool, ConfigIO::Result> has_data() {
        return GenericDatumFormatter<bool, Result>("has_data", "", [&](const Result & res) {
          return res.has_data;
        });
      }

      GenericDatumFormatter<bool, ConfigIO::Result> has_complete_data() {
        return GenericDatumFormatter<bool, Result>("has_complete_data", "", [&](const Result & res) {
          return res.has_data;
        });
      }

      GenericDatumFormatter<bool, ConfigIO::Result> preexisting_data(const std::map<std::string, ImportData> &data_results) {
        return GenericDatumFormatter<bool, Result>(
                 "preexisting_data",
                 "",
        [&](const Result & res) {
          return data_results.find(res.properties.from)->second.preexisting;
        },
        [&](const Result & res) {
          return data_results.find(res.properties.from) != data_results.end();
        });
      }

      GenericDatumFormatter<bool, ConfigIO::Result> import_data(const std::map<std::string, ImportData> &data_results) {
        return GenericDatumFormatter<bool, Result>(
                 "import_data",
                 "",
        [&](const Result & res) {
          return data_results.find(res.properties.from)->second.last_insert == res.properties.file_data.path();
        },
        [&](const Result & res) {
          return data_results.count(res.properties.from) != 0;
        });
      }

      GenericDatumFormatter<bool, ConfigIO::Result> import_additional_files(const std::map<std::string, ImportData> &data_results) {
        return GenericDatumFormatter<bool, Result>(
                 "import_additional_files",
                 "",
        [&](const Result & res) {
          auto it = data_results.find(res.properties.from);
          if(it != data_results.end()) {
            return it->second.copy_more;
          }
          return false;
        });
      }

      GenericDatumFormatter<double, ConfigIO::Result> lattice_deformation_cost() {
        return GenericDatumFormatter<double, Result>(
                 "lattice_deformation_cost", "",
        [&](const Result & res) {
          return res.properties.scalar("lattice_deformation_cost");
        },
        [&](const Result & res) {
          return res.properties.has_scalar("lattice_deformation_cost");
        });
      }

      GenericDatumFormatter<double, ConfigIO::Result> basis_deformation_cost() {
        return GenericDatumFormatter<double, Result>(
                 "basis_deformation_cost", "",
        [&](const Result & res) {
          return res.properties.scalar("basis_deformation_cost");
        },
        [&](const Result & res) {
          return res.properties.has_scalar("basis_deformation_cost");
        });
      }

      GenericDatumFormatter<double, ConfigIO::Result> relaxed_energy() {
        return GenericDatumFormatter<double, Result>(
                 "relaxed_energy", "",
        [&](const Result & res) {
          return res.properties.scalar("relaxed_energy");
        },
        [&](const Result & res) {
          return res.properties.has_scalar("relaxed_energy");
        });
      }

      GenericDatumFormatter<double, ConfigIO::Result> score(PropertiesDatabase &db_props) {
        return GenericDatumFormatter<double, Result>(
                 "score", "",
        [&](const Result & res) {
          return db_props.score(res.properties);
        },
        [&](const Result & res) {
          return res.has_data;
        });
      }

      GenericDatumFormatter<double, ConfigIO::Result> best_score(PropertiesDatabase &db_props) {
        return GenericDatumFormatter<double, Result>(
                 "best_score", "",
        [&](const Result & res) {
          return db_props.best_score(res.properties.to);
        },
        [&](const Result & res) {
          return db_props.find_via_to(res.properties.to) != db_props.end();
        });
      }

      GenericDatumFormatter<bool, ConfigIO::Result> is_best(PropertiesDatabase &db_props) {
        return GenericDatumFormatter<bool, Result>(
                 "is_best", "",
        [&](const Result & res) {
          return res.properties.from == db_props.relaxed_from(res.properties.to);
        },
        [&](const Result & res) {
          return db_props.find_via_to(res.properties.to) != db_props.end();
        });
      }

      /// Gives a 'selected' column, set all to false
      GenericDatumFormatter<bool, ConfigIO::Result> selected() {
        return GenericDatumFormatter<bool, Result>("selected", "", [&](const Result & res) {
          return false;
        });
      }

      /// Insert default formatters to dictionary, for 'casm import'
      void default_import_formatters(
        DataFormatterDictionary<Result> &dict,
        PropertiesDatabase &db_props,
        const std::map<std::string, ImportData> &data_results) {

        default_update_formatters(dict, db_props);

        dict.insert(
          preexisting_data(data_results),
          import_data(data_results),
          import_additional_files(data_results));
      }

      /// Insert default formatters to dictionary, for 'casm update'
      void default_update_formatters(
        DataFormatterDictionary<Result> &dict,
        PropertiesDatabase &db_props) {

        dict.insert(
          path(),
          fail_msg(),
          configname(),
          from_configname(),
          to_configname(),
          has_data(),
          has_complete_data(),
          lattice_deformation_cost(),
          basis_deformation_cost(),
          relaxed_energy(),
          score(db_props),
          best_score(db_props),
          is_best(db_props),
          selected());
      }
    }

    /// \brief Return path to properties.calc.json that will be imported
    ///        checking a couple possible locations relative to pos_path
    ///
    /// checks:
    /// 1) is a JSON file? is pos_path ends in ".json" or ".JSON", return pos_path
    /// 2) assume pos_path is /path/to/POS, checks for /path/to/calctype.current/properties.calc.json
    /// 3) assume pos_path is /path/to/POS, checks for /path/to/properties.calc.json
    /// else returns empty path
    ///
    fs::path ConfigData::calc_properties_path(fs::path pos_path, PrimClex const &_pclex) {

      // check 1: is a JSON file
      if(pos_path.extension() == ".json" || pos_path.extension() == ".JSON") {
        return pos_path;
      }

      // check 2: /path/to/POS -> /path/to/calctype.current/properties.calc.json
      {
        fs::path dft_path = pos_path;
        dft_path.remove_filename();
        (dft_path /= ("calctype." + _pclex.settings().default_clex().calctype)) /= "properties.calc.json";
        if(fs::exists(dft_path)) {
          return dft_path;
        }
      }

      // check 3: /path/to/POS -> /path/to/properties.calc.json
      {
        fs::path dft_path = pos_path;
        dft_path.remove_filename();
        dft_path /= "properties.calc.json";
        if(fs::exists(dft_path)) {
          return dft_path;
        }
      }

      // not found, return empty path
      return fs::path();
    }

    // If pos_path can be used to resolve a properties.calc.json, return its path.
    // Otherwise return pos_path
    fs::path ConfigData::resolve_struc_path(fs::path pos_path, PrimClex const &_pclex) {
      fs::path p = ConfigData::calc_properties_path(pos_path, _pclex);
      if(!p.empty())
        pos_path = p;
      return pos_path;
    }

    /// \brief Path to default calctype training_data directory for config
    Database<Supercell> &ConfigData::db_supercell() const {
      return primclex().db<Supercell>();
    }

    PropertiesDatabase &ConfigData::db_props() const {
      return m_db_props_func();
    }

    /// \brief Path to default calctype training_data directory for config
    fs::path ConfigData::calc_dir(const std::string configname) const {
      return primclex().dir().configuration_calc_dir(configname,
                                                     primclex().settings().default_clex().calctype);
    }

    /// \brief Return true if there are existing files in the traning_data directory
    ///        for a particular configuration
    bool ConfigData::has_existing_files(const std::string &from_configname) const {
      fs::path p = calc_dir(from_configname);
      if(!fs::exists(p)) {
        return false;
      }
      return std::distance(fs::directory_iterator(p), fs::directory_iterator());
    }

    /// \brief Return true if there are existing files in the traning_data directory
    ///        for a particular configuration
    bool ConfigData::has_existing_data(const std::string &from_configname) const {
      return db_props().find_via_from(from_configname) != db_props().end();
    }

    bool ConfigData::has_existing_data_or_files(const std::string &from_configname) const {
      return has_existing_data(from_configname) || has_existing_files(from_configname);
    }

    /// Check if 'properties.calc.json' file has not changed since last read
    ///
    /// - Compares 'data_timestamp' && fs::last_write_time
    bool ConfigData::no_change(const std::string &configname) const {
      fs::path prop_path = calc_properties_path(configname, primclex());
      if(!prop_path.empty()) {
        auto it = db_props().find_via_from(configname);
        if(it != db_props().end()) {
          return it->file_data == FileData(prop_path);
        }
      }
      return false;
    }

    /// \brief Remove existing files in the traning_data directory for a particular
    ///        configuration
    void ConfigData::rm_files(const std::string &configname, bool dry_run) const {
      fs::path p = calc_dir(configname);
      if(!fs::exists(p)) {
        return;
      }
      file_log().custom(std::string("Remove calculation files: ") + configname);
      recurs_rm_files(p, dry_run, file_log());
      file_log() << std::endl;
    }

    /// \brief Copy files in the same directory as properties.calc.json into the
    ///        traning_data directory for a particular configuration
    ///
    /// - First: calc_props_path = _calc_properties_path(pos_path) to get properties.calc.json location
    /// - If calc_props_path.empty(), return
    /// - else if !copy_additional_files copy properties.calc.json file only and return
    /// - else, recursively copy all files from calc_props_path.remove_filename()
    ///   to the training data directory for the current calctype
    ///
    /// \returns {did_cp, did_cp_more}:
    /// - did_cp: if properties.calc.json file was found and copied
    /// - did_cp_more: if additional files were found and copied
    ///
    std::pair<bool, bool> ConfigData::cp_files(
      const fs::path &pos_path,
      const std::string &configname,
      bool dry_run,
      bool copy_additional_files) const {

      bool did_cp(false);
      bool did_cp_more(false);

      fs::path p = calc_dir(configname);
      if(!fs::exists(p)) {
        if(!dry_run) {
          fs::create_directories(p);
        }
      }

      fs::path calc_props_path = calc_properties_path(pos_path, primclex());
      if(calc_props_path.empty()) {
        return std::make_pair(did_cp, did_cp_more);
      }

      file_log().custom(std::string("Copy calculation files: ") + configname);
      if(!copy_additional_files) {
        file_log() << "cp " << calc_props_path << " " << p / "properties.calc.json" << std::endl;
        did_cp = true;
        if(!dry_run) {
          fs::copy_file(calc_props_path, p / "properties.calc.json");
        }
      }
      else {
        Index count = recurs_cp_files(calc_props_path.remove_filename(), p, dry_run, file_log());
        if(count) {
          did_cp_more = true;
        }
      }
      file_log() << std::endl;
      return std::make_pair(did_cp, did_cp_more);
    }

  }
}
