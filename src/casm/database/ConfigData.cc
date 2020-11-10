#include "casm/database/ConfigData_impl.hh"

#include <boost/filesystem.hpp>
#include "casm/app/DirectoryStructure.hh"
#include "casm/app/ProjectSettings.hh"
#include "casm/database/PropertiesDatabase.hh"
#include "casm/database/DatabaseTypes_impl.hh"

namespace CASM {
  namespace Local {
    /// \brief Return path to properties.calc.json that will be imported
    ///        checking a couple possible locations relative to pos_path
    ///
    /// checks:
    /// 1) is a JSON file? is pos_path ends in ".json" or ".JSON", return pos_path
    /// 2) assume pos_path is /path/to/POS, checks for /path/to/calctype.current/properties.calc.json
    /// 3) assume pos_path is /path/to/POS, checks for /path/to/properties.calc.json
    /// else returns empty path
    ///
    static std::string _resolve_properties_path(std::string pos_path, PrimClex const &_pclex) {

      // check 1: is a JSON file
      if(pos_path.find(".json") != std::string::npos || pos_path.find(".JSON") != std::string::npos) {
        return pos_path;
      }

      // check 2: /path/to/POS -> /path/to/calctype.current/properties.calc.json
      {
        fs::path dft_path = pos_path;
        dft_path.remove_filename();
        (dft_path /= ("calctype." + _pclex.settings().default_clex().calctype)) /= "properties.calc.json";
        if(fs::exists(dft_path)) {
          return dft_path.string();
        }
      }

      // check 3: /path/to/POS -> /path/to/properties.calc.json
      {
        fs::path dft_path = pos_path;
        dft_path.remove_filename();
        dft_path /= "properties.calc.json";
        if(fs::exists(dft_path)) {
          return dft_path.string();
        }
      }

      // not found, return empty path
      return "";
    }
  }

  namespace DB {

    /// Create a new report directory to avoid overwriting existing results
    std::string create_report_dir(std::string report_dir) {
      Index i = 0;
      while(fs::exists(report_dir + "." + std::to_string(i))) {
        ++i;
      }
      report_dir += ("." + std::to_string(i));
      fs::create_directory(report_dir);
      return report_dir;
    }

    namespace ConfigIO {

      GenericDatumFormatter<std::string, ConfigIO::Result> initial_path() {
        return GenericDatumFormatter<std::string, Result>("initial_path", "",
        [](const Result & res) {
          return res.pos_path;
        });
      }

      GenericDatumFormatter<std::string, ConfigIO::Result> final_path() {
        return GenericDatumFormatter<std::string, Result>("final_path", "",
        [](const Result & res) {
          return res.properties.file_data.path();
        });
      }

      GenericDatumFormatter<std::string, ConfigIO::Result> fail_msg() {
        return GenericDatumFormatter<std::string, Result>("fail_msg", "",
        [](const Result & res) {
          return res.fail_msg;
        });
      }

      GenericDatumFormatter<std::string, ConfigIO::Result> data_origin() {
        return GenericDatumFormatter<std::string, Result>("data_origin", "",
        [](const Result & res) {
          return res.properties.origin;
        });
      }

      GenericDatumFormatter<std::string, ConfigIO::Result> to_configname() {
        return GenericDatumFormatter<std::string, Result>("to_configname", "",
        [](const Result & res) {
          return res.properties.to;
        });
      }

      GenericDatumFormatter<bool, ConfigIO::Result> has_data() {
        return GenericDatumFormatter<bool, Result>("has_data", "",
        [](const Result & res) {
          return res.has_data;
        });
      }

      GenericDatumFormatter<bool, ConfigIO::Result> has_complete_data() {
        return GenericDatumFormatter<bool, Result>("has_complete_data", "",
        [](const Result & res) {
          return res.has_data;
        });
      }

      GenericDatumFormatter<bool, ConfigIO::Result> preexisting_data() {
        return GenericDatumFormatter<bool, Result>(
                 "preexisting_data",
                 "",
        [&](const Result & res) {
          return res.import_data.preexisting;
        });
      }

      GenericDatumFormatter<bool, ConfigIO::Result> import_data() {
        return GenericDatumFormatter<bool, Result>(
                 "import_data",
                 "",
        [&](const Result & res) {
          return  res.import_data.is_best;// && it->second.import_data;
        });
      }

      GenericDatumFormatter<bool, ConfigIO::Result> import_additional_files() {
        return GenericDatumFormatter<bool, Result>(
                 "import_additional_files",
                 "",
        [&](const Result & res) {
          return res.import_data.copy_more;
        });
      }

      GenericDatumFormatter<double, ConfigIO::Result> lattice_deformation_cost() {
        return GenericDatumFormatter<double, Result>(
                 "lattice_deformation_cost", "",
        [](const Result & res) {
          return res.properties.scalar("lattice_deformation_cost");
        },
        [](const Result & res) {
          return res.properties.has_scalar("lattice_deformation_cost");
        });
      }

      GenericDatumFormatter<double, ConfigIO::Result> atomic_deformation_cost() {
        return GenericDatumFormatter<double, Result>(
                 "atomic_deformation_cost", "",
        [](const Result & res) {
          return res.properties.scalar("atomic_deformation_cost");
        },
        [](const Result & res) {
          return res.properties.has_scalar("atomic_deformation_cost");
        });
      }

      GenericDatumFormatter<double, ConfigIO::Result> relaxed_energy() {
        return GenericDatumFormatter<double, Result>(
                 "relaxed_energy", "",
        [](const Result & res) {
          return res.properties.scalar("relaxed_energy");
        },
        [](const Result & res) {
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
          return res.properties.origin == db_props.find_via_to(res.properties.to)->origin;
        },
        [&](const Result & res) {
          return db_props.find_via_to(res.properties.to) != db_props.end();
        });
      }

      /// Gives a 'selected' column, set all to false
      GenericDatumFormatter<bool, ConfigIO::Result> selected() {
        return GenericDatumFormatter<bool, Result>("selected", "",
        [](const Result & res) {
          return false;
        });
      }

      /// Gives a 'selected' column, set all to false
      GenericDatumFormatter<bool, ConfigIO::Result> is_new_config() {
        return GenericDatumFormatter<bool, Result>("is_new_config", "",
        [](const Result & res) {
          return res.is_new_config;
        });
      }

      /// Insert default formatters to dictionary, for 'casm import'
      void default_import_formatters(
        DataFormatterDictionary<Result> &dict,
        PropertiesDatabase &db_props) {

        default_update_formatters(dict, db_props);

        dict.insert(
          preexisting_data(),
          import_data(),
          import_additional_files());
      }

      /// Insert default formatters to dictionary, for 'casm update'
      void default_update_formatters(
        DataFormatterDictionary<Result> &dict,
        PropertiesDatabase &db_props) {

        dict.insert(
          initial_path(),
          final_path(),
          fail_msg(),
          //configname(),
          data_origin(),
          to_configname(),
          is_new_config(),
          has_data(),
          has_complete_data(),
          lattice_deformation_cost(),
          atomic_deformation_cost(),
          relaxed_energy(),
          score(db_props),
          best_score(db_props),
          is_best(db_props),
          selected());
      }
    }


    // If pos_path can be used to resolve a properties.calc.json, return its path.
    // Otherwise return pos_path
    std::string ConfigData::resolve_struc_path(std::string pos_path, PrimClex const &_pclex) {
      std::string p = Local::_resolve_properties_path(pos_path, _pclex);
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
    std::string ConfigData::calc_dir(const std::string configname) const {
      return primclex().dir().configuration_calc_dir(configname,
                                                     primclex().settings().default_clex().calctype).string();
    }

    /// \brief Return true if there are existing files in the traning_data directory
    ///        for a particular configuration
    bool ConfigData::has_existing_files(const std::string &to_configname) const {
      fs::path p = calc_dir(to_configname);
      if(!fs::exists(p)) {
        return false;
      }
      return std::distance(fs::directory_iterator(p), fs::directory_iterator());
    }

    /// \brief Return true if there is data already associated with a particular configuration
    bool ConfigData::has_existing_data(const std::string &to_configname) const {
      return db_props().find_via_to(to_configname) != db_props().end();
    }

    bool ConfigData::has_existing_data_or_files(const std::string &configname) const {
      return has_existing_data(configname) || has_existing_files(configname);
    }

    /// Check if 'properties.calc.json' file has not changed since last read
    ///
    /// - Compares 'data_timestamp' && fs::last_write_time
    bool ConfigData::no_change(const std::string &configname) const {
      FileData tfile_data(calc_properties_path(primclex(), configname));

      auto it = db_props().find_via_origin(tfile_data.path());
      if(it != db_props().end()) {
        return it->file_data == tfile_data;
      }

      return !tfile_data.exists();
    }

    /// \brief Remove existing files in the traning_data directory for a particular
    ///        configuration
    void ConfigData::rm_files(const std::string &configname, bool dry_run) const {
      fs::path p = calc_dir(configname);
      if(!fs::exists(p)) {
        return;
      }
      log().custom(std::string("Remove calculation files: ") + configname);
      recurs_rm_files(p, dry_run, log());
      log() << std::endl;
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
    void ConfigData::cp_files(
      ConfigIO::Result &res,
      bool dry_run,
      bool copy_additional_files) const {

      res.import_data.copy_data = false;
      res.import_data.copy_more = false;

      fs::path p = calc_dir(res.properties.to);
      if(!fs::exists(p)) {
        if(!dry_run) {
          fs::create_directories(p);
        }
      }

      fs::path origin_props_path = Local::_resolve_properties_path(res.properties.file_data.path(), primclex());
      if(origin_props_path.empty()) {
        return;
      }

      log().custom(std::string("Copy calculation files: ") + res.properties.to);
      if(!copy_additional_files) {
        log() << "cp " << origin_props_path << " " << p / "properties.calc.json" << std::endl;
        res.import_data.copy_data = true;
        if(!dry_run) {
          fs::copy_file(origin_props_path, p / "properties.calc.json");
        }
      }
      else {
        Index count = recurs_cp_files(origin_props_path.remove_filename(), p, dry_run, log());
        if(count) {
          res.import_data.copy_more = true;
        }
      }
      res.properties.file_data = FileData((p / "properties.calc.json").string());
      log() << std::endl;
      return;
    }

  }
}
