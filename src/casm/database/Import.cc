#include "casm/database/Import_impl.hh"

#include <ctime>
#include "casm/app/ProjectSettings.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/external/boost.hh"
#include "casm/database/ConfigTypeTraits.hh"

typedef std::back_insert_iterator<std::vector<CASM::fs::path> > vector_path_back_inserter;

// explicit template instantiations
#define INST_Import(r, data, type) \
template class ConfigData<type>; \
template class StructureMap<type>; \
template class Remove<type>; \

namespace CASM {
  namespace DB {

    template std::pair<vector_path_back_inserter, int>
    construct_pos_paths<vector_path_back_inserter>(
      const PrimClex &primclex,
      const Completer::ImportOption &import_opt,
      vector_path_back_inserter result);


    BOOST_PP_SEQ_FOR_EACH(INST_Import, _, CASM_DB_CONFIG_TYPES)


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

      GenericDatumFormatter<std::string, ConfigIO::Result> pos() {
        return GenericDatumFormatter<std::string, Result>("pos", "", [&](const Result & res) {
          return res.pos.string();
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
          return res.mapped_props.from;
        });
      }

      GenericDatumFormatter<std::string, ConfigIO::Result> from_configname() {
        return GenericDatumFormatter<std::string, Result>("from_configname", "", [&](const Result & res) {
          return res.mapped_props.from;
        });
      }

      GenericDatumFormatter<std::string, ConfigIO::Result> to_configname() {
        return GenericDatumFormatter<std::string, Result>("to_configname", "", [&](const Result & res) {
          return res.mapped_props.to;
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
          return data_results.find(res.mapped_props.from)->second.preexisting;
        },
        [&](const Result & res) {
          return data_results.find(res.mapped_props.from) != data_results.end();
        });
      }

      GenericDatumFormatter<bool, ConfigIO::Result> import_data(const std::map<std::string, ImportData> &data_results) {
        return GenericDatumFormatter<bool, Result>(
                 "import_data",
                 "",
        [&](const Result & res) {
          return data_results.find(res.mapped_props.from)->second.last_insert == res.pos;
        },
        [&](const Result & res) {
          return data_results.count(res.mapped_props.from) != 0;
        });
      }

      GenericDatumFormatter<bool, ConfigIO::Result> import_additional_files(const std::map<std::string, ImportData> &data_results) {
        return GenericDatumFormatter<bool, Result>(
                 "import_additional_files",
                 "",
        [&](const Result & res) {
          auto it = data_results.find(res.mapped_props.from);
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
          return res.mapped_props.mapped["lattice_deformation_cost"].get<double>();
        },
        [&](const Result & res) {
          return res.mapped_props.mapped.contains("lattice_deformation_cost");
        });
      }

      GenericDatumFormatter<double, ConfigIO::Result> basis_deformation_cost() {
        return GenericDatumFormatter<double, Result>(
                 "basis_deformation_cost", "",
        [&](const Result & res) {
          return res.mapped_props.mapped["basis_deformation_cost"].get<double>();
        },
        [&](const Result & res) {
          return res.mapped_props.mapped.contains("basis_deformation_cost");
        });
      }

      GenericDatumFormatter<double, ConfigIO::Result> relaxed_energy() {
        return GenericDatumFormatter<double, Result>(
                 "relaxed_energy", "",
        [&](const Result & res) {
          return res.mapped_props.mapped["relaxed_energy"].get<double>();
        },
        [&](const Result & res) {
          return res.mapped_props.mapped.contains("relaxed_energy");
        });
      }

      GenericDatumFormatter<double, ConfigIO::Result> score(PropertiesDatabase &db_props) {
        return GenericDatumFormatter<double, Result>(
                 "score", "",
        [&](const Result & res) {
          return db_props.score(res.mapped_props);
        },
        [&](const Result & res) {
          return res.has_data;
        });
      }

      GenericDatumFormatter<double, ConfigIO::Result> best_score(PropertiesDatabase &db_props) {
        return GenericDatumFormatter<double, Result>(
                 "score", "",
        [&](const Result & res) {
          return db_props.best_score(res.mapped_props.to);
        },
        [&](const Result & res) {
          return db_props.find_via_to(res.mapped_props.to) != db_props.end();
        });
      }

      GenericDatumFormatter<bool, ConfigIO::Result> is_best(PropertiesDatabase &db_props) {
        return GenericDatumFormatter<bool, Result>(
                 "is_best", "",
        [&](const Result & res) {
          return res.mapped_props.from == db_props.relaxed_from(res.mapped_props.to);
        },
        [&](const Result & res) {
          return db_props.find_via_to(res.mapped_props.to) != db_props.end();
        });
      }

      /// Gives a 'selected' column, set all to false
      GenericDatumFormatter<bool, ConfigIO::Result> selected() {
        return GenericDatumFormatter<bool, Result>("selected", "", [&](const Result & res) {
          return false;
        });
      }

      /// Insert default formatters to dictionary
      void default_formatters(
        DataFormatterDictionary<Result> &dict,
        PropertiesDatabase &db_props,
        const std::map<std::string, ImportData> &data_results) {

        dict.insert(
          pos(),
          fail_msg(),
          configname(),
          from_configname(),
          to_configname(),
          has_data(),
          has_complete_data(),
          preexisting_data(data_results),
          import_data(data_results),
          import_additional_files(data_results),
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
    fs::path ConfigDataGeneric::calc_properties_path(fs::path pos_path) const {

      // check 1: is a JSON file
      if(pos_path.extension() == ".json" || pos_path.extension() == ".JSON") {
        return pos_path;
      }

      // check 2: /path/to/POS -> /path/to/calctype.current/properties.calc.json
      {
        fs::path dft_path = pos_path;
        dft_path.remove_filename();
        (dft_path /= ("calctype." + primclex().settings().default_clex().calctype)) /= "properties.calc.json";
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

    /// \brief Return true if there are existing files in the traning_data directory
    ///        for a particular configuration
    bool ConfigDataGeneric::has_existing_files(const std::string &from_configname) const {
      fs::path p = calc_dir(from_configname);
      if(!fs::exists(p)) {
        return false;
      }
      return std::distance(fs::directory_iterator(p), fs::directory_iterator());
    }

    /// \brief Return true if there are existing files in the traning_data directory
    ///        for a particular configuration
    bool ConfigDataGeneric::has_existing_data(const std::string &from_configname) const {
      return db_props().find_via_from(from_configname) != db_props().end();
    }

    bool ConfigDataGeneric::has_existing_data_or_files(const std::string &from_configname) const {
      return has_existing_data(from_configname) || has_existing_files(from_configname);
    }

    /// Check if 'properties.calc.json' file has not changed since last read
    ///
    /// - Compares 'data_timestamp' && fs::last_write_time
    bool ConfigDataGeneric::no_change(const std::string &configname) const {
      fs::path prop_path = calc_properties_path(configname);
      if(!prop_path.empty()) {
        auto it = db_props().find_via_from(configname);
        if(it != db_props().end()) {
          auto json_it = it->unmapped.find("data_timestamp");
          if(json_it != it->unmapped.end() &&
             json_it->get<time_t>() == fs::last_write_time(prop_path)) {

            return true;
          }
        }
      }
      return false;
    }

    /// \brief Remove existing files in the traning_data directory for a particular
    ///        configuration
    void ConfigDataGeneric::rm_files(const std::string &configname, bool dry_run) const {
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
    void ConfigDataGeneric::cp_files(
      const fs::path &pos_path,
      const std::string &configname,
      bool dry_run,
      bool copy_additional_files) const {
      fs::path p = calc_dir(configname);
      if(!fs::exists(p)) {
        if(!dry_run) {
          fs::create_directory(p);
        }
      }

      fs::path calc_props_path = calc_properties_path(pos_path);
      if(calc_props_path.empty()) {
        return;
      }

      file_log().custom(std::string("Copy calculation files: ") + configname);
      if(!copy_additional_files) {
        file_log() << "cp " << calc_props_path << " " << p << std::endl;
        if(!dry_run) {
          fs::copy_file(calc_props_path, p);
        }
      }
      else {
        recurs_cp_files(calc_props_path.remove_filename(), p, dry_run, file_log());
      }
      file_log() << std::endl;
    }

  }
}
