#include "casm/database/Import.hh"

#include <ctime>
#include "casm/app/ProjectSettings.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/external/boost.hh"

namespace CASM {
  namespace DB {

    /// \brief Constructor
    ImportBase::ImportBase(
      const PrimClex &primclex,
      bool import_data,
      bool copy_additional_files,
      bool overwrite,
      fs::path report_dir) :

      Logging(primclex),
      m_report_dir(report_dir),
      m_file_str(m_report_dir / "files"),
      m_file_log(m_file_str),
      m_primclex(primclex),
      m_import_data(import_data),
      m_copy_additional_files(copy_additional_files),
      m_overwrite(overwrite) {}

    /// Create a new report directory to avoid overwriting existing results
    fs::path ImportBase::create_report_dir(fs::path report_dir) {
      Index i = 0;
      while(fs::exists(report_dir.string() + "." + std::to_string(i))) {
        ++i;
      }
      report_dir = report_dir.string() + "." + std::to_string(i);
      fs::create_directory(report_dir);
      return report_dir;
    }

    Database<Supercell> &ImportBase::db_supercell() const {
      return m_primclex.db<Supercell>();
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
    fs::path ImportBase::_calc_properties_path(fs::path pos_path) const {

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
    bool ImportBase::_has_existing_files(const std::string &from_configname) const {
      fs::path p = _calc_dir(from_configname);
      if(!fs::exists(p)) {
        return false;
      }
      return std::distance(fs::directory_iterator(p), fs::directory_iterator());
    }

    /// \brief Return true if there are existing files in the traning_data directory
    ///        for a particular configuration
    bool ImportBase::_has_existing_data(const std::string &from_configname) const {
      return db_props().find_via_from(from_configname) != db_props().end();
    }

    bool ImportBase::_has_existing_data_or_files(const std::string &from_configname) const {
      return _has_existing_data(from_configname) || _has_existing_files(from_configname);
    }

    /// Check if 'properties.calc.json' file has not changed since last read
    ///
    /// - Compares 'data_timestamp' && fs::last_write_time
    bool ImportBase::_no_change(const std::string &configname) const {
      fs::path prop_path = _calc_properties_path(configname);
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
    void ImportBase::_rm_files(const std::string &configname) const {
      fs::path p = _calc_dir(configname);
      if(!fs::exists(p)) {
        return;
      }
      m_file_log.custom(std::string("Remove calculation files: ") + configname);
      _recurs_rm_files(p);
      m_file_log << std::endl;
    }

    void ImportBase::_recurs_rm_files(const fs::path &p) const {
      auto it = fs::directory_iterator(p);
      auto end = fs::directory_iterator();
      for(; it != end; ++it) {
        if(fs::is_regular_file(*it)) {
          m_file_log << "rm " << *it << std::endl;
          fs::remove(*it);
        }
        else {
          _recurs_rm_files(*it);
          m_file_log << "rm " << *it << std::endl;
        }
      }
      m_file_log << "rm " << p << std::endl;
    }

    /// \brief Copy files in the same directory as properties.calc.json into the
    ///        traning_data directory for a particular configuration
    ///
    /// - First: calc_props_path = _calc_properties_path(pos_path) to get properties.calc.json location
    /// - If calc_props_path.empty(), return
    /// - else if !m_copy_additional_files copy properties.calc.json file only and return
    /// - else, recursively copy all files from calc_props_path.remove_filename()
    ///   to the training data directory for the current calctype
    void ImportBase::_cp_files(const fs::path &pos_path, const std::string &configname) const {
      fs::path p = _calc_dir(configname);
      if(!fs::exists(p)) {
        fs::create_directory(p);
      }

      fs::path calc_props_path = _calc_properties_path(pos_path);
      if(calc_props_path.empty()) {
        return;
      }

      m_file_log.custom(std::string("Copy calculation files: ") + configname);
      if(!m_copy_additional_files) {
        m_file_log << "cp " << calc_props_path << " " << p << std::endl;
        fs::copy_file(calc_props_path, p);
      }
      else {
        _recurs_cp_files(calc_props_path.remove_filename(), p);
      }
      m_file_log << std::endl;
    }

    void ImportBase::_recurs_cp_files(const fs::path &from_dir, const fs::path &to_dir) const {
      auto it = fs::directory_iterator(from_dir);
      auto end = fs::directory_iterator();
      for(; it != end; ++it) {
        if(fs::is_regular_file(*it)) {
          m_file_log << "cp " << *it << " " << to_dir << std::endl;
          fs::copy_file(*it, to_dir / it->path().filename());
        }
        else {
          fs::path new_to_dir = to_dir / it->path().filename();
          fs::create_directories(new_to_dir);
          _recurs_cp_files(*it, new_to_dir);
        }
      }
    }

    void ImportBase::_import_report(
      std::vector<Result> &results,
      const std::map<std::string, ImportData> &data_results) {

      // map_fail: could not map
      // map_success: could map
      // import_data_fail: would import but couldn't (score < best_score && !data_results.count(from))
      // import_data_conflicts: conflicts with other in import batch && preexisting
      // - pos, config, score_method, chosen?, overwrite?, import data?, import additional files?, score, best_score, is_preexisting?


      // list of structures that could not be mapped
      std::vector<Result> map_fail;

      // list of structures that could be mapped
      std::vector<Result> map_success;

      // list of structures that would be imported except preexisting data prevents it
      std::vector<Result> import_data_fail;

      for(long i = 0; i < results.size(); ++i) {
        const auto &res = results[i];
        if(res.mapped_props.to.empty()) {
          map_fail.push_back(res);
        }
        else {
          map_success.push_back(res);
          if(res.has_data && db_props().score(res.mapped_props) < db_props().best_score(res.mapped_props.to)) {
            import_data_fail.push_back(res);
          }
        }
      }

      // list of conflicts (multiple config with same 'from')
      std::map<std::string, std::vector<long> > conflict_count;
      std::vector<Result> conflict;

      for(long i = 0; i < results.size(); ++i) {
        const auto &res = results[i];
        auto it = conflict_count.find(res.mapped_props.from);
        if(it == conflict_count.end()) {
          conflict_count[res.mapped_props.from] = std::vector<long>(1, i);
        }
        else {
          it->second.push_back(i);
        }
      }
      for(const auto &val : conflict_count) {
        if(val.second.size() > 1) {
          for(const auto &i : val.second) {
            conflict.push_back(results[i]);
          }
        }
      }


      // output a 'batch' file with paths to structures that could not be imported
      if(map_fail.size()) {

        fs::path p = m_report_dir / "import_map_fail";
        fs::ofstream sout(p);

        log() << "WARNING: Could not import " << map_fail.size() << " structures." << std::endl;
        log() << "  See: " << p << " for details" << std::endl << std::endl;

        DataFormatterDictionary<Result> dict;
        dict.insert(_pos(), _fail_msg());
        auto formatter = dict.parse({"pos", "fail_msg"});
        sout << formatter(map_fail.begin(), map_fail.end());
      }

      // - pos, config, score_method, import data?, import additional files?, score, best_score, is_preexisting?
      auto formatter = this->_import_formatter(data_results);

      if(map_success.size()) {

        fs::path p = m_report_dir / "import_map_success";
        fs::ofstream sout(p);

        log() << "Successfully imported " << map_success.size() << " structures." << std::endl;
        log() << "  See: " << p << " for details" << std::endl << std::endl;

        sout << formatter(map_success.begin(), map_success.end());
      }

      if(import_data_fail.size()) {

        fs::path p = m_report_dir / "import_data_fail";
        fs::ofstream sout(p);

        log() << "WARNING: Did not import data from "
              << import_data_fail.size() << " structures which have are a mapping score"
              " better than the existing data." << std::endl;
        log() << "  You may wish to inspect these structures and allow overwriting "
              "or remove existing data manually." << std::endl;
        log() << "  See: " << p << " for details" << std::endl << std::endl;

        sout << formatter(import_data_fail.begin(), import_data_fail.end());
      }

      if(conflict.size()) {
        fs::path p = m_report_dir / "import_conflict";
        fs::ofstream sout(p);

        log() << "WARNING: Imported data from structures that mapped to the same configuration." << std::endl
              << "  Data can only be imported from one of the conflicting structures." << std::endl
              << "  Based on the current conflict resolution method the 'best' result was automatically chosen, " << std::endl
              << "  but you may wish to inspect these results and manually select which structures to import." << std::endl;
        log() << "  See: " << p << " for details" << std::endl << std::endl;

        sout << formatter(conflict.begin(), conflict.end());
      }
    }


    GenericDatumFormatter<std::string, ImportBase::Result> ImportBase::_pos() const {
      return GenericDatumFormatter<std::string, Result>("pos", "", [&](const Result & res) {
        return res.pos.string();
      });
    }

    GenericDatumFormatter<std::string, ImportBase::Result> ImportBase::_fail_msg() const {
      return GenericDatumFormatter<std::string, Result>("fail_msg", "", [&](const Result & res) {
        return res.fail_msg;
      });
    }

    /// Use 'from_configname' as 'configname'
    GenericDatumFormatter<std::string, ImportBase::Result> ImportBase::_configname() const {
      return GenericDatumFormatter<std::string, Result>("configname", "", [&](const Result & res) {
        return res.mapped_props.from;
      });
    }

    GenericDatumFormatter<std::string, ImportBase::Result> ImportBase::_from_configname() const {
      return GenericDatumFormatter<std::string, Result>("from_configname", "", [&](const Result & res) {
        return res.mapped_props.from;
      });
    }

    GenericDatumFormatter<std::string, ImportBase::Result> ImportBase::_to_configname() const {
      return GenericDatumFormatter<std::string, Result>("to_configname", "", [&](const Result & res) {
        return res.mapped_props.to;
      });
    }

    GenericDatumFormatter<bool, ImportBase::Result> ImportBase::_has_data() const {
      return GenericDatumFormatter<bool, Result>("has_data", "", [&](const Result & res) {
        return res.has_data;
      });
    }

    GenericDatumFormatter<bool, ImportBase::Result> ImportBase::_has_complete_data() const {
      return GenericDatumFormatter<bool, Result>("has_complete_data", "", [&](const Result & res) {
        return res.has_data;
      });
    }

    GenericDatumFormatter<bool, ImportBase::Result> ImportBase::_preexisting_data(const std::map<std::string, ImportData> &data_results) const {
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

    GenericDatumFormatter<bool, ImportBase::Result> ImportBase::_import_data(const std::map<std::string, ImportData> &data_results) const {
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

    GenericDatumFormatter<bool, ImportBase::Result> ImportBase::_import_additional_files(const std::map<std::string, ImportData> &data_results) const {
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

    GenericDatumFormatter<double, ImportBase::Result> ImportBase::_lattice_deformation_cost() const {
      return GenericDatumFormatter<double, Result>(
               "lattice_deformation_cost", "",
      [&](const Result & res) {
        return res.mapped_props.mapped["lattice_deformation_cost"].get<double>();
      },
      [&](const Result & res) {
        return res.mapped_props.mapped.contains("lattice_deformation_cost");
      });
    }

    GenericDatumFormatter<double, ImportBase::Result> ImportBase::_basis_deformation_cost() const {
      return GenericDatumFormatter<double, Result>(
               "basis_deformation_cost", "",
      [&](const Result & res) {
        return res.mapped_props.mapped["basis_deformation_cost"].get<double>();
      },
      [&](const Result & res) {
        return res.mapped_props.mapped.contains("basis_deformation_cost");
      });
    }

    GenericDatumFormatter<double, ImportBase::Result> ImportBase::_relaxed_energy() const {
      return GenericDatumFormatter<double, Result>(
               "relaxed_energy", "",
      [&](const Result & res) {
        return res.mapped_props.mapped["relaxed_energy"].get<double>();
      },
      [&](const Result & res) {
        return res.mapped_props.mapped.contains("relaxed_energy");
      });
    }

    GenericDatumFormatter<double, ImportBase::Result> ImportBase::_score() const {
      return GenericDatumFormatter<double, Result>(
               "score", "",
      [&](const Result & res) {
        return this->db_props().score(res.mapped_props);
      },
      [&](const Result & res) {
        return res.has_data;
      });
    }

    GenericDatumFormatter<double, ImportBase::Result> ImportBase::_best_score() const {
      return GenericDatumFormatter<double, Result>(
               "score", "",
      [&](const Result & res) {
        return this->db_props().best_score(res.mapped_props.to);
      },
      [&](const Result & res) {
        return this->db_props().find_via_to(res.mapped_props.to) != this->db_props().end();
      });
    }

    GenericDatumFormatter<bool, ImportBase::Result> ImportBase::_is_best() const {
      return GenericDatumFormatter<bool, Result>(
               "is_best", "",
      [&](const Result & res) {
        return res.mapped_props.from == this->db_props().relaxed_from(res.mapped_props.to);
      },
      [&](const Result & res) {
        return this->db_props().find_via_to(res.mapped_props.to) != this->db_props().end();
      });
    }

    /// Gives a 'selected' column, set all to false
    GenericDatumFormatter<bool, ImportBase::Result> ImportBase::_selected() const {
      return GenericDatumFormatter<bool, Result>("selected", "", [&](const Result & res) {
        return false;
      });
    }

    /// Insert default formatters to dictionary
    void ImportBase::_default_formatters(
      DataFormatterDictionary<Result> &dict,
      const std::map<std::string, ImportData> &data_results) const {

      dict.insert(
        _pos(),
        _fail_msg(),
        _configname(),
        _from_configname(),
        _to_configname(),
        _has_data(),
        _has_complete_data(),
        _preexisting_data(data_results),
        _import_data(data_results),
        _import_additional_files(data_results),
        _lattice_deformation_cost(),
        _basis_deformation_cost(),
        _relaxed_energy(),
        _score(),
        _best_score(),
        _is_best(),
        _selected());
    }

  }
}
