#ifndef CASM_DB_Import_impl
#define CASM_DB_Import_impl

#include "casm/database/Import.hh"
#include "casm/app/DirectoryStructure.hh"
#include "casm/app/import.hh"
#include "casm/clex/PrimClex.hh"

namespace CASM {
  namespace DB {

    /// Construct pos_paths from input args
    ///
    /// \return pair of OutputIterator, returncode
    ///
    template<typename OutputIterator>
    std::pair<OutputIterator, int> construct_pos_paths(
      const PrimClex &primclex,
      const Completer::ImportOption &import_opt,
      OutputIterator result) {

      const po::variables_map &vm = import_opt.vm();
      Index count = 0;
      bool missing_files(false);

      // count files, check if existing, print warning messages for missing files
      auto lambda = [&](const fs::path & p) {
        ++count;
        if(!fs::exists(p)) {
          if(!missing_files) {
            primclex.err_log() << "*** ERROR: Missing file(s):\n";
          }
          missing_files = true;
          primclex.err_log()   << "           " << fs::absolute(p) << "\n";
        }
        *result++ = p;
        return true;
      };

      //read all --pos paths
      std::for_each(import_opt.pos_vec().begin(), import_opt.pos_vec().end(), lambda);

      //read all the paths from a batch file
      if(vm.count("batch")) {

        if(!fs::exists(import_opt.batch_path())) {
          primclex.err_log() << "ERROR: Batch import file does not exist at " << import_opt.batch_path() << "\n";
          return std::make_pair<result, ERR_MISSING_INPUT_FILE>;
        }

        fs::ifstream batchfile(import_opt.batch_path());
        fs::path tpath;
        while(batchfile >> tpath) {
          lambda(tpath);
          ++count;
          batchfile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        }
        if(tpath != pos_paths.back() && tpath.string().size() != 0 && fs::exists(tpath)) {
          lambda(tpath);
        }
      }

      if(count == 0) {
        primclex.err_log() <<   "ERROR: No files specified for import.\n";
        if(vm.count("batch")) {
          primclex.err_log() << "       Check batch file for errors.\n";
        }
        return std::make_pair<result, ERR_INVALID_INPUT_FILE>;
      }

      if(missing_files) {
        return std::make_pair<result, ERR_INVALID_INPUT_FILE>;
      }

      return std::make_pair<result, 0>;
    }

    // --- class ConfigData ---

    template<typename _ConfigType>
    Database<_ConfigType> &ConfigData<_ConfigType>::db_config() const {
      return primclex().template db<ConfigType>();
    }

    template<typename _ConfigType>
    PropertiesDatabase &ConfigData<_ConfigType>::db_props() const {
      return primclex().template db_props<ConfigType>();
    }

    /// \brief Path to default calctype training_data directory for config
    template<typename _ConfigType>
    fs::path ConfigData<_ConfigType>::calc_dir(const std::string configname) const {
      return primclex().dir().template calc_dir<ConfigType>(configname);
    }


    // --- ImportT ---

    template<typename _ConfigType>
    template<typename PathIterator>
    void ImportT<_ConfigType>::import(PathIterator begin, PathIterator end) {

      // vector of Mapping results, may be >1 per input if primitive and non-primitive
      std::vector<ConfigIO::Result> results;

      // map of data import results
      //   'configname' -> 'preexisting?, copy_data?, copy_more?, last_import
      std::map<std::string, ConfigIO::ImportData> data_results;

      auto it = begin;
      for(; it != end; ++it) {

        std::vector<ConfigIO::Result> tvec;

        // Outputs one or more mapping results from the structure located at specied path
        //   See _import documentation for more.
        m_structure_mapper.map(*it, this->db_config().end(), std::back_inserter(tvec));

        for(auto &res : tvec) {

          results.push_back(res);

          // reasons to import data or not:
          //
          // could not map || no data || !m_import_data:
          //   do not import
          // !_has_existing_data_or_files:
          //   import
          // _has_existing_data_or_files:
          //   new_data:
          //     better_score:
          //       import
          //     !better_score:
          //       do not import
          //   !new_data:
          //      better_score:
          //        overwrite:
          //          import
          //        !overwrite:
          //          do not import
          //      !better_score:
          //        do not import

          // if could not map, no data, or do not import data, continue
          if(res.mapped_props.to.empty() || !res.has_data || !m_import_data) {
            continue;
          }

          // else we will try to import data

          // for import we currently don't have a way to get the initial config,
          // so set from = to, if not already so
          //   could support properties.calc.json: "ideal_" or "initial_" DoF
          if(res.mapped_props.to != res.mapped_props.from) {
            res.mapped_props.from = res.mapped_props.to;
          }
          std::string from = res.mapped_props.from;
          std::string to = res.mapped_props.to;

          // store info about data import, note if preexisting data before this batch
          auto data_res_it = data_results.find(from);
          if(data_res_it == data_results.end()) {
            data_res_it = data_results.insert(std::make_pair(from, ConfigIO::Result())).first;
            data_res_it->second.preexisting = has_existing_data_or_files(from);
          }
          ConfigIO::Result &data_res = data_res_it->second;

          // if preexisting data, do not import new data unless overwrite option set
          if(data_res.preexisting && !m_overwrite) {
            continue;
          }

          // if existing data (could be from this batch), check if score would
          //   be improved
          if(this->db_props().find_via_to(to) != db_props().end()) {
            double score = db_props().score(res.mapped_props);
            double best_score = db_props().best_score(to);
            if(score >= best_score) {
              continue;
            }
          }

          // note which structure is the latest import to configuration in case of conflicts
          data_res.last_insert = res.pos;

          // insert data (erasing any existing):
          auto pos = db_props().find_via_from(from);
          if(pos != db_props().end()) {
            db_props().erase(pos);
          }
          db_props().insert(res.mapped_props);

          // copy files:
          //   there might be existing files in cases of import conflicts
          rm_files(from, false);
          std::tie(data_res.copy_data, data_res.copy_more) =
            cp_files(res.pos, from, false, m_copy_additional_files);
        }

      }

      _import_report(results, data_results);

      db_supercell().commit();
      db_config().commit();
      db_props().commit();
    }

    template<typename _ConfigType>
    void ConfigDataGeneric::_import_report(
      std::vector<ConfigIO::Result> &results,
      const std::map<std::string, ConfigIO::ImportData> &data_results) {

      // map_fail: could not map
      // map_success: could map
      // import_data_fail: would import but couldn't (score < best_score && !data_results.count(from))
      // import_data_conflicts: conflicts with other in import batch && preexisting
      // - pos, config, score_method, chosen?, overwrite?, import data?, import additional files?, score, best_score, is_preexisting?


      // list of structures that could not be mapped
      std::vector<ConfigIO::Result> map_fail;

      // list of structures that could be mapped
      std::vector<ConfigIO::Result> map_success;

      // list of structures that would be imported except preexisting data prevents it
      std::vector<ConfigIO::Result> import_data_fail;

      std::string prefix = "import_";
      prefix += traits<ConfigType>::short_name;

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
      std::vector<ConfigIO::Result> conflict;

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

        fs::path p = m_report_dir / (prefix + "_map_fail");
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

        fs::path p = m_report_dir / (prefix + "_map_success");
        fs::ofstream sout(p);

        log() << "Successfully imported " << map_success.size() << " structures." << std::endl;
        log() << "  See: " << p << " for details" << std::endl << std::endl;

        sout << formatter(map_success.begin(), map_success.end());
      }

      if(import_data_fail.size()) {

        fs::path p = m_report_dir / (prefix + "_data_fail");
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
        fs::path p = m_report_dir / (prefix + "_conflict");
        fs::ofstream sout(p);

        log() << "WARNING: Imported data from structures that mapped to the same configuration." << std::endl
              << "  Data can only be imported from one of the conflicting structures." << std::endl
              << "  Based on the current conflict resolution method the 'best' result was automatically chosen, " << std::endl
              << "  but you may wish to inspect these results and manually select which structures to import." << std::endl;
        log() << "  See: " << p << " for details" << std::endl << std::endl;

        sout << formatter(conflict.begin(), conflict.end());
      }
    }

    // --- UpdateT ---

    /// \brief Re-parse calculations 'from' all selected configurations
    template<typename _ConfigType>
    void UpdateT<_ConfigType>::update(const DB::Selection<ConfigType> &selection, bool force) {

      // vector of Mapping results
      std::vector<Result> results;

      for(const auto &val : selection.data()) {

        // if not selected, skip
        if(!val.second) {
          continue;
        }
        std::string name = val.first;
        fs::path pos = calc_properties_path(primclex(), name);

        // reasons to update data or not:
        //
        // for each selected 'from' config:
        //
        //   if no existing files:
        //     continue
        //   if existing data or files:
        //     if !force && no change:
        //       continue
        //     erase existing props && read && map
        //     if could not map:
        //       'to' = ""; continue;
        //     else:
        //       update props
        //

        if(!_has_existing_files(name) || (!force && _no_change(name))) {
          continue;
        }

        // erase existing data (not files), unlinking relaxation mappings && resetting 'best' data
        db_props().erase_via_from(name);

        std::vector<Result> tvec;
        auto config_it = db_config().find(name);
        m_structure_mapper.map(pos, config_it, std::back_inserter(tvec));

        for(auto &res : tvec) {

          results.push_back(res);

          // if mapped && has data, insert
          if(!res.mapped_props.to.empty() && res.has_data) {
            // insert data:
            db_props().insert(res.mapped_props);
          }
        }

      }

      _update_report(result, selection);

      db_supercell().commit();
      db_config().commit();
      db_props().commit();
    }

    template<typename _ConfigType>
    void UpdateT<_ConfigType>::_update_report(std::vector<Result> &results, const DB::Selection<ConfigType> &selection) const {

      // report:
      //   update_map_fail: all 'from' config that were not successfully mapped
      //   update_map_success: all 'from' config that were successfully mapped
      //   update_map_conflict: all 'from' config that map to same 'to' config, sorted by 'to' config
      //   update_unstable: all 'from' config where 'from' != 'to', sorted by 'to' config
      //   update_unselected: update_unstable, where 'to' is not selected, sorted by 'to' config


      // list of structures that could (not) be mapped
      std::vector<Result> fail;
      std::vector<Result> success;
      std::vector<Result> conflict;
      std::vector<Result> unstable;
      std::vector<Result> unselected;

      std::string prefix = "update_";
      prefix += traits<ConfigType>::short_name;

      std::set<std::string> all_to;

      for(long i = 0; i < results.size(); ++i) {
        const auto &res = results[i];
        if(res.mapped_props.to.empty()) {
          fail.push_back(res);
        }
        else {
          all_to.insert(res.mapped_props.to);
          success.push_back(res);
        }
      }

      for(long i = 0; i < results.size(); ++i) {
        const auto &res = results[i];
        if(all_to.count(res.mapped_props.to)) {
          conflict.push_back(res);
          if(res.mapped_props.from != res.mapped_props.to) {
            unstable.push_back(res);
            if(!selection.is_selected(res.mapped_props.to)) {
              unselected.push_back(res);
            }
          }
        }
      }

      auto formatter = this->_update_formatter();

      if(fail.size()) {

        fs::path p = m_report_dir / (prefix + "_fail");
        fs::ostream sout(p);

        log() << "WARNING: Could not map " << fail.size() << " results." << std::endl;
        log() << "  See: " << p << " for details" << std::endl << std::endl;

        sout << formatter(fail.begin(), fail.end());
      }

      if(success.size()) {

        fs::path p = m_report_dir / (prefix + "_map_success");
        fs::ostream sout(p);

        log() << "Successfully mapped " << success.size() << " results." << std::endl;
        log() << "  See: " << p << " for details" << std::endl << std::endl;

        sout << formatter(success.begin(), success.end());
      }

      if(conflict.size()) {

        fs::path p = m_report_dir / (prefix + "_map_conflict");
        fs::ostream sout(p);

        log() << "WARNING: Found " << conflict.size() << " conflicting relaxation results." << std::endl;
        log() << "  See: " << p << " for details" << std::endl << std::endl;

        sout << formatter(conflict.begin(), conflict.end());
      }

      if(unstable.size()) {

        fs::path p = m_report_dir / (prefix + "_unstable");
        fs::ostream sout(p);

        log() << "WARNING: Found " << unstable.size() << " unstable relaxations." << std::endl;
        log() << "  See: " << p << " for details" << std::endl << std::endl;

        sout << formatter(unstable.begin(), unstable.end());
      }

      if(unselected.size()) {

        fs::path p = m_report_dir / (prefix + "_unselected");
        fs::ostream sout(p);

        log() << "WARNING: Found " << unselected.size() << " unstable relaxations to unselected configurations." << std::endl;
        log() << "  See: " << p << " for details" << std::endl << std::endl;

        sout << formatter(unselected.begin(), unselected.end());
      }


    }


    // --- RemoveT ---

    /// \brief Erase Configurations that have no data
    template<typename _ConfigType>
    void RemoveT<_ConfigType>::erase(const DB::Selection<ConfigType> &selection, bool dry_run) {
      std::vector<std::string> fail;
      for(const auto &val : selection.data()) {
        if(!_has_existing_data_or_files(val.first)) {
          db_config().erase(val.first);
        }
        else {
          log() << "skipping " << val.first << ": has existing data or files" << std::endl;
          m_fail.push_back(val.first);
        }
      }

      if(fail.size()) {
        _erase_report(fail);
        log() << "Skipped " << fail.size() << " " << ConfigType::name << std::endl;
        log() << "  See " << m_report_dir / "remove_fail" << std::endl;
        s
      }
      db_config().commit();
    }

    /// \brief Erase data and files (permanently), but not Configuration
    template<typename _ConfigType>
    void RemoveT<_ConfigType>::erase_data(const DB::Selection<ConfigType> &selection, bool dry_run) {
      // erase data
      for(const auto &val : selection.data()) {
        auto it = db_props().find_via_from(val.first);
        db_props().erase(it);
        rm_files(val.first, dry_run);
      }
      db_props().commit();
    }

    /// \brief Removes Configurations and data and files (permanently)
    ///
    /// - Data are always associated with one 'from' configuration, so the
    ///   selection here indicates 'from' configurations
    /// - The 'to' configurations are updated with the new best mapping properties
    template<typename _ConfigType>
    void RemoveT<_ConfigType>::erase_all(const DB::Selection<ConfigType> &selection, bool dry_run) {

      // erase data
      erase_data(selection);

      // erase configs
      erase(selection);
    }

    template<typename _ConfigType>
    void RemoveT<_ConfigType>::_erase_report(const std::vector<std::string> &fail) {
      std::string prefix {"remove_"};
      prefix += traits<ConfigType>::short_name;
      fs::ofstream file(m_report_dir / (prefix + "_fail"));
      for(const auto &val : fail) {
        file << val << std::endl;
      }
      file.close();
    }

    // --- Remove<Configuration> ---

    template<typename ConfigType>
    Remove<ConfigType>::Remove(
      const PrimClex &primclex,
      fs::path report_dir,
      Log &_file_log) :
      RemoveT(primclex, report_dir, _file_log) {}

    template<typename ConfigType>
    const std::string Remove<ConfigType>::desc = std::string("") +

                                                 "Remove enumerated configurations and calculation results: \n\n"

                                                 "  'casm remove --type " + traits<ConfigType>::short_name + "' options: \n\n"

                                                 "  - Configurations to be erased can be specified with the --names and \n"
                                                 "    --selection options.\n"
                                                 "  - Use without additional options to only remove enumerated configurations\n"
                                                 "    that do not have any associated files or data.\n"
                                                 "  - Use --data (-d) to remove data only, not enumerated configurations. \n"
                                                 "  - Use --force (-f) to remove data and enumerated configurations. \n"
                                                 "  - Use --dry-run (-n) to do a \"dry-run\". \n\n"

                                                 "  After removing a configuration it may be re-enumerated but will have a new\n"
                                                 "  index because indices will not be repeated.\n\n";

    template<typename ConfigType>
    int Remove<ConfigType>::run(
      const PrimClex &primclex,
      const Completer::RmOption &opt) {

      // -- read selection --
      DB::Selection<ConfigType> selection(primclex, opt.selection_path());
      for(const auto &name : opt.name_strs()) {
        if(primclex.db<ConfigType>().count(name)) {
          selection.data()[name] = true;
        }
        else {
          std::stringstream msg;
          msg << "Invalid Configuration name: " << name;
          throw CASM::runtime_error(msg.str(), ERR_INVALID_ARG);
        }
      }

      // get remove report_dir, check if exists, and create new report_dir.i if necessary
      fs::path report_dir = primclex.dir().root_dir() / "remove_report";
      report_dir = create_report_dir(report_dir);

      // -- erase --
      Remove<ConfigType> f(primclex, report_dir, primclex.log());

      if(opt.force()) {
        f.erase_all(selection, opt.dry_run());
      }
      else if(opt.data()) {
        f.erase_data(selection, opt.dry_run());
      }
      else {
        f.erase(selection, opt.dry_run());
      }
      return 0;
    }

  }
}

#endif
