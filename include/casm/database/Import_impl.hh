#include "casm/database/Import.hh"

namespace CASM {
  namespace DB {

    /// Construct pos_paths from input args
    ///
    /// \return pair of OutputIterator, returncode
    ///
    template<typename OutputIterator>
    std::pair<OutputIterator, int> Import::construct_pos_paths(
      const PrimClex &primclex,
      const Completer::ImportOption &import_opt,
      OutputIterator result) {

      po::variables_map &vm = import_opt.vm();
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

    template<typename _ConfigType>
    template<typename PathIterator>
    void ImportT<_ConfigType>::import(PathIterator begin, PathIterator end) {

      // vector of Mapping results
      std::vector<Result> results;

      // map of data import results
      //   'configname' -> 'preexisting?, copy_data?, copy_more?, last_import
      std::map<std::string, ImportData> data_results;

      auto it = begin;
      for(; it != end; ++it) {

        std::vector<Result> tvec;

        // Outputs one or more mapping results from the structure located at specied path
        //   See _import documentation for more.
        _import(*it, db_config().end(), std::back_inserter(tvec));

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
            data_res_it = data_res_it.insert({from, ImportResult()}).first;
            data_res_it->preexisting = _has_existing_data_or_files(from);
          }
          ImportResult &data_res = data_res_it->second;

          // if preexisting data, do not import new data unless overwrite option set
          if(data_res.preexisting && !m_overwrite) {
            continue;
          }

          // if existing data (could be from this batch), check if score would
          //   be improved
          if(db_props().find_via_to(to) != db_props().end()) {
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
          _rm_files(from);
          std::tie(data_res.copy_data, data_res.copy_more) = _cp_files(res.pos, from);
        }

      }

      _import_report(results, data_results);

      db_supercell().commit();
      db_config<ConfigType>().commit();
      db_props<ConfigType>().commit();
    }

    /// \brief Re-parse calculations 'from' all selected configurations
    template<typename _ConfigType>
    void ImportT<_ConfigType>::update(const DB::Selection<ConfigType> &selection, bool force) {

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
        _import(pos, config_it, std::back_inserter(tvec));

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

    /// \brief Erase Configurations that have no data
    template<typename _ConfigType>
    void ImportT<_ConfigType>::erase(const DB::Selection<ConfigType> &selection) {
      std::vector<std::string> fail;
      for(const auto &val : selection.data()) {
        if(!_has_existing_data_or_files(val.first)) {
          db_config<ConfigType>().erase(val.first);
        }
        else {
          m_fail.push_back(val.first);
        }
      }

      _erase_fail_report(fail);
      db_config().commit();
    }

    /// \brief Erase data and files (permanently), but not Configuration
    template<typename _ConfigType>
    void ImportT<_ConfigType>::erase_data(const DB::Selection<ConfigType> &selection) {
      // erase data
      for(const auto &val : selection.data()) {
        auto it = db_props<ConfigType>().find_via_from(val.first);
        db_props<ConfigType>().erase(it);
        _rm_files(val.first);
      }
      db_props().commit();
    }

    /// \brief Removes Configurations and data and files (permanently)
    ///
    /// - Data are always associated with one 'from' configuration, so the
    ///   selection here indicates 'from' configurations
    /// - The 'to' configurations are updated with the new best mapping properties
    template<typename _ConfigType>
    void ImportT<_ConfigType>::erase_all(const DB::Selection<ConfigType> &selection) {

      // erase data
      erase_data(selection);

      // erase configs
      erase(selection);
    }

    /// \brief Path to default calctype training_data directory for config
    template<typename _ConfigType>
    fs::path ImportT<_ConfigType>::calc_dir(const std::string configname) const {
      return primclex().dir().calc_dir<ConfigType>(configname);
    }

    // --- Specializations ---

    template<typename _ConfigType>
    void ImportT<_ConfigType>::_update_report(std::vector<Result> &results, const DB::Selection<ConfigType> &selection) const {

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

        fs::path p = m_report_dir / "update_map_fail";
        fs::ostream sout(p);

        log() << "WARNING: Could not map " << fail.size() << " results." << std::endl;
        log() << "  See: " << p << " for details" << std::endl << std::endl;

        sout << formatter(fail.begin(), fail.end());
      }

      if(success.size()) {

        fs::path p = m_report_dir / "update_map_success";
        fs::ostream sout(p);

        log() << "Successfully mapped " << success.size() << " results." << std::endl;
        log() << "  See: " << p << " for details" << std::endl << std::endl;

        sout << formatter(success.begin(), success.end());
      }

      if(conflict.size()) {

        fs::path p = m_report_dir / "update_map_conflict";
        fs::ostream sout(p);

        log() << "WARNING: Found " << conflict.size() << " conflicting relaxation results." << std::endl;
        log() << "  See: " << p << " for details" << std::endl << std::endl;

        sout << formatter(conflict.begin(), conflict.end());
      }

      if(unstable.size()) {

        fs::path p = m_report_dir / "update_unstable";
        fs::ostream sout(p);

        log() << "WARNING: Found " << unstable.size() << " unstable relaxations." << std::endl;
        log() << "  See: " << p << " for details" << std::endl << std::endl;

        sout << formatter(unstable.begin(), unstable.end());
      }

      if(unselected.size()) {

        fs::path p = m_report_dir / "update_unselected";
        fs::ostream sout(p);

        log() << "WARNING: Found " << unselected.size() << " unstable relaxations to unselected configurations." << std::endl;
        log() << "  See: " << p << " for details" << std::endl << std::endl;

        sout << formatter(unselected.begin(), unselected.end());
      }


    }


  }
}
