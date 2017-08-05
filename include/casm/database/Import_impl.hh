#ifndef CASM_DB_Import_impl
#define CASM_DB_Import_impl

#include <boost/filesystem/fstream.hpp>
#include "casm/database/Import.hh"
#include "casm/database/Selection_impl.hh"
#include "casm/app/DirectoryStructure.hh"
#include "casm/app/import.hh"
#include "casm/clex/PrimClex_impl.hh"

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
          return std::make_pair(result, ERR_MISSING_INPUT_FILE);
        }

        fs::ifstream batchfile(import_opt.batch_path());
        fs::path tpath;
        while(batchfile >> tpath) {
          lambda(tpath);
          ++count;
          batchfile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        }
        if(tpath.string().size() != 0 && fs::exists(tpath)) {
          lambda(tpath);
        }
      }

      if(count == 0) {
        primclex.err_log() <<   "ERROR: No files specified for import.\n";
        if(vm.count("batch")) {
          primclex.err_log() << "       Check batch file for errors.\n";
        }
        return std::make_pair(result, ERR_INVALID_INPUT_FILE);
      }

      if(missing_files) {
        return std::make_pair(result, ERR_INVALID_INPUT_FILE);
      }

      return std::make_pair(result, 0);
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
        m_structure_mapper.map(*it, db_config().end(), std::back_inserter(tvec));

        for(auto &res : tvec) {

          results.push_back(res);

          // reasons to import data or not:
          //
          // could not map || no data || !m_import_data:
          //   do not import
          // !has_existing_data_or_files:
          //   import
          // has_existing_data_or_files:
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
            data_res_it = data_results.insert({from, ConfigIO::ImportData()}).first;
            data_res_it->second.preexisting = has_existing_data_or_files(from);
          }
          ConfigIO::ImportData &data_res = data_res_it->second;

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
    void ImportT<_ConfigType>::_import_report(
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

        DataFormatterDictionary<ConfigIO::Result> dict;
        dict.insert(ConfigIO::pos(), ConfigIO::fail_msg());
        auto formatter = dict.parse({"pos", "fail_msg"});
        sout << formatter(map_fail.begin(), map_fail.end());
      }

      // - pos, config, score_method, import data?, import additional files?, score, best_score, is_preexisting?
      auto formatter = _import_formatter(data_results);

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

  }
}

#endif
