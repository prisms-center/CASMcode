#ifndef CASM_DB_Import_impl
#define CASM_DB_Import_impl

#include <boost/filesystem/fstream.hpp>

#include "casm/app/DirectoryStructure.hh"
#include "casm/app/import.hh"
#include "casm/clex/PrimClex_impl.hh"
#include "casm/database/ConfigData_impl.hh"
#include "casm/database/Import.hh"
#include "casm/database/PropertiesDatabase.hh"
#include "casm/database/Selection_impl.hh"

namespace CASM {
namespace DB {

/// Construct pos_paths from input args
///
/// \return pair of OutputIterator, returncode
///
template <typename OutputIterator>
std::pair<OutputIterator, int> construct_pos_paths(
    const PrimClex &primclex, const Completer::ImportOption &import_opt,
    OutputIterator result) {
  const po::variables_map &vm = import_opt.vm();
  Index count = 0;
  bool missing_files(false);

  // count files, check if existing, print warning messages for missing files
  auto lambda = [&](const fs::path &p) {
    ++count;
    if (!fs::exists(p)) {
      if (!missing_files) {
        err_log() << "*** ERROR: Missing file(s):\n";
      }
      missing_files = true;
      err_log() << "           " << fs::absolute(p) << "\n";
    }
    *result++ = fs::absolute(p);
    return true;
  };

  // read all --pos paths
  std::for_each(import_opt.pos_vec().begin(), import_opt.pos_vec().end(),
                lambda);

  // read all the paths from a batch file
  if (vm.count("batch")) {
    if (!fs::exists(import_opt.batch_path())) {
      err_log() << "ERROR: Batch import file does not exist at "
                << import_opt.batch_path() << "\n";
      return std::make_pair(result, ERR_MISSING_INPUT_FILE);
    }

    fs::ifstream batchfile(import_opt.batch_path());
    fs::path tpath;
    while (batchfile >> tpath) {
      lambda(tpath);
      ++count;
      batchfile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
    if (tpath.string().size() != 0 && fs::exists(tpath)) {
      lambda(tpath);
    }
  }

  if (count == 0) {
    err_log() << "ERROR: No files specified for import.\n";
    if (vm.count("batch")) {
      err_log() << "       Check batch file for errors.\n";
    }
    return std::make_pair(result, ERR_INVALID_INPUT_FILE);
  }

  if (missing_files) {
    return std::make_pair(result, ERR_INVALID_INPUT_FILE);
  }

  return std::make_pair(result, 0);
}

// --- ImportT ---

/// \brief Constructor
template <typename _ConfigType>
ImportT<_ConfigType>::ImportT(const PrimClex &primclex,
                              const StructureMap<ConfigType> &mapper,
                              ImportSettings const &_set,
                              std::string report_dir)
    : ConfigData(primclex, TypeTag<ConfigType>()),
      m_structure_mapper(mapper),
      m_set(_set),
      m_report_dir(report_dir) {}

template <typename _ConfigType>
template <typename PathIterator>
void ImportT<_ConfigType>::import(PathIterator begin, PathIterator end) {
  // vector of Mapping results, may be >1 per input if primitive and
  // non-primitive
  std::vector<ConfigIO::Result> results;

  // map of data import results
  //   'configname' -> 'preexisting?
  std::map<std::string, bool> preexisting;
  std::map<std::string, bool> preexisting_files;

  auto const &project_settings = this->primclex().settings();
  auto calctype = project_settings.default_clex().calctype;
  auto required_properties =
      project_settings.required_properties(traits<ConfigType>::name, calctype);

  Log &log = CASM::log();
  auto it = begin;
  for (; it != end; ++it) {
    log << "Importing " << it->string() << std::endl;

    std::vector<ConfigIO::Result> tvec;

    // Outputs one or more mapping results from the structure located at specied
    // path
    //   See _import documentation for more.
    m_structure_mapper.map(it->string(), required_properties, nullptr,
                           std::back_inserter(tvec));

    // if successfully mapped:
    // - check for preexisting properties and files
    for (auto &res : tvec) {
      if (!res.properties.to.empty()) {
        // note if preexisting properties before this batch
        auto p_it = preexisting.find(res.properties.to);
        if (p_it == preexisting.end()) {
          p_it = preexisting
                     .emplace(res.properties.to,
                              has_existing_data(res.properties.to))
                     .first;
        }
        res.import_data.preexisting = p_it->second;
      }

      if (!res.properties.to.empty()) {
        // note if preexisting files before this batch
        auto p_it = preexisting_files.find(res.properties.to);
        if (p_it == preexisting_files.end()) {
          p_it = preexisting_files
                     .emplace(res.properties.to,
                              has_existing_files(res.properties.to))
                     .first;
        }
        res.import_data.preexisting_files = p_it->second;
      }
    }

    // import properties if:
    // - could map structure
    // - import_properties == true
    for (auto &res : tvec) {
      if (!res.properties.to.empty() && settings().import_properties) {
        // we will try to import data
        // insert properties
        // - first erase in case properties from structure already inserted
        // - assume no "force" necessary
        db_props().erase_via_origin(res.properties.origin);
        db_props().insert(res.properties);
      }
    }

    // add individual structure mapping results to batch results map
    for (auto &res : tvec) {
      results.push_back(res);
    }
  }

  // Check if result is the new best conflict scoring mapping
  double tol = this->primclex().settings().lin_alg_tol();
  for (auto &res : results) {
    res.import_data.is_new_best = false;
    if (res.properties.to.empty()) continue;

    fs::path p = calc_dir(res.properties.to);
    std::string calc_dir_origin = (p / "properties.calc.json").string();
    if (res.properties.origin == calc_dir_origin) continue;

    auto all_origins = db_props().all_origins(res.properties.to);
    if (all_origins.size() == 0) {
      res.import_data.is_new_best = true;
      continue;
    }

    auto first_origin = *all_origins.begin();
    if (first_origin != res.properties.origin) continue;
    if (all_origins.size() == 1) {
      res.import_data.is_new_best = true;
      continue;
    }

    auto db_it = db_props().find_via_origin(calc_dir_origin);
    if (db_it == db_props().end()) {
      res.import_data.is_new_best = true;
      continue;
    }

    if (!all_origins.count(calc_dir_origin)) {
      res.import_data.is_new_best = true;
      continue;
    } else {
      // check for approximate ties -- not a new best if approx equal
      double preexisting_score = db_props().score(calc_dir_origin);
      double score = db_props().score(first_origin);
      if (score + tol < preexisting_score) {
        res.import_data.is_new_best = true;
        continue;
      }
    }
  }

  // Copy files as needed/requested
  this->_copy_files(results);

  this->_import_report(results);

  db_supercell().commit();
  db_config<ConfigType>().commit();
  db_props().commit();
}

// *********************************************************************************

template <typename _ConfigType>
void ImportT<_ConfigType>::_copy_files(
    std::vector<ConfigIO::Result> &results) const {
  if (!settings().copy_structure_files) {
    return;
  }

  // Check for new best scoring data
  for (auto &res : results) {
    if (!res.import_data.is_new_best) continue;
    if (res.import_data.preexisting || res.import_data.preexisting_files) {
      if (!settings().overwrite) continue;
    }
    if (!res.has_files) continue;

    // at this point we are:
    // - going to copy structure files and maybe additional files into
    //   the training data directory for the "res.properties.to" configuration
    // - going to update the properties database to reflect the changes

    fs::path p = calc_dir(res.properties.to);
    std::string calc_dir_origin = (p / "properties.calc.json").string();

    // check existing properties in training_data
    auto db_it = db_props().find_via_origin(calc_dir_origin);
    if (db_it != db_props().end()) {
      db_props().erase(db_it);
    }

    // also, erase properties that were inserted from res origin, because
    // after copy we will re-insert with new origin
    db_it = db_props().find_via_origin(res.properties.origin);
    if (db_it != db_props().end()) {
      db_props().erase(db_it);
    }

    // copy files:
    //   there might be existing files in cases of import conflicts
    //   cp_files() will update res.properties.file_data
    bool dry_run = false;
    rm_files(res.properties.to, dry_run);
    this->cp_files(res, dry_run, settings().copy_additional_files);

    // update origin to new location and re-insert
    res.properties.origin = calc_dir_origin;
    auto insert_result = db_props().insert(res.properties);
    if (!insert_result.second) {
      std::stringstream msg;
      msg << "Unknown properties insert error after copying data into "
             "training_data directory for structure imported from "
          << res.pos_path;
      throw std::runtime_error(msg.str());
    }
  }
}

// *********************************************************************************

template <typename _ConfigType>
void ImportT<_ConfigType>::_import_report(
    std::vector<ConfigIO::Result> &results) const {
  // map_fail: could not map
  // map_success: could map
  // import_data_fail: would import but couldn't (score < best_score &&
  // !data_results.count(from)) import_data_conflicts: conflicts with other in
  // import batch && preexisting
  // - pos, config, score_method, chosen?, overwrite?, import data?, import
  // additional files?, score, best_score, is_preexisting?

  // list of structures that could not be mapped
  std::vector<ConfigIO::Result> map_fail;

  // list of structures that could be mapped
  std::vector<ConfigIO::Result> map_success;

  // list of structures that would be imported except preexisting data prevents
  // it
  std::vector<ConfigIO::Result> import_data_fail;

  std::map<std::string, int> all_to;

  for (auto const &res : results) {
    if (res.properties.to.empty()) {
      map_fail.push_back(res);
    } else {
      auto it = all_to.find(res.properties.to);
      if (it == all_to.end()) {
        it = all_to.insert(std::make_pair(res.properties.to, 0)).first;
      }
      ++(it->second);

      map_success.push_back(res);
      if (res.has_data && settings().import_properties &&
          db_props().find_via_to(res.properties.to) != db_props().end() &&
          db_props().score(res.properties) <
              db_props().best_score(res.properties.to)) {
        import_data_fail.push_back(res);
      }
    }
  }

  // list of conflicts (multiple config with same 'to')
  std::vector<ConfigIO::Result> conflict;
  for (auto const &res : results) {
    if (!res.properties.to.empty() && all_to[res.properties.to] > 1) {
      conflict.push_back(res);
    }
  }

  // output a 'batch' file with paths to structures that could not be mapped
  if (map_fail.size()) {
    fs::path p = fs::path(m_report_dir) / ("map_fail");
    if (settings().output_as_json) {
      p += ".json";
    }
    fs::ofstream sout(p);

    log() << "WARNING: Could not map " << map_fail.size() << " structures."
          << std::endl;
    log() << "  See detailed report: " << p << std::endl << std::endl;

    DataFormatterDictionary<ConfigIO::Result> dict;
    dict.insert(ConfigIO::initial_path(), ConfigIO::fail_msg());
    auto formatter = dict.parse({"initial_path", "fail_msg"});
    if (settings().output_as_json) {
      jsonParser json;
      json = formatter(map_fail.begin(), map_fail.end());
      sout << json;
    } else {
      sout << formatter(map_fail.begin(), map_fail.end());
    }
  }

  // - pos, config, score_method, import data?, import additional files?, score,
  // best_score, is_preexisting?
  auto formatter = _import_formatter();
  auto write_report = [&](std::vector<ConfigIO::Result> const &results,
                          std::ostream &sout) {
    if (settings().output_as_json) {
      jsonParser json;
      json = formatter(results.begin(), results.end());
      sout << json;
    } else {
      sout << formatter(results.begin(), results.end());
    }
  };

  if (map_success.size()) {
    fs::path p = fs::path(m_report_dir) / ("map_success");
    if (settings().output_as_json) {
      p += ".json";
    }
    fs::ofstream sout(p);

    log() << "Successfully mapped " << map_success.size() << " structures."
          << std::endl;
    log() << "  See detailed report: " << p << std::endl << std::endl;

    write_report(map_success, sout);
  }

  if (import_data_fail.size()) {
    fs::path p = fs::path(m_report_dir) / ("import_data_fail");
    if (settings().output_as_json) {
      p += ".json";
    }
    fs::ofstream sout(p);

    log() << "WARNING: Did not import data or copy files from "
          << import_data_fail.size()
          << " structures which have a mapping score"
             " better than the existing data."
          << std::endl;
    log() << "  You may wish to inspect these structures and allow overwriting "
             "or remove existing data manually."
          << std::endl;
    log() << "  See detailed report: " << p << std::endl << std::endl;

    write_report(import_data_fail, sout);
  }

  if (conflict.size()) {
    fs::path p = fs::path(m_report_dir) / ("import_conflict");
    if (settings().output_as_json) {
      p += ".json";
    }
    fs::ofstream sout(p);

    log()
        << "WARNING: Imported data from structures that mapped to the same "
           "configuration."
        << std::endl
        << "  Data can only be imported from one of the conflicting structures."
        << std::endl
        << "  Based on the current conflict resolution method the 'best' "
           "result was automatically chosen, "
        << std::endl
        << "  but you may wish to inspect these results and manually select "
           "which structures to import."
        << std::endl;
    log() << "  See detailed report: " << p << std::endl << std::endl;

    write_report(conflict, sout);
  }
}

}  // namespace DB
}  // namespace CASM

#endif
