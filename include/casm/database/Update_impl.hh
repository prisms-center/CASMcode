#ifndef CASM_DB_Update_impl
#define CASM_DB_Update_impl

#include <boost/filesystem/fstream.hpp>

#include "casm/clex/PrimClex_impl.hh"
#include "casm/database/ConfigData_impl.hh"
#include "casm/database/Selection_impl.hh"
#include "casm/database/Update.hh"

namespace CASM {
namespace DB {

// --- UpdateT ---

/// \brief Constructor
template <typename _ConfigType>
UpdateT<_ConfigType>::UpdateT(const PrimClex &primclex,
                              const StructureMap<ConfigType> &mapper,
                              std::string report_dir)
    : ConfigData(primclex, TypeTag<ConfigType>()),
      m_structure_mapper(mapper),
      m_report_dir(report_dir) {}

/// \brief Re-parse calculations 'from' all selected configurations
template <typename _ConfigType>
void UpdateT<_ConfigType>::update(const DB::Selection<ConfigType> &selection,
                                  bool force) {
  Log &log = CASM::log();
  auto const &project_settings = this->primclex().settings();
  auto calctype = project_settings.default_clex().calctype;
  auto required_properties =
      project_settings.required_properties(traits<ConfigType>::name, calctype);

  // vector of Mapping results
  std::vector<ConfigIO::Result> results;
  for (const auto &val : selection.data()) {
    // if not selected, skip
    if (!val.second) {
      continue;
    }
    std::string name = val.first;
    std::string pos = CASM::calc_properties_path(this->primclex(), name);

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

    if ((!force && this->no_change(name))) {
      continue;
    }

    // erase existing data (not files), unlinking relaxation mappings &&
    // resetting 'best' data
    db_props().erase_via_origin(pos);

    if (!fs::exists(pos)) continue;

    log << "Updating data records for " << name << std::endl;
    std::vector<ConfigIO::Result> tvec;
    auto config_it = db_config<ConfigType>().find(name);
    if (config_it == db_config<ConfigType>().end())
      m_structure_mapper.map(resolve_struc_path(pos, primclex()),
                             required_properties, nullptr,
                             std::back_inserter(tvec));
    else
      m_structure_mapper.map(resolve_struc_path(pos, primclex()),
                             required_properties,
                             notstd::make_unique<ConfigType>(*config_it),
                             std::back_inserter(tvec));
    for (auto &res : tvec) {
      results.push_back(res);
      // if mapped && has data, insert
      if (!res.properties.to.empty() && res.has_data) {
        // insert data:
        db_props().insert(res.properties);
      }
    }
  }
  _update_report(results, selection);

  db_supercell().commit();
  db_config<ConfigType>().commit();
  db_props().commit();
}

template <typename _ConfigType>
void UpdateT<_ConfigType>::_update_report(
    std::vector<ConfigIO::Result> &results,
    const DB::Selection<ConfigType> &selection) const {
  // report:
  //   update_map_fail: all 'from' config that were not successfully mapped
  //   update_map_success: all 'from' config that were successfully mapped
  //   update_map_conflict: all 'from' config that map to same 'to' config,
  //   sorted by 'to' config update_unstable: all 'from' config where 'from' !=
  //   'to', sorted by 'to' config update_unselected: update_unstable, where
  //   'to' is not selected, sorted by 'to' config

  // list of structures that could (not) be mapped
  std::vector<ConfigIO::Result> fail;
  std::vector<ConfigIO::Result> success;
  std::vector<ConfigIO::Result> conflict;
  std::vector<ConfigIO::Result> unstable;
  std::vector<ConfigIO::Result> unselected;
  std::vector<ConfigIO::Result> new_config;

  std::string prefix = "update_";
  prefix += traits<ConfigType>::short_name;

  std::map<std::string, int> all_to;

  for (auto const &res : results) {
    if (res.properties.to.empty()) {
      fail.push_back(res);
    } else {
      auto it = all_to.find(res.properties.to);
      if (it == all_to.end()) {
        it = all_to.insert(std::make_pair(res.properties.to, 0)).first;
      }
      ++(it->second);
      if (res.properties.origin !=
          CASM::calc_properties_path(this->primclex(), res.properties.to)) {
        unstable.push_back(res);
        if (res.is_new_config) {
          new_config.push_back(res);
        }
        if (!selection.is_selected(res.properties.to)) {
          unselected.push_back(res);
        }
      }
      success.push_back(res);
    }
  }

  for (auto const &res : results) {
    if (all_to[res.properties.to] > 1) {
      conflict.push_back(res);
    }
  }

  auto formatter = _update_formatter();

  if (fail.size()) {
    fs::path p = fs::path(m_report_dir) / (prefix + "_fail");
    fs::ofstream sout(p);

    log() << "WARNING: Could not map " << fail.size() << " results."
          << std::endl;
    log() << "  See detailed report: " << p << std::endl << std::endl;

    sout << formatter(fail.begin(), fail.end());
  }

  if (success.size()) {
    fs::path p = fs::path(m_report_dir) / (prefix + "_success");
    fs::ofstream sout(p);

    log() << "Successfully mapped " << success.size() << " results."
          << std::endl;
    log() << "  See detailed report: " << p << std::endl << std::endl;

    sout << formatter(success.begin(), success.end());
  }

  if (conflict.size()) {
    fs::path p = fs::path(m_report_dir) / (prefix + "_conflict");
    fs::ofstream sout(p);

    log() << "WARNING: Found " << conflict.size()
          << " conflicting relaxation results." << std::endl;
    log() << "  See detailed report: " << p << std::endl << std::endl;

    sout << formatter(conflict.begin(), conflict.end());
  }

  if (unstable.size()) {
    fs::path p = fs::path(m_report_dir) / (prefix + "_unstable");
    fs::ofstream sout(p);

    log() << "WARNING: Found " << unstable.size() << " unstable relaxations."
          << std::endl;
    log() << "  See detailed report: " << p << std::endl << std::endl;

    sout << formatter(unstable.begin(), unstable.end());
  }

  if (unselected.size()) {
    fs::path p = fs::path(m_report_dir) / (prefix + "_unselected");
    fs::ofstream sout(p);

    log() << "WARNING: Found " << unselected.size()
          << " unstable relaxations to unselected configurations." << std::endl;
    log() << "  See detailed report: " << p << std::endl << std::endl;

    sout << formatter(unselected.begin(), unselected.end());
  }

  if (new_config.size()) {
    fs::path p = fs::path(m_report_dir) / (prefix + "_new");
    fs::ofstream sout(p);

    log() << "WARNING: Found " << new_config.size()
          << " unstable relaxations to new configurations." << std::endl
          << "  These configurations have been added to the project."
          << std::endl;
    log() << "  See detailed report: " << p << std::endl << std::endl;

    sout << formatter(unselected.begin(), unselected.end());
  }
}

}  // namespace DB
}  // namespace CASM

#endif
