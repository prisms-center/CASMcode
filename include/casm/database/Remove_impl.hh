#ifndef CASM_DB_Remove_impl
#define CASM_DB_Remove_impl

#include <boost/filesystem/fstream.hpp>

#include "casm/app/DirectoryStructure.hh"
#include "casm/app/import.hh"
#include "casm/app/rm.hh"
#include "casm/clex/Calculable.hh"
#include "casm/clex/PrimClex_impl.hh"
#include "casm/database/ConfigData_impl.hh"
#include "casm/database/PropertiesDatabase.hh"
#include "casm/database/Remove.hh"
#include "casm/database/Selection_impl.hh"

namespace CASM {
namespace DB {

// --- RemoveT ---

template <typename _ConfigType>
RemoveT<_ConfigType>::RemoveT(const PrimClex &primclex, std::string report_dir)
    : ConfigData(primclex, TypeTag<ConfigType>()), m_report_dir(report_dir) {}

/// \brief Erase Configurations that have no data
template <typename _ConfigType>
void RemoveT<_ConfigType>::erase(const DB::Selection<ConfigType> &selection,
                                 bool dry_run) {
  std::vector<std::string> fail;
  for (const auto &val : selection.data()) {
    if (!has_existing_data_or_files(val.first)) {
      db_config<ConfigType>().erase(val.first);
    } else {
      log() << "skipping " << val.first << ": has existing data or files"
            << std::endl;
      fail.push_back(val.first);
    }
  }

  if (fail.size()) {
    _erase_report(fail);
    log() << "Skipped " << fail.size() << " " << traits<ConfigType>::name
          << std::endl;
    log() << "  See " << fs::path(m_report_dir) / "remove_fail" << std::endl;
  }
  db_config<ConfigType>().commit();
}

/// \brief Erase properties via structure origin file path
template <typename _ConfigType>
void RemoveT<_ConfigType>::erase_properties_via_origin(
    std::vector<std::string> const &origins, bool dry_run) {
  std::vector<std::string> fail;
  for (auto origin : origins) {
    if (db_props().erase_via_origin(origin)) {
      log() << "Removed properties from: " << origin << std::endl;
    } else {
      log() << "Could not remove: " << origin << std::endl;
      fail.push_back(origin);
    }
  }

  if (fail.size()) {
    _erase_report(fail);
    log() << "Could not remove properties from " << fail.size() << " structures"
          << std::endl;
    log() << "  See " << fs::path(m_report_dir) / "remove_fail" << std::endl;
  }
  if (!dry_run) {
    db_props().commit();
  }
}

/// \brief Erase data and files (permanently), but not Configuration
template <typename _ConfigType>
void RemoveT<_ConfigType>::erase_data(
    const DB::Selection<ConfigType> &selection, bool dry_run) {
  // erase data
  for (const auto &val : selection.data()) {
    for (auto it = db_props().find_via_to(val.first); it != db_props().end();
         it = db_props().find_via_to(val.first)) {
      db_props().erase(it);
    }
    auto it = db_props().find_via_origin(
        CASM::calc_properties_path(this->primclex(), val.first));
    if (it != db_props().end()) {
      db_props().erase(it);
    }
    rm_files(val.first, dry_run);
  }
  if (!dry_run) {
    db_props().commit();
  }
}

/// \brief Removes Configurations and data and files (permanently)
///
/// - Data are always associated with one 'to' configuration, so the
///   selection here indicates 'to' configurations
/// - The 'from' configurations are updated with the new best mapping properties
template <typename _ConfigType>
void RemoveT<_ConfigType>::erase_all(const DB::Selection<ConfigType> &selection,
                                     bool dry_run) {
  // erase data
  erase_data(selection, dry_run);

  // erase configs
  erase(selection, dry_run);
}

template <typename _ConfigType>
void RemoveT<_ConfigType>::_erase_report(const std::vector<std::string> &fail) {
  fs::ofstream file(fs::path(m_report_dir) / ("remove_fail"));
  for (const auto &val : fail) {
    file << val << std::endl;
  }
  file.close();
}

// --- Remove<Configuration> ---

template <typename ConfigType>
Remove<ConfigType>::Remove(const PrimClex &primclex, std::string report_dir)
    : RemoveT<ConfigType>(primclex, report_dir) {}

template <typename ConfigType>
std::string Remove<ConfigType>::desc() {
  std::string res =
      "Remove enumerated configurations and calculation results: \n\n"

      "  'casm remove --type " +
      traits<ConfigType>::short_name +
      "' options: \n\n"

      "  - Configurations to be erased can be specified with the --names and \n"
      "    --selection options.\n"
      "  - Use without additional options to only remove enumerated \n"
      "    configurations that do not have any associated files or \n"
      "    properties.\n"
      "  - Use --data (-d) to remove data only (includes files and \n"
      "    properties), not enumerated configurations. Removes any properties\n"
      "    mapped to the specified configurations, properties mapped from the\n"
      "    `properties.calc.json` file in the configurations' training_data \n"
      "    directory, and all files from the configurations' training_data \n"
      "    directory. \n"
      "  - Use --force (-f) to remove data (includes files and \n"
      "    properties) and enumerated configurations. \n"
      "  - Use --dry-run (-n) to do a \"dry-run\". \n\n"

      "  - Properties imported by mapping from a particular structure file \n"
      "    may be removed from the properties database using \n"
      "    --structure-properties.\n\n"

      "  After removing a configuration it may be re-enumerated but will have "
      "  a new index because indices will not be repeated.\n\n";

  return res;
}

template <typename ConfigType>
int Remove<ConfigType>::run(const PrimClex &primclex,
                            const Completer::RmOption &opt) {
  // get remove report_dir, check if exists, and create new report_dir.i if
  // necessary
  std::string report_dir =
      (fs::path(primclex.dir().reports_dir()) / "remove_report").string();
  report_dir = create_report_dir(report_dir);

  // -- erase --
  Remove<ConfigType> f(primclex, report_dir);

  if (opt.vm().count("structure-properties")) {
    f.erase_properties_via_origin(opt.structure_properties(), opt.dry_run());
  }

  if (opt.vm().count("selection") || opt.vm().count("names")) {
    // -- read selection --
    fs::path selection_path;
    if (opt.vm().count("selection")) {
      selection_path = opt.selection_path();
    } else {
      selection_path = "NONE";
    }
    DB::Selection<ConfigType> selection(primclex, selection_path);

    // Add command-line options
    for (const auto &name : opt.name_strs()) {
      if (primclex.db<ConfigType>().count(name)) {
        selection.data()[name] = true;
      } else {
        std::stringstream msg;
        msg << "Invalid Configuration name: " << name;
        throw CASM::runtime_error(msg.str(), ERR_INVALID_ARG);
      }
    }

    if (opt.force()) {
      f.erase_all(selection, opt.dry_run());
    } else if (opt.data()) {
      f.erase_data(selection, opt.dry_run());
    } else {
      f.erase(selection, opt.dry_run());
    }
  }
  return 0;
}

}  // namespace DB
}  // namespace CASM

#endif
