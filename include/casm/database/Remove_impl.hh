#ifndef CASM_DB_Remove_impl
#define CASM_DB_Remove_impl

#include <boost/filesystem/fstream.hpp>
#include "casm/database/Remove.hh"
#include "casm/database/Selection_impl.hh"
#include "casm/database/PropertiesDatabase.hh"
#include "casm/app/DirectoryStructure.hh"
#include "casm/app/import.hh"
#include "casm/clex/PrimClex_impl.hh"
#include "casm/app/rm.hh"

namespace CASM {
  namespace DB {

    // --- RemoveT ---

    /// \brief Erase Configurations that have no data
    template<typename _ConfigType>
    void RemoveT<_ConfigType>::erase(const DB::Selection<ConfigType> &selection, bool dry_run) {
      std::vector<std::string> fail;
      for(const auto &val : selection.data()) {
        if(!has_existing_data_or_files(val.first)) {
          db_config().erase(val.first);
        }
        else {
          log() << "skipping " << val.first << ": has existing data or files" << std::endl;
          fail.push_back(val.first);
        }
      }

      if(fail.size()) {
        _erase_report(fail);
        log() << "Skipped " << fail.size() << " " << traits<ConfigType>::name << std::endl;
        log() << "  See " << m_report_dir / "remove_fail" << std::endl;
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
      erase_data(selection, dry_run);

      // erase configs
      erase(selection, dry_run);
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
      RemoveT<ConfigType>(primclex, report_dir, _file_log) {}

    template<typename ConfigType>
    std::string Remove<ConfigType>::desc() {
      std::string res =
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

      return res;
    }

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
