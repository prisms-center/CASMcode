#include "casm/kinetics/RemovePrimPeriodicDiffTransOrbit.hh"
#include "casm/app/rm.hh"
#include "casm/app/DirectoryStructure.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/database/Selection.hh"
#include "casm/database/Remove.hh"
#include "casm/database/DatabaseDefs.hh"
#include "casm/database/DatabaseTypes.hh"

namespace CASM {
  namespace DB {

    Remove<PrimPeriodicDiffTransOrbit>::Remove(const PrimClex &_primclex, fs::path report_dir, Log &_file_log) :
      m_primclex(_primclex), m_report_dir(report_dir), m_file_log(_file_log) {}

    const std::string Remove<PrimPeriodicDiffTransOrbit>::desc =

      "Remove a difftrans, including all enumerated difftransconfigurations and calculation results: \n\n"

      "  'casm remove --type diff_trans' options: \n\n"

      "  - Diffusion transformations to be erased can be specified with the --names and \n"
      "    --selection options.\n"
      "  - Use without additional options to remove all enumerated difftransconfigurations\n"
      "    that do not have any associated files or data for each specified diffusion transformation.\n"
      "    If no data or files, erase diffusion transformation.\n"
      "  - Use --data (-d) to remove all difftransconfiguration data, but not enumerated \n"
      "    configurations, for each specified diffusion transformation. \n"
      "  - Use --force (-f) to remove specified diffusion transformations including all data and \n"
      "    enumerated configurations. \n"
      "  - Use --dry-run (-n) to do a \"dry-run\". \n\n";
    /*
    /// Helper struct base class
    struct EraseDifftransConfigsBase {

      EraseDiffTransConfigsBase(Remove<PrimPeriodicDiffTransOrbit> &_remover, std::string _orbitname, bool _dry_run) :
        primclex(_remover.primclex()), remover(_remover), orbitname(_orbitname),
        dry_run(_dry_run), remaining(0) {}

      template<typename T>
      DB::Selection<T> make_selection() const {
        DB::Selection<T> selection(remover.primclex(), "NONE");
        auto it = primclex.db<T>().scel_range(scelname).begin();
        auto end = primclex.db<T>().scel_range(scelname).end();
        for(; it != end; ++it) {
          selection.data()[it.name()] = true;
        }
        return selection;
      }

      template<typename T>
      void count_remaining() {
        if(!dry_run) {
          remaining += boost::distance(remover.primclex().db<T>().scel_range(scelname));
        }
        else {
          ConfigData<T> data(remover.primclex(), null_log());
          auto it = primclex.db<T>().scel_range(scelname).begin();
          auto end = primclex.db<T>().scel_range(scelname).end();
          for(; it != end; ++it) {
            if(data.has_existing_data_or_files(it.name())) {
              remaining++;
            }
          }
        }
      }

      const PrimClex &primclex;
      const Remove<Supercell> &remover;
      std::string orbitname;
      bool dry_run;
      Index remaining;

    };
    */
    /*
    /// Helper struct that calls Remove<ConfigType>::erase for all configs in a supercell
    struct EraseScelConfigs : public EraseScelConfigsBase {
      EraseScelConfigs(Remove<Supercell> &_remover, std::string _scelname, bool _dry_run) :
        EraseScelConfigsBase(_remover, _scelname, _dry_run) {}

      template<typename T>
      void eval() {
        Remove<T> f(primclex, remover.report_dir(), remover.file_log());
        DB::Selection<T> selection = make_selection<T>();
        f.erase(selection, dry_run);
        count_remaining<T>();
      }
    };
    */
    /*
    /// Helper struct that calls Remove<ConfigType>::erase_data for all configs in a supercell
    struct EraseDataScelConfigs : public EraseScelConfigsBase {
      EraseDataScelConfigs(Remove<Supercell> &_remover, std::string _scelname, bool _dry_run) :
        EraseScelConfigsBase(_remover, _scelname, _dry_run) {}

      template<typename T>
      void eval() {
        Remove<T> f(primclex, remover.report_dir(), remover.file_log());
        DB::Selection<T> selection = make_selection<T>();
        f.erase_data(selection, dry_run);
        count_remaining<T>();
      }
    };
    */
    /*
    /// Helper struct that calls Remove<ConfigType>::erase_all for all configs in a supercell
    struct EraseAllScelConfigs : public EraseScelConfigsBase {

      EraseAllScelConfigs(Remove<Supercell> &_remover, std::string _scelname, bool _dry_run) :
        EraseScelConfigsBase(_remover, _scelname, _dry_run) {}

      template<typename T>
      void eval() {
        Remove<T> f(primclex, remover.report_dir(), remover.file_log());
        DB::Selection<T> selection = make_selection<T>();
        f.erase_all(selection, dry_run);
        count_remaining<T>();
      }

    };
    */

    /*
    /// \brief Erase all enumerated configs that have no data; if no data, erase supercell
    void Remove<Supercell>::erase(
      const DB::Selection<Supercell> &selection,
      bool dry_run) {

      bool did_erase = false;

      // call Remove<ConfigType>::erase for each ConfigType
      auto it = selection.selected().begin();
      auto end = selection.selected().end();
      for(; it != end; ++it) {
        EraseScelConfigs f(*this, it.name(), dry_run);
        for_each_config_type(f);

        // if no existing data or files, erase Supercell ...
        if(!f.remaining) {
          primclex().log() << "will erase " << it.name() << "\n";
          if(!dry_run) {
            primclex().db<Supercell>().erase(it.name());
            did_erase = true;
          }
        }
        else {
          primclex().log() << "skipping " << it.name() << ": has "
                           << f.remaining << " configurations remaining.\n";
        }
      }

      if(did_erase) {
        primclex().db<Supercell>().commit();
      }
    }
    */

    /*
    /// \brief Erase all data ; if no data, erase supercell
    void Remove<Supercell>::erase_data(
      const DB::Selection<Supercell> &selection,
      bool dry_run) {

      // call Remove<ConfigType>::erase for each ConfigType
      auto it = selection.selected().begin();
      auto end = selection.selected().end();
      for(; it != end; ++it) {
        EraseDataScelConfigs f(*this, it.name(), dry_run);
        for_each_config_type(f);
      }
    }
    */

    /*
    /// \brief Erase all enumerated configs that have no data; if no data, erase supercell
    void Remove<Supercell>::erase_all(
      const DB::Selection<Supercell> &selection,
      bool dry_run) {

      // call Remove<ConfigType>::erase for each ConfigType
      auto it = selection.selected().begin();
      auto end = selection.selected().end();
      for(; it != end; ++it) {
        EraseAllScelConfigs f(*this, it.name(), dry_run);
        for_each_config_type(f);

        // if no existing data or files, erase Supercell ...
        if(!f.remaining) {
          primclex().log() << "will erase " << it.name() << "\n";
          if(dry_run) {
            primclex().db<Supercell>().erase(it.name());
          }
        }
        else {
          primclex().log() << "unknown error: erase_all called, but " << it.name()
                           << " still has " << f.remaining << " configurations.\n";
          primclex().log() << "stopping..." << std::endl;
          return;
        }
      }

      primclex().db<Supercell>().commit();
    }
    */

    int Remove<PrimPeriodicDiffTransOrbit>::run(
      const PrimClex &primclex,
      const Completer::RmOption &opt) {

      // -- read selection --
      DB::Selection<PrimPeriodicDiffTransOrbit> selection(primclex, opt.selection_path());
      for(const auto &name : opt.name_strs()) {
        if(primclex.db<PrimPeriodicDiffTransOrbit>().count(name)) {
          selection.data()[name] = true;
        }
        else {
          std::stringstream msg;
          msg << "Invalid Diffusion transformation name: " << name;
          throw CASM::runtime_error(msg.str(), ERR_INVALID_ARG);
        }
      }

      // get remove report_dir, check if exists, and create new report_dir.i if necessary
      fs::path report_dir = primclex.dir().root_dir() / "remove_report";
      report_dir = create_report_dir(report_dir);

      // -- erase --
      Remove<PrimPeriodicDiffTransOrbit> f(primclex, report_dir, primclex.log());
      /*
      if(opt.force()) {
        f.erase_all(selection, opt.dry_run());
      }
      else if(opt.data()) {
        f.erase_data(selection, opt.dry_run());
      }
      else {
        f.erase(selection, opt.dry_run());
      }*/
      return 0;
    }

    const PrimClex &Remove<PrimPeriodicDiffTransOrbit>::primclex() const {
      return m_primclex;
    }

    fs::path Remove<PrimPeriodicDiffTransOrbit>::report_dir() const {
      return m_report_dir;
    }

    Log &Remove<PrimPeriodicDiffTransOrbit>::file_log() const {
      return m_file_log;
    }

  }
}
