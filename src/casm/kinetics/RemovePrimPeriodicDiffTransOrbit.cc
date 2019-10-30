#include "casm/kinetics/RemovePrimPeriodicDiffTransOrbit.hh"
#include "casm/app/rm.hh"
#include "casm/app/DirectoryStructure.hh"
#include "casm/database/Selection_impl.hh"
#include "casm/database/Remove_impl.hh"
#include "casm/database/DatabaseTypes_impl.hh"

namespace CASM {
  namespace DB {

    Remove<PrimPeriodicDiffTransOrbit>::Remove(const PrimClex &_primclex, fs::path report_dir, Log &_file_log) :
      m_primclex(_primclex), m_report_dir(report_dir), m_file_log(_file_log) {}

    std::string Remove<PrimPeriodicDiffTransOrbit>::desc() {

      std::string res =
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
      return res;
    }

    /// Helper struct base class
    struct ErasePrimPeriodicDiffTransOrbitConfigsBase {

      ErasePrimPeriodicDiffTransOrbitConfigsBase(Remove<Kinetics::PrimPeriodicDiffTransOrbit> &_remover, std::string _orbitname, bool _dry_run) :
        primclex(_remover.primclex()), remover(_remover), orbitname(_orbitname),
        dry_run(_dry_run), remaining(0) {}

      DB::Selection<Kinetics::DiffTransConfiguration> make_selection() const {
        DB::Selection<Kinetics::DiffTransConfiguration> selection(remover.primclex(), "NONE");
        auto it = primclex.db<Kinetics::DiffTransConfiguration>().orbit_range(orbitname).begin();
        auto end = primclex.db<Kinetics::DiffTransConfiguration>().orbit_range(orbitname).end();
        for(; it != end; ++it) {
          selection.data()[it.name()] = true;
        }
        return selection;
      }

      void count_remaining() {
        if(!dry_run) {
          remaining += boost::distance(remover.primclex().db<Kinetics::DiffTransConfiguration>().orbit_range(orbitname));
        }
        else {
          ConfigData data(remover.primclex(), null_log(), TypeTag<Kinetics::DiffTransConfiguration>());
          auto it = primclex.db<Kinetics::DiffTransConfiguration>().orbit_range(orbitname).begin();
          auto end = primclex.db<Kinetics::DiffTransConfiguration>().orbit_range(orbitname).end();
          for(; it != end; ++it) {
            if(data.has_existing_data_or_files(it.name())) {
              remaining++;
            }
          }
        }
      }

      const PrimClex &primclex;
      const Remove<Kinetics::PrimPeriodicDiffTransOrbit> &remover;
      std::string orbitname;
      bool dry_run;
      Index remaining;

    };


    /// Helper struct that calls Remove<ConfigType>::erase for all diff_trans_configs in a orbit
    struct ErasePrimPeriodicDiffTransOrbitConfigs : public ErasePrimPeriodicDiffTransOrbitConfigsBase {
      ErasePrimPeriodicDiffTransOrbitConfigs(Remove<Kinetics::PrimPeriodicDiffTransOrbit> &_remover, std::string _orbitname, bool _dry_run) :
        ErasePrimPeriodicDiffTransOrbitConfigsBase(_remover, _orbitname, _dry_run) {}

      void eval() {
        Remove<Kinetics::DiffTransConfiguration> f(primclex, remover.report_dir(), remover.file_log());
        DB::Selection<Kinetics::DiffTransConfiguration> selection = make_selection();
        f.erase(selection, dry_run);
        count_remaining();
      }
    };


    /// Helper struct that calls Remove<ConfigType>::erase_data for all configs in a supercell
    struct EraseDataPrimPeriodicDiffTransOrbitConfigs : public ErasePrimPeriodicDiffTransOrbitConfigsBase {
      EraseDataPrimPeriodicDiffTransOrbitConfigs(Remove<Kinetics::PrimPeriodicDiffTransOrbit> &_remover, std::string _orbitname, bool _dry_run) :
        ErasePrimPeriodicDiffTransOrbitConfigsBase(_remover, _orbitname, _dry_run) {}

      void eval() {
        Remove<Kinetics::DiffTransConfiguration> f(primclex, remover.report_dir(), remover.file_log());
        DB::Selection<Kinetics::DiffTransConfiguration> selection = make_selection();
        f.erase_data(selection, dry_run);
        count_remaining();
      }
    };


    /// Helper struct that calls Remove<ConfigType>::erase_all for all configs in a supercell
    struct EraseAllPrimPeriodicDiffTransOrbitConfigs : public ErasePrimPeriodicDiffTransOrbitConfigsBase {

      EraseAllPrimPeriodicDiffTransOrbitConfigs(Remove<Kinetics::PrimPeriodicDiffTransOrbit> &_remover, std::string _orbitname, bool _dry_run) :
        ErasePrimPeriodicDiffTransOrbitConfigsBase(_remover, _orbitname, _dry_run) {}

      void eval() {
        Remove<Kinetics::DiffTransConfiguration> f(primclex, remover.report_dir(), remover.file_log());
        DB::Selection<Kinetics::DiffTransConfiguration> selection = make_selection();
        f.erase_all(selection, dry_run);
        count_remaining();
      }

    };


    /// \brief Erase all enumerated diff_trans_configs that have no data; if no data, erase diff_trans
    void Remove<Kinetics::PrimPeriodicDiffTransOrbit>::erase(
      const DB::Selection<Kinetics::PrimPeriodicDiffTransOrbit> &selection,
      bool dry_run) {

      bool did_erase = false;

      // call Remove<Kinetics::DiffTransConfiguration>::erase
      auto it = selection.selected().begin();
      auto end = selection.selected().end();
      for(; it != end; ++it) {
        ErasePrimPeriodicDiffTransOrbitConfigs f(*this, it.name(), dry_run);
        f.eval();

        // if no existing data or files, erase PrimPeriodicDiffTransOrbit ...
        if(!f.remaining) {
          primclex().log() << "will erase " << it.name() << "\n";
          if(!dry_run) {
            primclex().db<Kinetics::PrimPeriodicDiffTransOrbit>().erase(it.name());
            did_erase = true;
          }
        }
        else {
          primclex().log() << "skipping " << it.name() << ": has "
                           << f.remaining << " configurations remaining.\n";
        }
      }

      if(did_erase) {
        primclex().db<Kinetics::PrimPeriodicDiffTransOrbit>().commit();
      }
    }



    /// \brief Erase all data from diff_trans_configs from orbitname
    void Remove<Kinetics::PrimPeriodicDiffTransOrbit>::erase_data(
      const DB::Selection<Kinetics::PrimPeriodicDiffTransOrbit> &selection,
      bool dry_run) {

      // call Remove<Kinetics::DiffTransConfiguration>::erase
      auto it = selection.selected().begin();
      auto end = selection.selected().end();
      for(; it != end; ++it) {
        EraseDataPrimPeriodicDiffTransOrbitConfigs f(*this, it.name(), dry_run);
        f.eval();
      }
    }


    /// \brief Erase all enumerated diff_trans_configs ; erase diff_trans
    void Remove<Kinetics::PrimPeriodicDiffTransOrbit>::erase_all(
      const DB::Selection<Kinetics::PrimPeriodicDiffTransOrbit> &selection,
      bool dry_run) {

      // call Remove<Kinetics::DiffTransConfiguration>::erase_all
      auto it = selection.selected().begin();
      auto end = selection.selected().end();
      for(; it != end; ++it) {
        EraseAllPrimPeriodicDiffTransOrbitConfigs f(*this, it.name(), dry_run);
        f.eval();

        // if no existing data or files, erase PrimPeriodicDiffTransOrbit ...
        if(!f.remaining) {
          primclex().log() << "will erase " << it.name() << "\n";
          if(dry_run) {
            primclex().db<Kinetics::PrimPeriodicDiffTransOrbit>().erase(it.name());
          }
        }
        else {
          primclex().log() << "unknown error: erase_all called, but " << it.name()
                           << " still has " << f.remaining << " configurations.\n";
          primclex().log() << "stopping..." << std::endl;
          return;
        }
      }

      primclex().db<Kinetics::PrimPeriodicDiffTransOrbit>().commit();
    }

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
      Remove<Kinetics::PrimPeriodicDiffTransOrbit> f(primclex, report_dir, primclex.log());

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
