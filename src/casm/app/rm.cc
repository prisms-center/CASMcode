#include "casm/crystallography/CoordinateSystems.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/app/ProjectSettings.hh"
#include "casm/app/DirectoryStructure.hh"
#include "casm/app/AppIO.hh"
#include "casm/app/casm_functions.hh"
#include "casm/clex/ConfigSelection.hh"

#include "casm/completer/Handlers.hh"

namespace CASM {

  namespace Completer {
    RmOption::RmOption(): OptionHandlerBase("rm") {}

    bool RmOption::force() const {
      return vm().count("force");
    }

    bool RmOption::dry_run() const {
      return vm().count("dry-run");
    }

    void RmOption::initialize() {
      add_help_suboption();
      add_scelnames_suboption();
      add_confignames_suboption();
      add_configlist_nodefault_suboption();

      m_desc.add_options()
      ("dry-run,n", "Dry run")
      ("force,f", "Force remove.");

      return;
    }
  }

  int _rm_configs(const CommandArgs &args, const Completer::RmOption &rm_opt);
  int _rm_scel(const CommandArgs &args, const Completer::RmOption &rm_opt);

  // ///////////////////////////////////////
  // 'rm' function for casm

  int rm_command(const CommandArgs &args) {


    /// Set command line options using boost program_options
    Completer::RmOption rm_opt;
    po::variables_map &vm = rm_opt.vm();

    bool rm_config = false;
    bool rm_scel = false;

    try {
      po::store(po::parse_command_line(args.argc, args.argv, rm_opt.desc()), vm); // can throw

      bool call_help = false;

      // allow --config && --confignames   or  --scelnames,  but not both
      if(!vm.count("help") && !vm.count("desc")) {

        rm_config = vm.count("config") + vm.count("confignames");
        rm_scel = vm.count("scelnames");

        if(!rm_config && !rm_scel) {
          args.log << "Error in 'casm rm': At least one of --config, "
                   "--confignames, --scelnames must be given." << std::endl;
          call_help = true;
        }

        if(rm_config + rm_scel > 1) {
          args.log << "Error in 'casm rm': May not delete multiple types at the same time." << std::endl;
          call_help = true;
        }
      }

      /** --help option
      */
      if(vm.count("help") || call_help) {
        args.log << "\n";
        args.log << rm_opt.desc() << std::endl;

        return 0;
      }

      if(vm.count("desc")) {
        args.log << "\n";
        args.log << rm_opt.desc() << std::endl;
        args.log << "DESCRIPTION" << std::endl;
        args.log << "  Remove configuration data and supercells.\n"

                 "  Can remove: \n"
                 "  - A supercell, including all enumerated configurations   \n"
                 "    and data (for all calctypes) via --scelnames.\n"
                 "  - Calculation data (for the current calctype) for \n"
                 "    individual configurations via --config or --confignames.\n\n"

                 "  Cannot remove: \n"
                 "  - Specific enumerated configuration. The current version \n"
                 "    of CASM is restricted to consequitively numbered       \n"
                 "    configurations in each supercell, so configurations may\n"
                 "    only be erased by erasing a supercell.\n\n";

        return 0;
      }

      po::notify(vm); // throws on error, so do after help in case
      // there are any problems

    }
    catch(po::error &e) {
      args.err_log << "ERROR: " << e.what() << std::endl << std::endl;
      args.err_log << rm_opt.desc() << std::endl;
      return ERR_INVALID_ARG;
    }
    catch(std::exception &e) {
      args.err_log << "Unhandled Exception reached the top of main: "
                   << e.what() << ", application will now exit" << std::endl;
      return ERR_UNKNOWN;

    }

    const fs::path &root = args.root;
    if(root.empty()) {
      args.err_log.error("No casm project found");
      args.err_log << std::endl;
      return ERR_NO_PROJ;
    }

    // Only allow erasing one type at a time (check is above to call help):
    // - configurations: --config / --confignames
    // - supercells: --scelnames
    if(rm_config) {
      return _rm_configs(args, rm_opt);
    }
    else if(rm_scel) {
      return _rm_scel(args, rm_opt);
    }
    return ERR_INVALID_ARG;
  }

  int _rm_configs(const CommandArgs &args, const Completer::RmOption &rm_opt) {

    // If 'args.primclex', use that, else construct PrimClex in 'uniq_primclex'
    // Then whichever exists, store reference in 'primclex'
    std::unique_ptr<PrimClex> uniq_primclex;
    PrimClex &primclex = make_primclex_if_not(args, uniq_primclex);

    // if selection of configurations
    std::unique_ptr<ConstConfigSelection> selection;
    if(rm_opt.vm().count("config")) {
      selection = notstd::make_unique<ConstConfigSelection>(primclex, rm_opt.selection_path());
    }
    else if(rm_opt.vm().count("confignames")) {
      selection = notstd::make_unique<ConstConfigSelection>(primclex, fs::path("NONE"));
      for(const auto &configname : rm_opt.config_strs()) {
        selection->set_selected(configname, true);
      }
    }

    if(rm_opt.dry_run()) {
      args.log << "dry-run...\n\n";
    }

    for(auto it = selection->selected_config_begin(); it != selection->selected_config_end(); ++it) {
      args.log.custom(std::string("Erase ") + it->name() + " data");
      if(rm_opt.force() && fs::exists(it->calc_dir())) {
        recurs_rm_files(it->calc_dir(), rm_opt.dry_run(), args.log);
      }
      else if(fs::exists(it->calc_dir())) {
        args.log << "skipping " << it->name() << " data: --force option not given\n";
      }
      else {
        args.log << "skipping " << it->name() << " data: no data\n";
      }
      args.log << "\n";
    }

    return 0;
  }

  int _rm_scel(const CommandArgs &args, const Completer::RmOption &rm_opt) {

    // If 'args.primclex', use that, else construct PrimClex in 'uniq_primclex'
    // Then whichever exists, store reference in 'primclex'
    std::unique_ptr<PrimClex> uniq_primclex;
    PrimClex &primclex = make_primclex_if_not(args, uniq_primclex);

    // collect supercells to delete
    std::set<std::string> scel_to_delete;

    args.log.custom("Check supercells");
    for(const auto &scelname : rm_opt.supercell_strs()) {
      Index index;
      if(!primclex.contains_supercell(scelname, index)) {
        args.log << "skipping " << scelname << ": does not exist.\n";
        continue;
      }

      const Supercell &scel = primclex.get_supercell(index);

      if(!rm_opt.force() && scel.get_config_list().size()) {
        args.log << "skipping " << scelname << ": has "
                 << scel.get_config_list().size() << " configurations and --force not given.\n";
      }
      else {
        args.log << "will erase " << scelname << "\n";
        scel_to_delete.insert(scelname);
      }
    }

    args.log << std::endl;

    // Erase supercells from SCEL and config_list.json
    if(!scel_to_delete.size()) {
      args.log << "No supercells to erase\n";
    }
    else {
      for(const auto &scelname : scel_to_delete) {
        args.log.custom(std::string("Erase ") + scelname);
        const Supercell &scel = primclex.get_supercell(scelname);
        recurs_rm_files(scel.get_path(), rm_opt.dry_run(), args.log);
        args.log << "\n";
      }

      if(!rm_opt.dry_run()) {
        primclex.print_supercells(scel_to_delete);
        primclex.write_config_list(scel_to_delete);
      }
    }
    return 0;

  };

}
