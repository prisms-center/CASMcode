#include <cstring>

#include "casm/app/casm_functions.hh"
#include "casm/misc/algorithm.hh"
#include "casm/casm_io/FileEnumerator.hh"
#include "casm/completer/Handlers.hh"

namespace CASM {

  namespace Completer {
    FilesOption::FilesOption(): OptionHandlerBase("files") {}

    const std::vector<std::string> &FilesOption::calc_vec() const {
      return m_calc_vec;
    }

    const std::string &FilesOption::settings_str() const {
      return m_settings_str;
    }

    void FilesOption::initialize() {
      add_help_suboption();
      add_output_suboption();
      add_gzip_suboption();

      m_desc.add_options()
      ("settings",
       po::value<std::string>(&m_settings_str)->default_value("curr"),
       "Which settings to include. "
       "One of: 'curr' (default), 'all'.")

      ("calc",
       po::value<std::vector<std::string> >(&m_calc_vec)->multitoken(),
       "Which calculation files to include. "
       "May be zero (default) or more of: 'settings', 'status', 'all'.") //TODO: Include these in ArgHandler

      ("relative,R",
       "Print relative path from project root directory.");
    }
  }


  // ///////////////////////////////////////
  // 'files' function for casm
  //    (add an 'if-else' statement in casm.cpp to call this)

  int files_command(const CommandArgs &args) {

    po::variables_map vm;
    fs::path out_path;
    std::string settings;
    std::vector<std::string> calc;  //These variable names are bad and you should feel bad
    bool gz_flag;

    /// Set command line options using boost program_options
    Completer::FilesOption files_opt;

    try {
      po::store(po::parse_command_line(args.argc, args.argv, files_opt.desc()), vm); // can throw

      if(!vm.count("help")) {

        std::vector<std::string> allowed;

        allowed = std::vector<std::string> {"", "curr", "all"};
        if(!contains(allowed, settings)) {
          args.err_log << "Error in 'casm files'. '" << settings << "' is not a valid option to --settings." << std::endl;
          args.err_log << "  Valid options are: 'curr', 'all'" << std::endl;
          return ERR_INVALID_ARG;
        }

        allowed = std::vector<std::string> {"settings", "status", "all"};
        for(auto it = calc.begin(); it != calc.end(); ++it) {
          if(!contains(allowed, *it)) {
            args.err_log << "Error in 'casm files'. '" << *it << "' is not a valid option to --calc." << std::endl;
            args.err_log << "  Valid options are: 'settings', 'status', 'all'" << std::endl;
            return ERR_INVALID_ARG;
          }
        }
      }

      /** --help option
       */
      if(vm.count("help")) {
        args.log << "\n";
        args.log << files_opt.desc() << std::endl;

        return 0;
      }

      if(vm.count("desc")) {
        args.log << "\n";
        args.log << files_opt.desc() << std::endl;

        args.log << "DESCRIPTION  \n"
                 "    Enumerate files used by this CASM project\n"
                 "    - If --all given, include files from all settings (bset, \n"
                 "      calctype, ref, clex, eci). Otherwise only include files\n"
                 "      relevant to the current settings.\n"
                 "    - If --calc given, include training data files. This option\n"
                 "      accepts one or more of: 'settings', 'status', or 'all'.\n"
                 "      'settings': include files from the 'training_data' \n"
                 "                  settings directories.\n"
                 "      'status':   include 'properties.calc.json' and     \n"
                 "                  'status.json' files.\n"
                 "      'all':      include all files in the 'training_data'\n"
                 "                  directory, recursively.\n\n"

                 "    Examples:\n"
                 "      casm files\n"
                 "      - Prints basic project files for the current settings.\n"
                 "      - Prints absolute paths.\n\n"

                 "      casm files -R\n"
                 "      - Prints basic project files for the current settings.\n"
                 "      - Prints relative paths from project root directory.\n\n"

                 "      casm files -o absfiles.txt\n"
                 "      - Write 'absfiles.txt' with basic project files for \n"
                 "        the current settings.\n"
                 "      - Prints absolute paths.\n\n"

                 "      casm files -R -o relfiles.txt\n"
                 "      - Write 'relfiles.txt' with basic project files for \n"
                 "        the current settings.\n"
                 "      - Prints relative paths from project root directory.\n\n"

                 "      tar -czvf proj.tar.gz `casm files -o STDOUT` \n"
                 "      tar -czvf proj.tar.gz -T absfiles.txt \n"
                 "      - Use tar to create an archive of your CASM project. \n"
                 "      - Extract with 'tar -xzvf proj.tar.gz'. \n\n"

                 "      rsync -avPR `casm files -R -o STDOUT` destination -n \n"
                 "      - Copy project files to the 'destination' directory, \n"
                 "        creating sub-directories as necessary.             \n"
                 "      - **Must run from project root directory**           \n"
                 "      - Replace 'destination' with '-e ssh user@host:destination' \n"
                 "        for a remote destination.\n"
                 "      - \"-n\" option is for a \"dry-run\". Remove it when ready \n"
                 "        to do the transfer.\n\n"

                 "      rsync -avPR --files-from=relfiles.txt src destination -n \n"
                 "      - Copy project files to the 'destination' directory, \n"
                 "        creating sub-directories as necessary.             \n"
                 "      - 'src' is the path to the root directory of the CASM\n"
                 "        project you are copying from.                      \n"
                 "      - Replace 'destination' with '-e ssh user@host:destination' \n"
                 "        for a remote destination.\n"
                 "      - \"-n\" option is for a \"dry-run\". Remove it when ready \n"
                 "        to do the transfer.\n";

        return 0;
      }


      po::notify(vm); // throws on error, so do after help in case
      // there are any problems

      out_path = files_opt.output_path();
      settings = files_opt.settings_str();
      calc = files_opt.calc_vec();
      gz_flag = files_opt.gzip_flag();

    }
    catch(po::error &e) {
      args.err_log << files_opt.desc() << std::endl;
      args.err_log << "ERROR: " << e.what() << std::endl << std::endl;
      return ERR_INVALID_ARG;
    }
    catch(std::exception &e) {
      args.err_log << files_opt.desc() << std::endl;
      args.err_log << "ERROR: " << e.what() << std::endl << std::endl;
      return ERR_UNKNOWN;
    }

    auto check_gz = [ = ](fs::path p) {
      if(p.extension() == ".gz" || p.extension() == ".GZ") {
        return true;
      }
      return false;
    };

    if(check_gz(out_path)) {
      gz_flag = true;
    }

    const fs::path &root = args.root;
    if(root.empty()) {
      args.err_log.error("No casm project found");
      args.err_log << std::endl;
      return ERR_NO_PROJ;
    }


    // set output_stream: where the file paths are written
    std::unique_ptr<std::ostream> uniq_fout;
    std::ostream &output_stream = make_ostream_if(vm.count("output"), args.log, uniq_fout, out_path, gz_flag);

    // set status_stream: where query settings and PrimClex initialization messages are sent
    Log &status_log = (out_path.string() == "STDOUT") ? args.err_log : args.log;

    // If '_primclex', use that, else construct PrimClex in 'uniq_primclex'
    // Then whichever exists, store reference in 'primclex'
    std::unique_ptr<PrimClex> uniq_primclex;
    PrimClex &primclex = make_primclex_if_not(args, uniq_primclex, status_log);

    std::vector<fs::path> files;
    auto result = std::back_inserter(files);
    FileEnumerator enumerate(primclex, settings == "all", vm.count("relative"));

    result = enumerate.basic_files(result);
    result = enumerate.bset_files(result);
    result = enumerate.reference_files(result);
    result = enumerate.eci_files(result);

    bool calc_settings = contains(calc, "settings");
    bool calc_status = contains(calc, "status");
    bool calc_all = contains(calc, "all");

    if(calc_settings && !calc_all) {
      result = enumerate.calc_settings_files(result);
    }
    if(calc_status && !calc_all) {
      result = enumerate.calc_status_files(result);
    }
    if(calc_all) {
      result = enumerate.all_calc_files(result);
    }

    for(auto f : files) {
      output_stream << f.string() << "\n";
    }

    return 0;
  };

}

