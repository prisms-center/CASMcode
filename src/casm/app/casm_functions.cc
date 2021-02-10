#include "casm/app/casm_functions.hh"

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

#include "casm/app/APICommand.hh"
#include "casm/app/APICommand_impl.hh"
#include "casm/app/DirectoryStructure.hh"
#include "casm/app/bset.hh"
#include "casm/app/composition.hh"
#include "casm/app/enum.hh"
#include "casm/app/files.hh"
#include "casm/app/format.hh"
#include "casm/app/import.hh"
#include "casm/app/info.hh"
#include "casm/app/init.hh"
#include "casm/app/monte.hh"
#include "casm/app/query.hh"
#include "casm/app/ref.hh"
#include "casm/app/rm.hh"
#include "casm/app/run.hh"
#include "casm/app/select.hh"
#include "casm/app/settings.hh"
#include "casm/app/status.hh"
#include "casm/app/super.hh"
#include "casm/app/sym.hh"
#include "casm/app/update.hh"
#include "casm/app/view.hh"
#include "casm/casm_io/Log.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/external/gzstream/gzstream.h"
#include "casm/misc/algorithm.hh"
#include "casm/version/version.hh"

namespace CASM {

void print_splash(std::ostream &out) {
  out << "       .::::::::.        .:::::.         .::::::.          .:.       "
         ".:.     \n"
      << "     .:::::::::::.      .:::::::.      .::::::::::.       .:::.     "
         ".:::.    \n"
      << "    .:::'     ':::.    .:::' ':::.    .:::'   '::::      .:::::.   "
         ".:::::.   \n"
      << "    ::::       ::::   .:::'   ':::.   ::::     '::'     .:::::::. "
         ".:::::::.  \n"
      << "    ::::              ::::     ::::   '::::.            "
         "::::'':::.:::''::::  \n"
      << "    ::::              ::::     ::::    '::::::.         ::::  "
         "':::::'  ::::  \n"
      << "    ::::              :::::::::::::      ''::::::.      ::::   ':::' "
         "  ::::  \n"
      << "    ::::              :::::::::::::         '::::::     ::::    ':'  "
         "  ::::  \n"
      << "    ::::      .::::   ::::     ::::            :::::    ::::         "
         "  ::::  \n"
      << "    ':::.    .::::'   ::::     ::::   .::.     .::::    ::::         "
         "  ::::  \n"
      << "     '::::::::::'     ::::     ::::   :::::...:::::'    ::::         "
         "  ::::  \n"
      << "        ':::::'       '::'     '::'   ':::::::::::'     '::'         "
         "  '::'  \n";
}

/// \brief CommandArgs constructor
///
/// \param _argc int, as from main
/// \param _argv char*[], as from main
/// \param _primclex pointer to PrimClex or nullptr
/// \param _root location of CASM project. If empty path, will use root of the
///        CASM project containing current working directory
///
CommandArgs::CommandArgs(int _argc, char *_argv[], PrimClex *_primclex,
                         fs::path _root)
    : CLIParse(_argc, _argv), primclex(_primclex), root(_root) {
  _init();
}

/// \brief CommandArgs constructor
///
/// \param _args std::string of form 'casm [subcommand] [opt...]'
/// \param _primclex pointer to PrimClex or nullptr
/// \param _root location of CASM project. If empty path, will use root of the
///        CASM project containing current working directory
///
CommandArgs::CommandArgs(std::string _args, PrimClex *_primclex, fs::path _root)
    : CLIParse(_args), primclex(_primclex), root(_root) {
  _init();
}

CommandArgs::~CommandArgs() {}

void CommandArgs::_init() {
  // set project 'root' if not already set
  if (root.empty()) {
    root = find_casmroot(fs::current_path());
  }

  // set 'command'
  command = (argc() > 1) ? std::string(argv()[1]) : "";

  if (command == "--version") {
    command = "version";
  }

  // check if 'help' command
  std::vector<std::string> help_commands{
      "help", "-h", "--help", "version", "--version", "status", "format"};
  is_help = contains(help_commands, command);

  // check if LOG should be written
  write_log = !(is_help || root.empty() || command == "init");
}

/// \brief Return static CommandMap containing all CASM API commands
CommandMap &command_map() {
  static CommandMap command_map{
      {"version", version_command},
      {"status", status_command},
      {"format", format_command},
      {"init", init_command},
      {"info", run_api_command<InfoCommand>},
      {"settings", settings_command},
      {SymCommand::name, run_api_command<SymCommand>},
      {"composition", composition_command},
      {"ref", ref_command},
      {UpdateCommand::name, run_api_command<UpdateCommand>},
      {EnumCommand::name, run_api_command<EnumCommand>},
      {"super", super_command},
      {SelectCommand::name, run_api_command<SelectCommand>},
      {"bset", run_api_command<BsetCommand>},
      {"run", run_command},
      {RmCommand::name, run_api_command<RmCommand>},
      {QueryCommand::name, run_api_command<QueryCommand>},
      {"files", files_command},
      {ImportCommand::name, run_api_command<ImportCommand>},
      {"monte", monte_command},
      {"view", view_command},
      {"help", help_command}};

  return command_map;
}

namespace api_impl {

std::string date_time() {
  auto t = std::time(nullptr);
  char str[255];
  strftime(str, sizeof(str), "%Y-%m-%d %H:%M:%S", std::localtime(&t));
  return std::string(str);
};

void write_LOG_begin(const CommandArgs &args) {
  // If not a 'version', or 'help' command, write to LOG
  fs::ofstream log(args.root / "LOG", std::ofstream::out | std::ofstream::app);
  log << "# " << date_time() << std::endl;

  // record whoami@hostname

  std::string whoami, hostname;
  {
    FILE *fp;
    char path[PATH_MAX];
    fp = popen("whoami", "r");
    while (fgets(path, PATH_MAX, fp) != NULL) {
      whoami = std::string(path);
      whoami.resize(whoami.size() - 1);
    }
  }
  {
    FILE *fp;
    char path[PATH_MAX];
    fp = popen("hostname", "r");
    while (fgets(path, PATH_MAX, fp) != NULL) {
      hostname = std::string(path);
      hostname.resize(hostname.size() - 1);
    }
  }
  log << "# " << whoami << "@" << hostname << "\n";

  // record CASM version number
  log << "# " << version() << "\n";

  // record command arguments
  for (int i = 0; i < args.argc(); i++) {
    log << args.argv()[i] << " ";
  }
  log << "\n";
  log.close();
}

void write_LOG_end(const CommandArgs &args, int retcode) {
  fs::ofstream log(args.root / "LOG", std::ofstream::out | std::ofstream::app);
  log << "# return: " << retcode << " runtime(s): " << CASM::log().time_s()
      << "\n\n";
  log.close();
}

}  // namespace api_impl

/// \brief Executes CASM commands specified by args
int casm_api(const CommandArgs &args) {
  if (args.argc() == 1) {
    help_command(args);
    return ERR_INVALID_ARG;
  }

  auto it = command_map().find(args.command);
  if (it != command_map().end()) {
    if (args.write_log) {
      api_impl::write_LOG_begin(args);
    }

    CASM::log().restart_clock();
    int retcode = it->second(args);
    if (args.write_log) {
      api_impl::write_LOG_end(args, retcode);
    }
    return retcode;
  } else {
    help_command(args);
    return ERR_INVALID_ARG;
  }
}

/// \brief If !_primclex, construct new PrimClex stored in uniq_primclex, then
///        return reference to existing or constructed PrimClex
///
/// \param args CommandArgs reference
/// \param uniq_primclex Reference to null std::unique_ptr<PrimClex> to manage
/// PrimClex, if it is constructed
///
/// \returns reference to PrimClex (either newly constructed managed by
///          uniq_primclex, or existing pointed at by args.primclex)
///
PrimClex &make_primclex_if_not(const CommandArgs &args,
                               std::unique_ptr<PrimClex> &uniq_primclex) {
  if (!args.primclex) {
    uniq_primclex.reset(new PrimClex(args.root));
    return *uniq_primclex;
  }
  return *args.primclex;
}

/// \brief Return a reference to proper std::ostream
///
/// \param output Output mode: False: use 'sout', True: check 'out_path' and
/// 'gzip' to decide \param sout stream to use if not writing to file \param
/// fout will be given an open file if writing to file \param out_path, where to
/// write if 'output': if "STDOUT", use 'sout'; otherwise filename \param gzip:
/// if true, write to gzip file
///
/// \return reference to stream to use
///
std::ostream &make_ostream_if(bool output, std::ostream &sout,
                              std::unique_ptr<std::ostream> &fout,
                              fs::path out_path, bool gzip) {
  if (output) {
    if (out_path.string() == "STDOUT") {
      return sout;
    }

    out_path = fs::absolute(out_path);

    if (gzip) {
      fout.reset(new gz::ogzstream(out_path.string().c_str()));
      return *fout;
    }

    fout.reset(new fs::ofstream(out_path));
    return *fout;
  } else {
    return sout;
  }
}

/// \brief Print CASM help info to log()
int help_command(const CommandArgs &args) {
  log().custom("casm usage");
  log() << "\n";

  log() << "casm [--version] <command> [options] [args]" << std::endl
        << std::endl;
  log() << "available commands:" << std::endl;

  std::vector<std::string> subcom;
  for (auto it = command_map().begin(); it != command_map().end(); ++it) {
    std::string s = it->first;
    // check for --version
    if (s[0] != '-') {
      subcom.push_back(std::string("  ") + s);
    }
  }

  std::sort(subcom.begin(), subcom.end());
  for (auto it = subcom.begin(); it != subcom.end(); ++it) {
    log() << *it << "\n";
  }
  log() << "\n";

  log() << "For help using a command: 'casm <command> --help'" << std::endl
        << std::endl;
  log() << "For step by step help use: 'casm status -n'" << std::endl
        << std::endl;

  return 0;
};

int version_command(const CommandArgs &args) {
  log() << "casm version: " << version() << std::endl;
  return 0;
}

}  // namespace CASM
