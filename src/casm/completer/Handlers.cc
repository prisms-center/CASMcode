#ifndef HANDLERS_CC
#define HANDLERS_CC

#include "casm/completer/Handlers.hh"

#include "casm/app/APICommand.hh"
#include "casm/app/DirectoryStructure.hh"
#include "casm/app/ProjectSettings.hh"
#include "casm/app/enum/EnumInterface.hh"
#include "casm/app/enum/standard_enumerator_interfaces.hh"
#include "casm/app/info/InfoInterface.hh"
#include "casm/app/info/standard_info_method_interfaces.hh"
#include "casm/casm_io/Log.hh"
#include "casm/casm_io/dataformatter/DataFormatter.hh"
#include "casm/casm_io/enum/stream_io.hh"
#include "casm/clex/Configuration.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/Supercell.hh"
#include "casm/database/ConfigDatabase.hh"
#include "casm/database/DatabaseTypes.hh"
#include "casm/database/ScelDatabase.hh"
#include "casm/global/enum/io_traits.hh"

namespace CASM {
namespace Completer {
typedef ArgHandler::ARG_TYPE ARG_TYPE;

ARG_TYPE ArgHandler::determine_type(
    const po::option_description &boost_option) {
  // This string will become something like "<path>", or "arg", or "<path>
  // (=MASTER)"
  std::string raw_boost_format = boost_option.format_parameter();
  // Sometimes boost option has default arguments. We don't want to include that
  std::string argtype_str;
  std::string::size_type space_pos = raw_boost_format.find(' ');

  // Spaces found, probably printing default value as well. Strip it off.
  if (space_pos != std::string::npos) {
    argtype_str = raw_boost_format.substr(0, space_pos);
  }

  // No spaces found, so format_parameter already returned what we wanted
  else {
    argtype_str = raw_boost_format;
  }

  for (auto it = m_argument_table.begin(); it != m_argument_table.end(); ++it) {
    if (it->first == argtype_str) {
      return it->second;
    }
  }

  return ARG_TYPE::VOID;
}

std::string ArgHandler::path() { return m_argument_table[0].first; }

std::string ArgHandler::command() { return m_argument_table[1].first; }

std::string ArgHandler::supercell() { return m_argument_table[2].first; }

std::string ArgHandler::query() { return m_argument_table[3].first; }

std::string ArgHandler::operation() { return m_argument_table[4].first; }

std::string ArgHandler::configname() { return m_argument_table[5].first; }

std::string ArgHandler::coordtype() { return m_argument_table[6].first; }

std::string ArgHandler::dbtype() { return m_argument_table[7].first; }

std::string ArgHandler::enummethod() { return m_argument_table[8].first; }

std::string ArgHandler::configtype() { return m_argument_table[9].first; }

std::string ArgHandler::calctype() { return m_argument_table[10].first; }

std::string ArgHandler::bset() { return m_argument_table[11].first; }

std::string ArgHandler::clex() { return m_argument_table[12].first; }

std::string ArgHandler::ref() { return m_argument_table[13].first; }

std::string ArgHandler::eci() { return m_argument_table[14].first; }

std::string ArgHandler::property() { return m_argument_table[15].first; }

std::string ArgHandler::dof() { return m_argument_table[16].first; }

std::string ArgHandler::infomethod() { return m_argument_table[17].first; }

std::string ArgHandler::montemethod() { return m_argument_table[18].first; }

void ArgHandler::void_to_bash(std::vector<std::string> &arguments) { return; }

void ArgHandler::path_to_bash(std::vector<std::string> &arguments) {
  arguments.push_back("BASH_COMP_PATH");
  return;
}

void ArgHandler::command_to_bash(std::vector<std::string> &arguments) {
  arguments.push_back("BASH_COMP_BIN");
  return;
}

void ArgHandler::scelname_to_bash(std::vector<std::string> &arguments) {
  if (!find_casmroot(boost::filesystem::current_path()).empty()) {
    ScopedNullLogging logging;
    const PrimClex &pclex =
        PrimClex(find_casmroot(boost::filesystem::current_path()));
    for (const auto &scel : pclex.const_db<Supercell>()) {
      arguments.push_back(scel.name());
    }
  }
  return;
}

void ArgHandler::configname_to_bash(std::vector<std::string> &arguments) {
  if (!find_casmroot(boost::filesystem::current_path()).empty()) {
    ScopedNullLogging logging;
    const PrimClex &pclex =
        PrimClex(find_casmroot(boost::filesystem::current_path()));
    for (const auto &config : pclex.const_db<Configuration>()) {
      arguments.push_back(config.name());
    }
  }
  return;
}

void ArgHandler::operator_to_bash(std::vector<std::string> &arguments) {
  DataFormatterDictionary<Configuration> dict =
      make_dictionary<Configuration>();

  for (auto it = dict.begin(); it != dict.cend(); ++it) {
    if (it->type() == DatumFormatterClass::Property) {
      arguments.push_back(it->name());
    }
  }
  return;
}

void ArgHandler::query_to_bash(std::vector<std::string> &arguments) {
  DataFormatterDictionary<Configuration> config_dict =
      make_dictionary<Configuration>();

  for (auto it = config_dict.begin(); it != config_dict.cend(); ++it) {
    if (it->type() == DatumFormatterClass::Property) {
      arguments.push_back("config:" + it->name());
    }
  }
  DataFormatterDictionary<Supercell> scel_dict = make_dictionary<Supercell>();

  for (auto it = scel_dict.begin(); it != scel_dict.cend(); ++it) {
    if (it->type() == DatumFormatterClass::Property) {
      arguments.push_back("scel:" + it->name());
    }
  }

  return;
}

void ArgHandler::dbtype_to_bash(std::vector<std::string> &arguments) {
  for (auto &item : DB::types_short()) {
    arguments.push_back(item);
  }
  return;
}

void ArgHandler::enummethod_to_bash(std::vector<std::string> &arguments) {
  auto enumerator_interfaces = make_standard_enumerator_interfaces();
  for (auto const &e : enumerator_interfaces) {
    arguments.push_back(e->name());
  }
  return;
}

void ArgHandler::configtype_to_bash(std::vector<std::string> &arguments) {
  for (auto &item : DB::config_types_short()) {
    arguments.push_back(item);
  }
  return;
}

void ArgHandler::calctype_to_bash(std::vector<std::string> &arguments) {
  auto checkroot = find_casmroot(boost::filesystem::current_path());
  if (!checkroot.empty()) {
    const DirectoryStructure dir{checkroot};
    for (auto &item : dir.all_calctype()) {
      arguments.push_back(item);
    }
  }
  return;
}

void ArgHandler::bset_to_bash(std::vector<std::string> &arguments) {
  auto checkroot = find_casmroot(boost::filesystem::current_path());
  if (!checkroot.empty()) {
    const DirectoryStructure dir{checkroot};
    for (auto &item : dir.all_bset()) {
      arguments.push_back(item);
    }
  }
  return;
}

void ArgHandler::clex_to_bash(std::vector<std::string> &arguments) {
  auto checkroot = find_casmroot(boost::filesystem::current_path());
  if (!checkroot.empty()) {
    const ProjectSettings set = open_project_settings(checkroot);
    for (auto &item : set.cluster_expansions()) {
      arguments.push_back(item.first);
    }
  }
  return;
}

void ArgHandler::ref_to_bash(std::vector<std::string> &arguments) {
  auto checkroot = find_casmroot(boost::filesystem::current_path());
  if (!checkroot.empty()) {
    const ProjectSettings set = open_project_settings(checkroot);
    for (auto &item : set.dir().all_ref(set.default_clex().calctype)) {
      arguments.push_back(item);
    }
  }
  return;
}

void ArgHandler::eci_to_bash(std::vector<std::string> &arguments) {
  auto checkroot = find_casmroot(boost::filesystem::current_path());
  if (!checkroot.empty()) {
    const ProjectSettings set = open_project_settings(checkroot);
    const ClexDescription d = set.default_clex();
    for (auto &item :
         set.dir().all_eci(d.property, d.calctype, d.ref, d.bset)) {
      arguments.push_back(item);
    }
  }
  return;
}

void ArgHandler::property_to_bash(std::vector<std::string> &arguments) {
  auto checkroot = find_casmroot(boost::filesystem::current_path());
  if (!checkroot.empty()) {
    const DirectoryStructure dir{checkroot};
    for (auto &item : dir.all_property()) {
      arguments.push_back(item);
    }
  }
  return;
}

void ArgHandler::dof_to_bash(std::vector<std::string> &arguments) {
  auto checkroot = find_casmroot(boost::filesystem::current_path());
  if (!checkroot.empty()) {
    const DirectoryStructure dir{checkroot};
    // for(auto &item : dir.all_property()) {
    // arguments.push_back(item);
    //}
  }
  return;
}

void ArgHandler::infomethod_to_bash(std::vector<std::string> &arguments) {
  auto info_method_interfaces = make_standard_info_method_interfaces();
  for (auto const &e : info_method_interfaces) {
    arguments.push_back(e->name());
  }
  return;
}

/**
 * This construction right here determines what the value_name of the boost
 * options should be named. It is through these strings that bash completion can
 * know which types of completions to suggest.
 */

const std::vector<std::pair<std::string, ARG_TYPE> >
    ArgHandler::m_argument_table(
        {std::make_pair("<path>", ARG_TYPE::PATH),
         std::make_pair("<command>", ARG_TYPE::COMMAND),
         std::make_pair("<supercell>", ARG_TYPE::SCELNAME),
         std::make_pair("<query>", ARG_TYPE::QUERY),
         std::make_pair("<operation>", ARG_TYPE::OPERATOR),
         std::make_pair("<configuration>", ARG_TYPE::CONFIGNAME),
         std::make_pair("<type>", ARG_TYPE::COORDTYPE),
         std::make_pair("<dbtype>", ARG_TYPE::DBTYPE),
         std::make_pair("<enummethod>", ARG_TYPE::ENUMMETHOD),
         std::make_pair("<configtype>", ARG_TYPE::CONFIGTYPE),
         std::make_pair("<calctype>", ARG_TYPE::CALCTYPE),
         std::make_pair("<bset>", ARG_TYPE::BSET),
         std::make_pair("<clex>", ARG_TYPE::CLEX),
         std::make_pair("<ref>", ARG_TYPE::REF),
         std::make_pair("<eci>", ARG_TYPE::ECI),
         std::make_pair("<property>", ARG_TYPE::PROPERTY),
         std::make_pair("<dof>", ARG_TYPE::DOF),
         std::make_pair("<infomethod>", ARG_TYPE::INFOMETHOD),
         std::make_pair("<montemethod>", ARG_TYPE::MONTEMETHOD)});

//*****************************************************************************************************//

OptionHandlerBase::OptionHandlerBase(const std::string &init_option_tag)
    : m_tag(init_option_tag),
      m_desc(std::string("'casm ") + init_option_tag + std::string("' usage")),
      m_gzip_flag(false) {}

OptionHandlerBase::OptionHandlerBase(const std::string &init_option_tag,
                                     const std::string &init_descriptor)
    : m_tag(init_option_tag), m_desc(init_descriptor), m_gzip_flag(false) {}

const std::string &OptionHandlerBase::tag() const { return m_tag; }

/// Get the variables map
po::variables_map &OptionHandlerBase::vm() { return m_vm; }

/// Get the variables map
const po::variables_map &OptionHandlerBase::vm() const { return m_vm; }

/**
 * Check if there are any program options in the options description. If there
 * aren't, then this is the first time someone is asking for those values, which
 * we set through the initialize routine. If there are values there already,
 * just hand them back.
 */
const po::options_description &OptionHandlerBase::desc() {
  if (m_desc.options().size() == 0) {
    initialize();
  }

  return m_desc;
}

const po::options_description &OptionHandlerBase::desc() const {
  if (m_desc.options().size() == 0) {
    throw std::logic_error(
        "code error: Options must be initialized by calling non-const desc()");
  }

  return m_desc;
}

const fs::path &OptionHandlerBase::selection_path() const {
  return m_selection_path;
}

const std::vector<fs::path> &OptionHandlerBase::selection_paths() const {
  return m_selection_paths;
}

const fs::path &OptionHandlerBase::file_path() const { return m_file_path; }

const fs::path &OptionHandlerBase::prim_path() const { return m_prim_path; }

const std::string &OptionHandlerBase::verbosity_str() const {
  return m_verbosity_str;
}

/// Will throw if not expected string or int in range [0, 100]
int OptionHandlerBase::verbosity() const {
  auto val = Log::verbosity_level(verbosity_str());
  if (val.first) {
    return val.second;
  }

  throw std::invalid_argument(Log::invalid_verbosity_msg(verbosity_str()));
}

const fs::path OptionHandlerBase::settings_path() const {
  return m_settings_path;
}

std::string OptionHandlerBase::input_str() const { return m_input_str; }

const fs::path OptionHandlerBase::output_path() const { return m_output_path; }

bool OptionHandlerBase::gzip_flag() const { return m_gzip_flag; }

const std::vector<std::string> &OptionHandlerBase::help_opt_vec() const {
  return m_help_opt_vec;
}

const std::string &OptionHandlerBase::config_str() const {
  return m_config_str;
}

const std::vector<std::string> &OptionHandlerBase::config_strs() const {
  return m_config_strs;
}

const std::string &OptionHandlerBase::supercell_str() const {
  return m_supercell_str;
}

const std::vector<std::string> &OptionHandlerBase::supercell_strs() const {
  return m_supercell_strs;
}

const std::string &OptionHandlerBase::coordtype_str() const {
  return m_coordtype_str;
}

COORD_TYPE OptionHandlerBase::coordtype_enum() const {
  return from_string<COORD_TYPE>(coordtype_str());
}

/// Returns the names of the DoF type names for add_dofs_suboption()
const std::vector<std::string> &OptionHandlerBase::dof_strs() const {
  return m_dof_strs;
}

void OptionHandlerBase::add_selection_suboption(const fs::path &_default) {
  m_desc.add_options()("selection,c",
                       po::value<fs::path>(&m_selection_path)
                           ->default_value(_default)
                           ->value_name(ArgHandler::path()),
                       (std::string("Only consider the selected objects from "
                                    "the given selection file. ") +
                        singleline_help<DB::SELECTION_TYPE>())
                           .c_str());
  return;
}

/// Add --selection suboption (no default)
void OptionHandlerBase::add_selection_no_default_suboption() {
  m_desc.add_options()(
      "selection,c",
      po::value<fs::path>(&m_selection_path)->value_name(ArgHandler::path()),
      (std::string("Only consider the selected objects from the given "
                   "selection file. ") +
       standard_singleline_enum_help<DB::SELECTION_TYPE>("nodefault"))
          .c_str());
  return;
}

void OptionHandlerBase::add_configlist_suboption(const fs::path &_default) {
  m_desc.add_options()("config,c",
                       po::value<fs::path>(&m_selection_path)
                           ->default_value(_default)
                           ->value_name(ArgHandler::path()),
                       (std::string("Only consider the selected configurations "
                                    "of the given selection file. ") +
                        standard_singleline_enum_help<DB::SELECTION_TYPE>(_default.string(), "filename"))
                           .c_str());
    return;
}

void OptionHandlerBase::add_selections_suboption(const fs::path &_default) {
  m_desc.add_options()("selections,c",
                       po::value<std::vector<fs::path> >(&m_selection_paths)
                           ->default_value(std::vector<fs::path>{_default})
                           ->value_name(ArgHandler::path()),
                       (std::string("Only consider the selected objects from "
                                    "the given selection file. ") +
                        singleline_help<DB::SELECTION_TYPE>())
                           .c_str());
  return;
}

void OptionHandlerBase::add_configlists_suboption(const fs::path &_default) {
  m_desc.add_options()("configs,c",
                       po::value<std::vector<fs::path> >(&m_selection_paths)
                           ->default_value(std::vector<fs::path>{_default})
                           ->value_name(ArgHandler::path()),
                       (std::string("Only consider the selected configurations "
                                    "of the given selection file. ") +
                        singleline_help<DB::SELECTION_TYPE>())
                           .c_str());
  return;
}

void OptionHandlerBase::add_configlist_nodefault_suboption() {
  m_desc.add_options()(
      "config,c",
      po::value<fs::path>(&m_selection_path)->value_name(ArgHandler::path()),
      (std::string("Only consider the selected configurations of the given "
                   "selection file. ") +
       standard_singleline_enum_help<DB::SELECTION_TYPE>("", "filename"))
          .c_str());
  return;
}

void OptionHandlerBase::add_configlists_nodefault_suboption() {
  m_desc.add_options()(
      "configs,c",
      po::value<std::vector<fs::path> >(&m_selection_paths)
          ->value_name(ArgHandler::path()),
      (std::string("Only consider the selected configurations of the given "
                   "selection files. ") +
       standard_singleline_enum_help<DB::SELECTION_TYPE>("", "filename"))
          .c_str());
  return;
}

void OptionHandlerBase::add_file_path_suboption(const fs::path &_default) {
  m_desc.add_options()(
      "path",
      po::value<fs::path>(&m_file_path)
          ->default_value(_default)
          ->value_name(ArgHandler::path()),
      std::string("Path to directory in which to run the command. ").c_str());
  return;
}

void OptionHandlerBase::add_prim_path_suboption(const fs::path &_default) {
  m_desc.add_options()("prim",
                       po::value<fs::path>(&m_prim_path)
                           ->default_value(_default)
                           ->value_name(ArgHandler::path()),
                       std::string("Path to prim.json or POSCAR-like file that "
                                   "defines project reference crystal. ")
                           .c_str());
  return;
}

void OptionHandlerBase::add_configtype_suboption(
    std::string _default, std::set<std::string> _configtype_opts) {
  if (!_configtype_opts.size()) {
    m_configtype_opts = DB::config_types_short();
  } else {
    m_configtype_opts = _configtype_opts;
  }

  std::stringstream help;
  help << "Type of configurations. "
       << standard_singleline_help(m_configtype_opts, _default) << ".";

  m_desc.add_options()("type,t",
                       po::value<std::string>(&m_configtype)
                           ->default_value(_default)
                           ->value_name(ArgHandler::configtype()),
                       help.str().c_str());
}

std::string OptionHandlerBase::configtype() const { return m_configtype; }

std::set<std::string> OptionHandlerBase::configtype_opts() const {
  return m_configtype_opts;
}

void OptionHandlerBase::add_db_type_suboption(
    std::string _default, std::set<std::string> _db_type_opts) {
  if (!_db_type_opts.size()) {
    m_db_type_opts = DB::types_short();
  } else {
    m_db_type_opts = _db_type_opts;
  }

  std::stringstream help;
  help << "Type of database objects. "
       << standard_singleline_help(m_db_type_opts, _default) << ".";

  m_desc.add_options()("type,t",
                       po::value<std::string>(&m_db_type)
                           ->default_value(_default)
                           ->value_name(ArgHandler::dbtype()),
                       help.str().c_str());
}

std::string OptionHandlerBase::db_type() const { return m_db_type; }

std::set<std::string> OptionHandlerBase::db_type_opts() const {
  return m_db_type_opts;
}

void OptionHandlerBase::add_help_suboption() {
  m_desc.add_options()("help,h", "Print help message")(
      "desc", "Print extended usage description");
  return;
}

void OptionHandlerBase::add_general_help_suboption() {
  m_desc.add_options()(
      "help,h",
      po::value<std::vector<std::string> >(&m_help_opt_vec)
          ->multitoken()
          ->zero_tokens(),
      "Print general help. Use '--help properties' for a list of query-able "
      "properties or '--help operators' for a list of query operators");
  return;
}

void OptionHandlerBase::add_verbosity_suboption() {
  // TODO: add ArgHandler for this
  m_desc.add_options()(
      "verbosity",
      po::value<std::string>(&m_verbosity_str)->default_value("standard"),
      "Verbosity of output. Options are 'none', 'quiet', 'standard', "
      "'verbose', 'debug', or an integer 0-100 (0: none, 100: debug).");
  return;
}

void OptionHandlerBase::add_settings_suboption(bool required) {
  std::string help_str =
      "Settings input file specifying which parameters should be used. See "
      "'casm format --" +
      tag() + "'.";

  if (required) {
    m_desc.add_options()("settings,s",
                         po::value<fs::path>(&m_settings_path)
                             ->required()
                             ->value_name(ArgHandler::path()),
                         help_str.c_str());
  } else {
    m_desc.add_options()(
        "settings,s",
        po::value<fs::path>(&m_settings_path)->value_name(ArgHandler::path()),
        help_str.c_str());
  }

  return;
}

void OptionHandlerBase::add_input_suboption(bool required) {
  std::string help_str =
      "String specifying input settings. See 'casm format --" + tag() + "'.";

  if (required) {
    m_desc.add_options()("input,i",
                         po::value<std::string>(&m_input_str)->required(),
                         help_str.c_str());
  } else {
    m_desc.add_options()("input,i", po::value<std::string>(&m_input_str),
                         help_str.c_str());
  }

  return;
}

void OptionHandlerBase::add_output_suboption() {
  m_desc.add_options()(
      "output,o",
      po::value<fs::path>(&m_output_path)->value_name(ArgHandler::path()),
      "Name for output file. Use STDOUT to print results without extra "
      "messages.");
  return;
}

void OptionHandlerBase::add_output_suboption(const fs::path &_default) {
  m_desc.add_options()("output,o",
                       po::value<fs::path>(&m_output_path)
                           ->default_value(_default)
                           ->value_name(ArgHandler::path()),
                       (std::string("Name for output file. ") +
                        "Use STDOUT to print results without extra messages."
                        "If not specified, '" +
                        _default.string() + "' will be used.")
                           .c_str());
  return;
}

void OptionHandlerBase::add_gzip_suboption() {
  m_desc.add_options()("gzip,z", po::value(&m_gzip_flag)->zero_tokens(),
                       "Write gzipped output file.");
  return;
}

void OptionHandlerBase::add_scelname_suboption() {
  std::string help_str;
  help_str = "Single supercell name to use casm " + m_tag +
             " with, such as 'SCEL4_2_2_1_0_0_0'";
  m_desc.add_options()("scelname",
                       po::value<std::string>(&m_supercell_str)
                           ->value_name(ArgHandler::supercell()),
                       help_str.c_str());
  return;
}

void OptionHandlerBase::add_scelnames_suboption() {
  std::string help_str;
  help_str = "One or more supercells to use casm " + m_tag +
             " with, such as 'SCEL4_2_2_1_0_0_0'";
  m_desc.add_options()("scelnames",
                       po::value<std::vector<std::string> >(&m_supercell_strs)
                           ->multitoken()
                           ->value_name(ArgHandler::supercell()),
                       help_str.c_str());
  return;
}

void OptionHandlerBase::add_configname_suboption() {
  std::string help_str;

  help_str = "Single configuration to use casm " + m_tag +
             " with, such as 'SCEL4_2_2_1_0_0_0/3'";

  m_desc.add_options()("configname",
                       po::value<std::string>(&m_config_str)
                           ->value_name(ArgHandler::configname()),
                       help_str.c_str());

  return;
}

void OptionHandlerBase::add_confignames_suboption() {
  std::string help_str;

  help_str = "One or more configurations to use casm " + m_tag +
             " with, such as 'SCEL4_2_2_1_0_0_0/3'";

  m_desc.add_options()("confignames",
                       po::value<std::vector<std::string> >(&m_config_strs)
                           ->multitoken()
                           ->value_name(ArgHandler::configname()),
                       help_str.c_str());

  return;
}

void OptionHandlerBase::add_names_suboption() {
  std::string help_str;

  help_str = "One or more object names to use casm " + m_tag +
             " with, such as 'SCEL4_2_2_1_0_0_0/3'";

  m_desc.add_options()("names",
                       po::value<std::vector<std::string> >(&m_name_strs)
                           ->multitoken()
                           ->value_name(ArgHandler::configname()),
                       help_str.c_str());

  return;
}

const std::vector<std::string> &OptionHandlerBase::name_strs() const {
  return m_name_strs;
}

void OptionHandlerBase::add_coordtype_suboption() {
  m_desc.add_options()("coord",
                       po::value<std::string>(&m_coordtype_str)
                           ->default_value("frac")
                           ->value_name(ArgHandler::coordtype()),
                       (std::string("Type of coordinate system to use. ") +
                        singleline_help<COORD_TYPE>())
                           .c_str());
  return;
}

/// Add a --dofs suboption to specify DoF Types
void OptionHandlerBase::add_dofs_suboption() {
  std::string help_str;
  help_str = "One or more DoF types to use casm " + m_tag +
             " with, such as 'disp' or 'EAstrain'";
  m_desc.add_options()("dofs",
                       po::value<std::vector<std::string> >(&m_dof_strs)
                           ->multitoken()
                           ->value_name(ArgHandler::dof()),
                       help_str.c_str());
  return;
}

void OptionHandlerBase::add_dry_run_suboption(std::string msg) {
  m_desc.add_options()("dry-run,n", msg.c_str());
}

std::string OptionHandlerBase::default_dry_run_msg() {
  return std::string("Dry run (if supported by method)");
}

bool OptionHandlerBase::dry_run() const { return vm().count("dry-run"); }

//*****************************************************************************************************//

}  // namespace Completer
}  // namespace CASM

#endif
