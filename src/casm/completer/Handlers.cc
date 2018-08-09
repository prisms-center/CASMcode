#ifndef HANDLERS_CC
#define HANDLERS_CC

#include "casm/completer/Handlers.hh"
#include "casm/casm_io/DataFormatter.hh"
#include "casm/clex/Configuration.hh"

namespace CASM {
  namespace Completer {
    typedef ArgHandler::ARG_TYPE ARG_TYPE;


    ARG_TYPE ArgHandler::determine_type(const po::option_description &boost_option) {
      //This string will become something like "<path>", or "arg", or "<path> (=MASTER)"
      std::string raw_boost_format = boost_option.format_parameter();
      //Sometimes boost option has default arguments. We don't want to include that
      std::string argtype_str;
      std::string::size_type space_pos = raw_boost_format.find(' ');

      //Spaces found, probably printing default value as well. Strip it off.
      if(space_pos != std::string::npos) {
        argtype_str = raw_boost_format.substr(0, space_pos);
      }

      //No spaces found, so format_parameter already returned what we wanted
      else {
        argtype_str = raw_boost_format;
      }

      for(auto it = m_argument_table.begin(); it != m_argument_table.end(); ++it) {
        if(it->first == argtype_str) {
          return it->second;
        }
      }

      return ARG_TYPE::VOID;
    }

    std::string ArgHandler::path() {
      return m_argument_table[0].first;
    }

    std::string ArgHandler::command() {
      return m_argument_table[1].first;
    }

    std::string ArgHandler::supercell() {
      return m_argument_table[2].first;
    }

    std::string ArgHandler::query() {
      return m_argument_table[3].first;
    }

    std::string ArgHandler::operation() {
      return m_argument_table[4].first;
    }

    std::string ArgHandler::configname() {
      return m_argument_table[5].first;
    }

    std::string ArgHandler::coordtype() {
      return m_argument_table[6].first;
    }


    void ArgHandler::void_to_bash(std::vector<std::string> &arguments) {
      return;
    }

    void ArgHandler::path_to_bash(std::vector<std::string> &arguments) {
      arguments.push_back("BASH_COMP_PATH");
      return;
    }

    void ArgHandler::command_to_bash(std::vector<std::string> &arguments) {
      arguments.push_back("BASH_COMP_BIN");
      return;
    }

    void scelname_to_bash(std::vector<std::string> &arguments) {
      return;
    }

    void ArgHandler::query_to_bash(std::vector<std::string> &arguments) {
      DataFormatterDictionary<Configuration> dict = make_dictionary<Configuration>();

      for(auto it = dict.begin(); it != dict.cend(); ++it) {
        if(it->type() ==  BaseDatumFormatter<Configuration>::Property) {
          arguments.push_back(it->name());
        }
      }
      return;
    }

    void ArgHandler::operator_to_bash(std::vector<std::string> &arguments) {
      DataFormatterDictionary<Configuration> dict = make_dictionary<Configuration>();

      for(auto it = dict.begin(); it != dict.cend(); ++it) {
        if(it->type() ==  BaseDatumFormatter<Configuration>::Property) {
          arguments.push_back(it->name());
        }
      }
      return;
    }


    /**
     * This construction right here determines what the value_name of the boost options
     * should be named. It is through these strings that bash completion can
     * know which types of completions to suggest.
     */

    const std::vector<std::pair<std::string, ARG_TYPE> > ArgHandler::m_argument_table({
      std::make_pair("<path>", ARG_TYPE::PATH),
      std::make_pair("<command>", ARG_TYPE::COMMAND),
      std::make_pair("<supercell>", ARG_TYPE::SCELNAME),
      std::make_pair("<query>", ARG_TYPE::QUERY),
      std::make_pair("<operation>", ARG_TYPE::OPERATOR),
      std::make_pair("<configuration>", ARG_TYPE::CONFIGNAME),
      std::make_pair("<type>", ARG_TYPE::COORDTYPE)
    });


    //*****************************************************************************************************//

    OptionHandlerBase::OptionHandlerBase(const std::string &init_option_tag):
      m_tag(init_option_tag),
      m_desc(std::string("'casm ") + init_option_tag + std::string("' usage")),
      m_gzip_flag(false) {
    }

    OptionHandlerBase::OptionHandlerBase(const std::string &init_option_tag, const std::string &init_descriptor):
      m_tag(init_option_tag),
      m_desc(init_descriptor),
      m_gzip_flag(false) {
    }

    const std::string &OptionHandlerBase::tag() const {
      return m_tag;
    }

    ///Get the variables map
    po::variables_map &OptionHandlerBase::vm() {
      return m_vm;
    }

    ///Get the variables map
    const po::variables_map &OptionHandlerBase::vm() const {
      return m_vm;
    }

    /**
     * Check if there are any program options in the options description. If there aren't, then
     * this is the first time someone is asking for those values, which we set through the
     * initialize routine. If there are values there already, just hand them back.
     */

    const po::options_description &OptionHandlerBase::desc() {
      if(m_desc.options().size() == 0) {
        initialize();
      }

      return m_desc;
    }

    const fs::path &OptionHandlerBase::selection_path() const {
      return m_selection_path;
    }

    const std::vector<fs::path> &OptionHandlerBase::selection_paths() const {
      return m_selection_paths;
    }

    const std::string &OptionHandlerBase::verbosity_str() const {
      return m_verbosity_str;
    }

    const fs::path OptionHandlerBase::settings_path() const {
      return m_settings_path;
    }

    std::string OptionHandlerBase::input_str() const {
      return m_input_str;
    }

    const fs::path OptionHandlerBase::output_path() const {
      return m_output_path;
    }

    bool OptionHandlerBase::gzip_flag() const {
      return m_gzip_flag;
    }

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
      COORD_TYPE selected_mode;

      if(m_coordtype_str[0] == 'F' || m_coordtype_str[0] == 'f') {
        selected_mode = COORD_TYPE::FRAC;
      }

      else if(m_coordtype_str[0] == 'C' || m_coordtype_str[0] == 'c') {
        selected_mode = COORD_TYPE::CART;
      }

      else {
        selected_mode = COORD_TYPE::COORD_DEFAULT;
      }

      return selected_mode;
    }

    void OptionHandlerBase::add_configlist_suboption(const fs::path &_default) {
      m_desc.add_options()
      ("config,c",
       po::value<fs::path>(&m_selection_path)->default_value(_default)->value_name(ArgHandler::path()),
       (std::string("Only consider the selected configurations of the given selection file. "
                    "Standard selections are 'MASTER', 'CALCULATED', 'ALL', or 'NONE'. "
                    "If not specified, '") + _default.string() + std::string("' will be used.")).c_str());
      return;
    }

    void OptionHandlerBase::add_configlists_suboption(const fs::path &_default) {
      m_desc.add_options()
      ("configs,c",
       po::value<std::vector<fs::path> >(&m_selection_paths)->multitoken()->default_value(std::vector<fs::path> {_default})->value_name(ArgHandler::path()),
       (std::string("Only consider the selected configurations of the given selection file. "
                    "Standard selections are 'MASTER', 'CALCULATED', 'ALL', or 'NONE'. "
                    "If not specified, '") + _default.string() + std::string("' will be used.")).c_str());
      return;
    }

    void OptionHandlerBase::add_configlist_nodefault_suboption() {
      m_desc.add_options()
      ("config,c",
       po::value<fs::path>(&m_selection_path)->value_name(ArgHandler::path()),
       "Only consider the selected configurations of the given selection file. "
       "Standard selections are 'MASTER', 'CALCULATED', 'ALL', or 'NONE'. ");
      return;
    }

    void OptionHandlerBase::add_configlists_nodefault_suboption() {
      m_desc.add_options()
      ("configs,c",
       po::value<std::vector<fs::path> >(&m_selection_paths)->multitoken()->value_name(ArgHandler::path()),
       "Only consider the selected configurations of the given selection files. "
       "Standard selections are 'MASTER', 'CALCULATED', 'ALL', or 'NONE'. ");
      return;
    }

    void OptionHandlerBase::add_help_suboption() {
      m_desc.add_options()
      ("help,h", "Print help message")
      ("desc", "Print extended usage description");
      return;
    }

    void OptionHandlerBase::add_general_help_suboption() {
      m_desc.add_options()
      ("help,h", po::value<std::vector<std::string> >(&m_help_opt_vec)->multitoken()->zero_tokens(), "Print general help. Use '--help properties' for a list of query-able properties or '--help operators' for a list of query operators");
      return;
    }

    void OptionHandlerBase::add_verbosity_suboption() {
      //TODO: add ArgHandler for this
      m_desc.add_options()
      ("verbosity", po::value<std::string>(&m_verbosity_str)->default_value("standard"), "Verbosity of output. Options are 'none', 'quiet', 'standard', 'verbose', 'debug', or an integer 0-100 (0: none, 100: all).");
      return;
    }

    void OptionHandlerBase::add_settings_suboption(bool required) {
      std::string help_str = "Settings input file specifying which parameters should be used. See 'casm format --" + m_tag + "'.";

      if(required) {
        m_desc.add_options()
        ("settings,s", po::value<fs::path>(&m_settings_path)->required()->value_name(ArgHandler::path()), help_str.c_str());
      }
      else {
        m_desc.add_options()
        ("settings,s", po::value<fs::path>(&m_settings_path)->value_name(ArgHandler::path()), help_str.c_str());
      }

      return;
    }

    void OptionHandlerBase::add_input_suboption(bool required) {
      std::string help_str = "String specifying input settings. See 'casm format --" + m_tag + "'.";

      if(required) {
        m_desc.add_options()
        ("input,i", po::value<std::string>(&m_input_str)->required(), help_str.c_str());
      }
      else {
        m_desc.add_options()
        ("input,i", po::value<std::string>(&m_input_str), help_str.c_str());
      }

      return;
    }

    void OptionHandlerBase::add_output_suboption() {
      m_desc.add_options()
      ("output,o", po::value<fs::path>(&m_output_path)->value_name(ArgHandler::path()), "Name for output file. Use STDOUT to print results without extra messages.");
      return;
    }

    void OptionHandlerBase::add_gzip_suboption() {
      m_desc.add_options()
      ("gzip,z", po::value(&m_gzip_flag)->zero_tokens(), "Write gzipped output file.");
      return;
    }

    void OptionHandlerBase::add_scelname_suboption() {
      std::string help_str;
      help_str = "Single supercell name to use casm " + m_tag + " with, such as 'SCEL4_2_2_1_0_0_0'";
      m_desc.add_options()
      ("scelname", po::value<std::string>(&m_supercell_str)->value_name(ArgHandler::supercell()), help_str.c_str());
      return;
    }

    void OptionHandlerBase::add_scelnames_suboption() {
      std::string help_str;
      help_str = "One or more supercells to use casm " + m_tag + " with, such as 'SCEL4_2_2_1_0_0_0'";
      m_desc.add_options()
      ("scelnames", po::value<std::vector<std::string> >(&m_supercell_strs)->multitoken()->value_name(ArgHandler::supercell()), help_str.c_str());
      return;
    }

    void OptionHandlerBase::add_configname_suboption() {
      std::string help_str;

      help_str = "Single configuration to use casm " + m_tag + " with, such as 'SCEL4_2_2_1_0_0_0/3'";

      m_desc.add_options()
      ("configname", po::value<std::string>(&m_config_str)->value_name(ArgHandler::configname()), help_str.c_str());

      return;
    }

    void OptionHandlerBase::add_confignames_suboption() {
      std::string help_str;

      help_str = "One or more configurations to use casm " + m_tag + " with, such as 'SCEL4_2_2_1_0_0_0/3'";

      m_desc.add_options()
      ("confignames", po::value<std::vector<std::string> >(&m_config_strs)->multitoken()->value_name(ArgHandler::configname()), help_str.c_str());

      return;
    }


    void OptionHandlerBase::add_coordtype_suboption() {
      m_desc.add_options()
      ("coord", po::value<std::string>(&m_coordtype_str)->default_value("frac")->value_name(ArgHandler::coordtype()), "Type of coordinate system to use. Use 'frac' for fractional (default) or 'cart' for Cartesian.");
      return;
    }

    //*****************************************************************************************************//

  }
}

#endif
