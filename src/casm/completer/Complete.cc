#ifndef COMPLETE_CC
#define COMPLETE_CC

#include "casm/completer/Complete.hh"
#include "casm/casm_io/DataFormatter.hh"
#include "casm/clex/Configuration.hh"

namespace CASM {
  namespace Completer {
    typedef ArgHandler::ARG_TYPE ARG_TYPE;

    std::string strip_argument(const std::string &raw_input) {
      std::string stripped = raw_input;

      if(stripped[0] == '-') {
        stripped.erase(stripped.begin() + 1, stripped.end());
      }

      if(stripped[0] == '-') {
        stripped.erase(stripped.begin() + 1, stripped.end());
      }

      return stripped;
    }

    //*****************************************************************************************************//

    namespace Suboption_impl {
      std::string pull_short(const po::option_description &single_boost_option) {
        std::string possible_short;
        possible_short = single_boost_option.canonical_display_name(po::command_line_style::allow_dash_for_short);
        if(possible_short.size() > 2 || possible_short[0] != '-' || possible_short[1] == '-') {
          return "-_";
        }
        else {
          return possible_short;
        }
      }

      std::string pull_long(const po::option_description &single_boost_option) {
        return single_boost_option.canonical_display_name(po::command_line_style::allow_long);
      }
    }

    Suboption::Suboption(const std::string &init_longname, std::string init_short, ARG_TYPE init_expected_types):
      m_long(init_longname),
      m_short(init_short),
      m_expected_arg(init_expected_types) {
      _sanity_throw();
    }


    Suboption::Suboption(const po::option_description &init_boost_option):
      m_long(Suboption_impl::pull_long(init_boost_option)),
      m_short(Suboption_impl::pull_short(init_boost_option)),
      m_expected_arg(ArgHandler::determine_type(init_boost_option)) {
      _sanity_throw();
    }

    void Suboption::_sanity_throw() const {
      if(m_long.size() < 3 || m_short.size() != 2) {
        throw std::runtime_error("--long option must be at least 3 characters long and -s(hort) must be exactly 2!");
      }
      if(m_long[0] != '-' || m_long[1] != '-' || m_short[0] != '-') {
        throw std::runtime_error("Suboption --long and -s(hort) tags must include leading '-' characters!");
      }
    }

    std::string Suboption::long_tag() const {
      return m_long;
    }

    std::string Suboption::short_tag() const {
      return m_short;
    }

    bool Suboption::matches(const std::string &test_tag) const {
      if(test_tag.size() == 2) {
        return (test_tag == m_short);
      }

      else {
        return (test_tag == m_long);
      }

      throw std::runtime_error("The impossible has occurred. Stepped through an if/else block and hit neither");
    }

    ARG_TYPE Suboption::argument_type() const {
      return m_expected_arg;
    }

    //***************************************************************************************************//

    Option::Option(const std::string &init_tag, const std::vector<Suboption> &init_allowed_subopts):
      m_tag(init_tag),
      m_avail_suboptions(init_allowed_subopts) {
    }

    Option::Option(const std::string &init_tag, const po::options_description &init_premade_descs):
      m_tag(init_tag) {
      m_avail_suboptions.reserve(init_premade_descs.options().size());

      for(auto it = init_premade_descs.options().begin(); it != init_premade_descs.options().end(); ++it) {
        m_avail_suboptions.push_back(Suboption(**it));
      }
    }

    std::string Option::tag() const {
      return m_tag;
    }

    /**
     * Runs through all the stored Suboptions and creates a list of each --long
     * tag they have. The list of strings is then returned.
     */

    std::vector<std::string> Option::probe_suboptions() const {
      std::vector<std::string> suboption_tag_list;
      for(auto it = m_avail_suboptions.begin(); it != m_avail_suboptions.end(); ++it) {
        suboption_tag_list.push_back(it->long_tag());
      }

      return suboption_tag_list;
    }

    /**
     * Runs through all the available Suboptions trying to match the provided suboption_tag.
     * If the suboption_tag matches, then return the expected type of arguments of that suboption,
     * if the provided tag doesn't match anything, then default to returning VOID, i.e.
     * no completion should be done by bash.
     */

    ARG_TYPE Option::probe_argument_type(const std::string &suboption_tag) const {
      for(auto it = m_avail_suboptions.begin(); it != m_avail_suboptions.end(); ++it) {
        if(it->matches(suboption_tag)) {
          return it->argument_type();
        }
      }

      return ARG_TYPE::VOID;
    }

    bool Option::matches(const std::string &test_tag) const {
      return (test_tag == m_tag);
    }

    //***************************************************************************************************//

    Engine::Engine(const std::vector<Option> &init_options):
      m_avail_options(init_options) {
    }

    /**
     * Runs through all the available Options and generates a list of all their tags.
     * This list is then returned.
     */

    std::vector<std::string> Engine::probe_options() const {
      std::vector<std::string> option_list;
      for(auto it = m_avail_options.begin(); it != m_avail_options.end(); ++it) {
        option_list.push_back(it->tag());
      }

      return option_list;
    }

    /**
     * Run through all the Options and try to match the provided option_tag. If one of the Options
     * has that tag, then return all the Suboptions it has as a list of --long formatted strings.
     * If the proved tag doesn't match any options, return an empty list.
     */

    std::vector<std::string> Engine::probe_suboptions(const std::string &option_tag) const {
      for(auto it = m_avail_options.begin(); it != m_avail_options.end(); ++it) {
        if(it->matches(option_tag)) {
          return it->probe_suboptions();
        }
      }

      // as fallback, show files/dirs
      return std::vector<std::string>({"BASH_COMP_PATH"});
    }

    /**
     * Determine what kind of argument should follow the provided suboption and return
     * either a bash completion suggestion, or a keyword the bash completion script
     * can identify to make the completion itself (e.g. completing paths or available executables)
     */

    std::vector<std::string> Engine::probe_arguments(const std::string &option_tag, const std::string &suboption_tag) const {
      std::vector<std::string> arguments;

      ARG_TYPE required_arg = _probe_argument_type(option_tag, suboption_tag);

      switch(required_arg) {
      case ARG_TYPE::VOID:
        ArgHandler::void_to_bash(arguments);
        break;

      case ARG_TYPE::PATH:
        ArgHandler::path_to_bash(arguments);
        break;

      case ARG_TYPE::COMMAND:
        ArgHandler::command_to_bash(arguments);
        break;

      case ARG_TYPE::SCELNAME:
        //Add supercell names from PrimClex pointer somehow
        break;

      case ARG_TYPE::QUERY:
        ArgHandler::query_to_bash(arguments);
        break;

      case ARG_TYPE::OPERATOR:
        ArgHandler::operator_to_bash(arguments);
        break;

      case ARG_TYPE::CONFIGNAME:
        //ArgHandler::configname_to_bash();
        break;

      case ARG_TYPE::COORDTYPE:
        //ArgHandler::coord_mode_to_bash();
        break;

      default:
        break;
      }

      // as fallback, show files/dirs
      if(!arguments.size()) {
        arguments.push_back("BASH_COMP_PATH");
      }

      return arguments;
    }



    /**
     * First run through all the Options and see if you can find one that matches the provided
     * option_tag. If one is found, then run through all the available Suboptions of said tag
     * and see if the provided Suboption_tag can be matched. If it can be matched, return the
     * expected argument type of that Suboption.
     *
     * If either tag can't be found the routine returns VOID, i.e. bash completion shouldn't make
     * any suggestions.
     */

    ARG_TYPE Engine::_probe_argument_type(const std::string &option_tag, const std::string &suboption_tag) const {
      for(auto it = m_avail_options.begin(); it != m_avail_options.end(); ++it) {
        if(it->matches(option_tag)) {
          return it->probe_argument_type(suboption_tag);
        }
      }

      return ARG_TYPE::VOID;
    }

    void Engine::push_back(const Option &new_option) {
      m_avail_options.push_back(new_option);
      return;
    }
  }
}

#endif
