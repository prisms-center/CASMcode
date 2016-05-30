#ifndef COMPLETE_HH
#define COMPLETE_HH

#include <string>
#include <vector>
#include <stdexcept>
#include <utility>
#include <boost/program_options.hpp>
#include "casm/CASM_global_definitions.hh"

namespace CASM {

  class PrimClex;

  namespace Completer {

    ///Remove "--" or "-" from beginning of string if it exists, and return as new string
    std::string strip_argument(const std::string &raw_input);

    /**
     * Handle the value type names of po::option_description. This class
     * determines what the keywords mean, and translates them into the
     * ARG_TYPE enum as appropriate.
     *
     * If you want bash completion for your boost program options, never specify
     * option_description::value_name manually (e.g. raw string), always request
     * the string through this class.
     *
     * When the Engine class isn't returning strings corresponding to options
     * or suboptions that bash should complete, then the user is in the
     * process of writing out an argument. By printing out keywords, the
     * bash completer script can know what kind of arguments are needed for the
     * casm command:
     *
     * VOID:        Don't complete anything
     * PATH:        Suggest path completions to a file or directory
     * COMMAND:     Suggest executables within the environment PATH
     * SCELNAME:    Run through the PrimClex and suggest the enumerated supercell names
     * QUERY:       Suggest the available query options
     * OPERATORS:   Suggest the available operators casm knows how to use (TODO: This one might be tricky to implement)
     */

    class ArgHandler {
    public:

      enum ARG_TYPE {VOID, PATH, COMMAND, SCELNAME, QUERY, OPERATOR};

      ///Translate the stored boost value_name into an ARG_TYPE for the completer engine
      static ARG_TYPE determine_type(const po::option_description &boost_option);

      ///Get value_type string for path completion
      static std::string path();

      ///Get value_type string for command completion (i.e. stuff in your $PATH)
      static std::string command();

      ///Get value_type string for supercell completion
      static std::string supercell();

      ///Get value_type string for query completion
      static std::string query();

      ///Get value_type string for operation completion
      static std::string operation();


      ///Fill the output strings with bash completion appropriate values for VOID (i.e. do nothing)
      static void void_to_bash(std::vector<std::string> &arguments);

      ///Fill the output strings with bash completion appropriate values for PATH
      static void path_to_bash(std::vector<std::string> &arguments);

      ///Fill the output strings with bash completion appropriate values for COMMAND
      static void command_to_bash(std::vector<std::string> &arguments);

      ///Fill the output strings with bash completion appropriate values for SCELNAME
      //static void scelname_to_bash(std::vector<std::string> &arguments);

      ///Fill the output strings with bash completion appropriate values for QUERY
      static void query_to_bash(std::vector<std::string> &arguments);

      ///Fill the output strings with bash completion appropriate values for OPERATOR
      static void operator_to_bash(std::vector<std::string> &arguments);


    private:

      ///List of pairs relating the value type name of po::option_description to its corresponding argument type
      static const std::vector<std::pair<std::string, ARG_TYPE> > m_argument_table;

    };

    //*****************************************************************************************************//

    /**
     * A typical casm command takes 3 types of parameters:
     *  Option: Things like "super", "enum", "monte", etc
     *  Suboption: Anything that starts with "--" or "-", such as "--max" or "-z"
     *  Argument: The argument you pass to the suboption, such as an integer or filename
     *
     * The Suboption class consists of a string to identify the --long option, a character
     * to identify the -s(hort) version, and the type of arguments that are expected.
     */

    class Suboption {

    public:
      typedef ArgHandler::ARG_TYPE ARG_TYPE;

      ///Explicit construction. Be sure to include "--" and '-' in the tags
      Suboption(const std::string &init_longname, std::string init_short, ARG_TYPE init_expected_types);

      ///Construct with boost objects
      Suboption(const po::option_description &init_boost_option);

      ///Return long name in string format
      std::string long_tag() const;

      ///Return short name as char
      std::string short_tag() const;

      ///See if a provided string matches either the --long or -s(hort) tags. Expects leading '-' characters.
      bool matches(const std::string &test_tag) const;

      ///Return the expected types of arguments that follow *this
      ARG_TYPE argument_type() const;

    private:

      ///--long identifier (includes leading "--")
      const std::string m_long;

      ///-s(hort) identifyer (includes leading '-')
      const std::string m_short;

      ///Type of arguments expected
      const ARG_TYPE m_expected_arg;

      ///Make sure values make sense
      bool _sanity_throw() const;
    };

    namespace Suboption_impl {
      ///Get the -s(hort) tag from boost, or make it "- " if it doesn't exist
      std::string pull_short(const po::option_description &single_boost_option);

      ///Get the --long tag from boost
      std::string pull_long(const po::option_description &single_boost_option);
    }

    //*****************************************************************************************************//


    /**
     * A typical casm command takes 3 types of arguments:
     *  Option: Things like "super", "enum", "monte", etc
     *  Suboption: Anything that starts with "--" or "-", such as "--max" or "-z"
     *  Argument: The argument you pass to the suboption, such as an integer or filename
     *
     * The Option class has a string that identifies it (e.g. "enum"), as well as a list
     * of Suboptions that it can take.
     */


    class Option {
      typedef std::size_t size_type;
      typedef ArgHandler::ARG_TYPE ARG_TYPE;

    public:

      ///Construct manually passing all values needed
      Option(const std::string &init_tag, const std::vector<Suboption> &init_allowed_subopts);

      ///Construct with program options (eventually preferred)
      Option(const std::string &init_tag, const po::options_description &init_premade_descs);

      ///Return the identifying name of *this (e.g. "super", "monte", etc)
      std::string tag() const;

      ///Return what the suboptions (--long format) for *this are
      std::vector<std::string> probe_suboptions() const;

      ///For a particular --suboption, get what kind of arguments are expected. suboption_tag should be pre-stripped.
      ARG_TYPE probe_argument_type(const std::string &suboption_tag) const;

      ///Append a new suboption to the option
      void push_back(const Suboption &new_suboption);

      ///Check if the given string corresponds to the tag of *this
      bool matches(const std::string &test_tag) const;

    private:

      ///Name that identifies this casm option (e.g. "monte", "init", etc)
      std::string m_tag;  //Make this const once the default constructors of Engine are gone?

      ///List of all the available --suboptions this option has
      std::vector<Suboption> m_avail_suboptions;

    };


    //*****************************************************************************************************//


    /**
     * A typical casm command takes 3 types of arguments:
     *  Option: Things like "super", "enum", "monte", etc
     *  Suboption: Anything that starts with "--" or "-", such as "--max" or "-z"
     *  Argument: The argument you pass to the suboption, such as an integer or filename
     *
     * The Engine class contains all the information you need to
     * pass to a bash completion script. However, it is currently limited to some
     * specific kind of formatting. Boost program option can do a bunch of fancy
     * things, but we don't seem to be using them much, so the completer is currently
     * limited to expecting only --long and -s(hort) (single character!) options.
     */

    class Engine {
      typedef std::size_t size_type;
      typedef ArgHandler::ARG_TYPE ARG_TYPE;

    public:

      ///Default constructor so you can push back your own things
      Engine() {};

      ///Construct by passing list of options
      Engine(const std::vector<Option> &init_options);

      ///Construct by passing boost program options (eventually preferred so that it can update itself)
      //Engine(const std::vector<CASM::po::options_description> &init_premade_descs);

      ///Get a list of all available options
      std::vector<std::string> probe_options() const;

      ///For a particular option, get the available suboptions
      std::vector<std::string> probe_suboptions(const std::string &option_tag) const;

      ///Return the arguments that should be bash completed
      std::vector<std::string> probe_arguments(const std::string &option_tag, const std::string &suboption_tag) const;

      ///Guess what should be returned based on the current word (probably not gonna make it that smart)
      //std::vector<std::string> probe_guess(const std::string current_word) const;

      ///Append a new option to the engine
      void push_back(const Option &new_option);

    private:

      ///List of all the available (presumably casm) options (e.g. "monte", "init", etc)
      std::vector<Option> m_avail_options;

      ///For a particular option with suboption, get what kind of arguments are expected
      ARG_TYPE _probe_argument_type(const std::string &option_tag, const std::string &suboption_tag) const;

      ///Pointer to the project PrimClex, so the Engine can know the existing SCEL names
      //const *CASM::PrimClex m_pclex;

    };
  }
}

#endif
