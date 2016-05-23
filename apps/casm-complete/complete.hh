#include <string>
#include <vector>

#include "casm/CASM_global_definitions.hh"

namespace Completer {
  enum class ARG_TYPE {VOID, PATH, COMMAND, SCELNAME, QUERY, OPERATOR};

  class ArgumentBase {
  public:

  private:
  };

  class SuboptionBase {
  public:

  private:
  };

  class Option {
  public:

  private:
  };

  /**
   * A typical casm command takes 3 types of arguments:
   *  Option: Things like "super", "enum", "monte", etc
   *  Suboption: Anything that starts with "--" or "-", such as "--max" or "-z"
   *  Argument: The argument you pass to the suboption, such as an integer or filename
   *
   * The Engine class contains all the information you need to
   * pass to a bash completion script.
   */

  class Engine {
    typedef std::size_t size_type;
  public:

    ///Default constructor so you can push back your own things
    Engine() {};

    ///Construct by passing list of options
    Engine(const std::vector<Option> &init_options);

    ///Construct by passing boost program options (eventually preferred so that it can update itself)
    Engine(const std::vector<CASM::po::options_description> &init_premade_descs);

    ///Get a list of all available options
    std::vector<std::string> probe_options() const;

    ///For a particular option, get the available suboptions
    std::vector<std::string> probe_suboptions(const std::string &option) const;

    ///For a particular option with suboption, get what kind of arguments are expected
    ARG_TYPE probe_argument_type(const std::string &option, const std::string &suboption) const;

    ///Return the arguments that should be bash completed
    std::vector<std::string> probe_arguments(ARG_TYPE argument_type) const;

    ///Guess what should be returned based on the current word (probably not gonna make it that smart)
    std::vector<std::string> probe_gess(const std::string current_word) const;

    ///Append a new option to the engine
    void push_back(const Option &new_option);

  private:

    std::vector<Option> m_avail_options;
  };
}
