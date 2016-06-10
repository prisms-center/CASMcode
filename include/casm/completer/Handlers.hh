#ifndef HANDLERS_HH
#define HANDLERS_HH

#include <string>
#include <vector>
#include <boost/program_options.hpp>
#include "casm/CASM_global_definitions.hh"

namespace CASM {
  namespace Completer {
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

      enum ARG_TYPE {VOID, PATH, COMMAND, SCELNAME, QUERY, OPERATOR, CONFIGNAME, COORDMODE};

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

      ///Get value_type string for configuration completion
      static std::string configname();

      ///Get value_type string for coordinate mode completion
      static std::string coordmode();


      ///Fill the output strings with bash completion appropriate values for VOID (i.e. do nothing)
      static void void_to_bash(std::vector<std::string> &arguments);

      ///Fill the output strings with bash completion appropriate values for PATH
      static void path_to_bash(std::vector<std::string> &arguments);

      ///Fill the output strings with bash completion appropriate values for COMMAND
      static void command_to_bash(std::vector<std::string> &arguments);

      ///Fill the output strings with bash completion appropriate values for SCELNAME   TODO: This routine is currently unimplemented
      static void scelname_to_bash(std::vector<std::string> &arguments);

      ///Fill the output strings with bash completion appropriate values for QUERY
      static void query_to_bash(std::vector<std::string> &arguments);

      ///Fill the output strings with bash completion appropriate values for OPERATOR
      static void operator_to_bash(std::vector<std::string> &arguments);


    private:

      ///List of pairs relating the value type name of po::option_description to its corresponding argument type
      static const std::vector<std::pair<std::string, ARG_TYPE> > m_argument_table;

    };

    /**
     * This is the base class to be used for every casm option. The purpose of this is to
     * provide a way to generate standard suboptions, such as requesting a supercell name,
     * specifying a coordinate mode, or asking for an output file.
     * Derived classes will hold all the required variables for po::options_description construction
     * that get passed by referenced and parsed from command line input.
     */

    class OptionHandlerBase {
    public:

      ///Define the name of the command during construction
      OptionHandlerBase(const std::string &init_option_tag);

      ///Get the program options, filled with the initialized values
      const po::options_description &desc();

      ///The desired name for the casm option (Perhaps this should get tied with CommandArg?)
      const std::string &tag() const;


    protected:

      ///name of the casm command
      std::string m_tag;

      ///Boost program options. All the derived classes have them, but will fill them up themselves
      po::options_description m_desc;

      ///Fill in the options descriptions accordingly
      virtual void initialize() = 0;

      //-------------------------------------------------------------------------------------//

      ///Add --config suboption (defaults to MASTER)
      void add_configlist_suboption();

      ///The selection string to go with add_config_suboption
      std::string m_selection_str;

      ///Returns the string corresponding to add_config_suboption()
      const std::string &selection_str() const;

      //-------------------------------------------------------------------------------------//

      ///Add a plain --help suboption
      void add_help_suboption();

      //-------------------------------------------------------------------------------------//

      ///Add a smart --help suboption that takes "properties" or "operators"
      void add_general_help_suboption();

      ///The list of strings to go with add_general_help_suboption()
      std::vector<std::string> m_help_opt_vec;

      ///Returns the list of strings corresponding to add_general_help_suboption()
      const std::vector<std::string> &help_opt_vec() const;

      //-------------------------------------------------------------------------------------//

      ///Add a --verbosity suboption. Default "standard" of "none", "quiet", "standard", "verbose", "debug" or an int 0-100
      void add_verbosity_suboption();

      ///The verbosity string to go with add_config_suboption
      std::string m_verbosity_str;

      ///Returns the string corresponding to add_verbosity_suboption()
      const std::string &verbosity_str() const;

      //-------------------------------------------------------------------------------------//

      ///Add a --settings suboption. Expects a corresponding `casm format` to go with it.
      void add_settings_suboption();

      ///The settings path to go with add_settings_suboption()
      fs::path m_settings_path;

      ///Returns the path corresponding to add_settings_suboption
      const fs::path settings_path() const;

      //-------------------------------------------------------------------------------------//

      ///Add a --output suboption. Expects to allow "STDOUT" to print to screen.
      void add_output_suboption();

      ///The path that goes with add_output_suboption
      fs::path m_output_path;

      ///Returns the path corresponding to add_output_suboption()
      const fs::path output_path() const;

      //-------------------------------------------------------------------------------------//

      ///Add a --gzip suboption. The value will default to false unless overridden by the derived class.
      void add_gzip_suboption();

      ///The bool that goes with add_gzip_suboption
      bool m_gzip_flag;

      ///Returns the value assigned for add_gzip_suboption()
      bool gzip_flag() const;

      //-------------------------------------------------------------------------------------//

      ///Add a --scelname suboption.
      void add_scelname_suboption();

      ///The string of the supercell name of add_scelname_suboption()
      std::string m_supercell_str;

      ///Returns the name of the supercell for add_scelname_suboption()
      const std::string &supercell_str() const;

      //----------------------------//

      ///Add a --scelnames suboption.
      void add_scelnames_suboption();

      ///The list of supercell names of add_scelnames_suboption()
      std::vector<std::string> m_supercell_vec;

      ///Returns the list of the supercells for add_scelnames_suboption()
      const std::vector<std::string> &supercell_vec() const;

      //-------------------------------------------------------------------------------------//

      ///Add a --configname suboption.
      void add_configname_suboption();

      ///The name of a single configname to go with add_configname_suboption()
      std::string m_config_str;

      ///Returns the name of the supercell for add_configname_suboption(), for when multiple=false
      const std::string &config_str() const;

      //----------------------------//

      ///Add a --confignames suboption.
      void add_confignames_suboption();

      ///The list of the supercell names of add_configname_suboption()
      std::vector<std::string> m_config_vec;

      ///Returns the names of the supercells for add_configname_suboption(), for when multiple=false
      const std::vector<std::string> &config_vec() const;

      //-------------------------------------------------------------------------------------//

      //Add more general suboption adding routines here//
    };

    /**
     * Options set for `casm monte`. Get your Monte Carlo completion here.
     */

    class MonteOption : public OptionHandlerBase {


    public:

      using OptionHandlerBase::verbosity_str;
      using OptionHandlerBase::settings_path;

      MonteOption();

      Index condition_index() const;

    private:

      void initialize() override;

      Index m_condition_index;
    };

    //*****************************************************************************************************//

    /**
     * Options set for `casm run`. Get casm to run your executables here.
     */

    class RunOption : public OptionHandlerBase {

    public:

      using OptionHandlerBase::selection_str;

      RunOption();

      const std::string &exec_str() const;

    private:

      void initialize() override;

      std::string m_exec_str;
    };

    //*****************************************************************************************************//

    /**
     * Options set for `casm query`. Get casm to query your things here.
     */

    class QueryOption : public OptionHandlerBase {

    public:

      using OptionHandlerBase::selection_str;
      using OptionHandlerBase::output_path;
      using OptionHandlerBase::gzip_flag;
      using OptionHandlerBase::help_opt_vec;

      QueryOption();

      const std::vector<std::string> &columns_vec() const;

      const std::vector<std::string> &new_alias_vec() const;

      bool json_flag() const;

      bool verbatim_flag() const;

      bool no_header_flag() const;

    private:

      void initialize() override;

      std::vector<std::string> m_columns_vec;

      std::vector<std::string> m_new_alias_vec;

      bool m_json_flag;

      bool m_verbatim_flag;

      bool m_no_header_flag;

    };

    //*****************************************************************************************************//

    /**
     * Options set for `casm bset`. Get your clusters here.
     */

    class BsetOption : public OptionHandlerBase {

    public:

      BsetOption();

    private:

      void initialize() override;

    };

    //*****************************************************************************************************//

    /**
     * Options set for `casm bset`. Get your clusters here.
     */

    class CompositionOption : public OptionHandlerBase {

    public:

      CompositionOption();

      const std::string &axis_choice_str() const;

    private:

      void initialize() override;

      std::string m_axis_choice_str;

    };

    //*****************************************************************************************************//

    /**
     * Options set for `casm ref`. Get your reference set here.
     */

    class RefOption : public OptionHandlerBase {

    public:

      using OptionHandlerBase::supercell_str;
      using OptionHandlerBase::config_str;

      RefOption();

      const std::string &set_str() const;

    private:

      void initialize() override;

      std::string m_set_str;

    };

    //*****************************************************************************************************//

    /**
     * Options set for `casm files`. Get your casm project packed up here.
     */

    class FilesOption : public OptionHandlerBase {

    public:

      using OptionHandlerBase::output_path;
      using OptionHandlerBase::gzip_flag;

      FilesOption();

      const std::vector<std::string> &calc_vec() const;

      const std::string &settings_str() const;

    private:

      void initialize() override;

      std::vector<std::string> m_calc_vec;

      std::string m_settings_str;

    };

    //*****************************************************************************************************//

    /**
     * Options set for `casm format`. Get your input files here.
     */

    class FormatOption : public OptionHandlerBase {

    public:

      FormatOption();

    private:

      void initialize() override;

    };

    //*****************************************************************************************************//

    /**
     * Options set for `casm import`. Add new structures to your project here.
     */

    class ImportOption : public OptionHandlerBase {

    public:

      ImportOption();

      double vol_tolerance() const;

      double lattice_weight() const;

      const std::vector<fs::path> &pos_vec() const;

      const fs::path &batch_path() const;


    private:

      void initialize() override;

      double m_vol_tolerance;

      double m_lattice_weight;

      std::vector<fs::path> m_pos_vec;

      fs::path m_batch_path;

    };

    //*****************************************************************************************************//

    /**
     * Options set for `casm format`. Get your input files here.
     */

    class InitOption : public OptionHandlerBase {

    public:

      InitOption();

    private:

      void initialize() override;

    };

    //*****************************************************************************************************//

  }
}


#endif
