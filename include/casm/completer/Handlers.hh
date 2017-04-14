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

      enum ARG_TYPE {VOID, PATH, COMMAND, SCELNAME, QUERY, OPERATOR, CONFIGNAME, COORDTYPE};

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
      static std::string coordtype();


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

      ///More explicit initialization
      OptionHandlerBase(const std::string &init_option_tag, const std::string &init_descriptor);

      ///Get the variables map
      po::variables_map &vm();

      ///Get the variables map
      const po::variables_map &vm() const;

      ///Get the program options, filled with the initialized values
      const po::options_description &desc();

      ///The desired name for the casm option (Perhaps this should get tied with CommandArg?)
      const std::string &tag() const;


    protected:

      ///name of the casm command
      std::string m_tag;

      ///Boost program options. All the derived classes have them, but will fill them up themselves
      po::options_description m_desc;

      ///Boost program options variable map
      po::variables_map m_vm;

      ///Fill in the options descriptions accordingly
      virtual void initialize() = 0;

      //-------------------------------------------------------------------------------------//

      ///Add --config suboption (defaults to MASTER)
      void add_configlist_suboption(const fs::path &_default = "MASTER");

      ///The selection string to go with add_config_suboption
      fs::path m_selection_path;

      ///Returns the string corresponding to add_config_suboption()
      const fs::path &selection_path() const;

      //----------------------------//

      ///Add --configs suboption (defaults to MASTER)
      void add_configlists_suboption(const fs::path &_default = "MASTER");

      ///The selection string to go with add_config_suboption
      std::vector<fs::path> m_selection_paths;

      ///Returns the string corresponding to add_config_suboption()
      const std::vector<fs::path> &selection_paths() const;

      //-------------------------------------------------------------------------------------//

      ///Add --config suboption (no default)
      void add_configlist_nodefault_suboption();

      //----------------------------//

      ///Add --configs suboption (no default)
      void add_configlists_nodefault_suboption();

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

      ///Add a --input suboption. Expects a corresponding `casm format` to go with it.
      void add_input_suboption(bool required = true);

      ///The settings path to go with add_input_suboption()
      std::string m_input_str;

      ///Returns the path corresponding to add_input_suboption
      std::string input_str() const;

      //-------------------------------------------------------------------------------------//

      ///Add a --settings suboption. Expects a corresponding `casm format` to go with it.
      void add_settings_suboption(bool required = true);

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
      std::vector<std::string> m_supercell_strs;

      ///Returns the list of the supercells for add_scelnames_suboption()
      const std::vector<std::string> &supercell_strs() const;

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
      std::vector<std::string> m_config_strs;

      ///Returns the names of the supercells for add_configname_suboption(), for when multiple=false
      const std::vector<std::string> &config_strs() const;

      //-------------------------------------------------------------------------------------//

      ///Add a --coord suboption to specify FRAC or CART
      void add_coordtype_suboption();

      ///The enum value in the form of a string for add_coordtype_suboption(). Only the first letter matters, but knock yourself out.
      std::string m_coordtype_str;

      ///Return the coordinate type in the form of a string
      const std::string &coordtype_str() const;

      ///Return the coordinate type recasted as the CASM defined enum
      COORD_TYPE coordtype_enum() const;

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

      using OptionHandlerBase::selection_path;

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

      using OptionHandlerBase::selection_path;
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

      double min_va_frac() const;

      double max_va_frac() const;

      const std::vector<fs::path> &pos_vec() const;

      const fs::path &batch_path() const;


    private:

      void initialize() override;

      double m_vol_tolerance;

      double m_lattice_weight;

      double m_min_va_frac;

      double m_max_va_frac;

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

    /**
     * Options set for `casm perturb`. Get your defects here.
     */

    class PerturbOption : public OptionHandlerBase {

    public:

      using OptionHandlerBase::selection_path;

      PerturbOption();

      const fs::path &cspecs_path() const;

    private:

      void initialize() override;

      fs::path m_cspecs_path;

    };

    //*****************************************************************************************************//

    /**
     * Options set for `casm perturb`. Get your defects here.
     */

    class SelectOption : public OptionHandlerBase {

    public:

      using OptionHandlerBase::help_opt_vec;
      using OptionHandlerBase::selection_paths;
      using OptionHandlerBase::output_path;

      SelectOption();

      const std::vector<std::string> &criteria_vec() const;

    private:

      void initialize() override;

      std::vector<std::string> m_criteria_vec;


    };

    //*****************************************************************************************************//

    /**
     * Options set for `casm settings`. Get your casm project configured here.
     */

    class SettingsOption : public OptionHandlerBase {

    public:

      SettingsOption();

      const std::string &input_str() const;

      const std::vector<std::string> &input_vec() const;

    private:

      void initialize() override;

      std::string m_input_str;

      std::vector<std::string> m_input_vec;


    };

    //*****************************************************************************************************//

    /**
     * Options set for `casm status`. Find out about your progress here.
     */

    class StatusOption : public OptionHandlerBase {

    public:

      StatusOption();

    private:

      void initialize() override;

    };

    //*****************************************************************************************************//

    /**
     * Options set for `casm super`. Get your superstructures here.
     */

    class SuperOption : public OptionHandlerBase {

    public:

      using OptionHandlerBase::coordtype_enum;
      using OptionHandlerBase::supercell_strs;
      using OptionHandlerBase::config_strs;
      using OptionHandlerBase::selection_paths;

      SuperOption();

      const std::vector<fs::path> &transf_mat_paths() const;

      const fs::path &struct_path() const;

      const std::string &unit_scel_str() const;

      Index min_vol() const;

      double tolerance() const;

    private:

      void initialize() override;

      std::vector<fs::path> m_transf_mat_paths;

      fs::path m_struct_path;

      std::string m_unit_scel_str;

      Index m_min_vol;

      double m_tolerance;
    };

    //*****************************************************************************************************//

    /**
     * Options set for `casm sym`. Get your point groups here.
     */

    class SymOption : public OptionHandlerBase {

    public:

      using OptionHandlerBase::coordtype_enum;

      SymOption();

    private:

      void initialize() override;

    };

    //*****************************************************************************************************//

    /**
     * Options set for `casm rm`. Get your point groups here.
     */

    class RmOption : public OptionHandlerBase {

    public:

      using OptionHandlerBase::selection_path;
      using OptionHandlerBase::config_strs;
      using OptionHandlerBase::supercell_strs;

      RmOption();

      bool force() const;

      bool dry_run() const;

    private:

      void initialize() override;

    };

    //*****************************************************************************************************//

    /**
     * Options set for `casm update`. Get your project up to date here.
     */

    class UpdateOption : public OptionHandlerBase {

    public:

      using OptionHandlerBase::coordtype_enum;

      using OptionHandlerBase::selection_path;

      UpdateOption();

      double vol_tolerance() const;

      double lattice_weight() const;

      double min_va_frac() const;

      double max_va_frac() const;

    private:

      void initialize() override;

      double m_vol_tolerance;

      double m_lattice_weight;  //TODO: Push to base? Other commands use this

      double m_min_va_frac;

      double m_max_va_frac;

    };

    //*****************************************************************************************************//

    /**
     * Options set for `casm view`. See what your casm landscape looks like here.
     */

    class ViewOption : public OptionHandlerBase {

    public:

      using OptionHandlerBase::config_strs;
      using OptionHandlerBase::selection_path;

      ViewOption();

    private:

      void initialize() override;

    };

    //*****************************************************************************************************//

    /**
     * Options set for `casm enum`. Enumerate configurations and supercells here.
     */

    class EnumOption : public OptionHandlerBase {

    public:

      EnumOption();

      using OptionHandlerBase::settings_path;
      using OptionHandlerBase::input_str;
      using OptionHandlerBase::supercell_strs;

      const std::vector<std::string> &desc_vec() const {
        return m_desc_vec;
      }

      std::string method() const {
        return m_method;
      }

      int min_volume() const {
        return m_min_volume;
      }

      int max_volume() const {
        return m_max_volume;
      }

      bool all_existing() const {
        return m_all_existing;
      }

      const std::vector<std::string> &filter_strs() const {
        return m_filter_strs;
      }

    private:

      void initialize() override;

      std::vector<std::string> m_desc_vec;

      std::string m_method;
      int m_min_volume;
      int m_max_volume;
      bool m_all_existing;
      std::vector<std::string> m_filter_strs;

    };

    //*****************************************************************************************************//

  }
}


#endif
