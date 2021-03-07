#ifndef HANDLERS_HH
#define HANDLERS_HH

#include <boost/filesystem/path.hpp>
#include <boost/program_options.hpp>
#include <string>
#include <vector>

#include "casm/casm_io/container/stream_io.hh"
#include "casm/global/definitions.hh"
#include "casm/global/enum.hh"

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
 * SCELNAME:    Run through the PrimClex and suggest the enumerated supercell
 * names QUERY:       Suggest the available query options OPERATORS:   Suggest
 * the available operators casm knows how to use (TODO: This one might be tricky
 * to implement)
 */

class ArgHandler {
 public:
  enum ARG_TYPE {
    VOID,
    PATH,
    COMMAND,
    SCELNAME,
    QUERY,
    OPERATOR,
    CONFIGNAME,
    COORDTYPE,
    DBTYPE,
    ENUMMETHOD,
    CONFIGTYPE,
    CALCTYPE,
    BSET,
    CLEX,
    REF,
    ECI,
    PROPERTY,
    DOF,
    INFOMETHOD
  };

  /// Translate the stored boost value_name into an ARG_TYPE for the completer
  /// engine
  static ARG_TYPE determine_type(const po::option_description &boost_option);

  /// Get value_type string for path completion
  static std::string path();

  /// Get value_type string for command completion (i.e. stuff in your $PATH)
  static std::string command();

  /// Get value_type string for supercell completion
  static std::string supercell();

  /// Get value_type string for query completion
  static std::string query();

  /// Get value_type string for operation completion
  static std::string operation();

  /// Get value_type string for configuration completion
  static std::string configname();

  /// Get value_type string for coordinate mode completion
  static std::string coordtype();

  /// Get value_type string for dbtype mode completion
  static std::string dbtype();

  /// Get value_type string for enummethod mode completion
  static std::string enummethod();

  /// Get value_type string for configtype mode completion
  static std::string configtype();

  /// Get value_type string for calctype mode completion
  static std::string calctype();

  /// Get value_type string for bset mode completion
  static std::string bset();

  /// Get value_type string for clex mode completion
  static std::string clex();

  /// Get value_type string for ref mode completion
  static std::string ref();

  /// Get value_type string for eci mode completion
  static std::string eci();

  /// Get value_type string for property mode completion
  static std::string property();

  /// Get value_type string for property mode completion
  static std::string dof();

  /// Get value_type string for infomethod mode completion
  static std::string infomethod();

  /// Fill the output strings with bash completion appropriate values for VOID
  /// (i.e. do nothing)
  static void void_to_bash(std::vector<std::string> &arguments);

  /// Fill the output strings with bash completion appropriate values for PATH
  static void path_to_bash(std::vector<std::string> &arguments);

  /// Fill the output strings with bash completion appropriate values for
  /// COMMAND
  static void command_to_bash(std::vector<std::string> &arguments);

  /// Fill the output strings with bash completion appropriate values for
  /// SCELNAME
  static void scelname_to_bash(std::vector<std::string> &arguments);

  /// Fill the output strings with bash completion appropriate values for
  /// CONFIGNAME
  static void configname_to_bash(std::vector<std::string> &arguments);

  /// Fill the output strings with bash completion appropriate values for QUERY
  static void query_to_bash(std::vector<std::string> &arguments);

  /// Fill the output strings with bash completion appropriate values for
  /// OPERATOR
  static void operator_to_bash(std::vector<std::string> &arguments);

  /// Fill the output strings with bash completion appropriate values for DBTYPE
  static void dbtype_to_bash(std::vector<std::string> &arguments);

  /// Fill the output strings with bash completion appropriate values for
  /// ENUMMETHOD
  static void enummethod_to_bash(std::vector<std::string> &arguments);

  /// Fill the output strings with bash completion appropriate values for
  /// CONFIGTYPE
  static void configtype_to_bash(std::vector<std::string> &arguments);

  /// Fill the output strings with bash completion appropriate values for
  /// CALCTYPE
  static void calctype_to_bash(std::vector<std::string> &arguments);

  /// Fill the output strings with bash completion appropriate values for BSET
  static void bset_to_bash(std::vector<std::string> &arguments);

  /// Fill the output strings with bash completion appropriate values for CLEX
  static void clex_to_bash(std::vector<std::string> &arguments);

  /// Fill the output strings with bash completion appropriate values for REF
  static void ref_to_bash(std::vector<std::string> &arguments);

  /// Fill the output strings with bash completion appropriate values for ECI
  static void eci_to_bash(std::vector<std::string> &arguments);

  /// Fill the output strings with bash completion appropriate values for
  /// PROPERTY
  static void property_to_bash(std::vector<std::string> &arguments);

  /// Fill the output strings with bash completion appropriate values for DOF
  static void dof_to_bash(std::vector<std::string> &arguments);

  /// Fill the output strings with bash completion appropriate values for
  /// INFOMETHOD
  static void infomethod_to_bash(std::vector<std::string> &arguments);

 private:
  /// List of pairs relating the value type name of po::option_description to
  /// its corresponding argument type
  static const std::vector<std::pair<std::string, ARG_TYPE> > m_argument_table;
};

/**
 * This is the base class to be used for every casm option. The purpose of this
 * is to provide a way to generate standard suboptions, such as requesting a
 * supercell name, specifying a coordinate mode, or asking for an output file.
 * Derived classes will hold all the required variables for
 * po::options_description construction that get passed by referenced and parsed
 * from command line input.
 *
 * Shortcuts:
 * \code
 * - '-a': --all, --add-canonical
 * - '-b': --batch
 * - '-c': --selection, --selections, --config, --configs
 * - '-d': --display, --details, --data
 * - '-e': --exec
 * - '-f': --force
 * - '-h': --help
 * - '-i': --input, --images
 * - '-j': --json
 * - '-k': --columns
 * - '-l': --list
 * - '-m': --method
 * - '-n': --dry-run, --next
 * - '-o': --output
 * - '-p': --pos
 * - '-P': --path
 * - '-r': --recursive
 * - '-s': --settings
 * - '-t': --type
 * - '-u': --update
 * - '-v': --verbatim
 * - '-w': --warning
 * - '-z': --gzip
 * - '-R': --relative
 *
 * \endcode
 */

class OptionHandlerBase {
 public:
  /// Define the name of the command during construction
  OptionHandlerBase(const std::string &init_option_tag);

  /// More explicit initialization
  OptionHandlerBase(const std::string &init_option_tag,
                    const std::string &init_descriptor);

  /// Get the variables map
  po::variables_map &vm();

  /// Get the variables map
  const po::variables_map &vm() const;

  /// Get the program options, filled with the initialized values
  const po::options_description &desc();

  /// Get the program options, filled with the initialized values
  const po::options_description &desc() const;

  /// The desired name for the casm option (Perhaps this should get tied with
  /// CommandArg?)
  const std::string &tag() const;

 protected:
  /// name of the casm command
  std::string m_tag;

  /// Boost program options. All the derived classes have them, but will fill
  /// them up themselves
  mutable po::options_description m_desc;

  /// Boost program options variable map
  po::variables_map m_vm;

  /// Fill in the options descriptions accordingly
  virtual void initialize() = 0;

  //-------------------------------------------------------------------------------------//

  /// Add --selection suboption (defaults to MASTER)
  void add_selection_suboption(const fs::path &_default = "MASTER");

  /// Add --selection suboption (no default)
  void add_selection_no_default_suboption();

  /// Add --config suboption (defaults to MASTER)
  void add_configlist_suboption(const fs::path &_default = "MASTER");

  /// The selection string to go with add_config_suboption
  fs::path m_selection_path;

  /// Returns the string corresponding to add_config_suboption()
  const fs::path &selection_path() const;

  //----------------------------//

  /// Add --selections suboption (defaults to MASTER)
  void add_selections_suboption(const fs::path &_default = "MASTER");

  /// Add --configs suboption (defaults to MASTER)
  void add_configlists_suboption(const fs::path &_default = "MASTER");

  /// The selection string to go with add_config_suboption
  std::vector<fs::path> m_selection_paths;

  /// Returns the string corresponding to add_config_suboption()
  const std::vector<fs::path> &selection_paths() const;

  //-------------------------------------------------------------------------------------//

  /// Add --prim suboption
  void add_prim_path_suboption(const fs::path &_default = "");

  /// The path string to go with add_prim_path_suboption
  fs::path m_prim_path;

  /// Returns the string corsresponding to add_prim_path_suboption()
  const fs::path &prim_path() const;

  //-------------------------------------------------------------------------------------//

  /// Add --path suboption (defaults to MASTER)
  void add_file_path_suboption(const fs::path &_default = "");

  /// The path string to go with add_file_path_suboption
  fs::path m_file_path;

  /// Returns the string corsresponding to add_file_path_suboption()
  const fs::path &file_path() const;

  //-------------------------------------------------------------------------------------//

  /// Add --config suboption (no default)
  void add_configlist_nodefault_suboption();

  //----------------------------//

  /// Add --configs suboption (no default)
  void add_configlists_nodefault_suboption();

  //-------------------------------------------------------------------------------------//

  /// Add --type suboption (default, set of short_name of allowed ConfigTypes)
  void add_configtype_suboption(std::string _default,
                                std::set<std::string> _configtype_opts);

  /// User-specified config type
  std::string m_configtype;

  /// Set of valid config types
  std::set<std::string> m_configtype_opts;

  std::string configtype() const;

  std::set<std::string> configtype_opts() const;

  //-------------------------------------------------------------------------------------//

  /// Add --type suboption (default, set of short_name of allowed DataObject
  /// Types)
  void add_db_type_suboption(std::string _default,
                             std::set<std::string> _configtype_opts);

  /// User-specified config type
  std::string m_db_type;

  /// Set of valid config types
  std::set<std::string> m_db_type_opts;

  std::string db_type() const;

  std::set<std::string> db_type_opts() const;

  //-------------------------------------------------------------------------------------//

  /// Add a plain --help and --desc suboptions
  void add_help_suboption();

  //-------------------------------------------------------------------------------------//

  /// Add a smart --help suboption that takes "properties" or "operators"
  void add_general_help_suboption();

  /// The list of strings to go with add_general_help_suboption()
  std::vector<std::string> m_help_opt_vec;

  /// Returns the list of strings corresponding to add_general_help_suboption()
  const std::vector<std::string> &help_opt_vec() const;

  //-------------------------------------------------------------------------------------//

  /// Add a --verbosity suboption. Default "standard" of "none", "quiet",
  /// "standard", "verbose", "debug" or an int 0-100
  void add_verbosity_suboption();

  /// The verbosity string to go with add_config_suboption
  std::string m_verbosity_str;

  /// Returns the string corresponding to add_verbosity_suboption()
  const std::string &verbosity_str() const;

  /// Will throw if not expected string or int in range [0, 100]
  int verbosity() const;

  //-------------------------------------------------------------------------------------//

  /// Add a --input suboption. Expects a corresponding `casm format` to go with
  /// it.
  void add_input_suboption(bool required = true);

  /// The settings path to go with add_input_suboption()
  std::string m_input_str;

  /// Returns the path corresponding to add_input_suboption
  std::string input_str() const;

  //-------------------------------------------------------------------------------------//

  /// Add a --settings suboption. Expects a corresponding `casm format` to go
  /// with it.
  void add_settings_suboption(bool required = true);

  /// The settings path to go with add_settings_suboption()
  fs::path m_settings_path;

  /// Returns the path corresponding to add_settings_suboption
  const fs::path settings_path() const;

  //-------------------------------------------------------------------------------------//

  /// Add a --output suboption. Expects to allow "STDOUT" to print to screen.
  void add_output_suboption();

  /// Add a --output suboption, with default value.
  /// Expects to allow "STDOUT" to print to screen.
  void add_output_suboption(const fs::path &_default);

  /// The path that goes with add_output_suboption
  fs::path m_output_path;

  /// Returns the path corresponding to add_output_suboption()
  const fs::path output_path() const;

  //-------------------------------------------------------------------------------------//

  /// Add a --gzip suboption. The value will default to false unless overridden
  /// by the derived class.
  void add_gzip_suboption();

  /// The bool that goes with add_gzip_suboption
  bool m_gzip_flag;

  /// Returns the value assigned for add_gzip_suboption()
  bool gzip_flag() const;

  //-------------------------------------------------------------------------------------//

  /// Add a --scelname suboption.
  void add_scelname_suboption();

  /// The string of the supercell name of add_scelname_suboption()
  std::string m_supercell_str;

  /// Returns the name of the supercell for add_scelname_suboption()
  const std::string &supercell_str() const;

  //----------------------------//

  /// Add a --scelnames suboption.
  void add_scelnames_suboption();

  /// The list of supercell names of add_scelnames_suboption()
  std::vector<std::string> m_supercell_strs;

  /// Returns the list of the supercells for add_scelnames_suboption()
  const std::vector<std::string> &supercell_strs() const;

  //-------------------------------------------------------------------------------------//

  /// Add a --configname suboption.
  void add_configname_suboption();

  /// The name of a single configname to go with add_configname_suboption()
  std::string m_config_str;

  /// Returns the name of the supercell for add_configname_suboption(), for when
  /// multiple=false
  const std::string &config_str() const;

  //----------------------------//

  /// Add a --confignames suboption.
  void add_confignames_suboption();

  /// The list of the supercell names of add_configname_suboption()
  std::vector<std::string> m_config_strs;

  /// Returns the names of the supercells for add_configname_suboption(), for
  /// when multiple=false
  const std::vector<std::string> &config_strs() const;

  //----------------------------//

  /// Add a --names suboption.
  void add_names_suboption();

  /// The list of the supercell names of add_configname_suboption()
  std::vector<std::string> m_name_strs;

  /// Returns the names of the supercells for add_configname_suboption(), for
  /// when multiple=false
  const std::vector<std::string> &name_strs() const;

  //-------------------------------------------------------------------------------------//

  /// Add a --coord suboption to specify FRAC or CART
  void add_coordtype_suboption();

  /// The enum value in the form of a string for add_coordtype_suboption(). Only
  /// the first letter matters, but knock yourself out.
  std::string m_coordtype_str;

  /// Return the coordinate type in the form of a string
  const std::string &coordtype_str() const;

  /// Return the coordinate type recasted as the CASM defined enum
  COORD_TYPE coordtype_enum() const;

  //-------------------------------------------------------------------------------------//

  /// Add a --dofs suboption to specify DoF Types
  void add_dofs_suboption();

  /// The list of DoF type names
  std::vector<std::string> m_dof_strs;

  /// Returns the names of the DoF type names for add_dofs_suboption()
  const std::vector<std::string> &dof_strs() const;

  //-------------------------------------------------------------------------------------//

  void add_dry_run_suboption(std::string msg = default_dry_run_msg());

  static std::string default_dry_run_msg();

  bool dry_run() const;

  //-------------------------------------------------------------------------------------//

  // Add more general suboption adding routines here//
};

/**
 * Options set for `casm monte`. Get your Monte Carlo completion here.
 */

class MonteOption : public OptionHandlerBase {
 public:
  using OptionHandlerBase::settings_path;
  using OptionHandlerBase::verbosity_str;

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
 * Options set for `casm bset`. Get your clusters here.
 */

class BsetOption : public OptionHandlerBase {
 public:
  BsetOption();

  using OptionHandlerBase::coordtype_enum;

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
  using OptionHandlerBase::config_str;
  using OptionHandlerBase::supercell_str;

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
  using OptionHandlerBase::gzip_flag;
  using OptionHandlerBase::output_path;

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
 * Options set for `casm init`. Set up your input here.
 */

class InitOption : public OptionHandlerBase {
 public:
  InitOption();

  using OptionHandlerBase::file_path;

  using OptionHandlerBase::prim_path;

  using OptionHandlerBase::selection_path;

  using OptionHandlerBase::config_strs;

  using OptionHandlerBase::dof_strs;

  using OptionHandlerBase::coordtype_enum;

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
  using OptionHandlerBase::config_strs;
  using OptionHandlerBase::coordtype_enum;
  using OptionHandlerBase::selection_paths;
  using OptionHandlerBase::supercell_strs;

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
 * Options set for `casm view`. See what your casm landscape looks like here.
 */

class ViewOption : public OptionHandlerBase {
 public:
  using OptionHandlerBase::config_strs;
  using OptionHandlerBase::configtype;
  using OptionHandlerBase::configtype_opts;
  using OptionHandlerBase::selection_path;

  ViewOption();
  int m_images;

 private:
  void initialize() override;
};

class EnumOptionBase : public OptionHandlerBase {
 public:
  EnumOptionBase(std::string const &_name)
      : OptionHandlerBase(_name), m_all_existing(false) {}

  using OptionHandlerBase::config_strs;
  using OptionHandlerBase::coordtype_enum;
  using OptionHandlerBase::input_str;
  using OptionHandlerBase::settings_path;
  using OptionHandlerBase::supercell_strs;

  virtual ~EnumOptionBase() {}

  int min_volume() const { return m_min_volume; }

  int max_volume() const { return m_max_volume; }

  bool all_existing() const { return m_all_existing; }

 protected:
  int m_min_volume;
  int m_max_volume;
  bool m_all_existing;

 private:
  virtual void initialize() override {}
};

//*****************************************************************************************************//

/**
 * Options set for `casm sym`. Get your point groups here.
 */

class SymOption : public EnumOptionBase {
 public:
  using OptionHandlerBase::coordtype_enum;
  using OptionHandlerBase::coordtype_str;
  using OptionHandlerBase::dof_strs;
  using OptionHandlerBase::selection_path;

  SymOption();

  bool wedges() const { return m_wedges; }

  double tol() const { return m_tol; }

  fs::path poscar_path() const { return m_poscar_path; }

 private:
  bool m_wedges;
  double m_tol;
  fs::path m_poscar_path;

  void initialize() override;
};

/**
 * Options set for `casm enum`. Enumerate configurations,supercells, diff_trans
 * and diff_trans_configs here.
 */

class EnumOption : public EnumOptionBase {
 public:
  EnumOption();

  using OptionHandlerBase::dry_run;
  using OptionHandlerBase::verbosity_str;

  const std::vector<std::string> &desc_vec() const { return m_desc_vec; }

  std::string method() const { return m_method; }

  const std::string &filter_str() const { return m_filter_str; }

 private:
  void initialize() override;

  std::vector<std::string> m_desc_vec;

  std::string m_method;
  std::string m_filter_str;
};

//*****************************************************************************************************//

/**
 * Options set for `casm info`. Get info about CASM objects.
 */

class InfoOption : public OptionHandlerBase {
 public:
  using OptionHandlerBase::input_str;
  using OptionHandlerBase::settings_path;

  InfoOption();

  const std::vector<std::string> &desc_vec() const { return m_desc_vec; }

  std::string method() const { return m_method; }

 private:
  void initialize() override;

  std::vector<std::string> m_desc_vec;
  std::string m_method;
};

//*****************************************************************************************************//

/**
 * Options set for `casm query`. Get casm to query your things here.
 */

class QueryOption : public OptionHandlerBase {
 public:
  using OptionHandlerBase::db_type;
  using OptionHandlerBase::db_type_opts;
  using OptionHandlerBase::gzip_flag;
  using OptionHandlerBase::help_opt_vec;
  using OptionHandlerBase::output_path;
  using OptionHandlerBase::selection_path;

  QueryOption();

  bool verbatim_flag() const;

  const std::vector<std::string> &columns_vec() const;

  const std::vector<std::string> &new_alias_vec() const;

 private:
  void initialize() override;

  std::vector<std::string> m_columns_vec;

  std::vector<std::string> m_new_alias_vec;
};
//*****************************************************************************************************//

/**
 * Options set for `casm import`. Add new structures to your project here.
 */

class ImportOption : public OptionHandlerBase {
 public:
  using OptionHandlerBase::configtype;
  using OptionHandlerBase::configtype_opts;
  using OptionHandlerBase::input_str;
  using OptionHandlerBase::settings_path;

  ImportOption();

  const std::vector<fs::path> &pos_vec() const;

  const fs::path &batch_path() const;

 private:
  void initialize() override;

  std::vector<fs::path> m_pos_vec;

  fs::path m_batch_path;
};
//*****************************************************************************************************//

/**
 * Options set for `casm select`. Get your organized project here.
 */

class SelectOption : public OptionHandlerBase {
 public:
  using OptionHandlerBase::db_type;
  using OptionHandlerBase::db_type_opts;
  using OptionHandlerBase::help_opt_vec;
  using OptionHandlerBase::output_path;
  using OptionHandlerBase::selection_paths;

  SelectOption();

  const std::vector<std::string> &criteria_vec() const;

 private:
  void initialize() override;

  // vector necessary to allow --set/--set-on/--set-off with or without an
  // argument
  std::vector<std::string> m_criteria_vec;
};
//*****************************************************************************************************//

/**
 * Options set for `casm rm`. Get your clean project here.
 */

class RmOption : public OptionHandlerBase {
 public:
  using OptionHandlerBase::db_type;
  using OptionHandlerBase::db_type_opts;
  using OptionHandlerBase::dry_run;
  using OptionHandlerBase::name_strs;
  using OptionHandlerBase::selection_path;

  RmOption();

  bool force() const;

  bool data() const;

 private:
  void initialize() override;
};

//*****************************************************************************************************//

/**
 * Options set for `casm update`. Get your project up to date here.
 */

class UpdateOption : public OptionHandlerBase {
 public:
  using OptionHandlerBase::configtype;
  using OptionHandlerBase::configtype_opts;
  using OptionHandlerBase::input_str;
  using OptionHandlerBase::selection_path;
  using OptionHandlerBase::settings_path;

  UpdateOption();

  double vol_tolerance() const;

  double lattice_weight() const;

  double min_va_frac() const;

  double max_va_frac() const;

 private:
  void initialize() override;

  double m_vol_tolerance;

  double m_lattice_weight;  // TODO: Push to base? Other commands use this

  double m_min_va_frac;

  double m_max_va_frac;
};

//*****************************************************************************************************//

}  // namespace Completer
}  // namespace CASM

#endif
