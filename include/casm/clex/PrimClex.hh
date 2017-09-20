#ifndef PRIMCLEX_HH
#define PRIMCLEX_HH

#include "casm/external/boost.hh"

#include "casm/misc/cloneable_ptr.hh"
#include "casm/casm_io/Log.hh"

#include "casm/crystallography/Structure.hh"
#include "casm/clex/DoFManager.hh"
#include "casm/clex/CompositionConverter.hh"
#include "casm/clex/Supercell.hh"
#include "casm/clex/Clexulator.hh"
#include "casm/clex/ChemicalReference.hh"
#include "casm/misc/cloneable_ptr.hh"
#include "casm/clex/NeighborList.hh"

#include "casm/app/DirectoryStructure.hh"
#include "casm/app/ProjectSettings.hh"
#include "casm/app/EnumeratorHandler.hh"
#include "casm/app/QueryHandler.hh"

/// Cluster expansion class
namespace CASM {

  class ECIContainer;

  template<typename T, typename U> class ConfigIterator;

  /** \defgroup Clex
   *
   * \brief The functions and classes related to evaluating cluster expansions
   */

  /** \defgroup PrimClex
   *  \ingroup Clex
   *  \ingroup Project
   *  \brief PrimClex is the top-level data structure for a CASM project
   *
   *  @{
   */

  /// \brief PrimClex is the top-level data structure for a CASM project
  ///
  /// The PrimClex provides access to:
  /// - the parent crystal structure, the `prim`
  /// - project files & settings
  /// - composition axes & references
  /// - supercells & configurations
  /// - clusters, basis sets, neighbor lists, & ECI
  /// -
  class PrimClex : public Logging {

    fs::path root;

    DirectoryStructure m_dir;
    ProjectSettings m_settings;

    Structure prim;
    bool m_vacancy_allowed;
    Index m_vacancy_index;

    mutable DoFManager m_dof_manager;

    /// Contains all the supercells that were involved in the enumeration.
    boost::container::stable_vector< Supercell > supercell_list;


    /// CompositionConverter specifies parameteric composition axes and converts between
    ///   parametric composition and mol composition
    bool m_has_composition_axes = false;
    CompositionConverter m_comp_converter;

    /// ChemicalReference specifies a reference for formation energies, chemical
    /// potentials, etc.
    notstd::cloneable_ptr<ChemicalReference> m_chem_ref;

    /// Stores the neighboring UnitCell and which sublattices to include in neighbor lists
    /// - mutable for lazy construction
    mutable notstd::cloneable_ptr<PrimNeighborList> m_nlist;


  public:

    typedef ConfigIterator<Configuration, PrimClex> config_iterator;
    typedef ConfigIterator<const Configuration, const PrimClex> config_const_iterator;

    // **** Constructors ****

    /// Initial construction of a PrimClex, from a primitive Structure
    PrimClex(const Structure &_prim, const Logging &logging = Logging());

    /// Construct PrimClex from existing CASM project directory
    ///  - read PrimClex and directory structure to generate all its Supercells and Configurations, etc.
    PrimClex(const fs::path &_root, const Logging &logging = Logging());

    /// Reload PrimClex data from settings
    void refresh(bool read_settings = false,
                 bool read_composition = false,
                 bool read_chem_ref = false,
                 bool read_configs = false,
                 bool clear_clex = false);

    // **** Accessors ****

    /// Return project name
    std::string name() const;


    // ** Directory path accessors **

    const DirectoryStructure &dir() const {
      return m_dir;
    }

    ProjectSettings &settings() {
      return m_settings;
    }

    const ProjectSettings &settings() const {
      return m_settings;
    }

    double crystallography_tol() const {
      return settings().crystallography_tol();
    }

    double lin_alg_tol() const {
      return settings().lin_alg_tol();
    }

    /// Return casm project directory path
    fs::path get_path() const;

    /// Return supercell directory path
    fs::path get_path(const Index &scel_index) const;

    /// Return configuration directory path
    fs::path get_path(const Index &scel_index, const Index &config_index) const;

    /// Return config_list.json file path
    fs::path get_config_list_path() const;


    // ** Composition accessors **

    /// check if CompositionConverter object initialized
    bool has_composition_axes() const;

    /// const Access CompositionConverter object
    const CompositionConverter &composition_axes() const;

    // ** Chemical reference **

    /// check if ChemicalReference object initialized
    bool has_chemical_reference() const;

    /// const Access ChemicalReference object
    const ChemicalReference &chemical_reference() const;


    // ** Prim and Orbitree accessors **

    /// const Access to primitive Structure
    const Structure &get_prim() const;

    ///Access to the primitive neighbor list
    PrimNeighborList &nlist() const;

    /// returns true if vacancy are an allowed species
    bool vacancy_allowed() const;

    /// returns the index of vacancies in composition vectors
    Index vacancy_index() const;

    // ** Supercell and Configuration accessors **

    /// Access entire supercell_list
    boost::container::stable_vector<Supercell> &get_supercell_list();

    /// const Access entire supercell_list
    const boost::container::stable_vector<Supercell> &get_supercell_list() const;

    /// const Access supercell by index
    const Supercell &get_supercell(Index i) const;

    /// Access supercell by index
    Supercell &get_supercell(Index i);

    /// const Access supercell by name
    const Supercell &get_supercell(std::string scellname) const;

    /// Access supercell by name
    Supercell &get_supercell(std::string scellname);

    /// Access supercell by Lattice, adding if necessary
    Supercell &get_supercell(const Lattice &lat);

    /// access configuration by name (of the form "scellname/[NUMBER]", e.g., ("SCEL1_1_1_1_0_0_0/0")
    const Configuration &configuration(const std::string &configname) const;
    Configuration &configuration(const std::string &configname);

    /// Configuration iterator: begin
    config_iterator config_begin();

    /// Configuration iterator: end
    config_iterator config_end();

    /// Configuration iterator: begin
    config_const_iterator config_begin() const;

    /// Configuration iterator: end
    config_const_iterator config_end() const;

    /// const Configuration iterator: begin
    config_const_iterator config_cbegin() const;

    /// const Configuration iterator: end
    config_const_iterator config_cend() const;

    /// Configuration iterator: begin
    config_iterator selected_config_begin();

    /// Configuration iterator: end
    config_iterator selected_config_end();

    /// const Configuration iterator: begin
    config_const_iterator selected_config_cbegin() const;

    /// const Configuration iterator: end
    config_const_iterator selected_config_cend() const;


    Eigen::MatrixXd shift_vectors() const;


    // **** Mutators ****

    /// Sets the composition axes, updates all configuration references,
    ///   and writes the updated configuration info
    void set_composition_axes(const CompositionConverter &_converter);

    // **** IO ****

    ///Call Configuration::write on every configuration to update files
    ///  - call update to also read all files
    void write_config_list(std::set<std::string> scel_to_delete = {});


    // **** Operators ****

    // **** Functions for preparing CLEXulators ****

    //Generate the global orbitree
    //John G 011013
    /// Use the given CSPECS

    /// \brief Generate supercells of a certain volume and shape and store them in the array of supercells
    void generate_supercells(const ScelEnumProps &enum_props);

    //Enumerate configurations for all the supercells that are stored in 'supercell_list'
    void print_enum_info(std::ostream &stream);
    void print_supercells(std::set<std::string> scel_to_delete = {}) const;
    void print_supercells(std::ostream &stream, std::set<std::string> scel_to_delete = {}) const;
    void read_supercells(std::istream &stream);
    void print_clex_configurations();


    //ParamComposition i/o and calculators in PrimClex

    void read_config_list();

    ///Fill up props of every configuration for a partucluar supercell. This will be deprecated when props disappears
    void read_scel_props(int scel_index, const std::string &JSON_output);
    ///Call read_config_props on every Supercell
    void read_all_scel_props(const std::string &JSON_output);

    ///Count over the number of configurations that are selected in all supercells
    int amount_selected() const;

    bool contains_supercell(std::string scellname, Index &index) const;

    bool contains_supercell(const Supercell &scel) const;
    bool contains_supercell(const Supercell &scel, Index &index) const;

    Index add_supercell(const Lattice &superlat);

    Index add_canonical_supercell(const Lattice &superlat);

    Eigen::Matrix3i calc_transf_mat(const Lattice &superlat) const;

    /// Set internal values of each DoFEnvironment
    void set_global_dof_state(const Configuration &curr_config)const {
      m_dof_manager.set_global_dof_state(curr_config);
    };

    void set_local_dof_state(const Configuration &curr_config, Index l)const {
      m_dof_manager.set_local_dof_state(curr_config, l);
    };

    /// Delete 'properties.ref_state.X.json' files,
    /// Then call 'clear_reference_properties'
    //void clear_reference_states();

    /// Sets the root reference state to be the calculated properties of the chosen config
    /// Calls 'clear_reference_properties'
    //void set_reference_state(int refid, const Configuration &config);

    /// Check that it is valid to use 'config' as reference state 'refid', returns bool and if false, sets 'reason_invalid'
    ///   Currently checks:
    ///     1) that the necessary properties have been calculated,
    ///     2) that the same Configuration is not being used twice
    ///   Needs to check that reference states span composition space
    //bool valid_reference_state(int refid, const Configuration &config, std::string &reason_invalid) const;

    /// find calculated configurations closest to
    /// [0, 0, 0, ...], [1, 0, 0, ...], [0, 1, 0, ...], [0, 0, 1, ...], ...
    /// and set them as the root reference states, also calls regenerate_references
    /// Clears reference states and properties whether or not it succeeds
    //void set_reference_state_auto();

    /// Clear 'reference' and 'delta' properties from all Configurations
    /// Re-write all Configurations, updating:
    ///   param_composition.json
    ///   properties.calc.json
    ///   properties.ref.json
    ///   properties.delta.json
    //void generate_references();

    bool has_orbitree(const ClexDescription &key) const;
    const SiteOrbitree &orbitree(const ClexDescription &key) const;

    bool has_clexulator(const ClexDescription &key) const;
    Clexulator clexulator(const ClexDescription &key) const;

    bool has_eci(const ClexDescription &key) const;
    const ECIContainer &eci(const ClexDescription &key) const;

  private:

    /// Initialization routines
    void _init();

    mutable std::map<ClexDescription, SiteOrbitree> m_orbitree;
    mutable std::map<ClexDescription, Clexulator> m_clexulator;
    mutable std::map<ClexDescription, ECIContainer> m_eci;

  };


  /// \brief Make orbitree. For now specifically global.
  SiteOrbitree make_orbitree(Structure &prim, const jsonParser &json, double _tol);

  /// \brief Print clexulator
  void print_clexulator(const Structure &prim,
                        SiteOrbitree &tree,
                        const PrimNeighborList &nlist,
                        std::string class_name,
                        std::ostream &stream,
                        double xtal_tol);

}
#endif
