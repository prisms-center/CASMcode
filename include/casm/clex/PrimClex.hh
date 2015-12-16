#ifndef PRIMCLEX_HH
#define PRIMCLEX_HH

#define BOOST_NO_SCOPED_ENUMS
#define BOOST_NO_CXX11_SCOPED_ENUMS
#include <boost/filesystem.hpp>

#include "casm/BP_C++/BP_Parse.hh"

#include "casm/crystallography/Structure.hh"
#include "casm/clex/DoFManager.hh"
#include "casm/clex/ParamComposition.hh"
#include "casm/clex/Supercell.hh"
#include "casm/clex/Clexulator.hh"

#include "casm/app/DirectoryStructure.hh"
#include "casm/app/ProjectSettings.hh"
#include "casm/app/AppIO.hh"

/// Cluster expansion class
namespace CASM {

  class ParamComposition;
  class ECIContainer;

  template<typename T, typename U> class ConfigIterator;
  
  
  /// \defgroup Clex
  ///    
  /// \brief A Configuration represents the values of all degrees of freedom in a Supercell
  ///    
  
  
  /// \brief PrimClex stores the primitive Structure and lots of related data
  ///
  /// \ingroup Clex
  ///
  class PrimClex {

    fs::path root;

    DirectoryStructure m_dir;
    ProjectSettings m_settings;

    std::string m_name;

    Structure prim;

    mutable DoFManager m_dof_manager;


    // CASM project current settings: used to determine where to write things
    std::vector<std::string> curr_property;
    std::string curr_clex;
    std::string curr_calctype;
    std::string curr_ref;
    std::string curr_bset;
    std::string curr_eci;

    // Runtime library compilation settings: compilation options
    std::string compile_options;
    std::string so_options;

    SiteOrbitree global_orbitree;

    //One flowertree for each site. The flowertrees contain all the clusters that have a particular site at the pivot, divided in branches of points, pairs etc. Probably what you
    //thought a bouquet was to begin with (but you'd be wrong because bouquets have a branch for each basis site)
    Array<SiteOrbitree> flowertrees;


    /// Contains all the supercells that were involved in the enumeration.
    boost::container::stable_vector< Supercell > supercell_list;


    /// CompositionConverter specifies parameteric composition axes and converts between
    ///   parametric composition and mol composition
    bool m_has_composition_axes = false;
    CompositionConverter m_comp_converter;


    /// Stores the 'delta' UnitCellCoord needed to determine all the
    ///   sites in the neighborhood of a given primitive cell according to:
    ///
    ///    neighbor_unit_cell_coord = UnitCellCoord(b,i,j,k) + prim_nlist[nlist_index]
    ///
    ///  'delta' UnitCellCoord = prim_nlist[nlist_index]
    ///

    Array<UnitCellCoord> prim_nlist;

  public:

    typedef ConfigIterator<Configuration, PrimClex> config_iterator;
    typedef ConfigIterator<const Configuration, const PrimClex> config_const_iterator;
    //typedef ConfigIterator<Transition> trans_iterator;
    //typedef ConfigIterator<const Transition> trans_const_iterator;

    // **** Constructors ****

    /// Initial construction of a PrimClex, from a primitive Structure
    PrimClex(const Structure &_prim);

    /// Construct PrimClex from existing CASM project directory
    ///  - read PrimClex and directory structure to generate all its Supercells and Configurations, etc.
    PrimClex(const fs::path &_root, std::ostream &sout);



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

    double tol() const {
      return settings().tol();
    }

    /// Return casm project directory path
    fs::path get_path() const;

    /// Return supercell directory path
    fs::path get_path(const Index &scel_index) const;

    /// Return configuration directory path
    fs::path get_path(const Index &scel_index, const Index &config_index) const;

    /// Return config_list.json file path
    fs::path get_config_list_path() const;

    // ** Current settings accessors **

    /// Return current property settings
    const std::vector<std::string> &get_curr_property() const;

    /// Return current clex settings
    std::string get_curr_clex() const;

    /// Return current calctype setting
    std::string get_curr_calctype() const;

    /// Return current reference setting
    std::string get_curr_ref() const;

    /// Return basis set settings
    std::string get_curr_bset() const;

    /// Return current global clexulator name
    std::string get_curr_clexulator() const;

    /// Return current eci settings
    std::string get_curr_eci() const;

    /// Return compiler options
    std::string get_compile_options() const;

    /// Return shared library options
    std::string get_so_options() const;

    // ** Composition accessors **

    /// const Access CompositionConverter object
    bool has_composition_axes() const;

    /// const Access CompositionConverter object
    const CompositionConverter &composition_axes() const;


    // ** Prim and Orbitree accessors **

    /// const Access to primitive Structure
    const Structure &get_prim() const;

    /// const Access to global orbitree
    const SiteOrbitree &get_global_orbitree() const;

    ///const access to the primitive neighbor list
    const Array<UnitCellCoord> &get_prim_nlist() const;


    // ** Supercell and Configuration accessors **

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

    /*
    /// Transition iterator: begin
    trans_iterator trans_begin();

    /// Transition iterator: end
    trans_iterator trans_end();

    /// const Transition iterator: begin
    trans_const_iterator trans_cbegin() const;

    /// const Transition iterator: end
    trans_const_iterator trans_cend() const;
    */


    // ** Neighbor list accessors **

    /// Returns number of neighbors
    Index get_nlist_size() const;

    /// Returns UnitCellCoord, 'delta', indicating where neighbor site 'nlist_index' is located in the neighborhood
    ///    neighbor_unit_cell_coord = UnitCellCoord(b,i,j,k) + this->get_nlist_uccoord(nlist_index)
    const UnitCellCoord &get_nlist_uccoord(Index nlist_index) const;

    Eigen::MatrixXd shift_vectors() const;

    // **** Mutators ****


    // **** IO ****

    ///Call Configuration::write on every configuration to update files
    ///  - call update to also read all files
    void write_config_list();

    /// \brief Set the primitive neighbor list explicitly, useful when it has been saved
    void set_prim_nlist(const Array<UnitCellCoord> &_prim_nlist) {
      prim_nlist = _prim_nlist;
    }


    // **** Operators ****

    // **** Functions for preparing CLEXulators ****

    //Generate the global orbitree
    //John G 011013
    /// Use the given CSPECS

    //Read the global Orbitree from a clust.json file
    void read_global_orbitree(const fs::path &fclust);


    // Prepare neighbor lists
    /// Add one orbitree of sites to the nearest neighbor list and update sites in tree with index, only considering the first site of each cluster
    void append_to_nlist(SiteOrbitree &new_tree);
    /// Add one orbitree of sites to the nearest neighbor list and update sites in tree with index
    void append_to_nlist_perm(SiteOrbitree &new_tree);
    /// Add global_orbitree to nearest neighbor list and update sites in tree with index
    void generate_full_nlist();
    /// Generate basis flowers from global orbitree and divide them nicely into flowertrees
    void populate_flowertrees();
    //\John G

    //Generate supercells of a certain volume and store them in the array of supercells
    void generate_supercells(int volStart, int volEnd, bool verbose);

    //Enumerate configurations for all the supercells that are stored in 'supercell_list'
    void print_enum_info(std::ostream &stream);
    void print_supercells() const;
    void print_supercells(std::ostream &stream) const;
    void read_supercells(std::istream &stream);
    void print_clex_configurations();

    void generate_supercell_nlists();

    //ParamComposition i/o and calculators in PrimClex

    /// Sets the composition axes, updates all configuration references,
    ///   and writes the updated configuration info
    void set_composition_axes(const CompositionConverter &_converter);

    void read_config_list();

    ///Fill up props of every configuration for a partucluar supercell. This will be deprecated when props disappears
    void read_scel_props(int scel_index, const std::string &JSON_output);
    ///Call read_config_props on every Supercell
    void read_all_scel_props(const std::string &JSON_output);

    ///Count over the number of configurations that are selected in all supercells
    int amount_selected() const;

    bool contains_supercell(std::string scellname, Index &index) const;

    Index add_supercell(const Lattice &superlat);

    Index add_canonical_supercell(const Lattice &superlat);

    Matrix3<int> calc_transf_mat(const Lattice &superlat) const;

    /// Set internal values of each DoFEnvironment
    void set_global_dof_state(const Configuration &curr_config)const {
      m_dof_manager.set_global_dof_state(curr_config);
    };

    void set_local_dof_state(const Configuration &curr_config, Index l)const {
      m_dof_manager.set_local_dof_state(curr_config, l);
    };

    /// Delete 'properties.ref_state.X.json' files,
    /// Then call 'clear_reference_properties'
    void clear_reference_states();

    /// Sets the root reference state to be the calculated properties of the chosen config
    /// Calls 'clear_reference_properties'
    void set_reference_state(int refid, const Configuration &config);

    /// Check that it is valid to use 'config' as reference state 'refid', returns bool and if false, sets 'reason_invalid'
    ///   Currently checks:
    ///     1) that the necessary properties have been calculated,
    ///     2) that the same Configuration is not being used twice
    ///   Needs to check that reference states span composition space
    bool valid_reference_state(int refid, const Configuration &config, std::string &reason_invalid) const;

    /// find calculated configurations closest to
    /// [0, 0, 0, ...], [1, 0, 0, ...], [0, 1, 0, ...], [0, 0, 1, ...], ...
    /// and set them as the root reference states, also calls regenerate_references
    /// Clears reference states and properties whether or not it succeeds
    void set_reference_state_auto();

    /// Clear 'reference' and 'delta' properties from all Configurations
    /// Re-write all Configurations, updating:
    ///   param_composition.json
    ///   properties.calc.json
    ///   properties.ref.json
    ///   properties.delta.json
    void generate_references();

    Clexulator global_clexulator() const;
    ECIContainer global_eci(std::string clex_name) const;
  private:

    /// Return the configuration closest in param_composition to the target_param_comp
    ///   Tie break returns configuration in smallest supercell (first found at that size)
    const Configuration &closest_calculated_config(const Eigen::VectorXd &target_param_comp) const;


    mutable Clexulator m_global_clexulator;
  };


  /// \brief Make orbitree. For now specifically global.
  SiteOrbitree make_orbitree(Structure &prim, const jsonParser &json);

  /// \brief Print clexulator
  void print_clexulator(const Structure &prim,
                        SiteOrbitree &tree,
                        const Array<UnitCellCoord> &nlist,
                        std::string class_name,
                        std::ostream &stream);


  /// \brief Expand a neighbor list to include neighborhood of another SiteOrbitree
  void expand_nlist(const Structure &prim,
                    SiteOrbitree &tree,
                    Array<UnitCellCoord> &nlist);
}
#endif
