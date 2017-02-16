#ifndef PRIMCLEX_HH
#define PRIMCLEX_HH

#include "casm/external/boost.hh"

#include "casm/misc/cloneable_ptr.hh"
#include "casm/casm_io/Log.hh"

#include "casm/crystallography/Structure.hh"
#include "casm/clex/CompositionConverter.hh"
#include "casm/clex/Supercell.hh"
#include "casm/clex/Clexulator.hh"
#include "casm/clex/ChemicalReference.hh"
#include "casm/clex/NeighborList.hh"
#include "casm/clex/ClexBasis.hh"

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

    DirectoryStructure m_dir;
    ProjectSettings m_settings;

    Structure m_prim;
    bool m_vacancy_allowed;
    Index m_vacancy_index;

    /// Contains all the supercells that were involved in the enumeration.
    boost::container::stable_vector< Supercell > m_supercell_list;


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
    explicit PrimClex(const Structure &_prim, const Logging &logging = Logging());

    /// Construct PrimClex from existing CASM project directory
    ///  - read PrimClex and directory structure to generate all its Supercells and Configurations, etc.
    explicit PrimClex(const fs::path &_root, const Logging &logging = Logging());

    /// Reload PrimClex data from settings
    void refresh(bool read_settings = false,
                 bool read_composition = false,
                 bool read_chem_ref = false,
                 bool read_configs = false,
                 bool clear_clex = false);

    // ** Directory path and settings accessors **

    const DirectoryStructure &dir() const {
      return m_dir;
    }

    ProjectSettings &settings() {
      return m_settings;
    }

    const ProjectSettings &settings() const {
      return m_settings;
    }

    /// \brief Get the crystallography_tol
    double crystallography_tol() const {
      return settings().crystallography_tol();
    }


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


    // ** Accessors **

    /// const Access to primitive Structure
    const Structure &prim() const;

    ///Access to the primitive neighbor list
    PrimNeighborList &nlist() const;

    /// returns true if vacancy are an allowed species
    bool vacancy_allowed() const;

    /// returns the index of vacancies in composition vectors
    Index vacancy_index() const;


    // ** Supercell, Configuration, etc. databases **

    template<typename T>
    Database<T>& db();
    
    template<typename T>
    const Database<T>& db() const;
    
    template<typename T>
    const Database<T>& const_db();
    
    
    template<typename T>
    Database<T>& db();
    
    template<typename T>
    const Database<T>& db() const;
    
    template<typename T>
    const Database<T>& const_db();
    
    
    // **** IO ****

    ///Call Configuration::write on every configuration to update files
    ///  - call update to also read all files
    void write_config_list();


    // **** Operators ****

    // **** Functions for preparing CLEXulators ****

    /// \brief Generate supercells of a certain volume and shape and store them in the array of supercells
    void generate_supercells(const ScelEnumProps &enum_props);

    //Enumerate configurations for all the supercells that are stored in 'supercell_list'
    void print_enum_info(std::ostream &stream);
    void print_supercells() const;
    void print_supercells(std::ostream &stream) const;
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


    bool has_orbits(const ClexDescription &key) const;
    template<typename OrbitOutputIterator, typename SymCompareType>
    OrbitOutputIterator orbits(const ClexDescription &key,
                               OrbitOutputIterator result,
                               const SymCompareType &sym_compare) const;

    bool has_clex_basis(const ClexDescription &key) const;
    const ClexBasis &clex_basis(const ClexDescription &key) const;

    bool has_clexulator(const ClexDescription &key) const;
    Clexulator clexulator(const ClexDescription &key) const;

    bool has_eci(const ClexDescription &key) const;
    const ECIContainer &eci(const ClexDescription &key) const;

  private:

    /// Initialization routines
    void _init();

    mutable std::map<ClexDescription, ClexBasis> m_clex_basis;
    mutable std::map<ClexDescription, Clexulator> m_clexulator;
    mutable std::map<ClexDescription, ECIContainer> m_eci;

  };

  //*******************************************************************************************
  template<typename OrbitOutputIterator, typename SymCompareType>
  OrbitOutputIterator PrimClex::orbits(
    const ClexDescription &key,
    OrbitOutputIterator result,
    const SymCompareType &sym_compare) const {

    return read_clust(
             result,
             jsonParser(dir().clust(key.bset)),
             prim(),
             prim().factor_group(),
             sym_compare,
             settings().crystallography_tol());
  }


}
#endif
