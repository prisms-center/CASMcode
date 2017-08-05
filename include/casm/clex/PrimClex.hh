#ifndef CASM_PrimClex
#define CASM_PrimClex

#include "casm/CASM_global_definitions.hh"
#include "casm/casm_io/Log.hh"

/// Cluster expansion class
namespace CASM {

  class DirectoryStructure;
  class ProjectSettings;
  class ClexDescription;
  class CompositionConverter;
  class ChemicalReference;
  class Structure;
  class PrimNeighborList;
  class ClexBasis;
  class Clexulator;
  class ECIContainer;

  namespace DB {
    template<typename T> class ValDatabase;
    template<typename T> class Database;
    class PropertiesDatabase;
    class DatabaseHandler;
  }

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

  public:

    // **** Constructors ****

    /// Initial construction of a PrimClex, from a primitive Structure
    explicit PrimClex(const Structure &_prim, const Logging &logging = Logging());

    /// Construct PrimClex from existing CASM project directory
    ///  - read PrimClex and directory structure to generate all its Supercells and Configurations, etc.
    explicit PrimClex(const fs::path &_root, const Logging &logging = Logging());

    PrimClex(const PrimClex &) = delete;

    ~PrimClex();

    /// Reload PrimClex data from settings
    void refresh(bool read_settings = false,
                 bool read_composition = false,
                 bool read_chem_ref = false,
                 bool read_configs = false,
                 bool clear_clex = false);

    // ** Directory path and settings accessors **

    const DirectoryStructure &dir() const;

    ProjectSettings &settings();

    const ProjectSettings &settings() const;

    /// \brief Get the crystallography_tol
    double crystallography_tol() const;


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
    DB::ValDatabase<T> &generic_db() const;

    template<typename T>
    const DB::ValDatabase<T> &const_generic_db() const;


    template<typename T>
    DB::Database<T> &db() const;

    template<typename T>
    const DB::Database<T> &const_db() const;


    template<typename T>
    DB::PropertiesDatabase &db_props(std::string calc_type) const;

    template<typename T>
    const DB::PropertiesDatabase &const_db_props(std::string calc_type) const;


    DB::DatabaseHandler &db_handler() const;

    const DB::DatabaseHandler &const_db_handler() const;


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

    /// To avoid excessive includes, use "pointer to data"
    /// - PrimClexData is defined in PrimClex.cc
    struct PrimClexData;

    /// Initialization routines
    void _init();

    std::unique_ptr<PrimClexData> m_data;

  };

}
#endif
