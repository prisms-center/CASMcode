#ifndef PRIMCLEX_HH
#define PRIMCLEX_HH

#include "casm/external/boost.hh"

#include "casm/misc/cloneable_ptr.hh"
#include "casm/casm_io/Log.hh"

#include "casm/crystallography/Structure.hh"
#include "casm/clex/CompositionConverter.hh"
#include "casm/clex/Clexulator.hh"
#include "casm/clex/ChemicalReference.hh"
#include "casm/clex/NeighborList.hh"
#include "casm/clex/ClexBasis.hh"

#include "casm/app/DirectoryStructure.hh"
#include "casm/app/ProjectSettings.hh"
#include "casm/app/QueryHandler.hh"

#include "casm/database/DatabaseHandler.hh"

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
    DB::Database<T> &db() const {
      return m_db_handler->db<T>();
    }

    template<typename T>
    const DB::Database<T> &const_db() const {
      return m_db_handler->const_db<T>();
    }


    template<typename T>
    DB::Database<T> &db(std::string db_name) const {
      return m_db_handler->db<T>(db_name);
    }

    template<typename T>
    const DB::Database<T> &const_db(std::string db_name) const {
      return m_db_handler->const_db<T>(db_name);
    }


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

    DirectoryStructure m_dir;
    ProjectSettings m_settings;

    Structure m_prim;
    bool m_vacancy_allowed;
    Index m_vacancy_index;

    std::unique_ptr<DB::DatabaseHandler> m_db_handler;

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
