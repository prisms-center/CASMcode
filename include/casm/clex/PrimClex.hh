#ifndef CASM_PrimClex
#define CASM_PrimClex

#include <memory>

#include "casm/global/definitions.hh"

/// Cluster expansion class
namespace CASM {

class DirectoryStructure;
class ProjectSettings;
struct ClexDescription;
class CompositionConverter;
class ChemicalReference;
class PrimNeighborList;
struct ClexBasisSpecs;
class ClusterSpecs;
class ClexBasis;
class Clexulator;
class ECIContainer;
struct NeighborhoodInfo;
class Structure;

namespace DB {
template <typename T>
class ValDatabase;
template <typename T>
class Database;
class PropertiesDatabase;
class DatabaseHandler;
}  // namespace DB

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
class PrimClex {
 public:
  typedef Structure PrimType;

  // **** Constructors ****

  /// Initial construction of a PrimClex, from ProjectSettings and shared prim
  explicit PrimClex(ProjectSettings const &_project_settings,
                    std::shared_ptr<PrimType const> _shared_prim);

  /// Construct PrimClex from existing CASM project directory
  ///  - read PrimClex and directory structure to generate all its Supercells
  ///  and Configurations, etc.
  explicit PrimClex(const fs::path &_root);

  PrimClex(const PrimClex &) = delete;

  ~PrimClex();

  /// Reload PrimClex data from settings
  void refresh(bool read_settings = false, bool read_composition = false,
               bool read_chem_ref = false, bool read_configs = false,
               bool clear_clex = false);

  // ** Directory path and settings accessors **

  ProjectSettings &settings();

  const ProjectSettings &settings() const;

  /// Check if DirectoryStructure exists
  bool has_dir() const;

  /// Access DirectoryStructure object. Throw if not set.
  const DirectoryStructure &dir() const;

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
  const PrimType &prim() const;

  /// Access to the primitive Structure as a shared resource
  std::shared_ptr<PrimType const> const &shared_prim() const;

  /// const Access to number of basis atoms
  Index n_basis() const;

  /// Access to the primitive neighbor list as a shared resource
  std::shared_ptr<PrimNeighborList> const &shared_nlist() const;

  /// Access to the primitive neighbor list
  PrimNeighborList &nlist() const;

  /// returns true if vacancy are an allowed species
  bool vacancy_allowed() const;

  /// returns the index of vacancies in composition vectors
  Index vacancy_index() const;

  // ** Supercell, Configuration, etc. databases **

  template <typename T>
  DB::ValDatabase<T> &generic_db() const;

  template <typename T>
  const DB::ValDatabase<T> &const_generic_db() const;

  template <typename T>
  DB::Database<T> &db() const;

  template <typename T>
  const DB::Database<T> &const_db() const;

  template <typename T>
  DB::PropertiesDatabase &db_props(std::string calc_type) const;

  template <typename T>
  const DB::PropertiesDatabase &const_db_props(std::string calc_type) const;

  DB::DatabaseHandler &db_handler() const;

  const DB::DatabaseHandler &const_db_handler() const;

  // ** Basis set specs, clexulators, eci ** ------------------

  bool has_basis_set_specs(std::string const &basis_set_name) const;
  ClexBasisSpecs const &basis_set_specs(
      std::string const &basis_set_name) const;
  NeighborhoodInfo const &neighborhood_info(
      std::string const &basis_set_name) const;
  Clexulator clexulator(std::string const &basis_set_name) const;

  bool has_eci(const ClexDescription &key) const;
  ECIContainer const &eci(const ClexDescription &key) const;

 private:
  /// To avoid excessive includes, use "pointer to data"
  /// - PrimClexData is defined in PrimClex.cc
  struct PrimClexData;

  /// Initialization routines
  void _init();

  std::unique_ptr<PrimClexData> m_data;
};

/// Write clust.json, basis.json, and clexulator source code, by constructing
/// orbits and ClexBasis
///
/// Notes:
/// - Overwrites any existing files
/// - Uses DoFType::traits_dict() for DoFTraits
void write_basis_set_data(std::shared_ptr<Structure const> shared_prim,
                          ProjectSettings const &settings,
                          std::string const &basis_set_name,
                          ClexBasisSpecs const &basis_set_specs,
                          PrimNeighborList &prim_neighbor_list);

/// Make Clexulator from existing source code
///
/// Notes:
/// - Use `write_basis_set_data` to write Clexulator source code prior to
/// calling this function.
/// - This function will compile the Clexulator source code if that has not yet
/// been done.
Clexulator make_clexulator(ProjectSettings const &settings,
                           std::string const &basis_set_name,
                           PrimNeighborList &prim_neighbor_list);

}  // namespace CASM
#endif
