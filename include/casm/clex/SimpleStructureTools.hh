#ifndef CLEX_SIMPLESTRUCTURETOOLS_HH
#define CLEX_SIMPLESTRUCTURETOOLS_HH

#include <set>
#include <vector>

#include "casm/crystallography/DoFDecl.hh"
#include "casm/global/definitions.hh"
#include "casm/basis_set/DoFTraits.hh"

namespace CASM {

namespace xtal {
class SimpleStructure;
}

class ConfigDoF;
class Configuration;
struct MappedProperties;
class Supercell;

/// \brief Construct from ConfigDoF _dof belonging to provided Supercell _scel
/// @param _which_dofs List of DoF type-names that specify which DoF values to
/// utilize from _config when building the result
///        if empty, all DoFs are used. To exclude all DoFs from conversion,
///        pass _which_dofs={"none"}
xtal::SimpleStructure make_simple_structure(
    Supercell const &_scel, ConfigDoF const &_dof,
    std::vector<DoFKey> const &_which_dofs = {});

/// \brief Construct SimpleStructure from Configuration
/// @param _config Configuration used as source data
/// @param _which_dofs List of DoF type-names that specify which DoF values to
/// utilize from
///        _config when building the result. If empty, all DoFs are used. To
///        exclude all DoFs from conversion, pass _which_dofs={"none"}.
/// @param relaxed Flag specifying to use relaxed coordinates and parameters (if
/// true) or to only
///        use the imposed DoFs (if false).
///
/// Note:
/// - If `(relaxed==true && is_calculated(_config)`:
///   - MappedProperties from `_config.calc_properties()` are applied.
///   - Relaxed coordinates are obtained from
///   `_config.calc_properties().site.at("coordinate")`
///   - Relaxed lattice vectors  are obtained from
///   `_config.calc_properties().global.at("latvec")`
///   - It is expected, but not checked, that there exist no DoF whose
///   application further affects
///     the coordinates or lattice vectors
/// - Otherwise: no MappedProperties are applied.
xtal::SimpleStructure make_simple_structure(
    Configuration const &_config, std::vector<DoFKey> const &_which_dofs = {},
    bool relaxed = false);

/// \brief Construct from ConfigDoF _dof belonging to provided Supercell _scel
/// and using calculated properties
/// @param _props Record of calculated properties to use during conversion
/// @param relaxed flag specifying to use _props argument during conversion (if
/// true)
///
/// Note:
/// - If relaxed==true:
///   - Relaxed coordinates are obtained from
///   `_config.calc_properties().site.at("coordinate")`
///   - Relaxed lattice vectors  are obtained from
///   `_config.calc_properties().global.at("latvec")`
///   - It is expected, but not checked, that there exist no DoF whose
///   application further affects
///     the coordinates or lattice vectors.
/// - Otherwise: no MappedProperties are applied.
xtal::SimpleStructure make_simple_structure(
    Supercell const &_scel, ConfigDoF const &_dof,
    MappedProperties const &_props, std::vector<DoFKey> const &_which_dofs = {},
    bool relaxed = false);

/// \brief Determine which sites of a Configuration can host each atom of a
/// SimpleStructure result[i] is set of site indices in @param _config that can
/// host atom 'i' of @param sstruc
std::vector<std::set<Index> > atom_site_compatibility(
    xtal::SimpleStructure const &sstruc, Configuration const &_config);

/// \brief Determine which sites of a Configuration can host each molecule of a
/// SimpleStructure result[i] is set of site indices in @param _config that can
/// host molecule 'i' of @param sstruc
std::vector<std::set<Index> > mol_site_compatibility(
    xtal::SimpleStructure const &sstruc, Configuration const &_config);

namespace clex_SimpleStructureTools_impl {

/// Class that helps manage which order DoF are applied to a SimpleStructure
/// when building it from a Configuration. Each TransformDirective applies one
/// DoF type or "atomizes" molecules (populating SimpleStructure::atom_info from
/// SimpleStructure::mol_info). It is meant to be stored in a
/// std::set<TransformDirective> and has a comparison operator and uses
/// information from AnisoValTraits to appropriately order TransformDirective in
/// the set. Then they can be applied sequentially (using `transform`) to build
/// the SimpleStructure.
class TransformDirective {
 public:
  /// \brief consturct from transformation or DoF type name
  TransformDirective(std::string const &_name);

  /// \brief Name of DoFType or transformation
  std::string const &name() const { return m_name; }

  /// \brief Compare with _other TransformDirective. Returns true if this
  /// TransformDirective has precedence
  bool operator<(TransformDirective const &_other) const;

  /// \brief Applies transformation to _struc using information contained in
  /// _config
  void transform(ConfigDoF const &_config,
                 xtal::BasicStructure const &_reference,
                 xtal::SimpleStructure &_struc) const;

 private:
  /// \brief Build m_before object by recursively traversing DoF dependencies
  void _accumulate_before(std::set<std::string> const &_queue,
                          std::set<std::string> &_result) const;

  /// \brief Build m_after object by recursively traversing DoF dependencies
  void _accumulate_after(std::set<std::string> const &_queue,
                         std::set<std::string> &_result) const;

  std::string m_name;
  std::set<std::string> m_before;
  std::set<std::string> m_after;

  DoFType::Traits const *m_traits_ptr;
};
}
}  // namespace CASM

#endif
