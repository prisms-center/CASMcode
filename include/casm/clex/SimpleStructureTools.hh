#ifndef CLEX_SIMPLESTRUCTURETOOLS_HH
#define CLEX_SIMPLESTRUCTURETOOLS_HH

#include "casm/basis_set/DoFTraits.hh"
#include <string>
namespace CASM {
  class Supercell;
  class ConfigDoF;
  class Configuration;
  class MappedProperties;

  namespace DoFType {
    class Traits;
  }

  namespace xtal {
    class BasicStructure;
  }

  //TODO: What is this?
  class TransformDirective {
  public:

    /// \brief consturct from transformation or DoF type name
    TransformDirective(std::string const &_name);

    /// \brief Name of DoFType or transformation
    std::string const &name() const {
      return m_name;
    }

    /// \brief Compare with _other TransformDirective. Returns true if this TransformDirective has precedence
    bool operator<(TransformDirective const &_other) const;

    /// \brief Applies transformation to _struc using information contained in _config
    void transform(ConfigDoF const  &_config, xtal::BasicStructure const &_reference, SimpleStructure &_struc) const;

  private:
    /// \brief Build m_before object by recursively traversing DoF dependencies
    void _accumulate_before(std::set<std::string>const &_queue, std::set<std::string> &_result) const;

    /// \brief Build m_after object by recursively traversing DoF dependencies
    void _accumulate_after(std::set<std::string>const &_queue, std::set<std::string> &_result) const;

    std::string m_name;
    std::set<std::string> m_before;
    std::set<std::string> m_after;

    DoFType::Traits const *m_traits_ptr;
  };

  /// \brief Construct SimpleStructure from Configuration
  /// @param _config Configuration used as source data
  /// @param _which_dofs List of DoF type-names that specify which DoF values to utilize from _config when building the result
  ///        if empty, all DoFs are used. To exclude all DoFs from conversion, pass _which_dofs={"none"}
  /// @param relaxed flag specifying to use relaxed coordinates and parameters (if true) or to only use the imposed DoFs (if false)
  SimpleStructure make_simple_structure(Configuration const &_config,
                                        std::vector<DoFKey> const &_which_dofs = {},
                                        bool relaxed = false);

  /// \brief Construct from ConfigDoF _dof belonging to provided Supercell _scel
  /// @param _which_dofs List of DoF type-names that specify which DoF values to utilize from _config when building the result
  ///        if empty, all DoFs are used. To exclude all DoFs from conversion, pass _which_dofs={"none"}
  SimpleStructure make_simple_structure(Supercell const &_scel,
                                        ConfigDoF const &_dof,
                                        std::vector<DoFKey> const &_which_dofs = {});

  /// \brief Construct from ConfigDoF _dof belonging to provided Supercell _scel and using calculated properties
  /// @param _props Record of calculated properties to use during conversion
  /// @param relaxed flag specifying to use _props argument during conversion (if true)
  SimpleStructure make_simple_structure(Supercell const &_scel,
                                        ConfigDoF const &_dof,
                                        MappedProperties const &_props,
                                        std::vector<DoFKey> const &_which_dofs = {},
                                        bool relaxed = false);

  /// \brief Determine which sites of a Configuration can host each atom of a SimpleStructure
  /// result[i] is set of site indices in @param _config that can host atom 'i' of @param sstruc
  std::vector<std::set<Index> > atom_site_compatibility(SimpleStructure const &sstruc, Configuration const &_config);

  /// \brief Determine which sites of a Configuration can host each molecule of a SimpleStructure
  /// result[i] is set of site indices in @param _config that can host molecule 'i' of @param sstruc
  std::vector<std::set<Index> > mol_site_compatibility(SimpleStructure const &sstruc, Configuration const &_config);

  /// \brief Imposes DoF values from ConfigDoF _config onto *this, using using any necessary information contained in _reference
  void _apply_dofs(SimpleStructure &_sstruc, ConfigDoF const &_config, BasicStructure const &_reference, std::vector<DoFKey> which_dofs);

}

#endif
