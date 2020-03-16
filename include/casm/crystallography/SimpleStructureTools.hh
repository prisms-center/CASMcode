#ifndef SIMPLESTRUCTURETOOLS_HH
#define SIMPLESTRUCTURETOOLS_HH

#include <vector>
#include <string>
#include <set>
#include "casm/external/Eigen/Dense"
#include "casm/global/definitions.hh"
#include "casm/crystallography/DoFDecl.hh"


namespace CASM {
  namespace SimpleStructureTools {
    // Defined in SimpleStructure.hh
    class SpeciesMode;
  }


  namespace xtal {

    /** \ingroup Structure
     *  @{
     */
    class SimpleStructure;
    class Site;
    class BasicStructure;

    SimpleStructure make_superstructure(Eigen::Ref<const Eigen::Matrix3i> const &_T, SimpleStructure const &_sstruc);

    /// \brief Construct from decorated structure
    SimpleStructure make_simple_structure(BasicStructure const &_struc);

    std::vector<std::set<Index> > atom_site_compatibility(SimpleStructure const &sstruc, BasicStructure const &_prim);
    std::vector<std::set<Index> > mol_site_compatibility(SimpleStructure const &sstruc, BasicStructure const &_prim);

    /// \brief use information in _reference to initialize atom_info from mol_info
    void _atomize(SimpleStructure &_sstruc, Eigen::Ref<const Eigen::VectorXi> const &_mol_occ, BasicStructure const &_reference);
    /* SimpleStructure make_simple_structure(BasicStructure<Site> const &_struc); */

    /* /// \brief Construct SimpleStructure from Configuration */
    /* /// @param _config Configuration used as source data */
    /* /// @param _which_dofs List of DoF type-names that specify which DoF values to utilize from _config when building the result */
    /* ///        if empty, all DoFs are used. To exclude all DoFs from conversion, pass _which_dofs={"none"} */
    /* /// @param relaxed flag specifying to use relaxed coordinates and parameters (if true) or to only use the imposed DoFs (if false) */
    /* SimpleStructure make_simple_structure(Configuration const &_config, */
    /*                                       std::vector<DoFKey> const &_which_dofs = {}, */
    /*                                       bool relaxed = false); */

    /* /// \brief Construct from ConfigDoF _dof belonging to provided Supercell _scel */
    /* /// @param _which_dofs List of DoF type-names that specify which DoF values to utilize from _config when building the result */
    /* ///        if empty, all DoFs are used. To exclude all DoFs from conversion, pass _which_dofs={"none"} */
    /* SimpleStructure make_simple_structure(Supercell const &_scel, */
    /*                                       ConfigDoF const &_dof, */
    /*                                       std::vector<DoFKey> const &_which_dofs = {}); */

    /* /// \brief Construct from ConfigDoF _dof belonging to provided Supercell _scel and using calculated properties */
    /* /// @param _props Record of calculated properties to use during conversion */
    /* /// @param relaxed flag specifying to use _props argument during conversion (if true) */
    /* SimpleStructure make_simple_structure(Supercell const &_scel, */
    /*                                       ConfigDoF const &_dof, */
    /*                                       MappedProperties const &_props, */
    /*                                       std::vector<DoFKey> const &_which_dofs = {}, */
    /*                                       bool relaxed = false); */

    /* /// \brief Construct BasicStructure<Site> from SimpleStructure. */
    /* /// @param _sstruc SimpleStructure used as source data for conversion */
    /* /// @param _all_dofs holds names of additional DOFs to initialize in structure */
    /* /// @param mode specifies whether ATOM or MOL info of _sstruc should be used to build sites of structure */
    /* /// @param _allowed_occupants List of allowed molecules at each site; if empty, occupants are assumed to be atoms */
    /* ///        having the species names and attributes indicated by _sstruc */
    /* BasicStructure<Site> make_basic_structure(SimpleStructure const &_sstruc, */
    /*                                           std::vector<DoFKey> const &_all_dofs, */
    /*                                           SimpleStructureTools::SpeciesMode mode, */
    /*                                           std::vector<std::vector<Molecule> > _allowed_occupants = {}); */

    /* /// \brief Determine which sites of a BasicStructure can host each atom of a SimpleStructure */
    /* /// result[i] is set of site indices in @param _prim that can host atom 'i' of @param sstruc */
    /* std::vector<std::set<Index> > atom_site_compatibility(SimpleStructure const &sstruc, BasicStructure<Site> const &_prim); */

    /* /// \brief Determine which sites of a BasicStructure can host each molecule of a SimpleStructure */
    /* /// result[i] is set of site indices in @param _prim that can host molecule 'i' of @param sstruc */
    /* std::vector<std::set<Index> > mol_site_compatibility(SimpleStructure const &sstruc, BasicStructure<Site> const &_prim); */

    /* /// \brief Determine which sites of a Configuration can host each atom of a SimpleStructure */
    /* /// result[i] is set of site indices in @param _config that can host atom 'i' of @param sstruc */
    /* std::vector<std::set<Index> > atom_site_compatibility(SimpleStructure const &sstruc, Configuration const &_config); */

    /* /// \brief Determine which sites of a Configuration can host each molecule of a SimpleStructure */
    /* /// result[i] is set of site indices in @param _config that can host molecule 'i' of @param sstruc */
    /* std::vector<std::set<Index> > mol_site_compatibility(SimpleStructure const &sstruc, Configuration const &_config); */

    /* /// \brief Imposes DoF values from ConfigDoF _config onto *this, using using any necessary information contained in _reference */
    /* void _apply_dofs(SimpleStructure &_sstruc, ConfigDoF const &_config, BasicStructure<Site> const &_reference, std::vector<DoFKey> which_dofs); */

    /* /// \brief use information in _reference to initialize atom_info from mol_info */
    /* void _atomize(SimpleStructure &_sstruc, Eigen::Ref<const Eigen::VectorXi> const &_mol_occ, BasicStructure<Site> const &_reference); */

    /* /// \brief Output to JSON, excluding any molecular or atomic species contained in 'excluded_species' */
    /* jsonParser &to_json(xtal::SimpleStructure const &_struc, */
    /*                     jsonParser &json_supplement, */
    /*                     std::set<std::string> const &excluded_species = {"Va", "VA", "va"}, */
    /*                     std::string prefix = "", */
    /*                     COORD_TYPE mode = CART); */

  }

}
#endif
