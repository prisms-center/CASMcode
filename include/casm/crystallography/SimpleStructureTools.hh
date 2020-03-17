#ifndef SIMPLESTRUCTURETOOLS_HH
#define SIMPLESTRUCTURETOOLS_HH

#include <vector>
#include <string>
#include <set>
#include "casm/external/Eigen/Dense"
#include "casm/global/definitions.hh"
#include "casm/crystallography/DoFDecl.hh"


namespace CASM {


  namespace xtal {

    /** \ingroup Structure
     *  @{
     */

    //TODO: Move this into crystallography declarations or something of the sort
    namespace SimpleStructureTools {
      // Defined in SimpleStructure.hh
      enum class SpeciesMode;
    }

    class SimpleStructure;
    class Site;
    class Molecule;
    class BasicStructure;

    SimpleStructure make_superstructure(Eigen::Ref<const Eigen::Matrix3i> const &_T, SimpleStructure const &_sstruc);

    /// \brief Construct from decorated structure
    SimpleStructure make_simple_structure(BasicStructure const &_struc);

    /// \brief Determine which sites of a BasicStructure can host each atom of a SimpleStructure
    /// result[i] is set of site indices in @param _prim that can host atom 'i' of @param sstruc
    std::vector<std::set<Index> > atom_site_compatibility(SimpleStructure const &sstruc, BasicStructure const &_prim);

    /// \brief Determine which sites of a BasicStructure can host each molecule of a SimpleStructure
    /// result[i] is set of site indices in @param _prim that can host molecule 'i' of @param sstruc
    std::vector<std::set<Index> > mol_site_compatibility(SimpleStructure const &sstruc, BasicStructure const &_prim);

    /// \brief use information in _reference to initialize atom_info from mol_info
    void _atomize(SimpleStructure &_sstruc, Eigen::Ref<const Eigen::VectorXi> const &_mol_occ, BasicStructure const &_reference);

    /// \brief Construct BasicStructure from SimpleStructure.
    /// @param _sstruc SimpleStructure used as source data for conversion
    /// @param _all_dofs holds names of additional DOFs to initialize in structure
    /// @param mode specifies whether ATOM or MOL info of _sstruc should be used to build sites of structure
    /// @param _allowed_occupants List of allowed molecules at each site; if empty, occupants are assumed to be atoms
    ///        having the species names and attributes indicated by _sstruc
    BasicStructure make_basic_structure(SimpleStructure const &_sstruc,
                                        std::vector<DoFKey> const &_all_dofs,
                                        SimpleStructureTools::SpeciesMode mode,
                                        std::vector<std::vector<Molecule> > _allowed_occupants = {});
    /** @} */
  }

}
#endif
