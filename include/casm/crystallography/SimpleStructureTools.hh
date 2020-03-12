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

  }

}
#endif
