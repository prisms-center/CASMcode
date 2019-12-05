#ifndef SIMPLESTRUCTURETOOLS_HH
#define SIMPLESTRUCTURETOOLS_HH

#include <vector>
#include <string>
#include <set>
#include "casm/external/Eigen/Dense"
#include "casm/global/definitions.hh"
#include "casm/basis_set/DoFDecl.hh"


namespace CASM {
  class jsonParser;
  class Supercell;
  class ConfigDoF;
  class Configuration;

  namespace DoFType {
    class Traits;
  }


  namespace xtal {

    /** \ingroup Structure
     *  @{
     */
    class SimpleStructure;

    class Site;

    template<typename CoordType>
    class BasicStructure;

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
      void transform(ConfigDoF const  &_config, BasicStructure<Site> const &_reference, SimpleStructure &_struc) const;

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


    SimpleStructure make_superstructure(Eigen::Ref<const Eigen::Matrix3i> const &_T, SimpleStructure const &_sstruc);

    /// \brief Construct from decorated structure
    SimpleStructure make_simple_structure(BasicStructure<Site> const &_struc);

    /// \brief Construct from Configuration
    SimpleStructure make_simple_structure(Configuration const &_config,
                                          std::vector<DoFKey> const &_which_dofs = {});

    /// \brief Construct from ConfigDoF _dof belonging to provided Supercell _scel
    SimpleStructure make_simple_structure(Supercell const &_scel,
                                          ConfigDoF const &_dof,
                                          std::vector<DoFKey> const &_which_dofs = {});


    std::vector<std::set<Index> > atom_site_compatibility(SimpleStructure const &sstruc, BasicStructure<Site> const &_prim);
    std::vector<std::set<Index> > mol_site_compatibility(SimpleStructure const &sstruc, BasicStructure<Site> const &_prim);

    std::vector<std::set<Index> > atom_site_compatibility(SimpleStructure const &sstruc, Configuration const &_config);
    std::vector<std::set<Index> > mol_site_compatibility(SimpleStructure const &sstruc, Configuration const &_config);

    /// \brief Imposes DoF values from ConfigDoF _config onto *this, using using any necessary information contained in _reference
    void _apply_dofs(SimpleStructure &_sstruc, ConfigDoF const &_config, BasicStructure<Site> const &_reference, std::vector<DoFKey> which_dofs);

    /// \brief use information in _reference to initialize atom_info from mol_info
    void _atomize(SimpleStructure &_sstruc, Eigen::Ref<const Eigen::VectorXi> const &_mol_occ, BasicStructure<Site> const &_reference);

    /// \brief Output to JSON, excluding any molecular or atomic species contained in 'excluded_species'
    jsonParser &to_json(xtal::SimpleStructure const &_struc,
                        jsonParser &json_supplement,
                        std::set<std::string> const &excluded_species = {"Va", "VA", "va"},
                        std::string prefix = "");

    /// \brief Read from JSON
    void from_json(xtal::SimpleStructure &_struc, const jsonParser &json, std::string prefix = "");
  }

}
#endif
