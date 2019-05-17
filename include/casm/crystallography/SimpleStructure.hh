#ifndef SIMPLESTRUCTURE_HH
#define SIMPLESTRUCTURE_HH

#include <vector>
#include <string>
#include <set>
#include "casm/external/Eigen/Dense"
#include "casm/casm_io/jsonParser.hh"
#include "casm/CASM_global_definitions.hh"



namespace CASM {

  /** \ingroup Structure
   *  @{
   */
  class SimpleStructure;
  class Supercell;
  class ConfigDoF;
  class Configuration;
  class Site;

  template<typename CoordType>
  class BasicStructure;

  namespace DoFType {
    class Traits;
  }

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



  /// Used to construct a BasicStructure from a 'properties.calc.json' object
  class SimpleStructure {
  public:
    using StrucType = BasicStructure<Site>;

    SimpleStructure(std::string const &_prefix = std::string());

    /// \brief Construct from decorated structure and specify prefix for output quantities
    SimpleStructure(BasicStructure<Site> const &_struc, const std::string &_prefix = "");

    /// \brief Construct from Configuration and specify prefix for output quantities
    SimpleStructure(Configuration const &_config, const std::string &_prefix = "");

    /// \brief Construct from ConfigDoF _dof belonging to provided Supercell _scel; specify prefix for output quantities
    SimpleStructure(Supercell const &_scel, ConfigDoF const &_dof, const std::string &_prefix = "");

    /// \brief Apply homogeneous deformation gradient tensor _F to lat_column_mat, mol_info, and atom_info
    void deform(Eigen::Ref<const Eigen::Matrix3d> const &_F);

    /// \brief Apply homogeneous rotation matrix _R to lat_column_mat, mol_info, and atom_info
    void rotate(Eigen::Ref<const Eigen::Matrix3d> const &_R);

    /// \brief use information in _reference to initialize atom_info from mol_info
    void atomize(BasicStructure<Site> const &_reference);

    Index n_mol() const {
      return mol_info.names.size();
    }

    Index n_atom() const {
      return atom_info.names.size();
    }

    std::string const &prefix() const {
      return m_prefix;
    }

    struct Info {
      std::vector<std::string> names;
      Eigen::MatrixXd coords;
      Eigen::MatrixXi SD;
      jsonParser dofs;
      jsonParser props;
      mutable std::vector<Index> permute;

      Index size() const {
        return names.size();
      }
    };

    Eigen::Matrix3d lat_column_mat;
    Eigen::Matrix3d cartesian_isometry;
    bool selective_dynamics;

    // Use occupation vector in order to avoid messy molecule-name aliasing issues
    Eigen::VectorXi mol_occ;
    Info mol_info;

    Info atom_info;

    jsonParser global_dofs;
    jsonParser global_props;
    //std::map<std::string,Tensor<double> > local_dofs;
    //std::map<std::string,Tensor<double> > global_dofs;

  private:

    /// \brief Imposes DoF values from ConfigDoF _config onto *this, using using any necessary information contained in _reference
    void _apply_dofs(ConfigDoF const &_config, BasicStructure<Site> const &_reference);

    std::string m_prefix;
    bool m_atomized;
  };

  std::vector<std::set<Index> > site_compatibility(SimpleStructure const &sstruc, BasicStructure<Site> const &_prim);

  /// \brief Output to JSON, excluding any molecular or atomic species contained in 'excluded_species'
  jsonParser &to_json(SimpleStructure const &_struc, jsonParser &json, std::set<std::string> const &excluded_species = {"Va", "VA", "va"});

  /// \brief Read from JSON
  void from_json(SimpleStructure &_struc, const jsonParser &json);

  std::string POS_file(SimpleStructure &_struc);

  /** @} */
}

#endif
