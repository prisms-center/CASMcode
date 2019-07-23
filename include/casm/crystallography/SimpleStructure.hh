#ifndef SIMPLESTRUCTURE_HH
#define SIMPLESTRUCTURE_HH

#include <vector>
#include <string>
#include <set>
#include <map>
#include "casm/external/Eigen/Dense"
#include "casm/CASM_global_definitions.hh"



namespace CASM {

  /** \ingroup Structure
   *  @{
   */

  /// Used to construct a BasicStructure from a 'properties.calc.json' object
  class SimpleStructure {
  public:
    enum class SpeciesMode {ATOM, MOL};

    struct Info {
      std::vector<std::string> names;
      Eigen::MatrixXd coords;
      Eigen::MatrixXi SD;

      std::map<std::string, Eigen::MatrixXd> dofs;
      std::vector<Index> permute;

      void sort_by_name();

      Index size() const {
        return names.size();
      }
    };

    SimpleStructure(std::string const &_prefix = std::string());

    /// \brief Apply homogeneous deformation gradient tensor _F to lat_column_mat, mol_info, and atom_info
    void deform(Eigen::Ref<const Eigen::Matrix3d> const &_F);

    /// \brief Apply homogeneous rotation matrix _R to lat_column_mat, mol_info, and atom_info
    void rotate(Eigen::Ref<const Eigen::Matrix3d> const &_R);

    Index n_mol() const {
      return mol_info.size();
    }

    Index n_atom() const {
      return atom_info.size();
    }

    Info &info(SpeciesMode _mode) {
      return _mode == SpeciesMode::MOL ? mol_info : atom_info;
    }

    Info const &info(SpeciesMode _mode) const {
      return _mode == SpeciesMode::MOL ? mol_info : atom_info;
    }

    std::string const &prefix() const {
      return m_prefix;
    }

    Eigen::Matrix3d lat_column_mat;
    Eigen::Matrix3d cartesian_isometry;
    bool selective_dynamics;

    // Use occupation vector in order to avoid messy molecule-name aliasing issues
    Info mol_info;

    Info atom_info;

    std::map<std::string, Eigen::VectorXd> dofs;

  private:
    std::string m_prefix;
  };

  /** @} */
}

#endif
