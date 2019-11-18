#ifndef SIMPLESTRUCTURE_HH
#define SIMPLESTRUCTURE_HH

#include <vector>
#include <string>
#include <set>
#include <map>
#include "casm/external/Eigen/Dense"
#include "casm/global/definitions.hh"

namespace CASM {
  namespace xtal {

    /** \ingroup Structure
     *  @{
     */

    /// Used to construct a BasicStructure from a 'properties.calc.json' object
    class SimpleStructure {
    public:
      enum class SpeciesMode {ATOM, MOL};

      struct Info {
        // Can treat as a Eigen::VectorXd
        using Coord = Eigen::MatrixXd::ColXpr;
        using ConstCoord = Eigen::MatrixXd::ConstColXpr;

        std::vector<std::string> names;
        Eigen::MatrixXd coords;

        std::map<std::string, Eigen::MatrixXd> properties;

        std::vector<Index> sort_by_name();

        Coord coord(Index i) {
          return coords.col(i);
        }

        ConstCoord coord(Index i) const {
          return coords.col(i);
        }

        void resize(Index N) {
          names.resize(N, "Va");
          coords.setZero(3, N);
        }

        Index size() const {
          return names.size();
        }
      };

      SimpleStructure(std::string const &_prefix = std::string());

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

      /// \brief Apply homogeneous deformation gradient tensor _F to lat_column_mat, mol_info, and atom_info
      void deform_coords(Eigen::Ref<const Eigen::Matrix3d> const &_F);

      /// \brief Apply homogeneous rotation matrix _R to lat_column_mat, mol_info, and atom_info
      void rotate_coords(Eigen::Ref<const Eigen::Matrix3d> const &_R);

      Eigen::Matrix3d lat_column_mat;

      // Use occupation vector in order to avoid messy molecule-name aliasing issues
      Info mol_info;

      Info atom_info;

      std::map<std::string, Eigen::VectorXd> properties;

    private:
      std::string m_prefix;
    };

    /** @} */
  }
}

#endif
