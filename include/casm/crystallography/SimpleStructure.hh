#ifndef SIMPLESTRUCTURE_HH
#define SIMPLESTRUCTURE_HH

#include <map>
#include <set>
#include <string>
#include <vector>

#include "casm/external/Eigen/Dense"
#include "casm/global/definitions.hh"

namespace CASM {
namespace xtal {

/** \ingroup Structure
 *  @{
 */

namespace SimpleStructureTools {
// TODO: Capitalize so it looks the same as COORD_MODE? Put it in some
// crystallography declarations file?
///\brief enum to refer to a particular representation of the occupants (atomic
/// or molecular)
enum class SpeciesMode { ATOM, MOL };
}  // namespace SimpleStructureTools

/// \brief Representation of a crystal of molecular and/or atomic occupants,
/// and any additional properties.
/// It does not require a primitive Structure or BasicStructure object to act
/// as a reference and can be thought of as the class that embodies all the
/// information of a config.json or properties.calc.json file
class SimpleStructure {
 public:
  using SpeciesMode = SimpleStructureTools::SpeciesMode;

  /// \brief Struct to encode all information about the crystal basis
  /// Info may describe the basis in a atomic context or molecular context
  /// Within a particular context, species are considered indivisible points
  /// If both contexts are present within the SimpleStructure, they are assumed
  /// to be different representation of the same crystal basis. The only
  /// difference is that the 'atomic' context separates each molecule into its
  /// constituent atoms. SimpleStructure DOES NOT TRACK how the molecular and
  /// atomic representations are related
  struct Info {
    /// Return types for accessing individual coordinates
    /// Can treat as a Eigen::VectorXd
    using Coord = Eigen::MatrixXd::ColXpr;
    using ConstCoord = Eigen::MatrixXd::ConstColXpr;

    /// \brief names[i] is name of species that occupies sites 'i'
    std::vector<std::string> names;

    /// \brief (3 x names.size()) matrix of coordinates.
    /// coords.col(i) is Cartesian coordinate of site 'i'
    Eigen::MatrixXd coords;

    /// \brief map of [property name, (m x names.size()) matrix] for all
    /// numerical site properties properties are assumed to be vectors of some
    /// property-specific dimension 'm'
    std::map<std::string, Eigen::MatrixXd> properties;

    /// \brief permutation that results in sites sorted alphabetically by
    /// species Guaranteed stable: will not change order for two sites with same
    /// species
    std::vector<Index> sort_by_name();

    /// \brief Access coordinate of site 'i'
    Coord cart_coord(Index i) { return coords.col(i); }

    /// \brief const access for coordinate of site 'i'
    ConstCoord cart_coord(Index i) const { return coords.col(i); }

    /// \brief Resize to hold N sites. All coordinates set to zero, all
    /// occupants set to "Va"
    void resize(Index N) {
      names.resize(N, "Va");
      coords.setZero(3, N);
      for (auto &el : properties) el.second.setZero(Eigen::NoChange, N);
    }

    /// \brief Number of sites is defined as names.size()
    Index size() const { return names.size(); }
  };

  Index n_mol() const { return mol_info.size(); }

  Index n_atom() const { return atom_info.size(); }

  ///\brief Info of specified context (atomic/molecular)
  Info &info(SpeciesMode _mode) {
    return _mode == SpeciesMode::MOL ? mol_info : atom_info;
  }

  ///\brief Const info of specified context (atomic/molecular)
  Info const &info(SpeciesMode _mode) const {
    return _mode == SpeciesMode::MOL ? mol_info : atom_info;
  }

  /// \brief Map all coordinates within the unit cell
  void within(double tol = TOL);

  /// \brief Apply homogeneous deformation gradient tensor _F to lat_column_mat,
  /// mol_info, and atom_info
  void deform_coords(Eigen::Ref<const Eigen::Matrix3d> const &_F);

  /// \brief Apply homogeneous rotation matrix _R to lat_column_mat, mol_info,
  /// and atom_info
  void rotate_coords(Eigen::Ref<const Eigen::Matrix3d> const &_R);

  Eigen::Matrix3d lat_column_mat;

  // Use occupation vector in order to avoid messy molecule-name aliasing issues
  Info mol_info;

  Info atom_info;

  std::map<std::string, Eigen::MatrixXd> properties;
};

/** @} */
}  // namespace xtal
}  // namespace CASM

#endif
