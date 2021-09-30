#ifndef CASM_CompositionCalculator
#define CASM_CompositionCalculator

#include <string>
#include <vector>

#include "casm/crystallography/LinearIndexConverter.hh"
#include "casm/global/eigen.hh"

namespace CASM {

/// \brief Calculate composition from occupation vectors
class CompositionCalculator {
 public:
  /// \brief Constructor
  ///
  /// \param unitcellcoord_index_converter Converts between site index (index
  ///     into occupation vector) and UnitCellCoord
  /// \param components Names of components corresponding to each position in
  ///     the result. Must be consistent with occ_to_component_index_converter.
  ///     Vacancy components are detected using `xtal::is_vacancy`.
  /// \param occ_to_component_index_converter Lookup table for sublattice index
  ///     and occupant index to component index. May be obtained from
  ///     `xtal::make_index_converter`.
  CompositionCalculator(
      xtal::UnitCellCoordIndexConverter const &_unitcellcoord_index_converter,
      std::vector<std::string> const &_components,
      std::vector<std::vector<Index>> const &_occ_to_component_index_converter);

  /// \brief The order of components in composition vector results
  std::vector<std::string> components() const;

  /// \brief Returns the composition as number per primitive cell, in the order
  ///     of components
  Eigen::VectorXd mean_num_each_component(
      Eigen::VectorXi const &occupation) const;

  /// \brief Returns the composition as total number, in the order of components
  Eigen::VectorXi num_each_component(Eigen::VectorXi const &occupation) const;

  /// \brief Returns the composition as species fraction, with [Va] = 0.0, in
  ///        the order of components
  Eigen::VectorXd species_frac(Eigen::VectorXi const &occupation) const;

 private:
  // Converts between site index (index into occupation vector) and
  // UnitCellCoord
  xtal::UnitCellCoordIndexConverter m_unitcellcoord_index_converter;

  // Names of components corresponding to each position in
  // the result. Must be consistent with occ_to_component_index_converter.
  // Vacancy components are detected using `xtal::is_vacancy`.
  std::vector<std::string> m_components;

  // Lookup table for sublattice index
  // and occupant index to component index. May be obtained from
  // `xtal::make_index_converter`.
  std::vector<std::vector<Index>> m_occ_to_component_index_converter;

  // Integer volume (number of prim per supercell), cast as double to use for
  // normalization
  double m_volume;

  // Whether one of the components is a vacancy, as detected by
  // `xtal::is_vacancy`
  bool m_vacancy_allowed;

  // Index in components of vacancies
  Index m_vacancy_index;
};

}  // namespace CASM

#endif
