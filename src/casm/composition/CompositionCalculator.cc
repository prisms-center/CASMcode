#include "casm/composition/CompositionCalculator.hh"

#include "casm/crystallography/Molecule.hh"

namespace CASM {

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
CompositionCalculator::CompositionCalculator(
    xtal::UnitCellCoordIndexConverter const &_unitcellcoord_index_converter,
    std::vector<std::string> const &_components,
    std::vector<std::vector<Index>> const &_occ_to_component_index_converter)
    : m_unitcellcoord_index_converter(_unitcellcoord_index_converter),
      m_components(_components),
      m_occ_to_component_index_converter(_occ_to_component_index_converter),
      m_volume(m_unitcellcoord_index_converter.total_sites() /
               m_occ_to_component_index_converter.size()),
      m_vacancy_allowed(false) {
  for (int i = 0; i < m_components.size(); ++i) {
    if (xtal::is_vacancy(m_components[i])) {
      m_vacancy_allowed = true;
      m_vacancy_index = i;
    }
  }
}

/// \brief The order of components in composition vector results
std::vector<std::string> CompositionCalculator::components() const {
  return m_components;
}

/// \brief Returns the composition as number per primitive cell, in the order
///     of components
Eigen::VectorXd CompositionCalculator::mean_num_each_component(
    Eigen::VectorXi const &occupation) const {
  return num_each_component(occupation).cast<double>() / m_volume;
}

/// \brief Returns the composition as total number, in the order of components
Eigen::VectorXi CompositionCalculator::num_each_component(
    Eigen::VectorXi const &occupation) const {
  Index n_sites = m_unitcellcoord_index_converter.total_sites();
  Index n_component = m_components.size();

  // initialize
  Eigen::VectorXi result = Eigen::VectorXi::Zero(n_component);

  // count the number of each component
  for (Index i = 0; i < n_sites; i++) {
    Index sublat = m_unitcellcoord_index_converter(i).sublattice();
    Index component_index =
        m_occ_to_component_index_converter[sublat][occupation[i]];
    result[component_index] += 1;
  }

  return result;
}

/// \brief Returns the composition as species fraction, with [Va] = 0.0, in
///        the order of components
Eigen::VectorXd CompositionCalculator::species_frac(
    Eigen::VectorXi const &occupation) const {
  Eigen::VectorXd result = this->mean_num_each_component(occupation);
  if (m_vacancy_allowed) {
    result(m_vacancy_index) = 0.0;
  }
  double sum = result.sum();
  result /= sum;
  return result;
}

}  // namespace CASM
