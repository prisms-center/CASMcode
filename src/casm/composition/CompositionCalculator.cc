#include "casm/composition/CompositionCalculator.hh"

namespace CASM {

/// \brief Constructor
///
/// \param _components Names of components corresponding to each position in
///     the result. Must be consistent with occ_to_component_index_converter.
///     Vacancy components are detected using `xtal::is_vacancy`.
/// \param _occ_to_component_index_converter Lookup table for sublattice index
///     and occupant index to component index. May be obtained from
///     `xtal::make_index_converter`.
/// \param _vacancy_names Set of component names that should be recognized as
///     vacancies. An error is throw if more than one component is a vacancy.
CompositionCalculator::CompositionCalculator(
    std::vector<std::string> const &_components,
    std::vector<std::vector<Index>> const &_occ_to_component_index_converter,
    std::set<std::string> const &_vacancy_names)
    : m_components(_components),
      m_occ_to_component_index_converter(_occ_to_component_index_converter),
      m_vacancy_allowed(false) {
  for (int i = 0; i < m_components.size(); ++i) {
    if (_vacancy_names.count(m_components[i])) {
      if (m_vacancy_allowed) {
        throw std::runtime_error(
            "Error in CompositionCalculator: components contains multiple "
            "vacancy species");
      }
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
  Index n_sites = occupation.size();
  Index n_sublat = m_occ_to_component_index_converter.size();
  Index volume = n_sites / n_sublat;
  return num_each_component(occupation).cast<double>() / volume;
}

/// \brief Returns the composition as total number, in the order of components
Eigen::VectorXi CompositionCalculator::num_each_component(
    Eigen::VectorXi const &occupation) const {
  Index n_sites = occupation.size();
  Index n_sublat = m_occ_to_component_index_converter.size();
  Index volume = n_sites / n_sublat;
  Index n_component = m_components.size();

  // initialize
  Eigen::VectorXi result = Eigen::VectorXi::Zero(n_component);

  // count the number of each component
  for (Index i = 0; i < n_sites; i++) {
    // assumes occuaption vector organized in sublattice blocks
    Index sublat = i / volume;
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
