#include "casm/monte_carlo/canonical/CanonicalEvent.hh"

namespace CASM {
namespace Monte {

/// \brief Constructor
///
/// \param Nspecies The number of different molecular species in this
///     calculation (use CompositionConverter::components().size())
/// \param Ncorr The total number of correlations that could be
///     calculated (use Clexulator::corr_size)
///
CanonicalEvent::CanonicalEvent(size_type Nspecies, size_type Ncorr,
                               size_type Neta)
    : m_dCorr(Eigen::VectorXd(Ncorr)),
      m_dN(Eigen::VectorXl::Zero(Nspecies)),
      m_deta(Eigen::VectorXd(Neta)) {}

/// \brief Set the change in total (formation) energy associated with this event
void CanonicalEvent::set_dEf(double dEf) { m_dEf = dEf; }

/// \brief Return change in total (formation) energy associated with this event
double CanonicalEvent::dEf() const { return m_dEf; }

/// \brief const Access change in number of all species (extensive). Order as in
/// CompositionConverter::components().
const Eigen::VectorXl &CanonicalEvent::dN() const { return m_dN; }

/// \brief Return change in number of species (extensive) described by
/// size_type. Order as in CompositionConverter::components().
long int CanonicalEvent::dN(size_type species_type_index) const {
  return m_dN(species_type_index);
}

/// \brief Access change in order parameter
Eigen::VectorXd &CanonicalEvent::deta() { return m_deta; }

/// \brief const Access change in order parameter
Eigen::VectorXd const &CanonicalEvent::deta() const { return m_deta; }

/// \brief Set change in (extensive) potential energy
void CanonicalEvent::set_dEpot(double dpot_nrg) { m_dEpot = dpot_nrg; }

/// \brief Return change in potential energy: dEpot = dEf
double CanonicalEvent::dEpot() const { return m_dEpot; }

/// \brief Access the changes in correlations associated with this event
Eigen::VectorXd &CanonicalEvent::dCorr() { return m_dCorr; }

/// \brief const Access the changes in correlations associated with this event
const Eigen::VectorXd &CanonicalEvent::dCorr() const { return m_dCorr; }

/// \brief Access the data describing this event
OccEvent &CanonicalEvent::occ_event() { return m_occ_event; }

/// \brief const Access the data describing this event
const OccEvent &CanonicalEvent::occ_event() const { return m_occ_event; }

}  // namespace Monte
}  // namespace CASM
