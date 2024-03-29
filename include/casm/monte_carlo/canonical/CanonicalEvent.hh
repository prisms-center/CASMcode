#ifndef CASM_Monte_CanonicalEvent
#define CASM_Monte_CanonicalEvent

#include "casm/external/Eigen/Dense"
#include "casm/global/definitions.hh"
#include "casm/monte_carlo/MonteDefinitions.hh"
#include "casm/monte_carlo/OccLocation.hh"

namespace CASM {
namespace Monte {

/// \brief Data structure for storing information regarding a proposed grand
/// canonical Monte Carlo event
class CanonicalEvent {
 public:
  /// \brief Default constructor
  CanonicalEvent() {}

  /// \brief Constructor
  ///
  /// \param Nspecies The number of different molecular species in this
  /// calculation (use CompositionConverter::components().size()) \param Ncorr
  /// The total number of correlations that could be calculated (use
  /// Clexulator::corr_size)
  ///
  CanonicalEvent(size_type Nspecies, size_type Ncorr, size_type Neta = 0);

  /// \brief Set the change in (extensive) formation energy associated with this
  /// event
  void set_dEf(double dE);

  /// \brief Return change in (extensive) formation energy associated with this
  /// event
  double dEf() const;

  /// \brief const Access change in number of species per supercell. Zeros, size
  /// of CompositionConverter::components().
  const Eigen::VectorXl &dN() const;

  /// \brief Access change in order parameter (intensive)
  Eigen::VectorXd &deta();

  /// \brief const Access change in order parameter (intensive)
  Eigen::VectorXd const &deta() const;

  /// \brief Return change in number of species in supercell. Zeros, size of
  /// CompositionConverter::components().
  long int dN(size_type species_type_index) const;

  /// \brief Set change in (extensive) potential energy
  void set_dEpot(double dpot_nrg);

  /// \brief Return change in (extensive) potential energy, dEpot = dEf
  double dEpot() const;

  /// \brief Access the changes in (extensive) correlations associated with this
  /// event
  Eigen::VectorXd &dCorr();

  /// \brief const Access the changes in (extensive) correlations associated
  /// with this event
  const Eigen::VectorXd &dCorr() const;

  /// \brief Access the data describing this event
  OccEvent &occ_event();

  /// \brief const Access the data describing this event
  const OccEvent &occ_event() const;

 private:
  /// \brief Change in (extensive) correlations due to this event
  Eigen::VectorXd m_dCorr;

  /// \brief Change in (extensive) formation energy due to this event
  double m_dEf;

  /// \brief Change in (extensive) potential energy due to this event
  double m_dEpot;

  /// \brief Change in number of each species in supercell due to this event.
  ///        Zeros, size of primclex.get_param_comp().get_components()
  Eigen::VectorXl m_dN;

  /// \brief Change in order parameter due to this event.
  Eigen::VectorXd m_deta;

  /// \brief The modifications performed by this event
  OccEvent m_occ_event;
};

}  // namespace Monte
}  // namespace CASM

#endif
