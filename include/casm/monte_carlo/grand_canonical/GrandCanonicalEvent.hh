#ifndef CASM_GrandCanonicalEvent_HH
#define CASM_GrandCanonicalEvent_HH

#include <vector>
#include "casm/external/Eigen/Dense"
#include "casm/CASM_global_definitions.hh"
#include "casm/monte_carlo/DoFMod.hh"

namespace CASM {


  /// \brief Data structure for storing information regarding a proposed grand canonical Monte Carlo event
  class GrandCanonicalEvent {

  public:

    typedef Index size_type;

    /// \brief Default constructor
    GrandCanonicalEvent() {}

    /// \brief Constructor
    ///
    /// \param Nspecies The number of different molecular species in this calculation (use CompositionConverter::components().size())
    /// \param Ncorr The total number of correlations that could be calculated (use Clexulator::corr_size)
    ///
    GrandCanonicalEvent(size_type Nspecies, size_type Ncorr);


    /// \brief Set the change in (extensive) formation energy associated with this event
    void set_dEf(double dE);

    /// \brief Return change in (extensive) formation energy associated with this event
    double dEf() const;


    /// \brief Access change in number of species per supercell. Order as in CompositionConverter::components().
    Eigen::VectorXl &dN();

    /// \brief const Access change in number of species per supercell. Order as in CompositionConverter::components().
    const Eigen::VectorXl &dN() const;

    /// \brief Set the change in number of species in supercell. Order as in CompositionConverter::components().
    void set_dN(size_type species_type_index, long int dn);

    /// \brief Return change in number of species in supercell. Order as in CompositionConverter::components().
    long int dN(size_type species_type_index) const;


    /// \brief Set change in (extensive) potential energy, dEpot = dEf - sum_i(Nunit * param_chem_pot_i * dcomp_x_i)
    void set_dEpot(double dpot_nrg);

    /// \brief Return change in (extensive) potential energy, dEpot = dEf - sum_i(Nunit * param_chem_pot_i * dcomp_x_i)
    double dEpot() const;

    /// \brief Access the changes in (extensive) correlations associated with this event
    Eigen::VectorXd &dCorr();

    /// \brief const Access the changes in (extensive) correlations associated with this event
    const Eigen::VectorXd &dCorr() const;


    /// \brief Access the occupational modification for this event
    OccMod &occupational_change();

    /// \brief const Access the occupational modification for this event
    const OccMod &occupational_change() const;


  private:

    /// \brief Change in (extensive) correlations due to this event
    Eigen::VectorXd m_dCorr;

    /// \brief Change in (extensive) formation energy due to this event
    double m_dEf;

    /// \brief Change in (extensive) potential energy, dEpot = dEf - sum_i(Nunit * param_chem_pot_i * dcomp_x_i)
    double m_dEpot;

    /// \brief Change in number of each species in supercell due to this event.
    ///        The order is determined by primclex.get_param_comp().get_components()
    Eigen::VectorXl m_dN;

    /// \brief The ConfigDoF modification performed by this event
    OccMod m_occ_mod;

  };


  /// \brief Constructor
  ///
  /// \param Nspecies The number of different molecular species in this calculation (use CompositionConverter::components().size())
  /// \param Ncorr The total number of correlations that could be calculated (use Clexulator::corr_size)
  ///
  inline GrandCanonicalEvent::GrandCanonicalEvent(size_type Nspecies, size_type Ncorr) :
    m_dCorr(Eigen::VectorXd(Ncorr)),
    m_dN(Eigen::VectorXl(Nspecies)) { }


  /// \brief Set the change in total (formation) energy associated with this event
  inline void GrandCanonicalEvent::set_dEf(double dEf) {
    m_dEf = dEf;
  }

  /// \brief Return change in total (formation) energy associated with this event
  inline double GrandCanonicalEvent::dEf() const {
    return m_dEf;
  }


  /// \brief Access change in number of all species (extensive). Order as in CompositionConverter::components().
  inline Eigen::VectorXl &GrandCanonicalEvent::dN() {
    return m_dN;
  }

  /// \brief const Access change in number of all species (extensive). Order as in CompositionConverter::components().
  inline const Eigen::VectorXl &GrandCanonicalEvent::dN() const {
    return m_dN;
  }

  /// \brief Set the change in number of species (extensive) described by size_type. Order as in CompositionConverter::components().
  inline void GrandCanonicalEvent::set_dN(size_type species_type_index, long int dNi) {
    m_dN(species_type_index) = dNi;
  }

  /// \brief Return change in number of species (extensive) described by size_type. Order as in CompositionConverter::components().
  inline long int GrandCanonicalEvent::dN(size_type species_type_index) const {
    return m_dN(species_type_index);
  }


  /// \brief Set the change in potential energy: dEpot = dEf - sum_i(Nunit * param_chem_pot_i * dcomp_x_i)
  inline void GrandCanonicalEvent::set_dEpot(double dEpot) {
    m_dEpot = dEpot;
  }

  /// \brief Return change in potential energy: dEpot = dEf - sum_i(Nunit * param_chem_pot_i * dcomp_x_i)
  inline double GrandCanonicalEvent::dEpot() const {
    return m_dEpot;
  }

  /// \brief Access the changes in correlations associated with this event
  inline Eigen::VectorXd &GrandCanonicalEvent::dCorr() {
    return m_dCorr;
  }

  /// \brief const Access the changes in correlations associated with this event
  inline const Eigen::VectorXd &GrandCanonicalEvent::dCorr() const {
    return m_dCorr;
  }


  /// \brief Access the occupational modification for this event
  inline OccMod &GrandCanonicalEvent::occupational_change() {
    return m_occ_mod;
  }

  /// \brief const Access the occupational modification for this event
  inline const OccMod &GrandCanonicalEvent::occupational_change() const {
    return m_occ_mod;
  }

}

#endif
