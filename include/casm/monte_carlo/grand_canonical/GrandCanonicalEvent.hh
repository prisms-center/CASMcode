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

    
    /// \brief Set the change in (intensive) formation energy associated with this event
    void set_dformation_energy(double dE);
    
    /// \brief Return change in (intensive) formation energy associated with this event
    double dformation_energy() const;

    
    /// \brief Access change in number of species per primitive cell. Order as in CompositionConverter::components().
    Eigen::VectorXd& dcomp_n();
    
    /// \brief const Access change in number of species per primitive cell. Order as in CompositionConverter::components().
    const Eigen::VectorXd& dcomp_n() const;
    
    /// \brief Set the change in number of species per primitive cell described by size_type. Order as in CompositionConverter::components().
    void set_dcomp_n(size_type species_type_index, double dn);
    
    /// \brief Return change in number of species per primitive cell described by size_type. Order as in CompositionConverter::components().
    double dcomp_n(size_type species_type_index) const;

    
    /// \brief Set change in (intensive) potential energy, dpotential_energy = dformation_energy - sum_i(param_chem_pot_i * dcomp_x_i)
    void set_dpotential_energy(double dpot_nrg);
    
    /// \brief Return change in (intensive) potential energy, dpotential_energy = dformation_energy - sum_i(param_chem_pot_i * dcomp_x_i)
    double dpotential_energy() const;

    /// \brief Access the changes in (intensive) correlations associated with this event
    Eigen::VectorXd& dcorr();
    
    /// \brief const Access the changes in (intensive) correlations associated with this event
    const Eigen::VectorXd& dcorr() const;

    
    /// \brief Access the occupational modification for this event
    OccMod& occupational_change();
    
    /// \brief const Access the occupational modification for this event
    const OccMod& occupational_change() const;

  
  private:
    
    /// \brief Change in (intensive) correlations due to this event
    Eigen::VectorXd m_dcorr;

    /// \brief Change in (intensive) formation energy due to this event
    double m_dformation_energy;
    
    /// \brief Change in (intensive) potential energy, dpotential_energy = dformation_energy - sum_i(param_chem_pot_i * dcomp_x_i) 
    double m_dpotential_energy;
    
    /// \brief Change in number of each species per primitive cell due to this event. 
    ///        The order is determined by primclex.get_param_comp().get_components()
    Eigen::VectorXd m_dcomp_n;
    
    /// \brief The ConfigDoF modification performed by this event
    OccMod m_occ_mod;
    
  };
  
  
  /// \brief Constructor
  ///
  /// \param Nspecies The number of different molecular species in this calculation (use CompositionConverter::components().size())
  /// \param Ncorr The total number of correlations that could be calculated (use Clexulator::corr_size)
  ///
  inline GrandCanonicalEvent::GrandCanonicalEvent(size_type Nspecies, size_type Ncorr) :
    m_dcomp_n(Eigen::VectorXd(Nspecies)),
    m_dcorr(Eigen::VectorXd(Ncorr)) { }

  
  /// \brief Set the change in total (formation) energy associated with this event
  inline void GrandCanonicalEvent::set_dformation_energy(double dE) {
    m_dformation_energy = dE;
  }
  
  /// \brief Return change in total (formation) energy associated with this event
  inline double GrandCanonicalEvent::dformation_energy() const {
    return m_dformation_energy;
  }

  
  /// \brief Access change in number of all species (extensive). Order as in CompositionConverter::components().
  inline Eigen::VectorXd& GrandCanonicalEvent::dcomp_n() {
    return m_dcomp_n;
  }
  
  /// \brief const Access change in number of all species (extensive). Order as in CompositionConverter::components().
  inline const Eigen::VectorXd& GrandCanonicalEvent::dcomp_n() const {
    return m_dcomp_n;
  }
  
  /// \brief Set the change in number of species (extensive) described by size_type. Order as in CompositionConverter::components().
  inline void GrandCanonicalEvent::set_dcomp_n(size_type species_type_index, double dn) {
    m_dcomp_n(species_type_index) = dn;
  }
  
  /// \brief Return change in number of species (extensive) described by size_type. Order as in CompositionConverter::components().
  inline double GrandCanonicalEvent::dcomp_n(size_type species_type_index) const {
    return m_dcomp_n(species_type_index);
  }

  
  /// \brief Set the change in potential energy: dpot_nrg = dform_nrg - sum_i(param_chem_pot_i * dcomp_x_i)
  inline void GrandCanonicalEvent::set_dpotential_energy(double dpot_nrg) {
    m_dpotential_energy = dpot_nrg;
  }
  
  /// \brief Return change in potential energy: dpot_nrg = dform_nrg - sum_i(param_chem_pot_i * dcomp_x_i)
  inline double GrandCanonicalEvent::dpotential_energy() const {
    return m_dpotential_energy;
  }

  /// \brief Access the changes in correlations associated with this event
  inline Eigen::VectorXd& GrandCanonicalEvent::dcorr() {
    return m_dcorr;
  }
  
  /// \brief const Access the changes in correlations associated with this event
  inline const Eigen::VectorXd& GrandCanonicalEvent::dcorr() const {
    return m_dcorr;
  }

  
  /// \brief Access the occupational modification for this event
  inline OccMod& GrandCanonicalEvent::occupational_change() {
    return m_occ_mod;
  }
  
  /// \brief const Access the occupational modification for this event
  inline const OccMod& GrandCanonicalEvent::occupational_change() const {
    return m_occ_mod;
  }
  
}

#endif
