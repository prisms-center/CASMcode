#include <functional>
#include "casm/clex/ConfigIO.hh"
#include "casm/clex/ConfigIOHull.hh"
#include "casm/clex/Configuration.hh"
#include "casm/clex/PrimClex.hh"

namespace CASM {

  namespace ConfigIO_impl {
    
    // --- OnHullConfigFormatter Definitions -------------------
    
    /// \brief Constructor
    OnHullConfigFormatter::OnHullConfigFormatter() :
      BaseHullConfigFormatter<bool>("on_hull", 
        std::string("Whether configuration is a vertex on the formation_energy convex hull (i.e., is a groundstate).")
                  + " Only one Configuration out of a set that have identical or almost identical points in" 
                  + " composition/energy space will return true."
                  + " Accepts arguments ($selection, $composition)."
                  + " ($selection may be one of: <filename>, 'ALL', 'CALCULATED', 'MASTER' <--default)"
                  + " ($composition may be one of: 'comp', 'atom_frac' <--default)"
                  + " For 'comp', 'formation_energy' is used. For 'atom_frac', 'formation_energy_per_atom' is used."
                  + " Ex: clex_hull_dist, clex_hull_dist(MASTER,comp).",
        "MASTER",
        "atom_frac",
        CASM::TOL) {
      m_calculator_map["atom_frac"] = std::make_pair(CASM::species_frac, CASM::formation_energy_per_species);
      m_calculator_map["comp"] = std::make_pair(CASM::comp, CASM::formation_energy);
    }
    
    /// \brief Validate that the Configuration has a formation energy per species
    bool OnHullConfigFormatter::_validate(const Configuration &_config) const {
      return _config.delta_properties().contains("relaxed_energy");
    }
    
    /// \brief Check if the Configuration is a hull vertex 
    ///
    /// - Only returns true for one Configuration out of a set that have identical or almost 
    ///   identical points in composition/energy space
    bool OnHullConfigFormatter::_evaluate(const Configuration &_config) const {
      orgQhull::QhullVertexList vertices = _hull().data().vertexList();
      for(auto it=vertices.begin(); it!=vertices.end(); ++it) {
        if(_hull().configuration(*it).name() == _config.name()) {
          return true;
        }
      }
      return false;
    }
    
    
    // --- HullDistConfigFormatter Definitions -------------------
    
    /// \brief Constructor
    HullDistConfigFormatter::HullDistConfigFormatter() :
      BaseHullConfigFormatter<double>("hull_dist", 
        std::string("Distance, in eV, of a configuration's formation_energy_per_atom above the convex hull.")
                  + " Accepts arguments ($selection,$composition)."
                  + " ($selection may be one of: <filename>, 'ALL', 'CALCULATED', 'MASTER' <--default)"
                  + " ($composition may be one of: 'comp', 'atom_frac' <--default)."
                  + " For 'comp', 'formation_energy' is used. For 'atom_frac', 'formation_energy_per_atom' is used."
                  + " Ex: clex_hull_dist, clex_hull_dist(MASTER,comp).",
        "MASTER",
        "atom_frac",
        CASM::TOL) {
      m_calculator_map["atom_frac"] = std::make_pair(CASM::species_frac, CASM::formation_energy_per_species);
      m_calculator_map["comp"] = std::make_pair(CASM::comp, CASM::formation_energy);
    }
    
    /// \brief Validate that the Configuration has a formation energy per species
    bool HullDistConfigFormatter::_validate(const Configuration &_config) const {
      return _config.delta_properties().contains("relaxed_energy");
    }
    
    /// \brief Return the distance to the hull
    double HullDistConfigFormatter::_evaluate(const Configuration &_config) const {
      double d = _hull().dist_to_hull(_config);
      d = (std::abs(d) < m_dist_to_hull_tol) ? 0.0 : d;
      return d;
    }
    
    
    // --- OnClexHullConfigFormatter Definitions -------------------
    
    /// \brief Constructor
    OnClexHullConfigFormatter::OnClexHullConfigFormatter() :
      BaseHullConfigFormatter<bool>("on_clex_hull", 
        std::string("Whether configuration is a vertex on the *cluster-expanded* formation_energy convex hull (i.e., is a *predicted* groundstate).")
                  + " Only one Configuration out of a set that have identical or almost identical points in" 
                  + " composition/energy space will return true."
                  + " Accepts arguments ($selection,$composition)."
                  + " ($selection may be one of: <filename>, 'ALL', 'CALCULATED', 'MASTER' <--default)"
                  + " ($composition may be one of: 'comp', 'atom_frac' <--default)"
                  + " For 'comp', 'clex(formation_energy)' is used. For 'atom_frac', 'clex(formation_energy_per_atom)' is used."
                  + " Ex: clex_hull_dist, clex_hull_dist(MASTER,comp).",
        "MASTER",
        "atom_frac",
        CASM::TOL) {
      m_calculator_map["atom_frac"] = std::make_pair(CASM::species_frac, CASM::clex_formation_energy_per_species);
      m_calculator_map["comp"] = std::make_pair(CASM::comp, CASM::clex_formation_energy);
    }
    
    
    /// \brief Validate that the Configuration has a cluster expanded formation energy per species
    ///
    /// - Currently always returns true
    bool OnClexHullConfigFormatter::_validate(const Configuration &_config) const {
      return true;
    }
    
    /// \brief Check if the Configuration is a hull vertex 
    ///
    /// - Only returns true for one Configuration out of a set that have identical or almost 
    ///   identical points in composition/energy space
    bool OnClexHullConfigFormatter::_evaluate(const Configuration &_config) const {
      orgQhull::QhullVertexList vertices = _hull().data().vertexList();
      for(auto it=vertices.begin(); it!=vertices.end(); ++it) {
        if(_hull().configuration(*it).name() == _config.name()) {
          return true;
        }
      }
      return false;
    }
    
    
    // --- ClexHullDistConfigFormatter Definitions -------------------
    
    /// \brief Constructor
    ClexHullDistConfigFormatter::ClexHullDistConfigFormatter() :
      BaseHullConfigFormatter<double>("clex_hull_dist", 
        std::string("Distance, in eV, of a configuration's *cluster-expanded* formation_energy_per_atom above the convex hull.")
                  + " Accepts arguments ($selection,$composition)."
                  + " ($selection may be one of: <filename>, 'ALL', 'CALCULATED', 'MASTER' <--default)"
                  + " ($composition may be one of: 'comp', 'atom_frac' <--default)"
                  + " For 'comp', 'clex(formation_energy)' is used. For 'atom_frac', 'clex(formation_energy_per_atom)' is used."
                  + " Ex: clex_hull_dist, clex_hull_dist(MASTER,comp).",
        "MASTER",
        "atom_frac",
        CASM::TOL) {
      m_calculator_map["atom_frac"] = std::make_pair(CASM::species_frac, CASM::clex_formation_energy_per_species);
      m_calculator_map["comp"] = std::make_pair(CASM::comp, CASM::clex_formation_energy);
    }
    
    /// \brief Validate that the Configuration has a cluster expanded formation energy per species
    ///
    /// - Currently always returns true
    bool ClexHullDistConfigFormatter::_validate(const Configuration &_config) const {
      return true;
    }
    
    /// \brief Return the distance to the hull
    double ClexHullDistConfigFormatter::_evaluate(const Configuration &_config) const {
      double d = _hull().dist_to_hull(_config);
      d = (std::abs(d) < m_dist_to_hull_tol) ? 0.0 : d;
      return d;
    }
    
  }
}

