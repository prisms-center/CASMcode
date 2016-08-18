#include <functional>
#include "casm/clex/ConfigIOHull.hh"
#include "casm/clex/Configuration.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/Norm.hh"

namespace CASM {

  namespace ConfigIO {

    /// Returns a map containing default hull calculators
    Hull::CalculatorOptions hull_calculator_options() {
      Hull::CalculatorOptions options;
      options["atom_frac"] = Hull::CalculatorPair(SpeciesFrac().clone(), formation_energy_per_species().clone());
      options["comp"] = Hull::CalculatorPair(Comp().clone(), formation_energy().clone());
      return options;
    }

    /// Returns a map containing default clex hull calculators
    Hull::CalculatorOptions clex_hull_calculator_options() {
      Hull::CalculatorOptions options;
      Clex clex;
      options["comp"] = Hull::CalculatorPair(Comp().clone(), clex.clone());
      clex.parse_args("formation_energy,per_species");
      options["atom_frac"] = Hull::CalculatorPair(SpeciesFrac().clone(), clex.clone());
      return options;
    }


    /// Returns "atom_frac"
    std::string default_hull_calculator() {
      return "atom_frac";
    }

    // --- OnHull Definitions -------------------

    const std::string OnHull::Name = "on_hull";

    const std::string OnHull::Desc =
      "Whether configuration is a vertex on the formation_energy convex hull (i.e., is a groundstate)."
      " Only one Configuration out of a set that have identical or almost identical points in"
      " composition/energy space will return true."
      " Accepts arguments ($selection,$composition,$dim_tol,$bottom_tol)."
      " ($selection may be one of: <filename>, 'ALL', 'CALCULATED', 'MASTER' <--default)"
      " ($composition may be one of: 'comp', 'atom_frac' <--default)"
      " ($dim_tol: tolerance for detecting composition dimensionality, default=1e-8)"
      " ($bottom_tol: tolerance for detecting which facets form the convex hull bottom, default=1e-8)"
      " For 'comp', 'formation_energy' is used. For 'atom_frac', 'formation_energy_per_atom' is used."
      " Ex: on_hull, on_hull(MASTER,comp).";

    /// \brief Constructor
    OnHull::OnHull() :
      BaseHull<bool>(Name, Desc) {}

    /// \brief Validate that the Configuration has a formation energy per species
    bool OnHull::validate(const Configuration &_config) const {
      return has_formation_energy(_config);
    }

    /// \brief Check if the Configuration is a hull vertex
    ///
    /// - Only returns true for one Configuration out of a set that have identical or almost
    ///   identical points in composition/energy space
    bool OnHull::evaluate(const Configuration &_config) const {
      orgQhull::QhullVertexList vertices = _hull().data().vertexList();
      for(auto it = vertices.begin(); it != vertices.end(); ++it) {
        if(_hull().configuration(*it).name() == _config.name()) {
          return true;
        }
      }
      return false;
    }


    // --- HullDist Definitions -------------------

    const std::string HullDist::Name = "hull_dist";

    const std::string HullDist::Desc =
      "Distance, in eV, of a configuration's formation_energy_per_atom above the convex hull."
      " Accepts arguments ($selection,$composition,$dim_tol,$bottom_tol)."
      " ($selection may be one of: <filename>, 'ALL', 'CALCULATED', 'MASTER' <--default)"
      " ($composition may be one of: 'comp', 'atom_frac' <--default)."
      " ($dim_tol: tolerance for detecting composition dimensionality, default=1e-8)"
      " ($bottom_tol: tolerance for detecting which facets form the convex hull bottom, default=1e-8)"
      " For 'comp', 'formation_energy' is used. For 'atom_frac', 'formation_energy_per_atom' is used."
      " Ex: hull_dist, hull_dist(MASTER,comp).";

    /// \brief Constructor
    HullDist::HullDist() :
      BaseHull<double>(Name, Desc) {}

    /// \brief Validate that the Configuration has a formation energy per species
    bool HullDist::validate(const Configuration &_config) const {
      return has_formation_energy(_config);
    }

    /// \brief Return the distance to the hull
    double HullDist::evaluate(const Configuration &_config) const {
      double d = _hull().dist_to_hull(_config);
      d = (std::abs(d) < m_dist_to_hull_tol) ? 0.0 : d;
      return d;
    }


    // --- OnClexHull Definitions -------------------

    const std::string OnClexHull::Name = "on_clex_hull";

    const std::string OnClexHull::Desc =
      "Whether configuration is a vertex on the *cluster-expanded* formation_energy "
      "convex hull (i.e., is a *predicted* groundstate)."
      " Only one Configuration out of a set that have identical or almost identical points in"
      " composition/energy space will return true."
      " Accepts arguments ($selection,$composition,$dim_tol,$bottom_tol)."
      " ($selection may be one of: <filename>, 'ALL', 'CALCULATED', 'MASTER' <--default)"
      " ($composition may be one of: 'comp', 'atom_frac' <--default)"
      " ($dim_tol: tolerance for detecting composition dimensionality, default=1e-8)"
      " ($bottom_tol: tolerance for detecting which facets form the convex hull bottom, default=1e-8)"
      " For 'comp', 'clex(formation_energy)' is used. For 'atom_frac', 'clex(formation_energy_per_atom)' is used."
      " Ex: clex_hull_dist, clex_hull_dist(MASTER,comp).";

    /// \brief Constructor
    OnClexHull::OnClexHull() :
      BaseHull<bool>(Name, Desc, "MASTER", default_hull_calculator(), clex_hull_calculator_options()) {}


    /// \brief Validate that the Configuration has a cluster expanded formation energy per species
    ///
    /// - Currently always returns true
    bool OnClexHull::validate(const Configuration &_config) const {
      return true;
    }

    /// \brief Check if the Configuration is a hull vertex
    ///
    /// - Only returns true for one Configuration out of a set that have identical or almost
    ///   identical points in composition/energy space
    bool OnClexHull::evaluate(const Configuration &_config) const {
      orgQhull::QhullVertexList vertices = _hull().data().vertexList();
      for(auto it = vertices.begin(); it != vertices.end(); ++it) {
        if(_hull().configuration(*it).name() == _config.name()) {
          return true;
        }
      }
      return false;
    }


    // --- ClexHullDist Definitions -------------------

    const std::string ClexHullDist::Name = "clex_hull_dist";

    const std::string ClexHullDist::Desc =
      "Distance, in eV, of a configuration's *cluster-expanded* "
      "formation_energy_per_atom above the convex hull."
      " Accepts arguments ($selection,$composition,$dim_tol,$bottom_tol)."
      " ($selection may be one of: <filename>, 'ALL', 'CALCULATED', 'MASTER' <--default)"
      " ($composition may be one of: 'comp', 'atom_frac' <--default)"
      " ($dim_tol: tolerance for detecting composition dimensionality, default=1e-8)"
      " ($bottom_tol: tolerance for detecting which facets form the convex hull bottom, default=1e-8)"
      " For 'comp', 'clex(formation_energy)' is used. For 'atom_frac', 'clex(formation_energy_per_atom)' is used."
      " Ex: clex_hull_dist, clex_hull_dist(MASTER,comp).";

    /// \brief Constructor
    ClexHullDist::ClexHullDist() :
      BaseHull<double>(Name, Desc, "MASTER", default_hull_calculator(), clex_hull_calculator_options()) {}

    /// \brief Validate that the Configuration has a cluster expanded formation energy per species
    ///
    /// - Currently always returns true
    bool ClexHullDist::validate(const Configuration &_config) const {
      return true;
    }

    /// \brief Return the distance to the hull
    double ClexHullDist::evaluate(const Configuration &_config) const {
      double d = _hull().dist_to_hull(_config);
      d = (std::abs(d) < m_dist_to_hull_tol) ? 0.0 : d;
      return d;
    }

  }
}

