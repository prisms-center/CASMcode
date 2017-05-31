#include <functional>
#include "casm/clex/ConfigIOHull.hh"
#include "casm/clex/Configuration.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/Norm.hh"
#include "casm/database/ConfigDatabase.hh"

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


    // --- BaseHull Template Definitions (only needed here) -------------------

    /// \brief Constructor
    ///
    /// \param _name Formatter name, e.g. "on_hull", "hull_dist", etc.
    /// \param _desc Formatter help description
    /// \param _comp_calculator Function or functor that return the composition of a Configuration as a Eigen::VectorXd
    /// \param _energy_calcuator Function or functor that return the energy of a Configuration as a double
    /// \param _singular_value_tol Used with CASM::almost_zero to detect zero-valued singular values
    ///
    /// - The selection of Configuration used to generate the hull is passed as an argument. See parse_args
    /// - The units of composition and energy are determined by the calculator parameters
    ///
    template<typename ValueType>
    BaseHull<ValueType>::BaseHull(const std::string &_name,
                                  const std::string &_desc,
                                  const std::string &_default_selection,
                                  const std::string &_default_composition_type,
                                  const Hull::CalculatorOptions &_calculator_map,
                                  double _singular_value_tol,
                                  double _bottom_facet_tol) :
      BaseValueFormatter<ValueType, Configuration>(_name, _desc),
      m_selection(_default_selection),
      m_composition_type(_default_composition_type),
      m_calculator_map(_calculator_map),
      m_singular_value_tol(_singular_value_tol),
      m_bottom_facet_tol(_bottom_facet_tol),
      m_initialized(false) {}

    /// \brief Calculates the convex hull
    ///
    /// - Uses the parsed args to determine the selection to use to calculate the hull
    /// - Calls 'init' on the CompCalculator and EnergyCalculator
    /// - Constructs the hull
    ///
    template<typename ValueType>
    void BaseHull<ValueType>::init(const Configuration &_tmplt) const {

      DB::Selection<Configuration> selection(
        _tmplt.primclex().db<Configuration>(),
        m_selection);

      // Hull typedefs:
      //typedef std::pair<notstd::cloneable_ptr<CompCalculator>,
      //                  notstd::cloneable_ptr<EnergyCalculator> > CalculatorPair;
      //typedef std::map<std::string, CalculatorPair> CalculatorOptions;

      auto res = m_calculator_map.find(m_composition_type)->second;
      res.first->init(_tmplt);
      res.second->init(_tmplt);

      m_hull = std::make_shared<Hull>(selection,
                                      *res.first,
                                      *res.second,
                                      m_singular_value_tol,
                                      m_bottom_facet_tol);

    }

    /// \brief column header to use
    ///
    /// \returns 'name() + '(' + args + ')'
    /// - ex: 'on_hull(MASTER)' if name() is 'on_hull', and args is "MASTER"
    /// - ex: 'on_clex_hull(ALL,comp) if name() is 'on_clex_hull', and args is "ALL,comp"
    ///
    template<typename ValueType>
    std::string BaseHull<ValueType>::short_header(const Configuration &_config) const {
      std::stringstream t_ss;
      t_ss << this->name() << "(" << m_selection << "," << m_composition_type << ")";
      return t_ss.str();
    }

    /// \brief Determine the selection to use to generate the hull
    ///
    /// Args are: ($selection,$composition,$dim_tol,$bottom_tol)
    ///
    /// Options for $selection are:
    /// - "ALL", use all configurations
    /// - "MASTER", use the current MASTER selection (default for no args)
    /// - "CALCULATED", use configurations for which is_calculated is true
    /// - other, assume the argument is the filename for a selection to use
    ///
    /// Options for $composition are:
    /// - "atom_frac", (default) use atom_frac for the composition and "formation_energy_per_species" for the energy
    /// - "comp", use parametric composition for the composition and "formation_energy" for the energy
    ///
    /// $dim_tol, default=1e-8
    /// - singular value tolerance used for detecting composition dimensions
    ///
    /// $bottom_tol, default=1e-8
    /// - tolerance used for detecting facets on the convex hull bottom
    ///
    template<typename ValueType>
    bool BaseHull<ValueType>::parse_args(const std::string &args) {

      if(m_initialized) {
        return false;
      }

      std::vector<std::string> splt_vec;
      boost::split(splt_vec, args, boost::is_any_of(","), boost::token_compress_on);

      if(splt_vec.size() > 4) {
        throw std::runtime_error("Attempted to initialize format tag " + this->name()
                                 + " with " + std::to_string(splt_vec.size()) + " arguments ("
                                 + args + "), but only up to 4 arguments allowed.\n");
        return false;
      }

      while(splt_vec.size() < 4) {
        splt_vec.push_back("");
      }

      m_selection = splt_vec[0].empty() ? m_selection : splt_vec[0];
      m_composition_type = splt_vec[1].empty() ? m_composition_type : splt_vec[1];
      m_singular_value_tol = splt_vec[2].empty() ? m_singular_value_tol : std::stod(splt_vec[2]);
      m_bottom_facet_tol = splt_vec[3].empty() ? m_bottom_facet_tol : std::stod(splt_vec[3]);
      m_initialized = true;

      auto it = m_calculator_map.find(m_composition_type);
      if(it == m_calculator_map.end()) {

        std::stringstream ss;
        ss << "Composition type '" << m_composition_type << "' is not recognized. Options are: ";
        for(auto it = m_calculator_map.begin(); it != m_calculator_map.end(); ++it) {
          std::cerr << it->first << " ";
          if(it != m_calculator_map.begin()) {
            ss << ", ";
          }
          ss << "'" << it->first << "'";
        }
        throw std::runtime_error(ss.str());

      }

      // prevent combining formatters
      return false;
    }

    /// \brief const Access the Hull object
    template<typename ValueType>
    const Hull &BaseHull<ValueType>::_hull() const {
      return *m_hull;
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

    template class BaseHull<bool>;
    template class BaseHull<double>;
  }
}

