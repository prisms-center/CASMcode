#ifndef CONFIGIOHULL_HH
#define CONFIGIOHULL_HH

#include "casm/casm_io/DataFormatter.hh"
#include "casm/hull/Hull.hh"
#include "casm/clex/ConfigIterator.hh"
#include "casm/clex/PrimClex.hh"

namespace CASM {

  class Configuration;

  namespace ConfigIO {

    /// Returns a map containing default hull calculators
    Hull::CalculatorOptions hull_calculator_options();

    /// Returns a map containing default clex hull calculators
    Hull::CalculatorOptions clex_hull_calculator_options();

    /// Returns "atom_frac"
    std::string default_hull_calculator();


    /// \brief Base class for hull info formatters
    ///
    /// \ingroup ConfigIO
    ///
    template<typename ValueType>
    class BaseHull: public BaseValueFormatter<ValueType, Configuration> {

    public:

      /// \brief Constructor
      BaseHull(const std::string &_name,
               const std::string &_desc,
               const std::string &_default_selection = "MASTER",
               const std::string &_default_composition_type = default_hull_calculator(),
               const Hull::CalculatorOptions &_calculator_map = hull_calculator_options(),
               double _singular_value_tol = 1e-8,
               double _bottom_facet_tol = 1e-8);

      /// \brief Calculates the convex hull
      void init(const Configuration &_tmplt) const override;

      /// \brief column header to use
      std::string short_header(const Configuration &_config) const override;

      /// \brief Determine the selection to use to generate the hull
      bool parse_args(const std::string &args) override;


    protected:

      /// \brief const Access the Hull object
      const Hull &_hull() const;

      // for parse_args: determine energy type based on composition type
      Hull::CalculatorOptions m_calculator_map;

      // detect and avoid printing -0.0
      constexpr static double m_dist_to_hull_tol = 1e-14;

    private:

      // tol used to detect zero singular values during principal component analysis preprocessing before hull calculation
      double m_singular_value_tol;

      // tol used to detect which facets are on the bottom of the convex hull
      double m_bottom_facet_tol;

      // the hull object
      mutable std::shared_ptr<Hull> m_hull;

      // Parsed arguments
      //  -- what selection to use for constructing the hull
      std::string m_selection;

      // Parsed arguments
      //  -- what composition to use for constructing the hull
      //  -- determines comp calculator and energy calculator, via m_calculator_map
      std::string m_composition_type;

      // Prevent re-parsing args
      bool m_initialized;

    };

    /// \brief Returns a boolean indicating if a Configuration is a convex hull vertex
    ///
    /// Whether configuration is a vertex on the formation_energy convex hull (i.e., is a groundstate).
    /// Only one Configuration out of a set that have identical or almost identical points in
    /// composition/energy space will return true.
    ///
    /// Accepts arguments ($selection, $composition):
    /// - $selection may be one of: <filename>, 'ALL', 'CALCULATED', 'MASTER' <--default
    /// - $composition may be one of: 'comp', 'atom_frac' <--default
    ///   - For 'comp', 'formation_energy' is used. For 'atom_frac', 'formation_energy_per_atom' is used.
    ///
    /// Ex: 'on_hull', 'on_hull(MASTER,comp)'
    ///
    /// \ingroup ConfigIO
    ///
    class OnHull: public BaseHull<bool> {

    public:

      static const std::string Name;
      static const std::string Desc;


      /// \brief Constructor
      OnHull();

      /// \brief Clone
      std::unique_ptr<OnHull> clone() const {
        return std::unique_ptr<OnHull>(this->_clone());
      }

      /// \brief Validate that the Configuration has a formation energy per species
      bool validate(const Configuration &_config) const override;

      /// \brief Check if the Configuration is a hull vertex
      bool evaluate(const Configuration &_config) const override;

    private:

      /// \brief Clone
      OnHull *_clone() const override {
        return new OnHull(*this);
      }
    };

    /// \brief Returns the distance, in eV, of a configuration's formation_energy_per_atom above the convex hull
    ///
    /// Accepts arguments ($selection,$composition):
    /// - $selection may be one of: <filename>, 'ALL', 'CALCULATED', 'MASTER' <--default
    /// - $composition may be one of: 'comp', 'atom_frac' <--default
    ///   - For 'comp', 'formation_energy' is used. For 'atom_frac', 'formation_energy_per_atom' is used.
    ///
    /// Ex: 'hull_dist', 'hull_dist(MASTER,comp)'
    ///
    /// \ingroup ConfigIO
    ///
    class HullDist: public BaseHull<double> {

    public:

      static const std::string Name;
      static const std::string Desc;


      /// \brief Constructor
      HullDist();

      /// \brief Clone
      std::unique_ptr<HullDist> clone() const {
        return std::unique_ptr<HullDist>(this->_clone());
      }

      /// \brief Validate that the Configuration has a formation energy per species
      bool validate(const Configuration &_config) const override;

      /// \brief Return the distance to the hull
      double evaluate(const Configuration &_config) const override;

    private:

      /// \brief Clone
      HullDist *_clone() const override {
        return new HullDist(*this);
      }
    };


    /// \brief Returns a boolean indicating if a Configuration is a predicted convex hull vertex
    ///
    /// Whether configuration is a vertex on the *cluster-expanded* formation_energy convex hull (i.e., is a *predicted* groundstate).
    /// Only one Configuration out of a set that have identical or almost identical points in
    /// composition/energy space will return true.
    ///
    /// Accepts arguments ($selection, $composition):
    /// - $selection may be one of: <filename>, 'ALL', 'CALCULATED', 'MASTER' <--default
    /// - $composition may be one of: 'comp', 'atom_frac' <--default
    ///   - For 'comp', 'clex(formation_energy)' is used. For 'atom_frac', 'clex(formation_energy_per_atom)' is used.
    ///
    /// Ex: 'on_clex_hull', 'on_clex_hull(MASTER,comp)'
    ///
    /// \ingroup ConfigIO
    ///
    class OnClexHull: public BaseHull<bool> {

    public:

      static const std::string Name;
      static const std::string Desc;


      /// \brief Constructor
      OnClexHull();

      /// \brief Clone
      std::unique_ptr<OnClexHull> clone() const {
        return std::unique_ptr<OnClexHull>(this->_clone());
      }

      /// \brief Validate that the Configuration has a cluster expanded formation energy per species
      virtual bool validate(const Configuration &_config) const override;

      /// \brief Check if the Configuration is a hull vertex
      virtual bool evaluate(const Configuration &_config) const override;

    private:

      /// \brief Clone
      OnClexHull *_clone() const override {
        return new OnClexHull(*this);
      }

    };

    /// \brief Returns the distance, in eV, of a configuration's clex(formation_energy_per_atom) above the predicted convex hull
    ///
    /// Accepts arguments ($selection,$composition):
    /// - $selection may be one of: <filename>, 'ALL', 'CALCULATED', 'MASTER' <--default
    /// - $composition may be one of: 'comp', 'atom_frac' <--default
    ///   - For 'comp', 'clex(formation_energy)' is used. For 'atom_frac', 'clex(formation_energy_per_atom)' is used.
    ///
    /// Ex: 'clex_hull_dist', 'clex_hull_dist(MASTER,comp)'
    ///
    /// \ingroup ConfigIO
    ///
    class ClexHullDist: public BaseHull<double> {

    public:

      static const std::string Name;
      static const std::string Desc;


      /// \brief Constructor
      ClexHullDist();

      /// \brief Clone
      std::unique_ptr<ClexHullDist> clone() const {
        return std::unique_ptr<ClexHullDist>(this->_clone());
      }

      /// \brief Validate that the Configuration has a cluster expanded formation energy per species
      virtual bool validate(const Configuration &_config) const override;

      /// \brief Return the distance to the hull
      virtual double evaluate(const Configuration &_config) const override;

    private:

      /// \brief Clone
      ClexHullDist *_clone() const override {
        return new ClexHullDist(*this);
      }

    };




    // --- BaseHull Definitions -------------------

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
      m_calculator_map(_calculator_map),
      m_singular_value_tol(_singular_value_tol),
      m_bottom_facet_tol(_bottom_facet_tol),
      m_selection(_default_selection),
      m_composition_type(_default_composition_type),
      m_initialized(false) {}

    /// \brief Calculates the convex hull
    ///
    /// - Uses the parsed args to determine the selection to use to calculate the hull
    /// - Calls 'init' on the CompCalculator and EnergyCalculator
    /// - Constructs the hull
    ///
    template<typename ValueType>
    void BaseHull<ValueType>::init(const Configuration &_tmplt) const {

      ConstConfigSelection selection;
      const PrimClex &primclex = _tmplt.get_primclex();

      if(m_selection == "ALL") {
        selection = ConstConfigSelection(primclex);
        for(auto it = primclex.config_cbegin(); it != primclex.config_cend(); ++it) {
          selection.set_selected(it->name(), true);
        }
      }
      else if(m_selection == "MASTER") {
        selection = ConstConfigSelection(primclex);
      }
      else if(m_selection == "CALCULATED") {
        selection = ConstConfigSelection(primclex);
        for(auto it = primclex.config_cbegin(); it != primclex.config_cend(); ++it) {
          selection.set_selected(it->name(), is_calculated(*it));
        }

      }
      else {
        if(!fs::exists(m_selection)) {
          throw std::runtime_error(
            std::string("ERROR in $selection argument for '") + short_header(_tmplt) + "'." +
            " Expected <filename>, 'ALL', 'CALCULATED', or 'MASTER' <--default" +
            " No file named '" + m_selection + "'.");
        }
        selection = ConstConfigSelection(primclex, m_selection);
      }

      // Hull typedefs:
      //typedef std::pair<notstd::cloneable_ptr<CompCalculator>,
      //                  notstd::cloneable_ptr<EnergyCalculator> > CalculatorPair;
      //typedef std::map<std::string, CalculatorPair> CalculatorOptions;

      m_calculator_map.find(m_composition_type)->second.first->init(_tmplt);
      m_calculator_map.find(m_composition_type)->second.second->init(_tmplt);

      m_hull = std::make_shared<Hull>(selection,
                                      *m_calculator_map.find(m_composition_type)->second.first,
                                      *m_calculator_map.find(m_composition_type)->second.second,
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



  }
}

#endif
