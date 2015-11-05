#ifndef CONFIGIOHULL_HH
#define CONFIGIOHULL_HH

#include "casm/casm_io/DataFormatterTools.hh"
#include "casm/hull/Hull.hh"
#include "casm/clex/ConfigIterator.hh"
#include "casm/clex/PrimClex.hh"

namespace CASM {

  class Configuration;

  namespace ConfigIO_impl {
    
    /// \brief Base class for hull info formatters
    template<typename ValueType>
    class BaseHullConfigFormatter: public BaseValueFormatter<ValueType, Configuration> {
    
      public:
      
      /// \brief Constructor
      BaseHullConfigFormatter(const std::string &_name, 
                              const std::string &_desc, 
                              const std::string &_default_selection,
                              const std::string &_default_composition_type,
                              double _singular_value_tol);
      
      /// \brief Calculates the convex hull
      void init(const Configuration &_tmplt) const override;

      /// \brief column header to use
      std::string short_header(const Configuration &_config) const override;
      
      /// \brief Determine the selection to use to generate the hull
      bool parse_args(const std::string &args) override;
    
    
    protected:
      
      /// \brief const Access the Hull object
      const Hull& _hull() const;
      
      // for parse_args: determine energy type based on composition type
      std::map<std::string, std::pair<Hull::CompCalculator,Hull::EnergyCalculator> > m_calculator_map;
      
      // detect and avoid printing -0.0
      constexpr static double m_dist_to_hull_tol = 1e-14;
      
    private:
      
      // tol used to detect zero singular values during principal component analysis preprocessing before hull calculation
      double m_singular_value_tol;
      
      // the hull object
      mutable std::shared_ptr<Hull> m_hull;
      
      // Parsed arguments
      //  -- what selection to use for constructing the hull
      std::string m_selection;
      
      // Parsed arguments
      //  -- what composition to use for constructing the hull
      //  -- determines comp calculator and energy calculator, via m_calculator_map
      std::string m_composition_type;
      
    };
    
    /// \brief Returns a boolean indicating if a Configuration is a convex hull vertex
    ///
    /// - Only one Configuration out of set that have identical or almost identical 
    ///   points in composition/energy space will return true 
    class OnHullConfigFormatter: public BaseHullConfigFormatter<bool> {
      
      public:
      
      /// \brief Constructor
      OnHullConfigFormatter();
      
      /// \brief Clone
      BaseDatumFormatter<Configuration>* clone() const override {
        return new OnHullConfigFormatter(*this);
      }
      
      
      protected:
      
      /// \brief Validate that the Configuration has a formation energy per species
      bool _validate(const Configuration &_config) const override;
      
      /// \brief Check if the Configuration is a hull vertex 
      bool _evaluate(const Configuration &_config) const override;
    };
    
    /// \brief Returns the distance of a Configuration above the convex hull along the energy axis
    class HullDistConfigFormatter: public BaseHullConfigFormatter<double> {
      
      public:
      
      /// \brief Constructor
      HullDistConfigFormatter();
      
      /// \brief Clone
      BaseDatumFormatter<Configuration>* clone() const override {
        return new HullDistConfigFormatter(*this);
      }
      
      protected:
      
      /// \brief Validate that the Configuration has a formation energy per species
      bool _validate(const Configuration &_config) const override;
      
      /// \brief Return the distance to the hull
      double _evaluate(const Configuration &_config) const override;
    };
    
    
    /// \brief Returns a boolean indicating if a Configuration is a convex hull vertex
    ///
    /// - Only one Configuration out of set that have identical or almost identical 
    ///   points in composition/energy space will return true 
    class OnClexHullConfigFormatter: public BaseHullConfigFormatter<bool> {
      
      public:
      
      /// \brief Constructor
      OnClexHullConfigFormatter();
      
      /// \brief Clone
      BaseDatumFormatter<Configuration>* clone() const override {
        return new OnClexHullConfigFormatter(*this);
      }
      
      
      protected:
      
      /// \brief Validate that the Configuration has a cluster expanded formation energy per species
      virtual bool _validate(const Configuration &_config) const override;
      
      /// \brief Check if the Configuration is a hull vertex 
      virtual bool _evaluate(const Configuration &_config) const override;
    };
    
    /// \brief Returns the distance of a Configuration above the convex hull along the energy axis
    class ClexHullDistConfigFormatter: public BaseHullConfigFormatter<double> {
      
      public:
      
      /// \brief Constructor
      ClexHullDistConfigFormatter();
      
      /// \brief Clone
      BaseDatumFormatter<Configuration>* clone() const override {
        return new ClexHullDistConfigFormatter(*this);
      }
      
      protected:
      
      /// \brief Validate that the Configuration has a cluster expanded formation energy per species
      virtual bool _validate(const Configuration &_config) const override;
      
      /// \brief Return the distance to the hull
      virtual double _evaluate(const Configuration &_config) const override;
    };
    
    
    
    
    // --- BaseHullConfigFormatter Definitions -------------------
    
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
    BaseHullConfigFormatter<ValueType>::BaseHullConfigFormatter(const std::string &_name, 
                                                     const std::string &_desc, 
                                                     const std::string &_default_selection,
                                                     const std::string &_default_composition_type,
                                                     double _singular_value_tol) :
      BaseValueFormatter<ValueType, Configuration>(_name, _desc), 
      m_selection(_default_selection),
      m_composition_type(_default_composition_type),
      m_singular_value_tol(_singular_value_tol) {}

    /// \brief Calculates the convex hull
    /// 
    /// Uses the parsed args to determine the selection to use to calculate the hull
    ///
    template<typename ValueType>
    void BaseHullConfigFormatter<ValueType>::init(const Configuration &_tmplt) const {
      
      ConstConfigSelection selection;
      const PrimClex& primclex = _tmplt.get_primclex();
      
      if(m_selection == "ALL") {
        selection = ConstConfigSelection(primclex);
        for(auto it=primclex.config_cbegin(); it!=primclex.config_cend(); ++it) {
          selection.set_selected(it->name(), true);
        }
      }
      else if(m_selection == "MASTER") {
        selection = ConstConfigSelection(primclex);
      }
      else if(m_selection == "CALCULATED") {
        selection = ConstConfigSelection(primclex);
        for(auto it=primclex.config_cbegin(); it!=primclex.config_cend(); ++it) {
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
      
      m_hull = std::make_shared<Hull>(selection, 
                                      m_calculator_map.find(m_composition_type)->second.first, 
                                      m_calculator_map.find(m_composition_type)->second.second, 
                                      m_singular_value_tol);
      
    }

    /// \brief column header to use
    ///
    /// \returns 'name() + '(' + args + ')'
    /// - ex: 'on_hull(MASTER)' if name() is 'on_hull', and args is "MASTER"
    /// - ex: 'on_clex_hull(all) if name() is 'on_clex_hull', and args is "all" or there are no args
    ///
    template<typename ValueType>
    std::string BaseHullConfigFormatter<ValueType>::short_header(const Configuration &_config) const {
      std::stringstream t_ss;
      t_ss << this->name() << "(" << m_selection << "," << m_composition_type << ")";
      return t_ss.str();
    }
    
    /// \brief Determine the selection to use to generate the hull
    ///
    /// Options are:
    /// - "all", use all configurations (default for no args)
    /// - "MASTER", use the current MASTER selection
    /// - other, assume the argument is the filename for a selection to use
    ///
    template<typename ValueType>
    bool BaseHullConfigFormatter<ValueType>::parse_args(const std::string &args) {
      
      std::vector<std::string> splt_vec;
      boost::split(splt_vec, args, boost::is_any_of(","), boost::token_compress_on);
      
      if(splt_vec.size() > 2) {
        throw std::runtime_error("Attempted to initialize format tag " + this->name()
                                 + " with " + std::to_string(splt_vec.size()) + " arguments ("
                                 + args + "), but only up to 2 arguments allowed.\n");
        return false;
      }
      
      while(splt_vec.size() < 2) {
        splt_vec.push_back("");
      }
      
      m_selection = splt_vec[0].empty() ? m_selection : splt_vec[0];
      m_composition_type = splt_vec[1].empty() ? m_composition_type : splt_vec[1];
      
      auto it = m_calculator_map.find(m_composition_type);
      if(it == m_calculator_map.end()) {
        
        std::stringstream ss;
        ss << "Composition type " << m_composition_type << " is not recognized. Options are: ";
        for(auto it=m_calculator_map.begin(); it!=m_calculator_map.end(); ++it) {
          std::cout << it->first << " ";
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
    const Hull& BaseHullConfigFormatter<ValueType>::_hull() const {
      return *m_hull;
    }
    
    
    
  }
}
#endif

