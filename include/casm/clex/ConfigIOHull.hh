#ifndef CONFIGIOHULL_HH
#define CONFIGIOHULL_HH

#include "casm/casm_io/DataFormatter.hh"
#include "casm/hull/Hull.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/database/Selection.hh"

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
      bool init(const Configuration &_tmplt) const override;

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

  }
}

#endif

