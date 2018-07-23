#ifndef CASM_Hull
#define CASM_Hull

#include "casm/external/qhull/libqhullcpp/PointCoordinates.h"
#include "casm/external/qhull/libqhullcpp/Qhull.h"

#include "casm/clex/ConfigSelection.hh"
#include "casm/clex/Configuration.hh"

namespace CASM {

  /// \brief Generate and inspect the convex hull generated from a selection of Configurations
  ///
  /// - The underlying Qhull object may not be copy constructed or assigned, so neither can this class. Instead use shared_ptr<Hull>.
  class Hull {

  public:

    typedef VectorXdAttribute<Configuration> CompCalculator;
    typedef ScalarAttribute<Configuration> EnergyCalculator;
    typedef std::pair<notstd::cloneable_ptr<CompCalculator>,
            notstd::cloneable_ptr<EnergyCalculator> > CalculatorPair;
    typedef std::map<std::string, CalculatorPair> CalculatorOptions;

    /// \brief Constructor for convex hull in composition/energy space
    Hull(const ConstConfigSelection &_selection,
         const CompCalculator &_comp_calculator = ConfigIO::SpeciesFrac(),
         const EnergyCalculator &_energy_calculator = ConfigIO::formation_energy_per_species(),
         double _singular_value_tol = 1e-14,
         double _bottom_facet_tol = 1e-14);

    /// \brief const Access the hull object directly
    const orgQhull::Qhull &data() const;

    /// \brief Orthogonal transformation matrix from a point in full comp/energy space to dimension-reduced comp/energy space
    const Eigen::MatrixXd &reduce() const;

    /// \brief Orthogonal transformation matrix from a point in the dimension-reduced comp/energy space to the full comp/energy space
    const Eigen::Transpose<const Eigen::MatrixXd> expand() const;

    /// \brief Return the configuration corresponding to any point
    const Configuration &configuration(const orgQhull::QhullPoint &point) const;

    /// \brief Return the configuration corresponding to a hull vertex
    const Configuration &configuration(const orgQhull::QhullVertex &vertex) const;

    /// \brief Return the chemical potential corresponding to a facet
    Eigen::VectorXd mu(const orgQhull::QhullFacet facet) const;

    /// \brief Return the 0K ground state corresponding to the input chemical potential
    const Configuration &groundstate(const Eigen::VectorXd &mu) const;

    /// \brief Use the EnergyCalculator to return the energy of a Configuration
    double energy(const Configuration &config) const;

    /// \brief Use the CompCalculator to return the composition of a Configuration
    Eigen::VectorXd composition(const Configuration &config) const;

    /// \brief Return a vector corresponding to the coordinate of a given configuration in full composition/energy space
    Eigen::VectorXd point(const Configuration &config) const;

    /// \brief Return a vector corresponding to the coordinate of a given configuration in the reduced composition/energy space
    Eigen::VectorXd reduced_point(const Configuration &config) const;

    /// \brief The distance a Configuration is above the hull along the energy axis
    double dist_to_hull(const Configuration &config) const;

    /// \brief The distance a point in the reduced composition/energy space is above the hull along the energy axis
    double dist_to_hull(Eigen::VectorXd _reduced_point) const;


  private:

    struct CompareVertex {
      CompareVertex() {}

      bool operator()(orgQhull::QhullVertex A, orgQhull::QhullVertex B) const {
        return A.id() < B.id();
      }
    };

    // the hull object
    orgQhull::Qhull m_hull;

    // the selection of Configurations used to generate the hull
    ConstConfigSelection m_selection;

    // get composition coordinates for Configuration
    notstd::cloneable_ptr<CompCalculator> m_comp_calculator;

    // get energy for a Configuration
    notstd::cloneable_ptr<EnergyCalculator> m_energy_calculator;

    // transform full dimension comp/energy vector onto subspace range(comp)/energy
    Eigen::MatrixXd m_reduce;

    // all 'bottom' facets, along with projection of 'm_bottom' along the unit outward normal
    std::vector< std::pair<orgQhull::QhullFacet, double> > m_bottom_facets;

    // the vertices on the bottom hull
    std::set<orgQhull::QhullVertex, CompareVertex> m_bottom_vertices;

  };

}

#endif
