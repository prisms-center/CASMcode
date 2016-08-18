#include "casm/hull/Hull.hh"

#include <iterator>

#include "casm/clex/PrimClex.hh"
#include "casm/misc/PCA.hh"
#include "casm/external/Eigen/Dense"
#include "casm/external/qhull/libqhullcpp/QhullFacetList.h"
#include "casm/external/qhull/libqhullcpp/QhullVertexSet.h"
#include "casm/casm_io/DataFormatter.hh"
#include "casm/casm_io/DataFormatterTools.hh"

namespace CASM {

  namespace Hull_impl {

    /// \brief Print informational message and throw exception if input data is not valid
    void _validate_input(const ConstConfigSelection &selection,
                         const Hull::CompCalculator &comp_calculator,
                         const Hull::EnergyCalculator &energy_calculator);
  }


  /// \brief Constructor for convex hull in atom_frac & Ef/atom space
  Hull::Hull(const ConstConfigSelection &_selection,
             const CompCalculator &_comp_calculator,
             const EnergyCalculator &_energy_calculator,
             double _singular_value_tol,
             double _bottom_facet_tol) :
    m_selection(_selection),
    m_comp_calculator(_comp_calculator),
    m_energy_calculator(_energy_calculator) {

    Hull_impl::_validate_input(m_selection, *m_comp_calculator, *m_energy_calculator);

    // get the number of configurations and compositions (full composition space)
    Index Nselected = std::distance(m_selection.selected_config_begin(), m_selection.selected_config_end());
    Index Ncomp = composition(*m_selection.selected_config_begin()).size();

    // generate initial set of points (col vector matrix)
    Eigen::MatrixXd mat(Ncomp + 1, Nselected);

    Index i = 0;
    for(auto it = m_selection.selected_config_begin(); it != m_selection.selected_config_end(); ++it) {
      mat.block(0, i, Ncomp, 1) = composition(*it);
      mat(Ncomp, i) = energy(*it);
      ++i;
    }

    // principal component analysis to get rotation matrix
    PCA pca(mat.topRows(Ncomp), _singular_value_tol);
    m_reduce = pad(pca.reduce(), 1);

    Eigen::MatrixXd reduced_mat = m_reduce * mat;

    // Construct Qhull PointCoordiantes, with correct dimension
    orgQhull::PointCoordinates points(reduced_mat.rows(), "");

    // use matrix data to set point coordinates
    points.append(reduced_mat.size(), reduced_mat.data());

    // calculate hull
    std::string qh_command = "";
    m_hull.runQhull(points.comment().c_str(), points.dimension(), points.count(), &*points.coordinates(), qh_command.c_str());

    // check for errors
    if(m_hull.hasQhullMessage()) {
      std::cerr << "\nQhull message:\n" << m_hull.qhullMessage();
      m_hull.clearQhullMessage();
      throw std::runtime_error("Qhull Error (see message)");
    }

    // want to find minimum distance from point to a bottom facet

    // V: any vector from point to facet
    // N: outward unit normal of facet
    // D: unit 'down' direction (by convention always [0, 0, ..., -1])
    //
    // distance to facet along normal: V.dot(N) = a (Qhull reports this as a negative number)
    // proj of D onto N: D.dot(N) = b (if b is positive, then this is 'bottom' facet)
    // distance along D to facet: a/b

    // collect bottom facets (along with 'b')
    {
      auto begin = m_hull.facetList().begin();
      auto end = m_hull.facetList().end();
      int dim = m_hull.dimension();
      double b;
      Eigen::Map<const Eigen::VectorXd> outnorm((*begin).hyperplane().begin(), dim);

      for(auto facet_it = begin; facet_it != end; ++facet_it) {

        // reset location of outnorm data
        new(&outnorm) Eigen::Map<const Eigen::VectorXd>((*facet_it).hyperplane().begin(), dim);
        b = -outnorm(dim - 1);

        if(b > _bottom_facet_tol) {
          m_bottom_facets.push_back(std::make_pair(*facet_it, b));
        }
      }
    }


    // find set of bottom vertices using bottom facets
    // use pointers for the set comparison
    for(auto facet_it = m_bottom_facets.begin(); facet_it != m_bottom_facets.end(); ++facet_it) {
      orgQhull::QhullVertexSet vertices = (*facet_it).first.vertices();
      for(auto vertex_it = vertices.begin(); vertex_it != vertices.end(); ++vertex_it) {
        m_bottom_vertices.insert(*vertex_it);
      }
    }

  }

  /// \brief const Access the hull object directly
  const orgQhull::Qhull &Hull::data() const {
    return m_hull;
  }

  /// \brief Orthogonal transformation matrix from a point in full comp/energy space to dimension-reduced comp/energy space
  ///
  /// - If the points in the provided configuration selection do not span the entire composition
  ///   space, a principle component analysis identifies the subspace that is spanned and this
  ///   matrix rotates points into that subspace
  const Eigen::MatrixXd &Hull::reduce() const {
    return m_reduce;
  }

  /// \brief Orthogonal transformation matrix from a point in the dimension-reduced comp/energy space to the full comp/energy space
  const Eigen::Transpose<const Eigen::MatrixXd> Hull::expand() const {
    return m_reduce.transpose();
  }

  /// \brief Return the configuration corresponding to any point
  const Configuration &Hull::configuration(const orgQhull::QhullPoint &point) const {
    auto it = m_selection.selected_config_cbegin();
    std::advance(it, point.id());
    return *it;
  }

  /// \brief Return the configuration corresponding to a hull vertex
  const Configuration &Hull::configuration(const orgQhull::QhullVertex &vertex) const {
    return configuration(vertex.point());
  }

  /// \brief Return the chemical potential corresponding to a facet
  ///
  /// - Calculated from the normal of the provided facet, so units depend on composition/energy space
  Eigen::VectorXd Hull::mu(const orgQhull::QhullFacet facet) const {

    // mu_i = dG/dx_i
    // hyperplane normal = [dx_i, dx_j, ... dG] (normalized)

    int dim = m_hull.dimension();
    Eigen::Map<const Eigen::VectorXd> reduced_hyperplane(facet.hyperplane().begin(), dim);

    Eigen::VectorXd _mu = Eigen::VectorXd::Constant(dim - 1, reduced_hyperplane(dim - 1));

    return _mu.cwiseQuotient(expand().leftCols(dim - 1) * reduced_hyperplane.head(dim - 1));
  }

  /// \brief Return the 0K ground state corresponding to the input chemical potential
  ///
  /// - Units must be appropriate given the composition/energy space
  const Configuration &Hull::groundstate(const Eigen::VectorXd &mu) const {

    // mu_i = dG/dx_i
    // hyperplane normal = [dx_i, dx_j, ... dG] (normalized)


    // normal vector to hyperplane
    Eigen::VectorXd hyperplane(mu.size() + 1);
    hyperplane.head(mu.size()) = mu.cwiseInverse();
    hyperplane(mu.size()) = 1.0;
    hyperplane.normalize();

    Eigen::VectorXd reduced_hyperplane = reduce() * hyperplane;

    double max_offset = std::numeric_limits<double>::min();
    orgQhull::QhullVertex _groundstate = *m_bottom_vertices.begin();

    // m_bottom_vertices: a std::set<orgQhull::QhullVertex*>
    for(auto it = m_bottom_vertices.begin(); it != m_bottom_vertices.end(); ++it) {

      Eigen::Map<const Eigen::VectorXd> reduced_vertex_point(it->point().begin(), it->point().size());

      double offset = reduced_hyperplane.dot(reduced_vertex_point);

      if(offset > max_offset) {
        max_offset = offset;
        _groundstate = *it;
      }
    }

    return configuration(_groundstate);
  }

  /// \brief Use the EnergyCalculator to return the energy of a Configuration
  double Hull::energy(const Configuration &config) const {
    return (*m_energy_calculator)(config);
  }

  /// \brief Use the CompCalculator to return the composition of a Configuration
  Eigen::VectorXd Hull::composition(const Configuration &config) const {
    return (*m_comp_calculator)(config);
  }


  /// \brief Return a vector corresponding to the coordinate of a given configuration in full comp/energy space
  Eigen::VectorXd Hull::point(const Configuration &config) const {
    Eigen::VectorXd _point(m_reduce.cols());
    _point.head(m_reduce.cols() - 1) = composition(config);
    _point.tail(1)(0) = energy(config);
    return _point;
  }

  /// \brief Return a vector corresponding to the coordinate of a given configuration in the reduced composition/energy space
  Eigen::VectorXd Hull::reduced_point(const Configuration &config) const {
    return reduce() * point(config);
  }

  /// \brief The distance a Configuration is above the hull along the energy axis
  double Hull::dist_to_hull(const Configuration &config) const {
    return dist_to_hull(reduced_point(config));
  }

  /// \brief The distance a point in the reduced composition/energy space is above the hull along the energy axis
  double Hull::dist_to_hull(Eigen::VectorXd _reduced_point) const {

    // want to find minimum distance from point to a bottom facet

    // V: any vector from point to facet
    // N: outward unit normal of facet
    // D: unit 'down' direction (by convention always [0, 0, ..., -1])
    //
    // distance to facet along normal: V.dot(N) = a (Qhull reports this as a negative number)
    // proj of D onto N: D.dot(N) = b (if b is positive, then this is 'bottom' facet)
    // distance along D to facet: a/b

    orgQhull::QhullPoint qpoint(_reduced_point.size(), _reduced_point.data());

    double a, b;
    double dist_to_hull = std::numeric_limits<double>::max();

    for(auto facet_it = m_bottom_facets.begin(); facet_it != m_bottom_facets.end(); ++facet_it) {

      // Qhull distance is negative for internal points
      a = -(facet_it->first).distance(qpoint);
      b = facet_it->second;
      if(a / b < dist_to_hull) {
        dist_to_hull = a / b;
      }
    }

    return dist_to_hull;
  }

  namespace Hull_impl {

    /// \brief Print informational message and throw exception if input data is not valid
    void _validate_input(const ConstConfigSelection &selection,
                         const Hull::CompCalculator &comp_calculator,
                         const Hull::EnergyCalculator &energy_calculator) {

      typedef std::map<std::string, std::pair<bool, bool> > CheckMap;

      CheckMap invalid_data;
      for(auto it = selection.selected_config_cbegin(); it != selection.selected_config_cend(); ++it) {
        if(!comp_calculator.validate(*it) || !energy_calculator.validate(*it)) {
          invalid_data[it.name()] = std::make_pair(comp_calculator.validate(*it), energy_calculator.validate(*it));
        }
      }

      if(invalid_data.size()) {
        std::cerr << "Invalid data for hull construction:\n";
        typedef GenericDatumFormatter<std::string, CheckMap::value_type> Formatter;

        DataFormatter<CheckMap::value_type> f(
        Formatter("configname", "", [](CheckMap::value_type v) {
          return v.first;
        }),
        Formatter("composition", "", [](CheckMap::value_type v) {
          return v.second.first ? "OK" : "invalid";
        }),
        Formatter("energy", "", [](CheckMap::value_type v) {
          return v.second.second ? "OK" : "invalid";
        })
        );

        std::cerr << f(invalid_data.begin(), invalid_data.end()) << "\n";

        std::stringstream ss;
        ss << "Error in Hull(): Invalid composition or energy data. \n"
           "Make sure you have set composition axes, all selected configurations\n"
           "have calculation results, and you have set your chemical reference.";
        throw std::runtime_error(ss.str());
      }
    }

  }


}


