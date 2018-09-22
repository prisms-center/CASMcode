#include "casm/crystallography/Lattice.hh"
#include "casm/crystallography/LatticeMap.hh"
namespace CASM {
  LatticeMap::LatticeMap(const Lattice &_ideal, const Lattice &_strained, Index num_atoms, double _tol/*=TOL*/, int _range/*=2*/) :
    m_L1(Eigen::Matrix3d(_ideal.get_reduced_cell().lat_column_mat())), m_L2(Eigen::Matrix3d(_strained.get_reduced_cell().lat_column_mat())),
    m_scale(pow(std::abs(m_L2.determinant() / m_L1.determinant()), 1.0 / 3.0)), m_atomic_vol(std::abs(m_L2.determinant() / (double)num_atoms)),  m_tol(_tol), m_cost(1e20),
    m_inv_count(-std::abs(_range) * IMatType::Ones(3, 3), std::abs(_range) * IMatType::Ones(3, 3), IMatType::Ones(3, 3)) {

    m_U = Eigen::Matrix3d(_ideal.inv_lat_column_mat()) * m_L1;
    m_V_inv = m_L2.inverse() * Eigen::Matrix3d(_strained.lat_column_mat());
    // Initialize to first valid mapping
    next_mapping_better_than(1e10);
  }

  //*******************************************************************************************
  /*
   *  For L_strained = (*this).lat_column_mat() and L_ideal = _ideal_lat.lat_column_mat(), we find the mapping:
   *
   *        L_strained = F * L_ideal * N      -- (Where 'N' is an integer matrix of determinant=1)
   *
   *  That minimizes the cost function:
   *
   *        C = pow((*this).volume(),2/3) * trace(E_dev.transpose()*E_dev.transpose()) / 4
   *
   *  Where E_dev is the non-volumetric component of the green-lagrange strain
   *
   *        E_dev = E - trace(E/3)*Identity   -- (such that E_dev is traceless)
   *
   *  where the Green-Lagrange strain is
   *
   *        E = (F.transpose()*F-Identity)/2
   *
   *  The cost function approximates the mean-square-displacement of a point in a cube of volume (*this).volume() when it is
   *  deformed by deformation matrix 'F', but neglecting volumetric effects
   *
   *  The algorithm proceeds by counting over 'N' matrices (integer matrices of determinant 1) with elements on the interval (-2,2).
   *  (we actually count over N.inverse(), because....)
   *
   *  The green-lagrange strain for that 'N' is then found using the relation
   *
   *        F.transpose()*F = L_ideal.inverse().transpose()*N.inverse().transpose()*L_strained.transpose()*L_strained*N.inverse()*L_ideal.inverse()
   *
   *  The minimal 'C' is returned at the end of the optimization.
   */
  //*******************************************************************************************

  const LatticeMap &LatticeMap::best_strain_mapping() const {
    m_inv_count.reset();

    // Get an upper bound on the best mapping by starting with no lattice equivalence
    m_N = DMatType::Identity(3, 3);
    // m_cache -> value of m_inv_count that gives m_N = identity;
    m_cache = m_V_inv * m_U;
    m_F = m_L2 * m_cache * m_L1.inverse();
    //std::cout << "starting m_F is \n" << m_F << "  det: " << m_F.determinant() << "\n";
    double best_cost = _calc_strain_cost();
    //std::cout << "starting cost is " << m_cost << "\n";
    //std::cout << "Starting cost is " << m_cost << ", starting N is \n" << m_N << "\nand starting F is \n" << m_F << "\n";
    //std::cout << "Best_cost progression: " << best_cost;
    while(next_mapping_better_than(best_cost - m_tol).strain_cost() < best_cost) {
      best_cost = strain_cost();
      //std::cout << "    " << best_cost;
    }

    m_cost = best_cost;
    return *this;
  }
  //*******************************************************************************************

  const LatticeMap &LatticeMap::next_mapping_better_than(double max_cost) const {
    m_cost = 1e20;
    return _next_mapping_better_than(max_cost);
  }
  //*******************************************************************************************
  // Implements the algorithm as above, with generalized inputs:
  //       -- m_inv_count saves the state between calls
  //       -- the search breaks when a mapping is found with cost < max_cost
  const LatticeMap &LatticeMap::_next_mapping_better_than(double max_cost) const {

    DMatType init_F(m_F);
    double tcost = max_cost;
    for(++m_inv_count; m_inv_count.valid(); ++m_inv_count) {
      //continue if determinant is not 1, because it doesn't preserve volume
      if(!almost_equal(std::abs(m_inv_count().determinant()), 1))
        continue;

      m_F = m_L2 * m_inv_count().cast<double>() * m_L1.inverse(); // -> F
      tcost = _calc_strain_cost();
      if(tcost < max_cost) {
        m_cost = tcost;

        // need to undo the effect of transformation to reduced cell on 'N'
        // Maybe better to get m_N from m_F instead?  m_U and m_V_inv depend on the lattice reduction
        // that was performed in the constructor, so we would need to store "non-reduced" L1 and L2
        m_N = m_U * m_inv_count().cast<double>().inverse() * m_V_inv;
        //  We already have:
        //        m_F = m_L2 * m_inv_count().cast<double>() * m_L1.inverse();
        break;
      }
    }
    if(!(tcost < max_cost)) {
      // If no good mappings were found, uncache the starting value of m_F
      m_F = init_F;
      // m_N hasn't changed if tcost>max_cost
      // m_cost hasn't changed either
    }
    // m_N, m_F, and m_cost will describe the best mapping encountered, even if nothing better than max_cost was encountered


    //std::cout << "Final N is:\n" << N << "\n\nAND final cost is " << min_cost << "\n\n";
    return *this;
  }

  //*******************************************************************************************

  // strain_cost is the mean-square displacement of a point on the surface of a sphere having volume = relaxed_atomic_vol
  // when it is deformed by the volume-preserving deformation F_deviatoric = F/det(F)^(1/3)
  double LatticeMap::calc_strain_cost(const Eigen::Matrix3d &F, double relaxed_atomic_vol) {
    // -> epsilon=(F_deviatoric-identity)
    Eigen::Matrix3d cache = 0.5 * (F.transpose() * F / pow(std::abs(F.determinant()), 2.0 / 3.0) - Eigen::Matrix3d::Identity(3, 3));

    // geometric factor: (3*V/(4*pi))^(2/3)/3 = V^(2/3)/7.795554179
    return std::pow(std::abs(relaxed_atomic_vol), 2.0 / 3.0) * cache.squaredNorm() / 7.795554179;
  }
  //*******************************************************************************************
  double LatticeMap::_calc_strain_cost() const {
    // -> epsilon=(F_deviatoric-identity)
    m_cache.noalias() = 0.5 * (m_F.transpose() * m_F / (m_scale * m_scale) - Eigen::Matrix3d::Identity(3, 3));

    // geometric factor: (3*V/(4*pi))^(2/3)/3 = V^(2/3)/7.795554179
    return std::pow(m_atomic_vol, 2.0 / 3.0) * m_cache.squaredNorm() / 7.795554179;
  }


}
