#include "casm/crystallography/SymType.hh"
#include "casm/crystallography/Lattice.hh"
#include "casm/crystallography/LatticeMap.hh"
#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/crystallography/LatticeIsEquivalent.hh"

namespace CASM {
  namespace xtal {
    StrainCostCalculator::StrainCostCalculator(Eigen::Ref<const Eigen::MatrixXd> const &strain_gram_mat /*= Eigen::MatrixXd::Identity(9,9)*/) {
      if(strain_gram_mat.size() == 0 || strain_gram_mat.isIdentity(1e-9)) {
        m_sym_cost = false;
      }
      else {
        m_sym_cost = true;
        m_gram_mat.resize(6, 6);
        if(strain_gram_mat.rows() == 6 && strain_gram_mat.cols() == 6) {
          std::vector<Index> map({0, 5, 4, 1, 3, 2});
          for(Index i = 0; i < 6; ++i) {
            for(Index j = 0; j < 6; ++j) {
              m_gram_mat(i, j) = strain_gram_mat(map[i], map[j]);
              if(i > 2)
                m_gram_mat(i, j) *= sqrt(2.);
              if(j > 2)
                m_gram_mat(i, j) *= sqrt(2.);

            }
          }
        }
        if(strain_gram_mat.rows() == 9 && strain_gram_mat.cols() == 9) {
          Index m = 0;
          for(Index i = 0; i < 3; ++i) {
            for(Index j = i; j < 3; ++j, ++m) {
              Index n = 0;
              for(Index k = 0; k < 3; ++k) {
                for(Index l = k; l < 3; ++l, ++n) {
                  m_gram_mat(m, n) = strain_gram_mat(i * 3 + j, k * 3 + l);
                  if(m > 2)
                    m_gram_mat(m, n) *= sqrt(2.);
                  if(n > 2)
                    m_gram_mat(m, n) *= sqrt(2.);
                }
              }
            }
          }
        }
      }


    }

    //*******************************************************************************************
    //static function
    double StrainCostCalculator::iso_strain_cost(const Eigen::Matrix3d &F, double _vol_factor) {
      Eigen::Matrix3d tmat = polar_decomposition(F / _vol_factor);

      // -> epsilon=(F_deviatoric-identity)
      return ((tmat - Eigen::Matrix3d::Identity(3, 3)).squaredNorm()
              + (tmat.inverse() - Eigen::Matrix3d::Identity(3, 3)).squaredNorm()) / 6.;
    }

    //*******************************************************************************************
    //static function
    double StrainCostCalculator::iso_strain_cost(const Eigen::Matrix3d &F) {
      return iso_strain_cost(F, vol_factor(F));
    }

    //*******************************************************************************************

    // strain_cost is the mean-square displacement of a point on the surface of a unit sphere
    // when it is deformed by the volume-preserving deformation F_deviatoric = F/det(F)^(1/3)
    double StrainCostCalculator::strain_cost(const Eigen::Matrix3d &F, double _vol_factor) const {

      if(m_sym_cost) {
        double cost = 0;
        m_cache = polar_decomposition(F / _vol_factor);
        m_cache_inv = m_cache.inverse() - Eigen::Matrix3d::Identity(3, 3);
        m_cache -= Eigen::Matrix3d::Identity(3, 3);
        Index m = 0;
        for(Index i = 0; i < 3; ++i) {
          for(Index j = i; j < 3; ++j, ++m) {
            Index n = 0;
            for(Index k = 0; k < 3; ++k) {
              for(Index l = k; l < 3; ++l, ++n) {
                cost += m_gram_mat(m, n) * (m_cache(i, j) * m_cache(j, k) + m_cache_inv(i, j) * m_cache_inv(j, k)) / 6.;
              }
            }
          }
        }
        // geometric factor: (3*V/(4*pi))^(2/3)/3 = V^(2/3)/7.795554179
        return cost;
      }

      return iso_strain_cost(F, _vol_factor);
    }

    //*******************************************************************************************

    double StrainCostCalculator::strain_cost(const Eigen::Matrix3d &F) const {
      return strain_cost(F, vol_factor(F));
    }

    //*******************************************************************************************

    LatticeMap::LatticeMap(const Lattice &_parent,
                           const Lattice &_child,
                           Index num_atoms,
                           int _range,
                           SymOpVector const &_parent_point_group,
                           SymOpVector const &_child_point_group,
                           Eigen::Ref<const Eigen::MatrixXd> const &strain_gram_mat,
                           double _init_better_than /* = 1e20 */) :
      m_calc(strain_gram_mat),
      m_vol_factor(pow(std::abs(volume(_child) / volume(_parent)), 1. / 3.)),
      m_range(_range),
      m_cost(1e20),
      m_currmat(0) {

      Lattice reduced_parent = _parent.reduced_cell();
      m_parent = reduced_parent.lat_column_mat();

      Lattice reduced_child = _child.reduced_cell();
      m_child = reduced_child.lat_column_mat();

      m_U = _parent.inv_lat_column_mat() * m_parent;
      m_V_inv = m_child.inverse() * _child.lat_column_mat();

      if(_range == 1)
        m_mvec_ptr = &unimodular_matrices<1>();
      else if(_range == 2)
        m_mvec_ptr = &unimodular_matrices<2>();
      else if(_range == 3)
        m_mvec_ptr = &unimodular_matrices<3>();
      else if(_range == 4)
        m_mvec_ptr = &unimodular_matrices<4>();
      else
        throw std::runtime_error("LatticeMap cannot currently be invoked for range>4");

      // Construct inverse fractional symops for parent
      {
        LatticeIsEquivalent symcheck(reduced_parent);
        m_parent_fsym_mats.reserve(_parent_point_group.size());
        for(auto const &op : _parent_point_group) {
          if(!symcheck(op))
            continue;
          //if(symcheck.U().isIdentity())
          //continue;
          //std::cout << "Testing point op:\n" << get_matrix(op) << "\n";
          m_parent_fsym_mats.push_back(iround(symcheck.U().inverse()));
          for(Index i = 0; i < (m_parent_fsym_mats.size() - 1); ++i) {
            if(m_parent_fsym_mats[i] == m_parent_fsym_mats.back()) {
              m_parent_fsym_mats.pop_back();
              break;
            }
          }
        }
      }

      // Construct fractional symops for child
      {
        LatticeIsEquivalent symcheck(reduced_child);
        m_child_fsym_mats.reserve(_child_point_group.size());
        for(auto const &op : _child_point_group) {
          if(!symcheck(op))
            continue;
          //if(symcheck.U().isIdentity())
          //continue;
          //std::cout << "Testing point op:\n" << get_matrix(op) << "\n";
          m_child_fsym_mats.push_back(iround(symcheck.U()));
          for(Index i = 0; i < (m_child_fsym_mats.size() - 1); ++i) {
            if(m_child_fsym_mats[i] == m_child_fsym_mats.back()) {
              m_child_fsym_mats.pop_back();
              break;
            }
          }
        }
      }


      reset(_init_better_than);
    }


    LatticeMap::LatticeMap(Eigen::Ref<const LatticeMap::DMatType> const &_parent,
                           Eigen::Ref<const LatticeMap::DMatType> const &_child,
                           Index _num_atoms,
                           int _range,
                           SymOpVector const &_parent_point_group,
                           SymOpVector const &_child_point_group,
                           Eigen::Ref<const Eigen::MatrixXd> const &strain_gram_mat,
                           double _init_better_than /* = 1e20 */) :
      LatticeMap(Lattice(_parent),
                 Lattice(_child),
                 _num_atoms,
                 _range,
                 _parent_point_group,
                 _child_point_group,
                 strain_gram_mat,
                 _init_better_than) {}

    //*******************************************************************************************
    void LatticeMap::reset(double _better_than) {
      m_currmat = 0;

      // From relation F * parent * inv_mat.inverse() = child
      m_F = m_child * inv_mat().cast<double>() * m_parent.inverse(); // -> F

      double tcost = m_calc.strain_cost(m_F, m_vol_factor);

      // Initialize to first valid mapping
      if(tcost <= _better_than && _check_canonical()) {

        m_cost = tcost;
        // reconstruct correct N for unreduced lattice
        m_N = m_U * inv_mat().cast<double>().inverse() * m_V_inv;
      }
      else
        next_mapping_better_than(_better_than);
    }
    //*******************************************************************************************
    /*
     *  For L_child = (*this).lat_column_mat() and L_parent = _parent_lat.lat_column_mat(), we find the mapping:
     *
     *        L_child = F * L_parent * N      -- (Where 'N' is an integer matrix of determinant=1)
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
     *        F.transpose()*F = L_parent.inverse().transpose()*N.inverse().transpose()*L_child.transpose()*L_child*N.inverse()*L_parent.inverse()
     *
     *  The minimal 'C' is returned at the end of the optimization.
     */
    //*******************************************************************************************

    const LatticeMap &LatticeMap::best_strain_mapping() const {
      m_currmat = 0;

      // Get an upper bound on the best mapping by starting with no lattice equivalence
      m_N = DMatType::Identity(3, 3);
      // m_dcache -> value of inv_mat() that gives m_N = identity;
      m_dcache = m_V_inv * m_U;
      m_F = m_child * m_dcache * m_parent.inverse();
      //std::cout << "starting m_F is \n" << m_F << "  det: " << m_F.determinant() << "\n";
      double best_cost = m_calc.strain_cost(m_F, m_vol_factor);
      //std::cout << "starting cost is " << m_cost << "\n";
      //std::cout << "Starting cost is " << m_cost << ", starting N is \n" << m_N << "\nand starting F is \n" << m_F << "\n";
      //std::cout << "Best_cost progression: " << best_cost;
      while(next_mapping_better_than(best_cost).strain_cost() < best_cost) {
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
      // tcost initial value shouldn't matter unles m_inv_count is invalid
      double tcost = max_cost;

      while(++m_currmat < n_mat()) {
        if(!_check_canonical()) {
          //std::cout << "Not canonical: \n" << inv_mat() << "\n";
          continue;
        }

        // From relation F * parent * inv_mat.inverse() = child
        m_F = m_child * inv_mat().cast<double>() * m_parent.inverse(); // -> F
        tcost = m_calc.strain_cost(m_F, m_vol_factor);
        //DEBUG
        //Eigen::Matrix3d M=m_parent.transpose()*m_parent;
        //Eigen::Matrix3d M2=inv_mat().cast<double>().inverse().transpose()*M*inv_mat().cast<double>().inverse();
        //if(almost_equal(M,M2,1e-4)){
        //std::cout << "***POINT OP: \n" << "N: \n" << (m_U * inv_mat().cast<double>().inverse() * m_V_inv) << "\nF: \n" << m_F<< "\ncost: " << tcost << "\n***\n";
        //}
        //\DEBUG
        if(tcost < max_cost) {
          //std::cout << "cost: " << max_cost << " -> " << tcost << "\n";
          m_cost = tcost;

          // need to undo the effect of transformation to reduced cell on 'N'
          // Maybe better to get m_N from m_F instead?  m_U and m_V_inv depend on the lattice reduction
          // that was performed in the constructor, so we would need to store "non-reduced" parent and child
          m_N = m_U * inv_mat().cast<double>().inverse() * m_V_inv;
          //std::cout << "N:\n" << m_N << "\n";
          //  We already have:
          //        m_F = m_child * inv_mat().cast<double>() * m_parent.inverse();
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

    bool LatticeMap::_check_canonical() const {
      // Purpose of jmin is to exclude (i,j)=(0,0) element
      // jmin is set to 0 at end of i=0 pass;
      Index jmin = 1;
      for(Index i = 0; i < m_parent_fsym_mats.size(); ++i, jmin = 0) {
        auto const &inv_parent_op = m_parent_fsym_mats[i];
        for(Index j = jmin; j < m_child_fsym_mats.size(); ++j) {
          auto const &child_op = m_child_fsym_mats[j];
          m_icache = child_op * inv_mat() * inv_parent_op;
          // Skip ops that transform matrix out of range; they won't be enumerated
          if(std::abs(m_icache(0, 0)) > m_range
             || std::abs(m_icache(0, 1)) > m_range
             || std::abs(m_icache(0, 2)) > m_range
             || std::abs(m_icache(1, 0)) > m_range
             || std::abs(m_icache(1, 1)) > m_range
             || std::abs(m_icache(1, 2)) > m_range
             || std::abs(m_icache(2, 0)) > m_range
             || std::abs(m_icache(2, 1)) > m_range
             || std::abs(m_icache(2, 2)) > m_range)
            continue;

          if(std::lexicographical_compare(m_icache.data(), m_icache.data() + 9, inv_mat().data(), inv_mat().data() + 9))
            return false;
        }
      }
      return true;
    }
  }
}
