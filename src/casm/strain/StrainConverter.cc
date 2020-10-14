#include "casm/crystallography/Strain.hh"
#include "casm/symmetry/SymGroup.hh"
#include "casm/symmetry/SymGroupRep.hh"
#include "casm/symmetry/SymRepTools.hh"
#include "casm/symmetry/SymMatrixXd.hh"
#include "casm/strain/StrainConverter.hh"

namespace CASM {

  //Calculates the metric tensor of the the deformation gradient as
  //F^{T}F
  Matrix3d StrainConverter::metric_tensor(const Eigen::Ref<const Matrix3d> &F) {
    return strain::metric_tensor(F);
  }

  //*******************************************************************************************
  //calculates and returns the value of U where F = R*U
  Matrix3d StrainConverter::right_stretch_tensor(Matrix3d &C, const Eigen::Ref<const Matrix3d> &F) {
    C = strain::metric_tensor(F);
    return strain::right_stretch_tensor(F);
  }

  //*******************************************************************************************
  //overloaded version of the above
  Matrix3d StrainConverter::right_stretch_tensor(const Eigen::Ref<const Matrix3d> &F) {
    return strain::right_stretch_tensor(F);
  }

  //*******************************************************************************************
  /// GREEN_LAGRANGE = 1/2 * (F^{T} F - I)
  Matrix3d StrainConverter::green_lagrange(const Eigen::Ref<const Matrix3d> &F) {
    return strain::deformation_tensor_to_metric<STRAIN_METRIC::GREEN_LAGRANGE>(F);
  }

  //*******************************************************************************************
  /// GREEN_LAGRANGE = 1/2 * (F^{T} F - I)
  Matrix3d StrainConverter::green_lagrange_to_F(const Eigen::Ref<const Matrix3d> &E) {
    return strain::metric_to_deformation_tensor<STRAIN_METRIC::GREEN_LAGRANGE>(E);
  }

  //*******************************************************************************************
  /// BIOT = (U-I)
  Matrix3d StrainConverter::biot(const Eigen::Ref<const Matrix3d> &F) {
    return strain::deformation_tensor_to_metric<STRAIN_METRIC::BIOT>(F);
  }

  //*******************************************************************************************
  /// BIOT = (U-I)
  Matrix3d StrainConverter::biot_to_F(const Eigen::Ref<const Matrix3d> &B) {
    return strain::metric_to_deformation_tensor<STRAIN_METRIC::BIOT>(B);
  }

  //*******************************************************************************************
  /// HENCKY = log(C)/2
  Matrix3d StrainConverter::hencky(const Eigen::Ref<const Matrix3d> &F) {
    return strain::deformation_tensor_to_metric<STRAIN_METRIC::HENCKY>(F);
  }

  //*******************************************************************************************
  /// HENCKY = log(C)/2
  Matrix3d StrainConverter::hencky_to_F(const Eigen::Ref<const Matrix3d> &H) {
    return strain::metric_to_deformation_tensor<STRAIN_METRIC::HENCKY>(H);
  }

  //*******************************************************************************************
  /// EULER_ALMANSI = (I-(F F^{T})^(-1))/2
  Matrix3d StrainConverter::euler_almansi(const Eigen::Ref<const Matrix3d> &F) {
    return strain::deformation_tensor_to_metric<STRAIN_METRIC::EULER_ALMANSI>(F);
  }

  //*******************************************************************************************
  /// EULER_ALMANSI = (I-(F F^{T})^(-1))/2
  Matrix3d StrainConverter::euler_almansi_to_F(const Eigen::Ref<const Matrix3d> &A) {
    return strain::metric_to_deformation_tensor<STRAIN_METRIC::EULER_ALMANSI>(A);
  }

  //*******************************************************************************************
  /// DISP_GRAD = F
  Matrix3d StrainConverter::disp_grad(Eigen::Ref<const Matrix3d> const &F) {
    return F;
  }


  //*******************************************************************************************
  //Calculates the strain metric based on what mode is passed
  //in. Allowed modes are listed in STRAIN_METRIC
  Matrix3d StrainConverter::strain_metric(Eigen::Ref<const Matrix3d> const &F, STRAIN_METRIC MODE) {
    /// GREEN_LAGRANGE = 1/2 * (F^{T} F - I)
    if(MODE == STRAIN_METRIC::GREEN_LAGRANGE) {
      return green_lagrange(F);
    }
    /// BIOT = (U-I)
    else if(MODE == STRAIN_METRIC::BIOT) {
      return biot(F);
    }
    /// HENCKY = log(C)/2
    else if(MODE == STRAIN_METRIC::HENCKY) {
      return hencky(F);
    }
    /// EULER_ALMANSI = 0.5 * (I-(F F^{T})^(-1))
    else if(MODE == STRAIN_METRIC::EULER_ALMANSI) {
      return euler_almansi(F);
    }
    /// STRETCH
    else if(MODE == STRAIN_METRIC::STRETCH) {
      return right_stretch_tensor(F);
    }
    /// DISP_GRAD
    else if(MODE == STRAIN_METRIC::DISP_GRAD) {
      return disp_grad(F);
    }
    else {
      std::cerr << "CRITICAL ERROR: In StrainConverter::strain_metric. Unrecognized strain metric" << std::endl;
      std::cerr << "                Your only options are GREEN_LAGRANGE, BIOT, HENCKY, EULER_ALMANSI,STRETCH and DISP_GRAD" << std::endl;
      std::cerr << "                Exiting..." << std::endl;
      exit(1);
    }
  }

  //*******************************************************************************************
  //Calculates the strain metric based on what mode is passed
  //in. Allowed modes are listed in STRAIN_METRIC
  Matrix3d StrainConverter::strain_metric(Eigen::Ref<const Matrix3d> const &F) const {
    assert(curr_metric_func && "StrainConverter object improperly initialized!");
    return (*curr_metric_func)(F);
  }

  //*******************************************************************************************
  //Calculates the strain metric based on what mode is passed
  //in. Allowed modes are listed in STRAIN_METRIC
  Matrix3d StrainConverter::strain_metric_to_F(Eigen::Ref<const Matrix3d> const &E) const {
    assert(curr_inv_metric_func && "StrainConverter object improperly initialized!");
    return (*curr_inv_metric_func)(E);
  }

  //*******************************************************************************************
  //Returns the symmetrically unique elements of E (assuming your
  //strain metric is symmetric) ordered in a manner decided by
  //m_order_strain
  VectorXd StrainConverter::unroll_E(Eigen::Ref<const Matrix3d> const &E) const {
    VectorXd _unrolled_E(m_order_strain.size());
    for(Index i = 0; i < m_order_strain.size(); i++) {
      _unrolled_E(i) = m_weight_strain[i] * E(m_order_strain[i][0], m_order_strain[i][1]);
    }
    return _unrolled_E;
  }

  //*******************************************************************************************
  //Returns the symmetrically unique elements of E (assuming your
  //strain metric is symmetric) ordered in a manner decided by
  //m_order_strain
  Matrix3d StrainConverter::rollup_E(Eigen::Ref<const VectorXd> const &_unrolled_E) const {
    Matrix3d _rolled_E;
    for(Index i = 0; i < m_order_strain.size(); i++) {
      _rolled_E(m_order_strain[i][0], m_order_strain[i][1]) = _unrolled_E(i) / m_weight_strain[i];
      _rolled_E(m_order_strain[i][1], m_order_strain[i][0]) = _unrolled_E(i) / m_weight_strain[i];
    }
    return _rolled_E;
  }

  //*******************************************************************************************
  VectorXd StrainConverter::unrolled_strain_metric(Eigen::Ref<const Matrix3d> const &F) const {
    return unroll_E(strain_metric(F));
  }

  //*******************************************************************************************
  Matrix3d StrainConverter::unrolled_strain_metric_to_F(Eigen::Ref<const VectorXd> const &E) const {
    return strain_metric_to_F(rollup_E(E));
  }

  //*******************************************************************************************
  //Calculates a linear combination of the components of unroll_E
  //using the sop_transf_mat
  VectorXd StrainConverter::sop(Matrix3d &E, Matrix3d &C, Matrix3d &U, Eigen::Ref<const Matrix3d> const &F, STRAIN_METRIC MODE) const {
    U = right_stretch_tensor(C, F);
    E = strain_metric(F, MODE);
    VectorXd _unroll_E = unroll_E(E);
    VectorXd _sop = m_sop_transf_mat * _unroll_E;
    return _sop;
  }

  //*******************************************************************************************

  VectorXd StrainConverter::sop(Matrix3d &E, Matrix3d &C, Matrix3d &U, Eigen::Ref<const Matrix3d> const &F) const {
    return sop(E, C, U, F, STRAIN_METRIC_MODE);
  }

  //************************************* SET routines ****************************************

  void StrainConverter::set_mode(const std::string &mode_name) {

    /// GREEN_LAGRANGE = 1/2 * (F^{T} F - I)
    if(mode_name == "STRAIN_GL" || mode_name == "GL") {
      STRAIN_METRIC_MODE = STRAIN_METRIC::GREEN_LAGRANGE;
      set_conventional_order_symmetric();
      curr_metric_func = &StrainConverter::green_lagrange;
      curr_inv_metric_func = &StrainConverter::green_lagrange_to_F;
    }
    /// BIOT = (U-I)
    else if(mode_name == "STRAIN_B" || mode_name == "B") {
      STRAIN_METRIC_MODE = STRAIN_METRIC::BIOT;
      set_conventional_order_symmetric();
      curr_metric_func = &StrainConverter::biot;
      curr_inv_metric_func = &StrainConverter::biot_to_F;
    }
    /// HENCKY = log(C)/2
    else if(mode_name == "STRAIN_H" || mode_name == "H") {
      STRAIN_METRIC_MODE = STRAIN_METRIC::HENCKY;
      set_conventional_order_symmetric();
      curr_metric_func = &StrainConverter::hencky;
      curr_inv_metric_func = &StrainConverter::hencky_to_F;
    }
    /// EULER_ALMANSI = 0.5 * (I-(F F^{T})^(-1))
    else if(mode_name == "STRAIN_EA" || mode_name == "EA") {
      STRAIN_METRIC_MODE = STRAIN_METRIC::EULER_ALMANSI;
      set_conventional_order_symmetric();
      curr_metric_func = &StrainConverter::euler_almansi;
      curr_inv_metric_func = &StrainConverter::euler_almansi_to_F;
    }
    /// STRETCH = U
    else if(mode_name == "STRAIN_U" || mode_name == "U") {
      STRAIN_METRIC_MODE = STRAIN_METRIC::STRETCH;
      set_conventional_order_symmetric();
      curr_metric_func = &StrainConverter::right_stretch_tensor;
      curr_inv_metric_func = &StrainConverter::disp_grad_to_F;
    }
    /// DISP_GRAD = F
    else if(mode_name == "STRAIN_F" || mode_name == "F") {
      STRAIN_METRIC_MODE = STRAIN_METRIC::DISP_GRAD;
      set_conventional_order_unsymmetric();
      curr_metric_func = &StrainConverter::disp_grad;
      curr_inv_metric_func = &StrainConverter::disp_grad_to_F;
    }
    else {
      std::cerr << "CRITICAL ERROR: In StrainConverter::set_mode(). Unrecognized strain metric '" << mode_name << "'" << std::endl;
      std::cerr << "                Your only options are STRAIN_GL, STRAIN_B, STRAIN_H, STRAIN_EA, STRAIN_U, and STRAIN_F" << std::endl;
      std::cerr << "                Exiting..." << std::endl;
      exit(1);
    }
  }


  //*******************************************************************************************

  void StrainConverter::set_symmetrized_sop(const SymGroup &pg) {
    Eigen::VectorXd tvec(Eigen::VectorXd::Zero(m_order_strain.size()));
    Eigen::Matrix3d tmat;
    SymGroupRep strain_rep(pg);

    for(Index g = 0; g < pg.size(); g++) {
      Eigen::MatrixXd trep(tvec.size(), tvec.size());
      for(Index i = 0; i < tvec.size(); i++) {
        tvec[i] = 1;
        tmat = rollup_E(tvec);
        tvec[i] = 0;
        tmat = pg[g].matrix() * tmat * pg[g].matrix().transpose();
        trep.col(i) = unroll_E(tmat);
      }
      strain_rep.set_rep(pg[g], SymMatrixXd(trep));
    }
    Eigen::VectorXd breathing = unroll_E(Eigen::Matrix3d::Identity());
    breathing.normalize();

    Eigen::MatrixXd bigmat(tvec.size(), tvec.size() + 1);
    bigmat.col(0) = breathing;
    bigmat.rightCols(tvec.size()) = irrep_trans_mat(strain_rep, pg).transpose();
    //m_sop_transf_mat=irrep_trans_mat(strain_rep,pg).transpose();
    Eigen::HouseholderQR<Eigen::MatrixXd> QR(bigmat);
    m_sop_transf_mat = QR.householderQ();
    m_symrep_ID = pg.master_group().add_representation(coord_transformed_copy(strain_rep, m_sop_transf_mat.transpose()));
  }

  //*******************************************************************************************
  //Conventional strain order parameters:
  //   e1 = (E11+E22+E33) / sqrt(3)
  //   e2 = (E11-E22) / sqrt(2)
  //   e3 = (2E33-E11+E22) / sqrt(6)
  //   e4 =  E12
  //   e5 =  E23
  //   e6 =  E13
  void StrainConverter::set_conventional_sop_transf_mat() {
    m_sop_transf_mat.resize(6, 6);
    m_sop_transf_mat << 1.0 / sqrt(3.0), 1.0 / sqrt(3.0), 1.0 / sqrt(3.0), 0, 0, 0,
                     1.0 / sqrt(2.0), -1.0 / sqrt(2.0), 0, 0, 0, 0,
                     -1.0 / sqrt(6.0), -1.0 / sqrt(6.0), 2.0 / sqrt(6.0), 0, 0, 0,
                     0, 0, 0, 0, 0, 1,
                     0, 0, 0, 0, 1, 0,
                     0, 0, 0, 1, 0, 0;
  }

  //*******************************************************************************************
  //Conventional order_strain:
  // unroll_E = (E11 E22 E33 E23 E13 E12)
  void StrainConverter::set_conventional_order_symmetric() {
    m_order_strain.resize(6, {0, 0});

    m_order_strain[0][0] = 0;
    m_order_strain[0][1] = 0;
    m_order_strain[1][0] = 1;
    m_order_strain[1][1] = 1;
    m_order_strain[2][0] = 2;
    m_order_strain[2][1] = 2;
    m_order_strain[3][0] = 1;
    m_order_strain[3][1] = 2;
    m_order_strain[4][0] = 0;
    m_order_strain[4][1] = 2;
    m_order_strain[5][0] = 0;
    m_order_strain[5][1] = 1;

    m_weight_strain.setOnes(6);
    m_weight_strain.tail(3) *= sqrt(2.0);
  }


  //*******************************************************************************************
  //Conventional order_strain:
  // unroll_E = (E11 E12 E13 E21 E22 E23 E31 E32 E33)
  void StrainConverter::set_conventional_order_unsymmetric() {
    m_order_strain.resize(9, {0, 0});
    m_weight_strain.setOnes(9);

    Index l = 0;
    for(Index i = 0; i < 3; i++) {
      for(Index j = 0; j < 3; j++) {
        m_order_strain[l][0] = i;
        m_order_strain[l][1] = j;
        l++;
      }
    }
  }

}
