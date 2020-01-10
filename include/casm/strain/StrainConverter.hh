#ifndef STRAINCONVERTER_HH
#define STRAINCONVERTER_HH

#include <vector>
#include "casm/external/Eigen/Dense"
#include "casm/global/definitions.hh"
#include "casm/symmetry/SymGroupRepID.hh"


//--------------------------------------------------
// STRAIN CLASS
// NAMING CONVENTION:
//       - F: deformation tensor (F = R * U)
//       - E: strain metric (calculated based on what STRAIN_METRIC_MODE is)
//       - R: rotation tensor
//       - U: stretch tensor
//       - C: F^{T} * F
//       - unrolled_E: the 6 unique components of E ordered in an order
//                     specified by order_strain
//       - sop: strain order parameter (sop = sop_transf_mat * unrolled_E)

namespace CASM {
  class SymGroup;
  namespace SymRepTools {
    struct IrrepWedge;
  }

  typedef Eigen::VectorXd VectorXd;
  typedef Eigen::MatrixXd MatrixXd;
  typedef Eigen::Matrix3d Matrix3d;

  //Enum type that will define how the strain metrics used to
  //calculate strain order parameters are going to be calculated
  enum STRAIN_METRIC {GREEN_LAGRANGE = 0, BIOT = 1, HENCKY = 2, EULER_ALMANSI = 3, STRETCH = 4, DISP_GRAD = 5};


  class StrainConverter {
  public:
    // Static methods for conversion from F to various strain metrics

    /// GREEN_LAGRANGE = (C-I)/2
    static Matrix3d green_lagrange(Eigen::Ref<const Matrix3d> const &F);
    static Matrix3d green_lagrange_to_F(Eigen::Ref<const Matrix3d> const &E);

    /// BIOT = (U-I)
    static Matrix3d biot(Eigen::Ref<const Matrix3d> const &F);
    static Matrix3d biot_to_F(Eigen::Ref<const Matrix3d> const &B);

    /// HENCKY = log(C)/2
    static Matrix3d hencky(Eigen::Ref<const Matrix3d> const &F);
    static Matrix3d hencky_to_F(Eigen::Ref<const Matrix3d> const &H);

    /// EULER_ALMANSI = (I-(F F^{T})^(-1))/2
    static Matrix3d euler_almansi(Eigen::Ref<const Matrix3d> const &F);
    static Matrix3d euler_almansi_to_F(Eigen::Ref<const Matrix3d> const &A);

    /// DISP_GRAD = F
    static Matrix3d disp_grad(Eigen::Ref<const Matrix3d> const &F);
    static Matrix3d disp_grad_to_F(Eigen::Ref<const Matrix3d> const &F) {
      return F;
    }

    //-------------------------------------------------
    // Routines that calculate derived quantities given the
    // deformation tensor
    static Matrix3d metric_tensor(Eigen::Ref<const Matrix3d> const &F);
    static Matrix3d right_stretch_tensor(Matrix3d &C, Eigen::Ref<const Matrix3d> const &F);
    static Matrix3d right_stretch_tensor(Eigen::Ref<const Matrix3d> const &F);

    static Matrix3d strain_metric(Eigen::Ref<const Matrix3d> const &F, STRAIN_METRIC MODE);
    //-------------------------------------------------

    StrainConverter(bool override_flag = false) {
      if(!override_flag) {
        std::cerr << "WARNING in CASM::StrainConverter you are calling the default constructor" << std::endl;
        std::cerr << "This is going to initialize a \"default\" sop_transf_mat and order_strain" << std::endl;
        std::cerr << "The Green-Lagrange strain metric will be used for the calculation of all deformation metrics" << std::endl;
        std::cerr << "PLEASE USE WITH CAUTION" << std::endl;
      }
      STRAIN_METRIC_MODE = GREEN_LAGRANGE;
      curr_metric_func = &StrainConverter::green_lagrange;
      set_conventional_sop_transf_mat();
      set_conventional_order_symmetric();
    };

    StrainConverter(const STRAIN_METRIC &_MODE, const MatrixXd &_sop_transf_mat,
                    const std::vector<std::vector<Index> >  &_order_strain) :
      STRAIN_METRIC_MODE(_MODE), m_sop_transf_mat(_sop_transf_mat),
      m_order_strain(_order_strain) {
      if(_MODE == GREEN_LAGRANGE)
        curr_metric_func = &StrainConverter::green_lagrange;
      else if(_MODE == BIOT)
        curr_metric_func = &StrainConverter::biot;
      else if(_MODE == HENCKY)
        curr_metric_func = &StrainConverter::hencky;
      else if(_MODE == EULER_ALMANSI)
        curr_metric_func = &StrainConverter::euler_almansi;
      else if(_MODE == STRETCH)
        curr_metric_func = &StrainConverter::right_stretch_tensor;
      else if(_MODE == DISP_GRAD)
        curr_metric_func = &StrainConverter::disp_grad;
    }

    StrainConverter(const std::string &mode_name) :
      StrainConverter(true) {
      set_mode(mode_name);
    }

    Index dim() const {
      return m_order_strain.size();
    }

    //-------------------------------------------------
    /// Get the strain metric in the current mode
    Matrix3d strain_metric(Eigen::Ref<const Matrix3d> const &F) const;
    Matrix3d strain_metric_to_F(Eigen::Ref<const Matrix3d> const &E) const;

    /// Unrolls the green-lagrange metric ( or any symmetric metric)
    VectorXd unroll_E(Eigen::Ref<const Matrix3d> const &E) const;
    Matrix3d rollup_E(Eigen::Ref<const VectorXd> const &_unrolled_E) const;

    VectorXd unrolled_strain_metric(Eigen::Ref<const Matrix3d> const &F) const;
    Matrix3d unrolled_strain_metric_to_F(Eigen::Ref<const VectorXd> const &E) const;

    VectorXd sop(Matrix3d &E, Matrix3d &C, Matrix3d &U, Eigen::Ref<const Matrix3d> const &F) const;
    VectorXd sop(Matrix3d &E, Matrix3d &C, Matrix3d &U, Eigen::Ref<const Matrix3d> const &F, STRAIN_METRIC MODE) const;

    //--------------------------------------------------
    // Routines that set the internal parameters of the
    // class
    void set_mode(const std::string &mode_name);

    const Eigen::MatrixXd &sop_transf_mat() const {
      return m_sop_transf_mat;
    }

    SymGroupRepID symrep_ID() const {
      return m_symrep_ID;
    }

    std::vector<SymRepTools::IrrepWedge> irreducible_sop_wedges(const SymGroup &pg);

    std::vector<SymRepTools::IrrepWedge> irreducible_wedges(const SymGroup &pg);

    void set_symmetrized_sop(const SymGroup &pg);

    void set_conventional_sop_transf_mat();
    void set_conventional_order_symmetric();
    void set_conventional_order_unsymmetric();

    // Need to rewrite these routines to write out only the settings
    // jsonParser &to_json(const StrainConverter &strain, jsonParser &json);
    // void from_json(StrainConverter &strain, const jsonParser &json);

  private:
    STRAIN_METRIC STRAIN_METRIC_MODE; //set the mode when you initialize your PRIM
    MatrixXd m_sop_transf_mat; //Use as sop = sop_transf_mat * unrolled_E
    // unrolled_E[i]= m_weight_strain[i]*strain(m_order_strain[i][0], m_order_strain[i][1])
    std::vector<std::vector<Index> > m_order_strain; //lists the order of strains to list unrolled_E
    Eigen::VectorXd m_weight_strain; // weights for the elements of unrolled_E

    SymGroupRepID m_symrep_ID;
    // typedef MetricFuncPtr for method pointers that take displacement gradient tensor as argument
    typedef Matrix3d(*MetricFuncPtr)(Eigen::Ref<const Matrix3d> const &);
    MetricFuncPtr curr_metric_func;
    MetricFuncPtr curr_inv_metric_func;

  };
}
#endif
