#ifndef CASM_SymRepBuilder
#define CASM_SymRepBuilder

#include <string>
#include <memory>
#include "casm/global/definitions.hh"
#include "casm/global/eigen.hh"

namespace CASM {

  /// \brief Abstract base class that provides interface for converting cartesian isometry to specialized transformation matrix
  /// Given a symmetry operation as a 3x3 orthogonal matrix acting on Cartesian space, a 3x1 Cartesian translation, and a time-reversal flag
  /// (true indicates +t goes to -t and false indicates +t goes to +t), returns a particular NxN matrix transformation whose dimension and
  /// value depends on the underlying implementation.
  class SymRepBuilderInterface {
  public:
    /// \brief Interace class manages the name of the SymRepBuilder and whether it is affected by time-reversal (i.e., time_reversal_active() is true)
    SymRepBuilderInterface(std::string const &_name,
                           bool _time_reversal_active) :
      m_name(_name),
      m_time_reversal_active(_time_reversal_active) {

    }

    /// \brief Returns name of SymRepBuilder type
    std::string const &name() const {
      return m_name;
    }

    /// \brief Returns true if symmetry representation is affected by time-reversal
    bool time_reversal_active() const {
      return m_time_reversal_active;
    }

    /// \brief Virtual destructor allows deletion of derived classes through pointer to interface
    virtual ~SymRepBuilderInterface() {};

    /// \brief Given the 3x3 rotation/rotoreflection matrix, translation vector 'tau', and time_reversal operator of Cartesian SymOp,
    /// constructs a dim x dim matrix representation of the symop
    /// Derived implementations may require that dim have a specific value, or fall within a range of allowed values.
    virtual Eigen::MatrixXd symop_to_matrix(Eigen::Ref<const Eigen::Matrix3d> const &_matrix,
                                            Eigen::Ref<const Eigen::Vector3d> const &_tau,
                                            bool time_reversal,
                                            Index dim) const = 0;

    /// \brief Given the 3x3 rotation/rotoreflection matrix, translation vector 'tau', and time_reversal operator of Cartesian SymOp,
    /// constructs a dim x dim complexmatrix representation of the symop
    /// Derived implementations may require that dim have a specific value, or fall within a range of allowed values.
    virtual Eigen::MatrixXcd symop_to_complex_matrix(Eigen::Ref<const Eigen::Matrix3d> const &_matrix,
                                                     Eigen::Ref<const Eigen::Vector3d> const &_tau,
                                                     bool time_reversal,
                                                     Index dim) const {
      return symop_to_matrix(_matrix, _tau, time_reversal, dim).cast<std::complex<double> >();
    }

    std::unique_ptr<SymRepBuilderInterface> clone() const {
      return std::unique_ptr<SymRepBuilderInterface>(_clone());
    }

  private:
    virtual SymRepBuilderInterface *_clone() const = 0;

    std::string m_name;
    bool m_time_reversal_active;
  };

  template<bool uses_time_reversal>
  class TemplateSymRepBuilderBase : public SymRepBuilderInterface {
  protected:
    TemplateSymRepBuilderBase(std::string const &_name) :
      SymRepBuilderInterface(_name, uses_time_reversal) {}

  };


  using SymRepBuilderBase = TemplateSymRepBuilderBase<false>;

  using TimeReversalSymRepBuilderBase = TemplateSymRepBuilderBase<true>;

  /// \brief Un-cloneable class for specifying absence of valid SymRepBuilder
  class NullSymRepBuilder : public SymRepBuilderBase {
  public:
    NullSymRepBuilder() :
      SymRepBuilderBase("NULL") {}

    Eigen::MatrixXd symop_to_matrix(Eigen::Ref<const Eigen::Matrix3d> const &_matrix,
                                    Eigen::Ref<const Eigen::Vector3d> const &_tau,
                                    bool time_reversal,
                                    Index dim) const override {
      return Eigen::MatrixXd();
    }
  private:
    SymRepBuilderInterface *_clone() const override {
      return nullptr;
    }
  };

  /// \brief Builds symmetry representation as the Cartesian matrix of povided SymOp
  class CartesianSymRepBuilder : public SymRepBuilderBase {
  public:
    CartesianSymRepBuilder() :
      SymRepBuilderBase("Cartesian") {}

    Eigen::MatrixXd symop_to_matrix(Eigen::Ref<const Eigen::Matrix3d> const &_matrix,
                                    Eigen::Ref<const Eigen::Vector3d> const &_tau,
                                    bool time_reversal,
                                    Index dim) const override {
      return _matrix;
    }
  private:
    SymRepBuilderInterface *_clone() const override {
      return new CartesianSymRepBuilder();
    }
  };

  /// \brief Builds symmetry representation as the 'dim' x 'dim' identity matrix, regardless of symop
  class IdentitySymRepBuilder : public SymRepBuilderBase {
  public:
    IdentitySymRepBuilder() :
      SymRepBuilderBase("Identity") {}

    Eigen::MatrixXd symop_to_matrix(Eigen::Ref<const Eigen::Matrix3d> const &_matrix,
                                    Eigen::Ref<const Eigen::Vector3d> const &_tau,
                                    bool time_reversal,
                                    Index dim) const override {
      return Eigen::MatrixXd::Identity(dim, dim);
    }
  private:
    SymRepBuilderInterface *_clone() const override {
      return new IdentitySymRepBuilder();
    }
  };

  /// \brief Builds symmetry representation that is the angular momentum symmetry representation of provided symop
  class AngularMomentumSymRepBuilder : public TimeReversalSymRepBuilderBase {
  public:
    AngularMomentumSymRepBuilder() :
      TimeReversalSymRepBuilderBase("AngularMomentum") {}

    Eigen::MatrixXd symop_to_matrix(Eigen::Ref<const Eigen::Matrix3d> const &_matrix,
                                    Eigen::Ref<const Eigen::Vector3d> const &_tau,
                                    bool time_reversal,
                                    Index dim) const override {
      return ((time_reversal ? -1. : 1.) * _matrix.determinant()) * _matrix;
    }
  private:
    SymRepBuilderInterface *_clone() const override {
      return new AngularMomentumSymRepBuilder();
    }
  };

  /// \brief Builds symmetry representation that is 'dim'x'dim' +Identity (-Identity) matrix if time_reversal is false (true)
  class TimeReversalSymRepBuilder : public TimeReversalSymRepBuilderBase {
  public:
    TimeReversalSymRepBuilder() :
      TimeReversalSymRepBuilderBase("TimeReversal") {}

    Eigen::MatrixXd symop_to_matrix(Eigen::Ref<const Eigen::Matrix3d> const &_matrix,
                                    Eigen::Ref<const Eigen::Vector3d> const &_tau,
                                    bool time_reversal,
                                    Index dim) const override {
      return (time_reversal ? -1. : 1.) * Eigen::MatrixXd::Identity(dim, dim);
    }
  private:
    SymRepBuilderInterface *_clone() const override {
      return new TimeReversalSymRepBuilder();
    }
  };

  /// \brief Build 6x6 symmetry representation for a rank 2 Cartesian tensor represented in Kelvin notation
  class Rank2TensorSymRepBuilder : public SymRepBuilderBase {
  public:
    Rank2TensorSymRepBuilder() :
      SymRepBuilderBase("Rank2Tensor") {}

    Eigen::MatrixXd symop_to_matrix(Eigen::Ref<const Eigen::Matrix3d> const &S,
                                    Eigen::Ref<const Eigen::Vector3d> const &_tau,
                                    bool time_reversal,
                                    Index dim) const override {
      Eigen::MatrixXd result(6, 6);

      result <<
             S(0, 0)*S(0, 0), S(0, 1)*S(0, 1), S(0, 2)*S(0, 2), sqrt(2)*S(0, 1)*S(0, 2), sqrt(2)*S(0, 2)*S(0, 0), sqrt(2)*S(0, 0)*S(0, 1),
             S(1, 0)*S(1, 0), S(1, 1)*S(1, 1), S(1, 2)*S(1, 2), sqrt(2)*S(1, 1)*S(1, 2), sqrt(2)*S(1, 2)*S(1, 0), sqrt(2)*S(1, 0)*S(1, 1),
             S(2, 0)*S(2, 0), S(2, 1)*S(2, 1), S(2, 2)*S(2, 2), sqrt(2)*S(2, 1)*S(2, 2), sqrt(2)*S(2, 2)*S(2, 0), sqrt(2)*S(2, 0)*S(2, 1),
             sqrt(2)*S(1, 0)*S(2, 0), sqrt(2)*S(1, 1)*S(2, 1), sqrt(2)*S(1, 2)*S(2, 2), S(1, 1)*S(2, 2) + S(1, 2)*S(2, 1), S(1, 0)*S(2, 2) + S(1, 2)*S(2, 0), S(1, 1)*S(2, 0) + S(1, 0)*S(2, 1),
             sqrt(2)*S(2, 0)*S(0, 0), sqrt(2)*S(2, 1)*S(0, 1), sqrt(2)*S(2, 2)*S(0, 2), S(2, 1)*S(0, 2) + S(2, 2)*S(0, 1), S(2, 0)*S(0, 2) + S(2, 2)*S(0, 0), S(2, 1)*S(0, 0) + S(2, 0)*S(0, 1),
             sqrt(2)*S(0, 0)*S(1, 0), sqrt(2)*S(0, 1)*S(1, 1), sqrt(2)*S(0, 2)*S(1, 2), S(0, 1)*S(1, 2) + S(0, 2)*S(1, 1), S(0, 0)*S(1, 2) + S(0, 2)*S(1, 0), S(0, 1)*S(1, 0) + S(0, 0)*S(1, 1);
      return result;

    }
  private:
    SymRepBuilderInterface *_clone() const override {
      return new Rank2TensorSymRepBuilder();
    }
  };

  // \brief Build 30x30 symmetry representation for an unrolled pair of d-orbital occupation matrices
  class dOrbitalOccupationSymRepBuilder : public TimeReversalSymRepBuilderBase {
    public:
        dOrbitalOccupationSymRepBuilder() : TimeReversalSymRepBuilderBase("dOrbitalOccupation"){}
    
        Eigen::MatrixXd symop_to_matrix(Eigen::Ref<const Eigen::Matrix3d> const &_matrix,
                                            Eigen::Ref<const Eigen::Vector3d> const &_tau,
                                            bool time_reversal,
                                            Index dim) const override
        {return Eigen::MatrixXd::Identity(30,30);}

    private:
        SymRepBuilderInterface *_clone() const override {
         return new dOrbitalOccupationSymRepBuilder();
        }
  };

  // Named constructors for all previously defined SymRepBuilders
  namespace SymRepBuilder {
    inline
    CartesianSymRepBuilder Identity() {
      return CartesianSymRepBuilder();
    }

    inline
    CartesianSymRepBuilder Cartesian() {
      return CartesianSymRepBuilder();
    }

    inline
    AngularMomentumSymRepBuilder AngularMomentum() {
      return AngularMomentumSymRepBuilder();
    }

    inline
    TimeReversalSymRepBuilder TimeReversal() {
      return TimeReversalSymRepBuilder();
    }

    inline
    Rank2TensorSymRepBuilder Rank2Tensor() {
      return Rank2TensorSymRepBuilder();
    }
  }

}
#endif
