#ifndef CASM_AnisoValTraits
#define CASM_AnisoValTraits

#include <set>
#include <string>
#include <memory>
#include <vector>
#include "include/casm/CASM_global_definitions.hh"
#include "include/casm/CASM_global_Eigen.hh"

//Defines notstd::make_unique<>()
#include "include/casm/misc/cloneable_ptr.hh"

namespace CASM {

  /// \brief Abstract base class that provides interface for converting cartesian isometry to specialized transformation matrix
  class SymRepBuilderInterface {
  public:
    SymRepBuilderInterface(std::string const &_name,
                           bool _time_reversal) :
      m_name(_name),
      m_time_reversal(_time_reversal) {

    }

    std::string const &name() const {
      return m_name;
    }

    bool time_reversal_active() const {
      return m_time_reversal;
    }

    virtual ~SymRepBuilderInterface() {};

    virtual Eigen::MatrixXd symop_to_matrix(Eigen::Ref<const Eigen::Matrix3d> const &_matrix, bool time_reversal, Index dim) const = 0;

    std::unique_ptr<SymRepBuilderInterface> clone() const {
      return std::unique_ptr<SymRepBuilderInterface>(_clone());
    }

  private:
    virtual SymRepBuilderInterface *_clone() const = 0;

    std::string m_name;
    bool m_time_reversal;
  };

  template<bool uses_time_reversal>
  class TemplateSymRepBuilderBase : public SymRepBuilderInterface {
  protected:
    TemplateSymRepBuilderBase(std::string const &_name) :
      SymRepBuilderInterface(_name, uses_time_reversal) {}

  };

  using SymRepBuilderBase = TemplateSymRepBuilderBase<false>;

  using TimeReversalSymRepBuilderBase = TemplateSymRepBuilderBase<true>;

  class NullSymRepBuilder : public SymRepBuilderBase {
  public:
    NullSymRepBuilder() :
      SymRepBuilderBase("Null") {}

    Eigen::MatrixXd symop_to_matrix(Eigen::Ref<const Eigen::Matrix3d> const &_matrix, bool time_reversal, Index dim) const override {
      return Eigen::MatrixXd();
    }
  private:
    SymRepBuilderInterface *_clone() const override {
      return nullptr;
    }
  };

  class CartesianSymRepBuilder : public SymRepBuilderBase {
  public:
    CartesianSymRepBuilder() :
      SymRepBuilderBase("Cartesian") {}

    Eigen::MatrixXd symop_to_matrix(Eigen::Ref<const Eigen::Matrix3d> const &_matrix, bool time_reversal, Index dim) const override {
      return _matrix;
    }
  private:
    SymRepBuilderInterface *_clone() const override {
      return new CartesianSymRepBuilder();
    }
  };

  class AngularMomentumSymRepBuilder : public TimeReversalSymRepBuilderBase {
  public:
    AngularMomentumSymRepBuilder() :
      TimeReversalSymRepBuilderBase("AngularMomentum") {}

    Eigen::MatrixXd symop_to_matrix(Eigen::Ref<const Eigen::Matrix3d> const &_matrix, bool time_reversal, Index dim) const override {
      return ((time_reversal ? -1. : 1.) * _matrix.determinant()) * _matrix;
    }
  private:
    SymRepBuilderInterface *_clone() const override {
      return new AngularMomentumSymRepBuilder();
    }
  };

  class TimeReversalSymRepBuilder : public TimeReversalSymRepBuilderBase {
  public:
    TimeReversalSymRepBuilder() :
      TimeReversalSymRepBuilderBase("TimeReversal") {}

    Eigen::MatrixXd symop_to_matrix(Eigen::Ref<const Eigen::Matrix3d> const &_matrix, bool time_reversal, Index dim) const override {
      return -Eigen::MatrixXd::Identity(dim, dim);
    }
  private:
    SymRepBuilderInterface *_clone() const override {
      return new TimeReversalSymRepBuilder();
    }
  };

  class Rank2TensorSymRepBuilder : public SymRepBuilderBase {
  public:
    Rank2TensorSymRepBuilder() :
      SymRepBuilderBase("Rank2Tensor") {}

    Eigen::MatrixXd symop_to_matrix(Eigen::Ref<const Eigen::Matrix3d> const &S, bool time_reversal, Index dim) const override {
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

  /// In future, may include function pointers (wrapped in std::function<>) for controlling certain parts
  /// of program execution
  class AnisoValTraits {
  public:
    static const unsigned char LOCAL = 0;
    static const unsigned char GLOBAL = (1u << 0);
    static const unsigned char UNIT_LENGTH = (1u << 1);
    static const unsigned char DESCRIBES_ORIENTATION = (1u << 2);

    static AnisoValTraits disp();
    static AnisoValTraits force();
    static AnisoValTraits strain(std::string const &_prefix);
    static AnisoValTraits magspin();
    static AnisoValTraits magmom();

    AnisoValTraits(std::string const &_name,
                   std::vector<std::string> const &_std_var_names,
                   unsigned char _options,
                   SymRepBuilderInterface const &_symrep_builder = NullSymRepBuilder(),
                   std::set<std::string> const &_incompatible = {},
                   std::set<std::string> const &_must_apply_before = {},
                   std::set<std::string> const &_must_apply_after = {}) :
      m_name(_name),
      m_standard_var_names(_std_var_names),
      m_opt(_options),
      m_symrep_builder(_symrep_builder.clone()),
      m_incompatible(_incompatible),
      m_apply_before(_must_apply_before),
      m_apply_after(_must_apply_after) {

    }

    static std::string class_desc() {
      return "AnisoValTraits";
    }

    /// \brief Generate a symmetry representation for the supporting vector space
    Eigen::MatrixXd symop_to_matrix(Eigen::Ref<const Eigen::Matrix3d> const &_matrix, bool time_reversal) const {
      if(m_symrep_builder != nullptr) {
        return m_symrep_builder->symop_to_matrix(_matrix, time_reversal, dim());
      }
      //else
      return Eigen::MatrixXd::Identity(dim(), dim());
    }

    /// \brief const access of name
    std::string const &name() const {
      return m_name;
    }

    // Please use !AnisoValTraits::global() instead of implementing AnisoValTraits::local()
    //bool local()const {
    //  return !global();
    //}

    /// \brief returns true if DoF is global
    bool global()const {
      return m_opt & GLOBAL;
    }

    /// \brief returns true if DoF must always have unit length
    bool unit_length()const {
      return m_opt & UNIT_LENGTH;
    }

    /// \brief returns true if time-reversal changes the DoF value
    bool time_reversal_active() const {
      return m_symrep_builder->time_reversal_active();
    }

    /// \brief returns true if value tracks the orientation of an occupying molecule (not typical)
    bool describes_occupant_orientation() const {
      return m_opt & DESCRIBES_ORIENTATION;
    }

    /// \brief conventional dimensionality of this DoF, returns -1 if always variable
    Index dim() const {
      return standard_var_names().size();
    }

    /// \brief equality comparison of name
    bool operator==(std::string const &other_name) const {
      return name() == other_name;
    }

    /// \brief lexicographic comparison of name
    bool operator<(std::string const &other_name) const {
      return name() < other_name;
    }

    /// \brief comparison of name, domain (discrete/continuous) and mode (local/global)
    bool identical(AnisoValTraits const &other) const {
      return name() == other.name()
             && m_opt == other.m_opt;
    }

    /// \brief allow implicit conversion to std::string (name)
    operator std::string const &() const {
      return name();
    }

    /// \brief return standard coordinate axes for continuous variable space
    std::vector<std::string> const &standard_var_names() const {
      return m_standard_var_names;
    }

    std::set<std::string> const &incompatible() const {
      return m_incompatible;
    }

    /// \brief Return list of DoFs that *must* be applied before this DoF is applied
    std::set<std::string> const &must_apply_before() const {
      return m_apply_before;
    }

    /// \brief Return list of DoFs that *must* be applied after this DoF is applied
    std::set<std::string> const &must_apply_after() const {
      return m_apply_after;
    }

    std::unique_ptr<AnisoValTraits> clone() const {
      return notstd::make_unique<AnisoValTraits>(*this);
    }

  protected:
    std::string m_name;
    std::vector<std::string> m_standard_var_names;
    unsigned char m_opt;
    notstd::cloneable_ptr<const SymRepBuilderInterface> m_symrep_builder;
    std::set<std::string> m_incompatible;
    std::set<std::string> m_apply_before;
    std::set<std::string> m_apply_after;
  };


  namespace AnisoVal_impl {
    /// \brief A class to manage dynamic evaluation of BasisFunctions
    //struct TraitsConverter {
    /*
    inline
    notstd::cloneable_ptr<DoFType::BasicTraits> traits2cloneable_ptr(const DoFType::BasicTraits &value) {
      return notstd::cloneable_ptr<DoFType::BasicTraits>(value.clone().release());
    }
    */
    //};
  }

}

#endif
