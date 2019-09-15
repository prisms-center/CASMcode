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
  class SymRepBuilder {
  public:
    SymRepBuilder(std::string const &_name) :
      m_name(_name) {}

    std::string const &name() {
      return m_name;
    }

    virtual ~SymRepBuilder() {};

    virtual Eigen::MatrixXd symop_to_matrix(Eigen::Ref<const Eigen::Matrix3d> const &_matrix, bool time_reversal) const = 0;

    std::unique_ptr<SymRepBuilder> clone() const {
      return std::unique_ptr<SymRepBuilder>(_clone());
    }

  private:
    virtual SymRepBuilder *_clone() const = 0;

    std::string m_name;
  };

  /// In future, may include function pointers (wrapped in std::function<>) for controlling certain parts
  /// of program execution
  class AnisoValTraits {
  public:
    static const unsigned char LOCAL = 0;
    static const unsigned char GLOBAL = (1u << 0);
    static const unsigned char UNIT_LENGTH = (1u << 1);
    static const unsigned char TIME_REVERSAL = (1u << 2);
    static const unsigned char DESCRIBES_ORIENTATION = (1u << 3);

    AnisoValTraits(std::string const &_type_name,
                   std::vector<std::string> const &_std_var_names,
                   unsigned char _options,
                   SymRepBuilder const *_symrep_builder = nullptr,
                   std::set<std::string> const &_incompatible = {},
                   std::set<std::string> const &_must_apply_before = {},
                   std::set<std::string> const &_must_apply_after = {}) :
      m_type_name(_type_name),
      m_standard_var_names(_std_var_names),
      m_opt(_options),
      m_symrep_builder(_symrep_builder),
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
        return m_symrep_builder->symop_to_matrix(_matrix, time_reversal);
      }
      //else
      return Eigen::MatrixXd::Identity(dim(), dim());
    }

    /// \brief const access of type_name
    std::string const &type_name() const {
      return m_type_name;
    }

    /// \brief const access of type_name
    std::string const &name() const {
      return m_type_name;
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
      return m_opt & TIME_REVERSAL;
    }

    /// \brief returns true if value tracks the orientation of an occupying molecule (not typical)
    bool describes_occupant_orientation() const {
      return m_opt & DESCRIBES_ORIENTATION;
    }

    /// \brief conventional dimensionality of this DoF, returns -1 if always variable
    Index dim() const {
      return standard_var_names().size();
    }

    /// \brief equality comparison of type_name
    bool operator==(std::string const &other_name) const {
      return type_name() == other_name;
    }

    /// \brief lexicographic comparison of type_name
    bool operator<(std::string const &other_name) const {
      return type_name() < other_name;
    }

    /// \brief comparison of type_name, domain (discrete/continuous) and mode (local/global)
    bool identical(AnisoValTraits const &other) const {
      return type_name() == other.type_name()
             && m_opt == other.m_opt;
    }

    /// \brief allow implicit conversion to std::string (type_name)
    operator std::string const &() const {
      return type_name();
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
    std::string m_type_name;
    std::vector<std::string> m_standard_var_names;
    unsigned char m_opt;
    SymRepBuilder const *m_symrep_builder;
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
