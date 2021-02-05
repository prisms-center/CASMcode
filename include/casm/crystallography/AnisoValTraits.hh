#ifndef CASM_AnisoValTraits
#define CASM_AnisoValTraits

#include <memory>
#include <set>
#include <string>
#include <vector>

#include "casm/crystallography/SymRepBuilder.hh"
#include "casm/global/definitions.hh"
#include "casm/global/eigen.hh"

// Defines notstd::make_unique<>()
#include "casm/misc/ParsingDictionary.hh"
#include "casm/misc/cloneable_ptr.hh"

namespace CASM {

class AnisoValTraits;

template <>
ParsingDictionary<AnisoValTraits> make_parsing_dictionary<AnisoValTraits>();

/// \brief Specifies traits of (possibly) anisotropic crystal properties
class AnisoValTraits {
 public:
  static const unsigned char LOCAL = 0;
  static const unsigned char GLOBAL = (1u << 0);
  static const unsigned char UNIT_LENGTH = (1u << 1);
  static const unsigned char DESCRIBES_ORIENTATION = (1u << 2);
  static const unsigned char EXTENSIVE = (1u << 3);

  /// Named constructor for uninitialized AnisoValTraits
  static AnisoValTraits null();

  /// Named constructor for discrete site occupation AnisoValTraits
  static AnisoValTraits occ();

  /// Named constructor for total energy AnisoValTraits
  static AnisoValTraits energy();

  /// Named constructor for mapping cost AnisoValTraits
  static AnisoValTraits cost();

  /// Named constructor for selective_dynamics AnisoValTraits
  static AnisoValTraits selective_dynamics();

  /// Named constructor for site displacement AnisoValTraits
  static AnisoValTraits disp();

  /// Named constructor for site coordinate AnisoValTraits
  static AnisoValTraits coordinate();

  /// Named constructor for lattice vector AnisoValTraits
  static AnisoValTraits latvec();

  /// Named constructor for site force AnisoValTraits
  static AnisoValTraits force();

  /// Named constructor for global 'isometry' AnisoValTraits
  static AnisoValTraits isometry();

  /// Named constructor for global strain AnisoValTraits
  static AnisoValTraits strain(std::string const &_metric);

  /// Non-collinear magnetic spin, with spin-orbit coupling
  static AnisoValTraits SOmagspin();

  /// \brief Non-collinear magnetic spin, with spin-orbit coupling, constrained
  /// to unit length
  static AnisoValTraits SOunitmagspin();

  /// Non-collinear magnetic spin, without spin-orbit coupling
  static AnisoValTraits NCmagspin();

  /// \brief Non-collinear magnetic spin, without spin-orbit coupling,
  /// constrained to unit length
  static AnisoValTraits NCunitmagspin();

  /// Collinear magnetic spin
  static AnisoValTraits Cmagspin();

  /// Collinear magnetic spin, constrained to unit length
  static AnisoValTraits Cunitmagspin();

  /// Named constructor for d-orbital occupation AnisoValTraits
  static AnisoValTraits d_orbital_occupation();

  /// Named constructor for spin-polarized d-orbital occupation AnisoValTraits
  static AnisoValTraits d_orbital_occupation_spin_polarized();

  /// Returns string after final `delim` character
  ///
  /// Given a string, returns string with all characters before the final
  /// `delim` character deleted. For example, if a trajectory has properties
  /// with `key1="step1_force"`, `key2="step2_force"`, etc., then
  /// AnisoValTraits::name_suffix(key1) will return "force" if `delim='_'`.
  static std::string name_suffix(std::string const &_name, char delim = '_') {
    std::string result;
    for (char c : _name) {
      if (c == delim)
        result.clear();
      else
        result.push_back(c);
    }
    return result;
  }

  /// Explicit constructor for AnisoValTraits
  AnisoValTraits(
      std::string const &_name, std::vector<std::string> const &_std_var_names,
      unsigned char _options,
      SymRepBuilderInterface const &_symrep_builder = NullSymRepBuilder(),
      std::set<std::string> const &_incompatible = {},
      std::set<std::string> const &_must_apply_before = {},
      std::set<std::string> const &_must_apply_after = {},
      std::vector<std::string> const &_variable_descriptions = {},
      bool _default = false);

  /// Construct a copy of an existing AnisoValTraits with matching name suffix
  AnisoValTraits(std::string const &_name);

  /// Returns "AnisoValTraits"
  static std::string class_desc() { return "AnisoValTraits"; }

  /// Generate a symmetry representation for the supporting vector space
  Eigen::MatrixXd symop_to_matrix(
      Eigen::Ref<const Eigen::Matrix3d> const &_matrix,
      Eigen::Ref<const Eigen::Vector3d> const &_tau, bool time_reversal) const {
    if (m_symrep_builder != nullptr) {
      return m_symrep_builder->symop_to_matrix(_matrix, _tau, time_reversal,
                                               dim());
    }
    // else
    return Eigen::MatrixXd::Identity(dim(), dim());
  }

  /// Const access of name
  std::string const &name() const { return m_name; }

  /// \brief Return true if *this has 'default' designation, meaning it can be
  /// overridden
  bool is_default() const { return m_default; }

  /// Return 'options' bitflag
  unsigned char options() const { return m_opt; }

  /// Returns true if type is extensive
  bool extensive() const { return m_opt & EXTENSIVE; }

  // Please use !AnisoValTraits::global() instead of implementing
  // AnisoValTraits::local()
  // bool local()const {
  //  return !global();
  //}

  /// Returns true if type is global
  bool global() const { return m_opt & GLOBAL; }

  /// Returns true if type must always have unit length
  bool unit_length() const { return m_opt & UNIT_LENGTH; }

  /// Returns true if time-reversal changes the value
  bool time_reversal_active() const {
    return m_symrep_builder->time_reversal_active();
  }

  /// Returns true if value tracks the orientation of an occupying
  /// molecule (not typical)
  bool describes_occupant_orientation() const {
    return m_opt & DESCRIBES_ORIENTATION;
  }

  /// Conventional dimensionality of this type, returns -1 if always variable
  Index dim() const { return standard_var_names().size(); }

  /// Equality comparison of name
  bool operator==(std::string const &other_name) const {
    return name() == other_name;
  }

  /// Lexicographic comparison of name
  bool operator<(std::string const &other_name) const {
    return name() < other_name;
  }

  /// Allow implicit conversion to std::string (name)
  operator std::string const &() const { return name(); }

  /// Return standard coordinate axes for continuous variable space
  std::vector<std::string> const &standard_var_names() const {
    return m_standard_var_names;
  }

  /// Returns expanded description of each standard_var_name
  std::vector<std::string> const &variable_descriptions() const {
    return m_variable_descriptions;
  }

  /// \brief When used as a DoF, this is a set of DoFs that are incompatible
  /// with this type of DoF
  std::set<std::string> const &incompatible() const { return m_incompatible; }

  /// \brief When used as a DoF, this is a set of DoFs that *must* be applied
  /// before this DoF is applied
  std::set<std::string> const &must_apply_before() const {
    return m_apply_before;
  }

  /// \brief When used as a DoF, this is a set of DoFs that *must* be applied
  /// after this DoF is applied
  std::set<std::string> const &must_apply_after() const {
    return m_apply_after;
  }

  /// Clone this object
  std::unique_ptr<AnisoValTraits> clone() const {
    return notstd::make_unique<AnisoValTraits>(*this);
  }

  /// Return name of the SymRepBuilderInterface
  std::string symrep_builder_name() const {
    if (m_symrep_builder)
      return m_symrep_builder->name();
    else
      return "NULL";
  }

 protected:
  std::string m_name;
  bool m_default;
  std::vector<std::string> m_standard_var_names;
  std::vector<std::string> m_variable_descriptions;
  unsigned char m_opt;
  notstd::cloneable_ptr<const SymRepBuilderInterface> m_symrep_builder;
  std::set<std::string> m_incompatible;
  std::set<std::string> m_apply_before;
  std::set<std::string> m_apply_after;
};

/// Check for equivalence of AnisoValTraits
///
/// Compares: name(), standard_var_names(), options(), symrep_builder_name(),
/// incompatible(), must_apply_before(), must_apply_after(), and is_default()
inline bool identical(AnisoValTraits const &A, AnisoValTraits const &B) {
  return (A.name() == B.name() && A.is_default() == B.is_default() &&
          A.symrep_builder_name() == B.symrep_builder_name() &&
          A.standard_var_names() == B.standard_var_names() &&
          A.options() == B.options() &&
          A.must_apply_before() == B.must_apply_before() &&
          A.must_apply_after() == B.must_apply_after() &&
          A.incompatible() == B.incompatible());
}

}  // namespace CASM

#endif
