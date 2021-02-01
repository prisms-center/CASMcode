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

/// Class for specifying essential traits of (possibly) anisotropic crystal
/// properties
///
/// CASM needs to know how to use crystal properties in several places:
/// - Site and global properties in xtal::SimpleStructure
/// - Degrees of freedoms in xtal::BasicStructure through xtal::DoFSet and
/// xtal::SiteDoFSet
/// - Molecule properties in xtal::SpeciesAttribute
/// - Calculated properties in MappedProperties
///
/// Examples of various property types, and their name strings, include:
/// - occupation ("occ"): which atom or molecule occupies a particular crystal
/// site
/// - displacement ("disp"): displacement of crystal sites from an ideal
/// location
/// - strain ("Ustrain", "GLstrain", etc.): how lattice vectors are transformed
/// from the ideal
///   lattice vectors, under various metrics ("Ustrain": stretch tensor,
///   "GLstrain": Green-Lagrange strain metric, etc.)
/// - energy ("energy"): the crystal energy, relative to user chosen reference
/// states
/// - force ("force"): forces on atoms present without/with allowing atomic and
/// lattice relaxation
///
/// The AnisoValTraits class contains all the necessary information for CASM to
/// deal with different types of properties. For each property type there must
/// exist one AnisoValTraits object which controls how CASM treats the property,
/// including:
/// - the property name
/// - whether the property is a local (site) property or global property
/// - whether the property is continuous or discrete
/// - the dimension of the property value vector
/// - the standard basis for the property value (what each component means)
/// - how symmetry operations transform the property value
/// - etc.
///
/// Note: Traits for all properties are defined using AnisoValTraits, even
/// isotropic properties.
///
/// The AnisoValTraits class contains two constructors. One allows defining a
/// new property type which gets placed in a single static container. The other,
/// with just a "name" argument, returns a copy of the already existing
/// AnisoValTraits object for the property type with the given name.
///
///
/// CASM property naming convention
/// -------------------------------
///
/// A particular instance of a property, such as in BasicStructure or
/// SimpleStructure, is labeled with a property name string. Property name
/// strings must end with the property type name (i.e. "occ", "disp", "Ustrain",
/// "energy", etc.) and may also include a modifier describing the particular
/// property (i.e. "formation" in "formation_energy"). If a modifier exists, it
/// must be separated from the property type by an underscore character ('_').
/// The name of a custom property type may not include an underscore.
///
/// Thus, whenever a particular property name is encountered, the property type
/// name can be determined and used to lookup the corresponding AnisoValTraits
/// object.
///
/// Examples:
/// - property name: "energy" -> property type (AnisoValTraits) name: "energy"
/// - property name: "relaxed_energy" -> property type (AnisoValTraits) name:
/// "energy"
/// - property name: "formation_energy" -> property type (AnisoValTraits) name:
/// "energy"
/// - property name: "relaxedenergy" -> Error (because there is no underscore)
/// - property name: "Ustrain" -> property type (AnisoValTraits) name: "Ustrain"
/// - property name: "strain" ->  Error (because there no strain type prefix)
///
class AnisoValTraits {
 public:
  static const unsigned char LOCAL = 0;
  static const unsigned char GLOBAL = (1u << 0);
  static const unsigned char UNIT_LENGTH = (1u << 1);
  static const unsigned char DESCRIBES_ORIENTATION = (1u << 2);
  static const unsigned char EXTENSIVE = (1u << 3);

  ///\brief Named constructor for uninitialized AnisoValTraits
  static AnisoValTraits null();

  ///\brief Named constructor for discrete site occupation AnisoValTraits
  static AnisoValTraits occ();

  ///\brief Named constructor for total energy AnisoValTraits
  static AnisoValTraits energy();

  ///\brief Named constructor for mapping cost
  ///(basis_cost/strain_cost/total_cost) AnisoValTraits
  static AnisoValTraits cost();

  ///\brief Named constructor for selective_dynamics AnisoValTraits
  static AnisoValTraits selective_dynamics();

  ///\brief Named constructor for site displacement AnisoValTraits
  static AnisoValTraits disp();

  ///\brief Named constructor for site coordinate AnisoValTraits
  static AnisoValTraits coordinate();

  ///\brief Named constructor for lattice vector AnisoValTraits
  static AnisoValTraits latvec();

  ///\brief Named constructor for site force AnisoValTraits
  static AnisoValTraits force();

  ///\brief Named constructor for global 'isometry' AnisoValTraits
  static AnisoValTraits isometry();

  ///\brief Named constructor for global strain AnisoValTraits
  /// @param _metric specifies which strain metric. Choices are:
  ///  - "U"  : Right-stretch tensor
  ///  - "B"  : Biot
  ///  - "GL" : Green-Lagrange
  ///  - "AE" : Almansi-Euler
  ///  - "H"  : Hencky
  static AnisoValTraits strain(std::string const &_metric);

  ///\brief Anisovaltraits named constructor for magnetic spin with spin-orbit
  /// coupling
  static AnisoValTraits SOmagspin();

  ///\brief Anisovaltraits named constructor for magnetic spin with spin-orbit
  /// coupling (unit length)
  static AnisoValTraits SOunitmagspin();

  ///\brief Anisovaltraits named constructor for non-collinear magnetic spin
  /// WITHOUT spin-orbit coupling
  static AnisoValTraits NCmagspin();

  ///\brief Anisovaltraits named constructor for non-collinear magnetic spin
  /// WITHOUT spin-orbit coupling (unit length)
  static AnisoValTraits NCunitmagspin();

  ///\brief Anisovaltraits named constructor for collinear magnetic spin
  static AnisoValTraits Cmagspin();

  ///\brief Anisovaltraits named constructor for collinear magnetic spin (unit
  /// length)
  static AnisoValTraits Cunitmagspin();

  ///\brief Named constructor for d-orbital occupation AnisoValTraits
  static AnisoValTraits d_orbital_occupation();

  ///\brief Named constructor for spin-polarized d-orbital occupation
  /// AnisoValTraits
  static AnisoValTraits d_orbital_occupation_spin_polarized();

  /// \brief Given a string, returns string with all characters before the final
  /// @delim character deleted For example, if a trajector has properties with
  /// key1="step1_force", key2="step2_force", etc, then
  /// AnisoValTraits::name_suffix(key1) will return "force" if delim='_'
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

  /// \brief Explicit constructor for AnisoValTraits must specify *all*
  /// attributes at construction There is no mutator for AnisoValTraits.
  /// Explicitly constructing two AnisoValTraits objects with the same name but
  /// non-identical data members will result in a thrown exception.
  AnisoValTraits(
      std::string const &_name, std::vector<std::string> const &_std_var_names,
      unsigned char _options,
      SymRepBuilderInterface const &_symrep_builder = NullSymRepBuilder(),
      std::set<std::string> const &_incompatible = {},
      std::set<std::string> const &_must_apply_before = {},
      std::set<std::string> const &_must_apply_after = {},
      std::vector<std::string> const &_variable_descriptions = {},
      bool _default = false);

  /// \brief Returns previously explicitly initialized AnisoValTraits with name
  /// AnisoValTraits::name_suffix(_name) If no AnisoValTraits with matching name
  /// has been initialized, throws an exception
  AnisoValTraits(std::string const &_name);

  static std::string class_desc() { return "AnisoValTraits"; }

  /// \brief Generate a symmetry representation for the supporting vector space
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

  /// \brief const access of name
  std::string const &name() const { return m_name; }

  /// \brief return true if *this has 'default' designation, meaning it can be
  /// overridden
  bool is_default() const { return m_default; }

  /// \brief return 'options' bitflag
  unsigned char options() const { return m_opt; }

  bool extensive() const { return m_opt & EXTENSIVE; }

  // Please use !AnisoValTraits::global() instead of implementing
  // AnisoValTraits::local()
  // bool local()const {
  //  return !global();
  //}

  /// \brief returns true if DoF is global
  bool global() const { return m_opt & GLOBAL; }

  /// \brief returns true if DoF must always have unit length
  bool unit_length() const { return m_opt & UNIT_LENGTH; }

  /// \brief returns true if time-reversal changes the DoF value
  bool time_reversal_active() const {
    return m_symrep_builder->time_reversal_active();
  }

  /// \brief returns true if value tracks the orientation of an occupying
  /// molecule (not typical)
  bool describes_occupant_orientation() const {
    return m_opt & DESCRIBES_ORIENTATION;
  }

  /// \brief conventional dimensionality of this DoF, returns -1 if always
  /// variable
  Index dim() const { return standard_var_names().size(); }

  /// \brief equality comparison of name
  bool operator==(std::string const &other_name) const {
    return name() == other_name;
  }

  /// \brief lexicographic comparison of name
  bool operator<(std::string const &other_name) const {
    return name() < other_name;
  }

  /// \brief allow implicit conversion to std::string (name)
  operator std::string const &() const { return name(); }

  /// \brief return standard coordinate axes for continuous variable space
  std::vector<std::string> const &standard_var_names() const {
    return m_standard_var_names;
  }

  /// \brief returns expanded description of each standard_var_name
  std::vector<std::string> const &variable_descriptions() const {
    return m_variable_descriptions;
  }

  std::set<std::string> const &incompatible() const { return m_incompatible; }

  /// \brief Return list of DoFs that *must* be applied before this DoF is
  /// applied
  std::set<std::string> const &must_apply_before() const {
    return m_apply_before;
  }

  /// \brief Return list of DoFs that *must* be applied after this DoF is
  /// applied
  std::set<std::string> const &must_apply_after() const {
    return m_apply_after;
  }

  std::unique_ptr<AnisoValTraits> clone() const {
    return notstd::make_unique<AnisoValTraits>(*this);
  }

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

/// \brief comparison of name, domain (discrete/continuous) and mode
/// (local/global)
inline bool identical(AnisoValTraits const &A, AnisoValTraits const &B) {
  return (A.name() == B.name() && A.is_default() == B.is_default() &&
          A.symrep_builder_name() == B.symrep_builder_name() &&
          A.standard_var_names() == B.standard_var_names() &&
          A.options() == B.options() &&
          A.must_apply_before() == B.must_apply_before() &&
          A.must_apply_after() == B.must_apply_after() &&
          A.incompatible() == B.incompatible());
}

namespace AnisoVal_impl {
/// \brief A class to manage dynamic evaluation of BasisFunctions
// struct TraitsConverter {
/*
inline
notstd::cloneable_ptr<DoFType::BasicTraits> traits2cloneable_ptr(const
DoFType::BasicTraits &value) { return
notstd::cloneable_ptr<DoFType::BasicTraits>(value.clone().release());
}
*/
//};
}  // namespace AnisoVal_impl

}  // namespace CASM

#endif
