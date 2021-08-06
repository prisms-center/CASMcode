#include "casm/crystallography/AnisoValTraits.hh"

#include <map>
#include <string>

#include "casm/misc/ParsingDictionary.hh"
namespace CASM {

/// \class AnisoValTraits
/// \brief Specifies traits of (possibly) anisotropic crystal properties
///
/// CASM needs to know how to use crystal properties in several places:
/// - Site and global properties in xtal::SimpleStructure
/// - Degrees of freedoms in xtal::BasicStructure through xtal::DoFSet and
/// xtal::SiteDoFSet
/// - Molecule properties in xtal::SpeciesProperty
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
/// Named constructors for standard types
/// -------------------------------------
///
/// For standard types, AnisoValTraits are pre-defined and can be accessed via
/// the following named constructors:
/// - AnisoValTraits::occ()
/// - AnisoValTraits::disp()
/// - AnisoValTraits::energy()
/// - AnisoValTraits::cost()
/// - AnisoValTraits::coordinate()
/// - AnisoValTraits::latvec()
/// - AnisoValTraits::selective_dynamics()
/// - AnisoValTraits::SOmagspin()
/// - AnisoValTraits::SOunitmagspin()
/// - AnisoValTraits::NCmagspin()
/// - AnisoValTraits::NCunitmagspin()
/// - AnisoValTraits::Cmagspin()
/// - AnisoValTraits::Cunitmagspin()
/// - AnisoValTraits::isometry()
/// - AnisoValTraits::strain("B")
/// - AnisoValTraits::strain("U")
/// - AnisoValTraits::strain("EA")
/// - AnisoValTraits::strain("GL")
/// - AnisoValTraits::strain("H")
/// - AnisoValTraits::force()
///

template <>
ParsingDictionary<AnisoValTraits> make_parsing_dictionary<AnisoValTraits>() {
  ParsingDictionary<AnisoValTraits> dict;

  dict.insert(AnisoValTraits::disp(), AnisoValTraits::energy(),
              AnisoValTraits::cost(), AnisoValTraits::coordinate(),
              AnisoValTraits::latvec(), AnisoValTraits::selective_dynamics(),
              AnisoValTraits::Cmagspin(), AnisoValTraits::Cunitmagspin(),
              AnisoValTraits::NCmagspin(), AnisoValTraits::NCunitmagspin(),
              AnisoValTraits::SOmagspin(), AnisoValTraits::SOunitmagspin(),
              AnisoValTraits::isometry(), AnisoValTraits::strain("B"),
              AnisoValTraits::strain("U"), AnisoValTraits::strain("EA"),
              AnisoValTraits::strain("GL"), AnisoValTraits::strain("H"),
              AnisoValTraits::force(), AnisoValTraits::d_orbital_occupation(),
              AnisoValTraits::d_orbital_occupation_spin_polarized());

  return dict;
}

namespace Local {

static std::map<std::string, AnisoValTraits> &traits_map() {
  static std::map<std::string, AnisoValTraits> map;
  return map;
}

static void register_traits(AnisoValTraits new_traits) {
  auto it = traits_map().find(new_traits.name());
  // Potential name collisiont if there is already entry with new_traits.name()
  // AND it does not carry 'default' designation. If it has 'default'
  // designation, it can be overridden
  if (it != traits_map().end() && !(it->second).is_default()) {
    // If new_traits has 'default' designation or if new_traits is identical to
    // the existing traits of the same name, then we do nothing. Otherwise,
    // throw an exception
    if (!(new_traits.is_default() || identical(new_traits, it->second))) {
      throw std::runtime_error(
          "Name collision: Attempting to register new AnisoValTraits '" +
          new_traits.name() +
          "' but an incompatible AnisoValTraits with that type already "
          "exists.");
    }
  } else {
    traits_map().emplace(
        std::make_pair(new_traits.name(), std::move(new_traits)));
  }
}

static int initialize() {
  AnisoValTraits::occ();
  AnisoValTraits::disp();
  AnisoValTraits::energy();
  AnisoValTraits::cost();
  AnisoValTraits::coordinate();
  AnisoValTraits::latvec();
  AnisoValTraits::selective_dynamics();
  AnisoValTraits::SOmagspin();
  AnisoValTraits::SOunitmagspin();
  AnisoValTraits::NCmagspin();
  AnisoValTraits::NCunitmagspin();
  AnisoValTraits::Cmagspin();
  AnisoValTraits::Cunitmagspin();
  AnisoValTraits::isometry();
  AnisoValTraits::strain("B");
  AnisoValTraits::strain("U");
  AnisoValTraits::strain("EA");
  AnisoValTraits::strain("GL");
  AnisoValTraits::strain("H");
  AnisoValTraits::force();
  return 1;
}

static AnisoValTraits const &traits(std::string const &name) {
  static int i = initialize();
  auto it = traits_map().find(name);
  if (it == traits_map().end()) {
    throw std::runtime_error(
        "Attempted to construct AnisoValTraits '" + name +
        "' but no AnisoValTraits with that name is known.");
  }
  return it->second;
}
}  // namespace Local

/// Construct a copy of an existing AnisoValTraits with matching name suffix
///
/// Returns a copy of a previously explicitly initialized AnisoValTraits with
/// name `AnisoValTraits::name_suffix(_name)`. If no AnisoValTraits with
/// matching name has been initialized, throws an exception.
AnisoValTraits::AnisoValTraits(std::string const &_name) {
  *this = Local::traits(name_suffix(_name));
}

/// Explicit constructor for AnisoValTraits
///
/// *All* attributes must be specified at construction. There is no mutator for
/// AnisoValTraits. Explicitly constructing two AnisoValTraits objects with the
/// same name but non-identical data members will result in a thrown exception,
/// unless `this->is_default()==true`.
///
/// @param _name Name for this type
/// @param _std_var_names Names for standard coordinate axes of continuous
///     variable space.
/// @param _options Bitflag describing setting properties (LOCAL, GLOBAL,
///     UNIT_LENGTH, DESCRIBES_ORIENTATION, EXTENSIVE)
/// @param _symrep_builder Method for constructing symmetry representations
/// @param _incompatible When used as a DoF, this is a set of DoFs that are
///     incompatible with this type of DoF
/// @param _must_apply_before When used as a DoF, this is a set of DoFs (by
///     AnisoValTraits name) that must be applied before this DoF is applied
///     when constructing a SimpleStructure from a Configuration.
/// @param _must_apply_after When used as a DoF, this is a set of DoFs (by
///     AnisoValTraits name) that must be applied after this DoF is applied
///     when constructing a SimpleStructure from a Configuration.
/// @param _variable_descriptions Expanded description of each standard_var_name
/// @param _default True if this has 'default' designation, meaning it can be
///     overridden
///
/// The special keyword "atomize" is also valid for `_must_apply_before` and
/// `_must_apply_after`, indicating the step when molecules are taken and
/// "atomized" (populating SimpleStructure::atom_info from
/// SimpleStructure::mol_info).
///
AnisoValTraits::AnisoValTraits(
    std::string const &_name, std::vector<std::string> const &_std_var_names,
    unsigned char _options,
    SymRepBuilderInterface const &_symrep_builder /*= NullSymRepBuilder()*/,
    std::set<std::string> const &_incompatible /*= {}*/,
    std::set<std::string> const &_must_apply_before /*= {}*/,
    std::set<std::string> const &_must_apply_after /*= {}*/,
    std::vector<std::string> const &_variable_descriptions /*= {}*/,
    bool _default /*= false*/)
    : m_name(name_suffix(_name)),
      m_default(_default),
      m_standard_var_names(_std_var_names),
      m_variable_descriptions(_variable_descriptions),
      m_opt(_options),
      m_symrep_builder(_symrep_builder.clone()),
      m_incompatible(_incompatible),
      m_apply_before(_must_apply_before),
      m_apply_after(_must_apply_after) {
  if (m_name != _name) {
    throw std::runtime_error(
        "Attempting to initialize AnisoValTraits object that does not satisfy "
        "naming constraints. Name '" +
        _name + "' was reduced to '" + m_name + "'." +
        " Note: The underscore character '_' is not allowed.\n");
  }
  if (m_variable_descriptions.empty())
    m_variable_descriptions = m_standard_var_names;

  Local::register_traits(*this);
}

/// \brief Named constructor for uninitialized AnisoValTraits
///
/// Required and non-default values:
/// - name(): "NULL"
/// - standard_var_names(): {}
/// - options(): LOCAL
AnisoValTraits AnisoValTraits::null() {
  return AnisoValTraits("NULL", {}, LOCAL, NullSymRepBuilder(), {}, {}, {});
}

/// \brief Named constructor for discrete site occupation AnisoValTraits
///
/// Required and non-default values:
/// - name(): "occ"
/// - standard_var_names(): {}
/// - options(): LOCAL
AnisoValTraits AnisoValTraits::occ() {
  return AnisoValTraits("occ", {}, LOCAL);
}

/// \brief Named constructor for total energy AnisoValTraits
///
/// Required and non-default values:
/// - name(): "energy"
/// - standard_var_names(): {"E"}
/// - options(): GLOBAL
/// - SymRepBuilderInterface: SymRepBuilder::Identity()
/// - is_default(): true
AnisoValTraits AnisoValTraits::energy() {
  return AnisoValTraits("energy", {"E"}, GLOBAL, SymRepBuilder::Identity(), {},
                        {}, {}, {}, true);
}

/// \brief Named constructor for mapping cost AnisoValTraits
///
/// This includes basis_cost, strain_cost, total_cost.
///
/// Required and non-default values:
/// - name(): "cost"
/// - standard_var_names(): {"C"}
/// - options(): GLOBAL
/// - SymRepBuilderInterface: SymRepBuilder::Identity()
/// - is_default(): true
AnisoValTraits AnisoValTraits::cost() {
  return AnisoValTraits("cost", {"C"}, GLOBAL, SymRepBuilder::Identity(), {},
                        {}, {}, {}, true);
}

/// \brief Named constructor for selective dynamics AnisoValTraits
///
/// Required and non-default values:
/// - name(): "selectivedynamics"
/// - standard_var_names(): {"aflag", "bflag", "cflag"}
/// - options(): LOCAL
/// - SymRepBuilderInterface: SymRepBuilder::Identity()
/// - is_default(): true
AnisoValTraits AnisoValTraits::selective_dynamics() {
  return AnisoValTraits("selectivedynamics", {"aflag", "bflag", "cflag"}, LOCAL,
                        SymRepBuilder::Identity(), {}, {}, {}, {}, true);
}

/// \brief Named constructor for site displacement AnisoValTraits
///
/// Required and non-default values:
/// - name(): "disp"
/// - standard_var_names(): {"dx", "dy", "dz"}
/// - options(): LOCAL
/// - SymRepBuilderInterface: CartesianSymRepBuilder()
/// - must_apply_after(): {"atomize"}
/// - is_default(): true
AnisoValTraits AnisoValTraits::disp() {
  return AnisoValTraits("disp", {"dx", "dy", "dz"}, LOCAL,
                        CartesianSymRepBuilder(), {}, {}, {"atomize"}, {},
                        true);
}

/// \brief Named constructor for site coordinate AnisoValTraits
///
/// Required and non-default values:
/// - name(): "coordinate"
/// - standard_var_names(): {"cx", "cy", "cz"}
/// - options(): LOCAL
/// - SymRepBuilderInterface: CartesianSymRepBuilder()
/// - must_apply_after(): {"atomize"}
/// - is_default(): true
AnisoValTraits AnisoValTraits::coordinate() {
  return AnisoValTraits("coordinate", {"cx", "cy", "cz"}, LOCAL,
                        CartesianSymRepBuilder(), {}, {}, {"atomize"}, {},
                        true);
}

/// \brief Named constructor for lattice vector AnisoValTraits
///
/// Required and non-default values:
/// - name(): "latvec"
/// - standard_var_names(): {"L1x", "L1y", "L1z", "L2x", "L2y", "L2z", "L3x",
///   "L3y", "L3z"}
/// - options(): GLOBAL
/// - SymRepBuilderInterface:
///       KroneckerSymRepBuilder<IdentitySymRepBuilder,
///                              CartesianSymRepBuilder,
///                              3, 3>("Rank2AsymTensor")
/// - must_apply_after(): {"atomize"}
/// - is_default(): true
AnisoValTraits AnisoValTraits::latvec() {
  return AnisoValTraits(
      "latvec", {"L1x", "L1y", "L1z", "L2x", "L2y", "L2z", "L3x", "L3y", "L3z"},
      GLOBAL,
      KroneckerSymRepBuilder<IdentitySymRepBuilder, CartesianSymRepBuilder, 3,
                             3>("Rank2AsymTensor"),
      {}, {}, {"atomize"}, {}, true);
}

/// \brief Named constructor for site force AnisoValTraits
///
/// Required and non-default values:
/// - name(): "force"
/// - standard_var_names(): {"fx", "fy", "fz"}
/// - options(): LOCAL | EXTENSIVE
/// - SymRepBuilderInterface: CartesianSymRepBuilder()
/// - must_apply_before(): {"atomize"}
/// - is_default(): true
AnisoValTraits AnisoValTraits::force() {
  return AnisoValTraits("force", {"fx", "fy", "fz"}, LOCAL | EXTENSIVE,
                        CartesianSymRepBuilder(), {}, {"atomize"}, {}, {},
                        true);
}

/// \brief Named constructor for global 'isometry' AnisoValTraits
///
/// Required and non-default values:
/// - name(): "isometry"
/// - standard_var_names(): {"S11", "S21", "S31", "S12", "S22", "S32", "S13",
///   "S23", "S33"}
/// - options(): GLOBAL
/// - SymRepBuilderInterface:
///       KroneckerSymRepBuilder<IdentitySymRepBuilder,
///                              CartesianSymRepBuilder,
///                              3, 3>("Rank2AsymTensor")
/// - must_apply_before(): {"atomize", "disp"}
/// - is_default(): true
AnisoValTraits AnisoValTraits::isometry() {
  return AnisoValTraits(
      "isometry",
      {"S11", "S21", "S31", "S12", "S22", "S32", "S13", "S23", "S33"}, GLOBAL,
      KroneckerSymRepBuilder<IdentitySymRepBuilder, CartesianSymRepBuilder, 3,
                             3>("Rank2AsymTensor"),
      {}, {"atomize", "disp"}, {}, {}, true);
}

/// \brief Named constructor for global strain AnisoValTraits
///
/// @param _metric specifies which strain metric. Choices are:
///  - "U"  : Right-stretch tensor
///  - "B"  : Biot
///  - "GL" : Green-Lagrange
///  - "AE" : Almansi-Euler
///  - "H"  : Hencky
///
/// Required and non-default values:
/// - name(): _prefix + "strain"
/// - standard_var_names(): {"e_1", "e_2", "e_3", "e_4", "e_5", "e_6"}
/// - options(): GLOBAL
/// - SymRepBuilderInterface: Rank2TensorSymRepBuilder()
/// - must_apply_before(): {"atomize", "disp"}
/// - variable_descriptions(): {"Exx", "Eyy", "Ezz", "sqrt(2)Eyz",
///   "sqrt(2)Exz", "sqrt(2)Exy"}
/// - is_default(): true
AnisoValTraits AnisoValTraits::strain(std::string const &_prefix) {
  return AnisoValTraits(
      _prefix + "strain", {"e_1", "e_2", "e_3", "e_4", "e_5", "e_6"}, GLOBAL,
      Rank2TensorSymRepBuilder(), {}, {"atomize", "disp"}, {},
      {"Exx", "Eyy", "Ezz", "sqrt(2)Eyz", "sqrt(2)Exz", "sqrt(2)Exy"}, true);
}

/// \brief Non-collinear magnetic spin, with spin-orbit coupling
///
/// Required and non-default values:
/// - name(): "SOmagspin"
/// - standard_var_names(): {"sx", "sy", "sz"}
/// - options(): LOCAL | EXTENSIVE
/// - SymRepBuilderInterface: AngularMomentumSymRepBuilder()
/// - incompatible(): {}
/// - must_apply_before(): {}
/// - must_apply_after(): {}
/// - variable_descriptions(): {}
/// - is_default(): false
AnisoValTraits AnisoValTraits::SOmagspin() {
  return AnisoValTraits("SOmagspin", {"sx", "sy", "sz"}, LOCAL | EXTENSIVE,
                        AngularMomentumSymRepBuilder(), {}, {"atomize"}, {}, {},
                        true);
}

/// \brief Non-collinear magnetic spin, with spin-orbit coupling, constrained to
/// unit length
///
/// Required and non-default values:
/// - name(): "SOunitmagspin"
/// - standard_var_names(): {"sx", "sy", "sz"}
/// - options(): LOCAL | UNIT_LENGTH | EXTENSIVE
/// - SymRepBuilderInterface: AngularMomentumSymRepBuilder()
/// - must_apply_before(): {"atomize"}
/// - is_default(): true
AnisoValTraits AnisoValTraits::SOunitmagspin() {
  return AnisoValTraits(
      "SOunitmagspin", {"sx", "sy", "sz"}, LOCAL | UNIT_LENGTH | EXTENSIVE,
      AngularMomentumSymRepBuilder(), {}, {"atomize"}, {}, {}, true);
}

/// \brief Non-collinear magnetic spin, without spin-orbit coupling
///
/// Required and non-default values:
/// - name(): "NCmagspin"
/// - standard_var_names(): {"sx", "sy", "sz"}
/// - options(): LOCAL | EXTENSIVE
/// - SymRepBuilderInterface: TimeReversalSymRepBuilder()
/// - must_apply_before(): {"atomize"}
/// - is_default(): true
AnisoValTraits AnisoValTraits::NCmagspin() {
  return AnisoValTraits("NCmagspin", {"sx", "sy", "sz"}, LOCAL | EXTENSIVE,
                        TimeReversalSymRepBuilder(), {}, {"atomize"}, {}, {},
                        true);
}

/// \brief
/// Non-collinear magnetic spin, without spin-orbit coupling, constrained to
/// unit length
///
/// Required and non-default values:
/// - name(): "NCunitmagspin"
/// - standard_var_names(): {"sx", "sy", "sz"}
/// - options(): LOCAL | UNIT_LENGTH | EXTENSIVE
/// - SymRepBuilderInterface: TimeReversalSymRepBuilder()
/// - must_apply_before(): {"atomize"}
/// - is_default(): true
AnisoValTraits AnisoValTraits::NCunitmagspin() {
  return AnisoValTraits(
      "NCunitmagspin", {"sx", "sy", "sz"}, LOCAL | UNIT_LENGTH | EXTENSIVE,
      TimeReversalSymRepBuilder(), {}, {"atomize"}, {}, {}, true);
}

/// \brief Collinear magnetic spin
///
/// Required and non-default values:
/// - name(): "Cmagspin"
/// - standard_var_names(): {"m"}
/// - options(): LOCAL | EXTENSIVE
/// - SymRepBuilderInterface: TimeReversalSymRepBuilder()
/// - must_apply_before(): {"atomize"}
/// - is_default(): true
AnisoValTraits AnisoValTraits::Cmagspin() {
  return AnisoValTraits("Cmagspin", {"m"}, LOCAL | EXTENSIVE,
                        TimeReversalSymRepBuilder(), {}, {"atomize"}, {}, {},
                        true);
}

/// \brief Collinear magnetic spin, constrained to unit length
///
/// Required and non-default values:
/// - name(): "Cunitmagspin"
/// - standard_var_names(): {"m"}
/// - options(): LOCAL | UNIT_LENGTH | EXTENSIVE
/// - SymRepBuilderInterface: TimeReversalSymRepBuilder()
/// - must_apply_before(): {"atomize"}
/// - is_default(): true
AnisoValTraits AnisoValTraits::Cunitmagspin() {
  return AnisoValTraits("Cunitmagspin", {"m"}, LOCAL | UNIT_LENGTH | EXTENSIVE,
                        TimeReversalSymRepBuilder(), {}, {"atomize"}, {}, {},
                        true);
}

/// \brief Named constructor for d-orbital occupation AnisoValTraits
///
/// Required and non-default values:
/// - name(): "dorbitaloccupation"
/// - standard_var_names(): {"v_1", "v_2", "v_3", "v_4", "v_5", "v_6", "v_7",
///   "v_8", "v_9", "v_10", "v_11", "v_12", "v_13", "v_14", "v_15"}
/// - options(): LOCAL
/// - SymRepBuilderInterface: dOrbitalOccupationSymRepBuilder()
/// - must_apply_before(): {"atomize"}
AnisoValTraits AnisoValTraits::d_orbital_occupation() {
  return AnisoValTraits("dorbitaloccupation",
                        {"v_1", "v_2", "v_3", "v_4", "v_5", "v_6", "v_7", "v_8",
                         "v_9", "v_10", "v_11", "v_12", "v_13", "v_14", "v_15"},
                        LOCAL, dOrbitalOccupationSymRepBuilder(), {},
                        {"atomize"});
}

/// \brief Named constructor for spin-polarized d-orbital occupation
/// AnisoValTraits
///
/// Required and non-default values:
/// - name(): "dorbitaloccupationspinpolarized"
/// - standard_var_names(): {"u_1",  "u_2",  "u_3",  "u_4",  "u_5",  "u_6",
///   "u_7",  "u_8", "u_9",  "u_10", "u_11", "u_12", "u_13", "u_14", "u_15",
///   "d_1", "d_2",  "d_3",  "d_4",  "d_5",  "d_6",  "d_7",  "d_8",  "d_9",
///   "d_10", "d_11", "d_12", "d_13", "d_14", "d_15"}
/// - options(): LOCAL
/// - SymRepBuilderInterface: dOrbitalOccupationSpinPolarizedSymRepBuilder()
/// - must_apply_before(): {"atomize"}
AnisoValTraits AnisoValTraits::d_orbital_occupation_spin_polarized() {
  return AnisoValTraits(
      "dorbitaloccupationspinpolarized",
      {"u_1",  "u_2",  "u_3",  "u_4",  "u_5",  "u_6",  "u_7",  "u_8",
       "u_9",  "u_10", "u_11", "u_12", "u_13", "u_14", "u_15", "d_1",
       "d_2",  "d_3",  "d_4",  "d_5",  "d_6",  "d_7",  "d_8",  "d_9",
       "d_10", "d_11", "d_12", "d_13", "d_14", "d_15"},
      LOCAL, dOrbitalOccupationSpinPolarizedSymRepBuilder(), {}, {"atomize"});
}

}  // namespace CASM
