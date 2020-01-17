#include "casm/crystallography/AnisoValTraits.hh"
#include <map>
#include <string>
namespace CASM {

  namespace Local {

    static std::map<std::string, AnisoValTraits> &traits_map() {
      static std::map<std::string, AnisoValTraits> map;
      return map;
    }

    static void register_traits(AnisoValTraits new_traits) {
      auto it = traits_map().find(new_traits.name());
      // Potential name collisiont if there is already entry with new_traits.name() AND
      // it does not carry 'default' designation. If it has 'default' designation, it can
      // be overridden
      if(it != traits_map().end() && !(it->second).is_default()) {
        // If new_traits has 'default' designation or if new_traits is identical to the existing
        // traits of the same name, then we do nothing. Otherwise, throw an exception
        if(!(new_traits.is_default() || identical(new_traits, it->second))) {
          throw std::runtime_error("Name collision: Attempting to register new AnisoValTraits '"
                                   + new_traits.name() + "' but an incompatible AnisoValTraits with that type already exists.");
        }
      }
      else {
        traits_map().emplace(std::make_pair(new_traits.name(), std::move(new_traits)));
      }
    }

    static AnisoValTraits const &traits(std::string const &name) {
      auto it = traits_map().find(name);
      if(it == traits_map().end()) {
        throw std::runtime_error("Attempted to construct AnisoValTraits '"
                                 + name + "' but no AnisoValTraits with that name is known.");
      }
      return it->second;
    }
  }

  AnisoValTraits::AnisoValTraits(std::string const &_name) {
    *this = Local::traits(name_suffix(_name));
  }

  AnisoValTraits::AnisoValTraits(std::string const &_name,
                                 std::vector<std::string> const &_std_var_names,
                                 unsigned char _options,
                                 SymRepBuilderInterface const &_symrep_builder /*= NullSymRepBuilder()*/,
                                 std::set<std::string> const &_incompatible /*= {}*/,
                                 std::set<std::string> const &_must_apply_before /*= {}*/,
                                 std::set<std::string> const &_must_apply_after /*= {}*/,
                                 bool _default /*= false*/) :
    m_name(name_suffix(_name)),
    m_default(_default),
    m_standard_var_names(_std_var_names),
    m_opt(_options),
    m_symrep_builder(_symrep_builder.clone()),
    m_incompatible(_incompatible),
    m_apply_before(_must_apply_before),
    m_apply_after(_must_apply_after) {

    if(m_name != _name) {
      throw std::runtime_error("Attempting to initialize AnisoValTraits object that does not satisfy naming constraints. Name '" + _name + "' was reduced to '" + m_name + "'."
                               + " Note: The underscore character '_' is not allowed.\n");
    }
    Local::register_traits(*this);
  }


  AnisoValTraits AnisoValTraits::null() {
    return AnisoValTraits("NULL",
                          {},
                          LOCAL,
                          NullSymRepBuilder(),
                          {},
                          {});
  }

  AnisoValTraits AnisoValTraits::energy() {
    return AnisoValTraits("energy",
    {"E"},
    GLOBAL,
    SymRepBuilder::Identity(),
    {},
    {},
    {},
    true);
  }

  AnisoValTraits AnisoValTraits::cost() {
    return AnisoValTraits("cost",
    {"C"},
    GLOBAL,
    SymRepBuilder::Identity(),
    {},
    {},
    {},
    true);
  }

  AnisoValTraits AnisoValTraits::selective_dynamics() {
    return AnisoValTraits("selectivedynamics",
    {"xflag", "yflag", "zflag"},
    LOCAL,
    SymRepBuilder::Identity(),
    {},
    {},
    {},
    true);
  }

  AnisoValTraits AnisoValTraits::disp() {
    return AnisoValTraits("disp",
    {"dx", "dy", "dz"},
    LOCAL,
    CartesianSymRepBuilder(),
    {},
    {},
    {"atomize"},
    true);
  }

  AnisoValTraits AnisoValTraits::coordinate() {
    return AnisoValTraits("coordinate",
    {"cx", "cy", "cz"},
    LOCAL,
    CartesianSymRepBuilder(),
    {},
    {},
    {"atomize"},
    true);
  }

  AnisoValTraits AnisoValTraits::latvec() {
    return AnisoValTraits("latvec",
    {"L1x", "L1y", "L1z", "L2x", "L2y", "L2z", "L3x", "L3y", "L3z"},
    GLOBAL,
    KroneckerSymRepBuilder<IdentitySymRepBuilder, CartesianSymRepBuilder, 3, 3>("Rank2AsymTensor"),
    {},
    {},
    {"atomize"},
    true);
  }

  AnisoValTraits AnisoValTraits::force() {
    return AnisoValTraits("force",
    {"fx", "fy", "fz"},
    LOCAL | EXTENSIVE,
    CartesianSymRepBuilder(),
    {},
    {"atomize"},
    {},
    true);
  }

  AnisoValTraits AnisoValTraits::isometry() {
    return AnisoValTraits("isometry",
    {"S11", "S21", "S31", "S12", "S22", "S32", "S13", "S23", "S33"},
    GLOBAL,
    KroneckerSymRepBuilder<IdentitySymRepBuilder, CartesianSymRepBuilder, 3, 3>("Rank2AsymTensor"),
    {},
    {"atomize", "disp"},
    {},
    true);
  }

  AnisoValTraits AnisoValTraits::strain(std::string const &_prefix) {
    return AnisoValTraits(_prefix + "strain",
    {"e_1", "e_2", "e_3", "e_4", "e_5", "e_6"},
    GLOBAL,
    Rank2TensorSymRepBuilder(),
    {},
    {"atomize", "disp"},
    {},
    true);
  }

  AnisoValTraits AnisoValTraits::magspin() {
    return AnisoValTraits("magspin",
    {"sx", "sy", "sz"},
    LOCAL | EXTENSIVE,
    AngularMomentumSymRepBuilder(),
    {},
    {"atomize"},
    {},
    true);
  }

  AnisoValTraits AnisoValTraits::magmom() {
    return AnisoValTraits("magmom",
    {"mx", "my", "mz"},
    LOCAL | UNIT_LENGTH | EXTENSIVE,
    TimeReversalSymRepBuilder(),
    {},
    {"atomize"},
    {},
    true);
  }

  AnisoValTraits AnisoValTraits::d_orbital_occupation() {
    return AnisoValTraits("dorbitaloccupation",
    {"v_1", "v_2", "v_3", "v_4", "v_5", "v_6", "v_7", "v_8", "v_9", "v_10", "v_11", "v_12", "v_13", "v_14", "v_15"},
    LOCAL,
    dOrbitalOccupationSymRepBuilder(),
    {},
    {"atomize"});
  }

  AnisoValTraits AnisoValTraits::d_orbital_occupation_spin_polarized() {
    return AnisoValTraits("dorbitaloccupationspinpolarized",
    {"u_1", "u_2", "u_3", "u_4", "u_5", "u_6", "u_7", "u_8", "u_9", "u_10", "u_11", "u_12", "u_13", "u_14", "u_15",
     "d_1", "d_2", "d_3", "d_4", "d_5", "d_6", "d_7", "d_8", "d_9", "d_10", "d_11", "d_12", "d_13", "d_14", "d_15"},
    LOCAL,
    dOrbitalOccupationSpinPolarizedSymRepBuilder(),
    {},
    {"atomize"});
  }

}
