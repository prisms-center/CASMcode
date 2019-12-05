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
    SymRepBuilder::Identity());
  }

  AnisoValTraits AnisoValTraits::selective_dynamics() {
    return AnisoValTraits("selectivedynamics",
    {"xflag", "yflag", "zflag"},
    LOCAL,
    SymRepBuilder::Identity());
  }

  AnisoValTraits AnisoValTraits::disp() {
    return AnisoValTraits("disp",
    {"dx", "dy", "dz"},
    LOCAL,
    CartesianSymRepBuilder(),
    {},
    {},
    {"atomize"});
  }

  AnisoValTraits AnisoValTraits::force() {
    return AnisoValTraits("force",
    {"fx", "fy", "fz"},
    LOCAL,
    CartesianSymRepBuilder(),
    {},
    {"atomize"});
  }

  AnisoValTraits AnisoValTraits::strain(std::string const &_prefix) {
    return AnisoValTraits(_prefix + "strain",
    {"e_1", "e_2", "e_3", "e_4", "e_5", "e_6"},
    GLOBAL,
    Rank2TensorSymRepBuilder(),
    {},
    {"atomize", "disp"});
  }

  AnisoValTraits AnisoValTraits::magspin() {
    return AnisoValTraits("magspin",
    {"sx", "sy", "sz"},
    LOCAL,
    AngularMomentumSymRepBuilder(),
    {},
    {"atomize"});
  }

  AnisoValTraits AnisoValTraits::magmom() {
    return AnisoValTraits("magmom",
    {"mx", "my", "mz"},
    LOCAL | UNIT_LENGTH,
    TimeReversalSymRepBuilder(),
    {},
    {"atomize"});
  }

}
