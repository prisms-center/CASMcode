#include "casm/enumerator/Enumerator_impl.hh"
#include "casm/crystallography/SupercellEnumerator.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/clex/ConfigEnumAllOccupations.hh"
#include "casm/clex/ConfigEnumStrain.hh"
#include "casm/clex/ConfigEnumNormalCoords.hh"
#include "casm/clex/SuperConfigEnum.hh"
#include "casm/clex/ScelEnum.hh"
#include "casm/kinetics/DiffusionTransformationEnum.hh"
#include "casm/kinetics/DiffTransConfigEnumOccPerturbations.hh"
#include "casm/kinetics/DiffTransConfigInterpolation.hh"
#include "casm/kinetics/EnumDiffTransConfigEndpoints.hh"
#include "casm/app/enum.hh"

namespace CASM {

  ConfigEnumInput::ConfigEnumInput(Configuration const &_config, std::set<Index> const &_sites_selection, std::string const &_name) :
    m_name(_name),
    m_sites_selection(_sites_selection),
    m_config(_config) {

    if(m_name.empty())
      m_name = _config.name();

    if(m_sites_selection.empty()) {
      for(Index i = 0; i < _config.size(); ++i)
        m_sites_selection.insert(i);
    }

  }

  ConfigEnumInput::ConfigEnumInput(Supercell const &_scel, std::set<Index> const &_sites_selection) :
    ConfigEnumInput(Configuration::zeros(_scel), _sites_selection, _scel.name()) {}

  void ConfigEnumInput::_generate_group() const {
    for(PermuteIterator const &perm_it : config().factor_group()) {
      bool add_it = true;

      for(Index s : sites()) {
        if(sites().count(perm_it.permute_ind(s)) == 0) {
          add_it = false;
          break;
        }
      }
      if(add_it)
        m_group.push_back(perm_it);
    }

  }

  void ConfigEnumInput::_add_site(Index b) {
    Index V = m_config.supercell().volume();
    for(Index i = b * V; i < (b + 1)*V; ++i) {
      m_sites_selection.insert(i);
    }
  }

  void ConfigEnumInput::_add_site(UnitCellCoord const &_ucc) {
    m_sites_selection.insert(m_config.supercell().prim_grid().find(_ucc));
  }



  /// \brief Use to construct an InterfaceMap
  std::unique_ptr<InterfaceMap<Completer::EnumOption> > make_enumerator_map() {
    return make_interface_map<Completer::EnumOption>();
  }

  /// \brief Use to construct an EnumeratorMap with standard Enumerators (not plugins)
  std::unique_ptr<EnumeratorMap> make_standard_enumerator_map() {
    std::unique_ptr<EnumeratorMap> emap = make_enumerator_map();

    emap->insert(
      EnumInterface<ScelEnum>(),
      EnumInterface<ConfigEnumAllOccupations>(),
      EnumInterface<ConfigEnumStrain>(),
      EnumInterface<ConfigEnumNormalCoords>(),
      EnumInterface<SuperConfigEnum>(),
      EnumInterface<Kinetics::DiffusionTransformationEnum>(),
      EnumInterface<Kinetics::DiffTransConfigEnumOccPerturbations>(),
      EnumInterface<Kinetics::DiffTransConfigInterpolation>(),
      EnumInterface<Kinetics::EnumDiffTransConfigEndpoints>()
    );

    return emap;
  }

  /// \brief Standardizes parsing casm enum input options to make ScelEnum JSON input
  jsonParser make_enumerator_scel_enum_input(
    const jsonParser &_kwargs,
    const Completer::EnumOption &enum_opt) {

    jsonParser kwargs {_kwargs};
    if(kwargs.is_null()) {
      kwargs = jsonParser::object();
    }

    jsonParser scel_input;
    kwargs.get_if(scel_input, "supercells");

    // check supercell shortcuts
    if(enum_opt.vm().count("min")) {
      scel_input["min"] = enum_opt.min_volume();
    }

    if(enum_opt.vm().count("max")) {
      scel_input["max"] = enum_opt.max_volume();
    }

    if(enum_opt.all_existing()) {
      scel_input.erase("min");
      scel_input.erase("max");
      scel_input["existing_only"] = true;
    }



    std::vector<std::string> scelnames;
    kwargs.get_if(scelnames, "scelnames");
    if(enum_opt.vm().count("scelnames")) {
      for(std::string const &scelname : enum_opt.supercell_strs())
        scelnames.push_back(scelname);
    }
    if(scelnames.size())
      scel_input["name"] = scelnames;

    if(scel_input.begin() == scel_input.end())
      scel_input["name"].put_array();
    scel_input["existing_only"] = true;

    return scel_input;
  }

  /// \brief Standardizes parsing casm enum input options to make an ScelEnumProps
  ScelEnumProps make_enumerator_scel_enum_props(
    const PrimClex &primclex,
    const jsonParser &_kwargs,
    const Completer::EnumOption &enum_opt) {

    return make_scel_enum_props(
             primclex,
             make_enumerator_scel_enum_input(_kwargs, enum_opt));
  }

  /// \brief Standardizes parsing casm enum input options to make an SupercellEnumerator<Lattice>
  ///
  /// See SuperConfigEnum for example documentation
  std::unique_ptr<SupercellEnumerator<Lattice> > make_enumerator_superlat_enum(
    const PrimClex &primclex,
    const jsonParser &_kwargs,
    const Completer::EnumOption &enum_opt) {

    ScelEnumProps enum_props = make_enumerator_scel_enum_props(
                                 primclex,
                                 _kwargs,
                                 enum_opt);

    Supercell unit_cell(&primclex, enum_props.generating_matrix());

    return notstd::make_unique<SupercellEnumerator<Lattice> >(
             primclex.prim().lattice(),
             unit_cell.factor_group(),
             enum_props);
  }

  /// \brief Standardizes parsing casm enum input options to make an ScelEnum
  ///
  /// See ConfigEnumAllOccupations for example documentation
  std::unique_ptr<ScelEnum> make_enumerator_scel_enum(
    const PrimClex &primclex,
    const jsonParser &_kwargs,
    const Completer::EnumOption &enum_opt) {

    return notstd::make_unique<ScelEnum>(
             primclex,
             make_enumerator_scel_enum_input(_kwargs, enum_opt));
  }

  std::vector<ConfigEnumInput> make_enumerator_input_configs(
    PrimClex const &primclex,
    jsonParser const &_kwargs,
    Completer::EnumOption const &enum_opt,
    EnumeratorMap const *interface_map) {
    std::vector<ConfigEnumInput> result;
    std::vector<std::string> confignames;
    try {
      auto scel_enum = make_enumerator_scel_enum(primclex, _kwargs, enum_opt);
      for(auto const &scel : *scel_enum)
        result.push_back(ConfigEnumInput(scel));

      _kwargs.get_if(confignames, "confignames");
      if(enum_opt.vm().count("confignames")) {
        for(std::string const &configname : enum_opt.config_strs())
          confignames.push_back(configname);
      }
      for(std::string const &configname : confignames) {
        auto it = primclex.const_db<Configuration>().find(configname);
        auto end = primclex.const_db<Configuration>().end();
        if(it == end)
          throw std::runtime_error("Attempting to parse enumeration input, but no valid config exists named '" + configname + "'.");
        result.push_back(ConfigEnumInput(*it));
      }

      if(_kwargs.contains("configs_from_enum")) {
        throw std::runtime_error("How do we extract configurations from enumerator? OutputIterator interface?");
        //interface_map->run(primclex, _kwargs["configs_from_enum"], Completer::EnumOption(), interface_map);
      }

      bool is_init = false;
      auto find_it = _kwargs.find("sublats");
      if(find_it != _kwargs.end()) {
        is_init = true;
        std::vector<Index> sublats;
        find_it->get(sublats);
        for(ConfigEnumInput &config : result) {
          config.set_sites(sublats);
        }
      }

      find_it = _kwargs.find("sites");
      if(find_it != _kwargs.end()) {
        std::vector<UnitCellCoord> sites;
        find_it->get(sites, primclex.prim());
        Index l = 0;
        if(is_init) {
          l = result.size();
          result.reserve(2 * result.size());

          for(Index i = 0; i < l; ++i)
            result.push_back(result[i]);
        }
        for(; l < result.size(); ++l)
          result[l].set_sites(sites);
      }
    }
    catch(std::exception const &e) {
      std::cout << "ERROR: " << e.what() << "\n";
      throw std::runtime_error(std::string("Unable to parse configurtion input arguments:\n") + e.what());
    }
    return result;
  }

  /// \brief Standardizes parsing casm enum filter expressions
  ///
  /// See ConfigEnumAllOccupations for example documentation
  std::vector<std::string> make_enumerator_filter_expr(
    const jsonParser &_kwargs,
    const Completer::EnumOption &enum_opt) {

    std::vector<std::string> filter_expr;
    // check shortcuts
    if(enum_opt.vm().count("filter")) {
      filter_expr = enum_opt.filter_strs();
    }
    else if(_kwargs.contains("filter")) {
      filter_expr.push_back(_kwargs["filter"].get<std::string>());
    };

    return filter_expr;
  }

  /// \brief Get dry-run value (default=false)
  bool dry_run(
    const jsonParser &_kwargs,
    const Completer::EnumOption &enum_opt) {
    if(enum_opt.vm().count("dry-run")) {
      return true;
    }
    else {
      return _kwargs.get_if_else<bool>("dry_run", false);
    }
  }
}
