#include "casm/container/Enumerator_impl.hh"
#include "casm/crystallography/SupercellEnumerator.hh"
#include "casm/clex/ScelEnum.hh"

namespace CASM {

  /// \brief Standardizes parsing casm enum input options to make ScelEnum JSON input
  jsonParser make_enumerator_scel_enum_input(
    const jsonParser &_kwargs,
    const Completer::EnumOption &enum_opt) {

    jsonParser kwargs {_kwargs};
    if(kwargs.is_null()) {
      kwargs = jsonParser::object();
    }

    jsonParser scel_input;
    scel_input["existing_only"] = true;
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

    if(enum_opt.vm().count("scelnames")) {
      scel_input["name"] = enum_opt.supercell_strs();
    }

    return scel_input;
  }

  /// \brief Standardizes parsing casm enum input options to make an ScelEnumProps
  ScelEnumProps make_enumerator_scel_enum_props(
    PrimClex &primclex,
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
    PrimClex &primclex,
    const jsonParser &_kwargs,
    const Completer::EnumOption &enum_opt) {

    ScelEnumProps enum_props = make_enumerator_scel_enum_props(
                                 primclex,
                                 _kwargs,
                                 enum_opt);

    Supercell unit_cell(&primclex, enum_props.generating_matrix());

    return notstd::make_unique<SupercellEnumerator<Lattice> >(
             primclex.get_prim().lattice(),
             unit_cell.factor_group(),
             enum_props);
  }

  /// \brief Standardizes parsing casm enum input options to make an ScelEnum
  ///
  /// See ConfigEnumAllOccupations for example documentation
  std::unique_ptr<ScelEnum> make_enumerator_scel_enum(
    PrimClex &primclex,
    const jsonParser &_kwargs,
    const Completer::EnumOption &enum_opt) {

    return notstd::make_unique<ScelEnum>(
             primclex,
             make_enumerator_scel_enum_input(_kwargs, enum_opt));
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
}
