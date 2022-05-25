#ifndef CASM_enumerator_io_formatter_ConfigEnumIO_impl
#define CASM_enumerator_io_formatter_ConfigEnumIO_impl

#include "casm/app/enum/dataformatter/ConfigEnumIO.hh"
#include "casm/clex/Configuration.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/io/json/Configuration_json_io.hh"
#include "casm/database/ConfigDatabaseTools_impl.hh"
#include "casm/enumerator/ConfigEnumInput_impl.hh"
#include "casm/misc/TypeInfo.hh"

namespace CASM {

// Note:
// - The following formatters are all templated to expect a ConfigEnumDataType
// that follows the
//   minimum interface defined by `ConfigEnumData`
// - Besides the functions included here, you can also use
// `make_datum_formatter_adapter` to adapt
//   DatumFormatters from ConfigIO or ScelIO for use with ConfigEnumDataType.
//   Examples:
//   - `make_datum_formatter_adapter<ConfigEnumDataType,
//   Configuration>(ConfigIO::scelname())`
//   - `make_datum_formatter_adapter<ConfigEnumDataType,
//   Supercell>(ScelIO::TransfMat())`
//   - Note that `make_datum_formatter_adapter` has overrides that can change
//   the formatter name
//     if that is appropriate.
//   - Using `make_datum_formatter_adapter` requires defining an
//   adapter::Adapter

namespace adapter {

template <typename EnumeratorType, typename InitialStateType>
Configuration const &
Adapter<Configuration, ConfigEnumData<EnumeratorType, InitialStateType>>::
operator()(
    ConfigEnumData<EnumeratorType, InitialStateType> const &adaptable) const {
  return adaptable.configuration;
}

template <typename EnumeratorType, typename InitialStateType>
Supercell const &
Adapter<Supercell, ConfigEnumData<EnumeratorType, InitialStateType>>::
operator()(
    ConfigEnumData<EnumeratorType, InitialStateType> const &adaptable) const {
  return adaptable.configuration.supercell();
}

}  // namespace adapter

namespace ConfigEnumIO {

template <typename ConfigEnumDataType>
GenericDatumFormatter<std::string, ConfigEnumDataType> name() {
  return GenericDatumFormatter<std::string, ConfigEnumDataType>(
      "name", "Configuration name, if inserted, else \"none\"",
      [](ConfigEnumDataType const &data) -> std::string {
        if (data.is_excluded_by_filter) {
          return "none";
        }
        if (data.insert_result.canonical_it ==
            data.primclex.template db<Configuration>().end()) {
          return "none";
        }
        return data.insert_result.canonical_it->name();
      });
}

template <typename ConfigEnumDataType>
GenericDatumFormatter<bool, ConfigEnumDataType> selected() {
  return GenericDatumFormatter<bool, ConfigEnumDataType>(
      "selected", "True if configuration exists in database, else false",
      [](ConfigEnumDataType const &data) {
        if (data.is_excluded_by_filter) {
          return false;
        }
        return data.insert_result.canonical_it !=
               data.primclex.template db<Configuration>().end();
      });
}

template <typename ConfigEnumDataType>
GenericDatumFormatter<bool, ConfigEnumDataType> is_new() {
  return GenericDatumFormatter<bool, ConfigEnumDataType>(
      "is_new",
      "True if canonical configuration did not exist before insertion. False "
      "if not checked or was existing before insertion.",
      [](ConfigEnumDataType const &data) {
        if (data.is_excluded_by_filter) {
          return false;
        }
        return data.insert_result.insert_canonical;
      });
}

template <typename ConfigEnumDataType>
GenericDatumFormatter<bool, ConfigEnumDataType> is_existing() {
  return GenericDatumFormatter<bool, ConfigEnumDataType>(
      "is_existing",
      "True if canonical configuration did exist before insertion. False if "
      "not checked or not existing before insertion.",
      [](ConfigEnumDataType const &data) {
        if (data.is_excluded_by_filter) {
          return false;
        }
        return !data.insert_result.insert_canonical;
      });
}

template <typename ConfigEnumDataType>
GenericDatumFormatter<bool, ConfigEnumDataType> is_excluded_by_filter() {
  return GenericDatumFormatter<bool, ConfigEnumDataType>(
      "is_excluded_by_filter",
      "True if configuration was excluded by the filter, else false.",
      [](ConfigEnumDataType const &data) {
        return data.is_excluded_by_filter;
      });
}

template <typename ConfigEnumDataType>
GenericDatumFormatter<Index, ConfigEnumDataType> initial_state_index() {
  return GenericDatumFormatter<Index, ConfigEnumDataType>(
      "initial_state_index", "Initial state index. From 0.",
      [](ConfigEnumDataType const &data) { return data.initial_state_index; });
}

template <typename ConfigEnumDataType>
GenericDatumFormatter<std::string, ConfigEnumDataType> initial_state_name() {
  return GenericDatumFormatter<std::string, ConfigEnumDataType>(
      "initial_state_name",
      "Initial state name. Typically the initital state configuration name "
      "plus description of selected sites.",
      [](ConfigEnumDataType const &data) { return data.initial_state_name; });
}

template <typename ConfigEnumDataType>
GenericDatumFormatter<std::string, ConfigEnumDataType>
initial_state_configname() {
  return GenericDatumFormatter<std::string, ConfigEnumDataType>(
      "initial_state_configname", "Initial state configuration name.",
      [](ConfigEnumDataType const &data) {
        return data.initial_state.configuration().name();
      });
}

template <typename ConfigEnumDataType>
GenericDatumFormatter<jsonParser, ConfigEnumDataType> initial_state() {
  return GenericDatumFormatter<jsonParser, ConfigEnumDataType>(
      "initial_state",
      "Output initial enumeration state, as configuration and selected sites.",
      [](ConfigEnumDataType const &data) -> jsonParser {
        jsonParser json = jsonParser::object();
        to_json(data.initial_state, json);
        return json;
      });
}

template <typename ConfigEnumDataType>
Generic1DDatumFormatter<std::vector<Index>, ConfigEnumDataType>
selected_sites() {
  return Generic1DDatumFormatter<std::vector<Index>, ConfigEnumDataType>(
      "selected_sites",
      "Indices of selected sites in the initial state configuration. JSON "
      "output only.",
      [](ConfigEnumDataType const &data) {
        auto const &sites = data.initial_state->sites();
        return std::vector<Index>{sites.begin(), sites.end()};
      });
}

template <typename ConfigEnumDataType>
GenericDatumFormatter<Index, ConfigEnumDataType> n_selected_sites() {
  return GenericDatumFormatter<Index, ConfigEnumDataType>(
      "n_selected_sites",
      "Number of selected sites in the initial state configuration.",
      [](ConfigEnumDataType const &data) {
        return data.initial_state.sites().size();
      });
}

/// Template that may be specialized by enumerator type to get the current
/// normal coordinate
template <typename EnumeratorType>
Eigen::VectorXd get_normal_coordinate(EnumeratorType const &enumerator) {
  std::stringstream msg;
  msg << "Error in get_normal_coordinate: Not supported for this enumerator '"
      << type_name<EnumeratorType>() << "'";
  throw std::runtime_error(msg.str());
}

template <typename ConfigEnumDataType>
Generic1DDatumFormatter<Eigen::VectorXd, ConfigEnumDataType>
normal_coordinate() {
  return Generic1DDatumFormatter<Eigen::VectorXd, ConfigEnumDataType>(
      "normal_coordinate",
      "Normal coordinate used to generate DoF value. Supported enumerators "
      "only.",
      [=](ConfigEnumDataType const &data) -> Eigen::VectorXd {
        return get_normal_coordinate(data.enumerator);
      });
}

template <typename ConfigEnumDataType>
GenericDatumFormatter<jsonParser, ConfigEnumDataType> config() {
  return GenericDatumFormatter<jsonParser, ConfigEnumDataType>(
      "config",
      "Output as JSON the as-enumerated configuration, before being made "
      "primitive and/or canonical for database insertion.",
      [=](ConfigEnumDataType const &data) -> jsonParser {
        jsonParser json = jsonParser::object();
        Configuration const &config = data.configuration;
        to_json(config, json);
        return json;
      });
}

}  // namespace ConfigEnumIO

}  // namespace CASM

#endif
