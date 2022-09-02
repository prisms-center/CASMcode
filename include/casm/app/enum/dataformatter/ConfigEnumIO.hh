#ifndef CASM_enumerator_io_formatter_ConfigEnumIO
#define CASM_enumerator_io_formatter_ConfigEnumIO

#include "casm/misc/TypeInfo.hh"

namespace CASM {

template <typename EnumeratorType, typename InitialStateType>
struct ConfigEnumData;
class Supercell;

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
//   the name if that
//     is appropriate.
//   - Using `make_datum_formatter_adapter` requires defining an
//   adapter::Adapter

namespace adapter {

template <typename ToType, typename FromType>
struct Adapter;

template <typename EnumeratorType, typename InitialStateType>
struct Adapter<Configuration,
               ConfigEnumData<EnumeratorType, InitialStateType>> {
  Configuration const &operator()(
      ConfigEnumData<EnumeratorType, InitialStateType> const &adaptable) const;
};

template <typename EnumeratorType, typename InitialStateType>
struct Adapter<Supercell, ConfigEnumData<EnumeratorType, InitialStateType>> {
  Supercell const &operator()(
      ConfigEnumData<EnumeratorType, InitialStateType> const &adaptable) const;
};
}  // namespace adapter

namespace ConfigEnumIO {

template <typename ConfigEnumDataType>
GenericDatumFormatter<std::string, ConfigEnumDataType> canonical_configname();

template <typename ConfigEnumDataType>
GenericDatumFormatter<bool, ConfigEnumDataType> selected();

template <typename ConfigEnumDataType>
GenericDatumFormatter<bool, ConfigEnumDataType> is_new();

template <typename ConfigEnumDataType>
GenericDatumFormatter<bool, ConfigEnumDataType> is_existing();

template <typename ConfigEnumDataType>
GenericDatumFormatter<bool, ConfigEnumDataType> is_excluded_by_filter();

template <typename ConfigEnumDataType>
GenericDatumFormatter<Index, ConfigEnumDataType> initial_state_index();

template <typename ConfigEnumDataType>
GenericDatumFormatter<std::string, ConfigEnumDataType> initial_state_name();

template <typename ConfigEnumDataType>
GenericDatumFormatter<std::string, ConfigEnumDataType>
initial_state_configname();

template <typename ConfigEnumDataType>
GenericDatumFormatter<jsonParser, ConfigEnumDataType> initial_state();

template <typename ConfigEnumDataType>
Generic1DDatumFormatter<std::vector<Index>, ConfigEnumDataType>
selected_sites();

template <typename ConfigEnumDataType>
GenericDatumFormatter<Index, ConfigEnumDataType> n_selected_sites();

/// Template that may be specialized by enumerator type to get the current
/// normal coordinate and enable the `normal_coordinate` DatumFormatter. Will
/// throw if no specialization exists.
template <typename EnumeratorType>
Eigen::VectorXd get_normal_coordinate(EnumeratorType const &enumerator);

template <typename ConfigEnumDataType>
Generic1DDatumFormatter<Eigen::VectorXd, ConfigEnumDataType>
normal_coordinate();

template <typename ConfigEnumDataType>
GenericDatumFormatter<jsonParser, ConfigEnumDataType> config();

template <typename ConfigEnumDataType>
GenericDatumFormatter<jsonParser, ConfigEnumDataType> canonical_config();

}  // namespace ConfigEnumIO
}  // namespace CASM

#endif
