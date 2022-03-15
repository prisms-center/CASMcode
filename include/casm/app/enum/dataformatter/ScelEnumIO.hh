#ifndef CASM_enumerator_io_formatter_ConfigEnumIO
#define CASM_enumerator_io_formatter_ConfigEnumIO

#include "casm/misc/TypeInfo.hh"

namespace CASM {

template <typename EnumeratorType>
struct ScelEnumData;
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

template <typename EnumeratorType>
struct Adapter<Supercell, ScelEnumData<EnumeratorType>> {
  Supercell const &operator()(
      ScelEnumData<EnumeratorType> const &adaptable) const;
};

}  // namespace adapter

namespace ScelEnumIO {

template <typename ScelEnumDataType>
GenericDatumFormatter<std::string, ScelEnumDataType> name();

template <typename ScelEnumDataType>
GenericDatumFormatter<bool, ScelEnumDataType> selected();

template <typename ScelEnumDataType>
GenericDatumFormatter<bool, ScelEnumDataType> is_new();

template <typename ScelEnumDataType>
GenericDatumFormatter<bool, ScelEnumDataType> is_existing();

template <typename ScelEnumDataType>
GenericDatumFormatter<bool, ScelEnumDataType> is_excluded_by_filter();

template <typename ScelEnumDataType>
Generic2DDatumFormatter<Eigen::MatrixXi, ScelEnumDataType>
transformation_matrix_to_super();

}  // namespace ScelEnumIO
}  // namespace CASM

#endif
