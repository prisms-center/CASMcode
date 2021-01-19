#ifndef CASM_DataFormatterFilter_impl
#define CASM_DataFormatterFilter_impl

#include "casm/casm_io/dataformatter/DataFormatterFilter.hh"
#include "casm/casm_io/dataformatter/DataFormatter_impl.hh"
#include "casm/casm_io/dataformatter/DataStream.hh"

namespace CASM {

template <typename DataObject>
DataFormatterFilter<DataObject>::DataFormatterFilter(
    DataFormatter<DataObject> const &_filter)
    : m_filter(_filter) {}

template <typename DataObject>
bool DataFormatterFilter<DataObject>::operator()(
    DataObject const &object) const {
  ValueDataStream<bool> _stream;
  _stream << m_filter(object);
  return _stream.value();
}

template <typename DataObject>
DataFormatterFilter<DataObject> make_data_formatter_filter(
    std::string const &filter_expr,
    DataFormatterDictionary<DataObject> const &_dict) {
  return DataFormatterFilter<DataObject>{_dict.parse(filter_expr)};
}

template <typename DataObject>
DataFormatterFilter<DataObject> make_data_formatter_filter(
    std::vector<std::string> const &filter_expr,
    DataFormatterDictionary<DataObject> const &_dict) {
  return DataFormatterFilter<DataObject>{_dict.parse(filter_expr)};
}

}  // namespace CASM

#endif
