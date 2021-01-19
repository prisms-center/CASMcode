#ifndef CASM_DataFormatterFilter
#define CASM_DataFormatterFilter

#include "casm/casm_io/dataformatter/DataFormatter.hh"

namespace CASM {

template <typename DataObject>
struct DataFormatterFilter {
  DataFormatterFilter(DataFormatter<DataObject> const &_filter);

  bool operator()(DataObject const &object) const;

  DataFormatter<DataObject> m_filter;
};

template <typename DataObject>
DataFormatterFilter<DataObject> make_data_formatter_filter(
    std::string const &filter_expr,
    DataFormatterDictionary<DataObject> const &_dict);

template <typename DataObject>
DataFormatterFilter<DataObject> make_data_formatter_filter(
    std::vector<std::string> const &filter_expr,
    DataFormatterDictionary<DataObject> const &_dict);

}  // namespace CASM

#endif
