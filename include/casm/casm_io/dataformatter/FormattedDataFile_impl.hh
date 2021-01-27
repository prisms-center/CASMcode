#ifndef CASM_FormattedDataFile_impl
#define CASM_FormattedDataFile_impl

#include <boost/filesystem/fstream.hpp>

#include "casm/casm_io/dataformatter/FormattedDataFile.hh"
#include "casm/external/gzstream/gzstream.h"

namespace CASM {

template <typename DataObject>
FormattedDataFile<DataObject>::FormattedDataFile(
    FormattedDataFileOptions const &options)
    : FormattedDataFile(options.file_path, options.json_output,
                        options.json_arrays, options.compress) {}

template <typename DataObject>
FormattedDataFile<DataObject>::FormattedDataFile(fs::path file_path,
                                                 bool json_output,
                                                 bool json_arrays,
                                                 bool compress)
    : m_initialized(false),
      m_json_output(json_output),
      m_json_arrays(json_arrays) {
  if (!compress) {
    m_ostream = notstd::make_unique<fs::ofstream>(file_path);
  } else {
    m_ostream = notstd::make_unique<gz::ogzstream>(file_path.string().c_str());
  }
}

template <typename DataObject>
FormattedDataFile<DataObject>::~FormattedDataFile() {
  if (m_json_output) {
    *m_ostream << m_json;
  }
}

template <typename DataObject>
void FormattedDataFile<DataObject>::operator()(
    DataFormatter<DataObject> const &formatter, DataObject const &object) {
  if (!m_initialized) {
    if (!m_json_output) {
      // do nothing. headers will be printed... somehow
    } else if (m_json_arrays) {
      m_json = jsonParser::object();
    } else {
      m_json = jsonParser::array();
    }
    m_initialized = true;
  }
  if (!m_json_output) {
    *m_ostream << formatter(object);
  } else if (m_json_arrays) {
    formatter.to_json_arrays(object, m_json);
  } else {
    // formatter.to_json(object, m_json); // TODO: why doesn't this work?
    m_json.push_back(formatter(object));
  }
}

}  // namespace CASM

#endif
