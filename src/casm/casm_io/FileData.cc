#include "casm/casm_io/FileData.hh"

#include <boost/filesystem.hpp>

#include "casm/global/definitions.hh"

namespace CASM {

bool FileData::exists() const { return fs::exists(path()); }

void FileData::refresh() {
  m_timestamp = 0;
  if (this->exists()) {
    m_timestamp = fs::last_write_time(path());
  }
}
}  // namespace CASM
