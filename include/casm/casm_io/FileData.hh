#include <ctime>
#include <string>

#ifndef CASM_FileData_HH
#define CASM_FileData_HH
namespace CASM {
/// \brief Interface class to check and/or store path and last_write_time of
/// file on disk. Used to determine if files have been updated or deleted
class FileData {
 public:
  /// \brief Construct with path and explicit timestamp
  FileData(std::string const &_path, std::time_t _timestamp)
      : m_path(_path), m_timestamp(_timestamp) {}

  ///\brief Construct with path, timestamp is set if file exists at path
  FileData(std::string const &_path) : m_path(_path), m_timestamp(0) {
    this->refresh();
  }

  ///\brief default constructor
  FileData() : FileData("", 0) {}

  /// \brief return stored path
  std::string const &path() const { return m_path; }

  /// \brief return stored timestamp. Timestamp corresponds to last_write_time
  /// of file timestamp defaults to 0 if file does not exist.
  std::time_t timestamp() const { return m_timestamp; }

  /// \brief checks timestamp of file at path. Returns true if last_write_time
  /// of file matches timestamp returns false otherwise, or if file doesn't
  /// exist and stored timestamp is not zero (indicated file has been deleted)
  bool up_to_date() const { return (*this) == FileData(path()); }

  /// \brief Returns true if path is empty
  bool empty() const { return m_path.empty(); }

  /// \brief checks path for existing file. Does not alter stored timestamp
  bool exists() const;

  /// \brief Updates timestamp using stored path
  void refresh();

  /// \brief returns true if other FileData has same stored path and timestamp
  bool operator==(FileData const &other) const {
    return this->timestamp() == other.timestamp() &&
           this->path() == other.path();
  }

 private:
  std::string m_path;
  std::time_t m_timestamp;
};
}  // namespace CASM
#endif
