#ifndef CASM_SafeOfstream
#define CASM_SafeOfstream

#include <string>
#include <vector>
#include "casm/external/boost.hh"
#include "casm/CASM_global_definitions.hh"

namespace CASM {

  /// \brief Write to a temporary file to ensure a good write, then rename
  ///
  /// \code
  /// SafeOfstream file;
  /// jsonParser json;
  ///
  /// // by default, write temporary file "something.json.tmp", or specify your own extension
  //  // throws if "something.json.tmp" already exists
  /// file.open("something.json", "tmp");
  ///
  /// // write to the underlying stream
  /// json.write(file.ofstream());
  ///
  /// // 'close' results in closing "something.json.tmp", checking fs::ofstream::fail(), and if not failed,
  /// // then removing "something.json", and renaming "something.json.tmp" -> "something.json"
  /// file.close();
  /// \endcode
  class SafeOfstream {

  public:

    SafeOfstream() {}

    /// \brief Opens "file.tmp" for writing, with intended final target "file"
    ///
    /// \param name Name of target file
    /// \param tmp_ext String to be used as an extension for the temporary file that is written
    ///
    /// \throws if name.tmp_ext already exists
    ///
    /// Example:
    /// \code
    /// SafeOfstream file;
    /// jsonParser json;
    ///
    /// // By default, open temporary file "something.json.tmp" for writing.
    /// // - throws if "something.json.tmp" already exists
    /// file.open("something.json");
    ///
    /// // or specify your own extension, to open temporary file "something.json.backup"
    /// // - throws if "something.json.backup" already exists
    /// // file.open("something.json", "backup");
    ///
    /// // print to the underlying stream
    /// json.print(file.ofstream());
    ///
    /// // 'close' results in closing ofstream to "something.json.tmp", checking fs::ofstream::fail(),
    /// // and if not failed then removing "something.json", and renaming "something.json.tmp" -> "something.json"
    /// file.close();
    /// \endcode
    ///
    void open(fs::path name, std::string tmp_ext = "tmp") {
      m_name = name;
      m_tmp_name = m_name.string() + "." + tmp_ext;

      if(fs::exists(m_tmp_name)) {
        throw std::runtime_error(
          std::string("Error in 'SafeOfstream::open(fs::path name, std::string tmp_ext)'.\n") +
          "  File: " + m_tmp_name.string() + " already exists");
      }

      m_sout.open(m_tmp_name);
    }

    /// \brief Access underlying stream
    fs::ofstream &ofstream() {
      return m_sout;
    }

    /// \brief Closes stream, and if not a failed write, removes "file" and renames "file.tmp" to "file"
    void close() {
      m_sout.close();
      if(!m_sout.fail()) {
        fs::remove(m_name);
        fs::rename(m_tmp_name, m_name);
      }
    }

  private:

    fs::path m_name;
    fs::path m_tmp_name;
    fs::ofstream m_sout;


  };


}

#endif
