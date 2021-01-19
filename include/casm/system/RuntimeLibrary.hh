#ifndef RuntimeLibrary_HH
#define RuntimeLibrary_HH

#include <dlfcn.h>

#include <cstdlib>
#include <fstream>
#include <functional>
#include <iostream>
#include <string>

#include "casm/global/definitions.hh"

namespace CASM {

class runtime_lib_compile_error : public std::runtime_error {
 public:
  runtime_lib_compile_error(std::string _filename_base, std::string _cmd,
                            std::string _result, std::string _what);

  std::string filename_base;
  std::string cmd;
  std::string result;

  void print(std::ostream &sout) const;
};

class runtime_lib_shared_error : public std::runtime_error {
 public:
  runtime_lib_shared_error(std::string _filename_base, std::string _cmd,
                           std::string _result, std::string _what);

  std::string filename_base;
  std::string cmd;
  std::string result;

  void print(std::ostream &sout) const;
};

/// \brief Write, compile, load and use code at runtime
class RuntimeLibrary {
 public:
  /// \brief Construct a RuntimeLibrary object, with the options to be used for
  /// compile
  ///        the '.o' file and the '.so' file
  RuntimeLibrary(std::string _filename_base, std::string _compile_options,
                 std::string _so_options);

  ~RuntimeLibrary();

  /// \brief Obtain a function from the current library
  ///
  /// Must be a C-style function to enable symbol lookup, i.e your source code
  /// should use extern "C". This means no member functions or overloaded
  /// functions.
  ///
  template <typename Signature>
  std::function<Signature> get_function(std::string function_name) const {
    std::function<Signature> func =
        reinterpret_cast<Signature *>(dlsym(m_handle, function_name.c_str()));

    const char *dlsym_error = dlerror();
    if (dlsym_error) {
      throw std::runtime_error(std::string(
          "Cannot load symbol " + function_name + " \n" + dlsym_error));
    }

    return func;
  }

  /// \brief Remove the current library and source code
  void rm();

  /// \brief Default c++ compiler options
  static std::pair<std::string, std::string> default_cxxflags();

  /// \brief Default c++ compiler options
  static std::pair<std::string, std::string> default_soflags();

  /// \brief Return default compiler
  static std::pair<std::string, std::string> default_cxx();

  /// \brief Return default includedir for CASM
  static std::pair<fs::path, std::string> default_casm_includedir();

  /// \brief Return default libdir for CASM
  static std::pair<fs::path, std::string> default_casm_libdir();

  /// \brief Return default includedir for boost
  static std::pair<fs::path, std::string> default_boost_includedir();

  /// \brief Return default libdir for boost
  static std::pair<fs::path, std::string> default_boost_libdir();

 private:
  /// \brief Compile a shared library
  void _compile();

  /// \brief Load a library with a given name
  void _load();

  /// \brief Close the current library
  void _close();

  std::string m_filename_base;
  std::string m_compile_options;
  std::string m_so_options;

  void *m_handle;
};

std::string include_path(const fs::path &dir);

std::string link_path(const fs::path &dir);

}  // namespace CASM

#endif
