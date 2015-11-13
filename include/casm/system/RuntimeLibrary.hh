#ifndef RuntimeLibrary_HH
#define RuntimeLibrary_HH

#include <iostream>
#include <fstream>
#include <string>
#include <functional>
#include <dlfcn.h>
#include <cstdlib>
#include <boost/filesystem.hpp>
#include "casm/system/Popen.hh"

namespace CASM {

  /// \brief Write, compile, load and use code at runtime
  class RuntimeLibrary {

  public:

    /// \brief Construct a RuntimeLibrary object, with the options to be used for compile
    ///        the '.o' file and the '.so' file
    RuntimeLibrary(std::string _compile_options = RuntimeLibrary::default_compile_options(),
                   std::string _so_options = RuntimeLibrary::default_so_options()) :
      m_compile_options(_compile_options),
      m_so_options(_so_options),
      m_filename_base(""),
      m_handle(nullptr) {}

    ~RuntimeLibrary() {
      if(m_handle != nullptr) {
        close();
      }
    }

    /// \brief Compile a shared library
    ///
    /// \param _filename_base Base name for the source code file. For example, "hello" results in writing "hello.cc",
    ///        and compiling "hello.o" and "hello.so" in the current working directory.
    /// \param _source A std::string containing the source code to be written. For example,
    /// \code
    /// std::string cc_file;
    ///
    /// cc_file = std::string("#include <iostream>\n") +
    ///           "extern \"C\" int hello() {\n" +
    ///           "   std::cout << \"Hello, my name is Ultron. I'm here to protect you.\" << '\\n';\n" +
    ///           "   return 42;\n" +
    ///           "}\n";
    /// \endcode
    ///
    /// \result Writes a file "example.cc", and compiles an object file and shared library using the options
    ///         provided when this RuntimeLibrary object was constructed. By default, "example.o" and "example.so".
    ///
    /// To enable runtime symbol lookup use C-style functions, i.e use extern "C" for functions you want to use
    /// via get_function.  This means no member functions or overloaded functions.
    ///
    void compile(std::string _filename_base,
                 std::string _source) {

      if(m_handle != nullptr) {
        close();
      }

      m_filename_base = _filename_base;

      // write the source code
      std::ofstream file(m_filename_base + ".cc");
      file << _source;
      file.close();

      // compile the source code into a dynamic library
      Popen p;
      p.popen(m_compile_options + " -o " + m_filename_base + ".o" + " -c " + m_filename_base + ".cc");
      p.popen(m_so_options + " -o " + m_filename_base + ".so" + " " + m_filename_base + ".o");
    }

    /// \brief Compile a shared library
    ///
    /// \param _filename_base Base name for the source code file. For example, "/path/to/hello" looks for "/path/to/hello.cc",
    ///        and compile "/path/to/hello.o" and "/path/to/hello.so".
    ///
    /// \result Compiles file "/path/to/hello.cc" into an object file and shared library using the options
    ///         provided when this RuntimeLibrary object was constructed. By default, "example.o" and "example.so".
    ///
    /// To enable runtime symbol lookup use C-style functions, i.e use extern "C" for functions you want to use
    /// via get_function.  This means no member functions or overloaded functions.
    ///
    void compile(std::string _filename_base) {
      if(m_handle != nullptr) {
        close();
      }

      m_filename_base = _filename_base;

      // compile the source code into a dynamic library
      Popen p;
      p.popen(m_compile_options + " -o " + m_filename_base + ".o" + " -c " + m_filename_base + ".cc");
      p.popen(m_so_options + " -o " + m_filename_base + ".so" + " " + m_filename_base + ".o");
    }


    /// \brief Load a library with a given name
    ///
    /// \param _filename_base For "hello", this loads "hello.so"
    ///
    void load(std::string _filename_base) {

      if(m_handle != nullptr) {
        close();
      }

      m_filename_base = _filename_base;

      m_handle = dlopen((m_filename_base + ".so").c_str(), RTLD_NOW);
      if(!m_handle) {
        throw std::runtime_error(std::string("Cannot open library: ") + m_filename_base + ".so");;
      }
    }

    /// \brief Obtain a function from the current library
    ///
    /// Must be a C-style function to enable symbol lookup, i.e your source code should use extern "C".
    /// This means no member functions or overloaded functions.
    ///
    template<typename Signature>
    std::function<Signature> get_function(std::string function_name) const {

      dlerror();

      std::function<Signature> func = reinterpret_cast<Signature *>(dlsym(m_handle, function_name.c_str()));

      const char *dlsym_error = dlerror();
      if(dlsym_error) {
        throw std::runtime_error(std::string("Cannot load symbol " + function_name + " \n" + dlsym_error));
      }

      return func;
    }

    /// \brief Close the current library
    ///
    /// This is also done on destruction.
    void close() const {
      // close
      dlclose(m_handle);
    }

    /// \brief Remove the current library and source code
    void rm() const {
      if(m_filename_base == "") {
        return;
      }

      // rm
      Popen p;
      p.popen(std::string("rm -f ") + m_filename_base + ".cc " + m_filename_base + ".o " + m_filename_base + ".so");
    }

    /// \brief Default compilation options
    ///
    /// \returns "g++ -O3 -Wall -fPIC"
    static std::string default_compile_options() {
      
      return cxx() + " -O3 -Wall -fPIC --std=c++11 " + include_casmprefix();
    }

    /// \brief Default shared library options
    ///
    /// \returns "g++ -shared"
    static std::string default_so_options() {
      return cxx() + " -shared";
    }
    
    /// \brief Return default compiler
    ///
    /// - if environment variable CXX exists, uses that, otherwise "g++"
    static std::string cxx() {
      std::string result = "g++";
      char* CXX = std::getenv("CXX");
      if(CXX != nullptr) {
        result = std::string(CXX);
      }
      return result;
    }
    
    /// \brief Return CASMPREFIX include path option
    ///
    /// - if environment variable CASMPREFIX exists, returns "-I$CASMPREFIX/include", otherwise ""
    static std::string include_casmprefix() {
      std::string result = "";
      char* CASMPREFIX = std::getenv("CASMPREFIX");
      if(CASMPREFIX != nullptr) {
        result = "-I" + (boost::filesystem::path(CASMPREFIX) / "include").string();
      }
      return result;
    }

  private:

    std::string m_compile_options;
    std::string m_so_options;

    std::string m_filename_base;

    void *m_handle;

  };
}

#endif
