#include "casm/system/RuntimeLibrary.hh"

namespace CASM {
  
  /// \brief Construct a RuntimeLibrary object, with the options to be used for compile
  ///        the '.o' file and the '.so' file
  RuntimeLibrary::RuntimeLibrary(std::string _compile_options,
                                 std::string _so_options) :
    m_compile_options(_compile_options),
    m_so_options(_so_options),
    m_filename_base(""),
    m_handle(nullptr) {}

  RuntimeLibrary::~RuntimeLibrary() {
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
  void RuntimeLibrary::compile(std::string _filename_base, std::string _source) {

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
  void RuntimeLibrary::compile(std::string _filename_base) {
    if(m_handle != nullptr) {
      close();
    }

    m_filename_base = _filename_base;

    // compile the source code into a dynamic library
    Popen p;
    std::string cmd = m_compile_options + " -o " + m_filename_base + ".o" + " -c " + m_filename_base + ".cc";
    p.popen(cmd);
    if(p.exit_code()) {
      std::cerr << "Error compiling: " << m_filename_base + ".cc" << std::endl;
      std::cerr << "Attempted: " << cmd << std::endl;
      std::cerr << p.gets() << std::endl;
      throw std::runtime_error("Can not compile " + m_filename_base + ".cc");
    }
    
    cmd = m_so_options + " -o " + m_filename_base + ".so" + " " + m_filename_base + ".o";
    p.popen(cmd);
    if(p.exit_code()) {
      std::cerr << "Error compiling shared object: " << m_filename_base + ".so" << std::endl;
      std::cerr << "Attempted: " << cmd << std::endl;
      std::cerr << p.gets() << std::endl;
      throw std::runtime_error("Can not compile " + m_filename_base + ".o");
    }
  }


  /// \brief Load a library with a given name
  ///
  /// \param _filename_base For "hello", this loads "hello.so"
  ///
  void RuntimeLibrary::load(std::string _filename_base) {

    if(m_handle != nullptr) {
      close();
    }

    m_filename_base = _filename_base;

    m_handle = dlopen((m_filename_base + ".so").c_str(), RTLD_NOW);
    if(!m_handle) {
      throw std::runtime_error(std::string("Cannot open library: ") + m_filename_base + ".so");;
    }
  }

  /// \brief Close the current library
  ///
  /// This is also done on destruction.
  void RuntimeLibrary::close() {
    // close
    if(m_handle != nullptr && m_filename_base != "") {
      dlclose(m_handle);
    }
  }

  /// \brief Remove the current library and source code
  void RuntimeLibrary::rm() {
    if(m_filename_base == "") {
      return;
    }

    // rm
    Popen p;
    p.popen(std::string("rm -f ") + m_filename_base + ".cc " + m_filename_base + ".o " + m_filename_base + ".so");
    
    m_filename_base = "";
  }

  /// \brief Default compilation command
  ///
  /// \returns "$CXX -O3 -Wall -fPIC --std=c++11 $CASM_INCLUDE"
  ///
  /// $CXX and $CASM_INCLUDE depend on current environment variables:
  /// - $CXX is replaced with "$CXX" if CXX exists and "g++" otherwise
  /// - $CASM_INCLUDE is replaced with "-I$CASMPREFIX/include" if CASMPREFIX exists
  std::string RuntimeLibrary::default_compile_options() {
    
    return cxx() + " " + default_cxxflags() + " " + casm_include();
  }
  
  /// \brief Default c++ compiler options
  ///
  /// \returns "-O3 -Wall -fPIC --std=c++11"
  std::string RuntimeLibrary::default_cxxflags() {
    return "-O3 -Wall -fPIC --std=c++11";
  }

  /// \brief Default shared library options
  ///
  /// \returns "$CXX -shared"
  std::string RuntimeLibrary::default_so_options() {
    return cxx() + " -shared" + " -lboost_system";
  }
  
  /// \brief Return default compiler
  ///
  /// - if environment variable CXX exists, uses that, otherwise "g++"
  std::string RuntimeLibrary::cxx() {
    std::string result = "g++";
    char* CXX = std::getenv("CXX");
    if(CXX != nullptr) {
      result = std::string(CXX);
    }
    return result;
  }
  
  /// \brief Return include path option for CASM
  ///
  /// \returns "-I$HHCASM/include" if environment variable HHCASM exists,
  ///          "-I$CASMPREFIX/include" if environment variable CASMPREFIX exists, 
  ///          otherwise an empty string
  std::string RuntimeLibrary::casm_include() {
    std::string result = "";
    
    char* HHCASM = std::getenv("HHCASM");
    if(HHCASM != nullptr) {
      result = "-I" + (boost::filesystem::path(HHCASM) / "include").string();
    }
    else {
      char* CASMPREFIX = std::getenv("CASMPREFIX");
      if(CASMPREFIX != nullptr) {
        result = "-I" + (boost::filesystem::path(CASMPREFIX) / "include").string();
      }
    }
    return result;
  }
  
}