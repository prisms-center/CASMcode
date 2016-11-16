#include "casm/system/RuntimeLibrary.hh"
#include "casm/casm_io/Log.hh"

namespace CASM {

  /// \brief Construct a RuntimeLibrary object, with the options to be used for compile
  ///        the '.o' file and the '.so' file
  RuntimeLibrary::RuntimeLibrary(std::string filename_base,
                                 std::string compile_options,
                                 std::string so_options,
                                 std::string compile_msg,
                                 const Logging &logging) :
    Logging(logging),
    m_filename_base(filename_base),
    m_compile_options(compile_options),
    m_so_options(so_options),
    m_handle(nullptr) {
    
    // If the shared library doesn't exist
    if(!fs::exists(m_filename_base + ".so")) {

      // But the library source code does
      if(fs::exists(m_filename_base + ".cc")) {

        // Compile it
        try {
          log().compiling<Log::standard>(m_filename_base + ".cc");
          log().begin_lap();
          log() << compile_msg << std::endl;
          _compile();
          log() << "compile time: " << log().lap_time() << " (s)\n" << std::endl;
        }
        catch(std::exception &e) {
          log() << "Error compiling clexulator. To fix: \n";
          log() << "  - Check compiler error messages.\n";
          log() << "  - Check compiler options with 'casm settings -l'\n";
          log() << "    - Update compiler options with 'casm settings --set-compile-options '...options...'\n";
          log() << "    - Make sure the casm headers can be found by including '-I/path/to/casm'\n";
          log() << "  - The default compiler is 'g++'. Override by setting the environment variable CXX\n" << std::endl;
          throw e;
        }
      }
      else {
        throw std::runtime_error(
          std::string("Error in RuntimeLibrary\n") +
          "  Could not find '" + m_filename_base + ".so' or '" + m_filename_base + ".cc'");
      }
    }

    // If the shared library exists
    if(fs::exists(m_filename_base + ".so")) {

      // Load the library with the Clexulator
      _load();

    }
    else {
      throw std::runtime_error(
        std::string("Error in Clexulator constructor\n") +
        "  Did not find '" + m_filename_base + ".so'");
    }

  }

  RuntimeLibrary::~RuntimeLibrary() {
    if(m_handle != nullptr) {
      _close();
    }
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
  void RuntimeLibrary::_compile() {

    // compile the source code into a dynamic library
    Popen p;
    std::string cmd = m_compile_options + " -o " + m_filename_base + ".o" + " -c " + m_filename_base + ".cc";
    p.popen(cmd);
    if(p.exit_code()) {
      err_log() << "Error compiling: " << m_filename_base + ".cc" << std::endl;
      err_log() << "Attempted: " << cmd << std::endl;
      err_log() << p.gets() << std::endl;
      throw std::runtime_error("Can not compile " + m_filename_base + ".cc");
    }

    cmd = m_so_options + " -o " + m_filename_base + ".so" + " " + m_filename_base + ".o";
    p.popen(cmd);
    if(p.exit_code()) {
      err_log() << "Error compiling shared object: " << m_filename_base + ".so" << std::endl;
      err_log() << "Attempted: " << cmd << std::endl;
      err_log() << p.gets() << std::endl;
      throw std::runtime_error("Can not compile " + m_filename_base + ".o");
    }
  }

  /// \brief Load a library with a given name
  ///
  /// \param _filename_base For "hello", this loads "hello.so"
  ///
  void RuntimeLibrary::_load() {

    m_handle = dlopen((m_filename_base + ".so").c_str(), RTLD_NOW);
    if(!m_handle) {
      fprintf(stderr, "dlopen failed: %s\n", dlerror());
      throw std::runtime_error(std::string("Cannot open library: ") + m_filename_base + ".so");
    }
  }

  /// \brief Close the current library
  ///
  /// This is also done on destruction.
  void RuntimeLibrary::_close() {
    // close
    if(m_handle != nullptr) {
      dlclose(m_handle);
    }
  }

  /// \brief Remove the current library and source code
  void RuntimeLibrary::rm() {
    _close();
    // rm
    Popen p;
    p.popen(std::string("rm -f ") + m_filename_base + ".cc " + m_filename_base + ".o " + m_filename_base + ".so");
  }

  namespace {

    std::vector<std::string> _cxx_env() {
      return std::vector<std::string> {
        "CASM_CXX",
        "CXX"
      };
    }

    std::vector<std::string> _cxxflags_env() {
      return std::vector<std::string> {
        "CASM_CXXFLAGS",
      };
    }

    std::vector<std::string> _soflags_env() {
      return std::vector<std::string> {
        "CASM_SOFLAGS"
      };
    }

    std::vector<std::string> _casm_env() {
      return std::vector<std::string> {
        "CASM_PREFIX"
      };
    }

    std::vector<std::string> _boost_env() {
      return std::vector<std::string> {
        "CASM_BOOST_PREFIX"
      };
    }

    /// \brief Some function of environment variables
    std::pair<std::string, std::string> _use_env(std::vector<std::string> var, std::string _default = "") {
      for(const auto &v : var) {
        char *_env = std::getenv(v.c_str());
        if(_env != nullptr) {
          return std::make_pair(std::string(_env), v);
        }
      }
      return std::make_pair(_default, "default");
    }
  }

  /// \brief Return default compiler and specifying variable
  ///
  /// \returns "$CASM_CXX" if environment variable CASM_CXX exists,
  ///          "$CXX" if environment variable CXX exists,
  ///          otherwise "g++"
  std::pair<std::string, std::string> RuntimeLibrary::default_cxx() {
    return _use_env(_cxx_env(), "g++");
  }

  /// \brief Default c++ compiler options
  ///
  /// \returns "-O3 -Wall -fPIC --std=c++11"
  std::pair<std::string, std::string> RuntimeLibrary::default_cxxflags() {
    return _use_env(_cxxflags_env(), "-O3 -Wall -fPIC --std=c++11");
  }

  /// \brief Default c++ shared library options
  ///
  /// \returns "-shared -lboost_system"
  std::pair<std::string, std::string> RuntimeLibrary::default_soflags() {
    return _use_env(_soflags_env(), "-shared -lboost_system");
  }

  /// \brief Return include path option for CASM
  ///
  /// \returns $CASM_PREFIX if environment variable CASM_PREFIX exists,
  ///          otherwise "/usr/local"
  std::pair<fs::path, std::string> RuntimeLibrary::default_casm_prefix() {
    auto res = _use_env(_casm_env(), "/usr/local");
    return std::make_pair(fs::path(res.first), res.second);
  }

  /// \brief Return include path option for boost
  ///
  /// \returns $CASM_BOOST_PREFIX if environment variable CASM_BOOST_PREFIX exists,
  ///          otherwise an empty string
  std::pair<fs::path, std::string> RuntimeLibrary::default_boost_prefix() {
    auto res = _use_env(_boost_env());
    return std::make_pair(fs::path(res.first), res.second);
  }

  std::string include_path(const fs::path &prefix) {
    if(!prefix.empty()) {
      return "-I" + (prefix / "include").string();
    }
    return "";
  };

  std::string link_path(const fs::path &prefix) {
    if(!prefix.empty()) {
      return "-L" + (prefix / "lib").string();
    }
    return "";
  };

}