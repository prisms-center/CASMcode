#include "ccasm/api.hh"
#include <wordexp.h>
#include "boost/iostreams/stream.hpp"
#include "boost/iostreams/device/null.hpp"
#include "casm/clex/PrimClex.hh"
#include "casm/external/boost.hh"
#include "casm/app/casm_functions.hh"

using namespace CASM;

extern "C" {

  costream *casm_STDOUT() {
    return reinterpret_cast<costream *>(&std::cout);
  }

  costream *casm_STDERR() {
    return reinterpret_cast<costream *>(&std::cerr);
  }


  costream *casm_nullstream_new() {
    using namespace boost::iostreams;
    return reinterpret_cast<costream *>(new stream<null_sink>(null_sink()));
  }

  void casm_nullstream_delete(costream *ptr) {
    using namespace boost::iostreams;
    delete reinterpret_cast<stream<null_sink>*>(ptr);
  }


  costream *casm_ostringstream_new() {
    return reinterpret_cast<costream *>(new std::ostringstream());
  }

  void casm_ostringstream_delete(costream *ptr) {
    delete reinterpret_cast<std::ostringstream *>(ptr);
  }

  unsigned long casm_ostringstream_size(costream *ptr) {
    return reinterpret_cast<std::ostringstream *>(ptr)->tellp();
  }

  char *casm_ostringstream_strcpy(costream *ptr, char *c_str) {
    auto str = reinterpret_cast<std::ostringstream *>(ptr)->str();
    std::strcpy(c_str, str.c_str());
    return c_str;
  }


  cPrimClex *casm_primclex_new(char *path, costream *log) {
    Log _log(*reinterpret_cast<std::ostream*>(log));
    PrimClex *ptr = new PrimClex(fs::path(path), _log);
    return reinterpret_cast<cPrimClex *>(ptr);
  }

  void casm_primclex_delete(cPrimClex *ptr) {
    delete reinterpret_cast<PrimClex *>(ptr);
  }

  
  int casm_capi(char *args, cPrimClex *primclex, costream *log, costream *err_log) {
    PrimClex *_primclex = reinterpret_cast<PrimClex *>(primclex);
    Log _log(*reinterpret_cast<std::ostream*>(log));
    Log _err_log(*reinterpret_cast<std::ostream*>(err_log));
    
    std::string s("casm ");
    s += std::string(args);
    
    // parse args -> argc, argv
    wordexp_t p;
    int res = wordexp(s.c_str(), &p, 0);
    if(res) {
      _err_log << "Error parsing query: '" << args << "'" << std::endl;
      _err_log << "wordexp() error: " << res << std::endl;
      switch(res) {
      case 1: {
        _err_log << "Check for illegal unescaped characters: |, &, ;, <, >, (, ), {, }" << std::endl;
        break;
      }
      default: {
        _err_log << "Check 'man wordexp' for error code meaning" << std::endl;
      }
      }
      return res;
    }
    
    CommandArgs command_args(p.we_wordc, p.we_wordv, _primclex, _log, _err_log);
    res = casm_api(command_args);
    
    wordfree(&p);

    return res;
  }

}