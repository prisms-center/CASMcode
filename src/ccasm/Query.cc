#include "ccasm/Query.hh"
#include <wordexp.h>
#include "boost/iostreams/stream.hpp"
#include "boost/iostreams/device/null.hpp"
#include "casm/clex/PrimClex.hh"
#include "casm/external/boost.hh"
#include "casm/app/query.hh"

using namespace CASM;

extern "C" {

  costream *STDOUT() {
    return reinterpret_cast<costream *>(&std::cout);
  }

  costream *STDERR() {
    return reinterpret_cast<costream *>(&std::cerr);
  }


  costream *nullstream_new() {
    using namespace boost::iostreams;
    return reinterpret_cast<costream *>(new stream<null_sink>(null_sink()));
  }

  void nullstream_delete(costream *ptr) {
    using namespace boost::iostreams;
    delete reinterpret_cast<stream<null_sink>*>(ptr);
  }


  costream *ostringstream_new() {
    return reinterpret_cast<costream *>(new std::ostringstream());
  }

  void ostringstream_delete(costream *ptr) {
    delete reinterpret_cast<std::ostringstream *>(ptr);
  }

  unsigned long ostringstream_size(costream *ptr) {
    return reinterpret_cast<std::ostringstream *>(ptr)->tellp();
  }

  char *ostringstream_strcpy(costream *ptr, char *c_str) {
    auto str = reinterpret_cast<std::ostringstream *>(ptr)->str();
    std::strcpy(c_str, str.c_str());
    return c_str;
  }


  cPrimClex *primclex_new(char *path, costream *sout) {
    std::ostream &_sout = *reinterpret_cast<std::ostream *>(sout);
    PrimClex *ptr = new PrimClex(fs::path(path), _sout);
    return reinterpret_cast<cPrimClex *>(ptr);
  }

  void primclex_delete(cPrimClex *ptr) {
    delete reinterpret_cast<PrimClex *>(ptr);
  }

  void primclex_check(cPrimClex *ptr) {
    std::cout << reinterpret_cast<PrimClex *>(ptr)->get_path() << std::endl;
  }

  int query(char *args, cPrimClex *_primclex, costream *sout, costream *serr) {
    PrimClex *ptr = reinterpret_cast<PrimClex *>(_primclex);
    std::ostream &_sout = *reinterpret_cast<std::ostream *>(sout);
    std::ostream &_serr = *reinterpret_cast<std::ostream *>(serr);

    fs::path curr = fs::current_path();

    // parse args -> argc, argv
    wordexp_t p;
    int res = wordexp(args, &p, 0);
    if(res) {
      _serr << "Error parsing query: '" << args << "'" << std::endl;
      _serr << "wordexp() error: " << res << std::endl;
      switch(res) {
      case 1: {
        _serr << "Check for illegal unescaped characters: |, &, ;, <, >, (, ), {, }" << std::endl;
        break;
      }
      default: {
        _serr << "Check 'man wordexp' for error code meaning" << std::endl;
      }
      }
      return res;
    }
    res = query_command(p.we_wordc, p.we_wordv, ptr, _sout, _serr);
    wordfree(&p);

    fs::current_path(curr);

    return res;
  }

}