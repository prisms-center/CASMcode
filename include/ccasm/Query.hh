#include <iostream>
#include <sstream>

/// For CASM::PrimClex*
typedef struct cPrimClex cPrimClex;

/// For std::ostream*
typedef struct costream costream;

extern "C" {

  costream *STDOUT();

  costream *STDERR();


  costream *nullstream_new();

  void nullstream_delete(costream *ptr);


  costream *ostringstream_new();

  void ostringstream_delete(costream *ptr);

  unsigned long ostringstream_size(costream *ptr);

  char *ostringstream_strcpy(costream *ptr, char *c_str);


  cPrimClex *primclex_new(char *path, costream *sout);

  void primclex_delete(cPrimClex *ptr);

  void primclex_check(cPrimClex *ptr);

  int query(char *args, cPrimClex *primclex, costream *sout, costream *serr);

}