#include <iostream>
#include <sstream>

/// For CASM::PrimClex*
typedef struct cPrimClex cPrimClex;

/// For std::ostream*
typedef struct costream costream;

/// For CASM::Log*
typedef struct cLog cLog;

extern "C" {

  cLog *log_new(costream *ptr, int verbosity, bool show_clock);
  
  void log_delete(cLog *ptr);
  
  
  costream *STDOUT();

  costream *STDERR();


  costream *nullstream_new();

  void nullstream_delete(costream *ptr);


  costream *ostringstream_new();

  void ostringstream_delete(costream *ptr);

  unsigned long ostringstream_size(costream *ptr);

  char *ostringstream_strcpy(costream *ptr, char *c_str);


  cPrimClex *primclex_new(char *path, cLog *log);

  void primclex_delete(cPrimClex *ptr);

  void primclex_check(cPrimClex *ptr);

  int query(char *args, cPrimClex *primclex, cLog *log, costream *serr);

}