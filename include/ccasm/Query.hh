#include <iostream>
#include <sstream>

/// For CASM::PrimClex*
typedef struct cPrimClex cPrimClex;

/// For std::ostream*
typedef struct costream costream;

/// For CASM::Log*
typedef struct cLog cLog;

extern "C" {

  cLog *casm_log_new(costream *ptr, int verbosity, bool show_clock);
  
  void casm_log_delete(cLog *ptr);
  
  
  costream *casm_STDOUT();

  costream *casm_STDERR();


  costream *casm_nullstream_new();

  void casm_nullstream_delete(costream *ptr);


  costream *casm_ostringstream_new();

  void casm_ostringstream_delete(costream *ptr);

  unsigned long casm_ostringstream_size(costream *ptr);

  char *casm_ostringstream_strcpy(costream *ptr, char *c_str);


  cPrimClex *casm_primclex_new(char *path, cLog *log);

  void casm_primclex_delete(cPrimClex *ptr);

  void casm_primclex_check(cPrimClex *ptr);

  
  int casm_cmain(char *args, cPrimClex *primclex, cLog *log, cLog *err_log);

}