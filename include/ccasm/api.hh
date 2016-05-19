#include <iostream>
#include <sstream>

/// For CASM::PrimClex*
typedef struct cPrimClex cPrimClex;

/// For std::ostream*
typedef struct costream costream;


extern "C" {

  costream *casm_STDOUT();

  costream *casm_STDERR();


  costream *casm_nullstream_new();

  void casm_nullstream_delete(costream *ptr);


  costream *casm_ostringstream_new();

  void casm_ostringstream_delete(costream *ptr);

  unsigned long casm_ostringstream_size(costream *ptr);

  char *casm_ostringstream_strcpy(costream *ptr, char *c_str);


  cPrimClex *casm_primclex_new(char *path, costream *log);

  void casm_primclex_delete(cPrimClex *ptr);
  
  
  int casm_capi(char *args, cPrimClex *primclex, costream *log, costream *err_log);

}