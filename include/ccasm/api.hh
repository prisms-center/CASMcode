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


  costream *casm_ostringstream_new();

  void casm_ostringstream_delete(costream *ptr);

  unsigned long casm_ostringstream_size(costream *ptr);

  char *casm_ostringstream_strcpy(costream *ptr, char *c_str);


  cPrimClex *casm_primclex_new(char *path, costream *log, costream *debug_log, costream *err_log);

  void casm_primclex_delete(cPrimClex *ptr);

  void casm_primclex_refresh(cPrimClex *ptr,
                             bool read_settings,
                             bool read_composition,
                             bool read_chem_ref,
                             bool read_configs,
                             bool clear_clex);


  int casm_capi(char *args, cPrimClex *primclex, costream *log, costream *debug_log, costream *err_log);

}