#include <iostream>
#include <sstream>

/** \ingroup API
 *  @{
 */

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


  cPrimClex *casm_primclex_null();

  cPrimClex *casm_primclex_new(char *path, costream *log, costream *debug_log, costream *err_log);

  void casm_primclex_delete(cPrimClex *ptr);

  void casm_primclex_refresh(cPrimClex *ptr,
                             bool read_settings,
                             bool read_composition,
                             bool read_chem_ref,
                             bool read_configs,
                             bool clear_clex);

  void casm_primclex_set_logging(cPrimClex *primclex, costream *log, costream *debug_log, costream *err_log);

  void casm_command_list(costream *ostringstream_log);

  int casm_capi(char *args, cPrimClex *primclex, char *root, costream *log, costream *debug_log, costream *err_log);

  int casm_capi_call(char *args, cPrimClex *primclex);
}

/** @} */
