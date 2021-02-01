#ifndef CASM_clex_ChemicalReference_file_io
#define CASM_clex_ChemicalReference_file_io

#include "casm/global/definitions.hh"

namespace CASM {

namespace xtal {
class BasicStructure;
}
class ChemicalReference;

ChemicalReference read_chemical_reference(fs::path filename,
                                          const xtal::BasicStructure &prim,
                                          double tol);

void write_chemical_reference(const ChemicalReference &chem_ref,
                              fs::path filename);

}  // namespace CASM

#endif
