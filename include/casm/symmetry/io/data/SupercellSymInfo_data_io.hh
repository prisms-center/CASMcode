#ifndef CASM_xtal_SupercellSymInfo_data_io
#define CASM_xtal_SupercellSymInfo_data_io

#include "casm/casm_io/dataformatter/DataFormatterDecl.hh"

namespace CASM {

class SupercellSymInfo;

template <>
DataFormatterDictionary<SupercellSymInfo>
make_attribute_dictionary<SupercellSymInfo>();

}  // namespace CASM

#endif
