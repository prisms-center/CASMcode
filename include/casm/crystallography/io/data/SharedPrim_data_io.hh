#ifndef CASM_xtal_SharedPrim_data_io
#define CASM_xtal_SharedPrim_data_io

#include <memory>

#include "casm/casm_io/dataformatter/DataFormatterDecl.hh"

namespace CASM {

class Structure;

template <>
DataFormatterDictionary<std::shared_ptr<const Structure>>
make_attribute_dictionary<std::shared_ptr<const Structure>>();

}  // namespace CASM

#endif
