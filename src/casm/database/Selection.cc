#include "casm/clex/PrimClex.hh"
#include "casm/database/Selection_impl.hh"
#include "casm/database/DatabaseTypes.hh"

// explicit template instantiations
#define INST_Selection(r, data, type) \
template class Selection<type>; \

namespace CASM {
  namespace DB {
    BOOST_PP_SEQ_FOR_EACH(INST_Selection, _, CASM_DB_TYPES)
  }
}
