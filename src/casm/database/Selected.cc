#include "casm/clex/PrimClex.hh"
#include "casm/database/Selected_impl.hh"

// explicit template instantiations
#define INST_Selected(r, data, type) \
template class Selected<type>; \

namespace CASM {
  namespace DB {
    BOOST_PP_SEQ_FOR_EACH(INST_Selected, _, CASM_DB_TYPES)
  }
}
