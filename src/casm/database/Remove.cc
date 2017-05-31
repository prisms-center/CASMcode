#include "casm/database/Remove_impl.hh"

#include "casm/database/ConfigTypeTraits.hh"

// explicit template instantiations
#define INST_Remove(r, data, type) \
template class Remove<type>; \
template class RemoveT<type>; \

namespace CASM {
  namespace DB {

    BOOST_PP_SEQ_FOR_EACH(INST_Remove, _, CASM_DB_CONFIG_TYPES)

  }
}
