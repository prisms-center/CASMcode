#include "casm/app/QueryHandler_impl.hh"
#include "casm/database/DatabaseTypeTraits.hh"

namespace CASM {

  // explicit template instantiations
#define INST_QueryHandler(r, data, type) template class QueryHandler<type>;
  BOOST_PP_SEQ_FOR_EACH(INST_QueryHandler, _, CASM_DB_TYPES)

}

