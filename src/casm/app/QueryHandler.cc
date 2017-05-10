#include "casm/app/QueryHandler_impl.hh"
#include "casm/database/DatabaseTypeTraits.hh"

namespace CASM {
  namespace {
    typedef std::insert_iterator<std::map<std::string, std::shared_ptr<RuntimeLibrary> > > runtimelib_it_type;
    \
  }

  // explicit template instantiations
#define INST_QueryHandler(r, data, type) \
template class QueryHandler<type>; \
template std::pair<std::insert_iterator<DataFormatterDictionary<type> >, runtimelib_it_type> load_query_plugins( \
  const ProjectSettings &set, \
  std::insert_iterator<DataFormatterDictionary<type> > dict_it, \
  runtimelib_it_type lib_it);

  BOOST_PP_SEQ_FOR_EACH(INST_QueryHandler, _, CASM_DB_TYPES)

}

