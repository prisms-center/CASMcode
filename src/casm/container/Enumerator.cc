#include "casm/container/Enumerator_impl.hh"

namespace CASM {

  namespace {
    typedef std::back_insert_iterator<std::vector<std::shared_ptr<RuntimeLibrary> > > runtimelib_it_type;
    typedef std::insert_iterator<EnumeratorMap> enum_it_type;
  }

  template std::pair<enum_it_type, runtimelib_it_type> load_enumerator_plugins(
    const PrimClex &primclex,
    enum_it_type enum_it,
    runtimelib_it_type lib_it);

}
