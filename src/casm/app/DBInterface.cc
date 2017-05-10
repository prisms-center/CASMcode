#include "casm/app/DBInterface_impl.hh"
#include "casm/app/select.hh"
#include "casm/app/query.hh"

namespace CASM {
  namespace DB {
    template InterfaceData<Configuration>::InterfaceData(const SelectCommand &);
    template InterfaceData<Configuration>::InterfaceData(const QueryCommand &);
  }
}

