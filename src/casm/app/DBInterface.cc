#include "casm/app/DBInterface_impl.hh"
#include "casm/app/select.hh"
#include "casm/app/query.hh"

namespace CASM {
  namespace DB {
    template InterfaceData<Configuration>::InterfaceData(const SelectCommand &);
    template InterfaceData<Configuration>::InterfaceData(const QueryCommand &);
    template InterfaceData<Supercell>::InterfaceData(const SelectCommand &);
    template InterfaceData<Supercell>::InterfaceData(const QueryCommand &);
    template InterfaceData<PrimPeriodicDiffTransOrbit>::InterfaceData(const SelectCommand &);
    template InterfaceData<PrimPeriodicDiffTransOrbit>::InterfaceData(const QueryCommand &);
  }
}

