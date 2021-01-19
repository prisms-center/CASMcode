#include "casm/app/DBInterface_impl.hh"
#include "casm/app/query.hh"
#include "casm/app/select.hh"
#include "casm/clex/Supercell.hh"

namespace CASM {
namespace DB {
template InterfaceData<Configuration>::InterfaceData(const SelectCommand &);
template InterfaceData<Configuration>::InterfaceData(const QueryCommand &);
template InterfaceData<Supercell>::InterfaceData(const SelectCommand &);
template InterfaceData<Supercell>::InterfaceData(const QueryCommand &);
}  // namespace DB
}  // namespace CASM
