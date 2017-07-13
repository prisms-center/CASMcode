#include "casm/app/DBInterface_impl.hh"
#include "casm/app/select.hh"
#include "casm/app/query.hh"
#include "casm/clex/Supercell.hh"
#include "casm/kinetics/DiffusionTransformation.hh"
#include "casm/kinetics/DiffTransConfiguration.hh"

namespace CASM {
  namespace DB {
    template InterfaceData<Configuration>::InterfaceData(const SelectCommand &);
    template InterfaceData<Configuration>::InterfaceData(const QueryCommand &);
    template InterfaceData<Supercell>::InterfaceData(const SelectCommand &);
    template InterfaceData<Supercell>::InterfaceData(const QueryCommand &);
    template InterfaceData<PrimPeriodicDiffTransOrbit>::InterfaceData(const SelectCommand &);
    template InterfaceData<PrimPeriodicDiffTransOrbit>::InterfaceData(const QueryCommand &);
    template InterfaceData<Kinetics::DiffTransConfiguration>::InterfaceData(const SelectCommand &);
    template InterfaceData<Kinetics::DiffTransConfiguration>::InterfaceData(const QueryCommand &);
  }
}

