#ifndef CASM_AppIO_impl
#define CASM_AppIO_impl

#include "casm/app/AppIO.hh"
#include "casm/casm_io/SafeOfstream.hh"
#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/clusterography/IntegralCluster.hh"
#include "casm/clusterography/io/json/IntegralCluster_json_io.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/database/Selection_impl.hh"

namespace CASM {

template <typename IterType>
void CompositionAxes::insert_enumerated(IterType begin, IterType end) {
  int i = 0;
  for (; begin != end; ++begin, ++i) {
    while (all_axes.count(std::to_string(i))) ++i;
    all_axes[std::to_string(i)] = *begin;
    enumerated.insert(std::to_string(i));
  }
}

}  // namespace CASM

#endif
