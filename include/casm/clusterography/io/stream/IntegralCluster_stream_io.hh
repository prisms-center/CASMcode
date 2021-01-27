#ifndef CASM_IntegralCluster_stream_io
#define CASM_IntegralCluster_stream_io

#include <iostream>

namespace CASM {

class IntegralCluster;

/// \brief Print IntegralCluster to stream, using default
/// Printer<IntegralCluster>
std::ostream &operator<<(std::ostream &sout, const IntegralCluster &clust);

}  // namespace CASM

#endif
