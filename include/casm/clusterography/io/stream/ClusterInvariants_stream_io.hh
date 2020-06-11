#ifndef CASM_ClusterInvariants_stream_io
#define CASM_ClusterInvariants_stream_io

#include <iostream>

namespace CASM {

  class ClusterInvariants;
  class WithinScelClusterInvariants;

  /// \brief Print ClusterInvariants
  std::ostream &operator<<(std::ostream &sout, ClusterInvariants const &invariants);

  /// \brief Print WithinScelClusterInvariants
  std::ostream &operator<<(std::ostream &sout, WithinScelClusterInvariants const &invariants);

}

#endif
