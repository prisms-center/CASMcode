#include "casm/clusterography/io/stream/ClusterInvariants_stream_io.hh"

#include <iomanip>

#include "casm/clusterography/ClusterInvariants.hh"

namespace CASM {

/// \brief Print ClusterInvariants
std::ostream &operator<<(std::ostream &sout,
                         ClusterInvariants const &invariants) {
  if (invariants.size() <= 1) {
    sout << "  #Points: " << invariants.size();
  } else {
    sout << "  #Points: " << invariants.size();
    sout << "  Site Distances: {";
    for (int i = 0; i < invariants.displacement().size(); i++) {
      if (i != 0) {
        sout << ", ";
      }
      sout << std::setprecision(5) << invariants.displacement()[i];
    }
    sout << "}";
  }
  return sout;
}

/// \brief Print WithinScelClusterInvariants
std::ostream &operator<<(std::ostream &sout,
                         WithinScelClusterInvariants const &invariants) {
  if (invariants.size() <= 1) {
    sout << "  #Points: " << invariants.size();
  } else {
    sout << "  #Points: " << invariants.size();
    sout << "  Site Distances: {";
    for (int i = 0; i < invariants.displacement().size(); i++) {
      if (i != 0) {
        sout << ", ";
      }
      sout << std::setprecision(5) << invariants.displacement()[i];
    }
    sout << "}";
  }
  return sout;
}

}  // namespace CASM
