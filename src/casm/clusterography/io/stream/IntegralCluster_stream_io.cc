#include "casm/app/AppIO.hh"
#include "casm/clusterography/io/stream/IntegralCluster_stream_io.hh"
#include "casm/clusterography/IntegralCluster.hh"

namespace CASM {

  /// \brief Print IntegralCluster to stream, using default Printer<IntegralCluster>
  std::ostream &operator<<(std::ostream &sout, const IntegralCluster &clust) {
    OrbitPrinterOptions opt;
    opt.coord_type = INTEGRAL;
    SitesPrinter printer {opt};
    Log log(sout);
    printer.print(clust, log);
    return sout;
  }

}
