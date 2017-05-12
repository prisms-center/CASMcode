#include "casm/monte_carlo/MonteDriver_impl.hh"
#include "casm/monte_carlo/canonical/Canonical.hh"
#include "casm/monte_carlo/canonical/CanonicalIO.hh"
#include "casm/monte_carlo/grand_canonical/GrandCanonical.hh"
#include "casm/monte_carlo/grand_canonical/GrandCanonicalIO.hh"

namespace CASM {
  namespace Monte {
    template class MonteDriver<Canonical>;
    template class MonteDriver<GrandCanonical>;

  }
}
