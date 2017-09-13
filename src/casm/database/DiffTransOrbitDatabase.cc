#include "casm/database/DiffTransOrbitDatabase.hh"
#include "casm/kinetics/DiffusionTransformation_impl.hh"

namespace CASM {
  namespace DB {

    /// Find PrimPeriodicDiffTransOrbit in database by comparing prototype
    typename Database<PrimPeriodicDiffTransOrbit>::iterator
    Database<PrimPeriodicDiffTransOrbit>::search(const PrimPeriodicDiffTransOrbit &orbit) const {
      return std::find(begin(), end(), orbit);
    }

    /// Find DiffusionTransformation in database by comparing to orbit prototypes
    typename Database<PrimPeriodicDiffTransOrbit>::iterator
    Database<PrimPeriodicDiffTransOrbit>::search(const Kinetics::DiffusionTransformation &diff_trans) const {

      const auto &g = prim().factor_group();
      PrimPeriodicDiffTransSymCompare sym_compare(crystallography_tol());
      return diff_trans.find_sym_equivalent(
               prototype_iterator(begin()),
               prototype_iterator(end()),
               g,
               sym_compare).base();

    }

  }
}
