#ifndef CASM_LatticeCanonicalForm
#define CASM_LatticeCanonicalForm

#include <vector>

#include "casm/crystallography/Adapter.hh"
#include "casm/crystallography/Lattice.hh"
#include "casm/global/definitions.hh"

namespace CASM {
namespace xtal {
struct SymOp;
typedef std::vector<SymOp> SymOpVector;
class Lattice;

namespace canonical {

/// True if lat_column_mat is approximately equal to the lat_column_mat of the
/// canonical_form, using the lattice point group to find the most canonical
/// form
bool check(const Lattice &lat);

/// True if lat_column_mat is approximately equal to the lat_column_mat of the
/// canonical_form, using the provided symmetry operations to find the most
/// canonical form
bool check(const Lattice &lat, SymOpVector const &g);

template <typename ExternSymOpVector>
bool check(const Lattice &lat, ExternSymOpVector const &g) {
  return check(lat, adapter::Adapter<SymOpVector, ExternSymOpVector>()(g));
}

//**************************************************************************************************************//

/**
 * Return a canonical Lattice by first converting the
 * given Lattice into the standard Niggli form,
 * followed by applying the point group of the Lattice
 * so that the one oriented in the most standard manner
 * is selected.
 */
Lattice equivalent(Lattice const &in_lat, SymOpVector const &point_grp,
                   double compare_tol);

template <typename ExternSymOpVector>
Lattice equivalent(Lattice const &in_lat, ExternSymOpVector const &point_grp,
                   double compare_tol) {
  return equivalent(
      in_lat, adapter::Adapter<SymOpVector, ExternSymOpVector>()(point_grp),
      compare_tol);
}

/// Canonical equivalent lattice, using the provided group
Lattice equivalent(const Lattice &lat, SymOpVector const &g);

template <typename ExternSymOpVector>
Lattice equivalent(Lattice const &in_lat, ExternSymOpVector const &g) {
  return equivalent(in_lat,
                    adapter::Adapter<SymOpVector, ExternSymOpVector>()(g));
}

/// Canonical equivalent lattice, using this lattice's point group
Lattice equivalent(const Lattice &lat);

/// Canonical equivalent lattice, using this lattice's point group and using
/// specified tolerance
Lattice equivalent(const Lattice &lat, double compare_tol);
//**************************************************************************************************************//

/// Uses provided group to find 'to_canonical' SymOp
///
/// - Returns first SymOp for which canonical_form.is_equivalent(apply(op,
/// *this))
/// - Note that that copy_apply(this->to_canonical(), *this).is_canonical()
///   may be false because they may be equivalent, but without identical
///   lat_column_mat().
Index operation_index(const Lattice &lat, SymOpVector const &g);

template <typename ExternSymOpVector>
Index operation_index(const Lattice &lat, ExternSymOpVector const &g) {
  return operation_index(lat,
                         adapter::Adapter<SymOpVector, ExternSymOpVector>()(g));
}

/// Return the index of the operation that makes the lattice canonical
Index operation_index(Lattice const &in_lat, SymOpVector const &point_grp,
                      double compare_tol);
}  // namespace canonical

}  // namespace xtal
}  // namespace CASM

#endif
