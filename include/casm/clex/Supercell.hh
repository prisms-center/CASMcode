#ifndef CASM_Supercell
#define CASM_Supercell

#include "casm/clex/HasCanonicalForm.hh"
#include "casm/clex/HasPrimClex.hh"
#include "casm/clex/SupercellTraits.hh"
#include "casm/crystallography/Lattice.hh"
#include "casm/database/Named.hh"
#include "casm/misc/Comparisons.hh"
#include "casm/misc/cloneable_ptr.hh"
#include "casm/symmetry/SupercellSymInfo.hh"
#include "casm/symmetry/SymGroup.hh"
#include "casm/symmetry/SymGroupRepID.hh"

namespace CASM {

namespace xtal {
class Site;
class Coordinate;
class UnitCellCoord;
class IntegralCoordinateWithin_f;
}  // namespace xtal
using xtal::Coordinate;
using xtal::Site;
using xtal::UnitCellCoord;

template <typename T, typename U>
class ConfigIterator;
class PermuteIterator;
class PrimClex;
class Clexulator;
class PrimNeighborList;
class SuperNeighborList;
class Structure;

namespace DB {
template <typename T>
class DatabaseIterator;
}

/** \defgroup Supercell
 *  \ingroup Clex
 *  \brief Represents a supercell of the primitive parent crystal structure
 *  @{
 */

/// \brief Represents a supercell of the primitive parent crystal structure
///
class Supercell
    : public DB::Named<
          Comparisons<SupercellCanonicalForm<CRTPBase<Supercell>>>> {
 public:
  using permute_const_iterator = SupercellSymInfo::permute_const_iterator;

  // **** Constructors ****

  Supercell(const Supercell &RHS);

  Supercell(std::shared_ptr<Structure const> const &_shared_prim,
            const Lattice &superlattice);
  Supercell(std::shared_ptr<Structure const> const &_shared_prim,
            Eigen::Matrix3l const &superlattice_matrix);

  Supercell(const PrimClex *_prim, const Lattice &superlattice);
  Supercell(const PrimClex *_prim,
            const Eigen::Ref<const Eigen::Matrix3l> &superlattice_matrix);

  ~Supercell();

  // **** Coordinates ****

  /// \brief Return the sublattice index for a linear index
  Index sublat(Index linear_index) const;

  /// \brief Given a Coordinate and tolerance, return linear index into
  /// Configuration
  Index linear_index(const Coordinate &coord, double tol = TOL) const;

  /// \brief Return the linear index corresponding to integral coordinates
  Index linear_index(const UnitCellCoord &bijk) const;

  /// \brief Return the linear index corresponding to integral coordinates
  Coordinate coord(Index linear_index) const;

  /// \brief Return the integral coordinates corresponding to a linear index
  UnitCellCoord uccoord(Index linear_index) const;

  /// \brief returns maximum allowed occupation bitstring -- used for
  /// initializing enumeration counters
  Eigen::VectorXi max_allowed_occupation() const;

  // **** Accessors ****

  const Structure &prim() const;

  std::shared_ptr<Structure const> const &shared_prim() const;

  double crystallography_tol() const;

  /// Use while transitioning Supercell to no longer need a `PrimClex const *`
  bool has_primclex() const;

  /// Use while transitioning Supercell to no longer need a `PrimClex const *`
  void set_primclex(PrimClex const *primclex_ptr) const;

  /// Use while transitioning Supercell to no longer need a `PrimClex const *`
  const PrimClex &primclex() const;

  /// Return number of primitive cells that fit inside of *this
  Index volume() const;

  Index basis_size() const;

  Index num_sites() const;

  SymGroupRep const &permutation_symrep() const;

  Eigen::Matrix3l transf_mat() const;

  /// \brief The super lattice
  const Lattice &lattice() const;

  /// \brief Returns the SuperNeighborList
  ///
  /// Requires that the prim_nlist has been set by one of:
  /// - constructing Supercell with a PrimClex const *
  /// - using set_primclex to set a PrimClex const *
  ///
  /// At each access, the underlying PrimNeighborList will be checked and if it
  /// has been expanded then the SuperNeighborList will be extended also.
  /// References obtained from this function will be out of date if the
  /// underlying PrimNeighborList has been expanded, so it is prudent to only
  /// access the SuperNeighborList for immediate use.
  const SuperNeighborList &nlist() const;

  // Factor group of this supercell
  const SymGroup &factor_group() const;

  // SymInfo object of this supercell
  const SupercellSymInfo &sym_info() const;

  bool operator<(const Supercell &B) const;

  /// \brief Insert the canonical form of this into the database
  std::pair<DB::DatabaseIterator<Supercell>, bool> insert() const;

 private:
  friend Comparisons<SupercellCanonicalForm<CRTPBase<Supercell>>>;
  friend DB::Named<Comparisons<SupercellCanonicalForm<CRTPBase<Supercell>>>>;

  bool eq_impl(const Supercell &B) const;

  // **** Generating functions ****

  std::string generate_name_impl() const;

  // Note:
  // - Prefer not to access PrimClex via Supercell. In future, PrimClex access
  // via Supercell will
  //   be removed completely.
  // - Until this is removed, it may be nullptr, in which case, some features
  // will throw. Only
  //   access via this->primclex() so that an error will be thrown if m_primclex
  //   is nullptr.
  // - Mutable as a temporary workaround
  mutable PrimClex const *m_primclex;

  std::shared_ptr<Structure const> m_shared_prim;

  SupercellSymInfo m_sym_info;

  /// SuperNeighborList, mutable for lazy construction
  mutable notstd::cloneable_ptr<SuperNeighborList> m_nlist;

  /// Store size of PrimNeighborList at time of construction of
  /// SuperNeighborList to enable checking if SuperNeighborList should be
  /// re-constructed
  mutable Index m_nlist_size_at_construction;
};

/// Make the supercell name from a Superlattice
std::string make_supercell_name(Structure const &prim,
                                xtal::Superlattice const &superlattice);

/// Make the canonical supercell name from a Superlattice
std::string make_canonical_supercell_name(
    Structure const &prim, xtal::Superlattice const &superlattice);

/// Construct a Superlattice from the supercell name
xtal::Superlattice make_superlattice_from_supercell_name(
    Structure const &prim, std::string supercell_name);

/// Apply symmetry operation to Supercell
Supercell &apply(const SymOp &op, Supercell &scel);

/// Copy and apply symmetry operation to Supercell
Supercell copy_apply(const SymOp &op, const Supercell &scel);

// --- The following are deprecated ----

const Supercell &make_supercell(const PrimClex &primclex, std::string name);

std::shared_ptr<Supercell> make_shared_supercell(const PrimClex &primclex,
                                                 std::string name);

Eigen::Matrix3l transf_mat(const Lattice &prim_lat, const Lattice &super_lat);

std::string generate_name(const Eigen::Matrix3l &transf_mat);

std::string scelname(const Structure &prim, const Lattice &superlat);

std::string canonical_scelname(const Structure &prim, const Lattice &superlat);

namespace xtal {
IntegralCoordinateWithin_f make_bring_within_f(const Supercell &scel);
}

/** @} */
}  // namespace CASM
#endif
