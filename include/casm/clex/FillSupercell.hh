#ifndef CASM_FillSupercell
#define CASM_FillSupercell

#include <memory>
#include <vector>

#include "casm/global/definitions.hh"

namespace CASM {

namespace xtal {
class Lattice;
}
class Configuration;
class Supercell;
class SymOp;

/// Return true if a configuration (in any equivalent orientation) can be tiled
/// into a supercell
///
/// Note:
/// - This will return true if, for any SymOp in the prim factor group,
/// `apply(symop, sub_configuration_lattice)`
///   can be used used to fill the Supercell. Otherwise it will return false.
bool is_valid_sub_configuration(xtal::Lattice const &sub_configuration_lattice,
                                Supercell const &supercell);

/// Create a super configuration by tiling the motif Configuration into the
/// supercell.
///
/// Note:
/// - This overload finds the first SymOp in the prim factor group such that
/// apply(symop, motif)
///   can be used to fill the Supercell. If none can be found, it throws.
Configuration fill_supercell(
    Configuration const &motif,
    std::shared_ptr<Supercell const> const &shared_supercell);

/// Create a super configuration by tiling the motif Configuration into the
/// supercell.
///
/// Note:
/// - Prefer to use `std::shared_ptr<Supercell const>` overload if possible
/// - This overload finds the first SymOp in the prim factor group such that
/// apply(symop, motif)
///   can be used to fill the Supercell. If none can be found, it throws.
Configuration fill_supercell(Configuration const &motif,
                             Supercell const &supercell);

/// Create a super configuration by tiling the apply(symop, motif) Configuration
/// into the supercell
///
/// Note:
/// - Will throw if apply(symop, motif) cannot be tiled into the supercell
Configuration fill_supercell(
    SymOp const &symop, Configuration const &motif,
    std::shared_ptr<Supercell const> const &shared_supercell);

/// Create a super configuration by tiling the apply(symop, motif) Configuration
/// into the supercell
///
/// Note:
/// - Prefer to use `std::shared_ptr<Supercell const>` overload if possible
/// - Will throw if apply(symop, motif) cannot be tiled into the supercell
Configuration fill_supercell(SymOp const &symop, Configuration const &motif,
                             Supercell const &supercell);

/// Make all super configurations that fill a particular supercell
template <typename ConfigOutputIterator>
ConfigOutputIterator make_all_super_configurations(
    Configuration const &configuration,
    std::shared_ptr<Supercell const> shared_supercell,
    ConfigOutputIterator result);

/// Functor to fill a Supercell with a tiling of a motif Configuration
///
/// Note:
/// - In general, the motif Configuration must be transformed by a prim factor
/// group operation
///   so that the motif's smaller supercell can be tiled into the larger
///   supercell.
/// - Construct with the Supercell to be filled and either the Lattice of the
/// motif configuration,
///   or the SymOp used to transform the motif before filling the larger
///   supercell.
/// - Using this functor directly will be more efficient if transforming many
/// motif Configuration
///   with the same Supercell into super configurations with the same Supercell
///   because will preserve some temporary data.
class FillSupercell {
 public:
  /// Constructor, so that the generated super configuration uses
  /// `std::shared_ptr<Supercell>`
  ///
  /// \param _shared_supercell Supercell to be filled
  FillSupercell(std::shared_ptr<Supercell const> const &_shared_supercell);

  /// Constructor, so that the generated super configuration uses `Supercell
  /// const &`
  ///
  /// \param _scel Supercell to be filled
  FillSupercell(Supercell const &_supercell);

  /// Create a super configuration by tiling the motif Configuration into the
  /// supercell.
  ///
  /// Note:
  /// - This overload finds the first SymOp in the prim factor group such that
  ///   apply(symop, motif) can be used to fill the Supercell. If none can be
  ///   found, it throws.
  /// - If it is called repeatedly on motif configurations that share the same
  /// Supercell it uses
  ///   the same SymOp.
  Configuration operator()(Configuration const &motif) const;

  /// Create a super configuration by tiling the apply(symop, motif)
  /// Configuration into the supercell
  Configuration operator()(SymOp const &_symop,
                           Configuration const &motif) const;

  /// Find first SymOp in the prim factor group such that apply(symop, motif)
  /// can be used to fill the Supercell. Returns nullptr if none can tile.
  const SymOp *find_symop(xtal::Lattice const &_motif_lattice) const;

  /// Returns the SymOp used by `operator()`
  const SymOp &symop() const { return *m_symop_ptr; }

 private:
  void _init(Supercell const &_motif_scel) const;

  std::shared_ptr<Supercell const> m_shared_supercell;
  Supercell const *m_supercell_ptr;

  // These mutable members are used by `operator()`. They are mutable members
  // instead of temporary variables so they can be reused if FillSupercell is
  // called repeatedly on motif Configuration from the same Supercell.

  mutable SymOp const *m_symop_ptr;
  mutable Supercell const *m_motif_supercell;
  mutable std::vector<std::vector<Index> > m_index_table;
};
}  // namespace CASM

#endif
