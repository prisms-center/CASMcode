#include "casm/clex/FillSupercell.hh"

#include "casm/clex/Configuration.hh"
#include "casm/clex/Supercell.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/crystallography/SymTools.hh"
#include "casm/crystallography/UnitCellCoord.hh"
#include "casm/symmetry/SymOp.hh"
#include "casm/symmetry/SymTools.hh"

namespace CASM {

/// Return true if a configuration (in any equivalent orientation) can be tiled
/// into a supercell
///
/// Note:
/// - This will return true if, any SymOp in the prim factor group,
/// `apply(symop, configuration)`
///   can be used used to fill the Supercell. Otherwise it will return false.
bool is_valid_sub_configuration(Lattice const &sub_configuration_lattice,
                                Supercell const &supercell) {
  FillSupercell f{supercell};
  return f.find_symop(sub_configuration_lattice) != nullptr;
}

/// Create a super configuration by tiling the motif Configuration into the
/// supercell.
///
/// Note:
/// - This overload finds the first SymOp in the prim factor group such that
/// apply(symop, motif)
///   can be used to fill the Supercell. If none can be found, it throws.
Configuration fill_supercell(
    Configuration const &motif,
    std::shared_ptr<Supercell const> const &shared_supercell) {
  FillSupercell f{shared_supercell};
  return f(motif);
}

/// Create a super configuration by tiling the motif Configuration into the
/// supercell.
///
/// Note:
/// - Prefer to use `std::shared_ptr<Supercell const>` overload if possible
/// - This overload finds the first SymOp in the prim factor group such that
/// apply(symop, motif)
///   can be used to fill the Supercell. If none can be found, it throws.
Configuration fill_supercell(Configuration const &motif,
                             Supercell const &supercell) {
  FillSupercell f{supercell};
  return f(motif);
}

/// Create a super configuration by tiling the apply(symop, motif) Configuration
/// into the supercell
///
/// Note:
/// - Will throw if apply(symop, motif) cannot be tiled into the supercell
Configuration fill_supercell(
    SymOp const &symop, Configuration const &motif,
    std::shared_ptr<Supercell const> const &shared_supercell) {
  FillSupercell f{shared_supercell};
  return f(symop, motif);
}

/// Create a super configuration by tiling the apply(symop, motif) Configuration
/// into the supercell
///
/// Note:
/// - Prefer to use `std::shared_ptr<Supercell const>` overload if possible
/// - Will throw if apply(symop, motif) cannot be tiled into the supercell
Configuration fill_supercell(SymOp const &symop, Configuration const &motif,
                             Supercell const &supercell) {
  FillSupercell f{supercell};
  return f(symop, motif);
}

/// \brief Constructor
///
/// \param _scel Supercell to be filled
/// \param _op SymOp that transforms the input motif before tiling into the
///        Supercell that is filled
FillSupercell::FillSupercell(
    std::shared_ptr<Supercell const> const &_shared_supercell)
    : m_shared_supercell(_shared_supercell),
      m_supercell_ptr(m_shared_supercell.get()),
      m_symop_ptr(nullptr),
      m_motif_supercell(nullptr) {}

/// \brief Constructor
///
/// \param _scel Supercell to be filled
/// \param _op SymOp that transforms the input motif before tiling into the
///        Supercell that is filled
FillSupercell::FillSupercell(Supercell const &_supercell)
    : m_supercell_ptr(&_supercell),
      m_symop_ptr(nullptr),
      m_motif_supercell(nullptr) {}

Configuration FillSupercell::operator()(Configuration const &motif) const {
  if (&motif.supercell() != m_motif_supercell) {
    m_symop_ptr = find_symop(motif.supercell().lattice());
    _init(motif.supercell());
  }
  return (*this)(*m_symop_ptr, motif);
}

Configuration FillSupercell::operator()(SymOp const &symop,
                                        Configuration const &motif) const {
  if (&motif.supercell() != m_motif_supercell || &symop != m_symop_ptr) {
    if (&symop != m_symop_ptr) {
      Lattice const &motif_lattice = motif.supercell().lattice();
      if (!is_superlattice(m_supercell_ptr->lattice(),
                           sym::copy_apply(symop, motif_lattice),
                           motif_lattice.tol())
               .first) {
        throw std::runtime_error(
            "Error in 'FillSupercell: super lattice != sym::copy_apply(symop, "
            "motif lattice)");
      }
      m_symop_ptr = &symop;
    }
    _init(motif.supercell());
  }

  std::unique_ptr<Configuration> result;
  if (m_shared_supercell) {
    result = notstd::make_unique<Configuration>(m_shared_supercell);
  } else {
    result = notstd::make_unique<Configuration>(*m_supercell_ptr);
  }

  // We reorien the starting configuration (motif) by symop in two steps.
  // In the first step, we transform the DoFs of motif by symop, without
  // permuting their site indices Then, we utilize m_index_table to decorate the
  // new cell with these transformed DoFs. m_index_table is built in
  // FillSupercel::_init() by explicitly transforming each UnitCellCoordinate of
  // the starting supercell (*m_motif_supercell) and finding it in the resultant
  // supercell (*m_supercell_ptr) This ensures that sublattice permutation and
  // microscopic translations are handled appropriately
  ConfigDoF trans_motif(motif.configdof());
  trans_motif.apply_sym_no_permute(symop);

  // copy transformed dof, as many times as necessary to fill the supercell
  for (auto const &dof : trans_motif.global_dofs())
    result->configdof().set_global_dof(dof.first, dof.second.values());

  for (Index s = 0; s < m_index_table.size(); ++s) {
    for (Index i = 0; i < m_index_table[s].size(); ++i) {
      Index scel_s = m_index_table[s][i];

      if (trans_motif.has_occupation()) {
        result->configdof().occ(scel_s) = trans_motif.occ(s);
      }
    }
  }

  for (auto const &dof : trans_motif.local_dofs()) {
    LocalContinuousConfigDoFValues &res_ref(
        result->configdof().local_dof(dof.first));
    for (Index s = 0; s < m_index_table.size(); ++s) {
      for (Index i = 0; i < m_index_table[s].size(); ++i) {
        Index scel_s = m_index_table[s][i];
        res_ref.site_value(scel_s) = dof.second.site_value(s);
      }
    }
  }
  return *result;
}

/// \brief Find first SymOp in the prim factor group such that apply(op, motif)
///        can be used to fill the Supercell
SymOp const *FillSupercell::find_symop(Lattice const &motif_lattice) const {
  const Lattice &scel_lattice = m_supercell_ptr->lattice();
  auto begin = m_supercell_ptr->prim().factor_group().begin();
  auto end = m_supercell_ptr->prim().factor_group().end();
  auto res = xtal::is_equivalent_superlattice(scel_lattice, motif_lattice,
                                              begin, end, motif_lattice.tol());
  if (res.first == end) {
    return nullptr;
  }
  return &(*res.first);
}

void FillSupercell::_init(Supercell const &_motif_supercell) const {
  m_motif_supercell = &_motif_supercell;

  // ------- site dof ----------
  Lattice oriented_motif_lat =
      sym::copy_apply(*m_symop_ptr, m_motif_supercell->lattice());

  auto oriended_motif_lattice_points = xtal::make_lattice_points(
      oriented_motif_lat, m_supercell_ptr->lattice(), TOL);

  const Structure &prim = m_supercell_ptr->prim();
  m_index_table.resize(m_motif_supercell->num_sites());

  // for each site in motif
  for (Index s = 0; s < m_motif_supercell->num_sites(); s++) {
    // apply symmetry to re-orient and find unit cell coord
    UnitCellCoord oriented_uccoord =
        sym::copy_apply(*m_symop_ptr, m_motif_supercell->uccoord(s),
                        prim.lattice(), prim.basis_permutation_symrep_ID());

    // for each unit cell of the oriented motif in the supercell, copy the
    // occupation
    for (const UnitCell &oriented_motif_uc : oriended_motif_lattice_points) {
      UnitCell oriented_motif_uc_relative_to_prim =
          oriented_motif_uc.reset_tiling_unit(oriented_motif_lat,
                                              prim.lattice());

      Index prim_motif_tile_ind =
          m_supercell_ptr->sym_info().unitcell_index_converter()(
              oriented_motif_uc_relative_to_prim);

      UnitCellCoord mc_uccoord(
          oriented_uccoord.sublattice(),
          m_supercell_ptr->sym_info().unitcell_index_converter()(
              prim_motif_tile_ind) +
              oriented_uccoord.unitcell());

      m_index_table[s].push_back(m_supercell_ptr->linear_index(mc_uccoord));
    }
  }
}

}  // namespace CASM
