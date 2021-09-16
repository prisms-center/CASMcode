#include "casm/crystallography/SimpleStructureTools.hh"

#include <stdexcept>

#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/IntegralCoordinateWithin.hh"
#include "casm/crystallography/SimpleStructure.hh"
#include "casm/crystallography/Site.hh"
#include "casm/crystallography/UnitCellCoord.hh"
#include "casm/external/Eigen/Core"

namespace CASM {
namespace xtal {
namespace Local {

static SimpleStructure::Info _replicate(SimpleStructure::Info const &_info,
                                        Index mult) {
  SimpleStructure::Info result;
  result.resize(_info.size() * mult);

  for (Index i = 0; i < _info.size(); ++i)
    result.coords.block(0, i * mult, 3, mult) =
        _info.cart_coord(i).replicate(1, mult);

  for (auto const &p : _info.properties) {
    result.properties.emplace(
        p.first, Eigen::MatrixXd(p.second.rows(), mult * p.second.cols()));
    for (Index i = 0; i < p.second.cols(); ++i)
      result.coords.block(0, i * mult, p.second.rows(), mult) =
          p.second.col(i).replicate(1, mult);
  }

  Index l = 0;
  for (Index b = 0; b < _info.size(); ++b) {
    for (Index g = 0; g < mult; ++g) {
      result.names[l++] = _info.names[b];
    }
  }
  return result;
}
}  // namespace Local

//***************************************************************************

SimpleStructure make_superstructure(Eigen::Ref<const Eigen::Matrix3i> const &_T,
                                    SimpleStructure const &_sstruc) {
  SimpleStructure superstructure;
  superstructure.lat_column_mat = _sstruc.lat_column_mat * _T.cast<double>();
  superstructure.properties = _sstruc.properties;

  auto all_lattice_points = make_lattice_points(_T.cast<long>());

  Index Nvol = all_lattice_points.size();

  superstructure.mol_info = Local::_replicate(_sstruc.mol_info, Nvol);
  superstructure.atom_info = Local::_replicate(_sstruc.atom_info, Nvol);

  Index nm = _sstruc.mol_info.size();
  Index na = _sstruc.atom_info.size();

  for (Index g = 0; g < Nvol; ++g) {
    Eigen::Vector3d lattice_point_vector =
        _sstruc.lat_column_mat * all_lattice_points[g].cast<double>();

    for (Index m = 0; m < nm; ++m) {
      superstructure.mol_info.cart_coord(g + m * Nvol) += lattice_point_vector;
    }
    for (Index a = 0; a < na; ++a) {
      superstructure.atom_info.cart_coord(g + a * Nvol) += lattice_point_vector;
    }
  }

  return superstructure;
}

//***************************************************************************************************

/// Calculates the parent basis index of each site in a supercell that is
/// generated with the make_superstructure method
/// @param _T is the transformation matrix linking `_sstruc` and the supercell
/// @param _sstruc is a SimpleStructure
/// @return std::vector<Index> with the Index in the ith entry corresponding to
/// the index of site i in _sstruc
std::vector<Index> superstructure_basis_idx(
    Eigen::Ref<const Eigen::Matrix3i> const &_T,
    SimpleStructure const &_sstruc) {
  auto all_lattice_points = make_lattice_points(_T.cast<long>());
  Index Nvol = all_lattice_points.size();
  std::vector<Index> basis_idx(_sstruc.atom_info.size() * Nvol, -1);
  for (Index grid_idx = 0; grid_idx < Nvol; ++grid_idx) {
    for (Index atom_idx = 0; atom_idx < _sstruc.atom_info.size(); ++atom_idx)
      basis_idx[grid_idx + atom_idx * Nvol] = atom_idx;
  }
  return basis_idx;
}

//***************************************************************************

SimpleStructure make_simple_structure(BasicStructure const &_struc) {
  SimpleStructure result;
  result.lat_column_mat = _struc.lattice().lat_column_mat();

  result.mol_info.coords.resize(3, _struc.basis().size());
  result.mol_info.names.reserve(_struc.basis().size());
  Eigen::VectorXi _mol_occ;
  // For now, default to first occupant. In future we may decide
  // to force user to pass mol_occ explicitly
  _mol_occ.setZero(_struc.basis().size());
  for (Index b = 0; b < _struc.basis().size(); ++b) {
    result.mol_info.cart_coord(b) = _struc.basis()[b].const_cart();
    result.mol_info.names.push_back(
        _struc.basis()[b].occupant_dof()[_mol_occ[b]].name());
  }
  _atomize(result, _mol_occ, _struc);
  return result;
}

BasicStructure make_basic_structure(
    SimpleStructure const &_sstruc, std::vector<DoFKey> const &_all_dofs,
    SimpleStructure::SpeciesMode mode,
    std::vector<std::vector<Molecule>> _allowed_occupants) {
  std::map<DoFKey, DoFSet> global_dof;
  std::map<DoFKey, SiteDoFSet> local_dof;
  for (DoFKey const &dof : _all_dofs) {
    if (AnisoValTraits(dof).global()) {
      global_dof.emplace(dof, AnisoValTraits(dof));
    } else {
      local_dof.emplace(dof, AnisoValTraits(dof));
    }
  }

  auto const &info = _sstruc.info(mode);
  if (_allowed_occupants.empty()) _allowed_occupants.resize(info.size());
  for (Index i = 0; i < info.size(); ++i) {
    if (_allowed_occupants[i].empty()) {
      _allowed_occupants[i].push_back(Molecule::make_atom(info.names[i]));
    }
    if (_allowed_occupants[i].size() == 1) {
      std::map<std::string, SpeciesProperty> attr_map =
          _allowed_occupants[i][0].properties();
      for (auto it = attr_map.begin(); it != attr_map.end(); ++it) {
        if (local_dof.count(it->first)) {
          auto er_it = it++;
          attr_map.erase(er_it);
        }
      }

      for (auto const &prop : info.properties) {
        if (local_dof.count(prop.first)) continue;

        if (prop.first == "disp") continue;

        if (!almost_zero(prop.second.col(i)))
          attr_map.emplace(prop.first,
                           SpeciesProperty(prop.first, prop.second.col(i)));
      }
      _allowed_occupants[i][0].set_properties(attr_map);
    }
  }

  BasicStructure result(Lattice(_sstruc.lat_column_mat));
  result.set_global_dofs(global_dof);
  std::vector<Site> tbasis(info.size(), Site(result.lattice()));

  for (Index i = 0; i < info.size(); ++i) {
    tbasis[i].cart() = info.cart_coord(i);
    tbasis[i].set_allowed_occupants(std::move(_allowed_occupants[i]));
    tbasis[i].set_dofs(local_dof);
  }

  result.set_basis(tbasis);
  return result;
}

//***************************************************************************

/**
** Calculates the invariant shuffle modes in the primitive unit cell. They
*symmetry preserving distortions are found by applying the Reynolds operator to
*the basis of displacement vectors. The average of the resulting basis vectors
*is used to form an orthonormal basis.
*
* @param factor_group use make_factor_group(struc) to obtain this group
* @param permute_group use make_permute_group(struc,factor_group) to obtain this
*group
*
* @return the vector of shuffle modes that are invariant under symmetry. Each
* element has a size `3 x basis_size`.
*
*/
std::vector<Eigen::MatrixXd> generate_invariant_shuffle_modes(
    const std::vector<xtal::SymOp> &factor_group,
    const std::vector<Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic,
                                               Index>> &permute_group) {
  if (factor_group.size() != permute_group.size()) {
    throw std::runtime_error(
        "error, the size of the symmetry operations in "
        "generate_invariant_shuffle_modes do not match");
  }
  int struc_basis_size = permute_group[0].indices().size();
  // Generate a basis consisting of individual shuffles of each atom in the
  // structure.
  std::vector<Eigen::MatrixXd> displacement_basis;
  for (int basis_idx = 0; basis_idx < struc_basis_size; ++basis_idx) {
    for (int dir_idx = 0; dir_idx < 3; ++dir_idx) {
      Eigen::MatrixXd single_shuffle =
          Eigen::MatrixXd::Zero(3, struc_basis_size);
      single_shuffle(dir_idx, basis_idx) = 1.0;
      displacement_basis.push_back(single_shuffle);
    }
  }
  std::vector<Eigen::VectorXd> displacement_aggregate(
      displacement_basis.size(),
      Eigen::VectorXd::Zero(displacement_basis[0].cols() *
                            displacement_basis[0].rows()));

  for (int idx = 0; idx < factor_group.size(); ++idx) {
    for (int disp_basis_idx = 0; disp_basis_idx < displacement_basis.size();
         ++disp_basis_idx) {
      Eigen::MatrixXd transformed_disp = factor_group[idx].matrix *
                                         displacement_basis[disp_basis_idx] *
                                         permute_group[idx];
      displacement_aggregate[disp_basis_idx] +=
          Eigen::VectorXd(Eigen::Map<Eigen::VectorXd>(
              transformed_disp.data(),
              transformed_disp.cols() * transformed_disp.rows()));
    }
  }
  Eigen::MatrixXd sym_disp_basis = Eigen::MatrixXd::Zero(
      displacement_aggregate.size(), displacement_aggregate[0].size());
  for (int disp_basis_idx = 0; disp_basis_idx < displacement_basis.size();
       ++disp_basis_idx) {
    displacement_aggregate[disp_basis_idx] =
        displacement_aggregate[disp_basis_idx] / double(factor_group.size());
    sym_disp_basis.row(disp_basis_idx) = displacement_aggregate[disp_basis_idx];
  }

  // Perform a SVD of the resulting aggregate matrix to obtain the rank and
  // space of the symmetry invariant basis vectors
  sym_disp_basis.transposeInPlace();
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(
      sym_disp_basis, Eigen::ComputeThinU | Eigen::ComputeThinV);
  int matrix_rank = (svd.singularValues().array().abs() >= CASM::TOL).sum();
  Eigen::MatrixXd sym_preserving_mode_matrix =
      svd.matrixV()(Eigen::all, Eigen::seq(0, matrix_rank - 1));
  std::vector<Eigen::MatrixXd> sym_preserving_modes;

  for (int sym_mode_idx = 0; sym_mode_idx < matrix_rank; ++sym_mode_idx) {
    Eigen::Map<Eigen::MatrixXd> _tmp_mode(
        sym_preserving_mode_matrix.col(sym_mode_idx).data(), 3,
        struc_basis_size);
    sym_preserving_modes.push_back(_tmp_mode);
  }
  return sym_preserving_modes;
}

//***************************************************************************
void _atomize(SimpleStructure &_sstruc,
              Eigen::Ref<const Eigen::VectorXi> const &_mol_occ,
              BasicStructure const &_reference) {
  Index N_atoms(0);

  Index nb = _reference.basis().size();
  Index nv = _mol_occ.size() / nb;
  Index s = 0;
  for (Index b = 0; b < nb; ++b) {
    for (Index v = 0; v < nv; ++v, ++s) {
      N_atoms += _reference.basis()[b].occupant_dof()[_mol_occ[s]].size();
    }
  }
  _sstruc.atom_info.coords.resize(3, N_atoms);
  _sstruc.atom_info.names.resize(N_atoms);

  // s indexes site (i.e., molecule), a is index of atom within the entire
  // structure
  Index a = 0;
  s = 0;
  for (Index b = 0; b < nb; ++b) {
    for (Index v = 0; v < nv; ++v, ++s) {
      Molecule const &molref =
          _reference.basis()[b].occupant_dof()[_mol_occ[s]];

      // Initialize atom_info.properties for *molecule* properties
      for (auto const &property : _sstruc.mol_info.properties) {
        auto it = _sstruc.atom_info.properties.find(property.first);
        if (it == _sstruc.atom_info.properties.end()) {
          Index dim = AnisoValTraits(property.first).dim();
          _sstruc.atom_info.properties.emplace(
              property.first, Eigen::MatrixXd::Zero(dim, N_atoms));
        }
      }

      // ma is index of atom within individual molecule
      for (Index ma = 0; ma < molref.size(); ++ma, ++a) {
        // Record position of atom
        _sstruc.atom_info.cart_coord(a) =
            _sstruc.mol_info.cart_coord(s) + molref.atom(ma).cart();
        // Record name of atom
        _sstruc.atom_info.names[a] = molref.atom(ma).name();

        // Initialize atom_info.properties for *atom* properties
        for (auto const &attr : molref.atom(ma).properties()) {
          auto it = _sstruc.atom_info.properties.find(attr.first);
          if (it == _sstruc.atom_info.properties.end()) {
            it = _sstruc.atom_info.properties
                     .emplace(attr.first,
                              Eigen::MatrixXd::Zero(attr.second.traits().dim(),
                                                    N_atoms))
                     .first;
          }
          // Record properties of atom
          it->second.col(a) = attr.second.value();
        }

        // Split molecule properties into atom properties using appropriate
        // extensivity rules If an property is specified both at the atom and
        // molecule levels then the two are added
        for (auto const &property : _sstruc.mol_info.properties) {
          auto it = _sstruc.atom_info.properties.find(property.first);
          if (AnisoValTraits(property.first).extensive()) {
            it->second.col(a) += property.second.col(s) / double(molref.size());
          } else {
            it->second.col(a) += property.second.col(s);
          }
        }
      }
    }
  }
}

//***************************************************************************

std::vector<std::set<Index>> mol_site_compatibility(
    SimpleStructure const &sstruc, BasicStructure const &_prim) {
  std::vector<std::set<Index>> result;
  result.reserve(sstruc.mol_info.names.size());
  for (std::string const &sp : sstruc.mol_info.names) {
    result.push_back({});
    for (Index b = 0; b < _prim.basis().size(); ++b) {
      if (_prim.basis()[b].contains(sp)) {
        result.back().insert(b);
      }
    }
  }
  return result;
}

std::vector<std::set<Index>> atom_site_compatibility(
    SimpleStructure const &sstruc, BasicStructure const &_prim) {
  std::vector<std::set<Index>> result;
  result.reserve(sstruc.atom_info.names.size());
  for (std::string const &sp : sstruc.atom_info.names) {
    result.push_back({});
    for (Index b = 0; b < _prim.basis().size(); ++b) {
      for (Molecule const &mol : _prim.basis()[b].occupant_dof()) {
        if (mol.contains(sp)) {
          result.back().insert(b);
          break;
        }
      }
    }
  }
  return result;
}
}  // namespace xtal
}  // namespace CASM
