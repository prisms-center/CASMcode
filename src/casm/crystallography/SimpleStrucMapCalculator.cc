#include "casm/crystallography/SimpleStrucMapCalculator.hh"

#include "casm/crystallography/AnisoValTraits.hh"
#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/BasicStructureTools.hh"
#include "casm/crystallography/Coordinate.hh"
#include "casm/crystallography/DoFDecl.hh"
#include "casm/crystallography/IntegralCoordinateWithin.hh"
#include "casm/crystallography/SimpleStructure.hh"
#include "casm/crystallography/SimpleStructureTools.hh"
#include "casm/crystallography/StrucMapping.hh"
#include "casm/crystallography/UnitCellCoord.hh"
#include "casm/crystallography/io/VaspIO.hh"
#include "casm/misc/algorithm.hh"

namespace CASM {
namespace xtal {
namespace Local {
static Coordinate _make_superlattice_coordinate(
    Index ijk_ix, const Superlattice &superlattice,
    impl::OrderedLatticePointGenerator index_to_ijk_f) {
  UnitCell ijk = index_to_ijk_f(ijk_ix);
  return make_superlattice_coordinate(ijk, superlattice);
}

}  // namespace Local

std::vector<Eigen::Vector3d> SimpleStrucMapCalculator::translations(
    MappingNode const &_node, SimpleStructure const &ichild_struc) const {
  SimpleStructure::Info const &p_info(this->struc_info(parent()));
  SimpleStructure::Info const &c_info(this->struc_info(ichild_struc));

  std::vector<Eigen::Vector3d> result;

  // Find a site in child whose occupant has the lowest number of
  // compatible sites in the parent. This will minimize translations
  Index i_trans(c_info.size()), n_best(p_info.size() + 2);
  for (Index i = 0; i < c_info.size() && n_best > 1; ++i) {
    auto it = _max_n_species().find(c_info.names[i]);
    if (it != _max_n_species().end()) {
      if (it->second < n_best) {
        n_best = it->second;
        i_trans = i;
      }
    } else {
      // child species not known in parent -- incompatible
      return result;
    }
  }

  std::string sp = c_info.names[i_trans];

  result.reserve(n_best);

  // Try translating child atom at i_trans onto each chemically compatible site
  // of parent
  Coordinate equiv_trans(_node.lattice_node.parent.superlattice());
  Coordinate t_trans(equiv_trans);
  Eigen::Vector3d best_equiv_trans;
  for (Index j = 0; j < _allowed_species().size(); ++j) {
    if (contains(_allowed_species()[j], sp)) {
      Eigen::Vector3d translation =
          p_info.cart_coord(j) - c_info.cart_coord(i_trans);
      best_equiv_trans = translation;
      bool unique_trans = true;
      if (this->internal_translations().size() > 1) {
        {
          double best_dist = 1e20;
          for (auto const &itrans : this->internal_translations()) {
            // equivalent translation is internal_translation + translation
            equiv_trans.cart() = itrans + translation;

            // See if equiv_trans is equivalent to any found translation, to
            // within a lattice translation
            for (auto const &rtrans : result) {
              t_trans.cart() = rtrans;
              if (equiv_trans.min_dist(t_trans) < this->xtal_tol()) {
                unique_trans = false;
                break;
              }
            }

            if (!unique_trans) break;

            // keep track of shortest equivalent translation
            equiv_trans.voronoi_within();
            double tdist = equiv_trans.const_cart().norm();
            if (tdist < best_dist) {
              best_dist = tdist;
              best_equiv_trans = equiv_trans.const_cart();
            }
          }
        }
      }
      if (unique_trans) {
        result.push_back(best_equiv_trans);
      }
    }
  }
  return result;
}

//*******************************************************************************************

/// \brief Creates copy of _child_struc by applying isometry, similarity,
/// translation, and site permutation of _node Result has all sites within the
/// unit cell
SimpleStructure SimpleStrucMapCalculator::resolve_setting(
    MappingNode const &_node, SimpleStructure const &_child_struc) const {
  SimpleStructure::Info const &c_info(_child_struc.atom_info);
  SimpleStructure::Info const &p_info((this->parent()).mol_info);
  auto const &cgrid = _node.lattice_node.child;
  auto const &pgrid = _node.lattice_node.parent;

  Index csize = c_info.size();

  SimpleStructure result;

  // Resolve setting of deformed lattice vectors
  // Symmetric deformation of parent lattice that takes it to de-rotated child
  // lattice:
  Eigen::Matrix3d U =
      (_node.lattice_node.stretch)
          .inverse();  // * _node.lattice_node.isometry).inverse();
  result.lat_column_mat = U * pgrid.superlattice().lat_column_mat();

  // Match number of sites in result to number of sites in supercell of parent
  SimpleStructure::Info &r_info(result.mol_info);
  r_info.resize(_node.mol_map.size());

  Eigen::MatrixXd mol_displacement;
  mol_displacement.setZero(3, _node.mol_map.size());

  impl::OrderedLatticePointGenerator parent_index_to_unitcell(
      pgrid.transformation_matrix_to_super());

  // transform and expand the coordinates of child to fill the resolved
  // structure
  {
    for (Index i = 0; i < _node.mol_map.size(); ++i) {
      std::set<Index> const &mol = _node.mol_map[i];

      for (Index j : mol) {
        if (_node.atom_permutation[j] < (csize * cgrid.size())) {
          mol_displacement.col(i) +=
              _node.atom_displacement.col(j) / double(mol.size());

          r_info.names[i] =
              _node.mol_labels[i]
                  .first;  //_child_struc.mol_info.names[j / cgrid.size()];
        } else if (mol.size() > 1)
          throw std::runtime_error(
              "SimpleStrucMapCalculator found multi-atom molecule with "
              "vacancy. This should not occur");
      }
      r_info.cart_coord(i) =
          U * (p_info.cart_coord(i / pgrid.size()) +
               Local::_make_superlattice_coordinate(i % pgrid.size(), pgrid,
                                                    parent_index_to_unitcell)
                   .const_cart() +
               mol_displacement.col(i));
    }
    result.within();
  }  // End coordinate transform

  // Transform and expand properties of child to fill resolved structure
  // global properties first
  for (auto const &el : _child_struc.properties) {
    Eigen::MatrixXd trans = AnisoValTraits(el.first).symop_to_matrix(
        _node.lattice_node.isometry, U * _node.atomic_node.translation,
        _node.atomic_node.time_reversal);
    result.properties.emplace(el.first, trans * el.second);
  }

  // site properties second
  for (auto const &el : c_info.properties) {
    AnisoValTraits ttraits(el.first);
    Eigen::MatrixXd trans = AnisoValTraits(el.first).symop_to_matrix(
        _node.lattice_node.isometry, U * _node.atomic_node.translation,
        _node.atomic_node.time_reversal);
    Eigen::MatrixXd tprop = trans * el.second;

    auto it = r_info.properties
                  .emplace(el.first, Eigen::MatrixXd::Zero(el.second.rows(),
                                                           r_info.size()))
                  .first;
    for (Index i = 0; i < _node.mol_map.size(); ++i) {
      for (Index j : _node.mol_map[i]) {
        Index k = _node.atom_permutation[j];
        if (k < (csize * cgrid.size())) {
          (it->second).col(i) += tprop.col(k / cgrid.size());
        }
        // else -- already checked for improper vacancy conditions above
      }
      if (ttraits.extensive())
        (it->second).col(i) / double(_node.mol_map[i].size());
    }
  }  // End properties transform

  // Add new local property -- site displacement of child relative to parent, at
  // parent strain state
  r_info.properties["disp"] = mol_displacement;

  // Add new global property -- right stretch tensor of child relative to
  // parent:
  {
    // Use AnisoValTraits constructor to get dimension, etc--will throw
    // exception if "Ustrain" is not recognized
    AnisoValTraits ttraits("Ustrain");
    Eigen::VectorXd Uvec(ttraits.dim());
    // Do conversion here for now, since StrainConverter is outside
    // crystallography module.
    Uvec << U(0, 0), U(1, 1), U(2, 2), sqrt(2.) * U(1, 2), sqrt(2.) * U(0, 2),
        sqrt(2.) * U(0, 1);
    result.properties[ttraits.name()] = Uvec;
  }

  // Add new global property -- isometry
  {
    AnisoValTraits ttraits("isometry");
    // unroll to 9-element vector
    result.properties[ttraits.name()] = Eigen::Map<const Eigen::VectorXd>(
        _node.lattice_node.isometry.data(), ttraits.dim());
  }

  return result;
}

//***************************************************************************************************

void SimpleStrucMapCalculator::finalize(
    MappingNode &node, SimpleStructure const &_child_struc,
    bool const &symmetrize_atomic_cost) const {
  populate_displacement(node, _child_struc);
  node.is_valid = this->_assign_molecules(node, _child_struc);
  if (!symmetrize_atomic_cost ||
      this->sym_invariant_displacement_modes().size() == 0)
    node.atomic_node.cost =
        StrucMapping::atomic_cost(node, this->struc_info(_child_struc).size());
  else {
    // Get the parent site indices for each site in the supercell
    auto basis_idx = superstructure_basis_idx(
        (node.lattice_node).parent.transformation_matrix_to_super().cast<int>(),
        parent());
    // Expand the sym invariant basis into this supercell
    std::vector<Eigen::MatrixXd> supercell_sym_invariant_modes(
        this->sym_invariant_displacement_modes().size(),
        Eigen::MatrixXd::Zero(3, basis_idx.size()));

    Eigen::Map<Eigen::VectorXd> disp_map(node.atom_displacement.data(),
                                         node.atom_displacement.size());
    Eigen::MatrixXd sym_removed_disp = node.atom_displacement;
    for (int sym_idx = 0;
         sym_idx < this->sym_invariant_displacement_modes().size(); ++sym_idx) {
      for (int scel_idx = 0; scel_idx < basis_idx.size(); ++scel_idx)
        supercell_sym_invariant_modes[sym_idx].col(scel_idx) =
            this->sym_invariant_displacement_modes()[sym_idx].col(
                basis_idx[scel_idx]);
      supercell_sym_invariant_modes[sym_idx] =
          supercell_sym_invariant_modes[sym_idx] /
          supercell_sym_invariant_modes[sym_idx].norm();
      // Project the displacement mode on to this symmetry invariant mode
      Eigen::VectorXd unrolled_sym_mode = Eigen::Map<Eigen::VectorXd>(
          supercell_sym_invariant_modes[sym_idx].data(),
          supercell_sym_invariant_modes[sym_idx].size());
      double projected_comp = disp_map.dot(unrolled_sym_mode);
      sym_removed_disp =
          sym_removed_disp -
          projected_comp * supercell_sym_invariant_modes[sym_idx];
    }
    node.atom_displacement = sym_removed_disp;
    node.atomic_node.cost =
        StrucMapping::atomic_cost(node, this->struc_info(_child_struc).size());
  }
  node.cost = node.atomic_weight * node.atomic_node.cost +
              node.lattice_weight * node.lattice_node.cost;

  return;
}

//****************************************************************************************************************
//            Assignment Problem methods
//****************************************************************************************************************

void SimpleStrucMapCalculator::populate_displacement(
    MappingNode &_node, SimpleStructure const &child_struc) const {
  const auto &pgrid = _node.lattice_node.parent;
  const auto &cgrid = _node.lattice_node.child;
  SimpleStructure::Info const &p_info(this->struc_info(parent()));
  SimpleStructure::Info const &c_info(this->struc_info(child_struc));
  // TODO: Just use linear index converter? could make things more obvious
  impl::OrderedLatticePointGenerator child_index_to_unitcell(
      cgrid.transformation_matrix_to_super());
  impl::OrderedLatticePointGenerator parent_index_to_unitcell(
      pgrid.transformation_matrix_to_super());

  _node.atom_permutation = _node.atomic_node.permutation();

  // initialize displacement matrix with all zeros
  _node.atom_displacement.setZero(3, pgrid.size() * p_info.size());

  Index cN = c_info.size() * cgrid.size();
  // Populate displacements given as the difference in the Coordinates
  // as described by node.permutation.
  for (Index i = 0; i < _node.atom_permutation.size(); i++) {
    // If we are dealing with a vacancy, its displacment must be zero.
    // if(_node.atom_permutation(i) >= child_struc.n_mol()) {
    //  --DO NOTHING--
    //}

    // Using min_dist routine to calculate the displacement vector that
    // corresponds to the distance used in the Cost Matrix and Hungarian
    // Algorithm The method returns the displacement vector pointing from the
    // IDEAL coordinate to the RELAXED coordinate
    Index j = _node.atom_permutation[i];
    if (j < cN) {
      // Initialize displacement as child coordinate
      Coordinate disp_coord(
          Local::_make_superlattice_coordinate(j % cgrid.size(), cgrid,
                                               child_index_to_unitcell)
              .const_cart(),
          pgrid.superlattice(), CART);

      disp_coord.cart() +=
          c_info.cart_coord(j / cgrid.size()) + _node.atomic_node.translation;

      Coordinate parent_coord = Local::_make_superlattice_coordinate(
          i % pgrid.size(), pgrid, parent_index_to_unitcell);
      parent_coord.cart() += p_info.cart_coord(i / pgrid.size());

      // subtract parent_coord to get displacement
      disp_coord -= parent_coord;
      disp_coord.voronoi_within();

      _node.atom_displacement.col(i) = disp_coord.const_cart();
    }
  }

  Eigen::Vector3d avg_disp =
      _node.atom_displacement.rowwise().sum() / max<double>(double(cN), 1.);

  for (Index i = 0; i < _node.atom_permutation.size(); i++) {
    if (_node.atom_permutation[i] < cN) {
      _node.atom_displacement.col(i) -= avg_disp;
    }
  }

  _node.atomic_node.translation -= avg_disp;
  // End of filling displacements
}

//****************************************************************************************
/*
 * Finding the cost_matrix given the relaxed structure
 * This will always return a square matrix with the extra elements
 * reflecting the vacancies specified in the ideal supercell.
 * Costs are calculated in context of the lattice.
 * cost_matrix(i,j) is cost of mapping child site 'j' onto parent site 'i'
 */
//****************************************************************************************
bool SimpleStrucMapCalculator::populate_cost_mat(
    MappingNode &_node, SimpleStructure const &child_struc) const {
  const auto &pgrid = _node.lattice_node.parent;
  const auto &cgrid = _node.lattice_node.child;

  Eigen::Vector3d const &translation(_node.atomic_node.translation);
  Eigen::MatrixXd &cost_matrix(_node.atomic_node.cost_mat);
  Eigen::Matrix3d metric =
      ((_node.lattice_node.stretch * _node.lattice_node.stretch).inverse() +
       Eigen::Matrix3d::Identity()) /
      2.;
  impl::OrderedLatticePointGenerator child_index_to_unitcell(
      cgrid.transformation_matrix_to_super());
  impl::OrderedLatticePointGenerator parent_index_to_unitcell(
      pgrid.transformation_matrix_to_super());

  SimpleStructure::Info const &p_info(this->struc_info(parent()));
  SimpleStructure::Info const &c_info(this->struc_info(child_struc));

  Index pN = p_info.size() * pgrid.size();
  Index cN = c_info.size() * cgrid.size();

  cost_matrix = Eigen::MatrixXd::Constant(pN, pN, StrucMapping::small_inf());

  if (pN < cN) {
    // std::cout << "Insufficient number of parent sites to accept child
    // atoms\n";
    return false;
  }
  Index residual = pN - cN;
  if (residual > this->max_n_va() * pgrid.size()) {
    // std::cout << "Parent cannot accommodate enough vacancies to make mapping
    // possible.\n";
    return false;
  }

  // index of basis atom in child_struc
  Index bc = 0;
  // index of atom in child supercell
  Index ac = 0;

  Coordinate shifted_parent_coord(pgrid.superlattice());
  Coordinate shifted_child_coord(pgrid.superlattice());

  // loop through all the sites of the child structure
  // As always, each sublattice is traversed contiguously
  for (; bc < c_info.size(); bc++) {
    auto it = this->_max_n_species().find(c_info.names[bc]);
    if (it == this->_max_n_species().end() ||
        (it->second * pgrid.size()) < cgrid.size()) {
      // std::cout << "Parent structure cannot accommodate enough atoms of type
      // '" << c_info.names[bc] << "'.\n"
      //          << "Must accommodate at least " <<cgrid.size()<<" but only
      //          accommodates " << (it->second*pgrid.size()) << "\n";
      return false;
    }

    // For each sublattice, loop over lattice points, 'n'. 'ac' tracks linear
    // index of atoms in child supercell
    for (Index lc = 0; lc < cgrid.size(); ++lc, ++ac) {
      shifted_child_coord.cart() = c_info.cart_coord(bc) +
                                   Local::_make_superlattice_coordinate(
                                       lc, cgrid, child_index_to_unitcell)
                                       .const_cart() +
                                   translation;

      // loop through all the sites in the parent supercell
      Index ap = 0;
      for (Index bp = 0; bp < p_info.size(); ++bp) {
        if (!contains(_allowed_species()[bp], c_info.names[bc])) {
          ap += pgrid.size();
          continue;
        }

        for (Index lp = 0; lp < pgrid.size(); ++lp, ++ap) {
          shifted_parent_coord.cart() =
              p_info.cart_coord(bp) + Local::_make_superlattice_coordinate(
                                          lp, pgrid, parent_index_to_unitcell)
                                          .const_cart();
          cost_matrix(ap, ac) =
              shifted_parent_coord.min_dist2(shifted_child_coord, metric);
        }
      }
    }
  }

  // If there are unvisited columns of cost_mat (because fewer sites in the
  // child structure than in parent), we will treat them as a vacant species in
  // the child struc and set them to zero cost for mapping onto parent sites
  // that allow vacancies

  if (residual) {
    // loop over parent sublats that allow Va, set corresponding block to zero
    // cost
    for (Index bp : this->va_allowed()) {
      cost_matrix.block(bp * pgrid.size(), ac, pgrid.size(), residual)
          .setZero();
    }
  }

  // JCT: I'm not sure if there's an easy way to check if the cost matrix is
  // viable in all cases
  //      Some of the simpler checks I could think of failed for edge cases with
  //      vacancies. If we return an invalid cost matrix, the Hungarian routines
  //      will detect that it is invalid, so maybe there's no point in doing
  //      additional checks here.
  return true;
}

//****************************************************************************************

bool SimpleStrucMapCalculator::_assign_molecules(
    MappingNode &node, SimpleStructure const &_child_struc) const {
  auto const &cgrid(node.lattice_node.child);
  auto const &pgrid(node.lattice_node.parent);
  node.mol_map.clear();
  node.mol_map.reserve(node.atom_permutation.size());
  node.mol_labels.clear();
  node.mol_labels.reserve(node.atom_permutation.size());
  Index j = 0;
  for (Index i : node.atom_permutation) {
    node.mol_map.emplace_back(MappingNode::AtomIndexSet({j}));
    Index bc = i / cgrid.size();
    std::string sp = "Va";
    if (bc < _child_struc.atom_info.size())
      sp = _child_struc.atom_info.names[bc];
    auto occ_itr =
        std::find(this->_allowed_species()[j / pgrid.size()].begin(),
                  this->_allowed_species()[j / pgrid.size()].end(), sp);
    // Index occ_i = find_index(this->_allowed_species()[j / pgrid.size()], sp);
    // if (occ_i == this->_allowed_species()[j / pgrid.size()].size())
    //   return false;
    if (occ_itr == this->_allowed_species()[j / pgrid.size()].end())
      return false;

    node.mol_labels.emplace_back(
        sp, occ_itr - this->_allowed_species()[j / pgrid.size()].begin());
    ++j;
  }

  return true;
}
}  // namespace xtal
}  // namespace CASM
