#include "casm/clex/ConfigMapping.hh"

#include "casm/crystallography/Adapter.hh"
#include "casm/crystallography/BasicStructureTools.hh"
#include "casm/crystallography/SymType.hh"
#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/casm_io/container/json_io.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/Configuration.hh"
#include "casm/clex/Supercell.hh"
#include "casm/clex/ParamComposition.hh"
#include "casm/strain/StrainConverter.hh"
#include "casm/crystallography/Lattice.hh"
#include "casm/crystallography/Niggli.hh"
#include "casm/crystallography/LatticeMap.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/crystallography/SimpleStructureTools.hh"
#include "casm/symmetry/PermuteIterator.hh"
#include "casm/completer/Handlers.hh"
#include "casm/database/ScelDatabase.hh"

namespace CASM {

  namespace Local {
    int permute_dist(std::vector<Index> const &_perm) {
      int result(0);
      for(int i = 0; i < _perm.size(); ++i) {
        result += std::abs(int(_perm[i] - i));
      }
      return result;
    }

    //*******************************************************************************************

    template<typename IterType>
    static IterType strictest_equivalent(IterType begin, IterType end, MappingNode const &_node) {
      SymOp op(_node.isometry(), _node.translation(), false, _node.tol());
      SymOp t_op;
      IterType best_it = begin;

      double tchar, tdist;
      int tdet, tpdist;
      double best_char = op.matrix().trace();
      double best_dist = op.tau().norm();
      int best_det = sgn(round(op.matrix().determinant()));
      int best_pdist = permute_dist(_node.permutation);

      Coordinate tau(_node.lat_node.parent.superlattice());
      while(begin != end) {
        t_op = begin->sym_op() * op;
        tdet = sgn(round(t_op.matrix().determinant()));
        bool skip_fg = false;
        if(tdet > best_det) {
          best_det = tdet;
          best_char = tdet * t_op.matrix().trace();
          best_pdist = permute_dist(begin->combined_permute() * _node.permutation);
          tau.cart() = t_op.tau();
          tau.voronoi_within();
          best_dist = tau.const_cart().norm();
          best_it = begin;
        }
        else if(tdet == best_det) {
          tchar = tdet * t_op.matrix().trace();
          if(almost_equal(tchar, best_char)) {
            tpdist = permute_dist(begin->combined_permute() * _node.permutation);
            if(tpdist > best_pdist) {
              best_det = tdet;
              best_char = tchar;
              best_pdist = tpdist;
              tau.cart() = t_op.tau();
              tau.voronoi_within();
              best_dist = tau.const_cart().norm();
              best_it = begin;
            }
            else if(tpdist == best_pdist) {
              tau.voronoi_within();
              tdist = tau.const_cart().norm();
              if(tdist < best_dist) {
                best_det = tdet;
                best_char = tchar;
                best_pdist = tpdist;
                best_dist = tdist;
                best_it = begin;
              }
            }
          }
          else if(tchar > best_char) {
            best_det = tdet;
            best_char = tchar;
            best_pdist = permute_dist(begin->combined_permute() * _node.permutation);
            tau.cart() = t_op.tau();
            tau.voronoi_within();
            best_dist = tau.const_cart().norm();
            best_it = begin;
          }
          else {
            skip_fg = true;
          }
        }
        else {
          skip_fg = true;
        }

        if(skip_fg) {
          begin = begin_next_fg_op(begin, end);
        }
        else {
          ++begin;
        }
      }
      return best_it;
    }
  }

  //*******************************************************************************************
  Index ConfigMapperResult::n_optimal(double tol/*=TOL*/) const {
    Index result = 0;
    auto it = maps.begin();
    while(it != maps.end() && almost_equal((it->first).cost, (maps.begin()->first).cost, tol))
      ++result;
    return result;
  }

  //*******************************************************************************************
  SimpleStructure PrimStrucMapCalculator::resolve_setting(MappingNode const &_node,
                                                          SimpleStructure const &_child_struc) const {

    throw std::runtime_error("ConfigMapping::resolve_setting is not implemented!");
    return _child_struc;
  }
  //*******************************************************************************************
  namespace ConfigMapping {

    xtal::StrucMapping::AllowedSpecies _allowed_species(BasicStructure const &_prim,
                                                        SimpleStructure::SpeciesMode _species_mode = SimpleStructure::SpeciesMode::ATOM) {
      std::vector<std::unordered_set<std::string> > result(_prim.basis().size());
      Index i = 0;
      for(Site const &site : _prim.basis()) {
        for(Molecule const &mol : site.occupant_dof().domain()) {
          if(_species_mode == SimpleStructure::SpeciesMode::MOL) {
            result[i].insert(mol.name());
          }
          else if(_species_mode == SimpleStructure::SpeciesMode::ATOM) {
            if(mol.size() != 1) {
              throw std::runtime_error("ConfigMapping::_allowed_species may only be called on structures with single-atom species.");
            }
            result[i].insert(mol.atom(0).name());
          }
        }
        ++i;
      }
      return result;
    }
  }

  //*******************************************************************************************

  /// \brief Reorders the permutation and compounds the spatial isometry (rotation + translation) of _node with that of _it
  MappingNode copy_apply(PermuteIterator const &_it, MappingNode const &_node, bool transform_cost_mat) {
    MappingNode result(_node);
    SymOp op = _it.sym_op();

    //Apply symmetry to LatticeNode:
    // lat_node.parent and lat_node.child remain unchanged, but isometry and stretch tensors are augmented
    //   parent superlattice Lp and child superlattice Lc are related via
    //      Lp = U * R * Lc  (where R is isometry and U is right stretch tensor of _node)
    //   'op' must be in point group of Lp, and invariance relation for Lp is
    //      op.matrix() * Lp * N = Lp  (where 'N' is integer fractional symop matrix)
    //   Substituting above expression and inserting Identity = op.matrix().transpose() * op.matrix() yields:
    //      Lp = (op.matrix() * U * op.matrix().transpose()) * (op.matrix() * R) * Lc * N
    //   where terms are grouped to reveal the transformation rules for matrices 'U' and 'R' for the symmetry-related mapping
    //   The fractional matrix 'N' is not used in this routine, since _node records the child lattice in its undeformed state
    //   and all atomic coordinates are recorded in undeformed cartesian coordinates
    result.lat_node.isometry = op.matrix() * _node.isometry();
    result.lat_node.stretch = op.matrix() * _node.stretch() * op.matrix().transpose();

    //Apply symmetry to HungarianNode:
    // parent coordinates Cp (3xN) and child coordinates Cc (3xN) are related via the mapping relation
    //    Cp = U * R * Cc * P.transpose() - D + T)]
    // where R is isometry, U is right stretch tensor, P is site permutation,
    // D is displacement field (3xN) and T is mapping translation (3x1 repeated for N columns)
    // 'op' must be in space group of Cp, and invariance relation for Cp is
    //    Cp = op.matrix() * Cp * Ps.transpose() + op.tau()
    // where Ps = _it.combined_permute().matrix() and op.tau() is added col-wise.
    // The invariance is only valid to within a lattice translation of the basis sites (which does not affect mapping score)
    // Inserting the mapping relation for Cp into the invariance relation yields
    //    Cp = [op.matrix() * U * op.matrix.transpose()] * [op.matrix * R] * Cc * [P.transpose() * Ps.transpose()] - [op.matrix() * D * Ps] + [op.matrix() * T + op.tau()]
    // where terms are grouped to reveal the transformation rules for mapping permutation 'P', displacement field 'D', and mapping translation 'T'
    // The transormation rules for 'U' and 'R' are identical to those above and are recorded in _node.lat_node.
    // These five translation rules specify all considerations necessary to describe how a symop in the space group of the parent crystal
    // relates a mapping of the child crystal onto the parent crystal to an equivalent mapping
    result.basis_node.translation = op.matrix() * _node.translation() + op.tau();
    Eigen::MatrixXd td = op.matrix() * result.displacement;
    for(Index i = 0; i < result.permutation.size(); ++i) {
      result.permutation[i] = _node.permutation[_it.permute_ind(i)];
      result.displacement.col(i) = td.col(_it.permute_ind(i));
    }

    //Attempt to transfrom the constrained mapping problem and cost matrix
    //This is distinct from the relations described above, as the asignments are not yet all known
    //Instead of permuting the child indices, we will permute the parent indices by the inverse
    //permutation, which should have the same effect in the end
    if(transform_cost_mat) {
      PermuteIterator inv_it = _it.inverse();
      for(Index i = 0; i < result.basis_node.irow.size(); ++i)
        result.basis_node.irow[i] = inv_it.permute_ind(_node.basis_node.irow[i]);

      result.basis_node.forced_on.clear();
      for(auto const &el : _node.basis_node.forced_on)
        result.basis_node.forced_on.emplace(inv_it.permute_ind(el.first), el.second);
      //No need to transform assignment vector or cost_mat, which are in terms of the nominal indexing
    }
    return result;
  }

  //*******************************************************************************************
  /// \brief Initializes configdof corresponding to a mapping (encoded by _node) of _child_struc onto _pclex
  ConfigDoF to_configdof(MappingNode const _node, SimpleStructure const &_child_struc, Supercell const  &_scel) {
    std::cerr << "to_configdof not fully implemented\n";
    exit(1);
    return _scel.zero_configdof(TOL);
  }
  //*******************************************************************************************

  PrimStrucMapCalculator::PrimStrucMapCalculator(BasicStructure const &_prim,
                                                 SimpleStructure::SpeciesMode _species_mode/*=StrucMapping::ATOM*/) :
    SimpleStrucMapCalculator(make_simple_structure(_prim),
                             std::vector<SymOp>({
    SymOp()
  }),
  _species_mode,
  ConfigMapping::_allowed_species(_prim)),

  m_prim(_prim) {
    xtal::SymOpVector factor_group_operations = xtal::make_factor_group(_prim);
    SymGroup fg = adapter::Adapter<SymGroup, xtal::SymOpVector>()(factor_group_operations, _prim.lattice());
    set_point_group(fg.copy_no_trans());

  }

  //*******************************************************************************************

  void ConfigMapper::add_allowed_lattices(std::vector<std::string> const &_lattice_names) {
    for(std::string const &_name : _lattice_names) {
      auto it = primclex().template db<Supercell>().find(_name);
      if(it == primclex().template db<Supercell>().end())
        throw std::runtime_error("Could not add mapping lattice constraint " + _name + " because no supercell having that name exists in database.\n");
      struc_mapper().add_allowed_lattice(it->lattice());
    }
  }


  //*******************************************************************************************

  ConfigMapper::ConfigMapper(PrimClex const &_pclex,
                             double _strain_weight,
                             double _max_volume_change/*=0.5*/,
                             int options/*=robust*/,
                             double _tol/*=-1.*/) :
    m_pclex(&_pclex),
    m_struc_mapper(PrimStrucMapCalculator(_pclex.prim()),
                   _strain_weight,
                   _max_volume_change,
                   options,
                   _tol > 0. ? _tol : _pclex.crystallography_tol()) {

  }

  //*******************************************************************************************

  ConfigMapperResult ConfigMapper::import_structure(SimpleStructure const &child_struc,
                                                    Configuration const *hint_ptr,
                                                    std::vector<DoFKey> const &_hint_dofs) const {
    ConfigMapperResult result;
    double best_cost = xtal::StrucMapping::big_inf();

    //bool is_new_config(true);

    if(hint_ptr != nullptr) {
      StrucMapper tmapper(*struc_mapper().calculator().quasi_clone(xtal::make_simple_structure(*hint_ptr, _hint_dofs),
                                                                   make_point_group(hint_ptr->point_group(), hint_ptr->supercell().sym_info().supercell_lattice()),
                                                                   SimpleStructure::SpeciesMode::ATOM),
                          struc_mapper().strain_weight(),
                          0.,
                          struc_mapper().options(),
                          struc_mapper().tol());

      auto config_maps = tmapper.map_deformed_struc_impose_lattice(child_struc,
                                                                   hint_ptr->ideal_lattice());

      // Refactor into external routine A. This is too annoying with the current way that supercells are managed
      if(!config_maps.empty()) {
        best_cost = config_maps.begin()->cost;
        const Supercell &scel(hint_ptr->supercell());
        for(auto const &map : config_maps) {
          SimpleStructure oriented_struc = struc_mapper().calculator().resolve_setting(map, child_struc);
          ConfigDoF tdof = to_configdof(map, oriented_struc, scel);
          Configuration tconfig(scel, jsonParser(), tdof);
          if(strict()) {
            // Strictness transformation reduces permutation swaps, translation magnitude, and isometry character
            PermuteIterator it_strict = Local::strictest_equivalent(scel.sym_info().permute_begin(), scel.sym_info().permute_end(), map);
            MappingNode tnode = copy_apply(it_strict, map);
            tconfig.apply_sym(it_strict);
            result.maps.emplace(std::make_pair(tnode, std::make_pair(ConfigMapperResult::MapData(""), std::move(tconfig))));
          }
          else {
            PermuteIterator it_canon = tconfig.to_canonical();
            MappingNode tnode = copy_apply(it_canon, map);
            tconfig.apply_sym(it_canon);
            result.maps.emplace(std::make_pair(tnode, std::make_pair(ConfigMapperResult::MapData(""), std::move(tconfig))));
          }
        }
      }
      //\End routine A
      //TODO: Check equivalence with hint_ptr
    }


    auto struc_maps = struc_mapper().map_deformed_struc(child_struc,
                                                        best_cost + struc_mapper().tol());
    // Refactor into external routine A. This is too annoying with the current way that supercells are managed
    for(auto const &map : struc_maps) {
      std::shared_ptr<Supercell> shared_scel = std::make_shared<Supercell>(&primclex(), map.lat_node.parent.superlattice());
      SimpleStructure oriented_struc = struc_mapper().calculator().resolve_setting(map, child_struc);
      ConfigDoF tdof = to_configdof(map, oriented_struc, *shared_scel);
      Configuration tconfig(shared_scel, jsonParser(), tdof);
      if(strict()) {
        // Strictness transformation reduces permutation swaps, translation magnitude, and isometry character
        PermuteIterator it_strict = Local::strictest_equivalent(shared_scel->sym_info().permute_begin(), shared_scel->sym_info().permute_end(), map);
        MappingNode tnode = copy_apply(it_strict, map);
        tconfig.apply_sym(it_strict);
        result.maps.emplace(std::make_pair(tnode, std::make_pair(ConfigMapperResult::MapData(""), std::move(tconfig))));
      }
      else {
        PermuteIterator it_canon = tconfig.to_canonical();
        MappingNode tnode = copy_apply(it_canon, map);
        tconfig.apply_sym(it_canon);
        result.maps.emplace(std::make_pair(tnode, std::make_pair(ConfigMapperResult::MapData(""), std::move(tconfig))));
      }
    }
    //\End routine A



    // calculate and store:
    // - 'relaxation_deformation'
    // - 'relaxation_displacement'
    // - cart_op
    // transform deformation tensor to match canonical form and apply operation to cart_op
    //result.success = true;

    return result;
  }
}
