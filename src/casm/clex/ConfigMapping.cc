#include "casm/clex/ConfigMapping.hh"

#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/casm_io/json_io/container.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/Configuration.hh"
#include "casm/clex/Supercell.hh"
#include "casm/clex/ParamComposition.hh"
#include "casm/strain/StrainConverter.hh"
#include "casm/crystallography/Lattice.hh"
#include "casm/crystallography/Niggli.hh"
#include "casm/crystallography/LatticeMap.hh"
#include "casm/crystallography/SupercellEnumerator.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/crystallography/SimpleStructureTools.hh"
#include "casm/symmetry/PermuteIterator.hh"
#include "casm/completer/Handlers.hh"
#include "casm/database/ScelDatabase.hh"

namespace CASM {


  namespace ConfigMapping {

    StrucMapping::AllowedSpecies _allowed_species(BasicStructure<Site> const &_prim,
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

  PrimStrucMapCalculator::PrimStrucMapCalculator(BasicStructure<Site> const &_prim,
                                                 SimpleStructure::SpeciesMode _species_mode/*=StrucMapping::ATOM*/) :
    SimpleStrucMapCalculator(to_simple_structure(_prim),
                             std::vector<SymOp>({SymOp()}),
  _species_mode,
  ConfigMapping::_allowed_species(_prim)),

  m_prim(_prim) {
    SymGroup fg;
    _prim.generate_factor_group(fg);
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
                             double _lattice_weight,
                             double _max_volume_change/*=0.5*/,
                             int options/*=robust*/,
                             double _tol/*=-1.*/) :
    m_pclex(&_pclex),
    m_struc_mapper(PrimStrucMapCalculator(_pclex.prim()),
                   _lattice_weight,
                   1,//_Nbest,
                   //std::vector<SymOp>({SymOp()}),
                   _max_volume_change,
                   options,
                   _tol > 0. ? _tol : _pclex.crystallography_tol()) {

  }

  //*******************************************************************************************

  ConfigMapperResult ConfigMapper::import_structure(SimpleStructure const &child_struc,
                                                    Configuration const *hint_ptr,
                                                    std::vector<DoFKey> const &_hint_dofs) const {
    ConfigMapperResult result;
    double best_cost = StrucMapping::big_inf();

    bool is_new_config(true);

    if(hint_ptr != nullptr) {
      StrucMapper tmapper(*struc_mapper().calculator().quasi_clone(to_simple_structure(*hint_ptr, "", _hint_dofs),
                                                                   make_point_group(hint_ptr->point_group()),
                                                                   SimpleStructure::SpeciesMode::ATOM),
                          struc_mapper().lattice_weight(),
                          1,
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
          SimpleStructure oriented_struc = resolve_setting(map, child_struc);
          ConfigDoF tdof = to_configdof(map, oriented_struc, primclex());
          std::unique_ptr<Configuration> tconfig = notstd::make_unique<Configuration>(scel, jsonParser(), tdof);
          if(strict()) {
            // Strictness transformation reduces permutation swaps, translation magnitude, and isometry character
            PermuteIterator it_strict = strictest_equivalent(scel.sym_info().permute_begin(), scel.sym_info().permute_end(), map);
            MappingNode tnode = copy_apply(it_strict, map);
            tconfig->apply_sym(it_strict);
            result.maps.emplace(tnode, ConfigMapperResult::Individual("", std::move(tconfig)));
          }
          else {
            PermuteIterator it_canon = tconfig->to_canonical();
            MappingNode tnode = copy_apply(it_canon, map);
            tconfig->apply_sym(it_canon);
            result.maps.emplace(tnode, ConfigMapperResult::Individual("", std::move(tconfig)));
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
      std::shared_ptr<Supercell> shared_scel = std::make_shared<Supercell>(&primclex(), map.lat_node.parent.scel_lattice());
      SimpleStructure oriented_struc = resolve_setting(map, child_struc);
      ConfigDoF tdof = to_configdof(map, oriented_struc, primclex());
      std::unique_ptr<Configuration> tconfig = notstd::make_unique<Configuration>(shared_scel, jsonParser(), tdof);
      if(strict()) {
        // Strictness transformation reduces permutation swaps, translation magnitude, and isometry character
        PermuteIterator it_strict = strictest_equivalent(shared_scel->sym_info().permute_begin(), shared_scel->sym_info().permute_end(), map);
        MappingNode tnode = copy_apply(it_strict, map);
        tconfig->apply_sym(it_strict);
        result.maps.emplace(tnode, ConfigMapperResult::Individual("", std::move(tconfig)));
      }
      else {
        PermuteIterator it_canon = tconfig->to_canonical();
        MappingNode tnode = copy_apply(it_canon, map);
        tconfig->apply_sym(it_canon);
        result.maps.emplace(tnode, ConfigMapperResult::Individual("", std::move(tconfig)));
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
