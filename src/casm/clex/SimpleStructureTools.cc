#include "casm/clex/SimpleStructureTools.hh"
#include "casm/crystallography/SimpleStructureTools.hh"
#include "casm/crystallography/SimpleStructure.hh"
#include "casm/crystallography/Coordinate.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/clex/Configuration.hh"
#include "casm/clex/Supercell.hh"
#include "casm/basis_set/DoFTraits.hh"

namespace CASM {
  xtal::SimpleStructure make_simple_structure(Supercell const &_scel,
                                              ConfigDoF const &_dof,
                                              std::vector<DoFKey> const &_which_dofs) {
    xtal::SimpleStructure result;
    result.lat_column_mat = _scel.lattice().lat_column_mat();


    result.mol_info.coords.resize(3, _dof.size());
    result.mol_info.names.reserve(_dof.size());

    for(Index b = 0, l = 0; b < _dof.n_sublat(); ++b) {
      for(Index v = 0; v < _dof.n_vol(); ++v, ++l) {
        result.mol_info.coord(l) = _scel.coord(l).const_cart();
        std::string mol_name = _scel.prim().basis()[ b ].occupant_dof()[_dof.occ(l)].name();
        result.mol_info.names.push_back(std::move(mol_name));
      }
    }

    _apply_dofs(result, _dof, _scel.prim(), _which_dofs);
    return result;
  }

  //***************************************************************************

  SimpleStructure make_simple_structure(Configuration const &_config, std::vector<DoFKey> const &_which_dofs) {
    return make_simple_structure(_config.supercell(), _config.configdof(), _which_dofs);
  }

  //***************************************************************************

  std::vector<std::set<Index> > mol_site_compatibility(SimpleStructure const &sstruc,
                                                       Configuration const &_config) {
    std::vector<std::set<Index> > result;
    result.reserve(sstruc.mol_info.names.size());
    for(std::string const &sp : sstruc.mol_info.names) {
      result.push_back({});
      for(Index l = 0; l < _config.size(); ++l) {
        if(_config.mol(l).name() == sp) {
          result.back().insert(l);
        }
      }
    }
    return result;
  }

  //***************************************************************************

  std::vector<std::set<Index> > atom_site_compatibility(SimpleStructure const &sstruc,
                                                        Configuration const &_config) {
    std::vector<std::set<Index> > result;
    result.reserve(sstruc.atom_info.names.size());
    for(std::string const &sp : sstruc.atom_info.names) {
      result.push_back({});
      for(Index l = 0; l < _config.size(); ++l) {
        if(_config.mol(l).contains(sp)) {
          result.back().insert(l);
        }
      }
    }
    return result;
  }

  //***************************************************************************

  void _apply_dofs(SimpleStructure &_sstruc, ConfigDoF const &_config, BasicStructure const &_reference, std::vector<DoFKey> which_dofs) {
    std::set<TransformDirective> tformers({TransformDirective("atomize")});
    if(which_dofs.empty()) {
      for(std::string const &dof : continuous_local_dof_types(_reference))
        which_dofs.push_back(dof);
      for(std::string const &dof : global_dof_types(_reference))
        which_dofs.push_back(dof);
    }

    for(DoFKey const &dof : which_dofs) {
      if(dof != "none" && dof != "occ")
        tformers.insert(dof);
    }

    //std::cout << "About to transform!!!\n";
    for(TransformDirective const &tformer : tformers) {
      tformer.transform(_config, _reference, _sstruc);
    }
  }

  //***************************************************************************

  TransformDirective::TransformDirective(std::string const &_name) :
    m_name(_name),
    m_traits_ptr(nullptr) {
    if(name() != "atomize") {
      m_traits_ptr = &DoFType::traits(name());
      _accumulate_before({_name}, m_before);
      _accumulate_after({_name}, m_after);
      if(m_after.count("atomize") == 0)
        m_before.insert("atomize");

    }
  }

  //***************************************************************************

  bool TransformDirective::operator<(TransformDirective const &_other) const {
    if(m_before.count(_other.name()) || _other.m_after.count(name())) {
      return false;
    }
    if(m_after.count(_other.name()) || _other.m_before.count(name())) {
      return true;
    }
    return name() < _other.name();
  }

  //***************************************************************************

  void TransformDirective::_accumulate_before(std::set<std::string>const &_queue, std::set<std::string> &_result) const {
    for(std::string const &el : _queue) {
      if(el != name())
        _result.insert(el);
      if(el != "atomize")
        _accumulate_before(AnisoValTraits(el).must_apply_before(), _result);
    }
  }

  //***************************************************************************

  void TransformDirective::_accumulate_after(std::set<std::string>const &_queue, std::set<std::string> &_result) const {
    for(std::string const &el : _queue) {
      if(el != name())
        _result.insert(el);
      if(el != "atomize")
        _accumulate_after(AnisoValTraits(el).must_apply_after(), _result);
    }
  }

  //***************************************************************************

  void TransformDirective::transform(ConfigDoF const  &_dof, BasicStructure const &_reference, SimpleStructure &_struc) const {
    if(m_traits_ptr) {
      if(m_traits_ptr->val_traits().global())
        _struc.properties[m_traits_ptr->name()] = _dof.global_dof(m_traits_ptr->name()).standard_values();
      else {
        _struc.mol_info.properties[m_traits_ptr->name()] = _dof.local_dof(m_traits_ptr->name()).standard_values();
      }
      m_traits_ptr->apply_dof(_dof, _reference, _struc);
    }
    else {
      xtal::_atomize(_struc, _dof.occupation(), _reference);
    }
  }

}
