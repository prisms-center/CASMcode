#include "casm/crystallography/SimpleStructure.hh"
#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/crystallography/Site.hh"
#include "casm/clex/Configuration.hh"
#include "casm/clex/Supercell.hh"
#include "casm/casm_io/jsonParser.hh"
#include "casm/basis_set/DoFTraits.hh"
#include "casm/casm_io/stream_io/container.hh"

namespace CASM {

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

  bool TransformDirective::operator<(TransformDirective const &_other) const {
    if(m_before.count(_other.name()) || _other.m_after.count(name()))
      return false;
    if(m_after.count(_other.name()) || _other.m_before.count(name()))
      return true;
    return name() < _other.name();
  }

  void TransformDirective::_accumulate_before(std::set<std::string>const &_queue, std::set<std::string> &_result) const {
    for(std::string const &el : _queue) {
      if(el == "atomize" || _result.count(el))
        continue;
      _result.insert(el);
      _accumulate_before(DoF::traits(el).before_dof_apply(), _result);
    }
  }

  void TransformDirective::_accumulate_after(std::set<std::string>const &_queue, std::set<std::string> &_result) const {
    for(std::string const &el : _queue) {
      if(el == "atomize" || _result.count(el))
        continue;
      _result.insert(el);
      _accumulate_after(DoF::traits(el).after_dof_apply(), _result);
    }
  }


  void TransformDirective::transform(ConfigDoF const  &_config, BasicStructure<Site> const &_reference, SimpleStructure &_struc) const {
    //std::cout << "Applying transformation: " << m_name << "\n";
    if(m_traits_ptr)
      m_traits_ptr->apply_dof(_config, _reference, _struc);
    else
      _struc.atomize(_reference);
  }

  SimpleStructure::SimpleStructure(const std::string &_prefix) :
    selective_dynamics(false),
    m_prefix(_prefix),
    m_atomized(false) {
  }

  SimpleStructure::SimpleStructure(BasicStructure<Site> const &_struc, const std::string &_prefix) :
    lat_column_mat(_struc.lattice().lat_column_mat()),
    selective_dynamics(_struc.selective_dynamics()),
    m_prefix(_prefix),
    m_atomized(false) {

    mol_info.SD = Eigen::MatrixXi::Zero(3, _struc.basis().size());

    mol_info.coords.resize(3, _struc.basis().size());
    mol_info.names.reserve(_struc.basis().size());
    mol_occ.resize(_struc.basis().size());
    for(Index b = 0; b < _struc.basis().size(); ++b) {
      mol_info.coords.col(b) = _struc.basis(b).const_cart();
      mol_info.names.push_back(_struc.basis(b).occ_name());
      mol_occ[b] = _struc.basis(b).site_occupant().value();
    }
    atomize(_struc);
  }

  SimpleStructure::SimpleStructure(Supercell const &_scel, ConfigDoF const &_dof, const std::string &_prefix) :
    lat_column_mat(_scel.lattice().lat_column_mat()),
    selective_dynamics(_scel.prim().selective_dynamics()),
    m_prefix(_prefix),
    m_atomized(false) {

    if(selective_dynamics)
      mol_info.SD = Eigen::MatrixXi::Zero(3, _dof.size());

    mol_info.coords.resize(3, _dof.size());
    mol_info.names.reserve(_dof.size());
    mol_occ = _dof.occupation();
    for(Index b = 0, l = 0; b < _dof.n_sublat(); ++b) {
      for(Index v = 0; v < _dof.n_vol(); ++v, ++l) {
        mol_info.coords.col(l) = _scel.coord(l).const_cart();
        std::string mol_name = _scel.prim().basis()[ b ].site_occupant()[_dof.occ(l)].name();
        mol_info.names.push_back(std::move(mol_name));
      }
    }

    _apply_dofs(_dof, _scel.prim());
  }

  SimpleStructure::SimpleStructure(Configuration const &_config, const std::string &_prefix) :
    SimpleStructure(_config.supercell(), _config.configdof(), _prefix) {}

  void SimpleStructure::deform(Eigen::Ref<const Eigen::Matrix3d> const &_F) {
    lat_column_mat = _F * lat_column_mat;
    if(mol_info.coords.rows() == 3)
      mol_info.coords = _F * mol_info.coords;
    if(atom_info.coords.rows() == 3)
      atom_info.coords = _F * atom_info.coords;

  }

  void SimpleStructure::atomize(BasicStructure<Site> const &_reference) {
    if(m_atomized)
      return;
    m_atomized = true;
    Index N_atoms(0);

    Index nb = _reference.basis().size();
    Index nv = mol_info.names.size() / nb;
    Index s = 0;
    for(Index b = 0; b < nb; ++b) {
      for(Index v = 0; v < nv; ++v, ++s) {
        N_atoms += _reference.basis(b).site_occupant()[mol_occ[s]].size();
      }
    }

    //std::cout << "Atomizing with N_atom = " << N_atoms << "; nv = " << nv << "; nb = " << nb << "\n";
    atom_info.coords.resize(3, N_atoms);
    atom_info.names.resize(N_atoms);

    atom_info.SD = Eigen::MatrixXi::Zero(3, N_atoms);

    // a indexes atom, s indexes site (molecule)
    Index a = 0;
    s = 0;
    for(Index b = 0; b < nb; ++b) {
      for(Index v = 0; v < nv; ++v, ++s) {
        Molecule const &molref = _reference.basis(b).site_occupant()[mol_occ[s]];
        //std::cout << "(b,v): (" << b << ", " << v << "); molref.size() = " << molref.size() << "\n";
        for(Index ms = 0; ms < molref.size(); ++ms, ++a) {
          atom_info.coords.col(a) = mol_info.coords.col(s) + molref.atom(ms).cart();
          atom_info.names[a] = molref.atom(ms).name();
          if(selective_dynamics) {
            atom_info.SD.col(a) = mol_info.SD.col(s);
            for(Index i = 0; i < 3; ++i) {
              if(molref.atom(i).sd_flag()[i])
                atom_info.SD(i, a) = 1;
            }
          }
        }
      }
    }
    //std::cout << "atom_info.names: " << atom_info.names << "\n"
    //        << "atom_info.coords: \n" << atom_info.coords << "\n";
  }

  jsonParser &to_json(SimpleStructure const &_struc, jsonParser &json, std::set<std::string> const &excluded_species) {

    std::string prefix = _struc.prefix();
    if(!prefix.empty())
      prefix.push_back('_');


    jsonParser &ajson = json["atom_type"].put_array();
    _struc.atom_info.permute.clear();
    for(Index i = 0; i < _struc.atom_info.names.size(); ++i) {
      if(excluded_species.count(_struc.atom_info.names[i]))
        continue;
      ajson.push_back(_struc.atom_info.names[i]);
      _struc.atom_info.permute.push_back(i);
    }

    jsonParser &mjson = json["mol_type"].put_array();
    _struc.mol_info.permute.clear();
    for(Index i = 0; i < _struc.mol_info.names.size(); ++i) {
      if(excluded_species.count(_struc.mol_info.names[i]))
        continue;
      mjson.push_back(_struc.mol_info.names[i]);
      _struc.mol_info.permute.push_back(i);
    }

    json[prefix + "lattice"] = _struc.lat_column_mat.transpose();
    json[prefix + "global_dofs"] = _struc.global_dofs;

    {
      auto
      it = _struc.atom_info.dofs.begin(),
      end_it = _struc.atom_info.dofs.end();

      for(; it != end_it; ++it) {
        auto
        it2 = it->begin(),
        end_it2 = it->end();
        for(; it2 != end_it2; ++it2) {
          if(it2->is_array()) {
            jsonParser &tjson = json[prefix + "atom_dofs"][it.name()][it2.name()].put_array();
            for(Index i : _struc.atom_info.permute)
              tjson.push_back((*it2)[i]);
          }
          else {
            json[prefix + "atom_dofs"][it.name()][it2.name()] = *it2;
          }
        }
      }
    }

    {
      auto
      it = _struc.mol_info.dofs.begin(),
      end_it = _struc.mol_info.dofs.end();

      for(; it != end_it; ++it) {
        auto
        it2 = it->begin(),
        end_it2 = it->end();
        for(; it2 != end_it2; ++it2) {
          if(it2->is_array()) {
            jsonParser &tjson = json[prefix + "mol_dofs"][it.name()][it2.name()].put_array();
            for(Index i : _struc.mol_info.permute)
              tjson.push_back((*it2)[i]);
          }
          else {
            json[prefix + "mol_dofs"][it.name()][it2.name()] = *it2;
          }
        }
      }
    }


    //Why does selective_dynamics default to ON?;
    if(_struc.selective_dynamics) {
      json["selective_dynamics"] = _struc.selective_dynamics;
      json["atom_selective_dynamics"].put_array();
      json["mol_selective_dynamics"].put_array();
    }


    {
      jsonParser &tjson = json[prefix + "atom_coords"].put_array();
      for(Index i : _struc.atom_info.permute) {
        tjson.push_back(_struc.atom_info.coords.col(i).transpose(), jsonParser::as_array());
        if(_struc.selective_dynamics)
          json["atom_selective_dynamics"].push_back(_struc.atom_info.SD.col(i).transpose(), jsonParser::as_array());
      }
    }

    {
      jsonParser &tjson = json[prefix + "mol_coords"].put_array();
      for(Index i : _struc.mol_info.permute) {
        tjson.push_back(_struc.mol_info.coords.col(i).transpose(), jsonParser::as_array());
        if(_struc.selective_dynamics)
          json["mol_selective_dynamics"].push_back(_struc.mol_info.SD.col(i).transpose(), jsonParser::as_array());
      }
    }
    return json;
  }

  void from_json(SimpleStructure &_struc, const jsonParser &json) {
    std::string prefix = _struc.prefix();
    if(!prefix.empty())
      prefix.push_back('_');

    try {
      std::string tstr;
      CASM::from_json(tstr, json["coord_mode"]);

      COORD_TYPE mode = CART;
      if(tstr == "direct" || tstr == "Direct")
        mode = FRAC;

      _struc.lat_column_mat = json[prefix + "lattice"].get<Eigen::Matrix3d>().transpose();
      if(json.contains(prefix + "global_dofs"))
        _struc.global_dofs = json[prefix + "global_dofs"];

      if(json.contains("atoms_per_type")) {
        std::vector<Index> ntype = json["atoms_per_type"].get<std::vector<Index> >();
        std::vector<std::string> type = json["atoms_type"].get<std::vector<std::string> >();

        for(Index i = 0; i < ntype.size(); ++i) {
          for(Index j = 0; j < ntype[i]; ++j) {
            _struc.atom_info.names.push_back(type[i]);
          }
        }

        if(mode == FRAC)
          _struc.atom_info.coords = _struc.lat_column_mat * json[prefix + "atom_coords"].get<Eigen::MatrixXd>().transpose();
        else
          _struc.atom_info.coords = json[prefix + "atom_coords"].get<Eigen::MatrixXd>().transpose();
      }

      if(json.contains("mols_per_type")) {
        std::vector<Index> ntype = json["mols_per_type"].get<std::vector<Index> >();
        std::vector<std::string> type = json["mols_type"].get<std::vector<std::string> >();

        for(Index i = 0; i < ntype.size(); ++i) {
          for(Index j = 0; j < ntype[i]; ++j) {
            _struc.mol_info.names.push_back(type[i]);
          }
        }

        if(mode == FRAC)
          _struc.mol_info.coords = _struc.lat_column_mat * json[prefix + "mol_coords"].get<Eigen::MatrixXd>().transpose();
        else
          _struc.mol_info.coords = _struc.lat_column_mat * json[prefix + "mol_coords"].get<Eigen::MatrixXd>().transpose();
      }

      if(json.contains(prefix + "atom_dofs"))
        _struc.atom_info.dofs = json[prefix + "atom_dofs"];

      if(json.contains(prefix + "mol_dofs"))
        _struc.mol_info.dofs = json[prefix + "mol_dofs"];

      json.get_if(_struc.selective_dynamics, "selective_dynamics");
      if(_struc.selective_dynamics) {
        if(json.contains("atom_selective_dynamics")) {
          _struc.atom_info.SD = json["atom_selective_dynamics"].get<Eigen::MatrixXi>().transpose();
        }
        if(json.contains("mol_selective_dynamics")) {
          _struc.mol_info.SD = json["mol_selective_dynamics"].get<Eigen::MatrixXi>().transpose();
        }
      }

    }
    catch(const std::exception &ex) {
      throw std::runtime_error(std::string("Unable to parse Structure from JSON object.  One or more tags were improperly specified:\n") + ex.what());
    }
  }


  void SimpleStructure::_apply_dofs(ConfigDoF const &_config, BasicStructure<Site> const &_reference) {
    std::set<TransformDirective> tformers({TransformDirective("atomize")});
    for(std::string const &dof : continuous_local_dof_types(_reference))
      tformers.insert(dof);
    for(std::string const &dof : global_dof_types(_reference))
      tformers.insert(dof);

    //std::cout << "About to transform!!!\n";
    for(TransformDirective const &tformer : tformers) {
      tformer.transform(_config, _reference, *this);
    }
  }

}
