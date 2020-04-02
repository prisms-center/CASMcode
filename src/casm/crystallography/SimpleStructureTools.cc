#include "casm/crystallography/SimpleStructure.hh"
#include "casm/crystallography/SimpleStructureTools.hh"
#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/IntegralCoordinateWithin.hh"
#include "casm/crystallography/Site.hh"
#include "casm/crystallography/UnitCellCoord.hh"
#include "casm/external/Eigen/Core"
#include <stdexcept>

namespace CASM {
  namespace xtal {
    namespace Local {

      SimpleStructure::Info _replicate(SimpleStructure::Info const &_info, Index mult) {
        SimpleStructure::Info result;
        result.resize(_info.size()*mult);

        for(Index i = 0; i < _info.size(); ++i)
          result.coords.block(0, i * mult, 3, mult) = _info.cart_coord(i).replicate(1, mult);

        for(auto const &p : _info.properties) {
          result.properties.emplace(p.first, Eigen::MatrixXd(p.second.rows(), mult * p.second.cols()));
          for(Index i = 0; i < p.second.cols(); ++i)
            result.coords.block(0, i * mult, p.second.rows(), mult) = p.second.col(i).replicate(1, mult);
        }

        Index l = 0;
        for(Index b = 0; b < _info.size(); ++b) {
          for(Index g = 0; g < mult; ++g) {
            result.names[l++] = _info.names[b];
          }
        }
        return result;
      }
    }

    //***************************************************************************

    SimpleStructure make_superstructure(Eigen::Ref<const Eigen::Matrix3i> const &_T, SimpleStructure const &_sstruc) {

      SimpleStructure superstructure;
      superstructure.lat_column_mat = _sstruc.lat_column_mat * _T.cast<double>();
      superstructure.properties = _sstruc.properties;

      Lattice sstruc_lattice(_sstruc.lat_column_mat);
      Lattice superstructure_lattice(superstructure.lat_column_mat);

      auto all_lattice_points = make_lattice_points(sstruc_lattice, superstructure_lattice, TOL);

      superstructure.mol_info = Local::_replicate(_sstruc.mol_info, all_lattice_points.size());
      superstructure.atom_info = Local::_replicate(_sstruc.atom_info, all_lattice_points.size());

      Index nm = _sstruc.mol_info.size();
      Index na = _sstruc.atom_info.size();

      Index Nvol = all_lattice_points.size();

      for(Index g = 0; g < Nvol; ++g) {
        Coordinate lattice_point_coordinate = make_superlattice_coordinate(all_lattice_points[g], sstruc_lattice, superstructure_lattice);
        for(Index m = 0; m < nm; ++m) {
          superstructure.mol_info.cart_coord(g + m * Nvol) += lattice_point_coordinate.const_cart();
        }
        for(Index a = 0; a < na; ++a) {
          superstructure.atom_info.cart_coord(g + a * Nvol) += lattice_point_coordinate.const_cart();
        }
      }

      return superstructure;
    }

    //***************************************************************************

    SimpleStructure make_simple_structure(BasicStructure const &_struc) {
      SimpleStructure result;
      result.lat_column_mat = _struc.lattice().lat_column_mat();

      result.mol_info.coords.resize(3, _struc.basis().size());
      result.mol_info.names.reserve(_struc.basis().size());
      Eigen::VectorXi _mol_occ;
      //For now, default to first occupant. In future we may decide
      //to force user to pass mol_occ explicitly
      _mol_occ.setZero(_struc.basis().size());
      for(Index b = 0; b < _struc.basis().size(); ++b) {
        result.mol_info.cart_coord(b) = _struc.basis()[b].const_cart();
        result.mol_info.names.push_back(_struc.basis()[b].occupant_dof()[_mol_occ[b]].name());
      }
      _atomize(result, _mol_occ, _struc);
      return result;
    }

    BasicStructure make_basic_structure(SimpleStructure const &_sstruc,
                                        std::vector<DoFKey> const &_all_dofs,
                                        SimpleStructure::SpeciesMode mode,
                                        std::vector<std::vector<Molecule> > _allowed_occupants) {

      std::map<DoFKey, DoFSet> global_dof;
      std::map<DoFKey, SiteDoFSet> local_dof;
      for(DoFKey const &dof : _all_dofs) {
        if(AnisoValTraits(dof).global()) {
          global_dof.emplace(dof, AnisoValTraits(dof));
        }
        else {
          local_dof.emplace(dof, AnisoValTraits(dof));
        }
      }

      auto const &info = _sstruc.info(mode);
      if(_allowed_occupants.empty())
        _allowed_occupants.resize(info.size());
      for(Index i = 0; i < info.size(); ++i) {
        if(_allowed_occupants[i].empty()) {
          _allowed_occupants[i].push_back(Molecule::make_atom(info.names[i]));
        }
        if(_allowed_occupants[i].size() == 1) {
          std::map<std::string, SpeciesAttribute> attr_map = _allowed_occupants[i][0].attributes();
          for(auto it = attr_map.begin(); it != attr_map.end(); ++it) {
            if(local_dof.count(it->first)) {
              auto er_it = it++;
              attr_map.erase(er_it);
            }
          }

          for(auto const &prop : info.properties) {
            if(local_dof.count(prop.first))
              continue;

            if(prop.first == "disp")
              continue;

            if(!almost_zero(prop.second.col(i)))
              attr_map.emplace(prop.first, SpeciesAttribute(prop.first, prop.second.col(i)));
          }
          _allowed_occupants[i][0].set_attributes(attr_map);
        }
      }

      BasicStructure result(Lattice(_sstruc.lat_column_mat));
      result.set_global_dofs(global_dof);
      std::vector<Site> tbasis(info.size(), Site(result.lattice()));

      for(Index i = 0; i < info.size(); ++i) {
        tbasis[i].cart() = info.cart_coord(i);
        tbasis[i].set_allowed_occupants(std::move(_allowed_occupants[i]));
        tbasis[i].set_dofs(local_dof);
      }

      result.set_basis(tbasis);
      return result;
    }

    //***************************************************************************

    void _atomize(SimpleStructure &_sstruc,
                  Eigen::Ref<const Eigen::VectorXi> const &_mol_occ,
                  BasicStructure const &_reference) {
      Index N_atoms(0);

      Index nb = _reference.basis().size();
      Index nv = _mol_occ.size() / nb;
      Index s = 0;
      for(Index b = 0; b < nb; ++b) {
        for(Index v = 0; v < nv; ++v, ++s) {
          N_atoms += _reference.basis()[b].occupant_dof()[_mol_occ[s]].size();
        }
      }
      _sstruc.atom_info.coords.resize(3, N_atoms);
      _sstruc.atom_info.names.resize(N_atoms);

      // a indexes atom, s indexes site (molecule)
      Index a = 0;
      s = 0;
      for(Index b = 0; b < nb; ++b) {
        for(Index v = 0; v < nv; ++v, ++s) {
          Molecule const &molref = _reference.basis()[b].occupant_dof()[_mol_occ[s]];
          for(Index ms = 0; ms < molref.size(); ++ms, ++a) {
            _sstruc.atom_info.cart_coord(a) = _sstruc.mol_info.cart_coord(s) + molref.atom(ms).cart();
            _sstruc.atom_info.names[a] = molref.atom(ms).name();
          }
        }
      }
    }

    //***************************************************************************

    std::vector<std::set<Index> > mol_site_compatibility(SimpleStructure const &sstruc,
                                                         BasicStructure const &_prim) {
      std::vector<std::set<Index> > result;
      result.reserve(sstruc.mol_info.names.size());
      for(std::string const &sp : sstruc.mol_info.names) {
        result.push_back({});
        for(Index b = 0; b < _prim.basis().size(); ++b) {
          if(_prim.basis()[b].contains(sp)) {
            result.back().insert(b);
          }
        }
      }
      return result;
    }

    std::vector<std::set<Index> > atom_site_compatibility(SimpleStructure const &sstruc,
                                                          BasicStructure const &_prim) {

      std::vector<std::set<Index> > result;
      result.reserve(sstruc.atom_info.names.size());
      for(std::string const &sp : sstruc.atom_info.names) {
        result.push_back({});
        for(Index b = 0; b < _prim.basis().size(); ++b) {
          for(Molecule const &mol : _prim.basis()[b].occupant_dof()) {
            if(mol.contains(sp)) {
              result.back().insert(b);
              break;
            }
          }
        }
      }
      return result;
    }
  }
}
