#include "casm/crystallography/SimpleStructure.hh"
#include "casm/crystallography/SimpleStructureTools.hh"
#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/LatticePointWithin.hh"
#include "casm/crystallography/Site.hh"
#include "casm/crystallography/UnitCellCoord.hh"
#include "casm/external/Eigen/Core"
#include <stdexcept>

namespace CASM {
  namespace xtal {
    namespace Local {

      //***************************************************************************

      SimpleStructure::Info _replicate(SimpleStructure::Info const &_info, Index mult) {
        SimpleStructure::Info result;
        result.resize(_info.size()*mult);

        result.coords = _info.coords.replicate(1, mult);

        for(auto const &p : _info.properties) {
          result.properties[p.first] = p.second.replicate(1, mult);
        }

        Index l = 0;
        for(Index g = 0; g < mult; ++g) {
          for(Index b = 0; b < _info.size(); ++b) {
            result.names[l++] = _info.names[b];
          }
        }
        return result;
      }
    }

    //***************************************************************************

    SimpleStructure make_superstructure(Eigen::Ref<const Eigen::Matrix3i> const &_T, SimpleStructure const &_sstruc) {

      SimpleStructure superstructure;
      superstructure.lat_column_mat = _sstruc.lat_column_mat * _T.cast<double>();;
      superstructure.properties = _sstruc.properties;

      Lattice sstruc_lattice(_sstruc.lat_column_mat);
      Lattice superstructure_lattice(superstructure.lat_column_mat);

      auto all_lattice_points = make_lattice_points(sstruc_lattice, superstructure_lattice, TOL);

      superstructure.mol_info = Local::_replicate(_sstruc.mol_info, all_lattice_points.size());
      superstructure.atom_info = Local::_replicate(_sstruc.atom_info, all_lattice_points.size());

      Index nm = _sstruc.mol_info.size();
      Index na = _sstruc.atom_info.size();

      for(Index g = 0; g < all_lattice_points.size(); ++g) {
        Coordinate lattice_point_coordinate = make_superlattice_coordinate(all_lattice_points[g], sstruc_lattice, superstructure_lattice);
        superstructure.mol_info.coords.block(0, g * nm, 3, nm).colwise() += lattice_point_coordinate.const_cart();
        superstructure.atom_info.coords.block(0, g * na, 3, na).colwise() += lattice_point_coordinate.const_cart();
      }

      return superstructure;
    }

    //***************************************************************************

    SimpleStructure make_simple_structure(BasicStructure const &_struc, const Eigen::VectorXi &current_basis_occupants) {
      assert(_struc.basis().size() == current_basis_occupants.size());
      SimpleStructure result;
      result.lat_column_mat = _struc.lattice().lat_column_mat();

      result.mol_info.coords.resize(3, _struc.basis().size());
      result.mol_info.names.reserve(_struc.basis().size());

      for(Index b = 0; b < _struc.basis().size(); ++b) {
        result.mol_info.coord(b) = _struc.basis()[b].const_cart();
        result.mol_info.names.push_back(_struc.basis()[b].occupant_dof()[current_basis_occupants[b]].name());
      }
      _atomize(result, current_basis_occupants, _struc);
      return result;
    }

    SimpleStructure make_simple_structure(BasicStructure const &_struc) {
      //Make sure all the sites only allow a single occupant
      for(const Site &s : _struc.basis()) {
        if(s.allowed_occupants().size() != 1) {
          throw std::runtime_error("Conversion to SimpleStructure is ambiguous. Basis site has more than one allowed occupant.");
        }
      }

      return make_simple_structure(_struc, Eigen::VectorXi::Zero(_struc.basis().size()));
    }

    //***************************************************************************


    void _atomize(SimpleStructure &_sstruc,
                  Eigen::Ref<const Eigen::VectorXi> const &_mol_occ,
                  BasicStructure const &_reference) {
      Index N_atoms(0);

      Index nb = _reference.basis().size();
      Index nv = _sstruc.mol_info.names.size() / nb;
      Index s = 0;
      for(Index b = 0; b < nb; ++b) {
        for(Index v = 0; v < nv; ++v, ++s) {
          N_atoms += _reference.basis()[b].occupant_dof()[_mol_occ[s]].size();
        }
      }

      //std::cout << "Atomizing with N_atom = " << N_atoms << "; nv = " << nv << "; nb = " << nb << "\n";
      _sstruc.atom_info.coords.resize(3, N_atoms);
      _sstruc.atom_info.names.resize(N_atoms);

      // a indexes atom, s indexes site (molecule)
      Index a = 0;
      s = 0;
      for(Index b = 0; b < nb; ++b) {
        for(Index v = 0; v < nv; ++v, ++s) {
          Molecule const &molref = _reference.basis()[b].occupant_dof()[_mol_occ[s]];
          //std::cout << "(b,v): (" << b << ", " << v << "); molref.size() = " << molref.size() << "\n";
          for(Index ms = 0; ms < molref.size(); ++ms, ++a) {
            _sstruc.atom_info.coord(a) = _sstruc.mol_info.coord(s) + molref.atom(ms).cart();
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

    //***************************************************************************

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
