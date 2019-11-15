#include "casm/crystallography/SimpleStructure.hh"


namespace CASM {
  namespace xtal {

    std::vector<Index> SimpleStructure::Info::sort_by_name() {
      std::vector<Index> permute;
      permute.clear();
      std::map<std::string, std::vector<Index> > smap;
      for(Index i = 0; i < names.size(); ++i) {
        smap[names[i]].push_back(i);
      }
      for(auto const &name : smap) {
        permute.insert(permute.end(), name.second.begin(), name.second.end());
      }
      return permute;
    }

    //***************************************************************************

    SimpleStructure::SimpleStructure(const std::string &_prefix) :
      m_prefix(_prefix) {
    }

    //***************************************************************************

    void SimpleStructure::deform_coords(Eigen::Ref<const Eigen::Matrix3d> const &_F) {
      lat_column_mat = _F * lat_column_mat;
      if(mol_info.coords.rows() == 3)
        mol_info.coords = _F * mol_info.coords;
      if(atom_info.coords.rows() == 3)
        atom_info.coords = _F * atom_info.coords;

    }

    //***************************************************************************

    void SimpleStructure::rotate_coords(Eigen::Ref<const Eigen::Matrix3d> const &_R) {
      lat_column_mat = _R * lat_column_mat;
      if(mol_info.coords.rows() == 3)
        mol_info.coords = _R * mol_info.coords;
      if(atom_info.coords.rows() == 3)
        atom_info.coords = _R * atom_info.coords;

    }

  }
}
