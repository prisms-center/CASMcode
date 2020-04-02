#include "casm/casm_io/json/jsonParser.hh"
#include "casm/casm_io/container/json_io.hh"
#include "casm/symmetry/SymRepTools.hh"

namespace CASM {
  jsonParser &to_json(SymRepTools::IrrepInfo const &irrep, jsonParser &json) {
    json["further_reducible_as_complex"] = irrep.pseudo_irrep;
    if(irrep.complex) {
      json["irrep_subspace_axes"]["imaginary"] = -irrep.trans_mat.imag();
    }

    json["irrep_subspace_axes"]["real"] = irrep.trans_mat.real();

    if(irrep.complex) {
      to_json_array(irrep.characters.imag(), json["characters"]["imaginary"]);
    }

    to_json_array(irrep.characters.transpose().real(), json["characters"]["real"]);


    json["high_symmetry_directions"].put_array(irrep.directions.size());

    for(Index i = 0; i < irrep.directions.size(); ++i) {
      json["high_symmetry_directions"][i].put_array(irrep.directions[i].size());
      for(Index j = 0; j < irrep.directions[i].size(); ++j) {
        to_json_array(irrep.directions[i][j], json["high_symmetry_directions"][i][j]);
      }
    }
    return json;
  }

  jsonParser &to_json(SymRepTools::SubWedge const &wedge, jsonParser &json) {
    json["full_wedge_axes"] = wedge.trans_mat().transpose();

    for(Index i = 0; i < wedge.irrep_wedges().size(); ++i) {
      json["irrep_wedge_axes"]["irrep_" + std::to_string(i + 1)] = wedge.irrep_wedges()[i].axes.transpose();
    }
    return json;
  }

  jsonParser &to_json(VectorSpaceSymReport const &obj, jsonParser &json) {
    json["matrix_representation"] = obj.symgroup_rep;

    std::vector<Index> mult(obj.irreps.size(), 0);
    Index n = 0;
    for(Index i = obj.irreps.size() - 1; i >= 0; --i) {
      if(n == 0)
        n = obj.irreps[i].index + 1;
      mult[i] = n;
      if(obj.irreps[i].index == 0)
        n = 0;
    }

    for(Index i = 0; i < obj.irreps.size(); ++i) {
      std::string ind_str = std::to_string(i + 1);
      json["irreps"][ind_str] = obj.irreps[i];
      json["irreps"][ind_str]["multiplicity"] = mult[i];
    }

    for(Index i = 0; i < obj.irreducible_wedge.size(); ++i) {
      json["irreducible_wedge"]["subwedge_" + std::to_string(i + 1)] = obj.irreducible_wedge[i];
    }

    json["adapted_axes"] = obj.symmetry_adapted_dof_subspace.transpose();
    return json;
  }
}

