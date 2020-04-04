#include "casm/casm_io/json/jsonParser.hh"
#include "casm/casm_io/container/json_io.hh"
#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/symmetry/SymRepTools.hh"

namespace CASM {
  namespace Local {
    static std::string _to_sequential_string(Index i, Index max_i) {
      max_i = max(i, max_i);
      Index length = 1;
      while(max_i /= 10)
        length++;

      std::string tresult = std::to_string(i);

      std::string result(length - tresult.size(), '0');
      return result.append(tresult);
    }

  }

  jsonParser &to_json(SymRepTools::IrrepInfo const &irrep, jsonParser &json) {
    if(!almost_zero(irrep.characters.imag())) {
      to_json_array(irrep.characters.imag(), json["characters"]["imaginary"]);
    }
    to_json_array(irrep.characters.transpose().real(), json["characters"]["real"]);

    json["reducible_as_complex"] = irrep.pseudo_irrep;
    if(irrep.pseudo_irrep) {
      to_json_array(2 * irrep.characters.real(), json["characters"]["reducible_real_characters"]);
    }

    if(!almost_zero(irrep.trans_mat.imag())) {
      json["axes"]["imaginary"] = -irrep.trans_mat.imag();
    }

    json["axes"]["real"] = irrep.trans_mat.real();


    if(irrep.directions.empty()) {
      json["high_symmetry_directions"] = "none";
    }
    else {
      json["high_symmetry_directions"].put_array(irrep.directions.size());
    }

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
      std::string irrep_name = "irrep_" + Local::_to_sequential_string(i + 1, wedge.irrep_wedges().size());
      json["irrep_wedge_axes"][irrep_name] = wedge.irrep_wedges()[i].axes.transpose();
    }
    return json;
  }

  jsonParser &to_json(VectorSpaceSymReport const &obj, jsonParser &json) {
    json["symmetry_representation"] = obj.symgroup_rep;

    std::vector<Index> mults;
    for(auto const &irrep : obj.irreps) {
      if(irrep.index == 0)
        mults.push_back(0);
      mults.back()++;
    }

    Index i(0), l(0);
    for(Index mult : mults) {
      ++i;
      for(Index m = 0; m < mult; ++m) {
        std::string irrep_name =
          "irrep_"
          + Local::_to_sequential_string(i, mults.size())
          + "_"
          + Local::_to_sequential_string(m + 1, mult);

        json["irreducible_representations"][irrep_name] = obj.irreps[l];
        json["irreducible_representations"][irrep_name]["multiplicity"] = mult;
        if(obj.irreducible_wedge.size()) {
          json["irreducible_representations"]
          [irrep_name]
          ["irreducible_wedge"] = obj.irreducible_wedge[0].irrep_wedges()[l].axes.transpose();
        }
        ++l;
      }
    }
    json["irreducible_representations"]["adapted_axes"] = obj.symmetry_adapted_dof_subspace.transpose();


    for(Index i = 0; i < obj.irreducible_wedge.size(); ++i) {
      std::string subwedge_name =
        "subwedge_axes_"
        + Local::_to_sequential_string(i + 1, obj.irreducible_wedge.size());
      json["irreducible_wedge"][subwedge_name] = obj.irreducible_wedge[i].trans_mat();
    }

    return json;
  }
}

