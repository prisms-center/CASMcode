#include "casm/casm_io/json/jsonParser.hh"
#include "casm/casm_io/container/json_io.hh"
#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/symmetry/SymRepTools.hh"

namespace CASM {

  jsonParser &to_json(SymRepTools::IrrepInfo const &irrep, jsonParser &json) {

    if(irrep.pseudo_irrep) {
      to_json_array(irrep.characters.imag(), json["is_direct_sum_of_complex_irreps"]["symop_characters_imag_component"]);
      to_json_array(irrep.characters.real(), json["is_direct_sum_of_complex_irreps"]["symop_characters_real_component"]);
    }


    /*
    if(!almost_zero(irrep.trans_mat.imag())) {
      json["axes"]["imaginary"] = -irrep.trans_mat.imag();
    }

    json["axes"]["real"] = irrep.trans_mat.real();
    */
    /*
    if(irrep.directions.empty()) {
      json["high_symmetry_directions"] = "none";
    }
    else {
      json["high_symmetry_directions"].put_array(irrep.directions.size());
    }

    for(Index i = 0; i < irrep.directions.size(); ++i) {
      json["high_symmetry_directions"][i].put_array(irrep.directions[i].size());
      for(Index j = 0; j < irrep.directions[i].size(); ++j) {
        to_json_array(Eigen::MatrixXd(irrep.trans_mat.real()*irrep.directions[i][j]), json["high_symmetry_directions"][i][j]);
      }
    }
    */
    return json;
  }

  jsonParser &to_json(SymRepTools::SubWedge const &wedge, jsonParser &json) {
    json["full_wedge_axes"] = wedge.trans_mat().transpose();

    for(Index i = 0; i < wedge.irrep_wedges().size(); ++i) {
      std::string irrep_name = "irrep_" + to_sequential_string(i + 1, wedge.irrep_wedges().size());
      json["irrep_wedge_axes"][irrep_name] = wedge.irrep_wedges()[i].axes.transpose();
    }
    return json;
  }

  jsonParser &to_json(VectorSpaceSymReport const &obj, jsonParser &json) {
    json["symmetry_representation"] = obj.symgroup_rep;

    json["glossary"] = obj.axis_glossary;

    std::vector<Index> mults;
    for(auto const &irrep : obj.irreps) {
      if(irrep.index == 0)
        mults.push_back(0);
      mults.back()++;
    }

    Index NQ = obj.symmetry_adapted_dof_subspace.cols();
    Index i(0), l(0), q(1);
    for(Index mult : mults) {
      ++i;
      for(Index m = 0; m < mult; ++m) {
        std::string irrep_name =
          "irrep_"
          + to_sequential_string(i, mults.size())
          + "_"
          + to_sequential_string(m + 1, mult);
        auto const &irrep = obj.irreps[l];
        if(irrep.pseudo_irrep) {
          json["irreducible_representations"]["pseudo_irrep"][irrep_name] = irrep;
        }
        //json["irreducible_representations"][irrep_name]["multiplicity"] = mult;
        if(obj.irreducible_wedge.size()) {
          json["irreducible_representations"]
          ["irreducible_wedge"]
          [irrep_name] = (irrep.trans_mat * obj.irreducible_wedge[0].irrep_wedges()[l].axes).real().transpose();
        }
        json["irreducible_representations"]["irrep_axes"][irrep_name].put_array();
        for(Index a = 0; a < irrep.irrep_dim(); ++a, ++q) {
          std::string axis_name = "q" + to_sequential_string(q, NQ);
          json["irreducible_representations"]["irrep_axes"][irrep_name].push_back(axis_name);
        }

        jsonParser &irrep_matrices = json["irreducible_representations"]["symop_matrices"][irrep_name];//.put_array();
        for(Index o = 0; o < obj.symgroup_rep.size(); ++o) {
          Eigen::MatrixXd const &op = obj.symgroup_rep[i];
          std::string op_name = "op_" + to_sequential_string(o + 1, obj.symgroup_rep.size());
          irrep_matrices[op_name] = (irrep.trans_mat * op * irrep.trans_mat.transpose()).real();
        }

        {
          jsonParser &djson = json["irreducible_representations"]["subgroup_invariant_directions"];
          if(irrep.directions.empty()) {
            djson[irrep_name] = "none";
          }
          else {
            for(Index d = 0; d < irrep.directions.size(); ++d) {
              std::string orbit_name = "direction_orbit_" + to_sequential_string(d + 1, irrep.directions.size());
              djson[irrep_name][orbit_name].put_array(irrep.directions[d].size());
              for(Index j = 0; j < irrep.directions[d].size(); ++j) {
                to_json_array(Eigen::MatrixXd(irrep.trans_mat.real()*irrep.directions[d][j]), djson[irrep_name][orbit_name][j]);
              }
            }
          }
        }

        ++l;
      }
    }

    for(Index q = 0; q < obj.symmetry_adapted_dof_subspace.cols(); ++q) {
      std::string axis_name = "q" + to_sequential_string(q + 1, NQ);
      to_json_array(obj.symmetry_adapted_dof_subspace.col(q),
                    json["irreducible_representations"]["adapted_axes"][axis_name]);
    }


    for(Index i = 0; i < obj.irreducible_wedge.size(); ++i) {
      std::string subwedge_name =
        "subwedge_axes_"
        + to_sequential_string(i + 1, obj.irreducible_wedge.size());
      json["irreducible_wedge"][subwedge_name] = obj.irreducible_wedge[i].trans_mat();
    }

    return json;
  }
}

