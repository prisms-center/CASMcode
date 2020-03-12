#include "casm/crystallography/io/SimpleStructureIO.hh"
#include "casm/crystallography/SimpleStructure.hh"
#include "casm/crystallography/CoordinateSystems.hh"
#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/jsonParser.hh"

namespace {
  using namespace CASM;
  /// Read SimpleStructure::Info for provided species type -- sp="mol" for molecule or sp="atom" for atom
  /// and having the provide prefix
  static void _info_from_json(xtal::SimpleStructure &_struc,
                              const jsonParser &json,
                              Eigen::Matrix3d const &f2c_mat,
                              std::string sp,
                              std::string prefix) {
    xtal::SimpleStructure::Info &sp_info = (sp == "atom" ? _struc.atom_info : _struc.mol_info);
    if(json.contains(sp + "s_per_type")) {
      std::vector<Index> ntype = json[sp + "s_per_type"].get<std::vector<Index> >();

      std::vector<std::string> type = json[sp + "s_type"].get<std::vector<std::string> >();

      for(Index i = 0; i < ntype.size(); ++i) {
        for(Index j = 0; j < ntype[i]; ++j) {
          sp_info.names.push_back(type[i]);
        }
      }
    }
    else if(json.contains(sp + "_types")) {
      from_json(sp_info.names, json[sp + "_types"]);
    }
    else
      return;

    // Remainder of loop body only evaluates if continue statement above is not triggered
    {
      std::vector<std::string> fields({prefix + "basis", "basis", sp + "_coords"});
      for(std::string const &field : fields) {
        auto it = json.find(field);
        if(it != json.end()) {
          sp_info.coords = f2c_mat * it->get<Eigen::MatrixXd>().transpose();
          break;
        }
      }
    }

    {
      std::vector<std::string> fields({sp + "_dofs", sp + "_vals"});
      for(std::string const &field : fields) {
        auto it = json.find(field);
        if(it != json.end()) {
          for(auto it2 = it->begin(); it2 != it->end(); ++it2) {
            sp_info.properties[it2.name()] = (*it2)["value"].get<Eigen::MatrixXd>().transpose();
          }
        }
      }
    }

  }
}

namespace CASM {
  jsonParser &to_json(xtal::SimpleStructure const &_struc,
                      jsonParser &supplement,
                      std::set<std::string> const &excluded_species,
                      std::string prefix) {

    if(!prefix.empty() && prefix.back() != '_')
      prefix.push_back('_');


    std::vector<Index> atom_permute, mol_permute;
    jsonParser &ajson = supplement["atom_type"].put_array();

    for(Index i = 0; i < _struc.atom_info.names.size(); ++i) {
      if(excluded_species.count(_struc.atom_info.names[i]))
        continue;
      ajson.push_back(_struc.atom_info.names[i]);
      atom_permute.push_back(i);
    }

    jsonParser &mjson = supplement["mol_type"].put_array();
    for(Index i = 0; i < _struc.mol_info.names.size(); ++i) {
      if(excluded_species.count(_struc.mol_info.names[i]))
        continue;
      mjson.push_back(_struc.mol_info.names[i]);
      mol_permute.push_back(i);
    }

    supplement[prefix + "lattice"] = _struc.lat_column_mat.transpose();

    for(auto const &dof : _struc.properties) {
      to_json_array(dof.second, supplement[prefix + "global_dofs"][dof.first]["value"]);
    }

    for(auto const &dof : _struc.atom_info.properties) {
      jsonParser &tjson = supplement[prefix + "atom_dofs"][dof.first]["value"].put_array();
      for(Index i : atom_permute)
        tjson.push_back(dof.second.col(i), jsonParser::as_array());
    }

    for(auto const &dof : _struc.mol_info.properties) {
      jsonParser &tjson = supplement[prefix + "mol_dofs"][dof.first]["value"].put_array();
      for(Index i : mol_permute)
        tjson.push_back(dof.second.col(i), jsonParser::as_array());
    }

    {
      jsonParser &tjson = supplement[prefix + "atom_coords"].put_array();
      for(Index i : atom_permute) {
        tjson.push_back(_struc.atom_info.coord(i), jsonParser::as_array());
      }
    }

    {
      jsonParser &tjson = supplement[prefix + "mol_coords"].put_array();
      for(Index i : mol_permute) {
        tjson.push_back(_struc.mol_info.coord(i), jsonParser::as_array());
      }
    }
    return supplement;
  }

  //***************************************************************************

  void from_json(xtal::SimpleStructure &_struc, const jsonParser &json, std::string prefix) {


    Eigen::Matrix3d f2c_mat;
    f2c_mat.setIdentity();

    if(!prefix.empty() && prefix.back() != '_')
      prefix.push_back('_');

    try {
      std::string tstr;
      CASM::from_json(tstr, json["coord_mode"]);

      if(json.contains("lattice")) {
        _struc.lat_column_mat = json["lattice"].get<Eigen::Matrix3d>().transpose();
      }
      else if(json.contains(prefix + "lattice")) {
        _struc.lat_column_mat = json[prefix + "lattice"].get<Eigen::Matrix3d>().transpose();
      }

      COORD_TYPE mode = CART;
      if(tstr == "direct" || tstr == "Direct") {
        mode = FRAC;
        f2c_mat = _struc.lat_column_mat;
      }

      {
        std::vector<std::string> fields({"global_vals", "global_dofs"});
        for(std::string const &field : fields) {
          auto it = json.find(field);
          if(it != json.end()) {
            for(auto it2 = it->begin(); it2 != it->end(); ++it2) {
              _struc.properties[it2.name()] = (*it2)["value"].get<Eigen::MatrixXd>().transpose();
            }
          }
        }
        if(json.contains(prefix + "energy")) {
          _struc.properties["energy"] = json[prefix + "energy"].get<Eigen::MatrixXd>();
        }
        else if(json.contains("energy")) {
          _struc.properties["energy"] = json["energy"].get<Eigen::MatrixXd>();
        }
      }

      if(json.contains(prefix + "forces")) {
        _struc.atom_info.properties["force"] = json[prefix + "forces"].get<Eigen::MatrixXd>().transpose();
      }
      else if(json.contains("forces")) {
        _struc.atom_info.properties["force"] = json["forces"].get<Eigen::MatrixXd>().transpose();
      }

      for(std::string sp : {
            "atom", "mol"
          }) {
        ::_info_from_json(_struc, json, f2c_mat, sp, prefix);
      }

    }
    catch(const std::exception &ex) {
      throw std::runtime_error(std::string("Unable to parse Structure from JSON object.  One or more tags were improperly specified:\n") + ex.what());
    }
  }
}
