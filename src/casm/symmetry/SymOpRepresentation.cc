#include "casm/symmetry/SymOpRepresentation.hh"
#include "casm/symmetry/SymGroupRep.hh"
#include "casm/symmetry/SymMatrixXd.hh"
#include "casm/symmetry/SymPermutation.hh"

namespace CASM {

  void SymOpRepresentation::set_rep(const MasterSymGroup &new_group, Index new_rep_ID) {
    m_master_group = &new_group;
    rep_ID = new_rep_ID;
    SymGroupRep const *trep(new_group.representation(rep_ID));
    Index i;
    for(i = 0; i < trep->size(); i++) {
      if(this == trep->at(i)) {
        op_index = i;
        break;
      }
    }
    if(i == new_group.size()) op_index = -1;

    return;
  }

  //**********************************************************

  void SymOpRepresentation::set_identifiers(const MasterSymGroup &new_group, Index new_rep_ID, Index new_op_index) {
    m_master_group = &new_group;
    rep_ID = new_rep_ID;
    op_index = new_op_index;

    return;
  }

  //**********************************************************


  /// creates jsonParser using polymorphism
  jsonParser &to_json(const SymOpRepresentation *rep, jsonParser &json) {
    return rep->to_json(json);
  }

  //**********************************************************

  /// This allocates a new object to 'rep'.
  ///   It might need a Lattice
  ///
  void from_json(SymOpRepresentation *rep, const jsonParser &json, const Lattice &lat) {
    try {
      if(json["SymOpRep_type"] == "SymPermutation") {

        // prepare a SymPermutation and then read from json
        Array<Index> perm;
        SymPermutation trep(perm);
        CASM::from_json(trep, json);

        // copy to rep
        rep = new SymPermutation(trep);

      }
      else if(json["SymOpRep_type"] == "SymMatrixXd") {

        // prepare a SymMatrixXd and then read from json
        Eigen::MatrixXd mat;
        SymMatrixXd trep(mat);
        CASM::from_json(trep, json);

        // copy to rep
        rep = new SymMatrixXd(trep);

      }
      else if(json["SymOpRep_type"] == "SymOp") {

        // prepare a SymMatrixXd and then read from json
        Coordinate coord(lat);
        SymOp trep(coord);
        CASM::from_json(trep, json);

        // copy to rep
        rep = new SymOp(trep);

      }
      else {
        std::cout << "Error in 'void from_json(SymOpRepresentation *rep, const jsonParser &json, const Lattice &lat)'" << std::endl;
        std::cout << "Unrecognized 'SymOpRep_type': '" << json["SymOpRep_type"] << "'." << std::endl;
        std::cout << "Options are: 'SymPermutation', 'SymMatrixXd', or 'SymOp'." << std::endl;
        exit(1);
      }
    }
    catch(...) {
      /// re-throw exceptions
      throw;
    }
  }
  //**********************************************************

  //enum symmetry_type {identity_op, mirror_op, glide_op, rotation_op, screw_op, inversion_op, rotoinversion_op, invalid_op};

  jsonParser &to_json(const SymOpRepresentation::symmetry_type &stype, jsonParser &json) {

    if(stype == SymOpRepresentation::identity_op) {
      json = "identity";
    }
    else if(stype == SymOpRepresentation::mirror_op) {
      json = "mirror";
    }
    else if(stype == SymOpRepresentation::glide_op) {
      json = "glide";
    }
    else if(stype == SymOpRepresentation::rotation_op) {
      json = "rotation";
    }
    else if(stype == SymOpRepresentation::screw_op) {
      json = "screw";
    }
    else if(stype == SymOpRepresentation::inversion_op) {
      json = "inversion";
    }
    else if(stype == SymOpRepresentation::rotoinversion_op) {
      json = "rotoinversion";
    }
    else if(stype == SymOpRepresentation::invalid_op) {
      json = "invalid";
    }
    return json;
  }

  //**********************************************************

  void from_json(SymOpRepresentation::symmetry_type &stype, const jsonParser &json) {
    try {
      std::cout << "Reading in type" << std::endl;
      std::string type = json.get<std::string>();
      if(type == "identity")
        stype = SymOpRepresentation::identity_op;
      else if(type == "mirror")
        stype = SymOpRepresentation::mirror_op;
      else if(type == "glide")
        stype = SymOpRepresentation::glide_op;
      else if(type == "rotation")
        stype = SymOpRepresentation::rotation_op;
      else if(type == "screw")
        stype = SymOpRepresentation::screw_op;
      else if(type == "inversion")
        stype = SymOpRepresentation::inversion_op;
      else if(type == "rotoinversion")
        stype = SymOpRepresentation::rotoinversion_op;
      else if(type == "invalid")
        stype = SymOpRepresentation::invalid_op;
      else {
        std::cout << "Error in 'void from_json(SymOpRepresentation::symmetry_type &stype, const jsonParser &json)'" << std::endl;
        std::cout << "Unrecognized SymOpRepresentation::symmetry_type: '" << type << "'" << std::endl;
        exit(1);
      }
    }
    catch(...) {
      /// re-throw exceptions
      throw;
    }
  }
};

