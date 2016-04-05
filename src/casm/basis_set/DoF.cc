#include "casm/basis_set/DoF.hh"

// needed for to/from json... this should be fixed
#include "casm/crystallography/Molecule.hh"

namespace CASM {

  // ** OccupantDoF **

  /// overload for each template type to be used
  //    because we want to be able to do: void from_json(DoF *dof, const jsonParser &json)
  //

  // int version
  template<> jsonParser &OccupantDoF<int>::to_json(jsonParser &json) const {
    json.put_obj();
    json["DoF_type"] = "OccupantDoF";
    json["DoF_template_type"] = "int";
    json["m_type_name"] = m_type_name;
    json["m_domain"] = m_domain;
    json["m_current_state"] = m_current_state;
    return json;
  }

  jsonParser &to_json(const OccupantDoF<int> &dof, jsonParser &json) {
    return dof.to_json(json);
  }

  void from_json(OccupantDoF<int> &dof, const jsonParser &json) {
    try {
      std::string name = json["m_type_name"].get<std::string>();
      Array<int> domain = json["m_domain"].get<Array<int> >();
      int current_state = json["m_current_state"].get<int>();

      dof = OccupantDoF<int>(name, domain, current_state);
    }
    catch(...) {
      /// re-throw exceptions
      throw;
    }
  }

  // molecule version
  template<> jsonParser &OccupantDoF<Molecule>::to_json(jsonParser &json) const {
    json.put_obj();
    json["DoF_type"] = "OccupantDoF";
    json["DoF_template_type"] = "Molecule";
    json["m_type_name"] = m_type_name;

    json["m_domain"] = m_domain;
    json["m_current_state"] = m_current_state;

    return json;
  }

  jsonParser &to_json(const OccupantDoF<Molecule> &dof, jsonParser &json) {
    return dof.to_json(json);
  }

  // Note: as a hack this expects dof.domain[0] to be present and have the right lattice!!!
  //   it's just used to set the lattice for all the Molecules
  void from_json(OccupantDoF<Molecule> &dof, const jsonParser &json) {
    try {
      std::string name = json["m_type_name"].get<std::string>();

      int current_state = json["m_current_state"].get<int>();

      Array<Molecule> domain;
      //std::cout<<"The size of the dof is "<<dof.get_domain().size()<<std::endl;
      ////std::cout<<"The dof "
      Molecule mol(*(dof[0].home()));
      //domain.clear();
      //std::cout<<"Done initializing molecule"<<std::endl;
      for(int i = 0; i < json["m_domain"].size(); i++) {
        from_json(mol, json["m_domain"][i]);
        domain.push_back(mol);
      }
      //std::cout<<"Done adding the molecules to the list, trying to create"
      //               <<" OccupantDoF<Molecule>"<<std::endl;

      dof = OccupantDoF<Molecule>(name, domain, current_state);
    }
    catch(...) {
      /// re-throw exceptions
      throw;
    }
  }


  // ** ContinuousDoF **

  jsonParser &to_json(const ContinuousDoF &dof, jsonParser &json) {
    return dof.to_json(json);
  }

  void from_json(ContinuousDoF &dof, const jsonParser &json) {
    try {
      dof = ContinuousDoF(json["m_type_name"].get<std::string>(), json["min"].get<double>(), json["max"].get<double>());
    }
    catch(...) {
      /// re-throw exceptions
      throw;
    }
  }


  // ** DoF **

  /// creates jsonParser using polymorphism
  jsonParser &to_json(const DoF *dof, jsonParser &json) {
    return dof->to_json(json);
  }

  /// This allocates a new object to 'dof'.
  ///   It needs a Lattice in case it is a OccupantDoF<Molecule>
  ///
  void from_json(DoF *dof, const jsonParser &json, const Lattice &lat) {
    try {
      if(json["DoF_type"] == "OccupantDoF") {
        if(json["DoF_template_type"] == "int") {

          // prepare a OccupantDoF<int> and then read from json
          OccupantDoF<int> tdof;
          CASM::from_json(tdof, json);

          // copy to dof
          dof = new OccupantDoF<int>(tdof);

        }
        else if(json["DoF_template_type"] == "Molecule") {

          // prepare a OccupantDoF<Molecule> and then read from json
          Array<Molecule> init_domain;
          init_domain.push_back(Molecule(lat));
          OccupantDoF<Molecule> tdof(json["m_type_name"].get<std::string>(), init_domain);

          // this expects dof.domain[0] has a Molecule with the right lattice
          CASM::from_json(tdof, json);

          // copy to dof
          dof = new OccupantDoF<Molecule>(tdof);

        }
        else {
          std::cerr << "Error in 'jsonParser from_json(DoF *dof, const jsonParser &json, const Lattice &lat)'" << std::endl;
          std::cerr << "Unrecognized 'DoF_template_type': '" << json["DoF_template_type"] << "'." << std::endl;
          exit(1);
        }
      }
      else if(json["DoF_type"] == "ContinuousDoF") {

        // prepare a OccupantDoF<Molecule> and then read from json
        ContinuousDoF tdof(json["m_type_name"].get<std::string>(), 0.0, 0.0);
        CASM::from_json(tdof, json);

        // copy to dof
        dof = new ContinuousDoF(tdof);

      }
      else {
        std::cerr << "Error in 'void from_json(DoF *dof, const jsonParser &json, const Lattice &lat)'" << std::endl;
        std::cerr << "Unrecognized 'DoF_type': '" << json["DoF_type"] << "'." << std::endl;
        std::cerr << "Options are: 'OccupantDoF', or 'ContinuousDoF'." << std::endl;
        exit(1);
      }
    }
    catch(...) {
      /// re-throw exceptions
      throw;
    }
  }

}

