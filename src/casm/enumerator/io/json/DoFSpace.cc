#include "casm/casm_io/json/jsonParser.hh"
#include "casm/casm_io/container/json_io.hh"
#include "casm/crystallography/SimpleStructure.hh"
#include "casm/clex/SimpleStructureTools.hh"
#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/enumerator/DoFSpace.hh"

namespace CASM {

  jsonParser &to_json(DoFSpace const &dofspace, jsonParser &json, std::string name) {
    json["dof"] = dofspace.dof_key;
    {
      jsonParser &cjson = json["initial_configuration"];

      cjson["identifier"] = name;

      SimpleStructure sstruc = make_simple_structure(dofspace.config_region.configuration());

      cjson["lattice_vectors"] = sstruc.lat_column_mat.transpose();

      jsonParser &sjson = cjson["sites"];

      for(Index i = 0; i < sstruc.mol_info.size(); ++i) {
        to_json_array(sstruc.mol_info.cart_coord(i), sjson[to_sequential_string(i + 1, sstruc.mol_info.size())][sstruc.mol_info.names[i]]);
      }

      cjson["selected_sites"].put_array();
      for(Index s : dofspace.config_region.sites()) {
        cjson["selected_sites"].push_back(s + 1);
      }

    }

    json["glossary"] = make_axis_glossary(dofspace.dof_key,
                                          dofspace.config_region.configuration(),
                                          dofspace.config_region.sites());

    return json;
  }


}
