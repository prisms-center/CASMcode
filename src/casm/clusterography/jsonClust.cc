#include "casm/clusterography/jsonClust.hh"

namespace CASM {

  void from_json(SiteClusterJsonHelper clust_helper, const jsonParser &json) {
    SiteCluster &clust = clust_helper.value();
    const Structure &struc = clust_helper.struc();
    UnitCellCoord ucc;
    // can throw
    for(Index i = 0; i < json["sites"].size(); i++) {
      from_json(ucc, json["sites"][i]);
      clust.push_back(struc.get_site(ucc));
    }

  }

  //*******************************************************************************************

  jsonParser &to_json(const ConstSiteClusterJsonHelper &clust_helper, jsonParser &json) {
    const SiteCluster &clust = clust_helper.value();
    const Structure &struc = clust_helper.struc();

    json["min_length"] = clust.min_length(); // <-- not used for input, but useful for users
    json["max_length"] = clust.max_length(); // <-- not used for input, but useful for users
    json["sites"].put_array();
    for(Index i = 0; i < clust.size(); i++) {
      json["sites"].push_back(struc.get_unit_cell_coord(clust[i]));
    }
    return json;
  }

  //*******************************************************************************************
  /// Lean printing of Orbits (only the prototype gets printed; full orbit can later be reconstructed using struc.factor_group();
  jsonParser &to_json(const ConstSiteOrbitJsonHelper &orbit_helper, jsonParser &json) {
    const SiteOrbit &orbit = orbit_helper.value();
    const Structure &struc = orbit_helper.struc();

    json["prototype"] = jsonHelper(orbit.prototype, struc, orbit_helper.tol());
    return json;

  }

  //*******************************************************************************************

  void from_json(SiteOrbitBranchJsonHelper branch_helper, const jsonParser &json) {
    SiteOrbitBranch &branch = branch_helper.value();
    const Structure &struc = branch_helper.struc();

    //can throw

    if(json.contains("pivot")) {
      SiteCluster new_pivot(branch.pivot.home());
      from_json(jsonHelper(new_pivot, struc), json["pivot"]);
      branch.set_pivot(new_pivot);
    }

    branch.reserve(json["orbits"].size());

    for(Index i = 0; i < json["orbits"].size(); i++) {
      SiteCluster tproto(branch.pivot.home());
      from_json(jsonHelper(tproto, struc), json["orbits"][i]["prototype"]);
      tproto.calc_properties();
      branch.push_back(SiteOrbit(tproto));
      branch.back().get_equivalent(struc.factor_group(), branch_helper.tol());
    }


  }

  //*******************************************************************************************

  jsonParser &to_json(const ConstSiteOrbitBranchJsonHelper &branch_helper, jsonParser &json) {
    const SiteOrbitBranch &branch = branch_helper.value();
    const Structure &struc = branch_helper.struc();

    if(branch.pivot.size())
      json["pivot"] = jsonHelper(branch.pivot, struc);

    json["orbits"].put_array();
    for(Index i = 0; i < branch.size(); i++) {
      json["orbits"].push_back(jsonHelper(branch[i], struc, branch_helper.tol()));
    }
    return json;
  }

  //*******************************************************************************************

  void from_json(SiteOrbitreeJsonHelper tree_helper, const jsonParser &json) {
    SiteOrbitree &tree = tree_helper.value();
    const Structure &struc = tree_helper.struc();
    //tree.set_name(json["name"].get<std::string>());
    //tree.set_cspecs(json["cspecs"]);
    tree.set_lattice(json["lattice"].get<Lattice>(), CART);

    tree.resize(json["branches"].size());
    for(Index i = 0; i < json["branches"].size(); i++) {
      from_json(jsonHelper(tree[i], struc, tree.tol()), json["branches"][i]);
    }
    if(json.contains("bspecs")) {
      tree.set_bspecs(json["bspecs"]);
    }
    tree.collect_basis_info(struc);
  }

  //*******************************************************************************************

  jsonParser &to_json(const ConstSiteOrbitreeJsonHelper &tree_helper, jsonParser &json) {
    const SiteOrbitree &tree = tree_helper.value();
    const Structure &struc = tree_helper.struc();

    //json["name"] = tree.name();
    //json["cspecs"] = tree.cspecs();
    json["lattice"] = tree.lattice;
    json["branches"].put_array();

    for(Index i = 0; i < tree.size(); i++) {
      json["branches"].push_back(jsonHelper(tree[i], struc, tree.tol()));
    }
    json["bspecs"] = tree.bspecs();
    return json;
  }

}

