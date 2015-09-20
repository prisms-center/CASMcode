#include "casm/hull/Hull.hh"

#include <vector>

#include "casm/hull/GeometryPieces.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/ConfigSelection.hh"
#include "casm/clex/ConfigIterator.hh"


namespace CASM {

  /**
  * Determine which configurations are groundstates and the distance from hull for everyone. We then
  * loop through all the configurations again and add/update these values.
  *
  * This probably means that curr_property should now also include dist_from_hull and is_groundstate, otherwise
  * generic_print will let configurations slide through that don't have these values set, which can mess up
  * your energy file output.
  *
   * Also determine what the hull is and then make a json object that contains all
   * the information you need about it. Information stored as list of facets
   * followed by list of vertexes:
   *
   * facets:
   *    -index
   *    -normal
   *    -area
   *    -intercepts
   *    -neighboring facets
   *    -vertices
   *
   * vertexes:
   *    -index
   *    -corresponding configuration
   *    -on facet
   *    -neighboring vertices
  */

  jsonParser update_hull_props(PrimClex &primclex,
                               ConfigSelectionIterator<false, false> begin,
                               ConfigSelectionIterator<false, false> end,
                               double geo_tol) {

    //Return this
    jsonParser hulljson;
    //We're only going to go through a subset of the configurations that we don't know in advance, since
    //some energies might not be set. We keep track of the configurations that have relaxed_energy
    std::vector<ConfigSelectionIterator<false, false> > valid_config;

    BP::Geo hull;

    //We should have the tolerance set in .casmroot
    //std::cerr << "WARNING in PrimClex::update_hull_props" << std::endl;
    //std::cerr << "Setting hull tolerance to 0.5meV. This might not be what you want!" << std::endl;

    //hull.set_Geo_tol(geo_tol);
    bool hull_found = populate_convex_hull(hull, begin, end, std::back_inserter(valid_config));

    //All configurations have the same size parametric composition vector (right?). +1 for the energy row.
    Index rows = begin->get_param_composition().size() + 1;
    //The columns correspond to the number of configurations;
    Index cols = valid_config.size();

    BP::BP_Vec<int> hull_indices;     //Perhaps we can generalize BP_Vec into a CASM::stable_array?

    if(hull_found) {


      // Clear current hull data (clear from all configurations, not just selected)
      for(auto it = primclex.config_begin(); it != primclex.config_end(); ++it) {
        it->clear_hull_data();
      }

      //Set value to -1 to only get bottom of hull. It's what happens in EnergySet::calc_hull
      Eigen::VectorXd bottom_vector = Eigen::VectorXd::Zero(rows);
      bottom_vector(0) = -1;
      hull.CH_bottom(bottom_vector);
      hull_indices = hull.CH_verts_indices();
      hulljson = hull_data(hull);

      //Go through the configurations and set delta properties for is_groundstate and dist_from_hull
      for(Index i = 0; i < cols; i++) {

        //Determine if configuration i is groundstate or not
        bool is_groundstate = false;
        ulint vertex_col_index;

        if(hull_indices.find_first(i, vertex_col_index)) {
          is_groundstate = true;
          //Add the config name to the jsonParser that has the hull info
          hulljson["vertices"][vertex_col_index]["name"] = valid_config[i].name();

        }

        //Calculate distance from the hull for configuration i
        double dist_from_hull = hull.CH_dist_to_hull(i);

        //set distance from hull and whether it's a groundstate or not
        valid_config[i]->set_hull_data(is_groundstate, dist_from_hull);
      }
    }

    else {
      std::cerr << "ERROR in PrimClex::update_hull" << std::endl
                << "Could not generate the hull for your configurations." << std::endl;
      exit(324);
    }

    return hulljson;
  }


}


