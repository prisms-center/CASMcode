#include "casm/hull/GeometryPieces.hh"

#include "casm/BP_C++/BP_Geo.hh"
#include "casm/BP_C++/BP_Vec.hh"

#include "casm/container/Array.hh"
#include "casm/casm_io/jsonParser.hh"

namespace CASM {

  /**
   * The following is a horrendous piece of code I'm using to avoid #include nightmares
   * that show up trying to define the appropriate to_json/from_json for BP::BP_Vec.
   * The plan right now is to eventually move away from BP::Geo to qHull, so this
   * will disappear on the long run.
   */

  Array<int> BP_to_CASM(const BP::BP_Vec<int> &value) {
    Array<int> casmfied;
    for(Index i = 0; i < value.size(); i++) {
      casmfied.push_back(value[i]);
    }
    return casmfied;
  }

  /**
     * Given a Geo object (presumably a hull for your groundstates) this
     * will return a jsonParser object that has the relevant information
     * to plot 0K phase diagrams in composition and chemical potential space.
     * Includes
     *    rank: equivalent to number of components in your system
     *    facet x:    normal: normal vector to your facet with dimensionality of rank
     *                vertex a: index of one of the corners of facet x
     *                vertex b: ...
     *                ...
     *                intercept 1: intercept of facet x on energy axis at a composition of pure component 1
     *                intercept 2: ...
     *                ...
     *                neighbors: list of neighboring facets to facet x
     *    facet y:    ...
     *    ...
     *
     *    vertex x:   index
     *                position
     *                neighboring vertices
     *                facets it's on
     *
     *
     * The intercepts are related to the chemical potential. Since everything
     * is done relative to one of the components this routine can be changed
     * to account for this. For now, what this routine returns is
     * purely geometrical information.
     */

  jsonParser hull_data(BP::Geo &hull) {
    //accessors in Geo aren't const
    //BP::Geo hull = ghull;

    //There's probably some sort of check that should go here that
    //makes sure your hull object isn't messed up.

    jsonParser hull_info;

    //First figure out the normal and rank
    hull_info["rank"] = hull.get_dim();

    //State tolerance
    hull_info["tolerance"] = hull.get_Geo_tol();

    //How many facets and points?
    hull_info["num_facets"] = hull.CH_facets_size();
    hull_info["num_vertices"] = hull.CH_verts_size();

    //Place to store all facets
    Array<jsonParser> facet_info_list;
    //loop over facets
    for(int f = 0; f < hull.CH_facets_size(); f++) {
      //make sublist of properties
      jsonParser facet_info = facet_data(hull, f);
      facet_info_list.push_back(facet_info);

      //store all the info you just found
      //hull_info[facetname]=facet_info;
    }

    //Place to store all vertices
    Array<jsonParser> vertex_info_list;
    //BP::BP_Vec<int> vertex_indices=hull.CH_verts_indices();
    //loop over vertices
    for(int v = 0; v < hull.CH_verts_size(); v++) {
      //jsonParser vertex_info=vertex_data(hull,vertex_indices[v]);
      jsonParser vertex_info = vertex_data(hull, v);
      vertex_info_list.push_back(vertex_info);
    }

    hull_info["facets"] = facet_info_list;
    hull_info["vertices"] = vertex_info_list;
    return hull_info;
  }

  /**
   * Returns info of just one facet of the hull Geo object specified by index.
   * (See hull data for how the layout is)
   */

  jsonParser facet_data(BP::Geo &hull, int findex) {
    jsonParser facet_info;

    //Store index
    facet_info["findex"] = findex;
    //Get the normal vector
    Eigen::VectorXd normal = hull.CH_facets_norm(findex);
    normal.normalize(); //Not done by default for some reason
    facet_info["normal"] = normal;
    //Get the area
    facet_info["area"] = hull.CH_facets_area(findex);

    //Include vertex info
    facet_info["vertices"] = BP_to_CASM(hull.CH_facets_nborverts(findex));

    //Get indices for neighboring facets FIX THIS
    facet_info["neighboring_facets"] = BP_to_CASM(hull.CH_facets_nborfacets(findex));

    //Figure out intercepts
    //for any point v on a facet with normal n we have n.v=const
    BP::BP_Vec<int> vertex_indices = hull.CH_facets_nborverts(findex);
    Eigen::VectorXd v = hull.CH_verts_pos(vertex_indices[0]); //I can't imagine not having a vertex, but maybe add a check?
    double k = v.dot(normal);
    int rank = hull.get_dim();

    //If the intercept occurs at a point on the facet x, we have x=0,0,1,0,x_E (example) and we want x_E
    double normE = normal(0);
    Array<double> intercepts;

    //first member of vector is the energy, so stop one before that. This loop finds intercepts for all
    //the specified (independent) components
    for(int i = 1; i < rank; i++) {
      double intercept = (k - normal(i)) / normE;
      intercepts.push_back(intercept);
    }

    //lastly get the dependent component intercept
    intercepts.push_back(k / normE);

    facet_info["intercepts"] = intercepts;
    return facet_info;
  }

  /**
   * Returns info of just one vertex of the hull Geo object specified by index.
   * (See hull data for how the layout is)
   */

  jsonParser vertex_data(BP::Geo &hull, int vindex) {
    jsonParser vertex_info;

    //Store index
    vertex_info["vindex"] = vindex;
    //Get position of vertex
    vertex_info["position"] = hull.CH_verts_pos(vindex);

    //Get indices for neighboring vertices
    vertex_info["neighboring_vertices"] = BP_to_CASM(hull.CH_verts_nborverts(vindex));

    //Get indices for neighboring facets
    vertex_info["neighboring_facets"] = BP_to_CASM(hull.CH_verts_nborfacets(vindex));

    //Check to see what facets this vertex is on
    //vertex_info["facets"]=BP_to_CASM(hull.CH_verts_nborfacets(index));
    return vertex_info;
  }
}

