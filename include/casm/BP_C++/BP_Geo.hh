
// To Do:
// 1) Test DT_calc_data()
// 2) Test DT_simps_circumcenter()



#ifndef BP_Geo_HH
#define BP_Geo_HH

#include <sstream>
#include <iostream>
#include <fstream>
#include "casm/external/Eigen/Dense"
//#include "casm/BP_C++/Array.hpp"
#include <sys/time.h>
#include <iomanip>
#include "casm/BP_C++/BP_Vec.hh"
#include "casm/BP_C++/BP_GVec.hh"
#include "casm/BP_C++/BP_basic.hh"


namespace BP {
  using namespace std;
  using namespace Eigen;
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  class V_Element;
  class Geo;

  bool vector_is(const Eigen::VectorXd &v, double val, double tol);
  bool vector_is(const Eigen::VectorXd &v, const Eigen::VectorXd &v2, double tol);
  void print_bins(BP_Vec< BP_Vec< BP_Gen_GVec_Member *> > &point_bins, BP_GVec<V_Element> &pts, int curr_facet);


  class V_Element {
  public:

    int dim;				//dimension of all space
    Eigen::VectorXd pos;			//Vector from origin to the vector space (might not be on the face)
    int rank;				//rank of this element's vector space
    //Eigen::MatrixXd Q;				//Unitary complete basis set, first 'rank' columns being vector space of the element
    //Eigen::MatrixXd P;				//Projector onto vector space of the element
    Eigen::VectorXd vertex_mean;
    Eigen::VectorXd out_norm;		// unit outward normal from element, depends on context
    double nvolume;
    int cut;

    V_Element intersect(const V_Element &E1);
    //Eigen::MatrixXd space()
    //{
    //	return Q.leftCols(rank);
    //};
    //Eigen::MatrixXd nullspace()
    //{
    //	return Q.rightCols(dim-rank);
    //};

    V_Element() {
      cut = 0;
      nvolume = 0;
      rank = -1;
    };

    /*bool operator==(const V_Element &E1)
    {
    	//std::cout << "checking rank" << std::endl;
    	// if rank and pos and Q.leftCols(rank) are equivalent, then the V_Element are equivalent
    	if( rank != E1.rank) return 0;
    	//if( (pos).isApprox(E1.pos, 1e-14)) return 1;
    	//else return 0;

    	if( rank == 0) std::cout << "rank are the same, checking pos" << std::endl;
    	bool equiv = 0;
    	double tol = 1e-8;		//this is kind of large, but was finding errors of >1e-14 in pos of line intersections...
    	double d;
    	//std::cout << "   pos: " << pos.transpose() << std::endl;
    	//std::cout << "E1.pos: " << E1.pos.transpose() << std::endl;
    	for( int i=0; i<pos.size(); i++)
    	{
    		d = pos(i) - E1.pos(i);
    		if( rank == 0) std::cout << d << " ";
    		if( d < -tol) return 0;
    		if( d > tol) return 0;
    	}
    	if( rank == 0) std::cout << std::endl;

    	if( rank == 0){ std::cout << " --- are equivalent" << std::endl; return 1;}

    	//std::cout << "pos and rank are the same, but planes may be different" << std::endl;
    	// if pos and rank are the same, we need to check if the planes are the same
    	//  check if null space of E1 has any overlap with space of this
    	for( int j=E1.rank; j<E1.dim; j++)
    	for( int i=0; i<rank; i++)
    	{
    		d = Q.col(i).dot( E1.Q.col(j));
    		if( d < -tol) return 0;
    		if( d > tol) return 0;
    	}
    	return 1;


    };
    */
    void write() {
      std::cout << "rank: " << rank << "  pos: " << pos.transpose() << " cut: " << cut << " vm: " << vertex_mean.transpose() << std::endl;
      //std::cout << "  Q:\n" << Q << std::endl ;
    };

    // rank, pos of element
    V_Element(int i1, Eigen::VectorXd &i2) { //, Eigen::MatrixXd &i3)
      rank = i1;
      pos = i2;
      if(rank == 0) vertex_mean = i2;
      dim = pos.size();
      cut = 0;
      nvolume = 0;


      // if skipping tests
      /*if( 1 == 1) return;

      // Make Q unitary
      FullPivHouseholderQR< Eigen::MatrixXd > tQR;
      tQR.compute(i3.leftCols(rank));
      Q = tQR.matrixQ();

      // Projector onto element space
      //P = Q.leftCols(rank)*Q.leftCols(rank).transpose();

      // check that pos and Q.leftCols(rank) are orhtogonal
      Eigen::MatrixXd m = Q.leftCols(rank).adjoint()*pos;
      bool are_orth = m.isZero();
      if( !are_orth) std::cout << " Possible Error in V_Element() Constructor: pos and the vector space are not orthogonal" << std::endl;
      */
    };

    // rank, pos, vertex_mean of element
    V_Element(int i1, Eigen::VectorXd &i2, Eigen::VectorXd &i3) { //, Eigen::MatrixXd &i3)
      rank = i1;
      pos = i2;
      vertex_mean = i3;
      dim = pos.size();
      cut = 0;
      nvolume = 0;
    };

    // rank, pos, vertex_mean, out_norm, of element
    V_Element(int i1, Eigen::VectorXd &i2, Eigen::VectorXd &i3, Eigen::VectorXd &i4) { //, Eigen::MatrixXd &i5)
      rank = i1;
      pos = i2;
      vertex_mean = i3;
      out_norm = i4;

      dim = pos.size();
      cut = 0;
      nvolume = 0;

      // if skipping tests
      if(1 == 1) return;


      /*
      // Make Q unitary
      FullPivHouseholderQR< Eigen::MatrixXd > tQR;
      tQR.compute(i5.leftCols(rank));
      Q = tQR.matrixQ();

      // Projector onto element space
      // P = Q.leftCols(rank)*Q.leftCols(rank).transpose();

      // check that pos and Q.leftCols(rank) are orhtogonal
      Eigen::MatrixXd m = Q.leftCols(rank).adjoint()*pos;
      bool are_orth = m.isZero();
      if( !are_orth) std::cout << " Possible Error in V_Element() Constructor: pos and the vector space are not orthogonal" << std::endl;
      */
    };

  };


  class CH_data_class {
  public:
    bool is_hull_point;
    BP_Vec<int> closest_facet;
    double dist_to_hull;
    Eigen::VectorXd vec_to_hull;

    CH_data_class() {
    };

    CH_data_class(bool i1, double i2, const Eigen::VectorXd &i3) {
      is_hull_point = i1;
      dist_to_hull = i2;
      vec_to_hull = i3;
    };
  };

  class DT_data_class {
  public:
    Eigen::VectorXd ccenter;

    DT_data_class() {

    };

    DT_data_class(const Eigen::VectorXd &i1) {
      ccenter = i1;
    };
  };

  class Geo {
  private:
    Eigen::VectorXd scale;
    BP_GVec< V_Element> points;
    BP_GVec< bool> connections;
    int dim;
    bool use_bottom;
    Eigen::VectorXd bottom;
    double Geo_tol;									// tolerance for whether points and facets are co-planar
    int verbosity;


    BP_Vec< BP_Vec<int> > equivalent_points;

    bool CH_exists;
    double ch_area;
    double ch_volume;
    double CH_tol;									// tolerance for checking if hull points are equivalent
    BP_Group< V_Element> CH_points;					// points on the hull (not filtered to get 'bottom' hull)
    BP_GVec< V_Element> CH_facets;
    BP_Group< V_Element> CH_filter_points;			// points on the full (including filtering to get 'bottom' hull)
    BP_Group< V_Element> CH_filter_facets;
    BP_Vec< CH_data_class> CH_data;
    bool CH_data_exists;

    bool DT_exists;
    double dt_volume;
    BP_GVec< V_Element> DT_simplexes;
    BP_Vec< DT_data_class> DT_data;
    bool DT_data_exists;

    bool VD_exists;
    double VD_tol;									// tolerance for checking if circumcenters are equivalent
    BP_GVec< V_Element> *VD_elements;
    bool VD_data_exists;

    BP_Graph< V_Element, bool> CH_graph;
    BP_Graph< V_Element, bool> DT_graph;
    BP_Graph< V_Element, bool> VD_graph;


  public:
    Geo();
    Geo(const Geo &geo);
    Geo(const Eigen::MatrixXd &);							// initialize a Geo object with the a matrix of points (each col is a point's position)
    Geo(const Eigen::MatrixXd &m, bool check_for_repeats);	// initialize a Geo object, and tell whether to check for repeated points
    ~Geo();										// Destructor

    void reset_points(const Eigen::MatrixXd &m, bool check_for_repeats);	// reset the Geo object to a new set of points, and tell whether to check for repeated points
    int get_dim();									// the dimension of space
    int size();										// size of points
    Eigen::VectorXd pos(int);								// position of point
    void set_verbosity(int);						// control how much is printed, 10 == print more, other == print less

    BP_Vec< BP_Vec<int> > get_equivalent_points();
    void write_equivalent_points(std::ostream &sout);

    double get_Geo_tol() {
      return Geo_tol;
    };
    double get_CH_tol() {
      return CH_tol;
    };
    double get_VD_tol() {
      return VD_tol;
    };
    void set_Geo_tol(double i1) {
      Geo_tol = i1;
    };
    void set_CH_tol(double i1) {
      CH_tol = i1;
    };
    void set_VD_tol(double i1) {
      VD_tol = i1;
    };

    // Convex Hull functions
    bool	calc_CH();								// find the convex hull of the points
    // if there are non-simplical hull facets, the hull if found correctly, but some hull points may not be found
    // these can be found afterwards by checking CH_dist_to_hull(i) < get_Geo_tol()
    void	CH_bottom(const Eigen::VectorXd &);			// filter the convex hull so that we only consider the 'bottom', which is defined by the vector
    void	CH_wholehull();							// use the whole hull
    double	CH_volume();							// volume of the whole hull
    double	CH_area();								// area of the (filtered) hull
    void    CH_write(std::ostream &);					// write out the facets and vertices of the hull

    // functions to get data on any point
    BP_Vec<int>	CH_closest_facet(int);				// Vector of indices of the closest facets to a point (if filtered, this is the facet along 'bottom' vector)
    double		CH_dist_to_hull(int);				// Distance to 'closest_facet'
    Eigen::VectorXd	CH_vec_to_hull(int);				// Vector to 'closest_facet'
    double		CH_dist_to_hull(const Eigen::VectorXd &)const;		// Distance to 'closest_facet' of given point
    //Eigen::VectorXd	CH_vec_to_hull( const Eigen::VectorXd &);

    // functions to get data on any facets
    int			CH_facets_size();					// # of facets in (filtered) hull
    double		CH_facets_area(int);				// area of this facet
    Eigen::VectorXd    CH_facets_norm(int);				// outward unit normal of this facet
    int			CH_facets_nborverts_size(int);		// # of neighboring vertices to this facet
    BP_Vec<int>	CH_facets_nborverts(int);			// Vector of indices (into CH_filter_points) of neighboring vertices to this facet
    int			CH_facets_nborfacets_size(int);		// # of neighboring facets to this facet
    BP_Vec<int>	CH_facets_nborfacets(int);			// Vector of indices (into CH_filter_facets) of neighboring facets to this facet

    // functions to get data on any point on the hull
    int			CH_verts_size();					// # of vertices on the (filtered) hull
    Eigen::VectorXd	CH_verts_pos(int);					// position of this vertex
    BP_Vec<int>	CH_verts_indices();					// indices of the vertices on the hull (i.e., the columns in original data that correspond to hull points)
    int			CH_verts_nborverts_size(int);		// # of neighboring vertices to this vertex
    BP_Vec<int>	CH_verts_nborverts(int);			// # of
    int			CH_verts_nborfacets_size(int);
    BP_Vec<int>	CH_verts_nborfacets(int);

    // Delaunay Triangulation functions
    bool		calc_DT();
    void		DT_write(std::ostream &);
    double		DT_volume();

    int			DT_simps_size();
    double		DT_simps_volume(int);
    Eigen::VectorXd	DT_simps_circumcenter(int);
    int			DT_simps_nborverts_size(int);
    BP_Vec<int>	DT_simps_nborverts(int);
    int			DT_simps_nborsimps_size(int);
    BP_Vec<int>	DT_simps_nborsimps(int);

    int			DT_verts_size();
    Eigen::VectorXd	DT_verts_pos(int);
    int			DT_verts_nborverts_size(int);
    BP_Vec<int>	DT_verts_nborverts(int);
    int			DT_verts_nborsimps_size(int);
    BP_Vec<int>	DT_verts_nborsimps(int);

    // Voronoi Diagram functions
    bool calc_VD();

    int			VD_elements_size(int);
    BP_Vec<int>	VD_elements_verts(int, int);
    double		VD_elements_volume(int, int);
    Eigen::VectorXd	VD_vertex_pos(int i1);
    Eigen::VectorXd	VD_elements_vertex_mean(int, int);		// dimension of element, index of element

  private:

    void init(const Eigen::MatrixXd &m, bool check_for_repeats);
    void set_points(const Eigen::MatrixXd &m, bool check_for_repeats);

    // CH functions
    bool generate_CH(BP_GVec< V_Element> &, BP_Group< V_Element> &, BP_GVec< V_Element> &, BP_Graph< V_Element, bool> &);
    void form_simplex_from_points(BP_GVec<V_Element> &, BP_Group<V_Element> &, BP_GVec<V_Element> &, BP_Graph<V_Element, bool> &, BP_Group<V_Element> &);
    bool check_hull_point(const Eigen::VectorXd &, BP_GVec<V_Element> &, BP_Graph<V_Element, bool> &, BP_Vec<BP_Gen_Vertex *> &);
    bool check_hull_point_BFS(const Eigen::VectorXd &, BP_GVec<V_Element> &, BP_Graph<V_Element, bool> &, BP_Vec<BP_Gen_Vertex *> &, int);
    BP_Gen_Vertex *form_facet_from_points(BP_Gen_Vertex *, BP_Vec<BP_Gen_Vertex *> &, BP_GVec<V_Element> &, BP_Graph<V_Element, bool> &, const Eigen::VectorXd &);
    void write_hull(int, int, BP_GVec<V_Element> &, BP_Group<V_Element> &);
    bool use_facet(int);
    void CH_calc_distances();

    // DT functions
    void DT_calc_data();

    // VD functions
    void VD_calc_data();
    void VD_face_finder(int, int, BP_Gen_GVec_Member *, BP_Vec<int> &,  BP_Vec<int> &, BP_Vec<int> &);
    BP_Gen_GVec_Member *VD_form_facet_from_ccenters(int, BP_GVec<V_Element> &, BP_Vec<int> &);

    int compare(double, double, double) const;



  };

};

#endif  //  BP_Geo_HH

