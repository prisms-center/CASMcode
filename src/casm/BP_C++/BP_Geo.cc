
// To Do:
// 1) Test DT_calc_data()
// 2) Test DT_simps_circumcenter()



#ifndef BP_Geo_CC
#define BP_Geo_CC

#include "casm/BP_C++/BP_Geo.hh"

namespace BP {
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  bool vector_is(const Eigen::VectorXd &v, double val, double tol) {
    //double tol = 1e-14;
    double valminus = val - tol;
    double valplus = val + tol;
    for(int i = 0; i < v.size(); i++) {
      if(v(i) < valminus) return 0;
      if(v(i) > valplus) return 0;
    }
    return 1;
  };

  bool vector_is(const Eigen::VectorXd &v, const Eigen::VectorXd &v2, double tol) {
    //double tol = 1e-12;
    double d;
    if(v.size() != v2.size()) return 0;
    for(int i = 0; i < v.size(); i++) {
      d = v(i) - v2(i);
      if(d < -tol) {
        //std::cout <<" d: " << d << " tol: " << tol << std::endl;
        return 0;
      }
      if(d > tol) {
        //std::cout <<" d: " << d << " tol: " << tol << std::endl;
        return 0;
      }
    }
    return 1;
  };

  void print_bins(BP_Vec< BP_Vec< BP_Gen_GVec_Member *> > &point_bins, BP_GVec<V_Element> &pts, int curr_facet) {
    std::cout << "PRINT BINS  curr_facet: " << curr_facet << std::endl;
    int i, j, k;
    for(i = 0; i < point_bins.size(); i++) {
      std::cout << i << " size: " << point_bins[i].size() << std::endl << "  pts: ";
      for(j = 0; j < point_bins[i].size(); j++) {
        std::cout << point_bins[i][j] << " " ;
        for(k = j + 1; k < point_bins[i].size(); k++) {
          if(point_bins[i][j] == point_bins[i][k]) {
            cout << "ERROR in point_bins" << endl;
            //BP_pause();
            exit(1);
          }
        }
      }
      std::cout << std::endl;
    }
  };


  /////////////////////////////////////////////////////////
  // V_Element functions

  // Returns V_Element formed by intersection of this V_Element and E1 (without connectivity)
  /*
  V_Element V_Element::intersect( const V_Element &E1)
  {
  	//std::cout << "begin intersect()" << std::endl;
  	//std::cout << " this->rank: " << rank << "  E1.rank: " << E1.rank << std::endl;
  	//if( rank == 1 & E1.rank == 1)
  	//{
  	//	std::cout << "this->pos: " << pos.transpose() << std::endl;
  	//	std::cout << "this->Q: " << Q.leftCols(rank).transpose() << std::endl;
  	//	std::cout << "E1.pos: " << E1.pos.transpose() << std::endl;
  	//	std::cout << "E1.Q: " << E1.Q.leftCols(E1.rank).transpose() << std::endl;
  	//}
  	int trank, tE_rank;
  	ColPivHouseholderQR<Eigen::MatrixXd> tcQR;
  	FullPivHouseholderQR< Eigen::MatrixXd > tQR;
  	Eigen::MatrixXd tQ1, tQ2;
  	Eigen::VectorXd tpos;

  	/////////////////////////////////////////////////
  	// find Q, rank: vector space of the intersection
  	//	the intersection is the null( union( null(this), null(E1)) )

  	// A is: union( null(this), null(E1) )
  	Eigen::MatrixXd A(dim, (dim-rank) + (dim - E1.rank));
  	A << Q.rightCols(dim-rank), E1.Q.rightCols(dim - E1.rank);

  	// find QR decomposition
  	tQR.compute(A);
  	trank = tQR.rank();

  	//Q is: the null of A (also, we include the range of A to have the complete basis)
  	tE_rank = dim - trank;
  	tQ1 = tQR.matrixQ();
  	tQ2.resize(tQ1.rows(), tQ1.cols());
  	tQ2 << tQ1.rightCols(dim-trank), tQ1.leftCols(trank);


  	/////////////////////////////////////////////////
  	// find tpos: the position of the space
  	//  solve the system:
  	//   space(tE.Q)' *     tpos        = 0
  	//	 nullspace(Q)'    * (tpos - pos)    = 0
  	//   nullspace(E1.Q)' * (tpos - E1.pos) = 0
  	//

  	A.resize(tE_rank+(dim - rank) + (dim - E1.rank),dim);
  	A.block(0,0,tE_rank,dim) = tQ2.block(0,0,dim,tE_rank).adjoint();
  	A.block(tE_rank,0,(dim-rank),dim) = Q.rightCols(dim-rank).adjoint();
  	A.block(tE_rank+(dim-rank),0,dim-E1.rank,dim) = E1.Q.rightCols(dim-E1.rank).adjoint();

  	Eigen::VectorXd b(tE_rank+(dim - rank) + (dim - E1.rank));
  	b.setZero();
  	b.segment(tE_rank,dim-rank) = Q.rightCols(dim-rank).adjoint()*pos;
  	b.segment(tE_rank+dim-rank,dim-E1.rank) = E1.Q.rightCols(dim-E1.rank).adjoint()*E1.pos;

  	// solve Ax=b
  	tcQR.compute(A);			//tQR sometimes fails here :(
  	tpos = tcQR.solve(b);


  	//if( vector_is( tpos, 0))
  	//{
  	//	std::cout << "this->pos: " << pos.transpose() << std::endl;
  	//	std::cout << "E1.pos: " << E1.pos.transpose() << std::endl;
  	//	std::cout << "A:\n" << A << std::endl;
  	//	std::cout << "b:\n" << b.transpose() << std::endl;
  	//	std::cout << "solve: " <<  tQR.solve(b) << std::endl;
  	//	std::cout << "pos:\n" << tpos.transpose() << std::endl;
  	//	std::cout << "tE.Q:\n" << tQ2.leftCols(tE_rank).transpose() << std::endl;
  	//	std::cout << "vector is zero" << std::endl;
  	//	BP_pause();
  	//}

  	//if the rank of the intersection is 0, check if the solution even exists
  	if( tE_rank == 0)
  	{

  		//if( tE.pos.isConstant(0)) std::cout << "isConstant(0)!!!" << std::endl;
  		//else std::cout << "tE.pos: " << tE.pos.transpose() << std::endl;

  		bool a_solution_exists = (A*tpos).isApprox(b, 1e-12);
  		if( !a_solution_exists)
  			std::cout << "The Solution Does not Exist! (tol=1e-12)" << std::endl;
  		//else
  		//	std::cout << "The Solution Exists!" << std::endl;
  		//tE.vertex_mean = tE.pos;

  	}

  	V_Element tE( tE_rank, tpos, tQ2);


  	//std::cout << "finish intersect()" << std::endl;
  	return tE;
  };
  */

  /////////////////////////////////////////////////////////
  // Geo functions

  Geo::Geo() {
    CH_exists = 0;
    DT_exists = 0;
    VD_exists = 0;
    ch_area = -1;
    ch_volume = -1;
    CH_data_exists = 0;

    DT_data_exists = 0;
    dt_volume = -1;

    use_bottom = 0;
    CH_tol = 1e-14;
    VD_tol = 1e-14;
    Geo_tol = 1e-14;
    verbosity = 10;
  };

  Geo::Geo(const Geo &geo) {
    CH_exists = 0;
    DT_exists = 0;
    VD_exists = 0;
    ch_area = -1;
    ch_volume = -1;
    CH_data_exists = 0;

    DT_data_exists = 0;
    dt_volume = -1;

    use_bottom = 0;
    CH_tol = 1e-14;
    VD_tol = 1e-14;
    Geo_tol = 1e-14;
    verbosity = 10;
  };

  Geo::Geo(const Eigen::MatrixXd &m) {
    init(m, true);
  };

  /// use this with Geo(m,false) to not check for repeats
  Geo::Geo(const Eigen::MatrixXd &m, bool i1) {
    init(m, i1);
  };

  void Geo::init(const Eigen::MatrixXd &m, bool check_for_repeats) {
    int i, j, k;

    CH_exists = 0;
    DT_exists = 0;
    VD_exists = 0;
    ch_area = -1;
    ch_volume = -1;
    CH_data_exists = 0;

    DT_data_exists = 0;
    dt_volume = -1;

    dim = m.rows();
    use_bottom = 0;
    CH_tol = 1e-14;
    VD_tol = 1e-14;
    Geo_tol = 1e-14;
    verbosity = 10;

    set_points(m, check_for_repeats);
  }

  void Geo::set_points(const MatrixXd &m, bool check_for_repeats) {
    int i, j, k;

    //scale points to range 1
    Eigen::VectorXd range_min, range_max;
    range_min.resize(dim);
    range_min.setZero();
    range_max.resize(dim);
    range_max.setZero();

    for(i = 0; i < m.cols(); i++)
      for(j = 0; j < m.rows(); j++) {
        if(i == 0) {
          range_min(j) = m(j, i);
          range_max(j) = m(j, i);
        }
        else {
          if(m(j, i) < range_min(j))
            range_min(j) = m(j, i);
          if(m(j, i) > range_max(j))
            range_max(j) = m(j, i);
        }
      }

    scale = range_max - range_min;

    // scaling messes up DT and VD, temporary solution:
    scale.setOnes();
    //cout << "init scale: " << scale << endl;



    if(check_for_repeats) {
      //////////////////////////////////////////////////
      // check for repeat points

      Eigen::MatrixXd tmp_m = m;
      Eigen::VectorXd mean_point;
      mean_point = m.rowwise().mean();

      BP_Vec< int> equiv_bin;
      //BP_Vec< BP_Vec<int> > equivalent_points;
      for(i = 0; i < m.cols(); i++) {
        equiv_bin.add(i);
        equivalent_points.add();
        equivalent_points.last().add(i);
      }

      bool is_same;
      for(i = 0; i < m.cols(); i++) {
        for(j = i + 1; j < m.cols(); j++) {
          if(equiv_bin[i] == equiv_bin[j])
            continue;

          is_same = true;
          for(k = 0; k < m.rows(); k++) {
            if(std::fabs(m(k, i) - m(k, j)) / scale(k) > CH_tol) {
              is_same = false;
              break;
            }

          }

          if(is_same) {
            int bin1 = equiv_bin[i];
            int bin2 =  equiv_bin[j];
            //std::cout << "is same " << i << " " << j << " bins: " << bin1 << " " << bin2 << std::endl;
            for(k = 0; k < equivalent_points[bin2].size(); k++) {
              equiv_bin[ equivalent_points[bin2][k]] = bin1;

            }

            equivalent_points[bin1].append(equivalent_points[bin2]);
            equivalent_points[bin2].clear();
            //std::cout << "  done" << std::endl;
          }
          else {
            //std::cout << "different " << i << " " << j << std::endl;
          }
        }

      }


      int inc = 1;
      for(i = 0; i < equivalent_points.size(); i++) {
        if(equivalent_points[i].size() > 1) {
          //std::cout << "Warning, equivalent point indices: ";
          //for( j=0; j<equivalent_points[i].size(); j++)
          //	std::cout << equivalent_points[i][j] << " ";
          //std::cout << std::endl;
          //std::cout << "   point: " ;
          //for( k=0; k<m.rows(); k++)
          //	std::cout << " " << m(k, equivalent_points[i][0]);
          //std::cout << std::endl;

          for(j = 1; j < equivalent_points[i].size(); j++) {
            for(k = 0; k < m.rows(); k++)
              tmp_m(k, equivalent_points[i][j]) = mean_point(k) + (inc * 2.0) * CH_tol;
            inc++;
          }
        }
      }

      //////////////////////////////////////////////////

      Eigen::VectorXd v;
      //Eigen::MatrixXd Midentity = Eigen::MatrixXd::Identity(dim,dim);
      points.capacity(tmp_m.cols());
      for(i = 0; i < tmp_m.cols(); i++) {
        v = tmp_m.col(i).cwiseQuotient(scale);
        //points.add( V_Element( 0, v, Midentity) );
        points.add(V_Element(0, v));
      }
    }
    else {
      //////////////////////////////////////////////////

      Eigen::VectorXd v;
      //Eigen::MatrixXd Midentity = Eigen::MatrixXd::Identity(dim,dim);
      points.capacity(m.cols());
      for(i = 0; i < m.cols(); i++) {
        v = m.col(i).cwiseQuotient(scale);
        //points.add( V_Element( 0, v, Midentity) );
        points.add(V_Element(0, v));
      }
    }



  };



  BP_Vec< BP_Vec<int> > Geo::get_equivalent_points() {
    return equivalent_points;
  };

  void Geo::write_equivalent_points(std::ostream &sout) {

    int i, j, k;
    Eigen::VectorXd v;
    for(i = 0; i < equivalent_points.size(); i++) {
      if(equivalent_points[i].size() > 1) {
        std::cout << "Warning, equivalent point indices: ";
        for(j = 0; j < equivalent_points[i].size(); j++)
          std::cout << equivalent_points[i][j] << " ";

        v = pos(equivalent_points[i][0]);

        std::cout << "   point: " ;
        for(k = 0; k < v.rows(); k++)
          std::cout << " " << v(k);
        std::cout << std::endl;
      }
    }
  };


  Geo::~Geo() {
    if(VD_exists) delete [] VD_elements;
  }

  int Geo::get_dim() {
    return dim;
  };

  void Geo::reset_points(const Eigen::MatrixXd &m, bool check_for_repeats) {
    int i, j, k;

    CH_graph.clear();
    CH_graph.free();
    DT_graph.clear();
    DT_graph.free();
    VD_graph.clear();
    VD_graph.free();

    CH_facets.clear();
    CH_facets.free();
    CH_data.clear();
    CH_data.free();
    DT_simplexes.clear();
    DT_simplexes.free();
    DT_data.clear();
    DT_data.free();
    if(VD_exists) delete [] VD_elements;

    VD_exists = 0;
    CH_exists = 0;
    DT_exists = 0;
    ch_area = -1;
    ch_volume = -1;
    CH_data_exists = 0;

    DT_data_exists = 0;
    dt_volume = -1;

    dim = m.rows();

    points.clear();
    equivalent_points.clear();
    connections.clear();

    set_points(m, check_for_repeats);
  };

  int Geo::size() {
    return points.size();
  };

  Eigen::VectorXd Geo::pos(int i1) {
    return points[i1].pos.cwiseProduct(scale);
  };

  void Geo::set_verbosity(int i1) {
    verbosity = i1;
  };

  /////////////////////////////////////////////////////////
  // Convex hull functions

  // return 1 if calculated ok, returns 0 if not calculated
  //   currently, points on the hull may not be listed as hull points if they are part of a non-simplical hull facet
  //	 these can be identified afterwards by checking for any CH_dist_to_hull(i) < get_Geo_tol()
  bool Geo::calc_CH() {
    if(CH_exists) return 1;

    CH_facets.capacity(1000);		//CONSTANT
    CH_facets.cap_increment(1000);		//CONSTANT
    connections.cap_increment(10000);	//CONSTANT
    //connections.capacity(10000);
    CH_data.clear();
    ch_area = -1;
    ch_volume = -1;
    CH_data_exists = 0;

    //std::cout << "points: " << std::endl;
    //for( int i=0; i<points.size(); i++)
    //{
    //	std::cout << pos(i).transpose() << std::endl;
    //}

    if(generate_CH(points, CH_points, CH_facets, CH_graph)) {
      CH_exists = 1;

      CH_filter_points.clear();
      for(int i = 0; i < CH_points.size(); i++)
        CH_filter_points.add(CH_points.member(i));

      CH_filter_facets.clear();
      for(int i = 0; i < CH_facets.size(); i++)
        CH_filter_facets.add(CH_facets.member(i));
    }
    return CH_exists;

  };


  // find convex hull in space of dimension i1 with matrix(dim, #points) of points
  // return 1 if found, 0 if not found
  //   currently, points on the hull may not be listed as hull points if they are part of a non-simplical hull facet
  //	 these can be identified afterwards by checking for any CH_dist_to_hull(i) < get_Geo_tol()
  bool Geo::generate_CH(BP_GVec< V_Element> &pts, BP_Group< V_Element> &hull_pts, BP_GVec< V_Element> &hull_fcts, BP_Graph< V_Element, bool> &graph) {
    // based on Quickhull algorithm:
    // Barber, C.B., Dobkin, D.P., and Huhdanpaa, H.T.,
    // "The Quickhull algorithm for convex hulls,"
    // ACM Trans. on Mathematical Software, 22(4):469-483, Dec 1996,
    //    http://www.qhull.org

    // The generate_CH() function works by creating a starting simplical hull,
    //	 binning each point to a facet it is 'outside' of,
    //	 then building the hull outwards by checking if each point is outside of the current hull.
    //   Any facets that can 'see' the point being checked are deleted,
    //   and new hull facets created from the new point and the 'ridge' of edges surrounding the
    //   facets that can see the new point.
    //   When a facet is deleted the points binned with it are re-binned with the new facets.
    //


    //std::cout << "begin generate_CH()" << std::endl;
    int i, j, k, ii, jj, kk, c, max_j;
    int dim = pts[(unsigned long int) 0].pos.rows();
    int curr_facet = 0;
    timeval tim;
    double start_time;
    double curr_time;
    double walltime;
    int lasttime = 1;
    gettimeofday(&tim, NULL);
    start_time = tim.tv_sec + (tim.tv_usec / 1000000.0);

    {
      // initialize variables
      Eigen::FullPivHouseholderQR< Eigen::MatrixXd > tQR;
      Eigen::VectorXi checked = Eigen::VectorXi::Zero(points.size());
      Eigen::MatrixXd start(dim, dim + 1);
      Eigen::MatrixXd startQ(dim, dim);
      Eigen::MatrixXd Midentity = Eigen::MatrixXd::Identity(dim, dim);
      Eigen::VectorXd v1, v2;//, order;
      Eigen::VectorXd v, center;
      double d, max_d;
      bool cont = 1;
      //BP_Vec< BP_Gen_GVec_Member*> start_list;
      BP_Group< V_Element> start_list;
      BP_Vec<BP_Gen_Vertex *> interior;
      BP_Vec< BP_Vec<BP_Gen_Vertex *> > perimeter;
      BP_Vec<BP_Gen_Vertex *> neighbors;
      BP_Vec<BP_Gen_Vertex *> new_elements;
      BP_Vec<BP_Gen_Vertex *> erase_list;
      BP_Vec< BP_Vec< BP_Gen_Vertex *> > tverts;
      BP_Vec< BP_Vec<int> > point_bins;
      BP_Gen_Vertex *n1, *n2, *n3, *p1, *p2, *cell;
      bool is_neighbor, not_sharing, found_vertex, add_point, placed;

      hull_pts.clear();
      hull_fcts.clear();
      graph.clear();
      hull_pts.capacity(pts.size());
      hull_fcts.capacity(10000);
      hull_fcts.cap_increment(10000);
      //connections.capacity(10000);
      connections.cap_increment(10000);
      graph.capacity(10000);
      graph.cap_increment(10000);




      //std::cout << "points: " << points.size() << std::endl;
      //for( i=0;i<pts.size(); i++)
      //{
      //	std::cout << "pos: " << pts[i].pos.transpose() << std::endl;
      //}

      // find dim+1 non-planar points, save in 'start_list'
      //std::cout << "find start_list" << std::endl;
      {
        start.col(0) = pts[(unsigned long int) 0].pos;
        start_list.add(pts.member((unsigned long int) 0));
        checked(0) = 1;
        i = 0;
        j = 1;
        tQR.setThreshold(1e-12);
        do {
          //std::cout << "i: " << i << std::endl;
          v1 = pts[j].pos - pts[(unsigned long int) 0].pos;
          startQ.col(i) = v1;
          tQR.compute(startQ.leftCols(i + 1));
          if(tQR.rank() == i + 1) {
            //cout << "  add: " << pts[j].pos.transpose() << std::endl;
            start.col(i + 1) = pts[j].pos;
            start_list.add(pts.member(j));
            checked(j) = 1;
            i++;
          }

          j++;
          if((i < dim) && (j == pts.size())) {
            std::cerr << "Warning in generate_CH(): range of points is less than dimension of space" << std::endl;
            //for( k=0; k<pts.size(); k++)
            //	std::cout << "pt: " << k << " pos: " << pts[k].pos.transpose() << std::endl;
            std::cerr << "Returning without finding hull" << std::endl;
            return 0;
          }

        }
        while(i < dim);
      }

      //std::cout << "start:\n" << start << std::endl;
      //std::cout << "startQ:\n" << startQ << std::endl;
      //std::cout << "  R:\n" << tQR.matrixQR() << std::endl;

      // form simplex from 'start' points
      {
        form_simplex_from_points(pts, hull_pts, hull_fcts, graph, start_list);
      }

      // find mean of 'start' points, this is inside the hull
      center.resize(dim);
      center.setZero();
      for(i = 0; i < start_list.size(); i++)
        //	center += pts[start_list[i]].pos;
        center += start_list[i].pos;
      center /= start_list.size();

      // create bins for all the points
      point_bins.capacity(pts.size());
      for(i = 0; i < hull_fcts.size(); i++) {
        point_bins.add();
        ii = pts.size() / hull_fcts.size();
        if(ii < 1000) ii = 1000;
        jj = ii / 10;
        point_bins[i].capacity(ii);
        point_bins[i].cap_increment(jj);
      }

      int points_left = 0;

      // bin the points according to the facets they are 'outside' of
      for(i = 0; i < pts.size(); i++)
        if(checked(i) == 0) {
          points_left++;
          for(j = 0; j < hull_fcts.size(); j++) {
            d = hull_fcts[j].out_norm.dot(pts[i].pos - hull_fcts[j].vertex_mean);
            //std::cout << "d: " << d << std::endl;
            c = compare(0, d, Geo_tol);
            // if pts[i] is 'on' or 'outside' hull_fcts[j]
            if(c != -1) {
              // add pts[i] to point_bins[j] & break;
              point_bins[j].add(i);
              break;
            }
          }
        }

      // find the first facet that has points in its bin
      bool curr_facet_ok = 0;
      while(!curr_facet_ok) {
        if(curr_facet < hull_fcts.size()) {
          if(point_bins[curr_facet].size() > 0)
            curr_facet_ok = 1;
          else
            curr_facet++;
        }
        else
          curr_facet_ok = 1;
      }

      //print_bins(point_bins, pts, curr_facet);

      if(verbosity >= 10) {
        std::cout << "Starting hull*********************************" << std::endl;
        std::cout << "#point: " << hull_pts.size() << std::endl;
        std::cout << "#facets: " << hull_fcts.size() << std::endl;
        std::cout << std::endl << "BEGIN LOOP OVER POINTS &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << std::endl;
      }

      // loop over all remaining points
      int check_count = 0;
      while(curr_facet < hull_fcts.size()) {
        // pick a trial point from the curr_facet's bin
        //  make the trial point the farthest point from the center
        max_j = 0;
        max_d = (pts[ point_bins[curr_facet][0]].pos - center).norm();
        kk = point_bins[curr_facet].size();
        for(j = 1; j < kk; j++) {
          d = (pts[ point_bins[curr_facet][j]].pos - center).norm();
          if(d > max_d) {
            max_d = d;
            max_j = j;
          }
        }

        // set i to index of trial point
        i = point_bins[curr_facet][max_j]; //pts.get_index( point_bins[curr_facet][max_j]);
        point_bins[curr_facet].remove(max_j);
        points_left--;

        //for(i=0;i<pts.size();i++)
        {
          //if( i % 1 == 0)
          //{
          //	std::cout << "i: " << i << std::endl;
          //	std::cout << "#point: " << hull_pts.size() << std::endl;
          //	std::cout << "#facets: " << hull_fcts.size() << std::endl;
          //
          //}

          //write_hull(check_count, i, pts, hull_pts);
          check_count++;

          gettimeofday(&tim, NULL);
          curr_time = tim.tv_sec + (tim.tv_usec / 1000000.0);
          walltime = curr_time - start_time;

          if(verbosity >= 10)
            if(walltime >= 5.0 * lasttime) {
              std::cout << "Runtime(s): ";
              std::cout << walltime;
              std::cout << "  points remaining: " << points_left << "  #facets: " << hull_fcts.size() << std::endl;
              lasttime++;

            }

          if(checked(i) == 0) {
            interior.clear();
            perimeter.clear();
            new_elements.clear();

            if(verbosity >= 20) std::cout << std::endl << "check point: " << pts[i].pos.transpose() << " unscaled: " << pts[i].pos.cwiseProduct(scale).transpose() << std::endl;

            // identify 'interior' facets that 'see' the new point
            add_point = check_hull_point_BFS(pts[i].pos, hull_fcts, graph, interior, curr_facet);
            //add_point = check_hull_point(pts[i].pos, hull_fcts, graph, interior);
            checked(i) = 1;
            if(add_point == 0) {
              //if( verbosity >= 20)
              std::cout << "-- do not add point" << std::endl;
              continue;
            }

            if(verbosity >= 20)
              std::cout << "++ add point" << std::endl;

            // add new point

            n1 = graph.add_vertex(hull_pts.add(pts.member(i)));

            //   'interior' facets will be 'cut' and new facets formed from the 'perimeter' of all
            //   the interior facets to the new point

            // for each 'interior' element, get list of vertices
            if(verbosity >= 20) std::cout << "get lists of vertices" << std::endl;
            tverts.clear();
            for(j = 0; j < interior.size(); j++) {
              tverts.add();
              tverts[j].clear();
              //std::cout << "tverts.size(): " << tverts.size() << std::endl;
              for(k = 0; k < graph.num_incident_edges(interior[j]); k++) {
                n2 = graph.get_neighbor(interior[j], k);
                if(graph.vert_val(n2).rank == 0) {
                  tverts[j].add(n2);
                  //std::cout << "  tverts[" << j << "].size(): " << tverts[j].size() << std::endl;
                }
              }

              //std::cout << "tverts[" << j << "].size(): " << tverts[j].size() << std::endl;
            }

            //std::cout << "output tverts: " << std::endl;
            //std::cout << "tverts.size(): " << tverts.size() << std::endl;
            //std::cout << "tverts[0].size(): " << tverts[0].size() << std::endl;
            //for( j=0; j<tverts.size(); j++)
            //{
            //	std::cout << "facet " << j << ", vertices: " ;
            //	for( k=0; k<tverts[j].size(); k++)
            //		std::cout << vert_to_index(tverts[j][k]) << " " ;
            //	std::cout << std::endl;
            //}

            // find perimeter segments
            if(verbosity >= 20) std::cout << "find perimeter segments" << std::endl;
            neighbors.clear();
            perimeter.clear();
            // loop over interior elements
            ///std::cout << "#interior: " << interior.size() << std::endl;
            for(j = 0; j < interior.size(); j++) {
              //std::cout << "   j: " << j << " nie: " << graph.num_incident_edges(interior[j])<< std::endl;
              // loop over their neighboring facets
              for(k = 0; k < graph.num_incident_edges(interior[j]); k++) {
                //std::cout << " k: " << k << std::endl;
                n2 = graph.get_neighbor(interior[j], k);
                //std::cout << " n2: " << vert_to_index(n2) << std::endl;
                // if facet, and not being cut
                if(graph.vert_val(n2).rank != 0 & graph.vert_val(n2).cut == 0) {
                  //std::cout << " add" << std::endl;
                  // we will make new element connected to this one
                  perimeter.add();
                  perimeter[perimeter.size() - 1].clear();
                  neighbors.add(n2);

                  // n2 is neighbor facet to interior[j]

                  // find shared vertices
                  //   loop over interior vertices and include all in both
                  //for( ii=0;ii<tverts[j].size();ii++)
                  //{
                  //	//std::cout << "  ii: " << ii << std::endl;
                  //	if( graph.are_1NN( tverts[j][ii], n2)) perimeter[perimeter.size()-1].add(tverts[j][ii]);
                  //	else
                  //	{	// only 1 unshared vertex for neighboring simplexes so add the rest
                  //		for( kk=ii+1;kk<tverts[j].size();kk++)
                  //			perimeter[perimeter.size()-1].add(tverts[j][kk]);
                  //		break;
                  //	}
                  //}


                  // mark all vertices of interior facet
                  for(ii = 0; ii < tverts[j].size(); ii++)
                    graph.vert_val(tverts[j][ii]).cut = 1;

                  int nborcount = 0;
                  // check vertices of neighbor facet
                  for(ii = 0; ii < graph.num_incident_edges(n2); ii++) {
                    n3 = graph.get_neighbor(n2, ii);
                    if(graph.vert_val(n3).rank == 0) {
                      nborcount++;
                      if(graph.vert_val(n3).cut == 1)	// if also neighbor with interior facet
                        perimeter[perimeter.size() - 1].add(n3);	// add to perimeter
                    }
                  }

                  if(nborcount != dim) {
                    std::cout << "nborcount = " << nborcount << std::endl;
                    std::cout << "error, nborcount != dim.  Maybe try increasing tol.  Geo_tol: " << Geo_tol << std::endl;
                    abort();
                    //BP_pause();
                  }

                  // unmark vertices of interior facet
                  for(ii = 0; ii < tverts[j].size(); ii++)
                    graph.vert_val(tverts[j][ii]).cut = 0;


                }
              }
            }

            //std::cout << " cc: " << check_count << " perimeter.size(): " << perimeter.size() << std::endl;
            //for( j=0; j<perimeter.size(); j++)
            //{
            //	std::cout << " j: " << j << " size: " << perimeter[j].size() << std::endl;
            //	if( perimeter[j].size() != dim-1)
            //		BP_pause();
            //	std::cout << "segment " << j << ", vertices: " ;
            //	for( k=0; k<perimeter[j].size(); k++)
            //		std::cout << vert_to_index(perimeter[j][k]) << " " ;
            //	std::cout << std::endl;
            //}


            // add new_elements for each of the perimeter segments, and add connections
            if(verbosity >= 20) std::cout << " add new elements" << std::endl;
            new_elements.clear();
            for(j = 0; j < perimeter.size(); j++) {
              // add new facet and connect with vertices
              //  returns 0 if rank of facet < dim-1 (vertices in a line)
              //std::cout << " add facet" << std::endl;
              n2 = form_facet_from_points(n1, perimeter[j], hull_fcts, graph, center);

              // connect to neighbor
              if(n2 != 0) {
                new_elements.add(n2);
                k = point_bins.size();
                point_bins.add();
                ii = point_bins[curr_facet].size() / interior.size();
                if(ii < 1000) ii = 1000;
                jj = ii / 10;
                //std::cout << " cap: " << ii << " incr: " << jj << std::endl;

                point_bins[k].capacity(ii);
                point_bins[k].cap_increment(jj);


                //std::cout << "connect to neighbor" << std::endl;
                graph.connect_vertices(n2, neighbors[j], connections.add(0), 0);
                //std::cout << "#Connections: " << edge_size() << std::endl;

              }


            }

            if(verbosity >= 20) std::cout << "try connect new elements  new_elements.size(): " << new_elements.size() << std::endl;
            // connect new_elements which share all but one vertex
            for(j = 0; j < new_elements.size(); j++) {
              if(new_elements[j] == 0) continue;
              // don't double count
              for(k = j + 1; k < new_elements.size(); k++) {
                if(new_elements[k] == 0) continue;
                // connect if sharing all but one vertex
                //if( verbosity >= 20) std::cout << " check " << j << " " << k << " perimeter[j].size(): " << perimeter.size() << "  #fcts: " << hull_fcts.size() << std::endl;
                not_sharing = 0;
                is_neighbor = 1;
                for(jj = 0; jj < perimeter[j].size(); jj++) {
                  // check if vertex perimeter[j][jj] is part of perimeter[k]
                  found_vertex = 0;
                  for(kk = 0; kk < perimeter[k].size(); kk++) {
                    if(perimeter[j][jj] == perimeter[k][kk]) {
                      found_vertex = 1;
                      break;
                    }

                  }

                  // if the vertex is not shared
                  if(found_vertex == 0) {
                    if(not_sharing == 0) not_sharing = 1;	// if not sharing 1, note
                    else if(not_sharing == 1) {				// if not sharing >1, these are not neighboring facets
                      is_neighbor = 0;
                      jj = perimeter[j].size();
                      kk = perimeter[k].size();
                    }
                  }

                }

                if(is_neighbor == 1) {

                  for(jj = 0; jj < graph.num_incident_edges(new_elements[j]); jj++) {
                    if(graph.get_neighbor(new_elements[j], jj) == new_elements[k]) {
                      std::cout << "double connection! j: " << j << " k: " << k << std::endl;
                      //BP_pause();
                      exit(1);
                    }

                  }

                  //if( verbosity >= 20) std::cout << "connect new elements" << std::endl;
                  graph.connect_vertices(new_elements[j], new_elements[k], connections.add(0), 0);
                  //if( verbosity >= 20) std::cout << "#Connections: " << graph.edge_size() << std::endl;

                }
              }

            }



            if(verbosity >= 20) std::cout << "erase interior elements: " << interior.size() << std::endl;
            //print_bins(point_bins, pts, curr_facet);

            // erase the interior elements that 'see' the new point
            while(interior.size() > 0) {
              erase_list.clear();
              //std::cout << "interior.size(): " << interior.size() << std::endl;
              for(k = 0; k < graph.num_incident_edges(interior[0]); k++) {
                //std::cout << " k: " << k << std::endl;
                n2 = graph.get_neighbor(interior[0], k);
                if(graph.vert_val(n2).rank == 0) {
                  //std::cout << "  is vert" << std::endl;
                  if(graph.num_incident_edges(n2) == 1) {
                    //std::cout << "   add erase vertex" << std::endl;
                    erase_list.add(n2);
                  }
                  else {
                    //std::cout << "  nie: " << num_incident_edges(n2) << std::endl;
                  }
                }
              }

              while(erase_list.size() > 0) {
                //std::cout << "erase point" << std::endl;
                graph.remove_vertex(erase_list[0]);
                hull_pts.remove(graph.vert_to_member_vert(erase_list[0]));
                erase_list.remove(0);
              }

              //std::cout << "  - erase facet" << std::endl;
              // before removing facet, must re-apportion that facet's outside points among the new_elements
              j = hull_fcts.get_index(graph.vert_to_member_vert(interior[0]));
              kk = point_bins[j].size();
              double tmp_min = 1;
              for(k = 0; k < kk; k++) {
                placed = 0;
                for(ii = 0; ii < new_elements.size(); ii++) {
                  d = graph.vert_val(new_elements[ii]).out_norm.dot(pts[point_bins[j][k]].pos - graph.vert_val(new_elements[ii]).vertex_mean);
                  //std::cout << "d: " << d << std::endl;
                  c = compare(0, d, Geo_tol);

                  //if( c == 0)
                  //{
                  //	std::cout << "c== 0" << std::endl;
                  //	std::cout << " on== " << graph.vert_val(new_elements[ii]).out_norm.transpose() << std::endl;
                  //	std::cout << "pos== " << pts[point_bins[j][k]].pos.transpose() << std::endl;
                  //	std::cout << " mv== " << graph.vert_val(new_elements[ii]).vertex_mean.transpose() << std::endl;
                  //	std::cout << " new_elements nbors: " << std::endl;
                  //	for( int j1=0; j1<graph.num_incident_edges( new_elements[ii]); j1++)
                  //	{
                  //		if( graph.vert_val( graph.get_neighbor( new_elements[ii], j1)).rank == 0)
                  //		{
                  //			std::cout << " nbr== " << graph.vert_val( graph.get_neighbor( new_elements[ii], j1)).pos.transpose() << std::endl;
                  //		}
                  //	}
                  //
                  //	//BP_pause();
                  //}

                  //std::cout << "c: " << c << std::endl;
                  // if point_bins[j][k] is 'on' or 'outside' hull_fcts[j]
                  if(c != -1) {
                    jj = hull_fcts.get_index(graph.vert_to_member_vert(new_elements[ii]));
                    // add pts[i] to point_bins[j] & break;
                    point_bins[jj].add(point_bins[j][k]);
                    placed = 1;
                    break;
                  }
                  else {
                    //std::cout << " point inside hull by distance, d: " << d << std::endl;
                    if(std::fabs(d) > tmp_min)
                      tmp_min = d;
                  }

                }

                if(placed == 0) {
                  //std::cout << " ^^^ point skipped: tmp_min: " << tmp_min << std::endl;
                  //BP_pause();
                  points_left--;
                }
              }

              if(j < curr_facet) {
                curr_facet = j;
                //std::cout << " j < curr_facet" << std::endl;
                //BP_pause();
              }

              point_bins[j].clear();
              point_bins[j].free();
              point_bins.remove(j);
              graph.erase(interior[0]);
              interior.remove(0);
              //print_bins(point_bins, pts, curr_facet);

            }

            kk = new_elements.size();
            for(i = 0; i < kk; i++) {
              jj = hull_fcts.get_index(graph.vert_to_member_vert(new_elements[i]));
              point_bins[jj].free();

            }

            curr_facet_ok = 0;

            while(!curr_facet_ok) {
              if(curr_facet < hull_fcts.size()) {
                if(point_bins[curr_facet].size() > 0)
                  curr_facet_ok = 1;
                else {
                  //std::cout << "N_alloc: " << point_bins[curr_facet].size_alloc() << std::endl;
                  curr_facet++;
                }
              }
              else
                curr_facet_ok = 1;
            }

            //std::cout << "hull_pts: " << std::endl;
            //for(j=0;j<hull_pts.size();j++)
            //{
            //	std::cout << "j: " << j << " pos: " << hull_pts[j].pos.transpose() << std::endl;
            //}


            //std::cout << "check" << std::endl << std::endl;
            for(j = 0; j < hull_fcts.size(); j++) {
              //	std::cout << " j: " << j << " nie: " << graph.num_incident_edges(hull_fcts.member(j)) << std::endl;
              //
              if(graph.num_incident_edges(hull_fcts.member(j)) != dim * 2) {
                std::cout << " hull connection error!" << std::endl;
                std::cout << "  nie: " << graph.num_incident_edges(hull_fcts.member(j)) << std::endl;
                for(k = 0; k < graph.num_incident_edges(hull_fcts.member(j)); k++)
                  std::cout << "   k: " << k << " rank: " << graph.vert_val(graph.get_neighbor(hull_fcts.member(j), k)).rank << std::endl;
                std::cout << "  Geo_tol: " << get_Geo_tol() << std::endl;
                //BP_pause();
                exit(1);
              }
            }


          }
        }
      }

      //write_hull(check_count, pts.size(), pts, hull_pts);

      //calc_volume(dim, CH_elements);


      if(verbosity >= 10) {
        std::cout << std::endl << "FINAL HULL *****************" << std::endl;
        std::cout << "#point: " << hull_pts.size() << std::endl;
        std::cout << "#facets: " << hull_fcts.size() << std::endl;

        gettimeofday(&tim, NULL);
        curr_time = tim.tv_sec + (tim.tv_usec / 1000000.0);
        walltime = curr_time - start_time;

        std::cout << "Runtime(s): " << walltime << std::endl;
      }

      //for(i=dim;i>=0;i--)
      //	for(j=0;j<CH_elements[i].size();j++)
      //	{
      //		CH_elements[i][j].write();
      //	}
      //for(i=0;i<dim+1;i++)
      //	for(j=0;j<CH_elements[i].size();j++)
      //		std::cout << "rank: " << i << " element: " << j <<  "  volume: " << CH_elements[i][j].nvolume << std::endl;
      //std::cout << "volume: " << CH_elements[dim][0].nvolume << std::endl;

      for(i = 0; i < point_bins.size(); i++)
        if(point_bins[i].size() > 0) {
          std::cout << "probable point_bins error" << std::endl;
          //BP_pause();
          exit(1);
        }
    }



    //std::cout << "finish generate_CH()" << std::endl;
    return 1;

  };

  void Geo::write_hull(int index, int atom, BP_GVec<V_Element> &pts,  BP_Group<V_Element> &hull_pts) {
    int i, j;
    int dm = pts[(unsigned long int) 0].pos.rows();
    std::string filename_base = "hull_";
    std::string filename;

    std::ofstream neb;
    neb.open("hull.neb");
    neb << index + 1 << std::endl;
    for(j = 0; j <= index; j++) {
      std::stringstream ss;
      ss << j;
      filename = filename_base + ss.str();
      neb << filename << std::endl;
    }
    neb.close();

    std::ofstream file;
    std::stringstream ss;
    ss << index;
    filename = filename_base + ss.str();
    file.open(filename.c_str());
    file << "comment" << std::endl;
    file << "1.0" << std::endl;
    file << "100 0 0" << std::endl;
    file << "0 100 0" << std::endl;
    file << "0 0 100" << std::endl;
    //if( atom < m.cols())
    {
      file << "O Sb P" << std::endl;
      file << pts.size() << " " << pts.size() << " 1" << std::endl;
    }
    //else
    //{
    //	file << "O Sb" << std::endl;
    //	file << m.cols() << " " << m.cols() << std::endl;
    //}
    file << "Direct" << std::endl;
    if(dm == 2) {
      for(i = 0; i < pts.size(); i++)
        file << pts[i].pos(0)*scale(0) / 100.0 << " 0.021 " << pts[i].pos(1)*scale(1) / 100.0 << std::endl;
      for(i = 0; i < pts.size(); i++)
        if(i < hull_pts.size())
          file << hull_pts[i].pos(0)*scale(0) / 100.0 << " 0.02 " << hull_pts[i].pos(1)*scale(1) / 100.0 << std::endl;
        else
          file << "0 0.9 0" << std::endl;
      if(atom < pts.size())
        file << pts[atom].pos(0)*scale(0) / 100.0 << " 0 " << pts[atom].pos(1)*scale(1) / 100.0 << std::endl;
      else
        file << "0 0.9 0" << std::endl;
    }
    else if(dm == 3) {
      for(i = 0; i < pts.size(); i++)
        file << pts[i].pos.transpose().cwiseProduct(scale) / 100.0 << std::endl;
      for(i = 0; i < pts.size(); i++)
        if(i < hull_pts.size())
          file << hull_pts[i].pos.transpose().cwiseProduct(scale) / 100.0 << std::endl;
        else
          file << "0 0 0" << std::endl;
      if(atom < pts.size())
        file << pts[atom].pos.transpose().cwiseProduct(scale) / 100.0 << std::endl;
      else
        file << "0 0 0" << std::endl;
    }
    file.close();
  };

  bool Geo::check_hull_point(const Eigen::VectorXd &v, BP_GVec<V_Element> &hull, BP_Graph<V_Element, bool> &graph, BP_Vec<BP_Gen_Vertex *> &interior) {
    //std::cout << "begin check_hull_point()" << std::endl;

    int i, j, c;
    BP_Gen_Vertex *n2, *n3;
    double d;


    // check elements with rank dim-1 for elements that 'see' the new point
    //   these will become 'interior' elements after adding the point and must be cut
    for(i = 0; i < hull.size(); i++) {
      //std::cout << "facet:  " << i << std::endl;
      n2 = graph.member_vert_to_vert(hull.member(i));

      //std::cout << "out_norm: " << vert_val(n2).out_norm.transpose() << std::endl;
      //std::cout << "v-vm: " << (v - vert_val(n2).vertex_mean).transpose() << std::endl;

      d = hull[i].out_norm.dot(v - hull[i].vertex_mean);
      //std::cout << "d: " << d << std::endl;
      c = compare(0, d, Geo_tol);
      if(c == 0) {
        // point is in plane of n2
        //std::cout << "point in plane" << std::endl;
        //check if point already exists
        for(j = 0; j < graph.num_incident_edges(n2); j++) {
          n3 = graph.get_neighbor(n2, j);
          if(graph.vert_val(n3).rank == 0) {
            // if points are equivalent, ignore this point
            if(vector_is(v, graph.vert_val(n3).pos, CH_tol) == 1) {
              //std::cout << "finish check_hull_point() equivalent point" << std::endl;
              return 0;
            }
          }
        }
        // otherwise, this point is on this element's space, so we'll end up deleting it
        hull[i].cut = 1;
        interior.add(n2);

      }
      else if(c == 1) {
        // test point is outside the hull (relative to n2), so n2 will be cut
        //std::cout << "point outside facet" << std::endl;
        hull[i].cut = 1;
        interior.add(n2);

      }
      else if(c == -1) {
        // test point is inside the hull (relative to n2)
        //std::cout << "point inside facet" << std::endl;

      }



    }

    //std::cout << "interior.size(): " << interior.size() << std::endl;
    if(interior.size() > 0) {
      // new hull point
      //std::cout << "finish check_hull_point() add" << std::endl;
      return 1;
    }
    else {
      // not on hull
      //std::cout << "finish check_hull_point() do not add" << std::endl;
      return 0;
    }


  };

  // Breadth First Search of hull_facets, starting with curr_facet
  //   (not exactly a breadth first search as is...)
  bool Geo::check_hull_point_BFS(const Eigen::VectorXd &v, BP_GVec<V_Element> &hull, BP_Graph<V_Element, bool> &graph, BP_Vec<BP_Gen_Vertex *> &interior, int curr_facet) {
    //std::cout << "begin check_hull_point()" << std::endl;

    // The generate_CH() function works by creating a starting simplical hull,
    //	 binning each point to a facet it is 'outside' of,
    //	 then building the hull outwards by checking if each point is outside of the current hull.
    //   Any facets that can 'see' the point being checked are deleted,
    //   and new hull facets created from the new point and the 'ridge' of edges surrounding the
    //   facets that can see the new point.
    //   When a facet is deleted the points binned with it are re-binned with the new facets.
    //
    // This function checks if point 'v' should be added to the current hull
    //	 Returns True if the point is a current hull point, False if it is interior to the current hull
    //   'interior' returns with the set of hull facets that can 'see' the point 'v'
    //   'curr_facet' is a facet point 'v' is outside of

    int i, j, c;
    BP_Gen_Vertex *n1, *n2, *n3;
    double d;

    BP_Vec< BP_Gen_Vertex *> queued;
    //BP_Vec< BP_Gen_Vertex*> inplane;
    BP_Vec< BP_Gen_Vertex *> inside;

    //std::cout << "point: " << v.transpose() << std::endl;

    //n1 = graph.member_vert_to_vert( hull.member(curr_facet));
    //d = graph.vert_val(n1).out_norm.dot(v - graph.vert_val(n1).vertex_mean);
    //std::cout << " curr_facet d: " << d << std::endl;
    //for( int k=0; k<graph.num_incident_edges(n1);k++)
    //{
    //	if( graph.vert_val( graph.get_neighbor( n1, k)).rank == 0)
    //	{
    //		std::cout << "  vert: " << graph.vert_to_index( graph.get_neighbor( n1, k)) << "  pos: " << graph.vert_val( graph.get_neighbor( n1, k)).pos.transpose() << std::endl;
    //	}
    //}

    //std::cout << "add curr_facet" << std::endl;
    hull[curr_facet].cut = 1;
    interior.add(graph.member_vert_to_vert(hull.member(curr_facet)));
    queued.add(graph.member_vert_to_vert(hull.member(curr_facet)));

    while(queued.size() > 0) {
      //loop over neighbor facets to queued[0]
      for(i = 0; i < graph.num_incident_edges(queued[0]); i++) {
        n1 = graph.get_neighbor(queued[0], i);

        // if facet, and not already added
        if(graph.vert_val(n1).rank != 0)
          if(graph.vert_val(n1).cut == 0) {

            // if interior
            d = graph.vert_val(n1).out_norm.dot(v - graph.vert_val(n1).vertex_mean);
            //std::cout << "d: " << d << std::endl;
            c = compare(0, d, Geo_tol);
            if(c == 0) {
              // point is in plane of n1
              //std::cout << "++d: " << d << " norm: " << graph.vert_val(n1).out_norm.transpose() ;
              //std::cout << "  point in plane" << std::endl;
              //check if point already exists
              for(j = 0; j < graph.num_incident_edges(n1); j++) {
                n3 = graph.get_neighbor(n1, j);
                if(graph.vert_val(n3).rank == 0) {
                  // if points are equivalent, ignore this point
                  if(vector_is(v, graph.vert_val(n3).pos, CH_tol) == 1) {
                    //std::cout << "finish check_hull_point() equivalent point" << std::endl;
                    return 0;
                  }
                }
              }

              // otherwise, this point is on this element's space, so we'll end up deleting it
              graph.vert_val(n1).cut = 1;
              interior.add(n1);
              //inplane.add(n1);
              queued.add(n1);

              //for( int k=0; k<graph.num_incident_edges(n1);k++)
              //	{
              //		if( graph.vert_val( graph.get_neighbor( n1, k)).rank == 0)
              //		{
              //			std::cout << "  vert: " << graph.vert_to_index( graph.get_neighbor( n1, k)) << "  pos: " << graph.vert_val( graph.get_neighbor( n1, k)).pos.transpose() << std::endl;
              //		}
              //	}

            }
            else if(c == 1) {
              // test point is outside the hull (relative to n1), so n1 will be cut
              //std::cout << "++d: " << d << " norm: " << graph.vert_val(n1).out_norm.transpose() ;
              //std::cout << "  point outside facet" << std::endl;
              graph.vert_val(n1).cut = 1;
              interior.add(n1);
              queued.add(n1);

              //for( int k=0; k<graph.num_incident_edges(n1);k++)
              //	{
              //		if( graph.vert_val( graph.get_neighbor( n1, k)).rank == 0)
              //		{
              //			std::cout << "  vert: " << graph.vert_to_index( graph.get_neighbor( n1, k)) << "  pos: " << graph.vert_val( graph.get_neighbor( n1, k)).pos.transpose() << std::endl;
              //		}
              //	}
            }
            else if(c == -1) {
              graph.vert_val(n1).cut = 1;
              inside.add(n1);
              // test point is inside the hull (relative to n2)
              //std::cout << "d: " << d << " norm: " << graph.vert_val(n1).out_norm.transpose() ;
              //std::cout << "  point inside facet" << std::endl;

              //if( std::fabs(d) < 1e-5)
              //{
              //	for( int k=0; k<graph.num_incident_edges(n1);k++)
              //	{
              //		if( graph.vert_val( graph.get_neighbor( n1, k)).rank == 0)
              //		{
              //			std::cout << "  vert: " << graph.vert_to_index( graph.get_neighbor( n1, k)) << "  pos: " << graph.vert_val( graph.get_neighbor( n1, k)).pos.transpose() << std::endl;
              //		}
              //	}
              //
              //}

            }
          }
      }

      // remove this one from queue
      queued.remove(0);

    }

    for(i = 0; i < inside.size(); i++)
      graph.vert_val(inside[i]).cut = 0;

    //std::cout << "interior.size(): " << interior.size() << std::endl;
    if(interior.size() > 0) {
      // // since there are 'interior' facets, we won't end up cutting the inplane facets
      //for( i=0; i<inplane.size(); i++)
      //	graph.vert_val(inplane[i]).cut = 0;

      // new hull point
      //std::cout << "finish check_hull_point() add" << std::endl;
      return 1;
    }
    else {
      //if( inplane.size() == 0)
      //{
      // not on hull
      //std::cout << "finish check_hull_point() do not add" << std::endl;
      return 0;
      //}
      //else
      //{
      //	// no 'interior' planes, so the 'inplane' facets will be removed
      //	for( i=0;i<inplane.size();i++)
      //		interior.add( inplane[i]);
      //
      //	return 1;
      //}
    }


  };


  // Form a simplex from the non-planar points in point_list
  void Geo::form_simplex_from_points(BP_GVec<V_Element> &pts, BP_Group<V_Element> &grp_pts, BP_GVec<V_Element> &simplexes, BP_Graph<V_Element, bool> &graph, BP_Group<V_Element> &point_list) {
    //std::cout << "begin form_simplex_from_points()" << std::endl;

    int i, j, k, ref;
    double d;
    //int dim = pts[point_list[0]].pos.rows();
    int dim = point_list[0].pos.rows();
    Eigen::VectorXd simplex_center, pos, vertex_mean, out_norm, vtmp1, vtmp2;
    Eigen::MatrixXd m2(dim, dim - 1);
    Eigen::MatrixXd tQ = Eigen::MatrixXd::Identity(dim, dim);
    Eigen::FullPivHouseholderQR< Eigen::MatrixXd > tQR;
    BP_Gen_Vertex *cell, *n1;
    BP_Group<V_Element> member_list;
    //BP_Vec<BP_Gen_GVec_Member*> member_list;
    BP_Gen_GVec_Member *m1;

    // add an element for the entire cell
    //cell = add_vertex( simplex[dim].add(V_Element(dim, v, tQ)));

    simplex_center.resize(dim);
    simplex_center.setZero();

    // add an element for each vertex
    for(i = 0; i < point_list.size(); i++) {
      //simplex_center += pts[ point_list[i]].pos;
      //grp_pts.add(point_list[i]);
      //graph.add_vertex( point_list[i]);
      simplex_center += point_list[i].pos;
      grp_pts.add(point_list.member(i));
      graph.add_vertex(point_list.member(i));
    }

    simplex_center /= point_list.size();

    //std::cout << "simplex_center: " << simplex_center.transpose() << std::endl;

    //for each point
    for(i = 0; i < point_list.size(); i++) {
      //std::cout << "i: " << i << std::endl;
      // find the plane defined by the other points (!= i)
      // create matrix m2(dim, dim-1) giving position of dim-1 vertices relative to
      //   the 'ref' vertex, which is also in the plane
      k = 0;
      if(i == 0) ref = 1;
      else ref = 0;

      for(j = 0; j < point_list.size(); j++) {
        if(j == i || j == ref) continue;
        //m2.col(k) = pts[point_list[j]].pos - pts[point_list[ref]].pos;
        m2.col(k) = point_list[j].pos - point_list[ref].pos;
        k++;
      }
      //std::cout << "m2:\n" << m2 << std::endl;

      // matrixQ().leftCols(dim-1) gives the basis of the plane
      tQR.compute(m2);
      tQ = tQR.matrixQ();

      //std::cout << "tQ:\n" << tQ << std::endl;
      // vertex_mean
      //std::cout << "m.col(ref): " << m.col(ref).transpose() << std::endl;
      //std::cout << "m2.rowwise().mean(): " << m2.rowwise().mean() << std::endl;

      //vertex_mean = pts[point_list[ref]].pos + m2.rowwise().mean();
      /////vertex_mean = point_list[ref].pos + m2.rowwise().mean();
      vertex_mean = m2.rowwise().sum() / dim + point_list[ref].pos;
      //std::cout << "vertex_mean:" << vertex_mean.transpose() << std::endl;

      // set out_norm to outward normal direction to plane
      vtmp1 = vertex_mean - simplex_center;
      vtmp2 = tQ.rightCols(1);
      d = vtmp1.dot(vtmp2);
      //std::cout << "d: " << d << std::endl;
      if(compare(0, d, Geo_tol) == 1) out_norm = tQ.rightCols(1);
      else if(compare(0, d, Geo_tol) == -1) out_norm = -tQ.rightCols(1);
      else {
        std::cout << "error calculating out_norm A, consider increasing rank determining threshold" << std::endl;
        std::cout << "point_list.size(): " << point_list.size() << std::endl;
        std::cout << "  m2: " << std::endl;
        std::cout << m2 << std::endl;
        std::cout << "rank: " << tQR.rank() << std::endl;
        std::cout << "dim: " << dim << std::endl;
        std::cout << "  R: " << std::endl;
        std::cout << tQR.matrixQR() << std::endl;
      }

      // find the vector normal to the plane giving the plane's position
      //   from the component of one of the vertices along the nullspace of the plane
      //pos = out_norm*out_norm.dot(pts[point_list[ref]].pos);
      pos = out_norm * out_norm.dot(point_list[ref].pos);
      //std::cout << "pos: " << pos.transpose() << std::endl;

      // add plane to simplex
      //m1 = simplexes.add( V_Element(dim-1, pos, vertex_mean, out_norm, tQ ) );
      m1 = simplexes.add(V_Element(dim - 1, pos, vertex_mean, out_norm));
      member_list.add(m1);
      n1 = graph.add_vertex(m1);

      //std::cout << "connect plane and points" << std::endl;
      // connect the plane with each of the other points
      for(j = 0; j < point_list.size(); j++)
        if(j != i)
          graph.connect_vertices(n1, graph.member_vert_to_vert(point_list.member(j)), connections.add(0), 0);
    }

    //std::cout << "connect planes" << std::endl;
    // connect each plane with the other planes
    for(i = 0; i < member_list.size(); i++)
      for(j = i + 1; j < member_list.size(); j++)
        graph.connect_vertices(member_list.member(i), member_list.member(j), connections.add(0), 0);

    //std::cout << "finish form_simplex_from_points()" << std::endl;
  };


  BP_Gen_Vertex *Geo::form_facet_from_points(BP_Gen_Vertex *n1, BP_Vec<BP_Gen_Vertex *> &perimeter, BP_GVec<V_Element> &hull_fcts, BP_Graph<V_Element, bool> &graph, const Eigen::VectorXd &center) {
    //std::cout << "begin form_facet_from_points()" << std::endl;

    int i, j, k, ref;
    double d;
    Eigen::VectorXd simplex_center, pos, vertex_mean, out_norm, vtmp1, vtmp2;
    int dim = center.rows();

    //std::cout << "!!! dim: " << dim << std::endl;

    Eigen::MatrixXd m2(dim, dim - 1);
    Eigen::MatrixXd tQ;// = Eigen::MatrixXd::Identity(dim,dim);
    Eigen::FullPivHouseholderQR< Eigen::MatrixXd > tQR;
    BP_Gen_Vertex *n2;

    // find the plane defined by the  points
    // create matrix m2(dim, dim-1) giving position of dim-1 vertices relative to n1
    //std::cout << "perimeter.size(): " << perimeter.size() << std::endl;
    for(j = 0; j < perimeter.size(); j++) {
      //std::cout << "pos: " <<  graph.vert_val(perimeter[j]).pos.transpose() << std::endl;
      m2.col(j) = graph.vert_val(perimeter[j]).pos - graph.vert_val(n1).pos;
    }
    //std::cout << "n1: " << graph.vert_val(n1).pos.transpose() << std::endl;

    // matrixQ().leftCols(dim-1) gives the basis of the plane
    tQR.compute(m2);
    tQ = tQR.matrixQ();

    if(tQR.rank() < dim - 1) {
      std::cout << "error, too small rank: " << tQR.rank() << std::endl;
      return 0;
    }

    // vertex_mean
    /////vertex_mean = graph.vert_val(n1).pos + m2.rowwise().mean();
    vertex_mean = m2.rowwise().sum() / dim + graph.vert_val(n1).pos;
    //std::cout << "vertex_mean: " << vertex_mean.transpose() << std::endl;
    //std::cout << "center: " << center.transpose() << std::endl;

    // set out_norm to outward normal direction to plane
    vtmp1 = vertex_mean - center;
    vtmp2 = tQ.rightCols(1);
    d = vtmp1.dot(vtmp2);
    if(compare(0, d, Geo_tol) == 1) out_norm = tQ.rightCols(1);
    else if(compare(0, d, Geo_tol) == -1) out_norm = -tQ.rightCols(1);
    else std::cout << "error calculating out_norm B" << std::endl;
    //std::cout << "out_norm: " << out_norm.transpose() << std::endl;

    // find the vector normal to the plane giving the plane's position
    //   from the component of one of the vertices along the nullspace of the plane
    pos = out_norm * out_norm.dot(graph.vert_val(n1).pos);


    // add plane to hull
    //n2 = graph.add_vertex( hull_fcts.add( V_Element(dim-1, pos, vertex_mean, out_norm, tQ ) ));
    n2 = graph.add_vertex(hull_fcts.add(V_Element(dim - 1, pos, vertex_mean, out_norm)));


    // and connect vertices to facet
    for(k = 0; k < perimeter.size(); k++) {
      //std::cout << "connect to perimeter" << std::endl;
      graph.connect_vertices(n2, perimeter[k], connections.add(0), 0);
      //std::cout << "#Connections: " << edge_size() << std::endl;

    }
    //std::cout << "connect to new point" << std::endl;
    graph.connect_vertices(n1, n2, connections.add(0), 0);
    //std::cout << "#Connections: " << edge_size() << std::endl;


    //std::cout << "finish form_facet_from_points() " << std::endl;
    return n2;
  };

  // return 1: if b > a + tol; returns -1 if b < a - tol; return 0 if a - tol < b < a + tol
  int Geo::compare(double a, double b, double tol) const {
    double d = b - a;
    //double tol = 1e-8;
    if(d > tol) return 1;
    else if(d < -tol) return -1;
    else return 0;

  };

  //////////////////////////
  /// Geo CH public functions

  void Geo::CH_write(std::ostream &sout) {
    BP_Vec<int> a;
    int i, j;
    sout << "Facets: " << std::endl;
    for(i = 0; i < CH_facets_size(); i++) {
      sout << i << ": " << std::endl;
      sout << "      area: " << CH_facets_area(i)  << std::endl;
      sout << "      norm: " << CH_facets_norm(i).transpose() << std::endl;

      a = CH_facets_nborverts(i);
      sout << "     verts: " ;
      for(j = 0; j < a.size(); j++)
        sout << " " << a[j] ;
      sout << std::endl;

      a = CH_facets_nborfacets(i);
      sout << "    facets: " ;
      for(j = 0; j < a.size(); j++)
        sout << " " << a[j] ;
      sout << std::endl;
    }
    sout << std::endl;
    sout << "Verts: " << std::endl;
    for(i = 0; i < CH_verts_size(); i++) {
      sout << i << ": " << std::endl;
      sout << "       pos: " << CH_verts_pos(i).transpose()  << std::endl;

      a = CH_verts_nborverts(i);
      sout << "     verts: " ;
      for(j = 0; j < a.size(); j++)
        sout << " " << a[j] ;
      sout << std::endl;

      a = CH_verts_nborfacets(i);
      sout << "    facets: " ;
      for(j = 0; j < a.size(); j++)
        sout << " " << a[j] ;
      sout << std::endl;
    }
  };

  // filter so that queries only include facets/points on the bottom of the hull,
  //   with the "bottom" direction defined by v
  void Geo::CH_bottom(const Eigen::VectorXd &v) {
    CH_data.clear();
    CH_data_exists = 0;
    ch_volume = -1;
    ch_area = -1;
    use_bottom = 1;
    bottom = v.cwiseQuotient(scale).normalized();

    CH_filter_points.clear();
    CH_filter_facets.clear();

    int i, j, k;
    bool first_found;
    BP_Gen_Vertex *n1;
    BP_GVec_Member<V_Element> *m1, *m2;

    // find bottom hull facets
    for(i = 0; i < CH_facets.size(); i++)
      if(use_facet(i)) CH_filter_facets.add(CH_facets.member(i));


    // find bottom hull points
    //
    // for each facet
    for(i = 0; i < CH_filter_facets.size(); i++) {
      m1 = CH_filter_facets.member(i);
      // add each neighboring vertex
      for(j = 0; j < CH_graph.num_incident_edges(m1); j++) {
        m2 = CH_graph.get_neighbor(m1, j);
        if(m2->get_val().rank == 0)
          CH_filter_points.add(m2);
      }
    }



  };

  void Geo::CH_wholehull() {
    if(use_bottom == 1) {
      ch_volume = -1;
      ch_area = -1;
      CH_data.clear();
      CH_data_exists = 0;

      CH_filter_points.clear();
      for(int i = 0; i < CH_points.size(); i++)
        CH_filter_points.add(CH_points.member(i));

      CH_filter_facets.clear();
      for(int i = 0; i < CH_facets.size(); i++)
        CH_filter_facets.add(CH_facets.member(i));

    }
    use_bottom = 0;

  };

  //returns volume contained within convex hull in dimension dim
  //  use_bottom does not affect the returned volume
  double Geo::CH_volume() {
    //std::cout << "begin CH_volume()" << std::endl;

    if(ch_volume > 0) return ch_volume;
    if(ch_area < 0) CH_area();

    ch_volume = 0;
    double height;
    int i;

    Eigen::VectorXd v;
    Eigen::VectorXd center;
    center.resize(dim);
    center.setZero();
    //std::cout << "here 1" << std::endl;
    for(i = 0; i < CH_points.size(); i++) {
      center += CH_points[i].pos;
    }
    center /= CH_points.size();

    //std::cout << "here 2" << std::endl;
    for(i = 0; i < CH_filter_facets.size(); i++) {
      v = center - CH_filter_facets[i].vertex_mean;
      height = v.cwiseProduct(scale).dot(-CH_filter_facets[i].out_norm.cwiseProduct(scale).normalized());
      ch_volume += (1.0 / dim) * CH_filter_facets[i].nvolume * height;
    }

    //std::cout << "finish CH_volume()" << std::endl;
    return ch_volume;

  };

  //returns area of convex hull in dimension dim-1
  //  if use_bottom is set, the area only includes those facets with out_norm vector that has positive component along bottom vector
  double Geo::CH_area() {
    //std::cout << "begin CH_area()" << std::endl;
    if(ch_area > 0) return ch_area;

    ch_area = 0;

    int i, j, r;
    BP_Vec<int> nverts;
    Eigen::VectorXd ref;
    Eigen::MatrixXd m1, m2;
    m1.resize(dim, dim - 1);
    Eigen::FullPivHouseholderQR< Eigen::MatrixXd > tQR;
    bool include;

    double fact = 1;
    for(i = 2; i < dim; i++)
      fact *= i;

    for(i = 0; i < CH_filter_facets.size(); i++) {

      nverts = CH_facets_nborverts(i);
      ref = CH_verts_pos(nverts[0]); //CH_filter_points[nverts[0]].pos; //is scaled to full range
      for(j = 1; j < nverts.size(); j++) {
        m1.col(j - 1) =  CH_verts_pos(nverts[j]) - ref;   //CH_filter_points[nverts[j]].pos - ref; //is scaled to full range
      }

      tQR.compute(m1);
      m2 = tQR.matrixQR();
      r = tQR.rank();

      CH_filter_facets[i].nvolume = 1.0;
      for(j = 0; j < r; j++)
        CH_filter_facets[i].nvolume *= m2(j, j);
      CH_filter_facets[i].nvolume = std::fabs(CH_filter_facets[i].nvolume);
      CH_filter_facets[i].nvolume /= fact;

      ch_area += CH_filter_facets[i].nvolume;
    }

    //std::cout << "finish CH_area()" << std::endl;
    return ch_area;

  };

  bool Geo::use_facet(int i) {
    if(use_bottom == 0) return 1;

    if(compare(0, bottom.dot(CH_facets[i].out_norm), Geo_tol) == 1) return 1;
    else return 0;
  }

  // I should probably automatically run this after generate_CH() so that co-planar points get included in the hull
  //   could then also merge co-planar facets or triangulate them.
  //   But this could be very slow... would be better to identify during generate_CH()
  void Geo::CH_calc_distances() {
    if(CH_data_exists) return;
    if(CH_exists == 0)
      calc_CH();

    int i, j, k, c;
    BP_Vec<int> min_array;
    Eigen::VectorXd v;
    double a, b, d;
    double min_d;
    bool first_found;
    BP_Gen_Vertex *n1;

    if(CH_data.size() != 0) std::cout << "error calc_CH_distances 0" << std::endl;

    for(i = 0; i < points.size(); i++) {
      CH_data.add();
      CH_data[i].is_hull_point = 0;
    }
    // loop over hull points
    for(i = 0; i < CH_filter_points.size(); i++) {
      //j = CH_filter_points.member(i)->BP_GVec_index;
      //CH_data[j].is_hull_point = 1;
      CH_data[ points.get_index(CH_filter_points.member(i)) ].is_hull_point = 1;
    }

    for(i = 0; i < points.size(); i++) {
      if(CH_data[i].is_hull_point == 1) {
        CH_data[i].dist_to_hull = 0;
        CH_data[i].vec_to_hull.resize(dim);
        CH_data[i].vec_to_hull.setZero();
        //CH_data[i].closest_facet = -1;
        CH_data[i].closest_facet.clear();
        continue;
      }

      first_found = 0;
      min_array.clear();
      for(j = 0; j < CH_filter_facets.size(); j++) {
        if(use_bottom == 0) {
          v = CH_filter_facets[j].pos - points[i].pos;
          d = v.dot(CH_filter_facets[j].out_norm);
          if(compare(0, d, Geo_tol) == -1) {
            std::cout << "error, pos outside hull. point: " << i << " pos: " << points[i].pos.cwiseProduct(scale).transpose() << std::endl;
            exit(1);
          }
          if(compare(0, d, Geo_tol) == 0) {
            std::cout << "warning, point " << i << " ( " << points[i].pos.cwiseProduct(scale).transpose() << " ) is also on the hull." << std::endl;
            std::cout << "  point " << i << " is co-planar with the following points: " << std::endl;
            BP_Vec<int> nborverts = CH_facets_nborverts(j);
            BP_Vec<int> vert_indices = CH_verts_indices();
            for(k = 0; k < nborverts.size(); k++) {
              std::cout << "  point: " << vert_indices[nborverts[k]] << "  ( " << CH_verts_pos(nborverts[k]).transpose() << " )" << std::endl;
            }
          }

          if(first_found == 0) {
            first_found = 1;
            min_d = d;
            min_array.add(j);
          }
          else {
            c = compare(min_d, d, Geo_tol);
            if(c == -1) {
              min_d = d;
              min_array.clear();
              min_array.add(j);
            }
            else if(c == 0) {
              min_array.add(j);
            }
          }
        }
        else {
          a = (CH_filter_facets[j].pos - points[i].pos).dot(CH_filter_facets[j].out_norm);
          b = bottom.dot(CH_filter_facets[j].out_norm);
          d = a / b;

          if(compare(0, d, Geo_tol) == 0) {
            std::cout << "warning, point " << i << " ( " << points[i].pos.cwiseProduct(scale).transpose() << " ) is also on the bottom hull." << std::endl;
            std::cout << "  point " << i << " is co-planar with the following points: " << std::endl;
            BP_Vec<int> nborverts = CH_facets_nborverts(j);
            BP_Vec<int> vert_indices = CH_verts_indices();
            for(k = 0; k < nborverts.size(); k++) {
              std::cout << "  point: " << vert_indices[nborverts[k]] << "  ( " << CH_verts_pos(nborverts[k]).transpose() << " )" << std::endl;
            }
          }
          if(compare(0, d, Geo_tol) == -1) {
            std::cout << "error, pos outside bottom hull. point: " << i << " pos: " << points[i].pos.cwiseProduct(scale).transpose() << std::endl;
            exit(1);
          }

          if(first_found == 0) {
            first_found = 1;
            min_d = d;
            min_array.add(j);
          }
          else {
            c = compare(min_d, d, Geo_tol);
            if(c == -1) {
              min_d = d;
              min_array.clear();
              min_array.add(j);
            }
            else if(c == 0) {
              min_array.add(j);
            }
          }
        }
      }

      CH_data[i].dist_to_hull = min_d;
      if(!use_bottom) CH_data[i].vec_to_hull = min_d * CH_filter_facets[min_array[0]].out_norm;
      else CH_data[i].vec_to_hull = min_d * bottom;
      CH_data[i].closest_facet = min_array;
    }

    CH_data_exists = 1;

  };



  // returns an array of indices of the closest facets to point i1
  // array size is >1 if projection falls on hull vertex, edge, etc. within tol
  //  returns -1 if hull point, since there are multiple and they should be found with CH_verts_nborfacets()
  //  is !use_bottom, this it the closest facet
  //  if use_bottom, this is the facet which is intersected when moving along bottom vector from the point to the hull
  BP_Vec<int> Geo::CH_closest_facet(int i1) {
    if(CH_data_exists == 0)
      CH_calc_distances();

    return CH_data[i1].closest_facet;
  };

  // distance of point i1 (from points list) to hull
  //  if !use_bottom, this is the closest distance
  //  if use_bottom, this is the distance along the bottom vector from point to the hull
  double Geo::CH_dist_to_hull(int i1) {
    if(CH_data_exists == 0) {
      CH_calc_distances();
    }
    if(CH_data[i1].is_hull_point) return 0;


    //return CH_data[i1].dist_to_hull;
    return CH_vec_to_hull(i1).norm();



  };

  // vector from point i1 (from points list) to hull
  //  if !use_bottom, this is the closest vector to hull
  //  if use_bottom, this is the vector in the direction of the bottom vector from point to the hull
  Eigen::VectorXd Geo::CH_vec_to_hull(int i1) {
    if(CH_data_exists == 0)
      CH_calc_distances();

    return CH_data[i1].vec_to_hull.cwiseProduct(scale);

  };

  // return distance from point to hull.
  //		Positive distances are for points inside hull
  //		Negative distances are for points outside hull
  double Geo::CH_dist_to_hull(const Eigen::VectorXd &i1) const {

    CH_data_class tmp_CH_data;
    int i, j, k, c;
    BP_Vec<int> min_array;
    Eigen::VectorXd v, pos;
    double a, b, d;
    double min_d;
    bool first_found;
    BP_Gen_Vertex *n1;
    tmp_CH_data.is_hull_point = 0;

    pos = i1.cwiseQuotient(scale);

    first_found = 0;
    min_array.clear();
    for(j = 0; j < CH_filter_facets.size(); j++) {
      if(use_bottom == 0) {
        v = CH_filter_facets[j].pos - pos;
        d = v.dot(CH_filter_facets[j].out_norm);
        //if( compare(0,d,Geo_tol) == -1) std::cout << "error, pos outside of hull" << std::endl;
        //if( compare(0,d,Geo_tol) == 0) std::cout << "warning, pos on hull" << std::endl;
        if(first_found == 0) {
          first_found = 1;
          min_d = d;
          min_array.add(j);
        }
        else {
          c = compare(min_d, d, Geo_tol);
          if(c == -1) {
            min_d = d;
            min_array.clear();
            min_array.add(j);
          }
          else if(c == 0) {
            min_array.add(j);
          }
        }
      }
      else {
        a = (CH_filter_facets[j].pos - pos).dot(CH_filter_facets[j].out_norm);
        b = bottom.dot(CH_filter_facets[j].out_norm);
        d = a / b;

        //if( compare(0,d,Geo_tol) == 0) std::cout << "warning, pos on hull (bottom)" << std::endl;
        //if( compare(0,d,Geo_tol) == -1) std::cout << "warning, pos outside hull (bottom)" << std::endl;

        if(first_found == 0) {
          first_found = 1;
          min_d = d;
          min_array.add(j);
        }
        else {
          c = compare(min_d, d, Geo_tol);
          if(c == -1) {
            min_d = d;
            min_array.clear();
            min_array.add(j);
          }
          else if(c == 0) {
            min_array.add(j);
          }
        }
      }
    }

    tmp_CH_data.dist_to_hull = min_d;
    if(!use_bottom) tmp_CH_data.vec_to_hull = min_d * CH_filter_facets[min_array[0]].out_norm;
    else tmp_CH_data.vec_to_hull = min_d * bottom;
    tmp_CH_data.closest_facet = min_array;

    if(min_d == 0)
      return 0;
    else if(min_d > 0)
      return tmp_CH_data.vec_to_hull.cwiseProduct(scale).norm();
    else
      return -tmp_CH_data.vec_to_hull.cwiseProduct(scale).norm();

  }

  //Eigen::VectorXd CH_vec_to_hull( const Eigen::VectorXd &i1)

  int Geo::CH_facets_size() {
    return CH_filter_facets.size();
  };

  double	Geo::CH_facets_area(int i1) {
    if(ch_area < 0)
      CH_area();

    return CH_filter_facets[i1].nvolume;
  };

  Eigen::VectorXd Geo::CH_facets_norm(int i1) {
    return CH_filter_facets[i1].out_norm.cwiseProduct(scale);
  };

  int Geo::CH_facets_nborverts_size(int i1) {
    //std::cout << "begin CH_facets_nborverts_size()" << std::endl;
    BP_GVec_Member<V_Element> *m1, *m2;
    //std::cout << "dim: " << dim << std::endl;
    int i, j, k;
    m1 = CH_filter_facets.member(i1);
    k = 0;
    for(i = 0; i < CH_graph.num_incident_edges(m1); i++) {
      m2 = CH_graph.get_neighbor(m1, i);
      j = CH_filter_points.get_index(m2);
      if(j != -1) {
        k++;

      }
    }
    //std::cout << "finish CH_facets_nborverts_size()" << std::endl;
    return k;
  };

  BP_Vec<int> Geo::CH_facets_nborverts(int i1) {
    //std::cout << "begin CH_facets_nborverts()" << std::endl;
    BP_GVec_Member<V_Element> *m1, *m2;
    BP_Vec<int> v;
    //std::cout << "dim: " << dim << std::endl;
    int i, j;
    m1 = CH_filter_facets.member(i1);
    for(i = 0; i < CH_graph.num_incident_edges(m1); i++) {
      m2 = CH_graph.get_neighbor(m1, i);
      j = CH_filter_points.get_index(m2);
      if(j != -1) {
        v.add(j);			/// filtered index

      }
    }
    //std::cout << "finish CH_facets_nborverts()" << std::endl;

    return v;
  };

  int Geo::CH_facets_nborfacets_size(int i1) {
    //std::cout << "begin CH_facets_nborfacets_size()" << std::endl;
    BP_GVec_Member<V_Element> *m1, *m2;
    //std::cout << "dim: " << dim << std::endl;
    int i, j, k;
    m1 = CH_filter_facets.member(i1);
    k = 0;
    for(i = 0; i < CH_graph.num_incident_edges(m1); i++) {
      m2 = CH_graph.get_neighbor(m1, i);
      j = CH_filter_facets.get_index(m2);
      if(j != -1) {
        k++;			/// filtered index

      }
    }
    //std::cout << "finish CH_facets_nborfacets_size()" << std::endl;
    return k;
  };

  BP_Vec<int> Geo::CH_facets_nborfacets(int i1) {
    //std::cout << "begin CH_facets_nborfacets()" << std::endl;
    BP_GVec_Member<V_Element> *m1, *m2;
    BP_Vec<int> v;

    int i, j;
    m1 = CH_filter_facets.member(i1);
    for(i = 0; i < CH_graph.num_incident_edges(m1); i++) {
      m2 = CH_graph.get_neighbor(m1, i);
      j = CH_filter_facets.get_index(m2);
      if(j != -1) {
        v.add(j);			/// filtered index
      }
    }
    //std::cout << "finish CH_facets_nborfacets()" << std::endl;

    return v;
  };

  int Geo::CH_verts_size() {
    return CH_filter_points.size();
  };

  Eigen::VectorXd Geo::CH_verts_pos(int i1) {
    return CH_filter_points[i1].pos.cwiseProduct(scale);
  };

  BP_Vec<int> Geo::CH_verts_indices() {
    BP_Vec<int> index_list;

    for(int i = 0; i < CH_filter_points.size(); i++) {
      index_list.add(points.get_index(CH_filter_points.member(i)));
    }

    return index_list;
  };

  int Geo::CH_verts_nborverts_size(int i1) {
    //std::cout << "begin CH_verts_nborverts_size()" << std::endl;
    BP_GVec_Member<V_Element> *m1, *m2, *m3;
    BP_Vec<int> facet_list1;
    BP_Vec<int> facet_list2;
    int i, j, k, ii, i2, j2;
    k = 0;
    m1 = CH_filter_points.member(i1);
    facet_list1 = CH_verts_nborfacets(i1);
    for(i = 0; i < CH_graph.num_incident_edges(m1); i++) {
      m2 = CH_graph.get_neighbor(m1, i);
      if(m2->get_val().rank != 0) {
        for(j = 0; j < CH_graph.num_incident_edges(m2); j++) {
          m3 = CH_graph.get_neighbor(m2, j);
          ii = CH_filter_points.get_index(m3);
          if(ii != -1)
            if(m3 != m1) {
              // check they share a facet, necessary for not including nbors by non-bottom facets
              facet_list2 = CH_verts_nborfacets(ii);
              for(i2 = 0; i2 < facet_list1.size(); i2++)
                for(j2 = 0; j2 < facet_list2.size(); j2++) {
                  if(facet_list1[i2] == facet_list2[j2]) {
                    i2 = facet_list1.size();
                    j2 = facet_list2.size();
                    k++;
                  }
                }


            }
        }
      }
    }
    //std::cout << "finish CH_verts_nborverts_size()" << std::endl;
    return k;
  };


  // this is not perfect when using 'bottom'
  //  it will include nborverts that are not nbor by 'bottom' facets
  BP_Vec<int> Geo::CH_verts_nborverts(int i1) {
    //std::cout << "begin CH_verts_nborverts()" << std::endl;
    BP_Vec<int> v;
    BP_GVec_Member<V_Element> *m1, *m2, *m3;
    BP_Vec<int> facet_list1;
    BP_Vec<int> facet_list2;
    int i, j, k, i2, j2;
    m1 = CH_filter_points.member(i1);
    facet_list1 = CH_verts_nborfacets(i1);
    for(i = 0; i < CH_graph.num_incident_edges(m1); i++) {
      m2 = CH_graph.get_neighbor(m1, i);
      if(m2->get_val().rank != 0) {
        for(j = 0; j < CH_graph.num_incident_edges(m2); j++) {
          m3 = CH_graph.get_neighbor(m2, j);
          k = CH_filter_points.get_index(m3);
          if(k != -1)
            if(CH_filter_points[k].cut != 1)
              if(m3 != m1) {
                // check they share a facet, necessary for not including nbors by non-bottom facets
                facet_list2 = CH_verts_nborfacets(k);
                for(i2 = 0; i2 < facet_list1.size(); i2++)
                  for(j2 = 0; j2 < facet_list2.size(); j2++) {
                    if(facet_list1[i2] == facet_list2[j2]) {
                      i2 = facet_list1.size();
                      j2 = facet_list2.size();
                      v.add(k);
                      CH_filter_points[k].cut = 1;

                    }
                  }

              }
        }
      }
    }
    //std::cout << "finish CH_verts_nborverts()" << std::endl;

    for(i = 0; i < v.size(); i++)
      CH_filter_points[v[i]].cut = 0;

    return v;
  };

  int Geo::CH_verts_nborfacets_size(int i1) {
    //std::cout << "begin CH_verts_nborfacets_size()" << std::endl;
    BP_GVec_Member<V_Element> *m1, *m2;
    int i, j, k;
    m1 = CH_filter_points.member(i1);
    k = 0;
    for(i = 0; i < CH_graph.num_incident_edges(m1); i++) {
      m2 = CH_graph.get_neighbor(m1, i);
      j = CH_filter_facets.get_index(m2);
      if(j != -1) {
        k++;
      }
    }
    //std::cout << "finish CH_verts_nborfacets_size()" << std::endl;
    return j;
  };

  BP_Vec<int> Geo::CH_verts_nborfacets(int i1) {
    //std::cout << "begin CH_verts_nborverts()" << std::endl;
    BP_Vec<int> v;
    BP_GVec_Member<V_Element> *m1, *m2;
    int i, j;
    m1 = CH_filter_points.member(i1);
    for(i = 0; i < CH_graph.num_incident_edges(m1); i++) {
      m2 = CH_graph.get_neighbor(m1, i);
      j = CH_filter_facets.get_index(m2);
      if(j != -1) {
        v.add(j);			/// filtered index
      }
    }
    //std::cout << "finish CH_verts_nborverts()" << std::endl;
    return v;
  };


  ///////////////////////////////////////
  /// Delaunay Triangulation functions

  bool Geo::calc_DT() {
    if(DT_exists) return 1;



    // Create DT by finding hull in dim+1 of points with value of dimension dim+1 == sum(sqr(pos))
    //
    BP_GVec< V_Element> tmp_points;
    BP_GVec< V_Element> tmp_CH_facets;
    BP_Graph<V_Element, bool> tmp_CH_graph;
    BP_Group< V_Element> tmp_CH_points;

    int i, j, k;
    BP_Gen_GVec_Member *m1, *m2;
    double d;

    Eigen::VectorXd v;
    v.resize(dim + 1);

    Eigen::VectorXd vzero;
    vzero.resize(dim);
    vzero.setZero();

    Eigen::VectorXd tmp_bottom;
    tmp_bottom.resize(dim + 1);
    tmp_bottom.setZero();
    tmp_bottom(dim) = -1;
    //std::cout << "tmp_bottom: " << tmp_bottom.transpose() << std::endl;

    //Eigen::MatrixXd Midentity_A = Eigen::MatrixXd::Identity(dim,dim);
    //Eigen::MatrixXd Midentity_B = Eigen::MatrixXd::Identity(dim+1,dim+1);

    // add points to tmp_points, with last dimension = sum(sqr)
    //std::cout << "points.size(): " << points.size() << std::endl;
    tmp_points.capacity(points.size());
    for(i = 0; i < points.size(); i++) {
      v.head(dim) = points[i].pos;
      d = points[i].pos.norm();
      v(dim) = sqr(d);
      //tmp_points.add( V_Element( 0, v, Midentity_B) );
      tmp_points.add(V_Element(0, v));
      //std::cout << v.transpose() << std::endl;
    }

    //for( i=0; i<tmp_points.size(); i++)
    //	std::cout << tmp_points[i].pos.transpose() << std::endl;


    if(verbosity >= 10)
      std::cout << "begin hull finding in dimension " << dim + 1 << std::endl;
    if(generate_CH(tmp_points, tmp_CH_points, tmp_CH_facets, tmp_CH_graph)) {
      if(tmp_CH_points.size() != tmp_points.size()) {
        std::cout << "ERROR, not all points on hull.  Maybe try increasing tol." << std::endl;
        abort();
      }
      if(verbosity >= 10)
        std::cout << "finish DT hull finding" << std::endl;


      // now use this graph to set the DT_graph
      //
      // remove 'top' facets (include 'side' facets in 'top)
      //std::cout << "# facets: " << tmp_CH_facets.size() << std::endl;
      BP_Vec<BP_Gen_GVec_Member *> top;
      for(i = 0; i < tmp_CH_facets.size(); i++) {
        if(compare(0, tmp_bottom.dot(tmp_CH_facets[i].out_norm), Geo_tol) != 1) {
          top.add(tmp_CH_facets.member(i));
        }

      }

      //return 1;


      if(verbosity >= 10)
        std::cout << "# top facets: " << top.size() << std::endl;

      for(i = 0; i < top.size(); i++)
        tmp_CH_facets.remove(top[i]);

      i = tmp_CH_facets.size();
      if(i < 1000) i = 1000;
      j = i / 10;

      DT_simplexes.capacity(i);			//CONSTANT
      DT_simplexes.cap_increment(j);		//CONSTANT
      DT_graph.capacity(i);
      DT_graph.cap_increment(j);


      // add tmp_CH_facets as DT_graph vertices
      for(i = 0; i < tmp_CH_facets.size(); i++) {
        //DT_graph.add_vertex( DT_simplexes.add(V_Element(dim,vzero,Midentity_A)));
        DT_graph.add_vertex(DT_simplexes.add(V_Element(dim, vzero)));
      }

      // add points to DT_graph
      for(i = 0; i < points.size(); i++) {
        DT_graph.add_vertex(points.member(i));
      }

      // add connections
      for(i = 0; i < tmp_CH_facets.size(); i++) {
        m1 = tmp_CH_facets.member(i);

        for(j = 0; j < tmp_CH_graph.num_incident_edges(m1); j++) {
          m2 = tmp_CH_graph.get_neighbor(m1, j);
          k = tmp_points.get_index(m2);
          if(k != -1) {
            DT_graph.connect_vertices(points.member(k), DT_simplexes.member(i), connections.add(0), 0);
          }
          k = tmp_CH_facets.get_index(m2);
          if(k != -1) {
            if(k > i)
              DT_graph.connect_vertices(DT_simplexes.member(k), DT_simplexes.member(i), connections.add(0), 0);
          }
        }
      }

      if(verbosity >= 10)
        std::cout << "#DT Simps: " << DT_simplexes.size() << std::endl;

      DT_calc_data();

      DT_exists = 1;

    }
    else {
      if(verbosity >= 10)
        std::cout << "hull failed!" << std::endl;
    }
    return DT_exists;

  };

  void Geo::DT_write(std::ostream &sout) {
    BP_Vec<int> a;
    int i, j;
    sout << "Total volume: " << DT_volume() << std::endl;
    sout << "Simplexes: " << std::endl;
    for(i = 0; i < DT_simps_size(); i++) {
      sout << i << ": " << std::endl;
      sout << "      vol: " << DT_simps_volume(i)  << std::endl;
      sout << "      circenter: " << DT_simps_circumcenter(i).transpose() << std::endl;

      a = DT_simps_nborverts(i);
      sout << "     verts: " ;
      for(j = 0; j < a.size(); j++)
        sout << " " << a[j] ;
      sout << std::endl;

      a = DT_simps_nborsimps(i);
      sout << "     simps: " ;
      for(j = 0; j < a.size(); j++)
        sout << " " << a[j] ;
      sout << std::endl;
    }
    sout << std::endl;
    sout << "Verts: " << std::endl;
    for(i = 0; i < DT_verts_size(); i++) {
      sout << i << ": " << std::endl;
      sout << "       pos: " << DT_verts_pos(i).transpose()  << std::endl;

      a = DT_verts_nborverts(i);
      sout << "     verts: " ;
      for(j = 0; j < a.size(); j++)
        sout << " " << a[j] ;
      sout << std::endl;

      a = DT_verts_nborsimps(i);
      sout << "     simps: " ;
      for(j = 0; j < a.size(); j++)
        sout << " " << a[j] ;
      sout << std::endl;
    }
  };

  int	Geo::DT_simps_size() {
    return DT_simplexes.size();
  };

  double Geo::DT_volume() {
    if(DT_data_exists == 0)
      DT_calc_data();

    return dt_volume;
  };

  void Geo::DT_calc_data() {
    //std::cout << "begin DT_calc_data()" << std::endl;
    Eigen::VectorXd ref, ccenter, x, b, b2, tmp_v;
    Eigen::MatrixXd m1, m2, m3;
    int i, j, k, r;
    BP_Vec<int> nverts;
    Eigen::FullPivHouseholderQR< Eigen::MatrixXd > tQR;
    DT_data.capacity(DT_simplexes.size());
    DT_data.clear();

    //std::cout << "dim: " << dim << std::endl;

    double fact = 1;
    for(i = 2; i <= dim; i++)
      fact *= i;

    //std::cout << "fact: " << fact << std::endl;

    m1.resize(dim, dim);
    m3.resize(dim + 1, dim + 1);
    m3.col(0).setOnes();
    b.resize(dim + 1);

    dt_volume = 0;

    for(i = 0; i < DT_simplexes.size(); i++) {
      // calculate simplex volume

      nverts = DT_simps_nborverts(i);
      ref = DT_verts_pos(nverts[0]);	//is scaled to rull range						//CH_filter_points[nverts[0]].pos;

      for(j = 1; j < nverts.size(); j++) {
        m1.col(j - 1) = DT_verts_pos(nverts[j]) - ref;	//is scaled to full range	//CH_filter_points[nverts[j]].pos - ref;
      }

      tQR.compute(m1);
      m2 = tQR.matrixQR();
      r = tQR.rank();

      if(r < dim) {
        std::cout << "error, DT rank < dim.  rank == " << r << "  dim == " << dim << std::endl;
        for(j = 0; j < nverts.size(); j++)
          std::cout << "j: " << j << "  vert: " << DT_verts_pos(nverts[j]).transpose() << std::endl;
        //abort();
      }

      DT_simplexes[i].nvolume = 1.0;
      for(j = 0; j < r; j++)
        DT_simplexes[i].nvolume *= m2(j, j);
      DT_simplexes[i].nvolume = std::fabs(DT_simplexes[i].nvolume);
      DT_simplexes[i].nvolume /= fact;

      // add to overall volume
      dt_volume += DT_simplexes[i].nvolume;

      //std::cout << "simp: " << i << " volume: " << DT_simplexes[i].nvolume << "  rank: " << r << std::endl;

      // calculate simplex circumcenter

      // using method described here: http://www.ics.uci.edu/~eppstein/junkyard/circumcenter.html
      // for each point p (up to N points), sum_i[ ( v_i,p - c_i,p)^2 ] = R^2,  i is index over dimension D
      //   make substitution q = R^2 - sum_i[ (c_i)^2 ]
      //   then  [1 2v_1,1 2v_2,1 ... 2v_D,1]*[ q   ] = [  sum_i[ (v_i,1)^2]  ]
      //         [1 2v_1,2 2v_2,2 ... 2v_D,2]*[ c_1 ] = [  sum_i[ (v_i,2)^2]  ]
      //         [1 2v_1,3 2v_2,3 ... 2v_D,2]*[ c_2 ] = [  sum_i[ (v_i,3)^2]  ]
      //         [...                       ]*[ ... ] = [  ...                ]
      //         [1 2v_1,P 2v_2,P ... 2v_D,P]*[ c_D ] = [  sum_i[ (v_i,P)^2]  ]

      //std::cout << " here 1" << std::endl;

      for(j = 0; j < nverts.size(); j++) {
        tmp_v = DT_verts_pos(nverts[j]);				// is scaled to full range

        for(k = 0; k < dim; k++) {
          m3(j, k + 1) = 2 * tmp_v(k);
        }

        b(j) = tmp_v.squaredNorm();

      }

      tQR.compute(m3);
      x = tQR.solve(b);
      ccenter = x.tail(dim);	// is scaled to full range

      DT_data.add(ccenter);

    }





    DT_data_exists = 1;
  };

  double Geo::DT_simps_volume(int i1) {
    if(DT_data_exists == 0)
      DT_calc_data();

    return DT_simplexes[i1].nvolume;
  };

  Eigen::VectorXd Geo::DT_simps_circumcenter(int i1) {
    if(DT_data_exists == 0)
      DT_calc_data();

    return DT_data[i1].ccenter;  // is already scaled
  };

  int	Geo::DT_simps_nborverts_size(int i1) {
    //std::cout << "begin DT_simps_nborverts_size()" << std::endl;
    BP_GVec_Member<V_Element> *m1, *m2;
    //std::cout << "dim: " << dim << std::endl;
    int i, j, k;
    m1 = DT_simplexes.member(i1);
    k = 0;
    for(i = 0; i < DT_graph.num_incident_edges(m1); i++) {
      m2 = DT_graph.get_neighbor(m1, i);
      j = points.get_index(m2);
      if(j != -1) {
        k++;

      }
    }
    //std::cout << "finish DT_simps_nborverts_size()" << std::endl;
    return k;
  };

  BP_Vec<int> Geo::DT_simps_nborverts(int i1) {
    //std::cout << "begin DT_simps_nborverts()" << std::endl;
    BP_GVec_Member<V_Element> *m1, *m2;
    BP_Vec<int> v;
    //std::cout << "dim: " << dim << std::endl;
    int i, j, k;
    m1 = DT_simplexes.member(i1);
    k = 0;
    for(i = 0; i < DT_graph.num_incident_edges(m1); i++) {
      m2 = DT_graph.get_neighbor(m1, i);
      j = points.get_index(m2);
      if(j != -1) {
        v.add(j);

      }
    }
    //std::cout << "finish DT_simps_nborverts()" << std::endl;
    return v;
  };

  int	Geo::DT_simps_nborsimps_size(int i1) {
    //std::cout << "begin DT_simps_nborsimps_size()" << std::endl;
    BP_GVec_Member<V_Element> *m1, *m2;
    //std::cout << "dim: " << dim << std::endl;
    int i, j, k;
    m1 = DT_simplexes.member(i1);
    k = 0;
    for(i = 0; i < DT_graph.num_incident_edges(m1); i++) {
      m2 = DT_graph.get_neighbor(m1, i);
      j = DT_simplexes.get_index(m2);
      if(j != -1) {
        k++;

      }
    }
    //std::cout << "finish DT_simps_nborsimps_size()" << std::endl;
    return k;
  };

  BP_Vec<int> Geo::DT_simps_nborsimps(int i1) {
    //std::cout << "begin DT_simps_nborsimps_size()" << std::endl;
    BP_GVec_Member<V_Element> *m1, *m2;
    BP_Vec<int> v;
    //std::cout << "dim: " << dim << std::endl;
    int i, j, k;
    m1 = DT_simplexes.member(i1);
    k = 0;
    for(i = 0; i < DT_graph.num_incident_edges(m1); i++) {
      m2 = DT_graph.get_neighbor(m1, i);
      j = DT_simplexes.get_index(m2);
      if(j != -1) {
        v.add(j);

      }
    }
    //std::cout << "finish DT_simps_nborsimps_size()" << std::endl;
    return v;
  };

  int	Geo::DT_verts_size() {
    return points.size();
  };

  Eigen::VectorXd Geo::DT_verts_pos(int i1) {
    return points[i1].pos.cwiseProduct(scale);
  };

  int	Geo::DT_verts_nborverts_size(int i1) {
    //std::cout << "begin DT_verts_nborverts_size()" << std::endl;
    BP_GVec_Member<V_Element> *m1, *m2, *m3;
    BP_Vec<int> v;
    int i, j, k, ii;
    k = 0;
    m1 = points.member(i1);
    for(i = 0; i < DT_graph.num_incident_edges(m1); i++) {
      m2 = DT_graph.get_neighbor(m1, i);

      for(j = 0; j < DT_graph.num_incident_edges(m2); j++) {
        m3 = DT_graph.get_neighbor(m2, j);
        ii = points.get_index(m3);
        if(ii != -1)
          if(points[ii].cut == 0)
            if(m3 != m1) {
              v.add(ii);
              points[ii].cut = 1;
              k++;
            }
      }

    }

    for(i = 0; i < v.size(); i++)
      points[v[i]].cut = 0;

    //std::cout << "finish DT_verts_nborverts_size()" << std::endl;
    return k;
  };

  BP_Vec<int> Geo::DT_verts_nborverts(int i1) {
    //std::cout << "begin DT_verts_nborverts_size()" << std::endl;
    BP_GVec_Member<V_Element> *m1, *m2, *m3;
    int i, j, k, ii;
    BP_Vec<int> v;
    k = 0;
    m1 = points.member(i1);
    for(i = 0; i < DT_graph.num_incident_edges(m1); i++) {
      m2 = DT_graph.get_neighbor(m1, i);

      for(j = 0; j < DT_graph.num_incident_edges(m2); j++) {
        m3 = DT_graph.get_neighbor(m2, j);
        ii = points.get_index(m3);
        if(ii != -1)
          if(points[ii].cut == 0)
            if(m3 != m1) {
              points[ii].cut = 1;
              v.add(ii);
            }
      }

    }

    for(i = 0; i < v.size(); i++)
      points[v[i]].cut = 0;

    //std::cout << "finish DT_verts_nborverts_size()" << std::endl;
    return v;
  };

  int	Geo::DT_verts_nborsimps_size(int i1) {
    //std::cout << "begin DT_verts_nborsimps_size()" << std::endl;
    BP_GVec_Member<V_Element> *m1, *m2;
    //std::cout << "dim: " << dim << std::endl;
    int i, j, k;
    m1 = points.member(i1);
    k = 0;
    for(i = 0; i < DT_graph.num_incident_edges(m1); i++) {
      m2 = DT_graph.get_neighbor(m1, i);
      j = DT_simplexes.get_index(m2);
      if(j != -1) {
        k++;

      }
    }
    //std::cout << "finish DT_verts_nborsimps_size()" << std::endl;
    return k;
  };

  BP_Vec<int> Geo::DT_verts_nborsimps(int i1) {
    //std::cout << "begin DT_verts_nborsimps()" << std::endl;
    BP_GVec_Member<V_Element> *m1, *m2;
    BP_Vec<int> v;
    //std::cout << "dim: " << dim << std::endl;
    int i, j, k;
    m1 = points.member(i1);
    k = 0;
    for(i = 0; i < DT_graph.num_incident_edges(m1); i++) {
      m2 = DT_graph.get_neighbor(m1, i);
      j = DT_simplexes.get_index(m2);
      if(j != -1) {
        v.add(j);

      }
    }
    //std::cout << "finish DT_verts_nborsimps()" << std::endl;
    return v;
  };


  //////////////////////////////////////////
  // Voronoi Diagram

  bool Geo::calc_VD() {
    if(VD_exists) return 1;

    // Create VD by finding Delaunay Triangulation

    if(verbosity >= 10) {
      std::cout << "find Voronoi Diagram:" << std::endl;
      std::cout << "1) find Delaunay Triangulation" << std::endl;
    }
    if(calc_DT()) {
      if(verbosity >= 10)
        std::cout << "   DT complete" << std::endl;

      // now use the Delaunay Triangulation to set the VD_graph
      int i, j, k, ii, d, nie1;
      BP_Gen_GVec_Member *tmp_member, *element, *m1, *m2;
      BP_Gen_Vertex *cell, *n1;
      BP_Vec<int> curr_ccenters;
      BP_Vec<int> possible_nbor_points;
      BP_Vec<int> curr_points;
      d = get_dim();


      // create memory for saving faces
      VD_elements = new BP_GVec< V_Element>[d + 1];
      VD_exists = 1;

      VD_elements[0].capacity(DT_simps_size());
      for(i = 1; i <= d; i++) {
        VD_elements[i].capacity(size());
        VD_elements[i].cap_increment(size());
      }


      //	for d=dim:-1:0
      //		if( d == dim)
      //			faces of dimension dim are all the ccenters connected to each point
      //		else
      //			faces of dimension d<dim are a subset of ccenters connected to face of dimension d+1
      //
      //	VD_Elements are connected (in VD_graph) such that voronoi cells are connected to neighbors, other elements only connect to their lower dimension faces
      Eigen::VectorXd VZero(d);
      Eigen::VectorXd cc, ps;
      Eigen::VectorXd vertex_mean(d);
      VZero.setZero();
      BP_Vec<int> a;

      bool eq_ccenter_found;
      BP_Vec<int> cc_index_map;

      // add ccenters to VD_Elements[0], without allowing repeats
      //  should modify this to check only neighboring simps for repeats
      if(verbosity >= 10)
        std::cout << "add elements of dimension 0 (vertices) to VD_Elements" << std::endl;
      for(i = 0; i < DT_simps_size(); i++) {
        //std::cout << std::endl << "---------" << std::endl << "i: " << i << std::endl;
        cc = DT_simps_circumcenter(i).cwiseQuotient(scale);	// is scaled to full range, must divide by scale
        //std::cout << "cc: " << cc.transpose() << std::endl;

        eq_ccenter_found = 0;
        for(j = 0; j < VD_elements[0].size(); j++) {
          //std::cout << "ps: " << VD_elements[0][j].pos.transpose() ;
          if(vector_is(cc, VD_elements[0][j].pos, VD_tol)) {
            //std::cout << "  EQUIV" << std::endl;
            eq_ccenter_found = 1;
            cc_index_map.add(j);
            break;
          }
          else {
            //std::cout << "  NOTEQ" << std::endl;
          }
        }
        if(!eq_ccenter_found) {
          // V_Elment( rank, vertex_mean)
          cc_index_map.add(VD_elements[0].size());
          VD_graph.add_vertex(VD_elements[0].add(V_Element(0, cc)));

        }


      }

      //for( i=0; i<cc_index_map.size(); i++)
      //{
      //	std::cout << "i: " << i << "  cc_index_map: " << cc_index_map[i] << "  pos: " << VD_elements[0][ cc_index_map[i]].pos.transpose() << std::endl;
      //}
      //
      //BP_pause();

      // add cells to VD_Elements[ get_dim()];
      //  and connect cells with vertices
      if(verbosity >= 10)
        std::cout << "add elements of dimension " << d << " (Voronoi cells) to VD_Elements" << std::endl;
      for(j = 0; j < size(); j++) {
        //get neighboring simplexes to this point
        //  the circumcenters of the neighboring simplexes are the vertices of the voronoi cell
        a = DT_verts_nborsimps(j);

        curr_ccenters.clear();
        for(k = 0; k < a.size(); k++) {
          eq_ccenter_found = 0;
          for(ii = 0; ii < curr_ccenters.size(); ii++) {
            if(curr_ccenters[ii] == cc_index_map[ a[k]]) {
              eq_ccenter_found = 1;
              ii = curr_ccenters.size();
            }
          }
          if(eq_ccenter_found == 0)
            curr_ccenters.add(cc_index_map[ a[k]]);
        }

        //find the mean of the neighboring simplexes ccenters, this is the vertex mean for the cell
        vertex_mean.setZero();
        for(k = 0; k < curr_ccenters.size(); k++)
          vertex_mean += VD_elements[0][ curr_ccenters[k]].pos;
        //vertex_mean += DT_simps_circumcenter(a[k]);
        vertex_mean /= 1.0 * curr_ccenters.size();

        //std::cout << " -a.size(): " << a.size() << std::endl;
        //std::cout << " -vertex_mean: " << vertex_mean.transpose() << std::endl;

        // add the cell to the VD_graph
        // V_Element( rank, pos, vertex_mean)
        ps = pos(j);
        tmp_member = VD_elements[d].add(V_Element(d, ps, vertex_mean));
        VD_graph.add_vertex(tmp_member);

        // connect the cell to the vertices
        for(k = 0; k < curr_ccenters.size(); k++) {
          VD_graph.connect_vertices(tmp_member, VD_elements[0].member(curr_ccenters[k]), connections.add(0), 0);
        }
      }

      //std::cout << "check connections: " << std::endl;
      //for( j=0; j<size(); j++)
      //for( i=0; i<VD_graph.num_incident_edges( VD_elements[d].member(j)); i++)
      //{
      //	if( VD_graph.vert_val( VD_graph.get_neighbor( VD_elements[d].member(j), i)).rank == 0)
      //	{
      //		std::cout << "  pt: " << j << "  nbor: " << i << "  pos: " << VD_graph.vert_val( VD_graph.get_neighbor( VD_elements[d].member(j), i)).pos.transpose() << std::endl;
      //	}
      //}

      // now we want information on all the facets of the voronoi cells

      // for each cell, create lists of ccenters and nborverts and curr_points = this cell
      // then call VD_face_finder
      if(verbosity >= 10)
        std::cout << "add all the other elements" << std::endl;
      for(i = 0; i < size(); i++) {
        curr_ccenters.clear();
        possible_nbor_points.clear();
        curr_points.clear();

        element = VD_elements[d].member(i);

        curr_points.add(i);

        a = DT_verts_nborsimps(i);
        for(j = 0; j < a.size(); j++) {
          eq_ccenter_found = 0;
          for(k = 0; k < curr_ccenters.size(); k++)
            if(curr_ccenters[k] == cc_index_map[ a[j]]) {
              eq_ccenter_found = 1;
              k = curr_ccenters.size();
            }
          if(!eq_ccenter_found)
            curr_ccenters.add(cc_index_map[ a[j]]);
        }

        for(j = 0; j < curr_ccenters.size(); j++) {
          m1 = VD_elements[0].member(curr_ccenters[j]);
          nie1 = VD_graph.num_incident_edges(m1);
          for(k = 0; k < nie1; k++) {
            m2 = VD_graph.get_neighbor(m1, k);
            if(m2 != element)
              if(VD_graph.vert_val(m2).rank == d) {
                if(VD_graph.vert_val(m2).cut == 0) {
                  VD_graph.vert_val(m2).cut = 1;
                  possible_nbor_points.add(VD_elements[d].get_index(m2));
                }
              }

          }
        }

        //std::cout << std::endl << "%%%%%%" << std::endl << "cell: " << i << "  pos: " << pos(i).transpose() << " curr_ccenters: ";
        //for( j=0; j<curr_ccenters.size(); j++)
        //{
        //	std::cout << " " << curr_ccenters[j] ;
        //}
        //std::cout << std::endl ;
        //std::cout << "cell: " << i << " possible_nbor_points: ";
        for(j = 0; j < possible_nbor_points.size(); j++) {
          //std::cout << " " << possible_nbor_points[j] ;
          VD_elements[d][possible_nbor_points[j]].cut = 0;
        }
        //std::cout << std::endl;


        VD_face_finder(d, d, element, curr_ccenters, possible_nbor_points, curr_points);

      }

      if(verbosity >= 10)
        std::cout << "Voronoi diagram has:" << std::endl;
      for(i = 0; i < d + 1; i++) {
        for(j = 0; j < VD_elements[i].size(); j++)
          VD_elements[i][j].nvolume = -1;


        if(verbosity >= 10) {
          if(i == 0) std::cout << "  " << VD_elements[i].size() << " vertices" << std::endl;
          else if(i == d) std::cout << "  " << VD_elements[i].size() << " cells of dimension " << i << std::endl;
          else std::cout << "  " << VD_elements[i].size() << " elements of dimension " << i << std::endl;
        }
      }

      //std::cout << "Calculate Voronoi diagram element volumes" << std::endl;
      //VD_calc_data();

      if(verbosity >= 10)
        std::cout << "Voronoi diagram complete" << std::endl;

      //for( i=0; i<VD_elements_size(2); i++)
      //{
      //	std::cout << " face: " << i << "  vertex_mean: " << VD_elements_vertex_mean(2,i).transpose() << std::endl;
      //	a = VD_elements_verts(2,i);
      //	for( j=0; j<a.size(); j++)
      //	{
      //		std::cout << "  j: " << j << "  vertex: " << a[j] << "  pos: " << VD_elements_vertex_mean(0,a[j]).transpose() << std::endl;
      //	}
      //
      //}

    }
    else {
      if(verbosity >= 10)
        std::cout << "Delaunay triangulation failed!" << std::endl;
    }
    return VD_exists;
  };


  // pass overall dimension d, current element, list of ccenters, list of nborverts, list of points making up this element
  void Geo::VD_face_finder(int curr_d, int dim, BP_Gen_GVec_Member *element, BP_Vec<int> &curr_ccenters, BP_Vec<int> &possible_nbor_points, BP_Vec<int> &curr_points) {
    //std::cout << "begin VD_face_finder(" << curr_d << ")" << std::endl;
    int i, j, k, nie1, nie2;
    BP_Vec<int> tmp_face;
    BP_Vec<int> tmp_points;
    BP_Vec<int> tmp_nbor_points;
    BP_Vec<int> a;
    BP_Gen_GVec_Member *tmp_member, *m1, *m2, *m3;
    Eigen::VectorXd out_norm, pos, vertex_mean;

    //for each neighboring point
    for(i = 0; i < possible_nbor_points.size(); i++) {
      //std::cout << "  i: " << i << "  possible_nbor_point: " << possible_nbor_points[i] << std::endl;

      //only consider combinations of curr_points that are in increasing order
      //if( possible_nbor_points[i] > curr_points[curr_points.size()-1])
      {
        //std::cout << "    in order" << std::endl;

        // face is all of curr_ccenters that are also 1NN with possible_nbor_points[i]
        tmp_face.clear();
        m1 = VD_elements[dim].member(possible_nbor_points[i]);
        for(j = 0; j < curr_ccenters.size(); j++) {

          if(VD_graph.are_1NN(VD_elements[0].member(curr_ccenters[j]), m1)) {
            tmp_face.add(curr_ccenters[j]);
            VD_elements[0][curr_ccenters[j]].cut = 1;
          }
          /*
          //if( curr_ccenters[j] is neighbor with nborverts[i])
          a = DT_simps_nborverts(curr_ccenters[j]);
          for( k=0; k<a.size(); k++)
          	if( a[k] == possible_nbor_points[i])
          	{
          		tmp_face.add( curr_ccenters[j]);
          		k = a.size();
          	}
          */
        }

        //std::cout << "    tmp_face.size()" << tmp_face.size() << std::endl;

        if(tmp_face.size() >= curr_d) {
          // check if this face already exists
          bool face_exists = 0;
          bool is_equiv;

          //std::cout << std::endl << "--------" << std::endl << "tmp_face: ";
          //for( j=0; j<tmp_face.size(); j++)
          //	std::cout << tmp_face[j] << " " ;
          //std::cout << std::endl;

          // check neighbor's of first vertex that are also of rank curr_d-1
          m1 = VD_elements[0].member(tmp_face[0]);
          nie1 = VD_graph.num_incident_edges(m1);
          for(j = 0; j < nie1; j++) {
            m2 = VD_graph.get_neighbor(m1, j);
            if(VD_graph.vert_val(m2).rank == curr_d - 1) {
              is_equiv = 1;

              //check if any of m2's vertices are different than this facet's vertices
              nie2 = VD_graph.num_incident_edges(m2);

              //std::cout << "trial facet: ";
              //for( k=0; k<nie2; k++)
              //{
              //	m3 = VD_graph.get_neighbor(m2,k);
              //	if( VD_graph.vert_val(m3).rank == 0)
              //	{
              //		std::cout << VD_elements[0].get_index(m3) << " " ;
              //	}
              //}
              //std::cout << std::endl;

              for(k = 0; k < nie2; k++) {
                m3 = VD_graph.get_neighbor(m2, k);
                if(VD_graph.vert_val(m3).rank == 0) {
                  // if not one of this facet's vertices
                  if(VD_graph.vert_val(m3).cut == 0) {
                    // set is_equiv to false
                    is_equiv = 0;
                    k = nie2;
                  }
                }

              }

              // if m2 is equivalent, then the face already exists
              if(is_equiv) {
                j = nie1;
                face_exists = 1;
              }
            }
          }

          for(j = 0; j < tmp_face.size(); j++)
            VD_elements[0][tmp_face[j]].cut = 0;


          if(!face_exists) {

            // add new face

            //std::cout << "      add new face, d: " << curr_d-1 << std::endl;
            //std::cout << std::endl << "--------" << std::endl << "tmp_face: ";
            //std::cout << "         tmp_face: " ;
            //for( j=0; j<tmp_face.size(); j++)
            //	std::cout << tmp_face[j] << " " ;
            //std::cout << std::endl;

            tmp_member = VD_form_facet_from_ccenters(curr_d - 1, VD_elements[curr_d - 1], tmp_face);
            VD_graph.add_vertex(tmp_member);

            //connect to element
            //std::cout << "      connect to element" << std::endl;
            VD_graph.connect_vertices(element, tmp_member, connections.add(0), 0);

            // maybe also want to connect to all higher dimension elements?

            //connect new face to vertices
            //std::cout << "      connect to vertices" << std::endl;

            if(curr_d - 1 == 1 & tmp_face.size() > 2) {
              std::cout << "warning: tmp_face.size(): " << tmp_face.size() << " d: " << curr_d - 1 << std::endl;
              //BP_pause();
            }

            for(k = 0; k < tmp_face.size(); k++)
              VD_graph.connect_vertices(tmp_member, VD_elements[0].member(tmp_face[k]), connections.add(0), 0);

            if(curr_d - 1 > 1) {
              // find faces of new face

              //maybe reset the list of nborverts?
              //tmp_nborverts = curr_points + (points that are connected to tmp_face ccenters)

              // set the list of points that are equidistant to new face
              //std::cout << "      create tmp_points" << std::endl;
              tmp_points.clear();
              tmp_points = curr_points;
              tmp_points.add(possible_nbor_points[i]);

              tmp_nbor_points.clear();
              for(j = 0; j < possible_nbor_points.size(); j++)
                if(j != i)
                  tmp_nbor_points.add(possible_nbor_points[j]);

              //std::cout << "        find faces of new face" << std::endl;
              VD_face_finder(curr_d - 1, dim, tmp_member, tmp_face, tmp_nbor_points, tmp_points);
            }
          }
          else {
            //std::cout << "  do not add, already exists" << std::endl;
          }

        }
        else {
          // unmark vertices of the facet
          for(j = 0; j < tmp_face.size(); j++)
            VD_elements[0][tmp_face[j]].cut = 0;
        }

      }

    }

    //std::cout << "finish VD_face_finder()" << std::endl;
    return;
  };

  // Form a facet from the planar points (circumcenters) in point_list, add to VD_elements_list, return BP_Gen_GVec_Member* ptr to new facet
  BP_Gen_GVec_Member *Geo::VD_form_facet_from_ccenters(int rank, BP_GVec<V_Element> &facet_list, BP_Vec<int> &ccenter_list) {
    //std::cout << " -begin VD_form_facet_from_points(" << rank << ")" << std::endl;

    int i;
    Eigen::VectorXd facet_center;

    facet_center.resize(get_dim());
    facet_center.setZero();

    // add an element for each vertex
    for(i = 0; i < ccenter_list.size(); i++) {
      facet_center += VD_elements[0][ ccenter_list[i]].pos;
    }

    facet_center /= ccenter_list.size();
    //std::cout << " -facet_center: " << facet_center.transpose() << std::endl;

    BP_Gen_GVec_Member *member = VD_elements[rank].add(V_Element(rank, facet_center, facet_center));
    //std::cout << " -finish VD_form_facet_from_points()" << std::endl;

    return member;
  };

  // returns the number of VD elements of rank i
  int Geo::VD_elements_size(int i) {
    if(VD_exists == 0)
      calc_VD();

    return VD_elements[i].size();
  };

  // returns the volume of the i2-th VD element of rank i1
  double Geo::VD_elements_volume(int i1, int i2) {
    if(VD_exists == 0)
      calc_VD();
    //if( VD_elements[i1][i2].nvolume != -1)
    //	return VD_elements[i1][i2].nvolume;

    //std::cout << "begin VD_elements_volume()" << std::endl;

    int i, j, k, vert_count;
    double height;
    Eigen::MatrixXd m1, m2, tQ;
    Eigen::VectorXd ref, center, v;
    Eigen::FullPivHouseholderQR< Eigen::MatrixXd > tQR;
    BP_Vec<int> a;

    if(i1 == 0) {
      VD_elements[0][i1].nvolume = 1;

    }
    else if(i1 == get_dim()) {

      a = VD_elements_verts(i1, i2);

      //std::cout << " a.size(): " << a.size() << std::endl;
      if(a.size() <= i1) {
        std::cout << "too few verts for element.  dim: " << i1 << " vert_count: " << a.size() << std::endl;
        VD_elements[i1][i2].nvolume = -1;
        return VD_elements[i1][i2].nvolume;
      }

      m2.resize(i1, a.size());

      for(k = 0; k < a.size(); k++)
        //m2.col(k) = VD_elements[0][ a[k]].pos;	// is scaled to range 1
        m2.col(k) = VD_vertex_pos(a[k]);			// is scaled to full range

      Geo tGeo(m2);
      tGeo.set_Geo_tol(get_Geo_tol());
      tGeo.set_verbosity(0);
      tGeo.calc_CH();
      VD_elements[i1][i2].nvolume = tGeo.CH_volume();

    }
    else {
      a = VD_elements_verts(i1, i2);
      m1.resize(get_dim(), a.size() - 1);
      //std::cout << " a.size(): " << a.size() << std::endl;
      if(a.size() <= i1) {
        std::cout << "too few verts for element.  dim: " << i1 << " vert_count: " << a.size() << std::endl;
        VD_elements[i1][i2].nvolume = -1;
        return VD_elements[i1][i2].nvolume;
      }
      m2.resize(i1, a.size());

      //ref = VD_elements[0][ a[0]].pos;	// is scaled to range 1
      ref = VD_vertex_pos(a[0]);		// is scaled to full range
      for(k = 1; k < a.size(); k++) {
        //m1.col(k-1) = VD_elements[0][ a[k]].pos - ref;	// is scaled to range 1
        m1.col(k - 1) = VD_vertex_pos(a[k]) - ref;		// is scaled to full range
      }

      //std::cout << "m1: " << std::endl;
      //std::cout << m1 << std::endl;


      // get basis of plane
      tQR.compute(m1);
      tQ = tQR.matrixQ();
      //std::cout << "tQ: " << std::endl;
      //std::cout << tQ << std::endl;
      //std::cout << "  rank: " << tQR.rank() << std::endl;

      //if( tQR.rank() != i)
      //{
      //	std::cout << "  ERROR in VD_calc_data(), rank: " << tQR.rank() << "  i: " << i << std::endl;
      //	std::cout << "  R: " << std::endl;
      //	std::cout << tQR.matrixQR() << std::endl;
      //	BP_pause();
      //}


      for(k = 0; k < a.size(); k++)
        //m2.col(k) = tQ.leftCols(i1).transpose()*VD_elements[0][ a[k]].pos;	// is scaled to range 1
        m2.col(k) = tQ.leftCols(i1).transpose() * VD_vertex_pos(a[k]);			// is scaled to full range

      Geo tGeo(m2);
      tGeo.set_Geo_tol(get_Geo_tol());
      tGeo.set_verbosity(0);
      tGeo.calc_CH();
      VD_elements[i1][i2].nvolume = tGeo.CH_volume();

    }

    if(VD_elements[i1][i2].nvolume == 0) {
      std::cout << "error, VD volume == 0" << std::endl;
      std::cout << " i1: " << i1 << "   i2: " << i2 << std::endl;
      std::cout << "  m2: " << std::endl;
      std::cout << m2.transpose() << std::endl;
      std::cout << " a.size(): " << a.size() << std::endl;
      std::cout << "  verts: " << std::endl;
      for(k = 0; k < a.size(); k++) {
        std::cout << "   k: " << k << "  pos: " << VD_vertex_pos(a[k]).transpose() << std::endl;
      }
      /*
      std::cout << "  m1: " << std::endl;
      std::cout << m1 << std::endl;
      std::cout << "   Q: " << std::endl;
      std::cout << tQ << std::endl;
      std::cout << "   R: " << std::endl;
      std::cout << tQR.matrixQR() << std::endl;
      std::cout << "rank: " << tQR.rank() << std::endl;
      */
      //abort();
    }

    //std::cout << "finish VD_elements_volume()" << std::endl;


    return VD_elements[i1][i2].nvolume;
  };

  // returns the position of the i1 vertex
  Eigen::VectorXd Geo::VD_vertex_pos(int i1) {
    if(VD_exists == 0)
      calc_VD();

    return VD_elements[0][i1].vertex_mean.cwiseProduct(scale);
  };

  // returns the position of the vertex mean of the i2-th element of rank i1
  Eigen::VectorXd Geo::VD_elements_vertex_mean(int i, int j) {
    if(VD_exists == 0)
      calc_VD();

    return VD_elements[i][j].vertex_mean.cwiseProduct(scale);
  };

  // returns the vertices defining the i2-th element of rank i1
  BP_Vec<int> Geo::VD_elements_verts(int i1, int i2) {
    if(VD_exists == 0)
      calc_VD();

    //std::cout << "begin DT_simps_nborsimps_size()" << std::endl;
    BP_GVec_Member<V_Element> *m1, *m2;
    BP_Vec<int> v;
    //std::cout << "dim: " << dim << std::endl;
    int i, j;
    m1 = VD_elements[i1].member(i2);
    for(i = 0; i < VD_graph.num_incident_edges(m1); i++) {
      m2 = VD_graph.get_neighbor(m1, i);
      j = VD_elements[0].get_index(m2);
      if(j != -1) {
        v.add(j);

      }
    }
    //std::cout << "finish DT_simps_nborsimps_size()" << std::endl;
    return v;
  };


};

#endif  //  BP_Geo_H

