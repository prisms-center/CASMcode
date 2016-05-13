#include "casm/crystallography/Lattice.hh"

#include "casm/crystallography/SupercellEnumerator.hh"
#include "casm/symmetry/SymOp.hh"

namespace CASM {

  Lattice::Lattice(const Eigen::Vector3d &vec1,
                   const Eigen::Vector3d &vec2,
                   const Eigen::Vector3d &vec3) {
    m_lat_mat << vec1, vec2, vec3;
    m_inv_lat_mat = m_lat_mat.inverse();
  }

  //********************************************************************

  ///Construct Lattice from a matrix of lattice vectors, where lattice vectors are columns
  ///(e.g., lat_mat is equivalent to lat_column_mat())
  Lattice::Lattice(const Eigen::Ref<const Eigen::Matrix3d> &lat_mat) :
    m_lat_mat(lat_mat),
    m_inv_lat_mat(lat_mat.inverse()) {
  }

  //********************************************************************

  Lattice Lattice::fcc() {
    Eigen::Matrix3d latmat;
    latmat <<
           0, 1, 1,
           1, 0, 1,
           1, 1, 0;
    latmat /= pow(latmat.determinant(), 1.0 / 3.0);
    return Lattice(latmat);
  }

  //********************************************************************

  Lattice Lattice::bcc() {
    Eigen::Matrix3d latmat;
    latmat <<
           -1, 1, 1,
           1, -1, 1,
           1, 1, -1;
    latmat /= pow(latmat.determinant(), 1.0 / 3.0);
    return Lattice(latmat);
  }

  //********************************************************************

  Lattice Lattice::cubic() {
    return Lattice(Eigen::Matrix3d::Identity());
  }

  //********************************************************************

  Lattice Lattice::hexagonal() {
    Eigen::Matrix3d latmat;
    latmat <<
           1, -1.0 / sqrt(3.0), 0,
           1, 1.0 / sqrt(3.0),  0,
           0, 0, sqrt(3.0) / 2.0;

    return Lattice(latmat.transpose());
  }

  //********************************************************************

  Lattice Lattice::scaled_lattice(double scale) const {
    return Lattice(scale * lat_column_mat());
  }

  //********************************************************************
  //Calculate length of lattice vector 'i'
  double Lattice::length(Index i) const {
    // Calculates Lengths
    return m_lat_mat.col(i).norm();
  }

  //********************************************************************
  double Lattice::angle(Index i) const {
    double t_a = m_lat_mat.col((i + 1) % 3).dot(m_lat_mat.col((i + 2) % 3)) / (length((i + 1) % 3) * length((i + 2) % 3));
    //Make sure that cos(angle) is between 0 and 1
    if((t_a - 1.0) > 0.0) {
      t_a = 1.0;
    }

    if((t_a + 1.0) < 0.0) {
      t_a = -1.0;
    }

    return (180.0 / M_PI) * acos(t_a);
  }

  //********************************************************************

  void Lattice::read(std::istream &stream) {
    double scale;
    stream >> scale;
    stream >> m_lat_mat;
    m_lat_mat *= scale;
    m_lat_mat.transposeInPlace();
    m_inv_lat_mat = m_lat_mat.inverse();

    return;
  }

  //********************************************************************

  void Lattice::print(std::ostream &stream, int _prec) const {
    int tprec = stream.precision();
    std::ios::fmtflags tflags = stream.flags();
    stream.precision(_prec);
    stream.width(_prec + 3);
    stream.flags(std::ios::showpoint | std::ios::fixed | std::ios::right);
    stream  << 1.0 << '\n';

    stream << ' ' << std::setw(_prec + 8) << m_lat_mat.col(0).transpose() << '\n';
    stream << ' ' << std::setw(_prec + 8) << m_lat_mat.col(1).transpose() << '\n';
    stream << ' ' << std::setw(_prec + 8) << m_lat_mat.col(2).transpose() << '\n';

    stream.precision(tprec);
    stream.flags(tflags);

    return;
  }

  //********************************************************************
  //Gets the reciprocal lattice from the lattice vectors... (AAB)
  Lattice Lattice::get_reciprocal() const {
    /* Old Expression
       return Lattice(2 * M_PI * cross_prod(vecs[1], vecs[2])/vol,
       2 * M_PI * cross_prod(vecs[2], vecs[0]) / vol,
       2 * M_PI * cross_prod(vecs[0], vecs[1]) / vol);
       return recip_lat;
    */
    return Lattice(2 * M_PI * inv_lat_column_mat().transpose()); //equivalent expression
  }

  //********************************************************************

  Array<int> Lattice::calc_kpoints(Array<int> prim_kpoints, Lattice prim_lat) {
    Array<int> super_kpoints = prim_kpoints;
    //    Lattice prim_recip_lat = (lattice.primitive->get_reciprocal());
    Lattice prim_recip_lat = (prim_lat.get_reciprocal());
    Lattice recip_lat = (*this).get_reciprocal();
    double prim_density = (prim_kpoints[0] * prim_kpoints[1] * prim_kpoints[2]) / (prim_recip_lat.vol());
    double super_density = 0;

    //        std::cout << prim_recip_lat << std::endl;
    //std::cout << recip_lat << std::endl;

    Array<double> prim_vec_lengths;

    for(int i = 0; i < 3; i++) {
      prim_vec_lengths.push_back(prim_recip_lat.length(i));
    }
    //std::cout<<prim_vec_lengths<<"\n";
    double shortest = prim_vec_lengths.min();
    int short_ind = prim_vec_lengths.find(shortest);

    //std::cout << "prim kpoints "<< prim_kpoints << std::endl;
    //    std::cout << "prim recip vol " << prim_recip_lat.vol << std::endl;
    //std::cout << "prim kpoint density " << prim_density << std::endl;

    double scale = (prim_kpoints[short_ind] / shortest);

    for(int i = 0; i < 3; i++) {
      super_kpoints[i] = int(ceil(scale * recip_lat.length(i)));
    }

    super_density = (super_kpoints[0] * super_kpoints[1] * super_kpoints[2]) / (recip_lat.vol());

    //std::cout << "super kpoint density " << super_density << std::endl;
    //std::cout << super_kpoints << std::endl;

    while(super_density < prim_density) {
      //  std::cout << prim_kpoints[short_ind] << std::endl;
      prim_kpoints[short_ind]++;
      //std::cout << prim_kpoints[short_ind] << std::endl;
      scale = (prim_kpoints[short_ind] / shortest);
      //std::cout << scale << std::endl;

      for(int i = 0; i < 3; i++) {
        super_kpoints[i] = int(ceil(scale * recip_lat.length(i)));
      }

      //std::cout << super_kpoints << std::endl;

      super_density = (super_kpoints[0] * super_kpoints[1] * super_kpoints[2]) / (recip_lat.vol());

      //std::cout << super_density << std::endl;
    }
    //std::cout << "---------------\n";

    return super_kpoints;
  }

  //********************************************************************
  // Finds subgroup of super_group that leaves this lattice invariant.
  // note that this routine is written so that pg_tol has exactly the same meaning
  // as in Lattice::generate_point_group
  void Lattice::find_invariant_subgroup(const SymGroup &super_group, SymGroup &sub_group, double pg_tol) const {
    if(sub_group.size() != 0) {
      std::cerr << "WARNING in Lattice::find_invariant_group" << std::endl;
      std::cerr << "The subgroup isn't empty and it's about to be rewritten!" << std::endl;
      sub_group.clear();
    }
    Eigen::Matrix3d tfrac_op, tMat;
    for(Index ng = 0; ng < super_group.size(); ng++) {
      tfrac_op = lat_column_mat().inverse() * super_group[ng].matrix() * lat_column_mat();

      //Use a soft tolerance of 1% to see if further screening should be performed
      if(!almost_equal(1.0, std::abs(tfrac_op.determinant()), 0.01) || !is_integer(tfrac_op, 0.01)) {
        continue;
      }

      //make tfrac_op integer.
      for(int i = 0; i < 3; i++) {
        for(int j = 0; j < 3; j++) {
          tfrac_op(i, j) = round(tfrac_op(i, j));
        }
      }

      // If symmetry is perfect, then ->  cart_op * lat_column_mat() == lat_column_mat() * frac_op  by definition
      // If we assum symmetry is imperfect, then ->   cart_op * lat_column_mat() == F * lat_column_mat() * frac_op
      // where 'F' is the displacement gradient tensor imposed by frac_op

      // tMat uses some matrix math to get F.transpose()*F*lat_column_mat();
      tMat = inv_lat_column_mat().transpose() * (tfrac_op.transpose() * lat_column_mat().transpose() * lat_column_mat() * tfrac_op);

      // Subtract lat_column_mat() from tMat, leaving us with (F.transpose()*F - Identity)*lat_column_mat().
      // This is 2*E*lat_column_mat(), where E is the green-lagrange strain
      tMat = (tMat - lat_column_mat()) / 2.0;

      //... and then multiplying by the transpose...
      tMat = tMat * tMat.transpose();

      // The diagonal elements of tMat describe the square of the distance by which the transformed vectors 'miss' the original vectors
      if(tMat(0, 0) < pg_tol * pg_tol && tMat(1, 1) < pg_tol * pg_tol && tMat(2, 2) < pg_tol * pg_tol) {
        sub_group.push_back(super_group[ng]);
      }

    }

  }
  //********************************************************************

  void Lattice::generate_point_group(SymGroup &point_group, double pg_tol) const {

    if(point_group.size() != 0) {
      std::cerr << "WARNING in Lattice::generate_point_group" << std::endl;
      std::cerr << "The point group for your lattice isn't empty and it's about to be rewritten!" << std::endl;
      point_group.clear();
    }

    point_group.set_lattice(*this);

    //Enumerate all possible matrices with elements equal to -1, 0, or 1
    //These represent operations that reorder lattice vectors or replace one
    //or more lattice vectors with a face or body diagonal.
    Eigen::Matrix3d tMat, tOp_cart;
    EigenCounter<Eigen::Matrix3i> pg_count(Eigen::Matrix3i::Constant(-1),
                                           Eigen::Matrix3i::Constant(1),
                                           Eigen::Matrix3i::Constant(1));

    //For this algorithm to work, lattice needs to be in reduced form.
    Lattice tlat_reduced(get_reduced_cell());
    do {

      //continue if determinant is not 1, because it doesn't preserve volume
      if(std::abs(pg_count().determinant()) != 1) {
        continue;
      }

      tOp_cart = tlat_reduced.lat_column_mat() * pg_count().cast<double>() * tlat_reduced.inv_lat_column_mat();

      //Find the effect of applying symmetry to the lattice vectors
      //The following is equivalent to point_group[i].matrix().transpose()*tlat_reduced.lat_column_mat()*point_group[i].matrix()
      tMat = tOp_cart.transpose() * tlat_reduced.lat_column_mat() * pg_count().cast<double>();

      //If pg_count() is a point_group operation, tMat should be equal to tlat_reduced.lat_column_mat().  We check by first taking the difference...
      tMat = (tMat - tlat_reduced.lat_column_mat()) / 2.0;

      //... and then multiplying by the transpose...
      tMat = tMat * tMat.transpose();

      // The diagonal elements are square of the distances by which the transformed lattice vectors "miss" the original lattice vectors
      // If they are less than the square of the tolerance, we add the operation to the point group
      if(tMat(0, 0) < pg_tol * pg_tol && tMat(1, 1) < pg_tol * pg_tol && tMat(2, 2) < pg_tol * pg_tol) {
        point_group.push_back(SymOp(tOp_cart, sqrt(tMat.diagonal().maxCoeff())));
      }

    }
    while(++pg_count);

    if(!point_group.is_group(pg_tol)) {
      std::cerr << "*** WARNING *** \n"
                << "    Lattice::generate_point_group() has been called on an ill-conditioned lattice \n"
                << "    (i.e., a well-defined point group could not be found with the supplied tolerance of " << pg_tol << ").\n"
                << "    CASM will use the group closure of the symmetry operations that were found.  Please consider using the \n"
                << "    CASM symmetrization tool on your input files.\n";
      std::cout << "Lat_column_mat:\n" << lat_column_mat() << "\n\n";

      point_group.enforce_group(pg_tol);

    }
    //Sort point_group by trace/conjugacy class
    point_group.sort_by_class();

    return;
  }

  //********************************************************************

  Array<double> Lattice::pg_converge(double large_tol) {
    Array<double> tarray;
    SymGroup point_group;
    generate_point_group(point_group, large_tol);
    if(!point_group.is_group(large_tol)) {
      std::cout << "This is not a group. It is being enforced...\n";
      point_group.enforce_group(large_tol);
    }
    for(Index i = 0; i < point_group.size(); i++) {
      tarray.push_back(point_group[i].map_error());
    }
    return tarray;
  }



  //********************************************************************

  void Lattice::pg_converge(double small_tol, double large_tol, double increment) {
    Array<double> tols;
    Array<bool> is_group, is_group_now;
    Array<int> num_ops, num_enforced_ops;
    Array<std::string> old_name, new_name;

    for(double i = small_tol; i <= large_tol; i += increment) {
      SymGroup point_group;

      tols.push_back(i);
      generate_point_group(point_group, i);
      point_group.character_table();
      old_name.push_back(point_group.get_name());
      num_ops.push_back(point_group.size());
      is_group.push_back(point_group.is_group(i));
      point_group.enforce_group(i);
      is_group_now.push_back(point_group.is_group(i));
      num_enforced_ops.push_back(point_group.size());
      point_group.character_table();
      new_name.push_back(point_group.get_name());
    }

    for(Index i = 0; i < tols.size(); i++) {
      std::cout << tols[i] << "\t" << num_ops[i] << "\t" << is_group[i] << "\t" << old_name[i] << "\t" << num_enforced_ops[i] << "\t" << is_group_now[i] << "\t" << new_name[i] << "\n";
    }

    return;
  }


  /// \brief Generate super Lattice
  ///
  /// Use SupercellEnumerator to enumerate possible HNF transformation matrices. Unique supercells
  /// are identified by applying point group operations and keeping the supercell if the HNF is 'canonical',
  /// meaning that the HNF indices in order H00, H11, H22, H12, H02, H01 are the lexicographically greatest.
  ///
  /// The supercell that is inserted in the 'supercell' container is the niggli cell, rotated to a
  /// standard orientation (see standard_orientation function).
  ///
  void Lattice::generate_supercells(Array<Lattice> &supercell,
                                    const SymGroup &effective_pg,
                                    int max_prim_vol,
                                    int min_prim_vol) const {
    std::cout << "Inside generate_supercells" << std::endl;
    SupercellEnumerator<Lattice> enumerator(*this, effective_pg, min_prim_vol, max_prim_vol + 1);
    supercell.clear();
    for(auto it = enumerator.begin(); it != enumerator.end(); ++it) {
      supercell.push_back(niggli(*it, effective_pg, TOL));
    }
    return;
  }


  //********************************************************************
  /**This function finds the reduced cell from the given primitive cell.
   *
   *
   * First, the translation vectors of the primitive cell are determined, then reduced.
   * The Bravais lattice is determined using Niggli's transformations. The method used to
   * determine the basis vectors are to find body and face diagonals across the cell.
   * The vectors satisfying these main conditions can be found by ensuring:
   * 		 i) the shortest face diagonal is not shorter than the longest edge of the same face
   * 		 ii) the shortest body diagonal is not shorter than the longest edge of the cell
   *
   * This is accomplished by enumerating the skew matrices (identity matrix with one of the
   * off-diagonal terms =1 or =-1).  This will replace one of the lattice vectors with the
   * face diagonal corresponding to that linear combination.  If the length of any of the
   * new skewed lattice vectors is shorter than the length of the original cell, it is replaced.
   * All skew matrices have determinants equal to 1, which preserves the volume of the original
   * cell.
   *
   */
  Lattice Lattice::get_reduced_cell() const {

    int i, j, k, nv;
    Array<Eigen::Matrix3d > skew;
    Eigen::Matrix3d tskew(Eigen::Matrix3d::Identity());
    bool minimized = false;

    //std::cout << "Before reduction: \n";
    //print(std::cout);


    //Creates 12 skew matrices
    //Do we also need the 12 "double skew" matrices corresponding to body diagonal substitution? YES
    skew.reserve(24);
    for(i = 0; i < 3; i++) {
      for(j = 0; j < 3; j++) {
        if(i != j) { //Checks to make sure we're not at a diagonal
          tskew(i, j) = 1;
          skew.push_back(tskew);
          tskew(i, j) = -1;
          skew.push_back(tskew);
          tskew(i, j) = 0;
        }
      }
    }

    //the 12 "double skew" matrices corresponding to body diagonal substitution
    for(i = 0; i < 3; i++) {		// column

      j = (i + 1) % 3;
      k = (i + 2) % 3;

      tskew(j, i) = 1;
      tskew(k, i) = 1;
      skew.push_back(tskew);

      tskew(j, i) = 1;
      tskew(k, i) = -1;
      skew.push_back(tskew);

      tskew(j, i) = -1;
      tskew(k, i) = 1;
      skew.push_back(tskew);

      tskew(j, i) = -1;
      tskew(k, i) = -1;
      skew.push_back(tskew);

      tskew(j, i) = 0;
      tskew(k, i) = 0;

    }


    Eigen::Matrix3d tMat, reduced_lat_mat;
    Eigen::Matrix3d tMat_mags, reduced_mags;
    reduced_lat_mat = lat_column_mat();  //Matrix3
    reduced_mags = reduced_lat_mat.transpose() * reduced_lat_mat;

    while(!minimized) {
      minimized = true;
      for(Index s = 0; s < skew.size(); s++) {

        //Multiply the original lattice (columns are lattice vectors) with the skew matrix
        //Right multiply does elementary column operation on lattice vectors
        tMat = (reduced_lat_mat * skew[s]);

        //Lattice vectors times their transpose give length and angle info
        tMat_mags = tMat.transpose() * tMat;

        for(nv = 0; nv < 3; nv++) {
          if(tMat_mags(nv, nv) < reduced_mags(nv, nv)) {
            reduced_lat_mat = tMat;
            reduced_mags = reduced_lat_mat.transpose() * reduced_lat_mat;
            minimized = false;
            s--; //try the same skew operator again.
            break;
          }
        }
      }
    }

    Lattice reduced_lat(reduced_lat_mat);
    return reduced_lat;
  }

  //********************************************************************

  Eigen::Vector3d Lattice::max_voronoi_vector(const Eigen::Vector3d &pos)const {

    double tproj(-1), maxproj(-1);
    int maxv;
    if(!voronoi_table.size()) {
      generate_voronoi_table();
    }

    for(Index nv = 0; nv < voronoi_table.size(); nv++) {
      tproj = pos.dot(voronoi_table[nv]);
      if(tproj > maxproj) {
        maxproj = tproj;
        maxv = nv;
      }
    }


    return voronoi_table[maxv];

  }

  //********************************************************************
  int Lattice::voronoi_number(const Eigen::Vector3d &pos)const {

    int tnum = 0;
    double tproj = 0;

    if(!voronoi_table.size()) {
      generate_voronoi_table();
    }

    for(Index nv = 0; nv < voronoi_table.size(); nv++) {
      tproj = pos.dot(voronoi_table[nv]);
      if(almost_equal(tproj, 1.0) < TOL) {
        tnum++;
      }
      else if(tproj > 1.0) {
        return -1;
      }
    }


    return tnum;
  }
  //********************************************************************

  void Lattice::generate_voronoi_table() const {
    voronoi_table.clear();
    //There are no fewer than 12 points in the voronoi table
    voronoi_table.reserve(12);

    Eigen::Vector3d tpoint;
    int i;

    Lattice tlat_reduced(get_reduced_cell());
    //Count over all lattice vectors, face diagonals, and body diagonals
    //originating from origin;
    EigenCounter<Eigen::Vector3i > combo_count(Eigen::Vector3i(-1, -1, -1),
                                               Eigen::Vector3i(1, 1, 1),
                                               Eigen::Vector3i(1, 1, 1));

    //std::cout << "For angles " << angles[0] << ", " << angles[1] << ", " << angles[2] << ", Voronoi table is: \n";
    //For each linear combination, check to see if it is on a face, edge, or vertex of the voronoi cell
    do {


      if(!(combo_count[0] || combo_count[1] || combo_count[2])) {
        continue;
      }

      //A linear combination does not fall on the voronoi boundary if the angle between
      //any two of the vectors forming that combination are obtuse
      for(i = 0; i < 3; i++) {
        if((combo_count[(i + 1) % 3] && combo_count[(i + 2) % 3])
           && std::abs(90.0 * abs(combo_count[(i + 1) % 3] - combo_count[(i + 2) % 3]) - angle(i)) + TOL < 90) {
          break;
        }
      }

      if(i == 3) {
        tpoint = tlat_reduced.lat_column_mat() * combo_count().cast<double>();
        double tval(tpoint.norm());

        tpoint /= tval * tval / 2;
        voronoi_table.push_back(tpoint);
        //std::cout << voronoi_table.back();
      }

    }
    while(++combo_count);
    return;
  }

  //********************************************************************
  /**
   * This function, given a linearly independent set of lattice vectors, finds the
   * dimensions along the unit cell vectors such that a sphere of given radius
   * fits within a uniform grid of 2dim[1]x2dim[2]x2dim[3] lattice points
   * centered at the origin.
   *
   * The algorithm works by getting the normal (e.g. n1) to each pair of lattice
   * vectors (e.g. a2, a3), scaling this normal to have length radius and
   * then projecting this normal parallel to the a2,a3 plane onto the
   * remaining lattice vector a1. This will tell us the number of a1 vectors
   * needed to make a grid to encompass the sphere.
   */
  //********************************************************************
  Eigen::Vector3i Lattice::enclose_sphere(double radius) const {

    // reciprocal vectors
    Eigen::Matrix3d recip(inv_lat_column_mat().transpose());

    for(int i = 0; i < 3; i++) {
      recip.col(i) *= radius / recip.col(i).norm();
    }
    //recip contains three column vectors of length 'radius' pointed along plane normals

    return (inv_lat_column_mat() * recip).cwiseAbs().unaryExpr(std::ptr_fun(ceil)).colwise().maxCoeff().cast<int>();
  }

  //********************************************************************
  //Change bool to an array of SymOps you want to use, default point group
  //Overload to only use identity
  //Return N matrix
  bool Lattice::is_supercell_of(const Lattice &tile, const Array<SymOp> &symoplist, Eigen::Matrix3d &multimat, double _tol) const {
    auto result = is_supercell(*this, tile, symoplist.begin(), symoplist.end(), _tol);
    multimat = result.second.cast<double>();
    return result.first != symoplist.end();
  }

  //********************************************************************
  //Change bool to an array of SymOps you want to use, default point group
  //Overload to only use identity
  //Return N matrix
  bool Lattice::is_supercell_of(const Lattice &tile, Eigen::Matrix3d &multimat, double _tol) const {
    auto result = is_supercell(*this, tile, _tol);
    multimat = result.second.cast<double>();
    return result.first;
  }

  //********************************************************************
  //Overload to only use identity
  //Return N matrix
  bool Lattice::is_supercell_of(const Lattice &tile, double _tol) const {
    return is_supercell(*this, tile, _tol).first;
  }

  //********************************************************************

  bool Lattice::is_supercell_of(const Lattice &tile, const Array<SymOp> &symoplist, double _tol) const {
    return is_supercell(*this, tile, symoplist.begin(), symoplist.end(), _tol).first != symoplist.end();
  }

  //********************************************************************

  Lattice &Lattice::make_right_handed() {

    if(lat_column_mat().determinant() < 0) {
      m_lat_mat = -m_lat_mat;
      m_inv_lat_mat = m_lat_mat.inverse();
    }

    return *this;
  }

  //\John G 121212
  //********************************************************************************************************

  Eigen::Vector3i Lattice::get_millers(Eigen::Vector3d plane_normal, double tolerance) const {
    //Get fractional coordinates of plane_normal in recip_lattice
    //These are h, k, l
    //For miller indeces h, k and l    plane_normal[CART]=h*a.recip+k*b.recip+l*c.recip
    return scale_to_int(lat_column_mat().transpose() * plane_normal, tolerance);
  }

  //John G 121015
  //********************************************************************
  /**
   *  Using miller indeces, the intercept of the plane is calculated. All intercepts
   *  are multiplied by a factor such that their values are integers. This corresponds
   *  to the plane being shifted so that it lands on three lattice sites. Using these
   *  lattice sites new vectors A and B that lie on the desired plane can be constructed.
   *  The choice of vector C is somewhat arbitrary. get_lattice_in_plane will first
   *  construct C to be the shortest most orthogonal vector possible. The user is then
   *  given the option to make C more orthogonal to A and B in exchange for a greater
   *  length.
   *  The returned lattice will be the one containing the most orthogonal C
   *  that is still smaller than max_vol
   */
  //********************************************************************

  Lattice Lattice::get_lattice_in_plane(Eigen::Vector3i millers, int max_vol) const {  //John G 121030
    //Hold new lattice vectors in these. Then at the end we make an actual Lattice out of it
    Eigen::Matrix3d surface_cell, last_surface_cell;    //Holds new lattice vectors, two of which are in the surface plane

    //Miller indeces of 100, 010 or 001 mean you don't need a new cell to expose the plane, however
    //you may want to reorient the vectors so that the ab plane is exposed (useful for Structure::stitch)

    if(millers == Eigen::Vector3i(0, 1, 0)) {
      std::cout << "No chopping neccesary" << std::endl;
      std::cout << "Flipping your vectors to bring a and b into plane:" << std::endl;
      std::cout << "b --> c" << std::endl;
      std::cout << "a --> b" << std::endl;
      std::cout << "c --> a" << std::endl;
      return Lattice(lat_column_mat().col(2), lat_column_mat().col(0), lat_column_mat().col(1));
    }

    else if(millers == Eigen::Vector3i(1, 0, 0)) {
      std::cout << "No chopping neccesary" << std::endl;
      std::cout << "Flipping your vectors to bring a and b into plane:" << std::endl;
      std::cout << "a --> c" << std::endl;
      std::cout << "b --> a" << std::endl;
      std::cout << "c --> b" << std::endl;
      return Lattice(lat_column_mat().col(1), lat_column_mat().col(2), lat_column_mat().col(0));
    }

    else if(millers == Eigen::Vector3i(0, 0, 1)) { // || millers==Eigen::Vector3i(0,1,0) || millers==Eigen::Vector3i(1,0,0))
      return *this;
    }

    //Miller indeces of xx0, 0xx or x0x mean one of the lattice vectors is the same as
    //one of the primitive lattice vectors. This means we only need to get TWO
    //points that are on the plane, which we connect to get the second lattice vector
    else if(millers[0] == 0 || millers[1] == 0 || millers[2] == 0) {
      int zero = 0;

      //Find out which miller index is 0
      while(true) {
        if(millers[zero] == 0) {
          break;
        }

        else {
          zero++;
        }
      }

      surface_cell.col(0) = lat_column_mat().col(zero);

      Eigen::Vector3d H_miller_point, K_miller_point;
      Eigen::Vector3d HK;
      Eigen::Vector3d millers_dubs;

      //Turn integer millers into doubles for mathematical purposes (inverse)
      millers_dubs = millers.cast<double>();

      //std::cout<<millers_dubs<<std::endl;

      Eigen::Vector3d inv_miller_dubs;
      Eigen::Vector3i inv_miller;


      //In actualility, the inverse miller of 0 is infinity, we set it to 0 here so that we can use
      //scale_to_int without going crazy. Since we won't be using this inverse miller index it's OK
      inv_miller_dubs[zero] = 0;
      inv_miller_dubs[(zero + 1) % 3] = 1.0 / millers_dubs[(zero + 1) % 3];
      inv_miller_dubs[(zero + 2) % 3] = 1.0 / millers_dubs[(zero + 2) % 3];

      inv_miller = scale_to_int(inv_miller_dubs, TOL);
      H_miller_point = inv_miller[(zero + 1) % 3] * lat_column_mat().col((zero + 1) % 3);
      K_miller_point = inv_miller[(zero + 2) % 3] * lat_column_mat().col((zero + 2) % 3);

      std::cout << "inv millers dubs: " << inv_miller_dubs << std::endl;
      std::cout << "inv millers : " << inv_miller << std::endl;

      HK = K_miller_point - H_miller_point;
      surface_cell.col(1) = HK;
    }



    else {
      //Get three points that lie on the plane
      //We'll want to find points that lie on the plane AND land of lattice points. In order to do
      //this we need the miller inverses multiplied by a factor that makes them integers
      Eigen::Vector3d H_miller_point, K_miller_point, L_miller_point;
      Eigen::Vector3d inv_miller_dubs;
      Eigen::Vector3i inv_miller;
      //Turn integer millers into doubles for mathematical purposes (inverse)
      Eigen::Vector3d millers_dubs;
      millers_dubs = millers.cast<double>();

      inv_miller_dubs[0] = 1.0 / millers_dubs[0];
      inv_miller_dubs[1] = 1.0 / millers_dubs[1];
      inv_miller_dubs[2] = 1.0 / millers_dubs[2];

      inv_miller = scale_to_int(inv_miller_dubs, TOL);

      H_miller_point = inv_miller[0] * lat_column_mat().col(0);
      K_miller_point = inv_miller[1] * lat_column_mat().col(1);
      L_miller_point = inv_miller[2] * lat_column_mat().col(2);

      //Get three vectors that connect the three points on the plane to each other. Any two of the following
      //vectors could be used for constructing the new lattice, but it's convenient to pick the two
      //most orthogonal vectors
      Eigen::Vector3d HK, KL, LH;
      Eigen::Vector3d tangles;


      HK = K_miller_point - H_miller_point;
      KL = L_miller_point - K_miller_point;
      LH = H_miller_point - L_miller_point;

      //John G 121212
      //The vectors that we got at this point are valid, but sometimes larger than they need to be.
      Eigen::Matrix3d templat;
      templat.col(0) = HK;
      templat.col(1) = KL;
      templat.col(2) = LH;

      //Find shortest vector
      int s = 0;
      for(int i = 1; i < 3; i++) {
        if(lat_column_mat().col(i).norm() < lat_column_mat().col(s).norm()) {
          s = i;
        }
      }

      //Try dividing by integers and see if they're still lattice combinations. If they are
      //shorten and continue
      for(int i = 0; i < 3; i++) {
        int maxval = round(templat.col(i).norm() / lat_column_mat().col(i).norm() + 1);
        for(int j = maxval; j > 1; j--) {
          Eigen::Vector3d shortened = inv_lat_column_mat() * templat.col(i) / j; //Shorten and convert to fractional
          bool combo = true;

          for(int k = 0; k < 3; k++) {
            if(!almost_zero(std::abs(round(shortened[k]) - shortened[k]))) {
              combo = false;
              break;
            }
          }

          if(combo) {
            templat.col(i) = lat_column_mat() * shortened;
            break;
          }
        }
      }

      HK = templat.col(0);
      KL = templat.col(1);
      LH = templat.col(2);

      //We select the two vectors that spawn the smallest area

      double HKKL, KLLH, LHHK;
      HKKL = HK.cross(KL).norm();
      KLLH = KL.cross(LH).norm();
      LHHK = LH.cross(HK).norm();

      if(HKKL <= KLLH && HKKL <= LHHK) {
        surface_cell.col(0) = HK;
        surface_cell.col(1) = KL;
      }

      else if(KLLH <= HKKL && KLLH <= LHHK) {
        surface_cell.col(0) = KL;
        surface_cell.col(1) = LH;
      }

      else {
        surface_cell.col(0) = LH;
        surface_cell.col(1) = HK;
      }
    }
    //\John G 121212
    //We now have lattice vectors a and b. The c vector can be any vector that is a linear combination
    //of the original primitive cell lattice vectors. Ideally the vector will be short and orthogonal
    //to the plane.
    Eigen::Vector3d normal;
    Eigen::Vector3i L_combination;
    int factor;


    normal = inv_lat_column_mat() * surface_cell.col(0).cross(surface_cell.col(1)); //101112
    factor = 1;

    //Divide by largest value in normal vector. We'll do something similar to when finding the miller
    //indeces. We'll approximate the real normal vector with an integer normal vector, increasing
    //the "resolution" of the integer vector in steps, until the desired accuary is reached.
    if(fabs(normal[0]) >= fabs(normal[1]) && fabs(normal[0]) >= fabs(normal[2]) && normal[0] != 0) {
      normal = normal / normal[0];
    }


    else if(fabs(normal[1]) >= fabs(normal[2]) && fabs(normal[1]) >= fabs(normal[0]) && fabs(normal[1]) != 0) {
      normal = normal / normal[1];
    }

    else {
      normal = normal / normal[2];
    }

    std::cout << "New cell vectors a and b have been generated:" << std::endl;
    std::cout << "Vector A: <" << surface_cell.col(0) << ">" << std::endl;
    std::cout << "Vector B: <" << surface_cell.col(1) << ">" << std::endl;
    std::cout << "Gamma :" << (180 / M_PI)*CASM::angle(surface_cell.col(0), surface_cell.col(1)) << "\u00B0" << std::endl << std::endl;
    std::cout << "Ready to make C..." << std::endl;

    Eigen::Vector3d tnormal;
    //orthoscore represents how close the linear combination is to the plane normal. 1 is perfect, 0 is stupid.
    double orthoscore = 1;
    double torthoscore = 0;

    tnormal = normal;
    double new_vol;

    do {
      normal = tnormal * factor;

      //Get linear combinations by rounding the normal to integers
      L_combination[0] = int(floor(normal[0] + 0.5));
      L_combination[1] = int(floor(normal[1] + 0.5));
      L_combination[2] = int(floor(normal[2] + 0.5));

      //After getting the normal vector in terms of integers, make a linear combination of the initial
      //lattice vectors to get third new lattice vector
      surface_cell.col(2) = lat_column_mat() * L_combination.cast<double>();
      orthoscore = fabs(cos(CASM::angle(lat_column_mat() * normal, surface_cell.col(2))));
      //Only use new linear combination if it's more orthogonal than the previous one
      if(orthoscore > torthoscore + TOL) {
        torthoscore = orthoscore;
        std::cout << "Attempt No." << factor << " to get third lattice vector:" << std::endl;
        std::cout << "Combine: " << L_combination[0] << "*a+" << L_combination[1] << "*b+" << L_combination[2] << "*c" << std::endl << std::endl;
        std::cout << "Cell overview:" << std::endl;
        std::cout << "Orthogonality score: " << orthoscore << std::endl;
        std::cout << "Vector A: <" << surface_cell.col(0) << " >" << std::endl;
        std::cout << "Vector B: <" << surface_cell.col(1) << " >" << std::endl;
        std::cout << "Vector C: <" << surface_cell.col(2) << " >" << std::endl << std::endl;
        std::cout << "Alpha :" << (180 / M_PI)*CASM::angle(surface_cell.col(1), surface_cell.col(2)) << "\u00B0" << std::endl << std::endl;
        std::cout << "Beta :" << (180 / M_PI)*CASM::angle(surface_cell.col(2), surface_cell.col(0)) << "\u00B0" << std::endl << std::endl;
        std::cout << "Gamma :" << (180 / M_PI)*CASM::angle(surface_cell.col(0), surface_cell.col(1)) << "\u00B0" << std::endl << std::endl << std::endl;

        last_surface_cell = surface_cell; //Remember currect generated cell, in case we can't find anything better later

        new_vol = fabs(surface_cell.col(2).dot(surface_cell.col(0).cross(surface_cell.col(1))));

        std::cout << "Volume: " << new_vol << std::endl;
        std::cout << "Volume equivalent to " << new_vol / vol() << " primitive volumes" << std::endl << std::endl;
      }

      factor++;

      if(factor == 100) {
        std::cerr << "Reached an outrageous size. Returning last generated cell" << std::endl;
        surface_cell = last_surface_cell;
        break;
      }

    }
    while(new_vol / vol() < max_vol && orthoscore < 1); //John G 121030

    //std::cout<<"SURFACE:"<<surface_cell<<std::endl;

    Lattice surface_lat(surface_cell.col(0), surface_cell.col(1), surface_cell.col(2));
    surface_lat.make_right_handed();
    SymGroup surf_lat_pg;

    surface_lat.generate_point_group(surf_lat_pg);

    Eigen::Matrix3d transmat;
    surface_lat.is_supercell_of(*this, surf_lat_pg, transmat);

    std::cout << "Your conversion matrix was:" << std::endl << transmat << std::endl;

    return surface_lat;
  }

  //********************************************************************

  ///Are two lattices the same, even if they have different lattice vectors, uses CASM::TOL
  bool Lattice::is_equivalent(const Lattice &B, double tol) const {
    Eigen::Matrix3d T = lat_column_mat().inverse() * B.lat_column_mat();
    return is_unimodular(T, tol) && is_integer(T, tol);
  }

  //********************************************************************

  ///Are lattice vectors identical for two lattices
  bool Lattice:: operator==(const Lattice &RHS) const {
    return almost_equal(RHS.lat_column_mat(), lat_column_mat());
  }


  //John G 011013
  //********************************************************************
  /**
   * Applies all operations of given SymGroup to the lattice and averages
   * out the lattice vectors, changing your lattice to perfectly match
   * with the SymGroup.
   */
  //********************************************************************

  void Lattice::symmetrize(const SymGroup &relaxed_pg) {
    Eigen::Matrix3d tLat2(Eigen::Matrix3d::Zero());
    Eigen::Matrix3d frac_mat;
    for(Index ng = 0; ng < relaxed_pg.size(); ng++) {
      frac_mat = iround(inv_lat_column_mat() * relaxed_pg[ng].matrix() * lat_column_mat()).cast<double>();
      tLat2 += frac_mat.transpose() * lat_column_mat().transpose() * lat_column_mat() * frac_mat;
    }

    tLat2 /= double(relaxed_pg.size());

    // tLat2 has the symmetrized lengths and angles -- it is equal to L.transpose()*L, where L=lat_column_mat()
    // we will find the sqrt of tLat2 and then reorient it so that it matches the original lattice
    Eigen::Matrix3d tMat(tLat2), tMat2;

    Eigen::JacobiSVD<Eigen::Matrix3d> tSVD(tMat);
    tMat = Eigen::Matrix3d::Zero();
    for(int i = 0; i < 3; i++) {
      tMat(i, i) = tSVD.singularValues()[i];
    }

    tMat2 = tSVD.matrixU() * tMat * tSVD.matrixV().transpose();

    tMat = lat_column_mat();

    tSVD.compute(tMat2 * tMat.transpose());

    tMat = tSVD.matrixV() * tSVD.matrixU().transpose() * tMat2;

    for(int i = 0; i < 3; i++) {
      for(int j = 0; j < 3; j++) {
        tLat2(i, j) = tMat(i, j);
      }
    }

    (*this) = Lattice(tLat2);



    return;
  }

  //***********************************************************
  /**
   * Same as the other symmetrize routine, except this one assumes
   * that the symmetry group you mean to use is the point group
   * of your lattice within a certain tolerance.
   */
  //***********************************************************

  void Lattice::symmetrize(double tol) {
    SymGroup point_group;
    generate_point_group(point_group, tol);
    symmetrize(point_group);
    return;
  }

  //***********************************************************

  bool Lattice::is_right_handed() const {
    if(vol() < 0) {
      return false;
    }
    else {
      return true;
    }
  }

  //********************************************************************
  // write Lattice in json as array of vectors
  jsonParser &to_json(const Lattice &lat, jsonParser &json) {
    json.put_array();
    json.push_back(lat[0]);
    json.push_back(lat[1]);
    json.push_back(lat[2]);
    return json;
  };

  //********************************************************************
  // read Lattice from a json array of Eigen::Vector3d
  void from_json(Lattice &lat, const jsonParser &json) {
    try {
      lat = Lattice(json[0].get<Eigen::Vector3d >(), json[1].get<Eigen::Vector3d >(), json[2].get<Eigen::Vector3d >());
    }
    catch(...) {
      /// re-throw exceptions
      throw;
    }
  };

  /// \brief Apply SymOp to a Lattice
  Lattice &apply(const SymOp &op, Lattice &lat) {
    return lat = Lattice(op.matrix() * lat.lat_column_mat());
  }

  /// \brief Copy and apply SymOp to a Lattice
  Lattice copy_apply(const SymOp &op, const Lattice &lat) {
    return Lattice(op.matrix() * lat.lat_column_mat());
  }

  namespace niggli_impl {

    /// Check that x < y, given some tolerance
    ///   Returns x < (y - tol)
    ///   Helper function for niggli
    bool _lt(double x, double y, double tol) {
      return x < (y - tol);
    }

    /// Check that x > y, given some tolerance
    ///   returns y < x - tol
    ///   Helper function for niggli
    bool _gt(double x, double y, double tol) {
      return y < (x - tol);
    }

    /// Check that x == y, given some tolerance
    ///   returns !(_lt(x,y,tol) || _gt(x,y,tol))
    ///   Helper function for niggli
    bool _eq(double x, double y, double tol) {
      return !(_lt(x, y, tol) || _gt(x, y, tol));
    }

    /// Product of off-diagonal signs of S (= lat.transpose()*lat)
    ///   Helper function for niggli
    int _niggli_skew_product_step3(const Eigen::Matrix3d &S, double tol) {
      int S12 = 0;
      if(_gt(S(1, 2), 0.0, tol)) {
        S12 = 1;
      }
      else if(_lt(S(1, 2), 0.0, tol)) {
        S12 = -1;
      }

      int S02 = 0;
      if(_gt(S(0, 2), 0.0, tol)) {
        S02 = 1;
      }
      else if(_lt(S(0, 2), 0.0, tol)) {
        S02 = -1;
      }

      return S12 * S02 * S12;
    }

    /// Product of off-diagonal signs of S (= lat.transpose()*lat)
    ///   Helper function for niggli
    int _niggli_skew_product(const Eigen::Matrix3d &S, double tol) {
      int S12 = 0;
      if(_gt(S(1, 2), 0.0, tol)) {
        S12 = 1;
      }
      else if(_lt(S(1, 2), 0.0, tol)) {
        S12 = -1;
      }

      int S02 = 0;
      if(_gt(S(0, 2), 0.0, tol)) {
        S02 = 1;
      }
      else if(_lt(S(0, 2), 0.0, tol)) {
        S02 = -1;
      }

      int S01 = 0;
      if(_gt(S(0, 1), 0.0, tol)) {
        S01 = 1;
      }
      else if(_lt(S(0, 1), 0.0, tol)) {
        S01 = -1;
      }

      return S12 * S02 * S01;
    }

    /// \brief Returns an equivalent \ref Lattice in Niggli form
    ///
    /// \returns an equivalent \ref Lattice in Niggli form
    ///
    /// \param lat a \ref Lattice
    /// \param tol tolerance for floating point comparisons
    ///
    /// The Niggli cell is a unique choice of lattice vectors for a particular lattice.
    /// It minimizes lattice vector lengths, and chooses a particular angular orientation.
    ///
    /// This implementation function does not set the standard spatial orientation.
    ///
    /// \see
    /// I. Krivy and B. Gruber, Acta Cryst. (1976). A32, 297.
    /// <a href="http://dx.doi.org/10.1107/S0567739476000636">[doi:10.1107/S0567739476000636]</a>
    /// R. W. Grosse-Kunstleve, N. K. Sauter and P. D. Adams, Acta Cryst. (2004). A60, 1.
    /// <a href="http://dx.doi.org/10.1107/S010876730302186X"> [doi:10.1107/S010876730302186X]</a>
    ///
    Lattice _niggli(const Lattice &lat, double tol) {

      //std::cout << "begin _niggli(const Lattice &lat, double tol)" << std::endl;

      // Get helper functions
      using namespace niggli_impl;

      //  S = lat.transpose()*lat, is a matrix of lattice vector scalar products,
      //    this is sometimes called the metric tensor
      //
      //    S(0,0) = a*a, S(0,1) = a*b, etc.
      //
      Eigen::Matrix3d reduced = lat.lat_column_mat();
      Eigen::Matrix3d S = reduced.transpose() * reduced;

      while(true) {

        // in notes: a = reduced.col(0), b = reduced.col(1), c = reduced.col(2)
        //   a*a is scalar product


        // 1)
        // if a*a > b*b or
        //    a*a == b*b and std::abs(b*c*2.0) > std::abs(a*c*2.0),
        // then permute a and b
        if(_gt(S(0, 0), S(1, 1), tol) ||
           (_eq(S(0, 0), S(1, 1), tol) && _gt(std::abs(S(1, 2)), std::abs(S(0, 2)), tol))) {
          reduced.col(0).swap(reduced.col(1));
          reduced *= -1.0;
          S = reduced.transpose() * reduced;
        }


        // 2)
        // if b*b > c*c or
        //    b*b == c*c and std::abs(a*c*2.0) > std::abs(a*b*2.0),
        // then permute b and c
        if(_gt(S(1, 1), S(2, 2), tol) ||
           (_eq(S(1, 1), S(2, 2), tol) && _gt(std::abs(S(0, 2)), std::abs(S(0, 1)), tol))) {
          reduced.col(1).swap(reduced.col(2));
          reduced *= -1.0;
          S = reduced.transpose() * reduced;

          continue;
        }

        // 3)
        // if (2*b*c)*(2*a*c)*(2*a*b) > 0,  *** is this right? **
        if(_niggli_skew_product(S, tol) > 0) {

          // if (2*b*c)*(2*a*c)*(2*b*c) > 0,
          // if(_niggli_skew_product_step3(S, tol) > 0) {

          // if (b*c) < 0.0, flip a
          if(_lt(S(1, 2), 0.0, tol)) {
            reduced.col(0) *= -1;
          }

          // if (a*c) < 0.0, flip b
          if(_lt(S(0, 2), 0.0, tol)) {
            reduced.col(1) *= -1;
          }

          // if (a*b) < 0.0, flip c
          if(_lt(S(0, 1), 0.0, tol)) {
            reduced.col(2) *= -1;
          }

          S = reduced.transpose() * reduced;

        }

        // 4)
        // if (2*b*c)*(2*a*c)*(2*a*b) <= 0,
        if(_niggli_skew_product(S, tol) <= 0) {

          int i = 1, j = 1, k = 1;
          int *p;


          // !! paper says (a*b), but I think they mean (b*c)
          // if (b*c) > 0.0, i = -1
          if(_gt(S(1, 2), 0.0, tol)) {
            i = -1;
          }
          // else if( !(b*c < 0)), p = &i;
          else if(!(_lt(S(1, 2), 0.0, tol))) {
            p = &i;
          }

          // if (a*c) > 0.0, j = -1
          if(_gt(S(0, 2), 0.0, tol)) {
            j = -1;
          }
          // else if( !(a*c < 0)), p = &j;
          else if(!(_lt(S(0, 2), 0.0, tol))) {
            p = &j;
          }

          // if (a*b) > 0.0, k = -1
          if(_gt(S(0, 1), 0.0, tol)) {
            k = -1;
          }
          // else if( !(a*b < 0)), p = &k;
          else if(!(_lt(S(0, 1), 0.0, tol))) {
            p = &k;
          }

          if(i * j * k < 0) {
            *p = -1;
          }

          reduced.col(0) *= i;
          reduced.col(1) *= j;
          reduced.col(2) *= k;

          S = reduced.transpose() * reduced;

        }

        // 5)
        // if std::abs(2.0*S(1,2)) > S(1,1) or
        //    (2.0*S(1,2) == S(1,1) and 2.0*2.0*S(0,2) < 2.0*S(0,1)) or
        //    (2.0*S(1,2) == -S(1,1) and 2.0*S(0,1) < 0.0)
        if(_gt(std::abs(2.0 * S(1, 2)), S(1, 1), tol) ||
           (_eq(2.0 * S(1, 2), S(1, 1), tol) && _lt(2.0 * S(0, 2), S(0, 1), tol)) ||
           (_eq(2.0 * S(1, 2), -S(1, 1), tol) && _lt(S(0, 1), 0.0, tol))) {

          // if (b*c) > 0.0, subtract b from c
          if(_gt(S(1, 2), 0.0, tol)) {
            reduced.col(2) -= reduced.col(1);
          }
          // else if (b*c) < 0.0, add b to c
          else if(_lt(S(1, 2), 0.0, tol)) {
            reduced.col(2) += reduced.col(1);
          }

          S = reduced.transpose() * reduced;

          continue;
        }

        // 6)
        // if std::abs(2.0*S(0,2)) > S(0,0) or
        //    (2.0*S(0,2) == S(0,0) and 2.0*2.0*S(1,2) < 2.0*S(0,1)) or
        //    (2.0*S(0,2) == -S(0,0) and 2.0*S(0,1) < 0.0)
        if(_gt(std::abs(2.0 * S(0, 2)), S(0, 0), tol) ||
           (_eq(2.0 * S(0, 2), S(0, 0), tol) && _lt(2.0 * S(1, 2), S(0, 1), tol)) ||
           (_eq(2.0 * S(0, 2), -S(0, 0), tol) && _lt(S(0, 1), 0.0, tol))) {

          // if (a*c) > 0.0, subtract a from c
          if(_gt(S(0, 2), 0.0, tol)) {
            reduced.col(2) -= reduced.col(0);
          }
          // else if (a*c) < 0.0, add a to c
          else if(_lt(S(0, 2), 0.0, tol)) {
            reduced.col(2) += reduced.col(0);
          }

          S = reduced.transpose() * reduced;

          continue;
        }

        // 7)
        // if std::abs(2.0*S(0,1)) > S(0,0) or
        //    (2.0*S(0,1) == S(0,0) and 2.0*2.0*S(1,2) < 2.0*S(0,2)) or
        //    (2.0*S(0,1) == -S(0,0) and 2.0*S(0,2) < 0.0)
        if(_gt(std::abs(2.0 * S(0, 1)), S(0, 0), tol) ||
           (_eq(2.0 * S(0, 1), S(0, 0), tol) && _lt(2.0 * S(1, 2), S(0, 2), tol)) ||
           (_eq(2.0 * S(0, 1), -S(0, 0), tol) && _lt(S(0, 2), 0.0, tol))) {

          // if (a*b) > 0.0, subtract a from b
          if(_gt(S(0, 1), 0.0, tol)) {
            reduced.col(1) -= reduced.col(0);
          }
          // else if (a*b) < 0.0, add a to b
          else if(_lt(S(0, 1), 0.0, tol)) {
            reduced.col(1) += reduced.col(0);
          }

          S = reduced.transpose() * reduced;

          continue;
        }

        // 8)
        // let tmp = 2*b*c + 2*a*c + 2*a*b + a*a + b*b
        // if  tmp < 0.0 or
        //     tmp == 0 and 2*(a*a + 2*a*c) + 2*a*b > 0
        double tmp = 2.0 * S(1, 2) + 2.0 * S(0, 2) + 2.0 * S(0, 1) + S(0, 0) + S(1, 1);
        if(_lt(tmp, 0.0, tol) ||
           (_eq(tmp, 0.0, tol) && _gt(2.0 * (S(0, 0) + 2.0 * S(0, 2)) + 2.0 * S(0, 1), 0.0, tol))) {

          // add a and b to c
          reduced.col(2) += reduced.col(0);
          reduced.col(2) += reduced.col(1);

          S = reduced.transpose() * reduced;

          continue;
        }

        break;

      } // end while

      return Lattice(reduced);
    }

  }

  /// \brief Returns an equivalent \ref Lattice in Niggli form
  ///
  /// \returns an equivalent \ref Lattice in Niggli form
  ///
  /// \param lat a \ref Lattice
  /// \param tol tolerance for floating point comparisons
  ///
  /// The Niggli cell is a unique choice of lattice vectors for a particular lattice.
  /// It minimizes lattice vector lengths, and chooses a particular angular orientation.
  ///
  /// With the angular orientation fixed, the final spatial orientation of the lattice is
  /// set using standard_orientation.
  ///
  /// \see
  /// I. Krivy and B. Gruber, Acta Cryst. (1976). A32, 297.
  /// <a href="http://dx.doi.org/10.1107/S0567739476000636">[doi:10.1107/S0567739476000636]</a>
  /// R. W. Grosse-Kunstleve, N. K. Sauter and P. D. Adams, Acta Cryst. (2004). A60, 1.
  /// <a href="http://dx.doi.org/10.1107/S010876730302186X"> [doi:10.1107/S010876730302186X]</a>
  ///
  Lattice niggli(const Lattice &lat, const SymGroup &point_grp, double tol) {
    Lattice reduced = niggli_impl::_niggli(lat, tol);
    return standard_orientation(reduced, point_grp, tol);
  }

  /// \brief Rotate the Lattice to a standard orientation using point group operations
  Lattice standard_orientation(const Lattice &lat, const SymGroup &point_grp, double tol) {

    Eigen::Matrix3d start = lat.lat_column_mat();
    Eigen::Matrix3d best = start;
    bool is_sym1 = almost_equal(best, best.transpose(), tol);
    Eigen::Matrix3d tmp;
    Eigen::Matrix3d ptmp = best;
    ptmp.col(0).swap(ptmp.col(2));
    ptmp.row(0).swap(ptmp.row(2));
    bool is_sym2 = almost_equal(ptmp, best, tol);
    // rotate cell so that the lattice vectors are mostly aligned
    //    with Cartesian coordinate axes... and the lattice matrix is symmetric (if possible)
    for(int i = 0; i < point_grp.size(); i++) {

      tmp = point_grp[i].matrix() * start;
      if(tmp.determinant() < 0.0) {
        continue;
      }

      bool better = false;

      bool tmp_sym1 = (almost_equal(tmp, tmp.transpose(), tol));
      if(is_sym1) {
        if(!tmp_sym1) {
          continue;
        }

        ptmp = tmp;
        ptmp.col(0).swap(ptmp.col(2));
        ptmp.row(0).swap(ptmp.row(2));
        bool tmp_sym2 = almost_equal(ptmp, tmp, tol);
        if(is_sym2) {
          if(!tmp_sym2) {
            continue;
          }
        }
        else if(tmp_sym2) {
          is_sym1 = tmp_sym1;
          is_sym2 = tmp_sym2;
          best = tmp;
          continue;
        }
      }
      else if(tmp_sym1) {
        ptmp = tmp;
        ptmp.col(0).swap(ptmp.col(2));
        ptmp.row(0).swap(ptmp.row(2));
        is_sym2 = almost_equal(ptmp, tmp, tol);
        is_sym1 = tmp_sym1;
        best = tmp;
        continue;
      }

      if(almost_equal(best(0, 0), tmp(0, 0), tol)) {
        if(almost_equal(best(1, 0), tmp(1, 0), tol)) {
          if(almost_equal(best(2, 0), tmp(2, 0), tol)) {
            if(almost_equal(best(1, 1), tmp(1, 1), tol)) {
              if(almost_equal(best(0, 1), tmp(0, 1), tol)) {
                if(almost_equal(best(2, 1), tmp(2, 1), tol)) {
                  if(almost_equal(best(2, 2), tmp(2, 2), tol)) {
                    if(almost_equal(best(0, 2), tmp(0, 2), tol)) {
                      if(almost_equal(best(1, 2), tmp(1, 2), tol)) {
                        better = false;
                      }
                      else if(tmp(1, 2) > best(1, 2)) {
                        better = true;
                      }
                    }
                    else if(tmp(0, 2) > best(0, 2)) {
                      better = true;
                    }
                  }
                  else if(tmp(2, 2) > best(2, 2)) {
                    better = true;
                  }
                }
                else if(tmp(2, 1) > best(2, 1)) {
                  better = true;
                }
              }
              else if(tmp(0, 1) > best(0, 1)) {
                better = true;
              }
            }
            else if(tmp(1, 1) > best(1, 1)) {
              better = true;
            }
          }
          else if(tmp(2, 0) > best(2, 0)) {
            better = true;
          }
        }
        else if(tmp(1, 0) > best(1, 0)) {
          better = true;
        }
      }
      else if(tmp(0, 0) > best(0, 0)) {
        better = true;
      }

      if(better) {
        best = tmp;
      }

    }

    return Lattice(best);
  }


  ///\brief returns Lattice that is smallest possible supercell of both input Lattice
  ///
  //*******************************************************************************************
  //
  //  Finds "superduper" Lattice L_{sd} (represented as a matrix with lattice vectors as its columns
  //  such that L_{sd} satisfies
  //        L_{sd} = L_1 * N_1 = L_2 * N_2,     (*1*)
  //  where N_1 and N_2 are integer matrices such that Eq.(*1*) is satisfied and det(N_1) and det(N_2) are minimized.
  //
  //  It is assumed that L_1 = L * M_1 and L_2 = L * M_2  (i.e., L_1 and L_2 are supercells of PRIM lattice L having
  //  integer transformation matrices M_1 and M_2, respectively).
  //
  //  Algorithm proceeds by noting inv(L_2)*L_1 = N_2*inv(N_1) = inv(M_2)*inv(M_1) = A/n, where A is an integer matrix and n = det(M_2). Assuming that
  //  'n' is small (n<10000), we can attempt to find 'n' and 'A'.
  //
  //  Solution: N_2 = A*N_1/n, s.t. det(N_1) is minimized and N_2 is integer
  //
  //  To minimize det(N_1), find smith normal form A = U*S*V, where U and V have det(1), S is diagonal,
  //  and all entries are integer.  Then choose N_1 = inv(V)*R, where R is a diagonal integer matrix with entries
  //             R(i,i)=n/gcf(S(i,i),n)
  //  The resulting solution will have det(M_1*N_1)>=lcm(det(M_1),det(M_2))
  //
  //*******************************************************************************************
  Lattice superdupercell(const Lattice &lat1, const Lattice &lat2) {

    Eigen::Matrix3d dA(lat2.inv_lat_column_mat()*lat1.lat_column_mat());
    long N = 1, num, denom;
    for(Index i = 0; i < 3; i++) {
      for(Index j = 0; j < 3; j++) {
        nearest_rational_number(dA(i, j), num, denom);
        dA *= double(denom);
        N *= denom;
      }
    }

    //std::cout << "dA is:\n" << dA << "\n\n";
    //std::cout << "and N is:\n" << N << "\n\n";
    Eigen::Matrix3l U, S, V;
    smith_normal_form(lround(dA), U, S, V);
    //std::cout << "Smith U is:\n" << U << "\n\n";
    //std::cout << "Smith S is:\n" << S << "\n\n";
    //std::cout << "Smith V is:\n" << V << "\n\n";
    //std::cout << "and U*S*V is:\n" << U*S*V << "\n\n";
    denom = N;

    //reuse matrix S for matrix 'R', as above
    for(Index i = 0; i < 3; i++) {
      S(i, i) = N / gcf(S(i, i), N);
    }
    //Matrix 'N_1', as above is now equal to inverse(V) * S

    Lattice tlat(lat1.lat_column_mat() * (inverse(V)*S).cast<double>());
    return tlat.get_reduced_cell();
  }

  //*******************************************************************************************

  /// Check if scel is a supercell of unitcell unit and some integer transformation matrix T
  std::pair<bool, Eigen::MatrixXi> is_supercell(const Lattice &scel, const Lattice &unit, double tol) {
    // check scel = unit*T, with integer T
    Eigen::MatrixXd T = unit.inv_lat_column_mat() * scel.lat_column_mat();

    if(is_integer(T, tol) && !almost_zero(T)) {
      return std::make_pair(true, iround(T));
    }
    return std::make_pair(false, T.cast<int>());
  }


}



