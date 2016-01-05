#include "casm/crystallography/Lattice.hh"

#include "casm/crystallography/SupercellEnumerator.hh"

namespace CASM {



  //Lattice::Lattice() {}


  //********************************************************************

  Lattice::Lattice(const Vector3<double> &vec1, const Vector3<double> &vec2,
                   const Vector3<double> &vec3) : vecs(3) {
    vecs[0] = vec1;
    vecs[1] = vec2;
    vecs[2] = vec3;
    calc_conversions();
    calc_properties();
  }

  //********************************************************************

  ///Construct Lattice from a matrix of lattice vectors, where lattice vectors are columns
  ///(e.g., lat_mat is equivalent to coord_trans_mat[FRAC])
  Lattice::Lattice(const Matrix3<double> &lat_mat) : vecs(3) {
    for(int i = 0; i < 3; i++) {
      for(int j = 0; j < 3; j++)
        vecs[i][j] = lat_mat(j, i);
    }
    calc_conversions();
    calc_properties();
  }

  //********************************************************************

  ///Construct Lattice from a matrix of lattice vectors, where lattice vectors are columns
  ///(e.g., lat_mat is equivalent to coord_trans_mat[FRAC])
  Lattice::Lattice(const Eigen::Matrix3d &lat_mat) : vecs(3) {
    for(int i = 0; i < 3; i++) {
      for(int j = 0; j < 3; j++)
        vecs[i][j] = lat_mat(j, i);
    }
    calc_conversions();
    calc_properties();
  }

  //********************************************************************

  Lattice::Lattice(const Lattice &RHS) : vecs(RHS.vecs) {
    calc_conversions();
    calc_properties();
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

  Lattice &Lattice::operator=(const Lattice &RHS) {
    if(this == &RHS)
      return *this;
    vecs = RHS.vecs;

    voronoi_table.clear();

    calc_conversions();
    calc_properties();

    return *this;
  }

  //********************************************************************
  /**
   * Calculates conversion matrices.
   *
   * Calculates matrices for both F (fractional) to
   * C (cartesian) and C to F conversions.
   */
  //********************************************************************

  void Lattice::calc_conversions() {
    for(int i = 0; i < 3; i++) {
      for(int j = 0; j < 3; j++) {
        coord_trans_mat[FRAC](j, i) = vecs[i][j];
      }
    }
    coord_trans_mat[CART] = coord_trans_mat[FRAC].inverse();

    return;
  }

  //********************************************************************
  /**
   * Calculates length and angles.
   *
   * Calculates the lengths of each unit cell vector
   * and assigns it to the lengths array.  Also calculates
   * angles between the unit cell vectors and assigns it
   * to the angles array.
   */
  void Lattice::calc_properties() {
    // Calculates Lengths
    for(int i = 0; i < 3; i++)
      lengths[i] = vecs[i].length();

    // Calculates Angles
    for(int i = 0; i < 3; i++) {

      angles[i] = vecs[(i + 1) % 3].dot(vecs[(i + 2) % 3]) / (lengths[(i + 1) % 3] * lengths[(i + 2) % 3 ]);

      //Make sure that cos(angle) is between 0 and 1
      if((angles[i] - 1.0) > 0.0)
        angles[i] = 1.0;

      if((angles[i] + 1.0) < 0.0)
        angles[i] = -1.0;

      angles[i] = (180.0 / M_PI) * acos(angles[i]);
    }

  }

  //********************************************************************

  Lattice Lattice::scaled_lattice(double scale) {
    return Lattice(scale * coord_trans(FRAC));
  }

  //********************************************************************

  void Lattice::read(std::istream &stream) {
    double scale;
    stream >> scale;
    vecs.resize(3);
    for(int i = 0; i < 3; i++) {
      stream >> vecs[i];
      vecs[i] *= scale;
    }

    calc_conversions();
    calc_properties();
    return;
  }

  //********************************************************************

  void Lattice::print(std::ostream &stream, int _prec) const {
    int tprec = stream.precision();
    std::ios::fmtflags tflags = stream.flags();
    stream.precision(_prec);
    stream.width(_prec+3);
    stream.flags(std::ios::showpoint | std::ios::fixed | std::ios::right);
    stream  << 1.0 << '\n';

    stream << ' ' << std::setw(_prec+8) << vecs[0] << '\n';
    stream << ' ' << std::setw(_prec+8) << vecs[1] << '\n';
    stream << ' ' << std::setw(_prec+8) << vecs[2] << '\n';

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
    return Lattice(2 * M_PI * coord_trans_mat[CART].transpose()); //equivalent expression
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
      prim_vec_lengths.push_back(prim_recip_lat.lengths[i]);
    }
    //std::cout<<prim_vec_lengths<<"\n";
    double shortest = prim_vec_lengths.min();
    int short_ind = prim_vec_lengths.find(shortest);

    //std::cout << "prim kpoints "<< prim_kpoints << std::endl;
    //    std::cout << "prim recip vol " << prim_recip_lat.vol << std::endl;
    //std::cout << "prim kpoint density " << prim_density << std::endl;

    double scale = (prim_kpoints[short_ind] / shortest);

    for(int i = 0; i < 3; i++) {
      super_kpoints[i] = int(ceil(scale * recip_lat.lengths[i]));
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
        super_kpoints[i] = int(ceil(scale * recip_lat.lengths[i]));
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
    Matrix3<double> tfrac_op, tMat;
    for(Index ng = 0; ng < super_group.size(); ng++) {
      tfrac_op = lat_column_mat().inverse() * super_group[ng].get_matrix(CART) * lat_column_mat();

      //Use a soft tolerance of 1% to see if further screening should be performed
      if(!almost_equal(1.0, std::abs(tfrac_op.determinant()), 0.01) || !tfrac_op.is_integer(0.01))
        continue;

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
      tMat = coord_trans_mat[CART].transpose() * (tfrac_op.transpose() * coord_trans_mat[FRAC].transpose() * coord_trans_mat[FRAC] * tfrac_op);

      // Subtract lat_column_mat() from tMat, leaving us with (F.transpose()*F - Identity)*lat_column_mat().
      // This is 2*E*lat_column_mat(), where E is the green-lagrange strain
      tMat = (tMat - coord_trans_mat[FRAC]) / 2.0;

      //... and then multiplying by the transpose...
      tMat = tMat * tMat.transpose();

      // The diagonal elements of tMat describe the square of the distance by which the transformed vectors 'miss' the original vectors
      if(tMat(0, 0) < pg_tol * pg_tol && tMat(1, 1) < pg_tol * pg_tol && tMat(2, 2) < pg_tol * pg_tol)
        sub_group.push_back(super_group[ng]);

    }

  }
  //********************************************************************

  void Lattice::generate_point_group(SymGroup &point_group, double pg_tol) const {

    if(point_group.size() != 0) {
      std::cerr << "WARNING in Lattice::generate_point_group" << std::endl;
      std::cerr << "The point group for your lattice isn't empty and it's about to be rewritten!" << std::endl;
      point_group.clear();
    }

    //Enumerate all possible matrices with elements equal to -1, 0, or 1
    //These represent operations that reorder lattice vectors or replace one
    //or more lattice vectors with a face or body diagonal.
    Matrix3<double> tMat, tOp_cart;
    Counter<Matrix3<int> > pg_count(Matrix3<int>(-1),
                                    Matrix3<int>(1),
                                    Matrix3<int>(1));

    //For this algorithm to work, lattice needs to be in reduced form.
    Lattice tlat_reduced(get_reduced_cell());
    do {

      //continue if determinant is not 1, because it doesn't preserve volume
      if(std::abs(pg_count().determinant()) != 1) continue;

      tOp_cart = (tlat_reduced.coord_trans_mat[FRAC] * pg_count()) * tlat_reduced.coord_trans_mat[CART];

      //Find the effect of applying symmetry to the lattice vectors
      //The following is equivalent to point_group[i].get_matrix(CART).transpose()*tlat_reduced.coord_trans_mat[FRAC]*point_group[i].get_matrix(FRAC)
      tMat = tOp_cart.transpose() * (tlat_reduced.coord_trans_mat[FRAC] * pg_count());

      //If pg_count() is a point_group operation, tMat should be equal to tlat_reduced.coord_trans_mat[FRAC].  We check by first taking the difference...
      tMat = (tMat - tlat_reduced.coord_trans_mat[FRAC]) / 2.0;

      //... and then multiplying by the transpose...
      tMat = tMat * tMat.transpose();

      // The diagonal elements are square of the distances by which the transformed lattice vectors "miss" the original lattice vectors
      // If they are less than the square of the tolerance, we add the operation to the point group
      if(tMat(0, 0) < pg_tol * pg_tol && tMat(1, 1) < pg_tol * pg_tol && tMat(2, 2) < pg_tol * pg_tol) {

        Array<double> diags;
        diags.push_back(tMat(0, 0));
        diags.push_back(tMat(1, 1));
        diags.push_back(tMat(2, 2));

        point_group.push_back(SymOp(tOp_cart, *this, CART, sqrt(diags.max())));
      }

      /*** Old way of checking pg_count()

      if(tOp_cart.is_unitary(pg_tol)) {
      point_group.push_back(SymOp(tOp_cart, *this, CART));
      point_group.back().get_sym_type();
      }
      **/
    }
    while(++pg_count);

    if(!point_group.is_group(pg_tol)) {
      std::cerr << "*** WARNING *** \n"
                << "    Lattice::generate_point_group() has been called on an ill-conditioned lattice \n"
                << "    (i.e., a well-defined point group could not be found with the supplied tolerance of " << pg_tol << ").\n"
                << "    CASM will use the group closure of the symmetry operations that were found.  Please consider using the \n"
                << "    CASM symmetrization tool on your input files.\n";
      point_group.enforce_group(pg_tol);

    }
    //Sort point_group by trace/conjugacy class
    point_group.sort_by_class();

    for(Index i = 0; i < point_group.size(); i++) {
      point_group[i].get_sym_type();
    }

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
      tarray.push_back(point_group[i].get_map_error());
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
      point_group.get_character_table();
      old_name.push_back(point_group.get_name());
      num_ops.push_back(point_group.size());
      is_group.push_back(point_group.is_group(i));
      point_group.enforce_group(i);
      is_group_now.push_back(point_group.is_group(i));
      num_enforced_ops.push_back(point_group.size());
      point_group.get_character_table();
      new_name.push_back(point_group.get_name());
    }

    for(Index i = 0; i < tols.size(); i++) {
      std::cout << tols[i] << "\t" << num_ops[i] << "\t" << is_group[i] << "\t" << old_name[i] << "\t" << num_enforced_ops[i] << "\t" << is_group_now[i] << "\t" << new_name[i] << "\n";
    }

    return;
  }



  //********************************************************************
  // The mathematical description
  // A : lattice vectors [ a00 a01 a02  a10 a11 a12 a20a21a22]
  // M : transformation matrix
  // A': new lattice vectors (transformed lattice vectors, A)
  // A' = M*A
  // If the matrix elements M[i][j] are intergers and ||M|| =1, then the lattices A and A' coincide.
  // ||M|| > 1, then the lattice A' is superlattice of the lattice A', and the volume of the primitive
  // cell in A' is ||M|| times greater than the volume of the primitive cell in A.

  // Algorithm
  // 1. Make lattice vectors as a upper triangular matrix where the product of the diagonal
  //    elements equals the volume. for the elements above the diagonal choose all values less than
  //	the diagnoal element.
  //
  // 2. To check the lattice vectors found by looping is linearly independece from
  //	other previous lattice vectors A,
  //		2.1 apply point sysmetry to the lattice vectors Ai ==> point_sy*Ai = Ai'
  //		( M*B = Ai' ==> B is the lattice vectors previously found (stored in tsupercell))
  // 		2.2 M = Ai' * B.inverse(); and check the matrix elements N[i][j] are non-intergers and
  //			||M|| is great than 1.
  //4. Add it to the supercell list.
  /*
    void Lattice::generate_supercells(Array<Lattice> &supercell, const SymGroup &effective_pg, int max_prim_vol, int min_prim_vol) const {
      int vol;
      Index pg, ts;
      Matrix3<int> tslat(0);
      Matrix3<double> lin_com, tsup_lat_mat, tsym_lat_mat;
      Array<Lattice> tsupercell;
      supercell.clear();

      for(vol = min_prim_vol; vol <= max_prim_vol; vol++) {
        for(tslat(0, 0) = 1; tslat(0, 0) <= vol; tslat(0, 0)++) {
          if(vol % tslat(0, 0) != 0) continue; //Changed by John
          for(tslat(1, 1) = 1; tslat(1, 1) <= vol / tslat(0, 0); tslat(1, 1)++) {
            if((vol / tslat(0, 0)) % tslat(1, 1) != 0) continue; //Changed by John
            tslat(2, 2) = vol / (tslat(0, 0) * tslat(1, 1));

            for(tslat(0, 1) = 0; tslat(0, 1) < tslat(0, 0); tslat(0, 1)++) {
              for(tslat(0, 2) = 0; tslat(0, 2) < tslat(0, 0); tslat(0, 2)++) {
                for(tslat(1, 2) = 0; tslat(1, 2) < tslat(1, 1); tslat(1, 2)++) {

                  tsup_lat_mat = coord_trans_mat[FRAC] * tslat;

                  bool is_unique = true; //added by John
                  for(ts = 0; ts < tsupercell.size(); ts++) {
                    for(pg = 0; pg < effective_pg.size(); pg++) {

                      //tsym_lat_mat is supercell lattice vectors transformed by point group symmetry operation
                      tsym_lat_mat = effective_pg[pg].get_matrix(CART) * tsup_lat_mat;

                      //lin_com is matrix containing fractional coordinates of transformed candidate lattice vectors, in terms of supercell ts
                      lin_com = tsupercell[ts].coord_trans_mat[CART] * tsym_lat_mat;

                      if(lin_com.is_integer()) {
                        is_unique = false;
                        break;
                      }

                    }
                    if(!is_unique) break;  //Break out of the outside loop if is_unique is already false;
                  }

                  //We add the supercell to the list if, after all the checks, it is unique
                  if(is_unique) {

                    Lattice n = niggli(Lattice(tsup_lat_mat), TOL);
                    Lattice r = Lattice(tsup_lat_mat).get_reduced_cell();

                    if(!(n == n)) {
                      std::cout << "Niggli cell:\n";
                      n.print(std::cout);
                      std::cout << "\n\nReduced:\n";
                      r.print(std::cout);
                      std::cout << "\n\n";
                    }

                    // std::cout<<"----------\n"<<tslat<<"\n----------\n";
                    tsupercell.push_back(niggli(Lattice(tsup_lat_mat), TOL));  //Lattice constructor takes matrix where rows are lattice vectors

                  }


                }
              }
            }  //end loops over off-diagonals

          }  //end loop over second diagonal
        }    //end loop over third diagonal

        supercell.append(tsupercell);
        tsupercell.clear();
      }//end loop over allowed volumes

      return;
    }
  */

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
    SupercellEnumerator<Lattice> enumerator(*this, effective_pg, min_prim_vol, max_prim_vol + 1);
    supercell.clear();
    for(auto it = enumerator.begin(); it != enumerator.end(); ++it) {
      supercell.push_back(niggli(*it, effective_pg, TOL));
    }
    return;
  }


  //********************************************************************
  /*
  void Lattice::generate_supercells(Array<Lattice> &supercell, const MasterSymGroup &factor_group, int max_prim_vol, int min_prim_vol) const {
    generate_supercells(supercell, factor_group.point_group(), max_prim_vol, min_prim_vol);
    return;
  }
  */
  //********************************************************************
  // Find the smallest strain that takes (*this) lattice to a lattice that is symmetrically
  // equivalent to strained_lattice. "Equivalent" means that the resulting lattice may be rotated
  // by an arbitrary axis-angle relative to strained_lattice, and is determined by the crystallographic
  // setting of (*this) lattice.
  /*
  Eigen::Matrix3d Lattice::nearest_equivalent_strain(const Lattice &strained_lattice) const {
    Lattice thislat(get_reduced_cell()), strained_lat(strained_lattice.get_reduced_cell());

    Counter<Matrix3<int> > imat_count(Matrix3<int>(-2),
                                      Matrix3<int>(2),
                                      Matrix3<int>(1));

    Matrix3<double> FTF, L2TL2, L1Tinv, L1inv, FTFbest, identity(0);
    Matrix3<int> Nbest;
    identity(0, 0) = identity(1, 1) = identity(2, 2) = 1;
    double min_norm(1e20), tnorm;
    L1inv = thislat.lat_column_mat().inverse();
    L1Tinv = L1inv.transpose();
    std::cout << "Volume of lat1 is " << thislat.lat_column_mat().determinant() << "\n";
    std::cout << "Volume of lat2 is " << strained_lat.lat_column_mat().determinant() << "\n";
    L2TL2 = strained_lat.lat_column_mat().transpose() * strained_lat.lat_column_mat();
    //For this algorithm to work, lattice needs to be in reduced form.
    int n = 0;
    do {
      n++;
      //continue if determinant is not 1, because it doesn't preserve volume
      if(std::abs(imat_count().determinant()) != 1) continue;

      FTF = (L1Tinv * imat_count().transpose()) * (L2TL2 * imat_count()) * L1inv;

      tnorm = (FTF - identity).norm();
      if(tnorm < min_norm) {
        min_norm = tnorm;
        FTFbest = FTF;
        std::cout << "Found new FTFbest with n = " << n << " and min_norm = " << tnorm << '\n' << FTFbest << "\n\n";
      }
    }
    while(++imat_count);
    std::cout << "Nbest was \n" << Nbest << "\n\n";
    Eigen::Matrix3d tmat;
    for(int i = 0; i < 3; i++) {
      for(int j = 0; j < 3; j++) {
        tmat(i, j) = FTFbest(i, j);
      }
    }
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> sqrtsolver(tmat);
    return sqrtsolver.operatorSqrt() - Eigen::Matrix3d::Identity();
  }
  */
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
    Array<Matrix3<double> > skew;
    Matrix3<double> tskew(Matrix3<double>::identity());
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


    Matrix3<double> tMat, reduced_lat_mat;
    Matrix3<double> tMat_mags, reduced_mags;
    reduced_lat_mat = coord_trans_mat[FRAC];  //Matrix3
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

  Vector3<double> Lattice::max_voronoi_vector(const Vector3<double> &pos)const {

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
  int Lattice::voronoi_number(const Vector3<double> &pos)const {

    int tnum = 0;
    double tproj = 0;

    if(!voronoi_table.size()) {
      generate_voronoi_table();
    }

    for(Index nv = 0; nv < voronoi_table.size(); nv++) {
      tproj = pos.dot(voronoi_table[nv]);
      if(almost_equal(tproj, 1.0) < TOL)
        tnum++;
      else if(tproj > 1.0)
        return -1;
    }


    return tnum;
  }
  //********************************************************************

  void Lattice::generate_voronoi_table() const {
    voronoi_table.clear();
    //There are no fewer than 12 points in the voronoi table
    voronoi_table.reserve(12);

    Vector3<double> tpoint;
    int i;

    Lattice tlat_reduced(get_reduced_cell());
    //Count over all lattice vectors, face diagonals, and body diagonals
    //originating from origin;
    Counter<Vector3<int> > combo_count(Vector3<int>(-1, -1, -1),
                                       Vector3<int>(1, 1, 1),
                                       Vector3<int>(1, 1, 1));

    //std::cout << "For angles " << angles[0] << ", " << angles[1] << ", " << angles[2] << ", Voronoi table is: \n";
    //For each linear combination, check to see if it is on a face, edge, or vertex of the voronoi cell
    do {


      if(!(combo_count[0] || combo_count[1] || combo_count[2])) continue;

      //A linear combination does not fall on the voronoi boundary if the angle between
      //any two of the vectors forming that combination are obtuse
      for(i = 0; i < 3; i++) {
        if((combo_count[(i + 1) % 3] && combo_count[(i + 2) % 3])
           && std::abs(90.0 * abs(combo_count[(i + 1) % 3] - combo_count[(i + 2) % 3]) - angles[i]) + TOL < 90)
          break;
      }

      if(i == 3) {
        tpoint = tlat_reduced.coord_trans_mat[FRAC] * combo_count();
        double tval(tpoint.length());

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
  Vector3<int> Lattice::enclose_sphere(double radius) const {

    Vector3<double> normals[2];
    Vector3<int> dimension;

    for(int i = 0; i < 3; i++) {
      normals[CART] = vecs[(i + 1) % 3].cross(vecs[(i + 2) % 3]);
      normals[CART] = normals[CART] * radius / normals[CART].length();
      normals[FRAC] = coord_trans_mat[CART] * normals[CART];
      dimension[i] = (int)ceil(std::abs(normals[FRAC][i]));
    }
    return dimension;
  }

  //********************************************************************
  //Change bool to an array of SymOps you want to use, default point group
  //Overload to only use identity
  //Return N matrix
  bool Lattice::is_supercell_of(const Lattice &tile, const Array<SymOp> &symoplist, Matrix3<double> &multimat, double _tol) const {
    if(symoplist.size() == 0) { //John G 121212 extra error
      std::cerr << "ERROR: In Lattice::is_supercell_of. You've passed a point group with no elements!" << std::endl;
      std::cerr << "       You need at least identity in your group. Exiting..." << std::endl;
      exit(1);
    }

    Matrix3<double> tsym_lat_mat;
    for(Index pg = 0; pg < symoplist.size(); pg++) {

      tsym_lat_mat = symoplist[pg].get_matrix(CART) * coord_trans_mat[FRAC];
      multimat = tile.coord_trans_mat[CART] * tsym_lat_mat;

      if(multimat.is_integer(_tol) && !multimat.is_zero(_tol))
        return true;

    }
    return false;
  }

  //********************************************************************
  //Change bool to an array of SymOps you want to use, default point group
  //Overload to only use identity
  //Return N matrix
  bool Lattice::is_supercell_of(const Lattice &tile, Matrix3<double> &multimat, double _tol) const {

    multimat = tile.coord_trans_mat[CART] * coord_trans_mat[FRAC];

    if(multimat.is_integer(_tol) && !multimat.is_zero(_tol))
      return true;

    return false;
  }

  //********************************************************************
  //Overload to only use identity
  //Return N matrix
  bool Lattice::is_supercell_of(const Lattice &tile, double _tol) const {
    Matrix3<double> multimat;
    return is_supercell_of(tile, multimat, _tol);
  }

  //********************************************************************

  bool Lattice::is_supercell_of(const Lattice &tile, const Array<SymOp> &symoplist, double _tol) const {
    Matrix3<double> multimat;
    return is_supercell_of(tile, symoplist, multimat, _tol);
  }

  //********************************************************************
  /// Finds 'new_scel' equivalent to '*this' and 'new_prim' equivalent to 'prim', such that 'new_prim' perfectly tiles 'new_scel'
  /// Returns true if tessellation cannot be found.
  bool Lattice::find_tessellation(Lattice &prim, Lattice &new_scel, Lattice &new_prim) const {
    SymGroup prim_pg;
    prim.generate_point_group(prim_pg);
    Matrix3<double> dMat;
    Matrix3<int> iMat;
    if(!is_supercell_of(prim, prim_pg, dMat)) return false;
    for(int i = 0; i < 9; i++) {
      iMat[i] = round(dMat[i]);
    }

    Matrix3<int> U, S, V;
    iMat.smith_normal_form(U, S, V);
    if(!S.is_diagonal()) return false;
    new_prim = Lattice(prim.coord_trans_mat[FRAC] * U);
    new_scel = Lattice(new_prim.coord_trans_mat[FRAC] * S);
    return true;
  }

  //********************************************************************

  Lattice &Lattice::make_right_handed() {

    if(lat_column_mat().determinant() < 0) {
      swap(vecs[0], vecs[1]);
      calc_conversions();
      calc_properties();
    }

    return *this;
  }

  //\John G 121212
  //********************************************************************************************************

  Vector3< int > Lattice::get_millers(Vector3< double > plane_normal, double tolerance) const {
    Lattice recip_lattice = get_reciprocal();
    Coordinate tnorm(plane_normal, recip_lattice, CART);
    Vector3< double > double_millers;

    //Get fractional coordinates of plane_normal in recip_lattice
    //These are h, k, l
    //For miller indeces h, k and l    plane_normal[CART]=h*a.recip+k*b.recip+l*c.recip
    double_millers = tnorm(FRAC);
    return double_millers.scale_to_int(tolerance);
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

  Lattice Lattice::get_lattice_in_plane(Vector3< int > millers, int max_vol) const {  //John G 121030
    //Hold new lattice vectors in these. Then at the end we make an actual Lattice out of it
    Array<Vector3< double > > surface_cell(3);     //Holds new lattice vectors, two of which are in the surface plane
    Array<Vector3< double > > last_surface_cell(3);

    //Miller indeces of 100, 010 or 001 mean you don't need a new cell to expose the plane, however
    //you may want to reorient the vectors so that the ab plane is exposed (useful for Structure::stitch)

    if(millers == Vector3<int>(0, 1, 0)) {
      std::cout << "No chopping neccesary" << std::endl;
      std::cout << "Flipping your vectors to bring a and b into plane:" << std::endl;
      surface_cell[0] = vecs[2];
      surface_cell[1] = vecs[0];
      surface_cell[2] = vecs[1];
      std::cout << "b --> c" << std::endl;
      std::cout << "a --> b" << std::endl;
      std::cout << "c --> a" << std::endl;


      return Lattice(surface_cell[0], surface_cell[1], surface_cell[2]);
    }

    else if(millers == Vector3<int>(1, 0, 0)) {
      std::cout << "No chopping neccesary" << std::endl;
      std::cout << "Flipping your vectors to bring a and b into plane:" << std::endl;
      surface_cell[1] = vecs[2];
      surface_cell[0] = vecs[1];
      surface_cell[2] = vecs[0];
      std::cout << "a --> c" << std::endl;
      std::cout << "b --> a" << std::endl;
      std::cout << "c --> b" << std::endl;

      return Lattice(surface_cell[0], surface_cell[1], surface_cell[2]);
    }

    else if(millers == Vector3<int>(0, 0, 1)) { // || millers==Vector3<int>(0,1,0) || millers==Vector3<int>(1,0,0))
      std::cout << "Silly goose! You don't need a new lattice." << std::endl;
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

      surface_cell[0] = vecs[zero];

      Vector3 < double > H_miller_point, K_miller_point;
      Vector3 < double > HK;
      Vector3 < double > millers_dubs;

      //Turn integer millers into doubles for mathematical purposes (inverse)
      millers_dubs = millers;

      //std::cout<<millers_dubs<<std::endl;

      Vector3< double > inv_miller_dubs;
      Vector3< int > inv_miller;


      //In actualility, the inverse miller of 0 is infinity, we set it to 0 here so that we can use
      //scale_to_int without going crazy. Since we won't be using this inverse miller index it's OK
      inv_miller_dubs[zero] = 0;
      inv_miller_dubs[(zero + 1) % 3] = 1.0 / millers_dubs[(zero + 1) % 3];
      inv_miller_dubs[(zero + 2) % 3] = 1.0 / millers_dubs[(zero + 2) % 3];

      inv_miller = inv_miller_dubs.scale_to_int(TOL);
      H_miller_point = inv_miller[(zero + 1) % 3] * vecs[(zero + 1) % 3];
      K_miller_point = inv_miller[(zero + 2) % 3] * vecs[(zero + 2) % 3];

      std::cout << "inv millers dubs: " << inv_miller_dubs << std::endl;
      std::cout << "inv millers : " << inv_miller << std::endl;

      HK = K_miller_point - H_miller_point;
      surface_cell[1] = HK;
    }



    else {
      //Get three points that lie on the plane
      //We'll want to find points that lie on the plane AND land of lattice points. In order to do
      //this we need the miller inverses multiplied by a factor that makes them integers
      Vector3 < double > H_miller_point, K_miller_point, L_miller_point;
      Vector3< double > inv_miller_dubs;
      Vector3< int > inv_miller;
      //Turn integer millers into doubles for mathematical purposes (inverse)
      Vector3 < double > millers_dubs;
      millers_dubs = millers;

      inv_miller_dubs[0] = 1.0 / millers_dubs[0];
      inv_miller_dubs[1] = 1.0 / millers_dubs[1];
      inv_miller_dubs[2] = 1.0 / millers_dubs[2];

      inv_miller = inv_miller_dubs.scale_to_int(TOL);

      H_miller_point = inv_miller[0] * vecs[0];
      K_miller_point = inv_miller[1] * vecs[1];
      L_miller_point = inv_miller[2] * vecs[2];

      //Get three vectors that connect the three points on the plane to each other. Any two of the following
      //vectors could be used for constructing the new lattice, but it's convenient to pick the two
      //most orthogonal vectors
      Vector3 < double > HK, KL, LH;
      Vector3 < double > tangles;


      HK = K_miller_point - H_miller_point;
      KL = L_miller_point - K_miller_point;
      LH = H_miller_point - L_miller_point;

      //John G 121212
      //The vectors that we got at this point are valid, but sometimes larger than they need to be.
      Vector3<Vector3<double> > templat;
      templat[0] = HK;
      templat[1] = KL;
      templat[2] = LH;

      //Find shortest vector
      int s = 0;
      for(int i = 1; i < 3; i++) {
        if(vecs[i].norm() < vecs[s].norm()) {
          s = i;
        }
      }

      //Try dividing by integers and see if they're still lattice combinations. If they are
      //shorten and continue
      for(int i = 0; i < 3; i++) {
        int maxval = round(templat[i].norm() / vecs[i].norm() + 1);
        for(int j = maxval; j > 1; j--) {
          Vector3<double> shortened = coord_trans_mat[CART] * templat[i] / j; //Shorten and convert to fractional
          bool combo = true;

          for(int k = 0; k < 3; k++) {
            if(!almost_zero(std::abs(round(shortened[k]) - shortened[k]))) {
              combo = false;
              break;
            }
          }

          if(combo) {
            templat[i] = coord_trans_mat[FRAC] * shortened;
            break;
          }
        }
      }

      HK = templat[0];
      KL = templat[1];
      LH = templat[2];

      //We select the two vectors that spawn the smallest area

      double HKKL, KLLH, LHHK;
      HKKL = HK.cross(KL).norm();
      KLLH = KL.cross(LH).norm();
      LHHK = LH.cross(HK).norm();

      if(HKKL <= KLLH && HKKL <= LHHK) {
        surface_cell[0] = HK;
        surface_cell[1] = KL;
      }

      else if(KLLH <= HKKL && KLLH <= LHHK) {
        surface_cell[0] = KL;
        surface_cell[1] = LH;
      }

      else {
        surface_cell[0] = LH;
        surface_cell[1] = HK;
      }
    }
    //\John G 121212
    //We now have lattice vectors a and b. The c vector can be any vector that is a linear combination
    //of the original primitive cell lattice vectors. Ideally the vector will be short and orthogonal
    //to the plane.
    Vector3 < double > normal;
    Vector3 < int > L_combination;
    int factor;


    normal = coord_trans_mat[CART] * surface_cell[0].cross(surface_cell[1]); //101112
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
    std::cout << "Vector A: <" << surface_cell[0] << ">" << std::endl;
    std::cout << "Vector B: <" << surface_cell[1] << ">" << std::endl;
    std::cout << "Gamma :" << (180 / M_PI)*surface_cell[0].get_angle(surface_cell[1]) << "\u00B0" << std::endl << std::endl;
    std::cout << "Ready to make C..." << std::endl;

    Vector3< double > tnormal;
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
      surface_cell[2] = L_combination[0] * vecs[0] + L_combination[1] * vecs[1] + L_combination[2] * vecs[2];
      orthoscore = fabs((cos((coord_trans_mat[FRAC] * normal).get_angle(surface_cell[2]))));
      //Only use new linear combination if it's more orthogonal than the previous one
      if(orthoscore > torthoscore + TOL) {
        torthoscore = orthoscore;
        std::cout << "Attempt No." << factor << " to get third lattice vector:" << std::endl;
        std::cout << "Combine: " << L_combination[0] << "*a+" << L_combination[1] << "*b+" << L_combination[2] << "*c" << std::endl << std::endl;
        std::cout << "Cell overview:" << std::endl;
        std::cout << "Orthogonality score: " << orthoscore << std::endl;
        std::cout << "Vector A: <" << surface_cell[0] << " >" << std::endl;
        std::cout << "Vector B: <" << surface_cell[1] << " >" << std::endl;
        std::cout << "Vector C: <" << surface_cell[2] << " >" << std::endl << std::endl;
        std::cout << "Alpha :" << (180 / M_PI)*surface_cell[1].get_angle(surface_cell[2]) << "\u00B0" << std::endl << std::endl;
        std::cout << "Beta :" << (180 / M_PI)*surface_cell[2].get_angle(surface_cell[0]) << "\u00B0" << std::endl << std::endl;
        std::cout << "Gamma :" << (180 / M_PI)*surface_cell[0].get_angle(surface_cell[1]) << "\u00B0" << std::endl << std::endl << std::endl;

        last_surface_cell = surface_cell; //Remember currect generated cell, in case we can't find anything better later

        new_vol = fabs(surface_cell[2].dot(surface_cell[0].cross(surface_cell[1])));

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

    Lattice surface_lat(surface_cell[0], surface_cell[1], surface_cell[2]);
    surface_lat.make_right_handed();
    SymGroup surf_lat_pg;

    surface_lat.generate_point_group(surf_lat_pg);

    Matrix3<double> transmat;
    surface_lat.is_supercell_of(*this, surf_lat_pg, transmat);

    std::cout << "Your conversion matrix was:" << std::endl << transmat << std::endl;

    return surface_lat;
  }

  //********************************************************************

  ///Are two lattices the same, even if they have different lattice vectors, uses CASM::TOL
  bool Lattice::is_equivalent(const Lattice &B) const {

    Lattice niggli_A = niggli_impl::_niggli(*this, TOL);
    Lattice niggli_B = niggli_impl::_niggli(B, TOL);
    SymGroup point_grp_A;
    niggli_A.generate_point_group(point_grp_A, TOL);

    return std::find_if(point_grp_A.cbegin(),
                        point_grp_A.cend(),
    [&](const SymOp & op) {
      return niggli_B == Lattice(op.get_matrix(CART) * niggli_A.lat_column_mat());
    }) != point_grp_A.cend();
  }

  //********************************************************************

  ///Are lattice vectors identical for two lattices
  bool Lattice:: operator==(const Lattice &RHS) const {
    for(int i = 0; i < 3; i++) {
      if((vecs[i] - RHS[i]).is_zero(TOL))
        continue;
      else
        return false;
    }

    return true;
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
    Matrix3<double> tLat2(0);
    for(Index ng = 0; ng < relaxed_pg.size(); ng++) {
      tLat2 += relaxed_pg[ng].get_matrix(FRAC).transpose() * coord_trans_mat[FRAC].transpose() * coord_trans_mat[FRAC] * relaxed_pg[ng].get_matrix(FRAC);
    }

    tLat2 /= double(relaxed_pg.size());

    // tLat2 has the symmetrized lengths and angles -- it is equal to L.transpose()*L, where L=coord_trans_mat[FRAC]
    // we will find the sqrt of tLat2 and then reorient it so that it matches the original lattice
    Eigen::Matrix3d tMat(tLat2), tMat2;

    Eigen::JacobiSVD<Eigen::Matrix3d> tSVD(tMat);
    tMat = Eigen::Matrix3d::Zero();
    for(int i = 0; i < 3; i++)
      tMat(i, i) = tSVD.singularValues()[i];

    tMat2 = tSVD.matrixU() * tMat * tSVD.matrixV().transpose();

    tMat = coord_trans_mat[FRAC];

    tSVD.compute(tMat2 * tMat.transpose());

    tMat = tSVD.matrixV() * tSVD.matrixU().transpose() * tMat2;

    for(int i = 0; i < 3; i++) {
      for(int j = 0; j < 3; j++) {
        tLat2(i, j) = tMat(i, j);
      }
    }

    (*this) = Lattice(tLat2);

    calc_properties();
    calc_conversions();
    return;
  }

  //***********************************************************
  /**
   * Same as the other symmetrize routine, except this one assumes
   * that the symmetry group you mean to use is the point group
   * of your lattice within a certain tolerance.
   */
  //***********************************************************

  void Lattice::symmetrize(const double &tolerance) {
    SymGroup point_group;
    generate_point_group(point_group, tolerance);
    symmetrize(point_group);
    return;
  }

  //\John G

  //********************************************************************
  /**
     Linearly interpolates the lattices to trace the path between
     (*this) lattice and end_lattice.
  */
  //********************************************************************
  /*
  void Lattice::linear_interpolate(const Lattice &end_lattice, const int &num_images, Array<Lattice> &interp_lat) const {
    Matrix3<double> tlat_mat, inc_mat;
    tlat_mat = coord_trans_mat[FRAC];
    for(int i = 0; i < 3; i++) {
      for(int j = 0; j < 3; j++) {
        inc_mat(i, j) = (end_lattice.coord_trans_mat[FRAC](i, j) - coord_trans_mat[FRAC](i, j)) / (num_images + 1);
      }
    }
    interp_lat.push_back(Lattice(tlat_mat));
    for(int stepVal = 1; stepVal <= (num_images + 1); stepVal++) {
      tlat_mat = tlat_mat + inc_mat;
      interp_lat.push_back(Lattice(tlat_mat));
    }
  }
  */
  bool Lattice::is_right_handed() const {
    if(vol() < 0)
      return false;
    else
      return true;
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

  // read Lattice from a json array of Vector3<double>
  void from_json(Lattice &lat, const jsonParser &json) {
    try {
      lat = Lattice(json[0].get<Vector3<double> >(), json[1].get<Vector3<double> >(), json[2].get<Vector3<double> >());
    }
    catch(...) {
      /// re-throw exceptions
      throw;
    }
  };


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

    // rotate cell so that the lattice vectors are mostly aligned
    //    with Cartesian coordinate axes...
    for(int i = 0; i < point_grp.size(); i++) {

      Eigen::Matrix3d tmp = Eigen::Matrix3d(point_grp[i].get_matrix(CART)) * start;

      if(volume(Lattice(tmp)) < 0.0)
        continue;

      bool better = false;

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

    Matrix3<double> dA(lat2.inv_lat_column_mat()*lat1.lat_column_mat());
    Matrix3<long> iA;
    long N = 1, num, denom;
    //std::cout << "dA is:\n" << dA << "\n\n";
    for(Index i = 0; i < 3; i++) {
      for(Index j = 0; j < 3; j++) {
        nearest_rational_number(dA(i, j), num, denom);
        dA *= double(denom);
        N *= denom;
      }
    }
    for(Index i = 0; i < 3; i++) {
      for(Index j = 0; j < 3; j++) {
        iA(i, j) = round(dA(i, j));
      }
    }
    //std::cout << "iA is:\n" << iA << "\n\n";

    //std::cout << "and N is:\n" << N << "\n\n";
    Matrix3<long> U, S, V;
    iA.smith_normal_form(U, S, V);
    //std::cout << "Smith U is:\n" << U << "\n\n";
    //std::cout << "Smith S is:\n" << S << "\n\n";
    //std::cout << "Smith V is:\n" << V << "\n\n";
    //std::cout << "and U*S*V is:\n" << U*S*V << "\n\n";
    denom = N;

    //reuse matrix iA for matrix 'R', as above
    iA = 0;
    for(Index i = 0; i < 3; i++) {
      iA(i, i) = N / gcf(S(i, i), N);
    }
    //Reuse matrix iA for matrix 'N_1', as above
    iA = V.inverse() * iA;
    //std::cout << "N_1 is: \n" << iA << "\n\n";
    //Reuse matrix dA
    for(Index i = 0; i < 3; i++) {
      for(Index j = 0; j < 3; j++) {
        dA(i, j) = double(iA(i, j));
      }
    }

    Lattice tlat(lat1.lat_column_mat()*dA);
    return tlat.get_reduced_cell();

  }

  //*******************************************************************************************

  Lattice superdupercell(const std::vector<Lattice> &lat_list) {
    if(lat_list.size() == 0)
      return Lattice();

    Lattice tsupdup(lat_list[0]);
    for(Index i = 1; i < lat_list.size(); i++)
      tsupdup = superdupercell(tsupdup, lat_list[i]);
    return tsupdup;
  }




}



