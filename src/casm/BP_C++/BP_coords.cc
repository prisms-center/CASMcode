#ifndef BP_coords_CC
#define BP_coords_CC


#include "casm/BP_C++/BP_coords.hh"

// //////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////
//  Coordinate functions

namespace BP {

  ///		Reads a ucs_coord_class object from a std::string
  std::string &operator>>(std::string &s, ucs_coord_class &c) {
    s >> c[0] >> c[1] >> c[2] >> c[3];
    return s;
  };

  ///		Reads a position_class object from a std::string
  std::string &operator>>(std::string &s, position_class &p) {
    s >> p.site >> p.pbx[0] >> p.pbx[1] >> p.pbx[2];
    return s;
  };

  // //////////////////////////
  // Bonus functions

  // definitions:
  //rotations
  BP_Vec<cart_coord_class> &rotate_x(BP_Vec<cart_coord_class> &v_list, double theta) {
    BP_Vec< BP_Vec<double> > R = BP_Vec< BP_Vec<double> >(3, BP_Vec<double>(3, 0));
    R[0][0] = 1.0;
    R[1][1] = cos(theta);
    R[2][1] = sin(theta);
    R[1][2] = -sin(theta);
    R[2][2] = cos(theta);

    return rotate(v_list, R);

  };

  BP_Vec<cart_coord_class> &rotate_y(BP_Vec<cart_coord_class> &v_list, double theta) {
    BP_Vec< BP_Vec<double> > R = BP_Vec< BP_Vec<double> >(3, BP_Vec<double>(3, 0));
    R[1][1] = 1.0;
    R[0][0] = cos(theta);
    R[0][2] = sin(theta);
    R[2][0] = -sin(theta);
    R[2][2] = cos(theta);

    return rotate(v_list, R);

  };

  BP_Vec<cart_coord_class> &rotate_z(BP_Vec<cart_coord_class> &v_list, double theta) {
    BP_Vec< BP_Vec<double> > R = BP_Vec<BP_Vec<double> >(3, BP_Vec<double>(3, 0));
    R[2][2] = 1.0;
    R[0][0] = cos(theta);
    R[1][0] = sin(theta);
    R[0][1] = -sin(theta);
    R[1][1] = cos(theta);

    return rotate(v_list, R);

  };

  BP_Vec<cart_coord_class> &rotate(BP_Vec<cart_coord_class> &v_list, BP_Vec< BP_Vec<double> > &R) {

    cart_coord_class tmp_coord;

    for(int i = 0; i < v_list.size(); i++) {
      tmp_coord.set(0, 0, 0);
      tmp_coord[0] = R[0][0] * v_list[i][0] + R[0][1] * v_list[i][1] + R[0][2] * v_list[i][2];
      tmp_coord[1] = R[1][0] * v_list[i][0] + R[1][1] * v_list[i][1] + R[1][2] * v_list[i][2];
      tmp_coord[2] = R[2][0] * v_list[i][0] + R[2][1] * v_list[i][1] + R[2][2] * v_list[i][2];
      v_list[i] = tmp_coord;
    }

    return v_list;

  };

  BP_Vec< BP_Vec<double> >  rotation_matrix_x(double theta) {
    BP_Vec< BP_Vec<double> > R = BP_Vec< BP_Vec<double> >(3, BP_Vec<double>(3, 0));
    R[0][0] = 1.0;
    R[1][1] = cos(theta);
    R[2][1] = sin(theta);
    R[1][2] = -sin(theta);
    R[2][2] = cos(theta);

    return R;

  };

  BP_Vec< BP_Vec<double> >  rotation_matrix_y(double theta) {
    BP_Vec< BP_Vec<double> > R = BP_Vec< BP_Vec<double> >(3, BP_Vec<double>(3, 0));
    R[1][1] = 1.0;
    R[0][0] = cos(theta);
    R[0][2] = sin(theta);
    R[2][0] = -sin(theta);
    R[2][2] = cos(theta);

    return R;

  };

  BP_Vec< BP_Vec<double> >  rotation_matrix_z(double theta) {
    BP_Vec< BP_Vec<double> > R = BP_Vec<BP_Vec<double> >(3, BP_Vec<double>(3, 0));
    R[2][2] = 1.0;
    R[0][0] = cos(theta);
    R[1][0] = sin(theta);
    R[0][1] = -sin(theta);
    R[1][1] = cos(theta);

    return R;

  };

  //get reciprocal lattice vectors:
  BP_Vec<cart_coord_class> recip_vec(const BP_Vec<cart_coord_class> &unit_v) {
    BP_Vec<cart_coord_class> recip_v;
    recip_v.capacity(3);

    double f = 2 * PI / unit_v[0].dot(unit_v[1].cross(unit_v[2]));
    recip_v.add(unit_v[1].cross(unit_v[2])*f);
    recip_v.add(unit_v[2].cross(unit_v[0])*f);
    recip_v.add(unit_v[0].cross(unit_v[1])*f);

    return recip_v;
  };


  //conversion
  cart_coord_class frac_to_cart(BP_Vec<cart_coord_class> &sc_v, frac_coord_class f) {
    // sc_v is (supercell vectors) (vector of three cart_coord_class objects)

    double r_c[3];
    for(int k = 0; k < 3; k++)
      r_c[k] = f[0] * sc_v[0][k] + f[1] * sc_v[1][k] + f[2] * sc_v[2][k];	//cartesian coords
    return cart_coord_class(r_c[0], r_c[1], r_c[2]);

  };

  frac_coord_class cart_to_frac(BP_Vec<cart_coord_class> &sc_v, cart_coord_class c) {
    // sc_v is (supercell vectors) (vector of three cart_coord_class objects)

    int i, j;
    // supercell vectors
    Eigen::Matrix3d M, Minv;
    for(i = 0; i < 3; i++) {
      for(j = 0; j < 3; j++) {
        M(j, i) = sc_v[i][j];
      }
    }

    // inverse supercell vectors
    Minv = M.inverse();
    //for( i=0; i<3; i++)
    //{
    //	for( j=0; j<3; j++)
    //	{
    //		sc_inv[i][j] = Minv(j,i);
    //	}
    //}

    frac_coord_class f;
    for(int k = 0; k < 3; k++)
      //	f[k] = c[0]*sc_inv[0][k] + c[1]*sc_inv[1][k] + c[2]*sc_inv[2][k];
      f[k] = c[0] * Minv(k, 0) + c[1] * Minv(k, 1) + c[2] * Minv(k, 2);

    return f;


  };


}

#endif // BP_coords_CC
