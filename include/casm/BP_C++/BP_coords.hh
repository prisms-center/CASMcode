#ifndef BP_coords_HH
#define BP_coords_HH

#include "casm/BP_C++/BP_basic.hh"
#include "casm/BP_C++/BP_Vec.hh"
#include "casm/BP_C++/BP_Parse.hh"
#include "casm/external/Eigen/Dense"

// //////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////
//  Coordinate classes

namespace BP {

  /// \ingroup BP BP_coords
  class ucs_coord_class {
  public:
    long int coord[4];

    ucs_coord_class() {
      coord[0] = 0;
      coord[1] = 0;
      coord[2] = 0;
      coord[3] = 0;

    };

    /*ucs_coord_class( int i0, int i1, int i2, int i3)
    {
    	coord[0] = i0;
    	coord[1] = i1;
    	coord[2] = i2;
    	coord[3] = i3;

    };*/

    ucs_coord_class(long int i0, long int i1, long int i2, long int i3) {
      coord[0] = i0;
      coord[1] = i1;
      coord[2] = i2;
      coord[3] = i3;

    };

    long int &operator[](int i1) {
      return coord[i1];
    };

    const long int &operator[](int i1) const {
      return coord[i1];
    };

    bool operator==(ucs_coord_class B) const {
      if(coord[0] == B[0])
        if(coord[1] == B[1])
          if(coord[2] == B[2])
            if(coord[3] == B[3])
              return 1;

      return 0;
    };

    bool operator!=(ucs_coord_class B) const {
      if(coord[0] == B[0])
        if(coord[1] == B[1])
          if(coord[2] == B[2])
            if(coord[3] == B[3])
              return 0;

      return 1;
    };

    bool operator<(ucs_coord_class B) const {
      if(coord[0] < B[0])
        return true;
      else if(coord[0] == B[0]) {
        if(coord[1] < B[1])
          return true;
        else if(coord[1] == B[1]) {
          if(coord[2] < B[2])
            return true;

        }
      }

      return false;
    };

    ucs_coord_class operator+(ucs_coord_class B) {
      ucs_coord_class A(*this);

      A[0] += B[0];
      A[1] += B[1];
      A[2] += B[2];
      A[3] = B[3];

      return A;
    };

    ucs_coord_class operator+(ucs_coord_class B) const {
      ucs_coord_class A(*this);

      A[0] += B[0];
      A[1] += B[1];
      A[2] += B[2];
      A[3] = B[3];

      return A;
    };

    ucs_coord_class operator-(ucs_coord_class B) {
      ucs_coord_class A(*this);

      A[0] -= B[0];
      A[1] -= B[1];
      A[2] -= B[2];
      //A[3] unchanged

      return A;
    };

    ucs_coord_class operator-(ucs_coord_class B) const {
      ucs_coord_class A(*this);

      A[0] -= B[0];
      A[1] -= B[1];
      A[2] -= B[2];
      //A[3] unchanged

      return A;
    };

    ucs_coord_class operator*(long int i) const {
      ucs_coord_class A(*this);

      A[0] *= i;
      A[1] *= i;
      A[2] *= i;

      return A;
    };

    friend std::ostream &operator<<(std::ostream &outstream, const ucs_coord_class &c) {
      outstream << c[0] << " " << c[1] << " " << c[2] << " " << c[3] << " " ;
      return outstream;
    };

    friend std::istream &operator>>(std::istream &instream, ucs_coord_class &c) {
      instream >> c[0] >> c[1] >> c[2] >> c[3] ;
      return instream;
    };

    void write(BP_bWrite &fout) {
      fout.write_int(coord[0]);
      fout.write_int(coord[1]);
      fout.write_int(coord[2]);
      fout.write_int(coord[3]);

    };

    void read(BP_bParse &fin) {
      coord[0] = fin.next_int();
      coord[1] = fin.next_int();
      coord[2] = fin.next_int();
      coord[3] = fin.next_int();
    };

    void set(long int i0, long int i1, long int i2, long int i3) {
      coord[0] = i0;
      coord[1] = i1;
      coord[2] = i2;
      coord[3] = i3;

    };
  };


  /// \ingroup BP BP_coords
  class frac_coord_class {
  public:
    double coord[3];

    frac_coord_class() {
      coord[0] = 0;
      coord[1] = 0;
      coord[2] = 0;

    };

    frac_coord_class(double i0, double i1, double i2) {
      coord[0] = i0;
      coord[1] = i1;
      coord[2] = i2;

    };

    double &operator[](int i1) {
      if(i1 > 2) std::cout << "frac_coord index > 2" << std::endl;
      if(i1 < 0) std::cout << "frac_coord index < 0" << std::endl;
      return coord[i1];
    };

    const double &operator[](int i1) const {
      if(i1 > 2) std::cout << "frac_coord index > 2" << std::endl;
      if(i1 < 0) std::cout << "frac_coord index < 0" << std::endl;
      return coord[i1];
    };

    friend std::ostream &operator<<(std::ostream &outstream, const frac_coord_class &c) {
      outstream << c[0] << " " << c[1] << " " << c[2] << " " ;
      return outstream;
    }

    void set(double i0, double i1, double i2) {
      coord[0] = i0;
      coord[1] = i1;
      coord[2] = i2;

    };

    frac_coord_class operator+(frac_coord_class f) const {
      return frac_coord_class(coord[0] + f[0], coord[1] + f[1], coord[2] + f[2]);
    };

    frac_coord_class operator-(frac_coord_class f) const {
      return frac_coord_class(coord[0] - f[0], coord[1] - f[1], coord[2] - f[2]);
    };

    frac_coord_class operator*(double d) const {
      return frac_coord_class(coord[0] * d, coord[1] * d, coord[2] * d);
    };

    frac_coord_class operator/(double d) const {
      return frac_coord_class(coord[0] / d, coord[1] / d, coord[2] / d);
    };

    frac_coord_class &operator+=(frac_coord_class f) {
      coord[0] += f[0];
      coord[1] += f[1];
      coord[2] += f[2];

      return *this;
    };

    frac_coord_class &operator-=(frac_coord_class f) {
      coord[0] -= f[0];
      coord[1] -= f[1];
      coord[2] -= f[2];

      return *this;
    };

    bool equal(frac_coord_class f, double eps) {
      if(std::fabs((*this)[0] - f[0]) < eps)
        if(std::fabs((*this)[1] - f[1]) < eps)
          if(std::fabs((*this)[2] - f[2]) < eps) {
            return true;
          }

      return false;
    }

    bool operator<(frac_coord_class B) const {
      double eps = 1e-6;

      if(coord[0] < B[0] - eps)
        return true;
      else if(coord[0] < B[0] + eps) {
        if(coord[1] < B[1] - eps)
          return true;
        else if(coord[1] < B[1] + eps) {
          if(coord[2] < B[2] - eps)
            return true;

        }
      }

      return false;
    };

    double mag() {
      return sqrt(coord[0] * coord[0] + coord[1] * coord[1] + coord[2] * coord[2]);
    };

    void in_unit() {
      int i;
      for(i = 0; i < 3; i++)
        while(coord[i] < 0.0)
          coord[i] += 1.0;

      for(i = 0; i < 3; i++)
        while(coord[i] >= 1.0)
          coord[i] -= 1.0;
    };

    frac_coord_class in_unit() const {
      frac_coord_class c = *this;
      int i;
      for(i = 0; i < 3; i++)
        while(c[i] < 0.0)
          c[i] += 1.0;

      for(i = 0; i < 3; i++)
        while(c[i] >= 1.0)
          c[i] -= 1.0;
      return c;
    };

    ///		Returns frac_coord_class such that
    ///		- -0.5 <= c[i] <= 0.5
    ///		Warning, for skewed systems this may not be the shortest distance
    ///
    frac_coord_class pbc() const {	// this is for distances
      frac_coord_class c = *this;
      int i;
      for(i = 0; i < 3; i++)
        while(c[i] < -0.5)
          c[i] += 1.0;

      for(i = 0; i < 3; i++)
        while(c[i] > 0.5)
          c[i] -= 1.0;

      return c;
    };
  };

  /// \ingroup BP BP_coords
  class cart_coord_class {
  public:
    double coord[3];

    cart_coord_class() {
      coord[0] = 0;
      coord[1] = 0;
      coord[2] = 0;

    };

    cart_coord_class(double i0, double i1, double i2) {
      coord[0] = i0;
      coord[1] = i1;
      coord[2] = i2;

    };

    double &operator[](int i1) {
      if(i1 > 2) std::cout << "cart_coord index > 2" << std::endl;
      if(i1 < 0) std::cout << "cart_coord index < 0" << std::endl;
      return coord[i1];
    };

    const double &operator[](int i1) const {
      if(i1 > 2) std::cout << "cart_coord index > 2" << std::endl;
      if(i1 < 0) std::cout << "cart_coord index < 0" << std::endl;
      return coord[i1];
    };

    friend std::ostream &operator<<(std::ostream &outstream, const cart_coord_class &c) {
      outstream << c[0] << " " << c[1] << " " << c[2] << " " ;
      return outstream;
    }

    void set(double i0, double i1, double i2) {
      coord[0] = i0;
      coord[1] = i1;
      coord[2] = i2;

    };

    double mag() const {
      return sqrt(coord[0] * coord[0] + coord[1] * coord[1] + coord[2] * coord[2]);
    };

    cart_coord_class norm() const {
      cart_coord_class coord_b = *this;
      return coord_b / coord_b.mag();
    };

    double dot(const cart_coord_class &c) const {
      return coord[0] * c[0] + coord[1] * c[1] + coord[2] * c[2];
    };

    cart_coord_class cross(const cart_coord_class &b) const {
      cart_coord_class c;
      c[0] = coord[1] * b[2] - b[1] * coord[2];
      c[1] = -(coord[0] * b[2] - b[0] * coord[2]);
      c[2] = coord[0] * b[1] - b[0] * coord[1];

      return c;
    };

    cart_coord_class operator+(cart_coord_class c) const {
      return cart_coord_class(coord[0] + c[0], coord[1] + c[1], coord[2] + c[2]);
    };

    cart_coord_class operator-(cart_coord_class c) const {
      return cart_coord_class(coord[0] - c[0], coord[1] - c[1], coord[2] - c[2]);
    };

    cart_coord_class &operator+=(cart_coord_class c) {
      coord[0] += c[0];
      coord[1] += c[1];
      coord[2] += c[2];

      return *this;

    };

    cart_coord_class &operator-=(cart_coord_class c) {
      coord[0] -= c[0];
      coord[1] -= c[1];
      coord[2] -= c[2];

      return *this;
    };

    cart_coord_class operator*(int x) const {
      return cart_coord_class(coord[0] * x, coord[1] * x, coord[2] * x);
    };

    cart_coord_class operator*(long int x) const {
      return cart_coord_class(coord[0] * x, coord[1] * x, coord[2] * x);
    };

    cart_coord_class operator*(unsigned long int x) const {
      return cart_coord_class(coord[0] * x, coord[1] * x, coord[2] * x);
    };

    cart_coord_class operator*(double x) const {
      return cart_coord_class(coord[0] * x, coord[1] * x, coord[2] * x);
    };

    cart_coord_class operator/(double x) const {
      return cart_coord_class(coord[0] / x, coord[1] / x, coord[2] / x);
    };
  };

  /// This class tracks the number of periodic boundary crossings.
  /**
  *	pbx_class::val is incremented when an atom crosses a periodic boundary in a positive direction, and
  *	decremented when a periodic boundary is crossed in the negative direction.
  *
  *	\ingroup BP BP_coords
  */

  class pbx_class {
  public:
    int val[3];			///< The number of periodic boundary crossings along the directions corresponding to prim_class::pr_v

    pbx_class() {
      val[0] = val[1] = val[2] = 0;
    };

    pbx_class(int i0, int i1, int i2) {
      val[0] = i0;
      val[1] = i1;
      val[2] = i2;
    };

    inline int &operator[](int i) {
      return val[i];
    };

    inline int operator[](int i) const {
      return val[i];
    };

    ///		Sets pbx_class::val[0] = i0; pbx_class::val[1] = i1; and pbx_class::val[0] = i1;
    void set(int i0, int i1, int i2) {
      val[0] = i0;
      val[1] = i1;
      val[2] = i2;

    };

    void set(pbx_class p2) {
      val[0] = p2[0];
      val[1] = p2[1];
      val[2] = p2[2];

    };

    inline pbx_class operator+(pbx_class p2) {
      pbx_class p1(*this);
      p1.val[0] += p2.val[0];
      p1.val[1] += p2.val[1];
      p1.val[2] += p2.val[2];
      return p1;
    };

    inline void operator+=(pbx_class p2) {
      val[0] += p2.val[0];
      val[1] += p2.val[1];
      val[2] += p2.val[2];
    };

    inline pbx_class operator-(pbx_class p2) {
      pbx_class p1(*this);
      p1.val[0] -= p2.val[0];
      p1.val[1] -= p2.val[1];
      p1.val[2] -= p2.val[2];
      return p1;
    };

    inline void operator-=(pbx_class p2) {
      val[0] -= p2.val[0];
      val[1] -= p2.val[1];
      val[2] -= p2.val[2];
    };

    bool operator==(pbx_class p2) const {
      if(val[0] == p2[0])
        if(val[1] == p2[1])
          if(val[2] == p2[2])
            return true;

      return false;
    };

    bool operator!=(pbx_class p2) const {
      if(val[0] == p2[0])
        if(val[1] == p2[1])
          if(val[2] == p2[2])
            return false;

      return true;
    };

    friend std::ostream &operator<<(std::ostream &outstream, const pbx_class &p) {
      outstream << p[0] << " " << p[1] << " " << p[2] << " " ;
      return outstream;
    };

  };




  /// This class gives the position of an atom by its site and number of periodic boundary crossings.
  /// \ingroup BP BP_coords
  ///
  class position_class {
  public:
    ucs_coord_class	site;
    pbx_class pbx;						//# times crossed the periodic boundaries

    friend std::ostream &operator<<(std::ostream &outstream, const position_class &p) {
      outstream << p.site << p.pbx[0] << " " << p.pbx[1] << " " << p.pbx[2] << " " ;
      return outstream;
    };

    friend std::istream &operator>>(std::istream &instream, position_class &p) {
      instream >> p.site >> p.pbx[0] >> p.pbx[1] >> p.pbx[2] ;
      return instream;
    }

    position_class() {
      site.set(-1, -1, -1, -1);
      pbx[0] = pbx[1] = pbx[2] = 0;
    };

    /**
    *	Constructor. Sets position_class::site = i1;
    */
    position_class(ucs_coord_class i1) {
      site = i1;
      pbx[0] = pbx[1] = pbx[2] = 0;
    };

    /**		Constructor. Sets position_class::site = i1;
    /		position_class::pbx[0] = i1; position_class::pbx[1] = i2; position_class::pbx[2] = i3;
    */
    position_class(ucs_coord_class i1, int i2, int i3, int i4) {
      site = i1;
      pbx[0] = i2;
      pbx[1] = i3;
      pbx[2] = i4;
    };

    /**		Constructor. Sets position_class::site = i1;
    /		position_class::pbx[0] = pbx;
    */
    position_class(ucs_coord_class i1, pbx_class i2) {
      site = i1;
      pbx = i2;
    };

    ///		Sets all values to 0.
    ///
    void clear() {
      site.set(0, 0, 0, 0);
      pbx[0] = pbx[1] = pbx[2] = 0;
    };

    void write(BP_bWrite &fout) {
      site.write(fout);
      fout.write_int(pbx[0]);
      fout.write_int(pbx[1]);
      fout.write_int(pbx[2]);

    };

    void read(BP_bParse &fin) {
      site.read(fin);

      pbx[0] = fin.next_int();
      pbx[1] = fin.next_int();
      pbx[2] = fin.next_int();
    };

    bool operator==(position_class &P) const {
      if(site == P.site)
        if(pbx == P.pbx)
          return true;

      return false;

    };

    bool operator!=(position_class &P) const {
      if(site == P.site)
        if(pbx == P.pbx)
          return false;

      return true;

    };
  };



  // //////////////////////////
  // Bonus functions

  // declarations:

  ///		Reads a ucs_coord_class object from a std::string
  std::string &operator>>(std::string &s, ucs_coord_class &c);

  ///		Reads a position_class object from a std::string
  std::string &operator>>(std::string &s, position_class &p);


  //rotations
  BP_Vec<cart_coord_class> &rotate_x(BP_Vec<cart_coord_class> &v_list, double theta);	///< \ingroup BP_coords
  BP_Vec<cart_coord_class> &rotate_y(BP_Vec<cart_coord_class> &v_list, double theta);	///< \ingroup BP_coords
  BP_Vec<cart_coord_class> &rotate_z(BP_Vec<cart_coord_class> &v_list, double theta);	///< \ingroup BP_coords
  BP_Vec<cart_coord_class> &rotate(BP_Vec<cart_coord_class> &v_list, BP_Vec< BP_Vec<double> > &R);	///< \ingroup BP_coords
  BP_Vec< BP_Vec<double> >  rotation_matrix_x(double theta);	///< \ingroup BP_coords
  BP_Vec< BP_Vec<double> >  rotation_matrix_y(double theta);	///< \ingroup BP_coords
  BP_Vec< BP_Vec<double> >  rotation_matrix_z(double theta);	///< \ingroup BP_coords
  template<class T> cart_coord_class operator*(const BP_Vec< BP_Vec<T> > &mA, const cart_coord_class &vB);

  //conversion
  cart_coord_class frac_to_cart(BP_Vec<cart_coord_class> &sc_v, frac_coord_class f);		///< \ingroup BP_coords
  frac_coord_class cart_to_frac(BP_Vec<cart_coord_class> &sc_v, cart_coord_class c);		///< \ingroup BP_coords

  //get reciprocal lattice vectors:
  BP_Vec<cart_coord_class> recip_vec(const BP_Vec<cart_coord_class> &unit_v);	///< \ingroup BP_coords

  // template definitions:

  template<class T> cart_coord_class operator*(const BP_Vec< BP_Vec<T> > &mA, const cart_coord_class &vB) {
    unsigned long int r, n, i, k;
    r = mA.size();
    n = mA[0].size();

    cart_coord_class v(0, 0, 0);

    for(i = 0; i < r; i++)
      for(k = 0; k < n; k++)
        v[i] += mA[i][k] * vB[k];

    return v;
  }

}

#endif // BP_coords_HH
