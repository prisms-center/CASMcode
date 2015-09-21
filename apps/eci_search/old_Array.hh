#ifndef ARRAY_CLASS
#define ARRAY_CLASS
//This code is for Matrix operation and is still under development. Qingchuan Xu
//when you claim a matrix m, use "array m(r,c)", r and c are number of rows and number of columns.
//there are two methods to initialize the matrix elements, the first is to use "cin>>m", input
//element value from screen, the 2nd is to set a array "double n[]={......}", then call member
//function "setArray(n)".
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <stack>
#include <queue>
#include <ctime>
# define TINY 1.0e-20;
#define SWAP(g,h) {y=(g);(g)=(h);(h)=y;}
#define SIGN(a,b) ((b) >=0.0 ? fabs(a) : -fabs(a))
//SIGN(a,b) makes the sign of a is the same as the sign of b.
static double maxarg1, maxarg2;
#define FMAX(a,b) (maxarg1=(a), maxarg2=(b),(maxarg1) > (maxarg2) ? (maxarg1) : (maxarg2))
static int iminarg1, iminarg2;
#define IMIN(a,b) (iminarg1=(a), iminarg2=(b), (iminarg1) < (iminarg2) ? (iminarg1) : (iminarg2))
#define TOL 2.220446049250313e-16
using namespace std;

namespace QXu {

  class Array {
    friend ostream &operator<<(ostream &, const Array &);  //overload operator <<
    friend istream &operator>>(istream &, Array &);        //overload operator >>
    friend Array operator * (double, Array);              //overload *, do number*matrix
    friend Array operator * (Array , double);             //overload *, do matrix*number
    friend double distance_pl(Array p3, Array p1, Array p2);//calculate the distance between point p3 and straight line
    //p1p2
    friend bool in_tri(Array p4, Array p1, Array p2, Array p3); //check whether point p4 is within triangle p1p2p3 or not
    friend bool on_plane(Array p4, Array p1, Array p2, Array p3);//check whether point p4 is on the plane p1-p2-p3 or not
    friend bool same_side(Array p3, Array p4, Array p1, Array p2);//check whether point p3 and p4 are to the same side of line p1-p2
    friend bool cut(Array p3, Array p4, Array p5, Array p1, Array p2);//check whether triangle p3-p4-p5 is cut by line p1-p2
  public:
    Array(int = 3, int = 3);             // default constructor
    Array(const Array &);                // copy constructor
    ~Array();                            // destructor
    int getSize() const;                 // return size
    const Array &operator=(const Array &);   //overload =, assign an known matrix to unknown matrix
    bool operator==(const Array &) const;   // overload ==, compare equal
    Array operator + (Array);              //overload +, do matrix+matrix
    Array operator - (Array);              //overload -, do matrix-matrix
    Array operator * (Array);              //overload *, do matrix*matrix
    double det();                           //calculate determinant using definition
    double detlu();                         //calculate determinant using LU decomposition
    Array inverse();                        //calculate inverse matrix using definition
    Array inverLU();                        //calculate inverse matrix using LU decomposition
    Array inverSVD();                       //calculate inverse matrix using SVD
    Array transpose();                      //transpose matrix
    double tr();                            //trace of a matrix
    Array LU();                             //do LU decomposion of a matrix
    Array eig();                            //calculate eigenvalue of a matrix
    Array solveleq();                       //slove linear equations,using Gaussian elimination with full pivot
    void  svd(Array &, Array &, Array &);   //Singular Value Decompose of a matrix
    Array svdeq();                          //Use SVD to solve AX=B for vector X
    Array svdls(Array &, Array &, double &, double &, int &);  //multidimensional least square fit using SVD, cv is calculated by analytical solution of Leave One Out
    double cv(Array &, int &);                             //calculate cv score using definition of Leave One Out cv score
    Array svdls_MCCV(Array &, Array &, double &, double &, int , int, int &);
    //multidimensional least square fit using SVD, cv is calculated by using Monte Carlo CV score
    double MCCV(Array &, int, int, int &);          //calculate cv using Monte Carlo cv method
    double MCCV_formula(Array &, int, int, Array &); //calculate cv using Monte Carlo method and formula
    Array svdls_core(Array &);               //do least square fit using svd, just return fitted values
    Array hull2d(int &);                         //find the hull of a set of points(2 dimension), int is the number of hull points
    void exrow(int i, int j);                      // exchange ith row and jth row, i and j are from 1
    void excol(int i, int j);                      // exchange ith column and jth column, i and j are from 1
    void assort_a(int i);                  // assort data according to the ith column, in ascending order, i is from 1
    void assort_d(int j);                  // assort data according to the ith column, in decending order, i is from 1
    double useful(double, double);          //a useful function used by svd()
    Array hull_base_p(int n);                 // find hull base of n-dimension case, n=2, or 3, only hull points
    Array hullnd_p(int n);                    //find hull of n-dimension case, n=2, or 3, only hull points
    Array hull_base(int n);                 // find hull base of n-dimension case, n=2, or 3
    Array hullnd(int n, Array &f);                    //find hull of n-dimension case, n=2, or 3
    void setArray(double *a);               //initialize matrix element, using array a
    bool operator!=(const Array &right) const { // Determine if two arrays are not equal and
      return !(*this == right);  // return true, otherwise return false (uses operator==).
    }
    double elem(int i, int j) {
      return ptr[(i - 1) * nc + j - 1]; // get the matrix element (i,j), note that i and j
    }
    // are from 1, not from 0.

    void set_elem(int i, int j, double val) {
      ptr[(i - 1)*nc + j - 1] = val; //note that i and j
    }
    // are from 1, not from 0.

    double zelem(int i, int j) {
      return ptr[(i) * nc + j]; // get the matrix element (i,j), note that i and j
    }
    // are from 0, not from 1.

    void set_zelem(int i, int j, double val) {
      ptr[(i)*nc + j] = val; //note that i and j
    }
    // are from 0, not from 1.



    Array row(int i);                          // get the ith row of matrix, i is from 1, not 0.
    Array col(int j);                          // get the jth row of matrix, j is from 1, not 0.
    int num_row() {
      return nl; // return number of rows
    }
    int num_col() {
      return nc; // return number of columns
    }
    double modu() {                            // compute module of a vector(nl,1), sum over squared elements then
      // square root of the sum.
      double modu_valu = 0.0;
      for(int i = 0; i < nl; i++) modu_valu += ptr[i] * ptr[i];
      modu_valu = sqrt(modu_valu);
      return modu_valu;
    }
    void identity() {                          // set identity matrix (nl, nl)
      int i;
      for(i = 0; i < nl * nc; i++) ptr[i] = 0.0;
      for(i = 0; i < nl; i++) ptr[i * nc + i] = 1.0;
    }

  private:
    int nl, nc; // size of the array, nl is number of row and nc is number of column
    double *ptr; // pointer to first element of array

  };


  // Member function definitions for class Array

  // Default constructor for class Array (default size 3x3)
  Array::Array(int arrayL, int arrayC) {
    nl = (arrayL > 0 ? arrayL : 3);
    nc = (arrayC > 0 ? arrayC : 3);
    ptr = new double[ nl * nc ]; // create space for array
    for(int i = 0; i < nl * nc; i++)
      ptr[ i ] = 0;           // initialize array
  }

  // Copy constructor for class Array
  // must receive a reference to prevent infinite recursion
  Array::Array(const Array &init) : nl(init.nl), nc(init.nc) {
    ptr = new double[ nl * nc ]; // create space for array
    for(int i = 0; i < nl * nc; i++)
      ptr[ i ] = init.ptr[ i ];  // copy init into object
  }

  // Destructor for class Array
  Array::~Array() {
    delete [] ptr;            // reclaim space for array
  }

  // Get the size of the array
  int Array::getSize() const {
    cout << "the dimension is " << nl << "x" << nc << endl;
    return nl * nc;
  }

  //matrix addition, matrix+matrix, return matrix
  Array Array::operator +(Array A2) {
    Array Asum(nl, nc);
    for(int i = 0; i < nl * nc; i++)
      Asum.ptr[i] = ptr[i] + A2.ptr[i];
    return Asum;
  }
  //matrix substraction, matrix-matrix, return matrix
  Array Array::operator -(Array A) {
    Array Asub(nl, nc);
    for(int i = 0; i < nl * nc; i++)
      Asub.ptr[i] = ptr[i] - A.ptr[i];
    return Asub;
  }
  //matrix product, matrix*matrix, return matrix
  Array Array::operator *(Array Ap) {
    Array Aprod(nl, Ap.nc);
    for(int i = 0; i < nl; i++)
      for(int j = 0; j < Ap.nc; j++) {
        double cc = 0.0;
        for(int k = 0; k < nc; k++)
          cc = cc + (ptr[i * nc + k]) * (Ap.ptr[k * (Ap.nc) + j]);
        Aprod.ptr[i * Ap.nc + j] = cc;
      }
    return Aprod;
  }
  //product of matrix and a scalar, retrun matrix
  Array operator *(double a, Array m) {
    Array ASP(m.nl, m.nc);
    for(int i = 0; i < m.nl; i++)
      for(int j = 0; j < m.nc; j++)
        ASP.ptr[i * m.nc + j] = a * m.ptr[i * m.nc + j];
    return ASP;
  }
  Array operator *(Array m, double a) {
    Array ASP(m.nl, m.nc);
    for(int i = 0; i < m.nl; i++)
      for(int j = 0; j < m.nc; j++)
        ASP.ptr[i * m.nc + j] = a * m.ptr[i * m.nc + j];
    return ASP;
  }


  //***********************************************
  //matrix determinant, using definition
  double Array::det() {
    int i, j, j1, j2;
    double deter;
    Array subm(nl - 1, nl - 1);
    if(nl == 1) {
      deter = ptr[0];
    }
    else if(nl == 2) {
      deter = ptr[0] * ptr[3] - ptr[2] * ptr[1];
    }
    else {
      deter = 0;
      for(j1 = 0; j1 < nl; j1++) {
        for(i = 1; i < nl; i++) {
          j2 = 0;
          for(j = 0; j < nl; j++) {
            if(j == j1) continue;
            subm.ptr[(i - 1)*subm.nc + j2] = ptr[i * nc + j];
            j2++;
          }
        }
        deter += pow(-1.0, j1 + 1.0 + 1.0) * ptr[j1] * (subm.det());
      }
    }
    return(deter);
  }
  //*****************************************************************
  //matrix determinant using LU method
  double Array::detlu() {
    Array ArrayLU(nl, nc);
    double pro = 1.0;
    ArrayLU = this->LU();
    for(int i = 0; i < nl; i++)
      pro *= ArrayLU.ptr[i * nc + i];
    return pro;
  }
  //*********************************************************************
  //matrix inverse, using definition, return matrix
  Array Array::inverse() {
    int ii, jj, ii1, jj1, ii2, jj2;
    Array isubm(nl - 1, nc - 1), invm(nl, nc), inv(nl, nc);
    double ib;
    if(det() == 0) {
      cout << "singular!" << endl;
      exit(1);
    }
    if(nl == 1) {
      inv.ptr[0] = 1.0 / ptr[0];
      return inv;
    }
    else {
      for(ii1 = 0; ii1 < nl; ii1++)
        for(jj1 = 0; jj1 < nl; jj1++) {
          ii2 = 0;
          for(ii = 0; ii < nl; ii++) {
            if(ii == ii1) continue;
            jj2 = 0;
            for(jj = 0; jj < nl; jj++) {
              if(jj == jj1) continue;
              isubm.ptr[ii2 * (nc - 1) + jj2] = ptr[ii * nc + jj];
              jj2++;
            }
            ii2++;
          }
          ib = isubm.det();
          invm.ptr[jj1 * nc + ii1] = pow(-1.0, ii1 + jj1) * ib;
        }
      inv = (1 / det()) * invm;
      return inv;
    }
  }
  //***************************************************************
  //matrix inverse using LU decompose
  Array Array::inverLU() {
    Array arrlu(nl, nc), inver(nl, nc);
    double *b = new double[nl];
    double *y = new double[nl];
    double *x = new double[nl];
    int i, j, k;
    double sum;
    arrlu = this->LU();
    for(j = 0; j < nc; j++) {
      for(i = 0; i < nl; i++) b[i] = 0.0;
      b[j] = 1.0;
      y[0] = b[0];
      for(i = 1; i < nl; i++) {
        sum = b[i];
        for(k = 0; k <= i - 1; k++)
          sum = sum - y[k] * arrlu.ptr[i * nc + k];
        y[i] = sum;
      }
      x[nl - 1] = y[nl - 1] / arrlu.ptr[(nl - 1) * nc + (nl - 1)];
      for(i = nl - 2; i >= 0; i--) {
        sum = y[i];
        for(k = i + 1; k < nl; k++)
          sum = sum - x[k] * arrlu.ptr[i * nc + k];
        x[i] = sum / arrlu.ptr[i * nc + i];
      }
      for(i = 0; i < nl; i++)
        inver.ptr[i * nc + j] = x[i];
    }

    return inver;
  }
  //****************************************************************

  Array Array::transpose() {
    Array transm(nc, nl);
    for(int i = 0; i < nc; i++)
      for(int j = 0; j < nl; j++)
        transm.ptr[i * transm.nc + j] = ptr[j * nc + i];
    return transm;
  }
  //***********************************************************
  //calculate the trace of a matrix
  double Array::tr() {
    double trace;
    trace = 0;
    for(int i = 0; i < nl; i++) trace += ptr[i * nc + i];
    return trace;
  }

  //***********************************************
  //LU decompose of matrix. The lower part is restored in the lower part of the returned matrix and the diagonal
  //elements are 1s, the upper part is restored in the upper part of the returned matrix. The diagonal elements
  //of the input matrix should not be zero
  Array Array::LU() {
    Array alu(nl, nc);
    int i, j, k;
    double  sum, dum;
    for(int i = 0; i < nl * nc; i++)
      alu.ptr[ i ] = ptr[ i ];
    for(j = 0; j < nl; j++) {
      for(i = 0; i < j; i++) {
        sum = alu.ptr[i * nc + j];
        for(k = 0; k < i; k++) sum -= (alu.ptr[i * nc + k]) * (alu.ptr[k * nc + j]);
        alu.ptr[i * nc + j] = sum;
      }
      for(i = j; i < nl; i++) {
        sum = alu.ptr[i * nc + j];
        for(k = 0; k < j; k++) sum -= (alu.ptr[i * nc + k]) * (alu.ptr[k * nc + j]);
        alu.ptr[i * nc + j] = sum;

      }
      if(alu.ptr[j * nc + j] == 0.0) {
        alu.ptr[j * nc + j] = TINY;
        cout << "diagonal element is zero" << endl;
      }
      if(j != nl) {
        dum = 1.0 / (alu.ptr[j * nc + j]);
        for(i = j + 1; i < nl; i++) alu.ptr[i * nc + j] *= dum;
      }

    }

    return alu;
  }

  //****************************************************************************************************
  //calculate eigenvalue of a matrix, return a 2-row matrix, the 1st row contains the real parts of the
  //eigenvalus and the 2nd row contains the imaginary parts of the eigenvalus
  Array Array::eig() {
    double *back = new double[nl * nc];
    double *a = new double[nl * nc + 1];
    double *wr = new double[nl + 1];
    double *wi = new double[nl + 1];
    Array eigenvalue(2, nc);
    wr[0] = wi[0] = 0.0;
    int last, j, i, m, nn, l, k, its, mmin;
    double s, r, g, f, c, sqrdx, y, x, z, w, v, u, t, q, p, anorm;
    for(i = 0; i < nl * nc; i++)
      back[i] = ptr[i];
    // this balance algorithm is from numerical recipes in C. William H. Press, the original matrix is balanced.
    sqrdx = 2.0 * 2.0;
    last = 0;
    while(last == 0) {
      last = 1;
      for(i = 0; i < nl; i++) {
        r = c = 0.0;
        for(j = 0; j < nl; j++)
          if(j != i) {
            c += fabs(ptr[j * nc + i]);
            r += fabs(ptr[i * nc + j]);
          }
        if(c && r) {
          g = r / 2.0;
          f = 1.0;
          s = c + r;
          while(c < g) {
            f *= 2.0;
            c *= sqrdx;
          }
          g = r * 2.0;
          while(c > g) {
            f /= 2.0;
            c /= sqrdx;
          }
          if((c + r) / f < 0.95 * s) {
            last = 0;
            g = 1.0 / f;
            for(j = 0; j < nl; j++) ptr[i * nc + j] *= g;
            for(j = 0; j < nl; j++) ptr[j * nc + i] *= f;
          }
        }
      }
    }

    //the following code reduce a balanced matrix to upper Hessenberg Form with identical eigenvalues. the Hessenberg
    //matrix is in elements ptr[i*nc+j] with i<=j+1. Elements with i>j+1 are to be thought of as zero, but are returned
    //with random values.We use an algorithm analogous to Gaussian elimination wit pivoting which is more efficient than
    //the Householder method.
    for(m = 1; m < nl - 1; m++) {
      x = 0.0;
      i = m;
      for(j = m; j < nl; j++) {
        if(fabs(ptr[j * nc + m - 1]) > fabs(x)) {
          x = ptr[j * nc + m - 1];
          i = j;
        }
      }
      if(i != m) {
        for(j = m - 1; j < nl; j++) SWAP(ptr[i * nc + j], ptr[m * nc + j])
          for(j = 0; j < nl; j++) SWAP(ptr[j * nc + i], ptr[j * nc + m])
          }
      if(x) {
        for(i = m + 1; i < nl; i++) {
          if((y = ptr[i * nc + m - 1]) != 0.0) {
            y /= x;
            ptr[i * nc + m - 1] = y;
            for(j = m; j < nl; j++)
              ptr[i * nc + j] -= y * ptr[m * nc + j];
            for(j = 0; j < nl; j++)
              ptr[j * nc + m] += y * ptr[j * nc + i];
          }
        }
      }
    }
    // the following code finds all eigenvalues of an upper Hessenberg matrix using QR algorithm with
    // shifts On input the matrix can be exactly as output from the above code, on output it is
    //destroyed. The real and imaginary parts of the eigenvalues are returned in wr[0...nl-1] and
    //wi[0...nl-1], respectively.
    a[0] = 1.0;
    for(i = 0; i < nl * nc; i++)
      a[i + 1] = ptr[i];
    anorm = fabs(a[1]);
    for(i = 2; i <= nl; i++)
      for(j = (i - 1); j <= nl; j++)
        anorm += fabs(a[i * nc + j - nc]);
    nn = nl;
    t = 0.0;
    while(nn >= 1) {
      its = 0;
      do {
        for(l = nn; l >= 2; l--) {
          s = fabs(a[(l - 1) * nc + l - 1 - nc]) + fabs(a[l * nc + l - nc]);
          if(s == 0.0) s = anorm;
          if((double)(fabs(a[l * nc + l - 1 - nc]) + s) == s) break;
        }
        x = a[nn * nc + nn - nc];
        if(l == nn) {
          wr[nn] = x + t;
          wi[nn--] = 0.0;
        }
        else {
          y = a[(nn - 1) * nc + nn - 1 - nc];
          w = a[nn * nc + nn - 1 - nc] * a[(nn - 1) * nc + nn - nc];
          if(l == (nn - 1)) {
            p = 0.5 * (y - x);
            q = p * p + w;
            z = sqrt(fabs(q));
            x += t;
            if(q >= 0.0) {
              z = p + SIGN(z, p);
              wr[nn - 1] = wr[nn] = x + z;
              if(z) wr[nn] = x - w / z;
              wi[nn - 1] = wi[nn] = 0.0;
            }
            else {
              wr[nn - 1] = wr[nn] = x + p;
              wi[nn - 1] = -(wi[nn] = z);
            }
            nn -= 2;
          }
          else {
            if(its == 30) cout << "Too many iterations in hqr" << endl;
            if(its == 10 || its == 20) {
              t  += x;
              for(i = 1; i <= nn; i++) a[i * nc + i - nc] -= x;
              s = fabs(a[nn * nc + nn - 1 - nc]) + fabs(a[(nn - 1) * nc + nn - 2 - nc]);
              y = x = 0.75 * s;
              w = -0.4375 * s * s;
            }
            ++its;
            for(m = (nn - 2); m >= 1; m--) {
              z = a[m * nc + m - nc];
              r = x - z;
              s = y - z;
              p = (r * s - w) / a[(m + 1) * nc + m - nc] + a[m * nc + m + 1 - nc];
              q = a[(m + 1) * nc + m + 1 - nc] - z - r - s;
              r = a[(m + 2) * nc + m + 1 - nc];
              s = fabs(p) + fabs(q) + fabs(r);
              p /= s;
              q /= s;
              r /= s;
              if(m == 1) break;
              u = fabs(a[m * nc + m - 1 - nc]) * (fabs(q) + fabs(r));
              v = fabs(p) * (fabs(a[(m - 1) * nc + m - 1 - nc]) + fabs(z) + fabs(a[(m + 1) * nc + m + 1 - nc]));
              if((double)(u + v) == v) break;
            }
            for(i = m + 2; i <= nn; i++) {
              a[i * nc + i - 2 - nc] = 0.0;
              if(i != (m + 2)) a[i * nc + i - 3 - nc] = 0.0;
            }
            for(k = m; k <= nn - 1; k++) {
              if(k != m) {
                p = a[k * nc + k - 1 - nc];
                q = a[(k + 1) * nc + k - 1 - nc];
                r = 0.0;
                if(k != (nn - 1)) r = a[(k + 2) * nc + k - 1 - nc];
                if((x = fabs(p) + fabs(q) + fabs(r)) != 0.0) {
                  p /= x;
                  q /= x;
                  r /= x;
                }
              }
              if((s = SIGN(sqrt(p * p + q * q + r * r), p)) != 0.0) {
                if(k == m) {
                  if(l != m)
                    a[k * nc + k - 1 - nc] = -a[k * nc + k - 1 - nc];
                }
                else
                  a[k * nc + k - 1 - nc] = -s * x;
                p += s;
                x = p / s;
                y = q / s;
                z = r / s;
                q /= p;
                r /= p;
                for(j = k; j <= nn; j++) {
                  p = a[k * nc + j - nc] + q * a[(k + 1) * nc + j - nc];
                  if(k != (nn - 1)) {
                    p += r * a[(k + 2) * nc + j - nc];
                    a[(k + 2)*nc + j - nc] -= p * z;
                  }
                  a[(k + 1)*nc + j - nc] -= p * y;
                  a[k * nc + j - nc] -= p * x;
                }
                mmin = nn < k + 3 ? nn : k + 3;
                for(i = l; i <= mmin; i++) {
                  p = x * a[i * nc + k - nc] + y * a[i * nc + k + 1 - nc];
                  if(k != (nn - 1)) {
                    p += z * a[i * nc + k + 2 - nc];
                    a[i * nc + k + 2 - nc] -= p * r;
                  }
                  a[i * nc + k + 1 - nc] -= p * q;
                  a[i * nc + k - nc] -= p;
                }
              }
            }
          }
        }
      }
      while(l < nn - 1);
    }

    for(i = 0; i < nl; i++) {
      eigenvalue.ptr[i] = wr[i + 1];
      eigenvalue.ptr[i + nc] = wi[i + 1];
    }
    for(i = 0; i < nl * nc; i++)
      ptr[i] = back[i];
    delete [] wr;
    delete [] wi;
    delete [] a;
    delete [] back;
    return eigenvalue;
  }

  // Overloaded assignment operator
  // const return avoides: ( a1 = a2 ) = a3
  const Array &Array::operator=(const Array &right) {
    if(&right != this) {    // check for self-assignment

      // for arrays of different sizes, deallocate original
      // left side array, then allocate new left side array.
      if((nl != right.nl) || (nc != right.nc)) {
        delete [] ptr;         // reclaim space
        nl = right.nl;     // resize this object
        nc = right.nc;
        ptr = new double[ nl * nc ]; // create space for array copy
        //assert( ptr != 0);     // terminate if not allocated
      }

      for(int i = 0; i < nl * nc; i++)
        ptr[ i ] = right.ptr[ i ];  // copy array into object
    }

    return *this;   // enables x = y = z;
  }

  // Determine if two arrays are equal and
  // return true, otherwise return false.
  bool Array::operator==(const Array &right) const {
    if((nl != right.nl) || (nc != right.nc))
      return false;   // arrays of different sizes

    for(int i = 0; i < nl * nc; i++)
      if(ptr[ i ] != right.ptr[ i ])
        return false; // arrays are not equal

    return true;        // arrays are equal
  }

  // Overloaded input operator for class Array;
  // inpurts values for entire array.
  istream &operator>>(istream &input, Array &a) {
    for(int i = 0; i < ((a.nl) * (a.nc)); i++)
      input >> a.ptr[ i ];
    cout << "INPUT IS FINISHED!" << endl;

    return input;   // enables cin >> x >> y;
  }

  void Array::setArray(double *a) {
    for(int i = 0; i < nl * nc; i++)
      ptr[i] = a[i];
  }

  // Overloaded out put operator for class Array
  ostream &operator<<(ostream &output, const Array &a) {
    int i, j, c;
    c = a.nc;
    for(i = 0; i < a.nl; i++) {
      for(j = 0; j < a.nc; j++)
        output << setw(12) << a.ptr[i * c + j];
      output << endl;
    }
    return output;   // enables count << x << y;
  }

  //*********************************************************************
  //solve linear equations Ax=b using Gaussian elimination with full pivots, the
  //input is a combined matrix of A and b, the return is the vector solu(nx1)
  Array Array::solveleq() {
    int *js, l, k, i, j, is, p, q;
    double d, t;
    js = new int[nl];
    Array a(nl, nc - 1);
    Array b(nl, 1);
    Array solu(nl, 1);
    for(i = 0; i < nl; i++) {
      for(j = 0; j < nc - 1; j++) a.ptr[i * (nc - 1) + j] = ptr[i * nc + j];
      b.ptr[i] = ptr[i * nc + j];
    }


    l = 1;
    for(k = 0; k <= nl - 2; k++) {
      d = 0.0;
      for(i = k; i <= nl - 1; i++)
        for(j = k; j <= nl - 1; j++) {
          t = fabs(a.ptr[i * nl + j]);
          if(t > d) {
            d = t;
            js[k] = j;
            is = i;
          }
        }
      if(d + 1.0 == 1.0) l = 0;
      else {
        if(js[k] != k)
          for(i = 0; i <= nl - 1; i++) {
            p = i * nl + k;
            q = i * nl + js[k];
            t = a.ptr[p];
            a.ptr[p] = a.ptr[q];
            a.ptr[q] = t;
          }
        if(is != k) {
          for(j = k; j <= nl - 1; j++) {
            p = k * nl + j;
            q = is * nl + j;
            t = a.ptr[p];
            a.ptr[p] = a.ptr[q];
            a.ptr[q] = t;
          }
          t = b.ptr[k];
          b.ptr[k] = b.ptr[is];
          b.ptr[is] = t;
        }
      }
      if(l == 0) {
        delete[] js;
        cout << "fail" << endl;
        exit(1);
      }
      d = a.ptr[k * nl + k];
      for(j = k + 1; j <= nl - 1; j++) {
        p = k * nl + j;
        a.ptr[p] = a.ptr[p] / d;
      }
      b.ptr[k] = b.ptr[k] / d;
      for(i = k + 1; i <= nl - 1; i++) {
        for(j = k + 1; j <= nl - 1; j++) {
          p = i * nl + j;
          a.ptr[p] = a.ptr[p] - a.ptr[i * nl + k] * a.ptr[k * nl + j];
        }
        b.ptr[i] = b.ptr[i] - a.ptr[i * nl + k] * b.ptr[k];
      }
    }
    d = a.ptr[(nl - 1) * nl + nl - 1];
    if(fabs(d) + 1.0 == 1.0) {
      delete[] js;
      cout << "fail" << endl;
      return(0);
    }
    solu.ptr[nl - 1] = b.ptr[nl - 1] / d;
    for(i = nl - 2; i >= 0; i--) {
      t = 0.0;
      for(j = i + 1; j <= nl - 1; j++)
        t = t + a.ptr[i * nl + j] * solu.ptr[j];
      solu.ptr[i] = b.ptr[i] - t;
    }
    js[nl - 1] = nl - 1;
    for(k = nl - 1; k >= 0; k--)
      if(js[k] != k) {
        t = solu.ptr[k];
        solu.ptr[k] = solu.ptr[js[k]];
        solu.ptr[js[k]] = t;
      }
    delete[] js;
    return solu;
  }


  //*************************************************************************************************
  //Given a matrix a(nlxnc), this routine computes its single value decomposition. A=U*W*Vt.
  //The matrix Ua(nlxnc) replaces A(nlxnc) on output. The diagonal matrix of singular values Wa
  //is a vector w(ncx1). The matrix Va(ncxnc) (not the transpose Vt) is as v(ncxnc). These three
  //parts are entered Ua, Wa and Va. See Ref. numerical recipes in C. William H. Press
  void Array::svd(Array &Ua, Array &Wa, Array &Va) {
    double *a = new double[nl * nc + 1];
    double *w = new double[nc + 1];
    double *v = new double[nc * nc + 1];
    double *rv1 = new double[nc + 1];
    int flag, i, its, j, jj, k, l, nm;
    double anorm, c, f, g, h, s, scale, x, y, z;
    rv1[0] = 1.0;
    a[0] = 1.0;
    w[0] = 1.0;
    v[0] = 1.0;
    for(i = 0; i < nl * nc; i++)
      a[i + 1] = ptr[i];
    g = scale = anorm = 0.0;
    for(i = 1; i <= nc; i++) {
      l = i + 1;
      rv1[i] = scale * g;
      g = s = scale = 0.0;
      if(i <= nl) {
        for(k = i; k <= nl; k++) scale += fabs(a[k * nc + i - nc]);
        if(scale) {
          for(k = i; k <= nl; k++) {
            a[k * nc + i - nc] /= scale;
            s += a[k * nc + i - nc] * a[k * nc + i - nc];
          }
          f = a[i * nc + i - nc];
          g = -SIGN(sqrt(s), f);
          h = f * g - s;
          a[i * nc + i - nc] = f - g;
          for(j = l; j <= nc; j++) {
            for(s = 0.0, k = i; k <= nl; k++) s += a[k * nc + i - nc] * a[k * nc + j - nc];
            f = s / h;
            for(k = i; k <= nl; k++) a[k * nc + j - nc] += f * a[k * nc + i - nc];
          }
          for(k = i; k <= nl; k++) a[k * nc + i - nc] *= scale;
        }
      }
      w[i] = scale * g;
      g = s = scale = 0.0;
      if(i <= nl && i != nc) {
        for(k = l; k <= nc; k++) scale += fabs(a[i * nc + k - nc]);
        if(scale) {
          for(k = l; k <= nc; k++) {
            a[i * nc + k - nc] /= scale;
            s += a[i * nc + k - nc] * a[i * nc + k - nc];
          }
          f = a[i * nc + l - nc];
          g = -SIGN(sqrt(s), f);
          h = f * g - s;
          a[i * nc + l - nc] = f - g;
          for(k = l; k <= nc; k++) rv1[k] = a[i * nc + k - nc] / h;
          for(j = l; j <= nl; j++) {
            for(s = 0.0, k = l; k <= nc; k++) s += a[j * nc + k - nc] * a[i * nc + k - nc];
            for(k = l; k <= nc; k++) a[j * nc + k - nc] += s * rv1[k];
          }
          for(k = l; k <= nc; k++) a[i * nc + k - nc] *= scale;
        }
      }
      anorm = FMAX(anorm, (fabs(w[i]) + fabs(rv1[i])));
    }
    for(i = nc; i >= 1; i--) {
      if(i < nc) {
        if(g) {
          for(j = l; j <= nc; j++)
            v[j * nc + i - nc] = (a[i * nc + j - nc] / a[i * nc + l - nc]) / g;
          for(j = l; j <= nc; j++) {
            for(s = 0.0, k = l; k <= nc; k++) s += a[i * nc + k - nc] * v[k * nc + j - nc];
            for(k = l; k <= nc; k++) v[k * nc + j - nc] += s * v[k * nc + i - nc];
          }
        }
        for(j = l; j <= nc; j++) v[i * nc + j - nc] = v[j * nc + i - nc] = 0.0;
      }
      v[i * nc + i - nc] = 1.0;
      g = rv1[i];
      l = i;
    }
    for(i = IMIN(nl, nc); i >= 1; i--) {
      l = i + 1;
      g = w[i];
      for(j = l; j <= nc; j++) a[i * nc + j - nc] = 0.0;
      if(g) {
        g = 1.0 / g;
        for(j = l; j <= nc; j++) {
          for(s = 0.0, k = l; k <= nl; k++) s += a[k * nc + i - nc] * a[k * nc + j - nc];
          f = (s / a[i * nc + i - nc]) * g;
          for(k = i; k <= nl; k++) a[k * nc + j - nc] += f * a[k * nc + i - nc];
        }
        for(j = i; j <= nl; j++) a[j * nc + i - nc] *= g;
      }
      else for(j = i; j <= nl; j++) a[j * nc + i - nc] = 0.0;
      ++a[i * nc + i - nc];
    }
    for(k = nc; k >= 1; k--) {
      for(its = 1; its <= 50; its++) {
        flag = 1;
        for(l = k; l >= 1; l--) {
          nm = l - 1;
          if((double)(fabs(rv1[l]) + anorm) == anorm) {
            flag = 0;
            break;
          }
          if((double)(fabs(w[nm]) + anorm) == anorm) break;
        }
        if(flag) {
          c = 0.0;
          s = 1.0;
          for(i = l; i <= k; i++) {
            f = s * rv1[i];
            rv1[i] = c * rv1[i];
            if((double)(fabs(f) + anorm) == anorm) break;
            g = w[i];
            h = useful(f, g);
            w[i] = h;
            h = 1.0 / h;
            c = g * h;
            s = -f * h;
            for(j = 1; j <= nl; j++) {
              y = a[j * nc + nm - nc];
              z = a[j * nc + i - nc];
              a[j * nc + nm - nc] = y * c + z * s;
              a[j * nc + i - nc] = z * c - y * s;
            }
          }
        }
        z = w[k];
        if(l == k) {
          if(z < 0.0) {
            w[k] = -z;
            for(j = 1; j <= nc; j++) v[j * nc + k - nc] = -v[j * nc + k - nc];
          }
          break;
        }
        if(its == 50) {
          cout << "no convergence in 30 svd iterations";
          exit(1);
        }
        x = w[l];
        nm = k - 1;
        y = w[nm];
        g = rv1[nm];
        h = rv1[k];
        f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
        g = useful(f, 1.0);
        f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g, f))) - h)) / x;
        c = s = 1.0;
        for(j = l; j <= nm; j++) {
          i = j + 1;
          g = rv1[i];
          y = w[i];
          h = s * g;
          g = c * g;
          z = useful(f, h);
          rv1[j] = z;
          c = f / z;
          s = h / z;
          f = x * c + g * s;
          g = g * c - x * s;
          h = y * s;
          y *= c;
          for(jj = 1; jj <= nc; jj++) {
            x = v[jj * nc + j - nc];
            z = v[jj * nc + i - nc];
            v[jj * nc + j - nc] = x * c + z * s;
            v[jj * nc + i - nc] = z * c - x * s;
          }
          z = useful(f, h);
          w[j] = z;
          if(z) {
            z = 1.0 / z;
            c = f * z;
            s = h * z;
          }
          f = c * g + s * y;
          x = c * y - s * g;
          for(jj = 1; jj <= nl; jj++) {
            y = a[jj * nc + j - nc];
            z = a[jj * nc + i - nc];
            a[jj * nc + j - nc] = y * c + z * s;
            a[jj * nc + i - nc] = z * c - y * s;
          }
        }
        rv1[l] = 0.0;
        rv1[k] = f;
        w[k] = x;
      }
    }
    for(i = 0; i < nl * nc; i++) Ua.ptr[i] = a[i + 1];
    for(i = 0; i < nc; i++) Wa.ptr[i] = w[i + 1];
    for(i = 0; i < nc * nc; i++) Va.ptr[i] = v[i + 1];
    delete [] a;
    delete [] w;
    delete [] v;
    delete [] rv1;
  }
  //*******************************************************************************************************
  //Use svd to solve AX=B for vector X, the input is the combined matrix of AB(nlxnc), return vector
  //solu(nlx1).A is specified by its svd components u(nlx(nc-1)), w[1...nc-1],v((nc-1)x(nc-1)) as
  //returned by svd(). solu((nc-1)x1) is the output solution vector.
  Array Array::svdeq() {
    Array a(nl, nc - 1);
    Array solu(nc - 1, 1);
    Array Ua(nl, nc - 1), Wa(nc - 1, 1), Va(nc - 1, nc - 1);
    double *b = new double [nl];
    double *tmp = new double [nc - 1];
    int jj, j, i;
    double s;

    for(i = 0; i < nl; i++) {
      for(j = 0; j < nc - 1; j++) a.ptr[i * (nc - 1) + j] = ptr[i * nc + j];
      b[i] = ptr[i * nc + j];
    }
    a.svd(Ua, Wa, Va);
    for(j = 0; j < nc - 1; j++) {
      s = 0.0;
      if(Wa.ptr[j]) {
        for(i = 0; i < nl; i++) s += Ua.ptr[i * (nc - 1) + j] * b[i];
        s /= Wa.ptr[j];
      }
      tmp[j] = s;
    }
    for(j = 0; j < nc - 1; j++) {
      s = 0.0;
      for(jj = 0; jj < nc - 1; jj++) s += Va.ptr[j * (nc - 1) + jj] * tmp[jj];
      solu.ptr[j] = s;
    }
    delete [] tmp;
    delete [] b;

    return solu;
  }
  //**************************************************************************
  //multidimensional least square fit using SVD. nc coefficients and nl functions.
  //nl>=nc. The input is vector Y(nl,1) and covariance matrix covar(nc,nc), if you have a
  //coefficient matrix coef(nl,nc), then coef.svdls(Y, covar) return the fitted
  //coefficient vectro solu(nc,1) and also get the covariance matrix covar(nc,nc). Assume standard
  //deviation of measured points is 1.0
  //cvn is the number of data used to calculate cv score, and Efitw stores the fitted data using
  //fitted coeficients solu
#define TOL 2.220446049250313e-16
  // TOL is the machine precision
  Array Array::svdls(Array &Y, Array &covar, double &rms, double &cv, int &cvn) {
    int i, j, jj, k, Stru_wcv, iii;
    double wmax, thresh, sum, s, a;
    double *tmp = new double [nc];
    double *wti = new double [nc];
    Array U(nl, nc), W(nc, 1), V(nc, nc), solu(nc, 1), U2(nc, nc), W2(nc, 1), V2(nc, nc);
    ((this->transpose()) * (*this)).svd(U2, W2, V2);
    for(iii = 0; iii < nc; iii++)
      if(W2.ptr[iii] < 1.0e-8) {
        //cout<<"when doing svdls(),singular matrix";
        rms = 1.0e20;
        cv = 1.0e20;
        delete [] tmp;
        delete [] wti;
        return solu;
      }
    this->svd(U, W, V);
    wmax = 0.0;
    for(j = 0; j < nc; j++)
      if(W.ptr[j] > wmax) wmax = W.ptr[j];
    thresh = nl * TOL * wmax;
    for(j = 0; j < nc; j++)
      if(W.ptr[j] < thresh) W.ptr[j] = 0.0;

    for(j = 0; j < nc; j++) {
      s = 0.0;
      if(W.ptr[j]) {
        for(i = 0; i < nl; i++) s += U.ptr[i * nc + j] * Y.ptr[i];
        s /= W.ptr[j];
      }
      tmp[j] = s;
    }
    for(j = 0; j < nc; j++) {
      s = 0.0;
      for(jj = 0; jj < nc; jj++) s += V.ptr[j * nc + jj] * tmp[jj];
      solu.ptr[j] = s;
    }                       //calculate fitted predictors stored in solu[], e.g. ECIs

    for(i = 0; i < nc; i++) {
      wti[i] = 0.0;
      if(W.ptr[i]) wti[i] = 1.0 / ((W.ptr[i]) * (W.ptr[i]));
    }
    for(i = 0; i < nc; i++) {
      for(j = 0; j <= i; j++) {
        for(sum = 0.0, k = 0; k < nc; k++) sum += V.ptr[i * nc + k] * V.ptr[j * nc + k] * wti[k];
        covar.ptr[j * nc + i] = covar.ptr[i * nc + j] = sum;
      }
    }		               //calculate covariance matrix covar(nc,nc) of fitted predictiors


    double *Efitw = new double [nl]; // store calculated E from ECIs, this E include structure weight factors
    cv = 0.0;                        //cross validation score
    rms = 0.0;                       //root mean square
    cvn = nl;
    for(i = 0; i < nl; i++) Efitw[i] = 0.0; // store fitted energy of structure i, include structure weight factor.
    for(i = 1; i <= nl; i++) {
      for(j = 1; j <= nc; j++) Efitw[i - 1] += solu.elem(j, 1) * elem(i, j); // calculate fitted energy of structure i, including structure weight factor.
      //cout<<((Y.elem(i,1)-Efitw[i-1]))*((Y.elem(i,1)-Efitw[i-1]))<<endl;
      rms += ((Y.elem(i, 1) - Efitw[i - 1])) * ((Y.elem(i, 1) - Efitw[i - 1])); // rms also includes the structure weight factor
      a = ((row(i)) * (((this->transpose()) * (*this)).inverSVD()) * ((row(i)).transpose())).elem(1, 1);
      if((a - 1.0 < 1.0e-11) && (a - 1.0 > -1.0e-11)) {
        cvn--;  // if a=1.0, then actual Ei=fitted Ei, don't include this
        continue;
      }
      // structure i in WCV calculation. Stru_wcv store the number
      // of structures included in wcv calculation.
      else cv += ((Y.elem(i, 1) - Efitw[i - 1]) / (1.0 - a)) * ((Y.elem(i, 1) - Efitw[i - 1]) / (1.0 - a));
    }
    cv = sqrt(cv / double(cvn));
    rms = sqrt(rms / double(nl));

    delete [] Efitw;
    delete [] tmp;
    delete [] wti;
    return solu;
  }
  //********************************************************************************************************************
  //multidimensional least square fit using SVD. nc coefficients and nl functions.
  //nl>=nc. The input is vector Y(nl,1) and covariance matrix covar(nc,nc), if you have a
  //coefficient matrix coef(nl,nc), then coef.svdls(Y, covar) return the fitted
  //coefficient vector solu(nc,1) and also get the covariance matrix covar(nc,nc). Assume standard
  //deviation of measured points is 1.0
  // Efitw stores the fitted data using fitted coeficients solu
  //calculate cv using Monte Carlo cv method, see ref. "Linear Model Selection by Cross-Validation"
  //Jun Shao, J. of the American Statistical Association, Jun 1993, 88, 422, pp486
#define TOL 2.220446049250313e-16
  // TOL is the machine precision
  //Nt is the number of data points used to fit
  //b is the number of subsets of data points that have size Nv=N-Nc
  //effect_b is the number of effective divisions
  Array Array::svdls_MCCV(Array &Y, Array &covar, double &rms, double &cv, int Nt, int b, int &effect_b) {
    //cout<<"svdls_MCCV is called"<<endl;
    int i, j, jj, k;
    double wmax, thresh, sum, s;
    double *tmp = new double [nc];
    double *wti = new double [nc];
    Array U(nl, nc), W(nc, 1), V(nc, nc), solu(nc, 1);
    this->svd(U, W, V);
    wmax = 0.0;
    for(j = 0; j < nc; j++)
      if(W.ptr[j] > wmax) wmax = W.ptr[j];
    thresh = nl * TOL * wmax;
    for(j = 0; j < nc; j++)
      if(W.ptr[j] < thresh) W.ptr[j] = 0.0;

    for(j = 0; j < nc; j++) {
      s = 0.0;
      if(W.ptr[j]) {
        for(i = 0; i < nl; i++) s += U.ptr[i * nc + j] * Y.ptr[i];
        s /= W.ptr[j];
      }
      tmp[j] = s;
    }
    for(j = 0; j < nc; j++) {
      s = 0.0;
      for(jj = 0; jj < nc; jj++) s += V.ptr[j * nc + jj] * tmp[jj];
      solu.ptr[j] = s;
    }               //calculate fitted predictors, e.g. ECIs, stored in solu[nc], using the whole data set

    for(i = 0; i < nc; i++) {
      wti[i] = 0.0;
      if(W.ptr[i]) wti[i] = 1.0 / ((W.ptr[i]) * (W.ptr[i]));
    }
    for(i = 0; i < nc; i++) {
      for(j = 0; j <= i; j++) {
        for(sum = 0.0, k = 0; k < nc; k++) sum += V.ptr[i * nc + k] * V.ptr[j * nc + k] * wti[k];
        covar.ptr[j * nc + i] = covar.ptr[i * nc + j] = sum;
      }
    }		      //calculate covariance matrix of predictors, stored in covar(nc,nc)


    double *Efitw = new double [nl]; // store calculated E from ECIs, this E include structure weight factors
    rms = 0.0;                       // root mean square
    for(i = 0; i < nl; i++) Efitw[i] = 0.0; // store fitted energy of structure i, include structure weight factor.
    for(i = 1; i <= nl; i++) {
      for(j = 1; j <= nc; j++) Efitw[i - 1] += solu.elem(j, 1) * elem(i, j); // calculate fitted energy of structure i, including structure weight factor.
      rms += ((Y.elem(i, 1) - Efitw[i - 1])) * ((Y.elem(i, 1) - Efitw[i - 1])); // rms also includes the structure weight factor
    }

    rms = sqrt(rms / double(nl));   //rms is based on the whole data set
    //cout<<"prepare to call MCCV()"<<endl;
    cv = MCCV(Y, Nt, b, effect_b);            //calculate cv using Monte Carlo Leave Many Out algorithm
    //cv = MCCV_formula(Y, Nt, b, solu);
    delete [] Efitw;
    delete [] tmp;
    delete [] wti;
    return solu;
  }
  ////////////////////////////////////////////////////////////////////////////////////
  //Nt is the number of data points used to fit
  //b is the number of subsets of data points that have size Nv=N-Nc, i.e. number of division
  //Y is the actual results, e.g. energies, of the whole data set
  //solu is the predicto values (e.g. ECIs) fitted based on the whole data set
  //MCCV_formula() return cv score using Monte Carlo Leave Many(Nv) Out and formula,
  //see ref. "Linear Model Selection by Cross-Validation", Jun Shao, J. of the American Statistical
  //Association, Jun 1993, 88, 422, pp486, formula (3.1)
  //Note: When running, should set random number seed by using "srand(static_cast<unsigned>(time(0)));"
  double Array::MCCV_formula(Array &Y, int Nt, int b, Array &solu) { //using formula is slower than using definetion
    cout << "MCCV_formula() is called" << endl;
    int i, j, k, m, fit_num, b_num, effect_div_num;
    double cv1;
    int Nv = nl - Nt;           //Nv is the number of data points in testing subset.
    Array YNt(Nt, 1), YNv(Nv, 1); //YNt stores the actual results in fitting subset, YNv stores the actual results in testing subset
    Array ANt(Nt, nc), ANv(Nv, nc); //ANt stores the correlation matrix in fitting subset, ANv stores the correlation matrix in testing subset
    Array solu2(nc, 1);            // stores the fitting coefficient, e.g. ECIs
    Array UNt(Nt, nc), WNt(nc, 1), VNt(nc, nc);
    Array X(nl, nc);
    for(i = 0; i < nl * nc; i++) X.ptr[i] = ptr[i];
    int *fit_array = new int [Nt];  //stores the index of data in fitting subset
    int *non_fit_array =  new int [Nv]; //stores the index of data in testing subset
    double *fit_E =  new double [Nv];   //stores the fitted energies in testing subset
    int fit_index;
    bool exist, in_fit, linear;
    cv1 = 0;
    effect_div_num = 0;
    for(b_num = 0; b_num < b; b_num++) {   //divide the whole data set b times
      // cout<<"Divide time "<<b_num<<endl;;
      for(i = 0; i < Nt; i++) fit_array[i] = 0;
      for(i = 0; i < Nv; i++) non_fit_array[i] = 0;
      for(i = 0; i < Nv; i++) fit_E[i] = 0;
      k = 0;                              //index of non_fit_array[Nv]
      linear =  false;
      for(fit_num = 0; fit_num < Nt; fit_num++) { //randomly select Nt data points to form fitting subset
        do {
          exist = false;
          fit_index = (int)((nl) * ((double)rand() / ((double)(RAND_MAX) + (double)(1)))) + 1; //randomly select a data point [1, nl]
          for(i = 0; i < Nt; i++)
            if(fit_array[i] == fit_index) {
              exist = true;
              break;
            }
        }
        while(exist == true);
        fit_array[fit_num] = fit_index;
      }
      for(i = 1; i <= nl; i++) {
        in_fit = false;
        for(j = 0; j < Nt; j++) if(i == fit_array[j]) {
            in_fit = true;
            break;
          }
        if(in_fit == false) {
          non_fit_array[k] = i;
          k++;
        }
      }
      for(i = 0; i < Nt; i++) { //store the actual results and correlation matrix in fitting subset in YNt and ANt
        YNt.ptr[i] = Y.ptr[fit_array[i] - 1];
        for(j = 0; j < nc; j++) ANt.ptr[i * nc + j] = ptr[(fit_array[i] - 1) * nc + j];
      }
      for(i = 0; i < Nv; i++) { //store the actual results and correlation matrix in testing subset in YNv and ANv
        YNv.ptr[i] = Y.ptr[non_fit_array[i] - 1];
        for(j = 0; j < nc; j++) ANv.ptr[i * nc + j] =  ptr[(non_fit_array[i] - 1) * nc + j];
      }
      ANt.svd(UNt, WNt, VNt);
      for(m = 0; m < nc; m++)
        if(WNt.ptr[m] < 1.0e-8) {
          //  cout<<"divide "<<b_num<<" is invalid"<<endl;
          linear = true;
          break;
        }
      if(linear) continue;    //this division induce linear dependence, is non-effective.
      else effect_div_num++;

      Array I_nv(Nv, Nv), Q_a_s(Nv, Nv);
      I_nv.identity();
      Q_a_s = ANv * ((X.transpose() * X).inverSVD()) * ANv.transpose();
      double Mod;
      Mod = ((I_nv - Q_a_s).inverSVD() * (YNv - ANv * solu)).modu();
      cv1 += Mod * Mod;

    }

    if(effect_div_num <= nl) cv1 = 1.0e20;           //if there are so many ineffective division, the cv is invalid.
    //we should insure that the effective division time is larger
    //than the number of the whole data points
    else cv1 = sqrt(cv1 / double(Nv * effect_div_num));

    //cout<<"The number of effective division is: "<<effect_div_num<<endl;
    //cout<<"nl is "<<nl<<endl;

    delete [] fit_array;
    delete [] non_fit_array;
    delete [] fit_E;
    return cv1;
  }

  /////////////////////////////////////////////////////////////////////////////////////////
  //Nt is the number of data points used to fit
  //b is the number of subsets of data points that have size Nv=N-Nc, i.e. number of division
  //Y is the actual results, e.g. energies
  //MCCV() return cv score using Monte Carlo Leave Many(Nv) Out
  //Note: When running, should set random number seed by using "srand(static_cast<unsigned>(time(0)));"
  double Array::MCCV(Array &Y, int Nt, int b, int &effect_div_num) {
    //cout<<"MCCV() is called"<<endl;
    int i, j, k, m, fit_num, b_num;
    double cv1;
    int Nv = nl - Nt;           //Nv is the number of data points in testing subset.
    Array YNt(Nt, 1), YNv(Nv, 1); //YNt stores the actual results in fitting subset, YNv stores the actual results in testing subset
    Array ANt(Nt, nc), ANv(Nv, nc); //ANt stores the correlation matrix in fitting subset, ANv stores the correlation matrix in testing subset
    Array solu2(nc, 1);            // stores the fitting coefficient, e.g. ECIs
    Array UNt(Nt, nc), WNt(nc, 1), VNt(nc, nc);
    int *fit_array = new int [Nt];  //stores the index of data in fitting subset
    int *non_fit_array =  new int [Nv]; //stores the index of data in testing subset
    double *fit_E =  new double [Nv];   //stores the fitted energies in testing subset
    int fit_index;
    bool exist, in_fit, linear;
    cv1 = 0;
    effect_div_num = 0;
    for(b_num = 0; b_num < b; b_num++) {   //divide the whole data set b times
      // cout<<"Divide time "<<b_num<<endl;;
      for(i = 0; i < Nt; i++) fit_array[i] = 0;
      for(i = 0; i < Nv; i++) non_fit_array[i] = 0;
      for(i = 0; i < Nv; i++) fit_E[i] = 0;
      k = 0;                              //index of non_fit_array[Nv]
      linear =  false;
      for(fit_num = 0; fit_num < Nt; fit_num++) { //randomly select Nt data points to form fitting subset
        do {
          exist = false;
          fit_index = (int)((nl) * ((double)rand() / ((double)(RAND_MAX) + (double)(1)))) + 1; //randomly select a data point [1, nl]
          for(i = 0; i < Nt; i++)
            if(fit_array[i] == fit_index) {
              exist = true;
              break;
            }
        }
        while(exist == true);
        fit_array[fit_num] = fit_index;
      }
      for(i = 1; i <= nl; i++) {
        in_fit = false;
        for(j = 0; j < Nt; j++) if(i == fit_array[j]) {
            in_fit = true;
            break;
          }
        if(in_fit == false) {
          non_fit_array[k] = i;
          k++;
        }
      }
      for(i = 0; i < Nt; i++) { //store the actual results and correlation matrix in fitting subset in YNt and ANt
        YNt.ptr[i] = Y.ptr[fit_array[i] - 1];
        for(j = 0; j < nc; j++) ANt.ptr[i * nc + j] = ptr[(fit_array[i] - 1) * nc + j];
      }
      for(i = 0; i < Nv; i++) { //store the actual results and correlation matrix in testing subset in YNv and ANv
        YNv.ptr[i] = Y.ptr[non_fit_array[i] - 1];
        for(j = 0; j < nc; j++) ANv.ptr[i * nc + j] =  ptr[(non_fit_array[i] - 1) * nc + j];
      }
      ANt.svd(UNt, WNt, VNt);
      for(m = 0; m < nc; m++)
        if(WNt.ptr[m] < 1.0e-8) {
          //  cout<<"divide "<<b_num<<" is invalid"<<endl;
          linear = true;
          break;
        }
      if(linear) continue;    //this division induce linear dependence, is non-effective.
      else effect_div_num++;

      solu2 = ANt.svdls_core(YNt);
      for(i = 0; i < Nv; i++) //calculate fitted energies for testing subset
        for(j = 0; j < nc; j++)
          fit_E[i] += solu2.ptr[j] * ANv.ptr[i * nc + j];
      //   cout<<"fit_E "<<i<<" =  "<<fit_E[i]<<" , "<<"actual E = "<<YNv.ptr[i]<<endl;
      for(i = 0; i < Nv; i++) cv1 = cv1 + (fit_E[i] - YNv.ptr[i]) * (fit_E[i] - YNv.ptr[i]);	// add cv1
    }

    if(effect_div_num <= (int(0.75 * nl))) cv1 = 1.0e20;         //if there are so many ineffective division, the cv is invalid.
    //we should insure that the effective division time is larger
    //than the number of the whole data points
    else cv1 = sqrt(cv1 / double(Nv * effect_div_num));

    //cout<<"The number of effective division is: "<<effect_div_num<<endl;
    //cout<<"nl is "<<nl<<endl;

    delete [] fit_array;
    delete [] non_fit_array;
    delete [] fit_E;
    return cv1;
  }


  ////////////////////////////////////////////////////////////////////////////////////
  //Array svdls_core(Array &)
  //multidimensional least square fit using SVD. nc coefficients and nl functions.
  //nl>=nc. The input is vector Y(nl,1). if you have a coefficient matrix coef(nl,nc),
  //then coef.svdls(Y1) return the fitted coefficient vector solu1(nc,1)
  Array Array::svdls_core(Array &Y1) {
    int i, j, jj;
    double wmax, thresh, s;
    double *tmp1 = new double [nc];
    Array U1(nl, nc), W1(nc, 1), V1(nc, nc), solu1(nc, 1);
    this->svd(U1, W1, V1);
    wmax = 0.0;
    for(j = 0; j < nc; j++)
      if(W1.ptr[j] > wmax) wmax = W1.ptr[j];
    thresh = nl * TOL * wmax;
    for(j = 0; j < nc; j++)
      if(W1.ptr[j] < thresh) W1.ptr[j] = 0.0;

    for(j = 0; j < nc; j++) {
      s = 0.0;
      if(W1.ptr[j]) {
        for(i = 0; i < nl; i++) s += U1.ptr[i * nc + j] * Y1.ptr[i];
        s /= W1.ptr[j];
      }
      tmp1[j] = s;
    }
    for(j = 0; j < nc; j++) {
      s = 0.0;
      for(jj = 0; jj < nc; jj++) s += V1.ptr[j * nc + jj] * tmp1[jj];
      solu1.ptr[j] = s;
    }
    delete [] tmp1;
    return solu1;
  }
  //****************************************************************************************************************
  //cv(Array &Y, int &num_cv) calculate cv by using Leave One Out cv definition.
  //num_cv is the number of structures used to calculate cv
  double Array::cv(Array &Y, int &num_cv) { // Y stores the actual energies
    int kk = 0;
    int ii, i, j, p, cvn_no_i, m, screen;
    num_cv = 0;
    bool linear;
    Array U_no_i(nl - 1, nc), W_no_i(nc, 1), V_no_i(nc, nc);
    Array A_no_i(nl - 1, nc);
    Array Y_no_i(nl - 1, 1);
    Array covar_no_i(nc, nc);
    Array solu_no_i(nc, 1);
    double rms_no_i, cv_no_i, Ei, cv_score;
    double Efit = 0.0;
    double sum = 0.0;
    cv_score = 0.0;
    double *X = new double [nc];
    Array ECI_no_i(nc, 1);
    for(ii = 0; ii < nl; ii++) {
      linear = false;
      kk = 0;
      Efit = 0.0;
      for(i = 0; i < nl; i++) {
        if(i == ii) {
          Ei = Y.ptr[i];
          for(p = 0; p < nc; p++)
            X[p] = ptr[i * nc + p];
          continue;
        }
        for(j = 0; j < nc; j++)
          A_no_i.ptr[kk * nc + j] = ptr[i * nc + j];
        Y_no_i.ptr[kk] = Y.ptr[i];
        kk++;
      }

      A_no_i.svd(U_no_i, W_no_i, V_no_i);
      for(m = 1; m <= nc; m++)
        if(W_no_i.elem(m, 1) < 1.0e-8) {
          linear = true;
          continue;
        }
      if(linear) continue;
      ECI_no_i = A_no_i.svdls(Y_no_i, covar_no_i, rms_no_i, cv_no_i, cvn_no_i);
      for(j = 0; j < nc; j++) Efit += ECI_no_i.ptr[j] * X[j];
      sum = sum + (Ei - Efit) * (Ei - Efit);
      num_cv++;
    }
    cv_score = sqrt(sum / double(num_cv));
    delete [] X;
    return cv_score;
  }

  //***************************************************************************************************
  double Array::useful(double a, double b) { //useful function used by svd()
    double absa, absb;
    absa = fabs(a);
    absb = fabs(b);
    if(absa > absb) return absa * sqrt(1.0 + (absb / absa) * (absb / absa));
    else return (absb == 0.0 ? 0.0 : absb * sqrt(1.0 + (absa / absb) * (absa / absb)));
  }

  //*************************************************************************************************
  //calculate inverse by using SVD method
  Array Array::inverSVD() {
    Array u(nl, nc), ut(nc, nl), w(nc, 1), v(nc, nc), wm(nc, nc);
    double wmax, thresh;
    this->svd(u, w, v);
    wmax = 0.0;
    for(int j = 0; j < nc; j++)
      if(w.ptr[j] > wmax) wmax = w.ptr[j];
    thresh = nl * TOL * wmax;
    for(int j = 0; j < nc; j++) {
      if(w.ptr[j] < thresh) {
        cout << "This matrix is singular";
        exit(1) ;
      }
      else  wm.ptr[j * nc + j] = 1.0 / w.ptr[j];
    }
    ut = u.transpose();
    return (v * wm) * ut;
  }

  //************************************
  // get the ith row of matrix, i is from 1.
  Array Array::row(int i) {
    Array row(1, nc);
    if(i > nl) {
      cout << "Error! The number of row is beyond the row range!";
      exit(1);
    }
    for(int j = 0; j < nc; j++) {
      row.ptr[j] = ptr[(i - 1) * nc + j];
    }
    return row;
  }
  //get the jth column of matrix, j is from 1.
  Array Array::col(int j) {
    Array column(nl, 1);
    if(j > nc) {
      cout << "Error! The number of column is beyond the column range!";
      exit(1);
    }
    for(int i = 0; i < nl; i++) {
      column.ptr[i] = ptr[i * nc + j - 1];
    }
    return column;
  }
  //******************* Find hull of a set of points
  // the input points are stored in a 3-column Array, the 1st column is the point label
  // the 2nd column is the x coordinate and the 3rd column is the y coordinate. This
  // code use gift wrapping algorithm. The order of the points in the input Array will
  // be changed. The code returns a matrix containing label and coordinates of hull points
  Array Array::hull2d(int &nh) { //nh is the number of hull points. This code may miss the points on the hull edge.
    Array areaA(3, 3);             //store the matrix for determinant calculation
    int i, j, k, m, x, y, z, ppp, pp, aa;
    bool h = true;
    stack <Array> hull_stack;
    queue <Array> edge;
    int num_edge = 0;
    Array edgeA(1, 2);
    Array transA(1, 3);
    queue <Array> remain_queue;
    int remain_size;
    areaA.ptr[6] = areaA.ptr[7] = areaA.ptr[8] = 1.0;
    k = 0;                          // k+1 is the number of hull points found
    for(i = 1; i < nl; i++) //find the starting point with lowest y, and put it into 1st row
      if(ptr[i * 3 + 2] < ptr[2]) exrow(1, i + 1);
    //for(jjj=0; jjj<3; jjj++)	transA.ptr[jjj]=ptr[jjj];
    //hull_stack.push(transA);
    for(m = 0; m < nl - 1; m++) {
      if(h == false) break;
      areaA.ptr[0] = ptr[m * 3 + 1]; //the 1st column of matrix area stores the determined hull point A
      areaA.ptr[3] = ptr[m * 3 + 2];
      edgeA.ptr[0] = ptr[m * 3 + 0]; //the label of hull point A is stored in edgeA[0]
      for(i = m + 1; i < nl; i++) {
        h = true;
        areaA.ptr[1] = ptr[i * 3 + 1]; //the 2nd column of matrix area stores the hull point candidate B
        areaA.ptr[4] = ptr[i * 3 + 2];
        for(j = 0; j < nl; j++) {
          if((j == i) || (j == m)) continue;
          areaA.ptr[2] = ptr[j * 3 + 1]; //the 3rd column of matrix area stores the third point C
          areaA.ptr[5] = ptr[j * 3 + 2];
          if(areaA.det() < 0) {
            h = false;  // if for all other points except A and B, the area is positive or zero
            break;
          }
        }                                     // that means all points C is to the left of line AB, then keep B as a
        // hull point
        if(h == true) {
          k++;    // the (k+1)th hull point is found
          exrow(i + 1, k + 1);
          edgeA.ptr[1] = ptr[k * 3 + 0]; // thte label of hull point B is stored in edgaA[1]
          edge.push(edgeA);
          num_edge++;
          break;
        }
      }
    }

    aa = k + 1;
    for(i = 0; i < aa; i++) { //push hull points into hull stack
      for(j = 0; j < 3; j++)
        transA.ptr[j] = ptr[i * 3 + j];
      hull_stack.push(transA);
    }
    for(i = aa; i < nl; i++) { //push remaining points into remain queue
      for(j = 0; j < 3; j++)
        transA.ptr[j] = ptr[i * 3 + j];
      remain_queue.push(transA);
    }

    for(ppp = 0; ppp < num_edge; ppp++) {
      if(remain_queue.empty()) break;
      remain_size = remain_queue.size();
      edgeA = edge.front();  // get a edge from the facet queue
      edge.push(edgeA);
      for(x = 0; x < 2; x++)
        for(y = 0; y < nl; y++)
          if(ptr[y * 3 + 0] == edgeA.ptr[x]) {
            for(z = 0; z < 2; z++)
              areaA.ptr[z * 3 + x] = ptr[y * 3 + z + 1];
            break;
          }                         //get one edge and put the corresponding points in the deter matrix
      for(pp = 0; pp < remain_size; pp++) {
        transA = remain_queue.front(); //get one point from remain queue.
        for(z = 0; z < 2; z++)                     // put the coordinates in the last column of deter matrix
          areaA.ptr[z * 3 + 2] = transA.ptr[z + 1];
        if(areaA.det() != 0) {          // the point is not on the edge
          remain_queue.push(transA);
          remain_queue.pop();
          continue;				// get the next point from remain queue
        }
        else {                         // the point is on the facet
          hull_stack.push(transA);          // push the point in the hull stack
          k++;
          remain_queue.pop();              // delete the point from the remain queue
        }
      }
      edge.pop();   //delete the old edge
    }

    Array hullA(k + 1, 3); // the 1st column is the label of the hull points, the 2nd and 3rd are coordinates of the hull points
    for(i = 0; i < (k + 1); i++) {
      transA = hull_stack.top();
      for(j = 0; j < 3; j++)
        hullA.ptr[i * 3 + j] = transA.ptr[j];
      hull_stack.pop();
    }

    nh = k + 1;

    return hullA;
  }
  //************************************************************************************************


  void Array::exrow(int ri, int rj) { // exchange the ri row and rj row, ri and rj are from 1
    int j;
    double *tempe = new double [nc];
    for(j = 0; j < nc; j++)
      tempe[j] = ptr[(ri - 1) * nc + j];
    for(j = 0; j < nc; j++)
      ptr[(ri - 1)*nc + j] = ptr[(rj - 1) * nc + j];
    for(j = 0; j < nc; j++)
      ptr[(rj - 1)*nc + j] = tempe[j];
    delete [] tempe;
  }

  void Array::excol(int ci, int cj) { // exchange the ci column and cj column, ci and cj are from 1.
    int i;
    double *tempe = new double [nl];
    for(i = 0; i < nl; i++)
      tempe[i] = ptr[i * nc + ci - 1];
    for(i = 0; i < nl; i++)
      ptr[i * nc + ci - 1] = ptr[i * nc + cj - 1];
    for(i = 0; i < nl; i++)
      ptr[i * nc + cj - 1] = tempe[i];
    delete [] tempe;
  }
  //****************
  void Array:: assort_a(int ci) { //assort data according to the ci column in acending order, ci is from 1
    int i, j;
    for(j = 0; j < nl - 1; j++)
      for(i = j + 1; i < nl; i++)
        if(ptr[i * nc + ci - 1] < ptr[j * nc + ci - 1]) exrow(i + 1, j + 1);
  }

  void Array:: assort_d(int ci) { //assort data according to the ci column in declining order, ci is from 1
    int i, j;
    for(j = 0; j < nl - 1; j++)
      for(i = j + 1; i < nl; i++)
        if(ptr[i * nc + ci - 1] > ptr[j * nc + ci - 1]) exrow(i + 1, j + 1);
  }

  //******************************************************************
  //find hull base for n dimensional case, n=2 or 3 , for the case that only find hull points, not edges or facets

  Array Array::hull_base_p(int n) {
    Array r_base(1, n);
    int i, j;
    bool h;
    if(n == 2) {
      double area_value;
      Array area(3, 3);
      area.ptr[6] = area.ptr[7] = area.ptr[8] = 1.0;
      r_base.ptr[0] = ptr[0];
      for(i = 1; i < nl; i++) //find the starting point with lowest y, and put it into 1st row
        if(ptr[i * 3 + 2] < ptr[2]) {
          r_base.ptr[0] = ptr[i * 3 + 0]; // label of point A is stored in p[0]
          exrow(1, i + 1);
        }
      area.ptr[0] = ptr[1];  //the 1st column of matrix area stores the determined hull point A
      area.ptr[3] = ptr[2];
      for(i = 1; i < nl; i++) {
        h = true;
        area.ptr[1] = ptr[i * 3 + 1]; //the 2nd column of matrix area stores the hull point candidate B
        area.ptr[4] = ptr[i * 3 + 2];
        for(j = 0; j < nl; j++) {
          if((j == i) || (j == 0)) continue;
          area.ptr[2] = ptr[j * 3 + 1]; //the 3rd column of matrix area stores the third point C
          area.ptr[5] = ptr[j * 3 + 2];
          area_value = area.det();
          if(area_value == 0) continue;
          else if(area_value < 0) {
            h = false;  // if for all other points except A and B, the area is positive or zero
            break;
          }
          else continue;
        }                                     // that means all points C is to the left of line AB, then keep B as a
        // hull point
        if(h == true) {
          r_base.ptr[1] = int(ptr[i * 3 + 0]);              // label of point B is stored in p[1]
          break;
        }
      }
      return r_base;
    }
    else {
      Array re_h_base(1, n - 1);
      double deter_value;
      int s, q, count;
      int sign;
      int pre_sign;
      Array deter((n + 1), (n + 1));
      bool same;
      Array reduce(nl, n);
      for(i = 0; i < nl; i++)
        for(j = 0; j < n; j++)
          reduce.ptr[i * n + j] = ptr[i * (n + 1) + j];
      re_h_base = reduce.hull_base(n - 1);
      for(j = 0; j < n - 1; j++)
        r_base.ptr[j] = re_h_base.ptr[j];
      for(j = 0; j < n - 1; j++)
        for(i = 0; i < nl; i++)
          if(r_base.ptr[j] == ptr[i * (n + 1) + 0]) exrow(i + 1, j + 1);
      for(j = 0; j < n + 1; j++)
        deter.ptr[n * (n + 1) + j] = 1.0;
      for(j = 0; j < n - 1; j++)
        for(i = 0; i < n; i++)
          deter.ptr[i * (n + 1) + j] = ptr[j * (n + 1) + i + 1];

      for(i = n - 1; i < nl; i++) {
        h = true;
        count = 0;
        for(j = 0; j < n; j++)
          deter.ptr[j * (n + 1) + n - 1] = ptr[i * (n + 1) + j + 1]; //the 2nd last column of matrix area stores the hull point candidate B
        for(j = 0; j < nl; j++) {
          same = false;
          for(q = 0; q < n - 1; q++)
            if((j == q) || (j == i))
              same = true;
          if(same) continue;
          count++;
          for(s = 0; s < n; s++)      // the last column of deter matrix
            deter.ptr[s * (n + 1) + n] = ptr[j * (n + 1) + s + 1];
          deter_value = deter.det();
          if(deter_value < 0) sign = -1;
          else if(deter_value > 0) sign = 1;
          else {
            count--;
            continue;
          }
          if(count == 1) {
            pre_sign = sign;
            continue;
          }
          if(sign != pre_sign) {
            h = false;
            break;
          }
          pre_sign = sign;
        }
        if(h == true) {
          r_base.ptr[n - 1] = int(ptr[i * (n + 1) + 0]);         // label of nth point of base is stored in p[n-1]
          break;
        }
      }
      return r_base;
    }

  }

  //***************************************************************
  // find hull of n-dimension case, for the case that only find hull points, not edges or facets
  // the code will return an matrix containing the label and coordinates of hull points
  // the 1st column is the hull point label.
  Array Array::hullnd_p(int n) { // n is the dimension = 2,or 3
    int screen;
    int *label = new int [n];
    bool h, same;
    Array tranA(1, nc);
    Array tran2A(1, nc);
    double deter2_value, theta;
    double *deter2 = new double [(n + 1) * (n + 1)];
    Array deter2A(n + 1, n + 1);
    double *sub_sys = new double [nc * (n - 1)];
    double *tran = new double [nc];
    double *tran2 = new double [nc];
    queue <Array>  remain_queue;
    stack <Array>  hull_stack;
    queue <Array>  new_elem;
    queue <Array>  facet;
    int i, j, r, s, x, y, z, zz, pp, ppp, remain_size, sign, pre_sign, count, remain_i;
    int num_hull = n; //initial number of hull points, hull base points
    int num_facet = 1;
    Array h_base(1, n);
    h_base = hull_base(n);
    facet.push(h_base);
    for(x = 0; x < n + 1; x++)
      deter2[n * (n + 1) + x] = 1.0;

    for(j = 0; j < n; j++)           // exchange base points into the first n rows
      for(i = 0; i < nl; i++)
        if(h_base.ptr[j] == ptr[i * (n + 1) + 0]) exrow(i + 1, j + 1);
    for(i = 0; i < n; i++) {
      for(j = 0; j < nc; j++)
        tran[j] = ptr[i * nc + j];
      tranA.setArray(tran);
      hull_stack.push(tranA);    // push hull base points to hull_stack
    }

    for(i = n; i < nl; i++) {
      for(j = 0; j < nc; j++)
        tran[j] = ptr[i * nc + j];
      tranA.setArray(tran);
      remain_queue.push(tranA);    // push remaining points to remain_queue
    }
    for(i = 0; i < n; i++)
      for(j = 0; j < n; j++) {
        if(j == i) continue;
        else {
          for(r = 0; r < nc; r++)
            tran[r] = ptr[j * nc + r];
          tranA.setArray(tran);
          new_elem.push(tranA);
        }  // push new sub_element from hull base to queue new_elem
      }


    while(!new_elem.empty()) {
      remain_size = remain_queue.size();
      if(remain_queue.empty()) break;

      for(s = 0; s < n - 1; s++) {
        tranA = new_elem.front();
        new_elem.pop();
        for(r = 0; r < nc; r++)
          sub_sys[s * nc + r] = tranA.ptr[r];
      }

      for(x = 0; x < n - 1; x++) {
        for(y = 0; y < n; y++)
          deter2[y * (n + 1) + x] = sub_sys[x * (n + 1) + y + 1];
        label[x] = int(sub_sys[x * (n + 1) + 0]);
      }

      for(remain_i = 0; remain_i < remain_size; remain_i++) { //1
        tranA = remain_queue.front(); // get a candidate from remain_queue, don't remove it from the queue at this time
        if(n == 3) {
          theta = ((sub_sys[1] - tranA.ptr[1]) * (sub_sys[5] - tranA.ptr[1]) + (sub_sys[2] - tranA.ptr[2]) * (sub_sys[6] - tranA.ptr[2]) + (sub_sys[7] - tranA.ptr[3]) * (sub_sys[3] - tranA.ptr[3])) /
                  (sqrt((sub_sys[1] - tranA.ptr[1]) * (sub_sys[1] - tranA.ptr[1]) + (sub_sys[2] - tranA.ptr[2]) * (sub_sys[2] - tranA.ptr[2]) + (sub_sys[3] - tranA.ptr[3]) * (sub_sys[3] - tranA.ptr[3])) *
                   sqrt((sub_sys[5] - tranA.ptr[1]) * (sub_sys[5] - tranA.ptr[1]) + (sub_sys[6] - tranA.ptr[2]) * (sub_sys[6] - tranA.ptr[2]) + (sub_sys[7] - tranA.ptr[3]) * (sub_sys[7] - tranA.ptr[3])));
          if((theta == 1.0) || (theta == -1.0)) {
            remain_queue.push(tranA);
            remain_queue.pop();
            continue;
          }
        }	 // make sure these three points are not on a line
        for(y = 0; y < n; y++)
          deter2[y * (n + 1) + n - 1] = tranA.ptr[0 * (n + 1) + y + 1]; // store hull candidate
        label[n - 1] = int(tranA.ptr[0 * (n + 1) + 0]);     //store label of candidate
        h = true;
        count = 0;
        for(j = 0; j < nl; j++) {
          same = false;
          sign = 0;
          for(int q = 0; q < n; q++)
            if(ptr[j * nc + 0] == label[q])
              same = true;
          if(same) continue;
          count++;
          for(s = 0; s < n; s++)      // the last column of deter matrix
            deter2[s * (n + 1) + n] = ptr[j * (n + 1) + s + 1];
          deter2A.setArray(deter2);
          deter2_value = deter2A.det();
          if(deter2_value < 0) sign = -1;
          else if(deter2_value > 0) sign = 1;
          else {
            count--;
            continue;
          }
          if(count == 1) {
            pre_sign = sign;
            continue;
          }
          if(sign != pre_sign) {
            h = false;
            tranA = remain_queue.front();
            remain_queue.push(tranA);
            remain_queue.pop();
            break;
          }
          pre_sign = sign;
        }
        if(h == true) {
          num_hull++;
          for(zz = 0; zz < n - 1; zz++)
            h_base.ptr[zz] = sub_sys[zz * nc + 0];
          tranA = remain_queue.front();
          h_base.ptr[n - 1] = tranA.ptr[0];
          facet.push(h_base);
          // cout<<h_base<<endl;
          num_facet++;
          hull_stack.push(tranA);
          remain_queue.pop();
          for(z = 0; z < n - 1; z++) {
            for(y = 0; y < n - 1; y++) {
              if(y == z) continue;
              for(x = 0; x < nc; x++)
                tran2[x] = sub_sys[y * nc + x];
              tran2A.setArray(tran2);
              new_elem.push(tran2A);
            }
            new_elem.push(tranA);
          }
          //  break;
        }
      }// 1
    }

    // find the in facet points

    for(ppp = 0; ppp < num_facet; ppp++) {
      if(remain_queue.empty()) break;
      remain_size = remain_queue.size();
      h_base = facet.front();  // get a facet from the facet queue
      facet.push(h_base);
      for(x = 0; x < n; x++)
        for(y = 0; y < nl; y++)
          if(ptr[y * nc + 0] == h_base.ptr[x]) {
            for(z = 0; z < n; z++)
              deter2A.ptr[z * (n + 1) + x] = ptr[y * nc + z + 1];
            break;
          }                         //get one facet and put the corresponding points in the deter matrix

      for(pp = 0; pp < remain_size; pp++) {
        tranA = remain_queue.front(); //get one point from remain queue.
        for(z = 0; z < n; z++)                     // put the coordinates in the last column of deter matrix
          deter2A.ptr[z * (n + 1) + n] = tranA.ptr[z + 1];
        if(deter2A.det() != 0) {          // the point is not on the facet
          remain_queue.push(tranA);
          remain_queue.pop();
          continue;				// get the next point from remain queue
        }
        else {                         // the point is on the facet

          hull_stack.push(tranA);          // push the point in the hull stack
          cout << tranA.ptr[0] << "   " << h_base << endl;
          num_hull++;
          remain_queue.pop();              // delete the point from the remain queue
        }
      }
      facet.pop();   //delete the old facet
    }

    Array hull_resultA(num_hull, nc);
    double *hull_result = new double [num_hull * nc];
    for(j = 0; j < num_hull; j++) {
      tranA = hull_stack.top();
      for(i = 0; i < nc; i++)
        hull_result[j * nc + i] = tranA.elem(1, i + 1);
      hull_stack.pop();
    }
    hull_resultA.setArray(hull_result);
    delete [] label;
    delete [] deter2;
    delete [] sub_sys;
    delete [] tran;
    delete [] tran2;
    delete [] hull_result;
    return hull_resultA;
  }

  //******************************************************************

  //find hull base for n dimensional case, n=2 or 3, for the case that find both hull points and facets or edges

  Array Array::hull_base(int n) {
    Array r_base(1, n);
    int i, j;
    bool h;
    if(n == 2) {
      double area_value, distance, pre_distance;
      int cand_size;
      queue <Array> cand;
      Array temp2d(1, 3);
      Array area(3, 3);
      area.ptr[6] = area.ptr[7] = area.ptr[8] = 1.0;
      r_base.ptr[0] = ptr[0];
      for(i = 1; i < nl; i++) //find the starting point with lowest y, and put it into 1st row
        if(ptr[i * 3 + 2] < ptr[2]) {
          r_base.ptr[0] = ptr[i * 3 + 0]; // label of point A is stored in p[0]
          exrow(1, i + 1);
        }
      area.ptr[0] = ptr[1];  //the 1st column of matrix area stores the determined hull point A
      area.ptr[3] = ptr[2];
      for(i = 1; i < nl; i++) {
        h = true;
        area.ptr[1] = ptr[i * 3 + 1]; //the 2nd column of matrix area stores the hull point candidate B
        area.ptr[4] = ptr[i * 3 + 2];
        for(j = 0; j < nl; j++) {
          if((j == i) || (j == 0)) continue;
          area.ptr[2] = ptr[j * 3 + 1]; //the 3rd column of matrix area stores the third point C
          area.ptr[5] = ptr[j * 3 + 2];
          area_value = area.det();
          if(area_value == 0) continue;
          else if(area_value < 0) {
            h = false;  // if for all other points except A and B, the area is positive or zero
            break;
          }
          else continue;
        }                                     // that means all points C is to the left of line AB, then keep B as a
        // hull point
        if(h == true) {
          temp2d.ptr[0] = ptr[i * 3 + 0];
          temp2d.ptr[1] = ptr[i * 3 + 1];
          temp2d.ptr[2] = ptr[i * 3 + 2];
          cand.push(temp2d);                // point B is the hull point associated with point A, continue to find such point of A
          // r_base.ptr[1]=int(ptr[i*3+0]);                    // label of point B is stored in p[1]
          //  break;
        }
      }
      cand_size = cand.size();                   //find the hull B that is farest from A
      if(cand_size == 1) {
        temp2d = cand.front();
        r_base.ptr[1] = int(temp2d.ptr[0]);
        return r_base;
      }
      else {
        pre_distance = 0.0;
        for(i = 0; i < cand_size; i++) {
          temp2d = cand.front();
          distance = sqrt((temp2d.ptr[1] - ptr[1]) * (temp2d.ptr[1] - ptr[1]) + (temp2d.ptr[2] - ptr[2]) * (temp2d.ptr[2] - ptr[2]));
          if(distance > pre_distance) {
            r_base.ptr[1] = temp2d.ptr[0];
            pre_distance = distance;
          }
          cand.pop();
        }
        return r_base;
      }
    }
    else {
      Array re_h_base(1, n - 1);
      queue <Array> cand_3d;
      Array temp3d(1, n + 1);
      Array p1(1, n), p2(1, n), p3(1, n);
      int cand_3d_size;
      double distance_3d, pre_distance_3d;
      double deter_value;
      int s, q, count, qq;
      int sign;
      int pre_sign;
      Array deter((n + 1), (n + 1));
      bool same;
      Array reduce(nl, n);
      for(i = 0; i < nl; i++)
        for(j = 0; j < n; j++)
          reduce.ptr[i * n + j] = ptr[i * (n + 1) + j];
      re_h_base = reduce.hull_base(n - 1);
      for(j = 0; j < n - 1; j++)
        r_base.ptr[j] = re_h_base.ptr[j];
      for(j = 0; j < n - 1; j++)
        for(i = 0; i < nl; i++)
          if(r_base.ptr[j] == ptr[i * (n + 1) + 0]) exrow(i + 1, j + 1);
      for(j = 0; j < n + 1; j++)
        deter.ptr[n * (n + 1) + j] = 1.0;
      for(j = 0; j < n - 1; j++)
        for(i = 0; i < n; i++)
          deter.ptr[i * (n + 1) + j] = ptr[j * (n + 1) + i + 1];

      for(i = n - 1; i < nl; i++) {
        h = true;
        count = 0;
        for(j = 0; j < n; j++)
          deter.ptr[j * (n + 1) + n - 1] = ptr[i * (n + 1) + j + 1]; //the 2nd last column of matrix area stores the hull point candidate B
        for(j = 0; j < nl; j++) {
          same = false;
          for(q = 0; q < n - 1; q++)
            if((j == q) || (j == i))
              same = true;
          if(same) continue;
          count++;
          for(s = 0; s < n; s++)      // the last column of deter matrix
            deter.ptr[s * (n + 1) + n] = ptr[j * (n + 1) + s + 1];
          deter_value = deter.det();
          if(deter_value < 0) sign = -1;
          else if(deter_value > 0) sign = 1;
          else {
            count--;
            continue;
          }
          if(count == 1) {
            pre_sign = sign;
            continue;
          }
          if(sign != pre_sign) {
            h = false;
            break;
          }
          pre_sign = sign;
        }
        if(h == true) {
          for(qq = 0; qq < n + 1; qq++)
            temp3d.ptr[qq] = ptr[i * (n + 1) + qq];
          //   r_base.ptr[n-1]=int(ptr[i*(n+1)+0]);                   // label of nth point of base is stored in p[n-1]
          //    break;
          cand_3d.push(temp3d);
        }
      }
      cand_3d_size = cand_3d.size();
      if(cand_3d_size == 1) {
        temp3d = cand_3d.front();
        r_base.ptr[n - 1] = int(temp3d.ptr[0]);
        return r_base;
      }
      else {
        pre_distance_3d = 0.0;
        for(i = 0; i < cand_3d_size; i++) {
          temp3d = cand_3d.front();
          for(qq = 0; qq < 3; qq++) {
            p1.ptr[qq] = ptr[qq + 1];
            p2.ptr[qq] = ptr[1 * 4 + qq + 1];
            p3.ptr[qq] = temp3d.ptr[qq + 1];
          }
          distance_3d = distance_pl(p3, p1, p2);
          if(distance_3d > pre_distance_3d) {
            r_base.ptr[2] = temp3d.ptr[0];
            pre_distance_3d = distance_3d;
          }
          cand_3d.pop();
        }
        return r_base;
      }
    }

  }

  //***************************************************************
  // find hull of n-dimension case, for the case that fint both hull points and facets or edges
  // the code will return an matrix containing the label and coordinates of hull points
  // the 1st column is the hull point label.
  Array Array::hullnd(int n, Array &f) { // n is the dimension = 2,or 3
    int  mark;
    int *label = new int [n];
    bool h, same, same_point, h_out;
    bool cut_value;
    Array tranA(1, nc), tranA_2(1, nc), keep(1, nc);
    Array tran2A(1, nc);
    double deter2_value, theta, dis, pre_dis;
    double *deter2 = new double [(n + 1) * (n + 1)];
    Array deter2A(n + 1, n + 1);
    double *sub_sys = new double [nc * (n - 1)];
    double *tran = new double [nc];
    double *tran2 = new double [nc];
    queue <Array>  remain_queue;
    queue <Array>  hull_queue;
    queue <Array>  new_elem;
    queue <Array>  facet;
    queue <Array> cut_facet;
    Array cut_facetA(1, 3);
    stack <Array>  facet_stack;
    Array facetA(1, 3);  //stores the label of three points
    Array edgeA(1, 2);    //stores the label of two points
    Array temp_edge(1, 2);
    queue <Array> edge;
    Array p1(1, 3), p2(1, 3), p3(1, 3), p4(1, 3);
    Array p11(1, 3), p22(1, 3), p33(1, 3), p44(1, 3), p55(1, 3);
    int keep_label, m, i, j, r, ss, s, x, y, z, zz, pp, ppp, remain_size, sign, pre_sign, count, remain_i, hull_queue_size, new_i, edge_size, num_edge, num_appear;
    int num_hull = n; //initial number of hull points, hull base points
    int num_facet = 1, cut_facetA_size;
    Array h_base(1, n);
    h_base = hull_base(n);
    facet.push(h_base);              //push hull base into facet queue.
    facet_stack.push(h_base);
    cut_facet.push(h_base);
    for(x = 0; x < n + 1; x++)
      deter2[n * (n + 1) + x] = 1.0;

    for(j = 0; j < n; j++)           // exchange base points into the first n rows
      for(i = 0; i < nl; i++)
        if(h_base.ptr[j] == ptr[i * (n + 1) + 0]) exrow(i + 1, j + 1);
    for(i = 0; i < n; i++) {
      for(j = 0; j < nc; j++)
        tran[j] = ptr[i * nc + j];
      tranA.setArray(tran);
      hull_queue.push(tranA);    // push hull base points to hull_queue
    }

    if(n == 3) {
      edgeA.ptr[0] = ptr[0 * 4 + 0];
      edgeA.ptr[1] = ptr[1 * 4 + 0];
      edge.push(edgeA);
      edgeA.ptr[0] = ptr[0 * 4 + 0];
      edgeA.ptr[1] = ptr[2 * 4 + 0];
      edge.push(edgeA);
      edgeA.ptr[0] = ptr[1 * 4 + 0];
      edgeA.ptr[1] = ptr[2 * 4 + 0];
      edge.push(edgeA);         //push initial three edges

      while(!facet.empty()) {
        facetA = facet.front(); //get one facet from the facet queue;
        // cout<<"the facet is "<<facetA;
        // cout<<"the edge in queue are"<<endl;
        edge_size = edge.size();

        num_edge = 0; //the edge of one facet, which is new or just appear one time

        for(i = 0; i < nl; i++)
          if(ptr[i * 4 + 0] == facetA.ptr[0]) {
            for(r = 0; r < nc; r++)
              tran[r] = ptr[i * nc + r];
            tranA.setArray(tran);
            break;
          }                                 // push new sub_element from hull base to queue new_elem
        edgeA.ptr[0] = tranA.ptr[0];
        for(i = 0; i < nl; i++)
          if(ptr[i * 4 + 0] == facetA.ptr[1]) {
            for(r = 0; r < nc; r++)
              tran[r] = ptr[i * nc + r];
            tranA_2.setArray(tran);
            break;
          }                                // push new sub_element from hull base to queue new_elem
        edgeA.ptr[1] = tranA_2.ptr[0];
        num_appear = 0;
        //cout<<"1st edge is "<<edgeA;
        for(i = 0; i < edge_size; i++) {
          temp_edge = edge.front();
          if(((edgeA.ptr[0] == temp_edge.ptr[0]) || (edgeA.ptr[0] == temp_edge.ptr[1])) && ((edgeA.ptr[1] == temp_edge.ptr[0]) || (edgeA.ptr[1] == temp_edge.ptr[1]))) num_appear++;
          edge.push(temp_edge);
          edge.pop();
        }
        if(num_appear >= 2) {
          ; //cout<<"1st edge is 2 times"<<endl
        }

        else {
          new_elem.push(tranA);  //push one edge into queue
          new_elem.push(tranA_2);
          num_edge++;
        }

        for(i = 0; i < nl; i++)
          if(ptr[i * 4 + 0] == facetA.ptr[0]) {
            for(r = 0; r < nc; r++)
              tran[r] = ptr[i * nc + r];
            tranA.setArray(tran);
            break;
          }                             // push new sub_element from hull base to queue new_elem
        edgeA.ptr[0] = tranA.ptr[0];
        for(i = 0; i < nl; i++)
          if(ptr[i * 4 + 0] == facetA.ptr[2]) {
            for(r = 0; r < nc; r++)
              tran[r] = ptr[i * nc + r];
            tranA_2.setArray(tran);
            break;
          }                            // push new sub_element from hull base to queue new_elem
        edgeA.ptr[1] = tranA_2.ptr[0];
        num_appear = 0;
        //cout<<"2nd edge is "<<edgeA;
        for(i = 0; i < edge_size; i++) {
          temp_edge = edge.front();
          if(((edgeA.ptr[0] == temp_edge.ptr[0]) || (edgeA.ptr[0] == temp_edge.ptr[1])) && ((edgeA.ptr[1] == temp_edge.ptr[0]) || (edgeA.ptr[1] == temp_edge.ptr[1]))) num_appear++;
          edge.push(temp_edge);
          edge.pop();
        }
        if(num_appear >= 2) {
          ;
        }//cout<<"2nd edge is 2 times"<<endl;}
        else {
          new_elem.push(tranA);
          new_elem.push(tranA_2);
          num_edge++;
        }

        for(i = 0; i < nl; i++)
          if(ptr[i * 4 + 0] == facetA.ptr[1]) {
            for(r = 0; r < nc; r++)
              tran[r] = ptr[i * nc + r];
            tranA.setArray(tran);
            break;
          }                              // push new sub_element from hull base to queue new_elem
        edgeA.ptr[0] = tranA.ptr[0];
        for(i = 0; i < nl; i++)
          if(ptr[i * 4 + 0] == facetA.ptr[2]) {
            for(r = 0; r < nc; r++)
              tran[r] = ptr[i * nc + r];
            tranA_2.setArray(tran);
            break;
          }                             // push new sub_element from hull base to queue new_elem
        edgeA.ptr[1] = tranA_2.ptr[0];
        num_appear = 0;
        //cout<<"3rd edge is "<<edgeA;
        for(i = 0; i < edge_size; i++) {
          temp_edge = edge.front();
          if(((edgeA.ptr[0] == temp_edge.ptr[0]) || (edgeA.ptr[0] == temp_edge.ptr[1])) && ((edgeA.ptr[1] == temp_edge.ptr[0]) || (edgeA.ptr[1] == temp_edge.ptr[1]))) num_appear++;
          edge.push(temp_edge);
          edge.pop();
        }
        if(num_appear >= 2) {
          ;
        }//cout<<"3rd edge is 2 times"<<endl;}
        else {
          new_elem.push(tranA);
          new_elem.push(tranA_2);
          num_edge++;
        }

        if(new_elem.empty()) {
          facet.pop();
          continue;
        }
        //cout<<num_edge<<endl;
        for(new_i = 0; new_i < num_edge; new_i++) { //get one edge each time
          for(s = 0; s < 2; s++) {
            tranA = new_elem.front();
            new_elem.pop();
            for(r = 0; r < nc; r++)
              sub_sys[s * nc + r] = tranA.ptr[r]; //get one edge to sub_sys
          }

          for(x = 0; x < 2; x++) {
            for(y = 0; y < 3; y++)
              deter2[y * 4 + x] = sub_sys[x * 4 + y + 1];
            label[x] = int(sub_sys[x * 4 + 0]);
          }
          pre_dis = 0.0;
          //cout<<"the edge is "<<sub_sys[0]<<"   "<<sub_sys[4]<<endl;

          for(i = 0; i < nl; i++) {
            for(j = 0; j < nc; j++)
              tran[j] = ptr[i * nc + j];
            tranA.setArray(tran);
            remain_queue.push(tranA);    // push all points to remain_queue
          }

          while(!remain_queue.empty()) {     //1
            tranA = remain_queue.front();                // get a candidate from remain_queue, don't remove it from the queue at this time
            if((tranA.ptr[0] == facetA.ptr[0]) || (tranA.ptr[0] == facetA.ptr[1]) || (tranA.ptr[0] == facetA.ptr[2])) {
              remain_queue.pop();
              continue;
            }
            theta = ((sub_sys[1] - tranA.ptr[1]) * (sub_sys[5] - tranA.ptr[1]) + (sub_sys[2] - tranA.ptr[2]) * (sub_sys[6] - tranA.ptr[2]) + (sub_sys[7] - tranA.ptr[3]) * (sub_sys[3] - tranA.ptr[3])) /
                    (sqrt((sub_sys[1] - tranA.ptr[1]) * (sub_sys[1] - tranA.ptr[1]) + (sub_sys[2] - tranA.ptr[2]) * (sub_sys[2] - tranA.ptr[2]) + (sub_sys[3] - tranA.ptr[3]) * (sub_sys[3] - tranA.ptr[3])) *
                     sqrt((sub_sys[5] - tranA.ptr[1]) * (sub_sys[5] - tranA.ptr[1]) + (sub_sys[6] - tranA.ptr[2]) * (sub_sys[6] - tranA.ptr[2]) + (sub_sys[7] - tranA.ptr[3]) * (sub_sys[7] - tranA.ptr[3])));
            if((theta == 1.0) || (theta == -1.0)) {
              remain_queue.pop();
              continue;
            }	// make sure these three points are not on a line

            for(y = 0; y < 3; y++)
              deter2[y * 4 + 2] = tranA.ptr[0 * 4 + y + 1]; // store hull candidate
            label[2] = int(tranA.ptr[0 * 4 + 0]);         //store label of candidate
            h = true;
            count = 0;

            for(j = 0; j < nl; j++) {
              same = false;
              sign = 0;
              for(int q = 0; q < 3; q++)
                if(ptr[j * nc + 0] == label[q])
                  same = true;
              if(same) continue;
              count++;
              for(s = 0; s < n; s++)      // the last column of deter matrix
                deter2[s * 4 + 3] = ptr[j * 4 + s + 1];
              deter2A.setArray(deter2);
              deter2_value = deter2A.det();
              if(deter2_value < 0) sign = -1;
              else if(deter2_value > 0) sign = 1;
              else {
                count--;
                continue;
              }
              if(count == 1) {
                pre_sign = sign;
                continue;
              }
              if(sign != pre_sign) {
                h = false;
                remain_queue.pop();
                break;
              }
              pre_sign = sign;
            }
            if(h == true) {
              p4.ptr[0] = tranA.ptr[1];
              p4.ptr[1] = tranA.ptr[2];
              p4.ptr[2] = tranA.ptr[3];
              p1.ptr[0] = sub_sys[1]; //p1 and p2 are two end points of the edge selected.
              p1.ptr[1] = sub_sys[2];
              p1.ptr[2] = sub_sys[3];
              p2.ptr[0] = sub_sys[5];
              p2.ptr[1] = sub_sys[6];
              p2.ptr[2] = sub_sys[7];
              for(i = 0; i < 3; i++) {
                if((facetA.ptr[i] == sub_sys[0]) || (facetA.ptr[i] == sub_sys[4])) continue;
                else mark = int(facetA.ptr[i]);
              }
              for(i = 0; i < nl; i++)
                if(mark == ptr[i * 4 + 0]) {
                  p3.ptr[0] = ptr[i * 4 + 1];      //p3 is another point of the selected facet
                  p3.ptr[1] = ptr[i * 4 + 2];
                  p3.ptr[2] = ptr[i * 4 + 3];
                  break;
                }

              if(on_plane(p4, p1, p2, p3)) {
                //	    cout<<"point "<<tranA.ptr[0]<<" on_plane"<<endl;            //if p4 is on the plane p1-p2-p3
                if(in_tri(p4, p1, p2, p3)) {
                  remain_queue.pop();  //if p4 is within triangle p1-p2-p3
                  continue;
                }
                else if(same_side(p4, p3, p1, p2)) {
                  remain_queue.pop();  //if p4 and p3 are to the same side of edge p1-p2
                  continue;
                }
                else {};
              }

              ////////////////
              cut_value = false;
              for(i = 0; i < 2; i++) {
                if(i == 0) {
                  p11 = p1;  //get one new edge;
                  p22 = p4;
                }
                else {
                  p11 = p2;  //get another new edge;
                  p22 = p4;
                }

                cut_facetA_size = cut_facet.size();
                for(pp = 0; pp < cut_facetA_size; pp++) {
                  cut_facetA = cut_facet.front();       // get one facet
                  for(s = 0; s < nl; s++)
                    if(ptr[s * nc + 0] == cut_facetA.ptr[0]) {
                      for(ss = 0; ss < nc; ss++)
                        tranA_2.ptr[ss] = ptr[s * nc + ss];
                    }
                  p33.ptr[0] = tranA_2.ptr[1];
                  p33.ptr[1] = tranA_2.ptr[2];
                  p33.ptr[2] = tranA_2.ptr[3];

                  for(s = 0; s < nl; s++)
                    if(ptr[s * nc + 0] == cut_facetA.ptr[1]) {
                      for(ss = 0; ss < nc; ss++)
                        tranA_2.ptr[ss] = ptr[s * nc + ss];
                    }
                  p44.ptr[0] = tranA_2.ptr[1];
                  p44.ptr[1] = tranA_2.ptr[2];
                  p44.ptr[2] = tranA_2.ptr[3];

                  for(s = 0; s < nl; s++)
                    if(ptr[s * nc + 0] == cut_facetA.ptr[2]) {
                      for(ss = 0; ss < nc; ss++)
                        tranA_2.ptr[ss] = ptr[s * nc + ss];
                    }
                  p55.ptr[0] = tranA_2.ptr[1];
                  p55.ptr[1] = tranA_2.ptr[2];
                  p55.ptr[2] = tranA_2.ptr[3];

                  if(cut(p33, p44, p55, p11, p22)) {
                    cut_value = true;  //if the new edge cut one existing facet
                    break;
                  }
                  cut_facet.push(cut_facetA);
                  cut_facet.pop();
                }
                if(cut_value == true) break;
              }
              if(cut_value == true) {
                remain_queue.pop();
                continue;
              }
              /////////////////

              dis = distance_pl(p4, p1, p2);        //calculate the distance between point p4 and edge p1-p2
              if(dis > pre_dis) {
                keep = tranA;
                pre_dis = dis;
                remain_queue.pop();
              }
              else {
                remain_queue.pop();
              }
            }
          }//1
          h_base.ptr[0] = sub_sys[0];
          h_base.ptr[1] = sub_sys[4];
          h_base.ptr[2] = keep.ptr[0];
          facet.push(h_base);
          facet_stack.push(h_base);
          cut_facet.push(h_base);
          num_facet++;                    //construct new facet
          //	cout<<"the obtained facet is "<<h_base;
          edgeA.ptr[0] = sub_sys[0];                   //create new edges
          edgeA.ptr[1] = sub_sys[4];
          edge.push(edgeA);
          edgeA.ptr[0] = sub_sys[0];
          edgeA.ptr[1] = keep.ptr[0];
          edge.push(edgeA);
          edgeA.ptr[0] = sub_sys[4];
          edgeA.ptr[1] = keep.ptr[0];
          edge.push(edgeA);

          same_point = false;
          hull_queue_size = hull_queue.size();
          for(i = 0; i < hull_queue_size; i++) {
            tranA = hull_queue.front();
            if(tranA.ptr[0] == keep.ptr[0]) {
              same_point = true;
              break;
            }
            else {
              hull_queue.push(tranA);
              hull_queue.pop();
            }
          }
          if(same_point == false) {
            hull_queue.push(keep);
            num_hull++;
          }
        }
        facet.pop();
      }

      Array facet_result(num_facet, 3);
      for(i = 0; i < num_facet; i++) {
        facetA = facet_stack.top();
        for(j = 0; j < 3; j++)
          facet_result.ptr[i * 3 + j] = facetA.ptr[j];
        facet_stack.pop();
      }

      // find the in facet points
      while(!remain_queue.empty()) remain_queue.pop();  //clear remain queue
      for(i = 0; i < nl; i++) {
        same = false;
        for(j = 0; j < num_hull; j++) {
          tranA = hull_queue.front();
          if(ptr[i * 4 + 0] == tranA.ptr[0]) {
            hull_queue.push(tranA);
            hull_queue.pop();
            same = true;
            break;
          }
          else {
            hull_queue.push(tranA);
            hull_queue.pop();
          }
        }
        if(same == false) {
          tranA.ptr[0] = ptr[i * 4 + 0];
          tranA.ptr[1] = ptr[i * 4 + 1];
          tranA.ptr[2] = ptr[i * 4 + 2];
          tranA.ptr[3] = ptr[i * 4 + 3];
          remain_queue.push(tranA);
        }
      }                                      //push points other than hull points into queue


      for(ppp = 0; ppp < num_facet; ppp++) {
        if(remain_queue.empty()) break;
        remain_size = remain_queue.size();
        h_base.ptr[0] = facet_result.ptr[ppp * 3 + 0];
        h_base.ptr[1] = facet_result.ptr[ppp * 3 + 1];
        h_base.ptr[2] = facet_result.ptr[ppp * 3 + 2]; // get a facet from the facet result
        for(x = 0; x < 3; x++)
          for(y = 0; y < nl; y++)
            if(ptr[y * nc + 0] == h_base.ptr[x]) {
              for(z = 0; z < n; z++)
                deter2A.ptr[z * 4 + x] = ptr[y * nc + z + 1];
              break;
            }                         //get one facet and put the corresponding points in the deter matrix
        p1.ptr[0] = deter2A.ptr[0];
        p1.ptr[1] = deter2A.ptr[4];
        p1.ptr[2] = deter2A.ptr[8];
        p2.ptr[0] = deter2A.ptr[1];
        p2.ptr[1] = deter2A.ptr[5];
        p2.ptr[2] = deter2A.ptr[9];
        p3.ptr[0] = deter2A.ptr[2];
        p3.ptr[1] = deter2A.ptr[6];
        p3.ptr[2] = deter2A.ptr[10];
        for(pp = 0; pp < remain_size; pp++) {
          tranA = remain_queue.front(); //get one point from remain queue.
          for(z = 0; z < 3; z++)
            p4.ptr[z] = tranA.ptr[z + 1];
          if(!in_tri(p4, p1, p2, p3)) {         // the point is not on the facet
            remain_queue.push(tranA);
            remain_queue.pop();
            continue;				// get the next point from remain queue
          }
          else {                         // the point is on the facet
            hull_queue.push(tranA);          // push the point in the hull stack
            cout << "Point " << setw(4) << tranA.ptr[0] << ", " << tranA.ptr[1] << "  is on facet  " << h_base;
            num_hull++;
            remain_queue.pop();              // delete the point from the remain queue
          }
        }
      }

      Array hull_resultA(num_hull, nc);
      for(j = 0; j < num_hull; j++) {
        tranA = hull_queue.front();
        for(i = 0; i < nc; i++)
          hull_resultA.ptr[j * nc + i] = tranA.ptr[i];
        hull_queue.push(tranA);
        hull_queue.pop();
      }

      delete [] label;
      delete [] deter2;
      delete [] sub_sys;
      delete [] tran;
      delete [] tran2;

      //cout<<endl<<"The hull facets are:"<<endl<< facet_result<<endl;
      f = facet_result;
      return hull_resultA;
    }
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    else if(n == 2) {
      edgeA.ptr[0] = ptr[0];
      edgeA.ptr[1] = ptr[3];
      edge.push(edgeA);
      num_edge = 1;
      deter2A.ptr[6] = deter2A.ptr[7] = deter2A.ptr[8] = 1.0;
      for(m = 1; m < nl; m++) {
        h_out = false;
        pre_dis = 0.0;
        deter2A.ptr[0] = ptr[m * 3 + 1]; //the 1st column of matrix area stores the determined hull point A
        deter2A.ptr[3] = ptr[m * 3 + 2];
        edgeA.ptr[0] = ptr[m * 3 + 0]; //the label of hull point A is stored in edgeA[0]
        for(i = 0; i < nl; i++) { //find the hull point B associated with A
          h = true;
          deter2A.ptr[1] = ptr[i * 3 + 1]; //the 2nd column of matrix area stores the hull point candidate B
          deter2A.ptr[4] = ptr[i * 3 + 2];
          for(j = 0; j < nl; j++) {
            if((j == i) || (j == m)) continue;
            deter2A.ptr[2] = ptr[j * 3 + 1]; //the 3rd column of matrix area stores the third point C
            deter2A.ptr[5] = ptr[j * 3 + 2];
            if(deter2A.det() < 0) {
              h = false;  // if for all other points except A and B, the area is positive or zero
              break;
            }
          }                                     // that means all points C is to the left of line AB, then keep B as a
          // hull point
          if(h == true) {
            h_out = true;
            dis = sqrt((ptr[i * 3 + 1] - ptr[m * 3 + 1]) * (ptr[i * 3 + 1] - ptr[m * 3 + 1]) + (ptr[i * 3 + 2] - ptr[m * 3 + 2]) * (ptr[i * 3 + 2] - ptr[m * 3 + 2]));
            if(dis > pre_dis) {
              pre_dis = dis;
              keep_label = i;
              keep.ptr[0] = ptr[i * 3 + 0];
              keep.ptr[1] = ptr[i * 3 + 1];
              keep.ptr[2] = ptr[i * 3 + 2];
              continue;
            }
            else {
              continue;
            }
          }
        }
        if(h_out == true) {
          if(keep.ptr[0] == ptr[0]) {
            edgeA.ptr[1] = ptr[0];
            edge.push(edgeA);
            num_edge++;
            break;
          }
          else {
            num_hull++;    // the (k+1)th hull point is found
            exrow(keep_label + 1, num_hull);
            edgeA.ptr[1] = ptr[(num_hull - 1) * 3 + 0]; // thte label of hull point B is stored in edgaA[1]
            edge.push(edgeA);
            num_edge++;
          }
        }
        else {
          edgeA.ptr[1] = ptr[0];
          edge.push(edgeA);
          num_edge++;
          break;
        }
      }

      while(!hull_queue.empty()) hull_queue.pop();
      for(i = 0; i < num_hull; i++) { //push hull points into hull stack
        for(j = 0; j < 3; j++)
          tranA.ptr[j] = ptr[i * 3 + j];
        hull_queue.push(tranA);
      }
      for(i = num_hull; i < nl; i++) { //push remaining points into remain queue
        for(j = 0; j < 3; j++)
          tranA.ptr[j] = ptr[i * 3 + j];
        remain_queue.push(tranA);
      }

      for(ppp = 0; ppp < num_edge; ppp++) {
        if(remain_queue.empty()) break;
        remain_size = remain_queue.size();
        edgeA = edge.front();  // get a edge from the edge queue
        edge.push(edgeA);
        for(x = 0; x < 2; x++)
          for(y = 0; y < nl; y++)
            if(ptr[y * 3 + 0] == edgeA.ptr[x]) {
              for(z = 0; z < 2; z++)
                deter2A.ptr[z * 3 + x] = ptr[y * 3 + z + 1];
              break;
            }                         //get one edge and put the corresponding points in the deter matrix
        for(pp = 0; pp < remain_size; pp++) {
          tranA = remain_queue.front(); //get one point from remain queue.
          for(z = 0; z < 2; z++)                     // put the coordinates in the last column of deter matrix
            deter2A.ptr[z * 3 + 2] = tranA.ptr[z + 1];
          if(deter2A.det() != 0) {          // the point is not on the edge
            remain_queue.push(tranA);
            remain_queue.pop();
            continue;				// get the next point from remain queue
          }
          else {                         // the point is on the facet
            cout << "Point " << setw(4) << tranA.ptr[0] << "   is on edge " << edgeA;
            hull_queue.push(tranA);          // push the point in the hull queue
            num_hull++;
            remain_queue.pop();              // delete the point from the remain queue
          }
        }
        edge.pop();   //delete the old edge
      }

      Array hullA(num_hull, 3);  // the 1st column is the label of the hull points, the 2nd and 3rd are coordinates of the hull points
      for(i = 0; i < num_hull; i++) {
        tranA = hull_queue.front();
        for(j = 0; j < 3; j++)
          hullA.ptr[i * 3 + j] = tranA.ptr[j];
        hull_queue.pop();
      }
      Array edge_result(num_edge, 2);
      for(i = 0; i < num_edge; i++) {
        edgeA = edge.front();
        edge_result.ptr[i * 2 + 0] = edgeA.ptr[0];
        edge_result.ptr[i * 2 + 1] = edgeA.ptr[1];
        edge.push(edgeA);
        edge.pop();
      }

      //cout<<endl<<"The hull edges are:"<<endl<<edge_result<<endl;
      f = edge_result;
      return hullA;

    }
    //////////////////////////////////////////////////////////////////////////////////
    else {
      cout << " the dimension is not 2 or 3!" << endl;
      exit(1);
    }
  }
  //************************************************************************************************************************
  //calculate distance between point p3 and straight line p1p2
  double distance_pl(Array p3, Array p1, Array p2) {
    double A, B, C, D, x, y, z;
    A = p1.ptr[0] - p2.ptr[0];
    B = p1.ptr[1] - p2.ptr[1];
    C = p1.ptr[2] - p2.ptr[2];
    D = -A * p3.ptr[0] - B * p3.ptr[1] - C * p3.ptr[2];
    x = (p1.ptr[0] * B * B + p1.ptr[0] * C * C - A * B * p1.ptr[1] - A * C * p1.ptr[2] - A * D) / (A * A + B * B + C * C);
    y = (p1.ptr[1] * A * A + p1.ptr[1] * C * C - B * A * p1.ptr[0] - B * C * p1.ptr[2] - B * D) / (A * A + B * B + C * C);
    z = (p1.ptr[2] * A * A + p1.ptr[2] * B * B - C * A * p1.ptr[0] - C * B * p1.ptr[1] - C * D) / (A * A + B * B + C * C);
    return sqrt((x - p3.ptr[0]) * (x - p3.ptr[0]) + (y - p3.ptr[1]) * (y - p3.ptr[1]) + (z - p3.ptr[2]) * (z - p3.ptr[2]));
  }

  //check whether point p4 is within triangle p1p2p3 or not
  bool in_tri(Array p4, Array p1, Array p2, Array p3) {
    Array p1p4(1, 3), p2p4(1, 3), p3p4(1, 3), p1p2(1, 3), p2p3(1, 3), p3p1(1, 3), cross1(1, 3), cross2(1, 3), cross3(1, 3);
    double norm1, norm2, norm3, theta12, theta23, theta31;

    p1p4.ptr[0] = p4.ptr[0] - p1.ptr[0];
    p1p4.ptr[1] = p4.ptr[1] - p1.ptr[1];
    p1p4.ptr[2] = p4.ptr[2] - p1.ptr[2];
    p1p2.ptr[0] = p2.ptr[0] - p1.ptr[0];
    p1p2.ptr[1] = p2.ptr[1] - p1.ptr[1];
    p1p2.ptr[2] = p2.ptr[2] - p1.ptr[2];

    p2p4.ptr[0] = p4.ptr[0] - p2.ptr[0];
    p2p4.ptr[1] = p4.ptr[1] - p2.ptr[1];
    p2p4.ptr[2] = p4.ptr[2] - p2.ptr[2];
    p2p3.ptr[0] = p3.ptr[0] - p2.ptr[0];
    p2p3.ptr[1] = p3.ptr[1] - p2.ptr[1];
    p2p3.ptr[2] = p3.ptr[2] - p2.ptr[2];

    p3p4.ptr[0] = p4.ptr[0] - p3.ptr[0];
    p3p4.ptr[1] = p4.ptr[1] - p3.ptr[1];
    p3p4.ptr[2] = p4.ptr[2] - p3.ptr[2];
    p3p1.ptr[0] = p1.ptr[0] - p3.ptr[0];
    p3p1.ptr[1] = p1.ptr[1] - p3.ptr[1];
    p3p1.ptr[2] = p1.ptr[2] - p3.ptr[2];

    cross1.ptr[0] = p1p4.ptr[1] * p1p2.ptr[2] - p1p4.ptr[2] * p1p2.ptr[1];
    cross1.ptr[1] = p1p4.ptr[2] * p1p2.ptr[0] - p1p4.ptr[0] * p1p2.ptr[2];
    cross1.ptr[2] = p1p4.ptr[0] * p1p2.ptr[1] - p1p4.ptr[1] * p1p2.ptr[0];

    cross2.ptr[0] = p2p4.ptr[1] * p2p3.ptr[2] - p2p4.ptr[2] * p2p3.ptr[1];
    cross2.ptr[1] = p2p4.ptr[2] * p2p3.ptr[0] - p2p4.ptr[0] * p2p3.ptr[2];
    cross2.ptr[2] = p2p4.ptr[0] * p2p3.ptr[1] - p2p4.ptr[1] * p2p3.ptr[0];

    cross3.ptr[0] = p3p4.ptr[1] * p3p1.ptr[2] - p3p4.ptr[2] * p3p1.ptr[1];
    cross3.ptr[1] = p3p4.ptr[2] * p3p1.ptr[0] - p3p4.ptr[0] * p3p1.ptr[2];
    cross3.ptr[2] = p3p4.ptr[0] * p3p1.ptr[1] - p3p4.ptr[1] * p3p1.ptr[0];

    norm1 = sqrt(cross1.ptr[0] * cross1.ptr[0] + cross1.ptr[1] * cross1.ptr[1] + cross1.ptr[2] * cross1.ptr[2]);
    norm2 = sqrt(cross2.ptr[0] * cross2.ptr[0] + cross2.ptr[1] * cross2.ptr[1] + cross2.ptr[2] * cross2.ptr[2]);
    norm3 = sqrt(cross3.ptr[0] * cross3.ptr[0] + cross3.ptr[1] * cross3.ptr[1] + cross3.ptr[2] * cross3.ptr[2]);
    if((norm1 == 0) || (norm2 == 0) || (norm3 == 0))  return true; //if p4 on one edge, it is also within the triangle.
    else {
      theta12 = (cross1.ptr[0] * cross2.ptr[0] + cross1.ptr[1] * cross2.ptr[1] + cross1.ptr[2] * cross2.ptr[2]) / (norm1 * norm2);
      theta23 = (cross2.ptr[0] * cross3.ptr[0] + cross2.ptr[1] * cross3.ptr[1] + cross2.ptr[2] * cross3.ptr[2]) / (norm2 * norm3);
      theta31 = (cross3.ptr[0] * cross1.ptr[0] + cross3.ptr[1] * cross1.ptr[1] + cross3.ptr[2] * cross1.ptr[2]) / (norm3 * norm1);
      if(((-1.0e-12 < (theta12 - 1.0)) && ((theta12 - 1.0) < 1.0e-12)) && ((-1.0e-12 < (theta23 - 1.0)) && ((theta23 - 1.0) < 1.0e-12)) && ((-1.0e-12 < (theta31 - 1.0)) && ((theta31 - 1.0) < 1.0e-12))) return true;
      else return false;
    }
  }
  //check whether point p4 is on the plane p1-p2-p3 or not
  bool on_plane(Array p4, Array p1, Array p2, Array p3) {
    Array area(4, 4);
    area.ptr[12] = area.ptr[13] = area.ptr[14] = area.ptr[15] = 1.0;
    area.ptr[0] = p1.ptr[0];
    area.ptr[4] = p1.ptr[1];
    area.ptr[8] = p1.ptr[2];
    area.ptr[1] = p2.ptr[0];
    area.ptr[5] = p2.ptr[1];
    area.ptr[9] = p2.ptr[2];
    area.ptr[2] = p3.ptr[0];
    area.ptr[6] = p3.ptr[1];
    area.ptr[10] = p3.ptr[2];
    area.ptr[3] = p4.ptr[0];
    area.ptr[7] = p4.ptr[1];
    area.ptr[11] = p4.ptr[2];
    if((-1.0e-12 <= area.det()) && (area.det() <= 1.0e-12)) return true;
    else return false;
  }
  //check whether point p3 and p4 are to the same side of line p1-p2
  bool same_side(Array p3, Array p4, Array p1, Array p2) {
    Array p1p2(1, 3), p1p3(1, 3), p1p4(1, 3), cross1(1, 3), cross2(1, 3);
    double norm1, norm2, theta;
    p1p2.ptr[0] = p2.ptr[0] - p1.ptr[0];
    p1p2.ptr[1] = p2.ptr[1] - p1.ptr[1];
    p1p2.ptr[2] = p2.ptr[2] - p1.ptr[2];

    p1p3.ptr[0] = p3.ptr[0] - p1.ptr[0];
    p1p3.ptr[1] = p3.ptr[1] - p1.ptr[1];
    p1p3.ptr[2] = p3.ptr[2] - p1.ptr[2];

    p1p4.ptr[0] = p4.ptr[0] - p1.ptr[0];
    p1p4.ptr[1] = p4.ptr[1] - p1.ptr[1];
    p1p4.ptr[2] = p4.ptr[2] - p1.ptr[2];

    cross1.ptr[0] = p1p2.ptr[1] * p1p4.ptr[2] - p1p2.ptr[2] * p1p4.ptr[1];
    cross1.ptr[1] = p1p2.ptr[2] * p1p4.ptr[0] - p1p2.ptr[0] * p1p4.ptr[2];
    cross1.ptr[2] = p1p2.ptr[0] * p1p4.ptr[1] - p1p2.ptr[1] * p1p4.ptr[0];

    cross2.ptr[0] = p1p2.ptr[1] * p1p3.ptr[2] - p1p2.ptr[2] * p1p3.ptr[1];
    cross2.ptr[1] = p1p2.ptr[2] * p1p3.ptr[0] - p1p2.ptr[0] * p1p3.ptr[2];
    cross2.ptr[2] = p1p2.ptr[0] * p1p3.ptr[1] - p1p2.ptr[1] * p1p3.ptr[0];

    norm1 = sqrt(cross1.ptr[0] * cross1.ptr[0] + cross1.ptr[1] * cross1.ptr[1] + cross1.ptr[2] * cross1.ptr[2]);
    norm2 = sqrt(cross2.ptr[0] * cross2.ptr[0] + cross2.ptr[1] * cross2.ptr[1] + cross2.ptr[2] * cross2.ptr[2]);
    if((norm1 == 0.0) || (norm2 == 0.0)) return true; // points on the line
    else {
      theta = (cross1.ptr[0] * cross2.ptr[0] + cross1.ptr[1] * cross2.ptr[1] + cross1.ptr[2] * cross2.ptr[2]) / (norm1 * norm2);

      if((-1.0e-12 < (theta - 1.0)) && ((theta - 1.0) < 1.0e-12))
        return true;
      else
        return false;
    }

  }
  //check whether triangle p3-p4-p5 is cut by line p1-p2, cut return true
  bool cut(Array p3, Array p4, Array p5, Array p1, Array p2) {

    if((!on_plane(p1, p3, p4, p5)) || (!on_plane(p2, p3, p4, p5))) return false; // line p1-p2 is not on the plane p3-p4-p5
    else if((!same_side(p3, p4, p1, p2)) || (!same_side(p3, p5, p1, p2)) || (!same_side(p4, p5, p1, p2))) return true;
    else return false;
  }

};


#endif
//*************************************** END of Array Class********************************

