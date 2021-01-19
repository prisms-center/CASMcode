#include "casm/misc/CASM_math.hh"

#include <map>

#include "casm/global/eigen.hh"

namespace CASM {
//*******************************************************************************************

int round(double val) {
  return int(val < 0 ? floor(val + 0.5) : ceil(val - 0.5));
};

//*******************************************************************************************

double ran0(int &idum) {
  int IA = 16807;
  int IM = 2147483647;
  int IQ = 127773;
  int IR = 2836;
  int MASK = 123459876;
  double AM = 1.0 / IM;

  // minimal random number generator of Park and Miller
  // returns uniform random deviate between 0.0 and 1.0
  // set or rest idum to any  integer value (except the
  // unlikely value MASK) to initialize the sequence: idum must
  // not be altered between calls for successive deviates
  // in a sequence

  int k;
  idum = idum ^ MASK;  // XOR the two integers
  k = idum / IQ;
  idum = IA * (idum - k * IQ) - IR * k;
  if (idum < 0) {
    idum = idum + IM;
  }
  double ran = AM * idum;
  idum = idum ^ MASK;
  return ran;
}

//*******************************************************************************************

int dl_string_dist(const std::string &a, const std::string &b) {
  // "infinite" distance is just the max possible distance
  int max_val = a.size() + b.size();

  // make and initialize the character array indices
  std::map<char, int> DA;

  // make the distance matrix H
  Eigen::MatrixXi H(a.size() + 2, b.size() + 2);

  // initialize the left and top edges of H
  H(0, 0) = max_val;
  for (int i = 0; i <= a.size(); ++i) {
    DA[a[i]] = 0;
    H(i + 1, 0) = max_val;
    H(i + 1, 1) = i;
  }
  for (int j = 0; j <= b.size(); ++j) {
    DA[b[j]] = 0;
    H(0, j + 1) = max_val;
    H(1, j + 1) = j;
  }

  // fill in the distance matrix H
  // look at each character in a
  for (int i = 1; i <= a.size(); ++i) {
    int DB = 0;
    // look at each character in b
    for (int j = 1; j <= b.size(); ++j) {
      int i1 = DA[b[j - 1]];
      int j1 = DB;
      int cost;
      if (a[i - 1] == b[j - 1]) {
        cost = 0;
        DB = j;
      } else
        cost = 1;
      H(i + 1, j + 1) = min(min(H(i, j) + cost,    // substitution
                                H(i + 1, j) + 1),  // insertion
                            min(H(i, j + 1) + 1,   // deletion
                                H(i1, j1) + (i - i1 - 1) + 1 + (j - j1 - 1)));
    }
    DA[a[i - 1]] = i;
  }
  return H(a.size() + 1, b.size() + 1);
}

//*******************************************************************************************
/// Find greatest common factor
int gcf(int i1, int i2) {
  i1 = std::abs(i1);
  i2 = std::abs(i2);
  while (i1 != i2 && i1 != 1 && i2 != 1) {
    if (i1 < i2) {
      i2 -= i1;
    } else {
      i1 -= i2;
    }
  }
  if (i1 == 1) {
    return i1;
  } else {
    return i2;
  }
}

//*******************************************************************************************
/// Find least common multiple
int lcm(int i1, int i2) { return std::abs(i1 * (i2 / gcf(i1, i2))); }

//*******************************************************************************************
// This evaluates Gaussians using the formula:
// f(x) = a*e^(-(x-b)^2/(c^2))
//*******************************************************************************************

double gaussian(double a, double x, double b, double c) {
  return a * exp(-((x - b) * (x - b)) / (c * c));
}

//*******************************************************************************************
// This calculates Gaussian moments given by the integral:
// m = \int_{\infty}^{\infty} dx
// x^pow*exp[-x^2/(2*sigma^2)]/(\sqrt(2*\pi)*sigma)
//*******************************************************************************************

double gaussian_moment(int expon, double sigma) {
  if (expon % 2) return 0.0;

  double m = pow(sigma, expon);

  expon -= 1;
  while (expon - 2 > 0) {
    m *= double(expon);
    expon -= 2;
  }
  return m;
}

//*******************************************************************************************
// This calculates Gaussian moments given by the integral:
// m = \int_{\infty}^{\infty} dx
// x^pow*exp[-(x-x0)^2/(2*sigma^2)]/(\sqrt(2*\pi)*sigma)
//*******************************************************************************************

double gaussian_moment(int expon, double sigma, double x0) {
  double m = 0;
  for (int i = 0; i <= expon; i++) {
    m += nchoosek(expon, i) * gaussian_moment(i, sigma) * pow(x0, expon - i);
  }

  return m;
}

//*******************************************************************************************
// finds rational number that approximates 'val' to within 'tol' --
//      --  almost_int(double(denominator)*val, tol) == true OR
//      almost_int(val/double(numerator), tol) == true OR
//      --  denominator is always positive
//      --  sign(numerator)=sign(val)
//      --  if(almost_zero(val,tol)) --> numerator = 0, denominator = 1
void nearest_rational_number(double val, long &numerator, long &denominator,
                             double tol) {
  if (almost_zero(val, tol)) {
    numerator = 0;
    denominator = 1;
    return;
  }
  long sgn(val < 0 ? -1 : 1);
  val = std::abs(val);
  double tdenom, tnum;
  long lim(max(long(100), long(1 / (10 * tol))));
  for (long i = 1; i < lim + 1; i++) {
    tdenom = double(i) / val;
    tnum = val / double(i);
    if (tdenom > 1 && almost_zero(tdenom - round(tdenom), tol)) {
      numerator = sgn * i;
      denominator = round(tdenom);
      return;
    } else if (tnum > 1 && almost_zero(tnum - round(tnum), tol)) {
      denominator = i;
      numerator = sgn * round(tnum);
      return;
    }
  }
}

//*******************************************************************************************
/* Finds best irrational number approximation of double 'val' and
 * returns tex-formated string that contains irrational approximation
 * searches numbers of the form (x/y)^(1/z), where x and y range from 1 to 'lim'
 * z ranges from 1 to 'max_pow'
 */
//*******************************************************************************************

std::string irrational_to_tex_string(double val, int lim, int max_pow) {
  std::stringstream tstr;
  if (almost_zero(round(val) - val)) {
    tstr << round(val);
    return tstr.str();
  }
  if (val < 0) {
    tstr << '-';
    val = std::abs(val);
  }
  double tval(val), tdenom, tnum;
  int idenom, inum;
  for (int ipow = 1; ipow < max_pow + 1; ipow++) {
    for (int i = 1; i < lim + 1; i++) {
      tdenom = double(i) / tval;
      tnum = tval / double(i);
      if (tdenom > 1 && almost_zero(std::abs(tdenom - round(tdenom)))) {
        inum = i;
        idenom = round(tdenom);
      } else if (tnum > 1 && almost_zero(std::abs(tnum - round(tnum)))) {
        idenom = i;
        inum = round(tnum);
      } else {
        continue;
      }

      if (ipow == 1) {
        tstr << inum << '/' << idenom;
        return tstr.str();
      }
      if (ipow == 2) {
        tstr << "\\sqrt{" << inum;
        if (idenom != 1) tstr << '/' << idenom;
        tstr << '}';
        return tstr.str();
      } else {
        tstr << '(' << inum;
        if (idenom != 1) tstr << '/' << idenom;
        tstr << ")^{1/" << ipow << "}";
        return tstr.str();
      }
    }
    tval *= val;
  }
  tstr << val;
  return tstr.str();
}

//*******************************************************************************************
std::string to_sequential_string(Index i, Index max_i,
                                 char prepend_char /*='0'*/) {
  max_i = max(i, max_i);
  Index length = 1;
  while (max_i /= 10) length++;

  std::string tresult = std::to_string(i);

  std::string result(length - tresult.size(), prepend_char);
  return result.append(tresult);
}

//*******************************************************************************************
// John G 010413
int mod(int a, int b) {
  if (b < 0) return mod(-a, -b);

  int ret = a % b;
  if (ret < 0) ret += b;
  return ret;
}

//*******************************************************************************************

double cuberoot(double number) {
  if (number < 0) {
    return -pow(-number, 1.0 / 3);
  }

  else {
    return pow(number, 1.0 / 3);
  }
}

}  // namespace CASM
