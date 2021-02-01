#include "casm/casm_io/dataformatter/DataStream.hh"

#include "casm/misc/CASM_math.hh"

namespace CASM {
namespace DataStream_impl {

/*template<typename OutType> template<>
OutType DataStreamPromoter<OutType>::promote(OutType a){
  return a;
  }*/

// PROMOTE TO DOUBLE
template <>
double DataStreamPromoter<double>::promote(std::string a) {
  return std::stod(a);
}
//\End of Double promoters

// PROMOTE TO LONG
template <>
long DataStreamPromoter<long>::promote(double a) {
  return round(a);
}

template <>
long DataStreamPromoter<long>::promote(std::string a) {
  return std::stol(a);
}
//\End of long promoters

// PROMOTE TO CHAR
template <>
char DataStreamPromoter<char>::promote(double a) {
  throw std::runtime_error(
      std::string("No viable conversion to type 'char' from double ") +
      std::to_string(a) + "\n");
  return round(a);
}

template <>
char DataStreamPromoter<char>::promote(std::string a) {
  if (a.size()) return a[0];
  return (char)0;
}

template <>
char DataStreamPromoter<char>::promote(bool a) {
  if (toupper(a))
    return 'T';
  else
    return 'F';
}
//\End of char promoters

// PROMOTE TO BOOL
template <>
bool DataStreamPromoter<bool>::promote(double a) {
  return !almost_zero(a);
}

template <>
bool DataStreamPromoter<bool>::promote(std::string a) {
  // in-place convert 'a' to lower case
  std::transform(a.begin(), a.end(), a.begin(), tolower);

  if (a == "true") return true;
  if (a == "false") return false;
  throw std::runtime_error("No viable convertion to type 'bool' from string '" +
                           a + "\n");
}

template <>
bool DataStreamPromoter<bool>::promote(char a) {
  if (toupper(a) == 'F') return false;
  if (toupper(a) == 'T') return true;
  throw std::runtime_error("No viable convertion to type 'bool' from char '" +
                           std::string(1, a) + "\n");
}
//\End of bool promoters

// PROMOTE TO STRING
template <>
std::string DataStreamPromoter<std::string>::promote(long a) {
  return std::to_string(a);
}

template <>
std::string DataStreamPromoter<std::string>::promote(double a) {
  return std::to_string(a);
}

template <>
std::string DataStreamPromoter<std::string>::promote(char a) {
  return std::string(1, a);
}

template <>
std::string DataStreamPromoter<std::string>::promote(bool a) {
  if (a)
    return "true";
  else
    return "false";
}
//\End of String promoters

}  // namespace DataStream_impl
}  // namespace CASM
