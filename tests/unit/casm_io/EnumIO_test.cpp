#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

/// What is being tested:
#include "casm/casm_io/EnumIO.hh"
#include "casm/CASM_global_enum.hh"

/// What is being used to test it:
#include <boost/filesystem.hpp>
#include "casm/casm_io/jsonParser.hh"


using namespace CASM;

template<typename EnumType>
void enum_io_test(EnumType TEST) {

  // 1. std::string to_string(ENUM val)
  // 2. ENUM from_string(const std::string &val)
  // 3. std::string multiline_help<ENUM>();
  // 4. std::string singleline_help<ENUM>();
  // 5. std::string help<ENUM>();
  // 6. std::ostream &operator<<(std::ostream &sout, const ENUM& val); \
  // 7. std::istream &operator>>(std::istream &sin, ENUM& val); \
  // 8. jsonParser &to_json(const ENUM &val, jsonParser &json); \
  // 9. void from_json(ENUM& val, const jsonParser& json);
  // 10. const std::string traits<ENUM>::name
  // 11. static const std::multimap<ENUM, std::vector<std::string> > traits<ENUM>::strval

  const std::vector<std::string> &opt = traits<EnumType>::strval.find(TEST)->second;
  std::string to_str = opt.front();

  // 1.
  BOOST_CHECK_EQUAL(to_string(TEST), to_str);

  // 2.
  BOOST_CHECK_THROW(from_string<EnumType>("mistake"), std::invalid_argument);

  // 3.
  BOOST_CHECK_EQUAL(!multiline_help<EnumType>().empty(), true);

  // 4.
  BOOST_CHECK_EQUAL(!singleline_help<EnumType>().empty(), true);

  // 5.
  BOOST_CHECK_EQUAL(help<EnumType>(), multiline_help<EnumType>());

  // 6.
  std::stringstream ss;
  ss << TEST;
  BOOST_CHECK_EQUAL(ss.str(), to_str);

  // check reading unrecognized strings
  {
    jsonParser json;
    EnumType tmp;

    // 7.
    std::stringstream ss;
    ss << "mistake";
    ss.seekg(0);
    BOOST_CHECK_THROW((ss >> tmp), std::invalid_argument);

    // 8. & 9. to_json  & from_json
    json.put_obj();
    json[traits<EnumType>::name] = "mistake";
    BOOST_CHECK_THROW(from_json(tmp, json[traits<EnumType>::name]), std::invalid_argument);
    BOOST_CHECK_THROW((tmp = from_json<EnumType>(json[traits<EnumType>::name])), std::invalid_argument);
    BOOST_CHECK_THROW((tmp = json[traits<EnumType>::name].template get<EnumType>()), std::invalid_argument);
  }

  // try reading all recognized strings
  for(const std::string &val : opt) {

    jsonParser json;
    EnumType tmp;

    // 2.
    BOOST_CHECK_EQUAL(TEST, from_string<EnumType>(val));

    // 7.
    std::stringstream ss;
    ss << val;
    ss.seekg(0);
    ss >> tmp;
    BOOST_CHECK_EQUAL(TEST, tmp);

    // 8. & 9. to_json  & from_json
    json.put_obj();
    json[traits<EnumType>::name] = val;
    from_json(tmp, json[traits<EnumType>::name]);
    BOOST_CHECK_EQUAL(TEST, tmp);
    BOOST_CHECK_EQUAL(TEST, from_json<EnumType>(json[traits<EnumType>::name]));
    BOOST_CHECK_EQUAL(TEST, json[traits<EnumType>::name].template get<EnumType>());
  }

}

BOOST_AUTO_TEST_SUITE(EnumIOTest)

// From CASM_global_enum

BOOST_AUTO_TEST_CASE(CoordTypeTest) {
  enum_io_test(COORD_TYPE::FRAC);
  enum_io_test(COORD_TYPE::CART);
  enum_io_test(COORD_TYPE::INTEGRAL);
  enum_io_test(COORD_TYPE::COORD_DEFAULT);
}

BOOST_AUTO_TEST_CASE(PeriodicityTypeTest) {
  enum_io_test(PERIODICITY_TYPE::PERIODIC);
  enum_io_test(PERIODICITY_TYPE::APERIODIC);
  enum_io_test(PERIODICITY_TYPE::LOCAL);
  BOOST_CHECK_EQUAL(PERIODICITY_TYPE::APERIODIC, PERIODICITY_TYPE::LOCAL);
  enum_io_test(PERIODICITY_TYPE::PERIODICITY_DEFAULT);
}

BOOST_AUTO_TEST_CASE(EquivalenceTypeTest) {
  enum_io_test(EQUIVALENCE_TYPE::PRIM);
  enum_io_test(EQUIVALENCE_TYPE::SCEL);
  enum_io_test(EQUIVALENCE_TYPE::CONFIG);
}

BOOST_AUTO_TEST_CASE(CellTypeTest) {
  enum_io_test(CELL_TYPE::PRIM);
  enum_io_test(CELL_TYPE::SCEL);
}

BOOST_AUTO_TEST_CASE(OnErrorTest) {
  enum_io_test(OnError::THROW);
  enum_io_test(OnError::WARN);
  enum_io_test(OnError::CONTINUE);
}

BOOST_AUTO_TEST_SUITE_END()
