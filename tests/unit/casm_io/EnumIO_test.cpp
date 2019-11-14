#include "gtest/gtest.h"

/// What is being tested:
#include "casm/casm_io/enum/json_io.hh"
#include "casm/casm_io/enum/stream_io.hh"
#include "casm/global/enum.hh"
#include "casm/global/enum/io_traits.hh"
#include "casm/global/enum/stream_io.hh"
#include "casm/global/enum/json_io.hh"

/// What is being used to test it:
#include <boost/filesystem.hpp>
#include "casm/casm_io/json/jsonParser.hh"


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
  ASSERT_EQ(to_string(TEST), to_str);

  // 2.
  ASSERT_THROW(from_string<EnumType>("mistake"), std::invalid_argument);

  // 3.
  ASSERT_EQ(!multiline_help<EnumType>().empty(), true);

  // 4.
  ASSERT_EQ(!singleline_help<EnumType>().empty(), true);

  // 5.
  ASSERT_EQ(help<EnumType>(), multiline_help<EnumType>());

  // 6.
  std::stringstream ss;
  ss << TEST;
  ASSERT_EQ(ss.str(), to_str);

  // check reading unrecognized strings
  {
    jsonParser json;
    EnumType tmp;

    // 7.
    std::stringstream ss;
    ss << "mistake";
    ss.seekg(0);
    ASSERT_THROW((ss >> tmp), std::invalid_argument);

    // 8. & 9. to_json  & from_json
    json.put_obj();
    json[traits<EnumType>::name] = "mistake";
    ASSERT_THROW(from_json(tmp, json[traits<EnumType>::name]), std::invalid_argument);
    ASSERT_THROW((tmp = from_json<EnumType>(json[traits<EnumType>::name])), std::invalid_argument);
    ASSERT_THROW((tmp = json[traits<EnumType>::name].template get<EnumType>()), std::invalid_argument);
  }

  // try reading all recognized strings
  for(const std::string &val : opt) {

    jsonParser json;
    EnumType tmp;

    // 2.
    ASSERT_EQ(TEST, from_string<EnumType>(val));

    // 7.
    std::stringstream ss;
    ss << val;
    ss.seekg(0);
    ss >> tmp;
    ASSERT_EQ(TEST, tmp);

    // 8. & 9. to_json  & from_json
    json.put_obj();
    json[traits<EnumType>::name] = val;
    from_json(tmp, json[traits<EnumType>::name]);
    ASSERT_EQ(TEST, tmp);
    ASSERT_EQ(TEST, from_json<EnumType>(json[traits<EnumType>::name]));
    ASSERT_EQ(TEST, json[traits<EnumType>::name].template get<EnumType>());
  }

}

// From CASM_global_enum

TEST(EnumIOTest, CoordTypeTest) {
  enum_io_test(COORD_TYPE::FRAC);
  enum_io_test(COORD_TYPE::CART);
  enum_io_test(COORD_TYPE::INTEGRAL);
  //enum_io_test(COORD_TYPE::COORD_DEFAULT);
}

TEST(EnumIOTest, PeriodicityTypeTest) {
  enum_io_test(PERIODICITY_TYPE::PERIODIC);
  enum_io_test(PERIODICITY_TYPE::APERIODIC);
  enum_io_test(PERIODICITY_TYPE::LOCAL);
  ASSERT_EQ(PERIODICITY_TYPE::APERIODIC, PERIODICITY_TYPE::LOCAL);
  //enum_io_test(PERIODICITY_TYPE::PERIODICITY_DEFAULT);
}

TEST(EnumIOTest, EquivalenceTypeTest) {
  enum_io_test(EQUIVALENCE_TYPE::PRIM);
  enum_io_test(EQUIVALENCE_TYPE::SCEL);
  enum_io_test(EQUIVALENCE_TYPE::CONFIG);
}

TEST(EnumIOTest, CellTypeTest) {
  enum_io_test(CELL_TYPE::PRIM);
  enum_io_test(CELL_TYPE::SCEL);
}

TEST(EnumIOTest, OnErrorTest) {
  enum_io_test(OnError::THROW);
  enum_io_test(OnError::WARN);
  enum_io_test(OnError::CONTINUE);
}
