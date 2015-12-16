#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

/// What is being tested:
#include "casm/casm_io/jsonParser.hh"

/// What is being used to test it:
#include <boost/filesystem.hpp>

using namespace CASM;

std::string json_str = 
R"({
  "int" : 34,
  "number" : 4.0023,
  "string" : "hello",
  "bool_true" : true,
  "bool_false" : false,
  "object" : {
    "int" : 34,
    "number" : 4.0023,
    "string" : "hello",
    "bool_true" : true,
    "bool_false" : false
  },
  "uniform_array" : [1, 2, 3, 4],
  "mixed_array" : [
    "hello",
    34,
    4.0023,
    {"int" : 34, "number" : 4.0023}
  ]
})";

BOOST_AUTO_TEST_SUITE(jsonParserTest)

BOOST_AUTO_TEST_CASE(Basic) {

  jsonParser json = jsonParser::parse(json_str);
  
  BOOST_CHECK_EQUAL(true, json.is_obj());
  
  // test jsonParser::get<T>()
  int i;
  from_json(i, json["int"]);
  BOOST_CHECK_EQUAL(i, 34);
  BOOST_CHECK_EQUAL(34, json["int"].get<int>());
  
  double d;
  from_json(d, json["number"]);
  BOOST_CHECK_EQUAL(d, 4.0023);
  BOOST_CHECK_EQUAL(4.0023, json["number"].get<double>());
}

BOOST_AUTO_TEST_SUITE_END()
 
#include "../src/casm/casm_io/jsonParser.cc"
