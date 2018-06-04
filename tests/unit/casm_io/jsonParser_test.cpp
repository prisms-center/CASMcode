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

template<typename T>
void test_at(T &json) {
  BOOST_CHECK_EQUAL(json.at("int").template get<int>(), 34);
  BOOST_CHECK_THROW(json.at("mistake"), std::invalid_argument);

  BOOST_CHECK_EQUAL(json.at(fs::path("object") / "int").template get<int>(), 34);
  BOOST_CHECK_THROW(json.at(fs::path("object") / "mistake"), std::invalid_argument);

  BOOST_CHECK_EQUAL(json.at(fs::path("mixed_array") / "1").template get<int>(), 34);
  BOOST_CHECK_THROW(json.at(fs::path("mixed_array") / "10"), std::invalid_argument);

  BOOST_CHECK_EQUAL(json["uniform_array"].at(0).template get<int>(), 1);
  BOOST_CHECK_THROW(json["object"].at(0), std::invalid_argument);
  BOOST_CHECK_THROW(json["uniform_array"].at(-1), std::out_of_range);
  BOOST_CHECK_THROW(json["uniform_array"].at(4), std::out_of_range);
  BOOST_CHECK_THROW(json["uniform_array"].at(100), std::out_of_range);
}

BOOST_AUTO_TEST_CASE(At) {

  jsonParser json = jsonParser::parse(json_str);
  test_at(json);
}

BOOST_AUTO_TEST_CASE(ConstAt) {

  const jsonParser json = jsonParser::parse(json_str);
  test_at(json);
}

template<typename T>
void test_find_at(T &json) {

  BOOST_CHECK_EQUAL(json.end() == json.find_at(fs::path()), true);
  BOOST_CHECK_EQUAL(json.end() == json.find_at(""), true);

  BOOST_CHECK_EQUAL(json.find_at("int")->template get<int>(), 34);
  BOOST_CHECK_EQUAL(json.find_at("mistake") == json.end(), true);

  BOOST_CHECK_EQUAL(json.find_at(fs::path("object") / "int")->template get<int>(), 34);
  BOOST_CHECK_EQUAL(json.find_at(fs::path("object") / "mistake") == json.end(), true);

  BOOST_CHECK_EQUAL(json.find_at(fs::path("mixed_array") / "1")->template get<int>(), 34);
  BOOST_CHECK_EQUAL(json.find_at(fs::path("mixed_array") / "10") == json.end(), true);
}

BOOST_AUTO_TEST_CASE(FindAt) {

  jsonParser json = jsonParser::parse(json_str);
  test_find_at(json);
}

BOOST_AUTO_TEST_CASE(ConstFindAt) {

  const jsonParser json = jsonParser::parse(json_str);

  test_find_at(json);
}

BOOST_AUTO_TEST_CASE(Get) {

  const jsonParser json = jsonParser::parse(json_str);

  BOOST_CHECK_EQUAL(json["int"].get<int>(), 34);
  BOOST_CHECK_THROW(json["int"].get<std::string>(), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(ArrayExtraTrailingComma) {

  std::string json_extra_trailing_comma =
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
"uniform_array" : [1, 2, 3, 4,],
"mixed_array" : [
"hello",
34,
4.0023,
{"int" : 34, "number" : 4.0023}
]
})";

  BOOST_CHECK_THROW(jsonParser::parse(json_extra_trailing_comma), std::runtime_error);

}

BOOST_AUTO_TEST_CASE(ArrayMissingComma) {

  std::string json_missing_comma =
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
"uniform_array" : [1, 2 3, 4],
"mixed_array" : [
"hello",
34,
4.0023,
{"int" : 34, "number" : 4.0023}
]
})";

  BOOST_CHECK_THROW(jsonParser::parse(json_missing_comma), std::runtime_error);

}

BOOST_AUTO_TEST_CASE(FindDiffTest) {

  jsonParser A = jsonParser::parse(json_str);
  jsonParser B {A};

  B["object"]["number"] = B["object"]["number"].get<double>() + 1e-8;

  BOOST_CHECK(A != B);
  BOOST_CHECK(A.almost_equal(B, 1e-5));

  fs::path diff_point = find_diff(A, B);
  BOOST_CHECK_EQUAL(diff_point, fs::path("object/number"));
  BOOST_CHECK_EQUAL(A.at(diff_point), 4.0023);

  diff_point = find_diff(A, B, 1e-5);
  BOOST_CHECK(diff_point.empty());

}

BOOST_AUTO_TEST_SUITE_END()

//#include "../src/casm/casm_io/jsonParser.cc"
