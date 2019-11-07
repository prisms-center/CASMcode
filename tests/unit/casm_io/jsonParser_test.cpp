#include "gtest/gtest.h"

/// What is being tested:
#include "casm/casm_io/json_io/jsonParser.hh"

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

TEST(jsonParserTest, Basic) {

  jsonParser json = jsonParser::parse(json_str);

  ASSERT_EQ(true, json.is_obj());

  // test jsonParser::get<T>()
  int i;
  from_json(i, json["int"]);
  ASSERT_EQ(i, 34);
  ASSERT_EQ(34, json["int"].get<int>());

  double d;
  from_json(d, json["number"]);
  ASSERT_EQ(d, 4.0023);
  ASSERT_EQ(4.0023, json["number"].get<double>());
}

template<typename T>
void test_at(T &json) {
  ASSERT_EQ(json.at("int").template get<int>(), 34);
  ASSERT_THROW(json.at("mistake"), std::invalid_argument);

  ASSERT_EQ(json.at(fs::path("object") / "int").template get<int>(), 34);
  ASSERT_THROW(json.at(fs::path("object") / "mistake"), std::invalid_argument);

  ASSERT_EQ(json.at(fs::path("mixed_array") / "1").template get<int>(), 34);
  ASSERT_THROW(json.at(fs::path("mixed_array") / "10"), std::invalid_argument);

  ASSERT_EQ(json["uniform_array"].at(0).template get<int>(), 1);
  ASSERT_THROW(json["object"].at(0), std::invalid_argument);
  ASSERT_THROW(json["uniform_array"].at(-1), std::out_of_range);
  ASSERT_THROW(json["uniform_array"].at(4), std::out_of_range);
  ASSERT_THROW(json["uniform_array"].at(100), std::out_of_range);
}

TEST(jsonParserTest, At) {

  jsonParser json = jsonParser::parse(json_str);
  test_at(json);
}

TEST(jsonParserTest, ConstAt) {

  const jsonParser json = jsonParser::parse(json_str);
  test_at(json);
}

template<typename T>
void test_find_at(T &json) {

  ASSERT_EQ(&json == &*json.find_at(fs::path()), true);
  ASSERT_EQ(&json == &*json.find_at(""), true);

  {
    auto it = json.find_at(fs::path());
    ASSERT_EQ(json.end() == it, false);
    it++;
    ASSERT_EQ(json.end() == it, true);
  }


  ASSERT_EQ(json.find_at("int")->template get<int>(), 34);
  ASSERT_EQ(json.find_at("mistake") == json.end(), true);

  ASSERT_EQ(json.find_at(fs::path("object") / "int")->template get<int>(), 34);
  ASSERT_EQ(json.find_at(fs::path("object") / "mistake") == json.end(), true);

  ASSERT_EQ(json.find_at(fs::path("mixed_array") / "1")->template get<int>(), 34);
  ASSERT_EQ(json.find_at(fs::path("mixed_array") / "10") == json.end(), true);
}

TEST(jsonParserTest, FindAt) {

  jsonParser json = jsonParser::parse(json_str);
  test_find_at(json);
}

TEST(jsonParserTest, ConstFindAt) {

  const jsonParser json = jsonParser::parse(json_str);

  test_find_at(json);
}

TEST(jsonParserTest, Get) {

  const jsonParser json = jsonParser::parse(json_str);

  ASSERT_EQ(json["int"].get<int>(), 34);
  ASSERT_THROW(json["int"].get<std::string>(), std::runtime_error);
}

TEST(jsonParserTest, ArrayExtraTrailingComma) {

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

  ASSERT_THROW(jsonParser::parse(json_extra_trailing_comma), std::runtime_error);

}

TEST(jsonParserTest, ArrayMissingComma) {

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

  ASSERT_THROW(jsonParser::parse(json_missing_comma), std::runtime_error);

}

TEST(jsonParserTest, FindDiffTest) {

  jsonParser A = jsonParser::parse(json_str);
  jsonParser B {A};

  B["object"]["number"] = B["object"]["number"].get<double>() + 1e-8;

  ASSERT_TRUE(A != B);
  ASSERT_TRUE(A.almost_equal(B, 1e-5));

  fs::path diff_point = find_diff(A, B);
  ASSERT_EQ(diff_point, fs::path("object/number"));
  ASSERT_EQ(A.at(diff_point), 4.0023);

  diff_point = find_diff(A, B, 1e-5);
  ASSERT_TRUE(diff_point.empty());

}

