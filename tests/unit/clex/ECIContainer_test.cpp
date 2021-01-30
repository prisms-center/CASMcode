#include "casm/clex/ECIContainer.hh"

#include "Common.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clex/io/json/ECIContainer_json_io.hh"
#include "gtest/gtest.h"

using namespace CASM;

class ECIContainerJSONTest : public testing::Test {
 protected:
  ECIContainerJSONTest()
      : eci_json(test::data_file("clex", "FCC_ternary_bspecs_ex0_basis.json")) {

  }

  void add_simple_eci() {
    // add 2 eci
    eci_json["orbits"][0]["eci"] = 1.0;
    eci_json["orbits"][3]["eci"] = 4.0;
  }

  jsonParser eci_json;
};

TEST_F(ECIContainerJSONTest, FromJSON) {
  // setup
  add_simple_eci();

  // read eci JSON
  InputParser<ECIContainer> eci_parser{eci_json};

  // assert
  EXPECT_TRUE(eci_parser.valid());

  ECIContainer const &eci = *eci_parser.value;
  EXPECT_EQ(eci.size(), 2);

  EXPECT_EQ(eci.index()[0], 0);
  EXPECT_TRUE(almost_equal(eci.value()[0], 1.0));

  EXPECT_EQ(eci.index()[1], 3);
  EXPECT_TRUE(almost_equal(eci.value()[1], 4.0));
}
