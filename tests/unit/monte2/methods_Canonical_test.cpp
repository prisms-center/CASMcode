#include "casm/crystallography/Structure.hh"
#include "casm/external/MersenneTwister/MersenneTwister.h"
#include "casm/monte2/Conversions.hh"
#include "casm/monte2/events/OccCandidate.hh"
#include "casm/monte2/events/OccLocation.hh"
#include "crystallography/TestStructures.hh"
#include "gtest/gtest.h"

using namespace CASM;

TEST(MethodsCanonicalTest, Test1) { EXPECT_EQ(1, 1); }
