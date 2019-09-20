#include "gtest/gtest.h"

/// Dependencies
#include "casm/container/Array.hh"

/// What is being tested:
#include "casm/container/Counter.hh"

/// What is being used to test it:
#include<vector>

using namespace CASM;

TEST(CounterTest, ArrayIntCounter) {
  {
    int size = 3;
    int count = 0;
    Array<int> initial(size, 0);
    Array<int> final(size, 1);
    Array<int> increment(size, 1);
    Counter<Array<int> > counter(initial, final, increment);

    EXPECT_EQ(counter.size(), size);

    do {
      count++;
    }
    while(++counter);

    EXPECT_EQ(count, 8);
  }

  {
    int size = 3;
    int count = 0;
    Array<int> initial(size, 0);
    Array<int> final(size, -2);
    Array<int> increment(size, -1);
    Counter<Array<int> > counter(initial, final, increment);

    do {
      count++;
    }
    while(++counter);

    EXPECT_EQ(count, 27);
  }

  {
    int size = 3;
    int count = 0;
    Array<int> initial(size, 0);
    Array<int> final(size, 2);
    Array<int> increment(size, 1);

    final[1] = -2;
    increment[1] = -1;

    Counter<Array<int> > counter(initial, final, increment);

    do {
      count++;
    }
    while(++counter);

    EXPECT_EQ(count, 27);
  }
}

TEST(CounterTest, stdVectorDoubleCounter) {
  {
    int size = 3;
    int count = 0;
    Array<double> initial(size, 0.0);
    Array<double> final(size, 1.0);
    Array<double> increment(size, 0.4);
    Counter<Array<double> > counter(initial, final, increment);

    EXPECT_EQ(counter.size(), size);

    do {
      count++;
    }
    while(++counter);

    EXPECT_EQ(count, 27);
  }

  {
    int size = 3;
    int count = 0;
    Array<double> initial(size, 0.0);
    Array<double> final(size, -1.0);
    Array<double> increment(size, -0.4);
    Counter<Array<double> > counter(initial, final, increment);

    do {
      count++;
    }
    while(++counter);

    EXPECT_EQ(count, 27);
  }

  {
    int size = 3;
    int count = 0;
    Array<double> initial(size, 0.0);
    Array<double> final(size, 1.0);
    Array<double> increment(size, 0.4);

    final[1] = -1.0;
    increment[1] = -0.4;

    Counter<Array<double> > counter(initial, final, increment);

    do {
      count++;
    }
    while(++counter);

    EXPECT_EQ(count, 27);
  }
}
