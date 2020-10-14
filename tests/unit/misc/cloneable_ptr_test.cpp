#include "gtest/gtest.h"

/// What is being tested:
#include "casm/misc/cloneable_ptr.hh"


class ATest {

public:

  ATest() :
    ATest(0) {}

  ATest(int _val) :
    m_value(_val) {};

  ATest(const ATest &other) :
    m_value(other.m_value) {};

  ATest(ATest &&other) {
    if(this != &other) {
      m_value = other.m_value;
      other.m_value = 0;
    }
  };

  ATest &operator=(const ATest &other) {
    if(this != &other) {
      m_value = other.m_value;
    }
    return *this;

  };

  ATest &operator=(ATest &&other) {
    if(this != &other) {
      m_value = other.m_value;
      other.m_value = 0;
    }
    return *this;

  };

  int value() const {
    return m_value;
  }

private:

  int m_value;
};

class BTest {

public:

  BTest() :
    BTest(-1) {}

  BTest(int _val, bool _cloned = false) :
    m_cloned(_cloned),
    m_value(_val) {};

  BTest(const BTest &other) :
    m_value(other.m_value) {};

  virtual ~BTest() {}

  BTest(BTest &&other) {
    if(this != &other) {
      m_value = other.m_value;
      other.m_value = -1;
    }
  };

  BTest &operator=(const BTest &other) {
    if(this != &other) {
      m_value = other.m_value;
    }
    return *this;

  };

  BTest &operator=(BTest &&other) {
    if(this != &other) {
      m_value = other.m_value;
      other.m_value = -1;
    }
    return *this;

  };

  virtual int value() const {
    return m_value;
  }

  std::unique_ptr<BTest> clone() const {
    return std::unique_ptr<BTest>(this->_clone());
  }

  bool cloned() const {
    return m_cloned;
  }

protected:

  bool m_cloned;

private:

  virtual BTest *_clone() const {
    return new BTest(m_value, true);
  }

  int m_value;
};

class CTest : public BTest {

public:

  CTest() :
    CTest(-2) {}

  CTest(int _val, bool _cloned = false) :
    BTest(_val + 10, _cloned),
    m_value(_val) {};

  CTest(const CTest &other) :
    m_value(other.m_value) {};

  CTest(CTest &&other) {
    if(this != &other) {
      m_value = other.m_value;
      other.m_value = -2;
    }
  };

  CTest &operator=(const CTest &other) {
    if(this != &other) {
      m_value = other.m_value;
    }
    return *this;

  };

  CTest &operator=(CTest &&other) {
    if(this != &other) {
      m_value = other.m_value;
      other.m_value = -2;
    }
    return *this;

  };

  int value() const override {
    return m_value;
  }

  std::unique_ptr<CTest> clone() const {
    return std::unique_ptr<CTest>(this->_clone());
  }

private:

  CTest *_clone() const override {
    return new CTest(m_value, true);
  }

  int m_value;
};


TEST(CloneablePtrTest, ConstructorTest) {

  ATest A {1};
  BTest B {2};
  CTest C {3};

  EXPECT_EQ(A.value(), 1);

  // default construct
  notstd::cloneable_ptr<ATest> Aptr {};
  EXPECT_TRUE(Aptr.unique().get() == nullptr);

  // copy A
  notstd::cloneable_ptr<ATest> Aptr2 {new ATest(A)};
  EXPECT_EQ(A.value(), 1);
  EXPECT_EQ(Aptr2->value(), 1);

  // move A
  notstd::cloneable_ptr<ATest> Aptr3 {new ATest(std::move(A))};
  EXPECT_EQ(A.value(), 0);
  EXPECT_EQ(Aptr3->value(), 1);

  // copy cloneable
  notstd::cloneable_ptr<ATest> Aptr4 {Aptr2};
  EXPECT_EQ(Aptr2->value(), 1);
  EXPECT_EQ(Aptr4->value(), 1);

  // move cloneable
  notstd::cloneable_ptr<ATest> Aptr5 {std::move(Aptr2)};
  EXPECT_EQ(Aptr5->value(), 1);
  EXPECT_TRUE(Aptr2.unique().get() == nullptr);

  // copy other
  notstd::cloneable_ptr<BTest> Bptr {new CTest(C)};
  EXPECT_EQ(Bptr->value(), 3);
  EXPECT_EQ(C.value(), 3);

  // move other
  notstd::cloneable_ptr<BTest> Bptr2 {new CTest(std::move(C))};
  EXPECT_EQ(Bptr2->value(), 3);
  EXPECT_EQ(C.value(), -2);


}

TEST(CloneablePtrTest, MakeCloneableTest) {
  notstd::cloneable_ptr<ATest> Aptr = notstd::make_cloneable<ATest>();
  EXPECT_EQ(Aptr->value(), 0);

  notstd::cloneable_ptr<ATest> Aptr2 = notstd::make_cloneable<ATest>(1);
  EXPECT_EQ(Aptr2->value(), 1);

  notstd::cloneable_ptr<BTest> Bptr = notstd::make_cloneable<BTest>(2);
  EXPECT_EQ(Bptr->value(), 2);

  notstd::cloneable_ptr<BTest> Bptr2 = notstd::make_cloneable<CTest>(3);
  EXPECT_EQ(Bptr2->value(), 3);

}

TEST(CloneablePtrTest, CloneTest) {

  EXPECT_TRUE(!notstd::has_clone<ATest>::value);
  EXPECT_TRUE(notstd::has_clone<BTest>::value);
  EXPECT_TRUE(notstd::has_clone<CTest>::value);

  ATest A {1};
  BTest B {2};
  CTest C {3};

  notstd::cloneable_ptr<ATest> Aptr = notstd::clone(A);
  EXPECT_EQ(Aptr->value(), 1);

  notstd::cloneable_ptr<BTest> Bptr = notstd::clone(B);
  EXPECT_EQ(B.value(), 2);
  EXPECT_EQ(Bptr->value(), 2);
  EXPECT_EQ(Bptr->cloned(), true);

  notstd::cloneable_ptr<BTest> Bptr2 = notstd::clone(C);
  EXPECT_EQ(C.value(), 3);
  EXPECT_EQ(Bptr2->value(), 3);
  EXPECT_EQ(Bptr2->cloned(), true);

}

TEST(CloneablePtrTest, CloneFunctionTest) {

  notstd::cloneable_ptr<BTest> Bptr = notstd::make_cloneable<BTest>(2);
  notstd::cloneable_ptr<BTest> Bptr2 = notstd::make_cloneable<CTest>(3);
  EXPECT_EQ(Bptr->value(), 2);
  EXPECT_EQ(Bptr2->value(), 3);

  using std::swap;
  swap(Bptr, Bptr2);

  EXPECT_EQ(Bptr->value(), 3);
  EXPECT_EQ(Bptr2->value(), 2);

}
