#ifndef DATASTREAM_HH
#define DATASTREAM_HH
#include <vector>
#include <string>
#include <iostream>
namespace CASM {

  class DataStream {
  public:
    enum DataStreamTraits {
      none = 0,
      skipfail = (1u << 0),
      failbit = (1u << 1)
    };

    static DataStream &endl(DataStream &_strm) {
      return _strm.newline();
    };

    DataStream(DataStreamTraits _traits = none) :
      m_traits(_traits) {}

    virtual ~DataStream() {}

    virtual DataStream &operator<<(const std::string &) {
      return *this;
    }

    virtual DataStream &operator<<(long) {
      return *this;
    }

    virtual DataStream &operator<<(double) {
      return *this;
    }

    virtual DataStream &operator<<(bool) {
      return *this;
    }

    virtual DataStream &operator<<(char) {
      return *this;
    }

    /*  Some day?  This seems annoying...
    virtual DataStream &operator<<(std::complex<double>) {
      return *this;
    }
    */

    DataStream &operator<<(DataStream & (*F)(DataStream &)) {
      return F(*this);
    }

    virtual DataStream &newline() {
      return *this;
    }

    DataStream &operator<<(DataStreamTraits set_bits) {
      m_traits |= set_bits;
      return *this;
    }

    bool fail() const {
      return m_traits & failbit;
    }

    void clear_fail() {
      m_traits &= ~failbit;
    }
  protected:

    bool _skipfail() {
      return m_traits & skipfail;
    }

  private:
    int m_traits;
  };

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Count the number of POD objects streamed
  class CountDataStream : public DataStream {
  public:
    CountDataStream(bool newline_reset = true) :
      DataStream(none), m_newline_reset(newline_reset), m_count(0) {}

    DataStream &operator<<(const std::string &) override {
      return increment();
    }

    DataStream &operator<<(long) override {
      return increment();
    }

    DataStream &operator<<(double) override {
      return increment();
    }

    DataStream &operator<<(bool) override {
      return increment();
    }

    DataStream &operator<<(char) override {
      return increment();
    }

    DataStream &newline() override {
      if(m_newline_reset)
        m_count = 0;
      return *this;
    }

    int count() const {
      return m_count;
    }
  private:
    DataStream &increment() {
      m_count++;
      return *this;
    }

    bool m_newline_reset;
    int m_count;

  };

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  namespace DataStream_impl {
    template<typename OutType>
    struct DataStreamPromoter {
      template<typename InType>
      static OutType promote(InType a) {
        return static_cast<OutType>(a);
      }
    };

    template<>
    struct DataStreamPromoter<std::string> {
      template<typename InType>
      static std::string promote(InType a) {
        return static_cast<std::string>(a);
      }

    };

    // Specialized String Promoters
    template<>
    std::string DataStreamPromoter<std::string>::promote(long a);

    template<>
    std::string DataStreamPromoter<std::string>::promote(double a);

    template<>
    std::string DataStreamPromoter<std::string>::promote(char a);

    template<>
    std::string DataStreamPromoter<std::string>::promote(bool a);
    //\End specialized string promoters

    template<>
    struct DataStreamPromoter<double> {
      template<typename InType>
      static double promote(InType a) {
        return static_cast<double>(a);
      }

    };

    //Specialized double promoter
    template<>
    double DataStreamPromoter<double>::promote(std::string a);

    template<>
    struct DataStreamPromoter<long> {
      template<typename InType>
      static long promote(InType a) {
        return static_cast<long>(a);
      }

    };

    //Specialized long promotion
    template<>
    long DataStreamPromoter<long>::promote(double a);

    template<>
    long DataStreamPromoter<long>::promote(std::string a);
    //\End Specialized long promotion

    template<>
    struct DataStreamPromoter<bool> {
      template<typename InType>
      static bool promote(InType a) {
        return static_cast<bool>(a);
      }

    };

    // Specialized bool promotion
    template<>
    bool DataStreamPromoter<bool>:: promote(double a);

    template<>
    bool DataStreamPromoter<bool>:: promote(std::string a);

    template<>
    bool DataStreamPromoter<bool>:: promote(char a);
    //\End Specialized bool promotion

    template<>
    struct DataStreamPromoter<char> {
      template<typename InType>
      static char promote(InType a) {
        return static_cast<char>(a);
      }
    };

    //Specialized char promotion
    template<>
    char DataStreamPromoter<char>::promote(double a);

    template<>
    char DataStreamPromoter<char>::promote(std::string a);

    template<>
    char DataStreamPromoter<char>::promote(bool a);
    //\End specialized char promotion

  }

  template <typename T, typename Promoter = DataStream_impl::DataStreamPromoter<T> >
  class ValueDataStream : public DataStream {
  public:

    ValueDataStream(DataStreamTraits _traits = none) :
      DataStream(_traits) {}

    DataStream &operator<<(const std::string &val) override {
      m_value = Promoter::promote(val);
      return *this;
    }

    DataStream &operator<<(long val) override {
      m_value = Promoter::promote(val);
      return *this;
    }

    DataStream &operator<<(double val) override {
      m_value = Promoter::promote(val);
      return *this;
    }

    DataStream &operator<<(bool val) override {
      m_value = Promoter::promote(val);
      return *this;
    }

    DataStream &operator<<(char val) override {
      m_value = Promoter::promote(val);
      return *this;
    }

    const T &value()const {
      return m_value;
    }
  private:
    T m_value;
  };


  template <typename T, typename Promoter = DataStream_impl::DataStreamPromoter<T> >
  class VectorDataStream : public DataStream {
  public:

    VectorDataStream(DataStreamTraits _traits = none) :
      DataStream(_traits) {}

    DataStream &operator<<(const std::string &val) override {
      m_vector.push_back(Promoter::promote(val));
      return *this;
    }

    DataStream &operator<<(long val) override {
      m_vector.push_back(Promoter::promote(val));
      return *this;
    }

    DataStream &operator<<(double val) override {
      m_vector.push_back(Promoter::promote(val));
      return *this;
    }

    DataStream &operator<<(bool val) override {
      m_vector.push_back(Promoter::promote(val));
      return *this;
    }

    DataStream &operator<<(char val) override {
      m_vector.push_back(Promoter::promote(val));
      return *this;
    }

    const std::vector<T> &vector()const {
      return m_vector;
    }
  private:
    std::vector<T> m_vector;
  };


  inline
  DataStream &operator<<(DataStream &_stream, int i) {
    return (_stream << (long)i);
  }

  inline
  DataStream &operator<<(DataStream &_stream, float f) {
    return (_stream << (double)f);
  }

  inline
  DataStream &operator<<(DataStream &_stream, unsigned int i) {
    return (_stream << (long)i);
  }

  inline
  DataStream &operator<<(DataStream &_stream, unsigned long i) {
    return (_stream << (long)i);
  }

  template<class T>
  inline
  DataStream &operator<<(DataStream &_stream, const std::vector<T> &vec) {
    for(auto it = vec.cbegin(); it != vec.cend(); ++it) {
      _stream << *it;
    }
    return _stream;
  }

}

#include "casm/external/Eigen/Dense"
namespace Eigen {
  template <typename Derived>
  CASM::DataStream &operator<<(CASM::DataStream &_stream, const MatrixBase<Derived> &value) {
    for(int i = 0; i < value.rows(); i++) {
      for(int j = 0; j < value.cols(); j++) {
        _stream << value(i, j);
      }
    }
    return _stream;
  }

}
#endif
