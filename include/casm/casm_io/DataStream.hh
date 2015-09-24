#ifndef DATASTREAM_HH
#define DATASTREAM_HH
#include <vector>
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

    DataStream &operator<<(const std::string &) {
      return increment();
    }

    DataStream &operator<<(long) {
      return increment();
    }

    DataStream &operator<<(double) {
      return increment();
    }

    DataStream &operator<<(bool) {
      return increment();
    }

    DataStream &operator<<(char) {
      return increment();
    }

    DataStream &newline() {
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
  namespace ArrayDataStream_impl {
    template<typename OutType>
    struct DataStreamPromoter {
      template<typename InType>
      static OutType promote(InType);
    };

    template<>
    struct DataStreamPromoter<std::string> {
      template<typename InType>
      static std::string promote(InType);

    };

    template<>
    struct DataStreamPromoter<double> {
      template<typename InType>
      static double promote(InType);

    };

    template<>
    struct DataStreamPromoter<long> {
      template<typename InType>
      static long promote(InType);

    };

    template<>
    struct DataStreamPromoter<bool> {
      template<typename InType>
      static bool promote(InType);

    };

    template<>
    struct DataStreamPromoter<char> {
      template<typename InType>
      static char promote(InType);

    };

  }

  template <typename T, typename Promoter = ArrayDataStream_impl::DataStreamPromoter<T> >
  class ArrayDataStream : public DataStream {
  public:

    ArrayDataStream(DataStreamTraits _traits = none) :
      DataStream(_traits) {}

    DataStream &operator<<(const std::string &val) {
      m_array.push_back(Promoter::promote(val));
      return *this;
    }

    DataStream &operator<<(long val) {
      m_array.push_back(Promoter::promote(val));
      return *this;
    }

    DataStream &operator<<(double val) {
      m_array.push_back(Promoter::promote(val));
      return *this;
    }

    DataStream &operator<<(bool val) {
      m_array.push_back(Promoter::promote(val));
      return *this;
    }

    DataStream &operator<<(char val) {
      m_array.push_back(Promoter::promote(val));
      return *this;
    }

    const std::vector<T> &array()const {
      return m_array;
    }
  private:
    std::vector<T> m_array;
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

#endif
