#ifndef FORMATFLAG_HH
#define FORMATFLAG_HH

#include <iostream>

namespace CASM {

  class FormatFlag {
  public:
    enum Value {
      cart = 0,
      va_off = 0,
      header_on = 0,
      frac = (1u << 0),
      va_on = (1u << 1),
      header_off = (1u << 2)
                   // add up to (1u << 7)
    };

    FormatFlag(int _value) : m_value(_value) {}

    FormatFlag(std::ostream &_stream) : m_value(_stream.iword(iword_index())) {}

    int value() const {
      return m_value;
    }

    std::ostream &set(std::ostream &_stream) const {
      _stream.iword(iword_index()) = value();
      return _stream;
    }

    bool print_header() const {
      return !bool(m_value & header_off);
    }

    FormatFlag &print_header(bool _set) {
      if(_set)
        m_value &= ~header_off;
      else
        m_value |= header_off;
      return *this;
    }

    bool print_va() const {
      return bool(m_value & va_on);
    }

    FormatFlag operator|(int RHS) {
      FormatFlag tflag(*this);
      return tflag |= RHS;
    }

    FormatFlag operator^(int RHS) {
      FormatFlag tflag(*this);
      return tflag ^= RHS;
    }

    FormatFlag operator&(int RHS) {
      FormatFlag tflag(*this);
      return tflag &= RHS;
    }

    FormatFlag &operator|=(int RHS) {
      m_value |= RHS;
      return *this;
    }

    FormatFlag &operator^=(int RHS) {
      m_value ^= RHS;
      return *this;
    }

    FormatFlag &operator&=(int RHS) {
      m_value &= RHS;
      return *this;
    }

  private:
    static int iword_index();
    int m_value;
  };

  std::ostream &operator<<(std::ostream &_stream, const FormatFlag &flag);

}
#endif

