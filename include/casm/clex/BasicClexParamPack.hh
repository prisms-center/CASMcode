#ifndef BASICCLEXPARAMPACK_HH
#define BASICCLEXPARAMPACK_HH
#include <cstddef>
#include <iostream>
#include <string>
#include <sstream>

#include "casm/CASM_global_definitions.hh"
#include "casm/clex/ClexParamPack.hh"
#include "casm/casm_io/EnumIO.hh"

namespace CASM {
  class BasicClexParamKey : public ClexParamPack_impl::BaseKey {
  public:
    typedef ClexParamPack::size_type size_type;


    BasicClexParamKey(std::string const &_name = "", size_type _ind = -1) :
      ClexParamPack_impl::BaseKey(_name),
      m_index(_ind) {}

    size_type index() const {
      return m_index;
    }
  protected:

    ClexParamPack_impl::BaseKey *_clone() const override {
      return new BasicClexParamKey(*this);
    }

    size_type m_index;
  };

  /// \brief Abstract base class for reading/writing clexulator parameters
  class BasicClexParamPack : public ClexParamPack {
  public:

    enum class EvalMode {
      DEFAULT, READ, DYNAM
    };

    static const EvalMode DEFAULT;
    static const EvalMode READ;
    static const EvalMode DYNAM;


    typedef BasicClexParamKey Key;

    size_type size(ClexParamKey const &_key) const override {
      return size(*static_cast<Key const *>(_key.ptr()));
    }

    size_type size(Key const &_key) const {
      return m_data[_key.index()].size();
    }

    std::string eval_mode(ClexParamKey const &_key) const override {
      return to_string(eval_mode(*static_cast<Key const *>(_key.ptr())));
    }

    EvalMode eval_mode(Key const &_key) const {
      return m_eval[_key.index()];
    }

    std::vector<double> const &read(ClexParamKey const &_key) const override {
      return read(*static_cast<Key const *>(_key.ptr()));
    }

    std::vector<double> const &read(Key const &_key) const {
      return m_data[_key.index()];
    }

    double const &read(ClexParamKey const &_key, size_type _ind) const override {
      return read(*static_cast<Key const *>(_key.ptr()), _ind);
    }

    double const &read(Key const &_key, size_type _ind) const {
      return m_data[_key.index()][_ind];
    }

    void set_eval_mode(ClexParamKey const &_key, std::string const &_mode) override {
      set_eval_mode(*static_cast<Key const *>(_key.ptr()), from_string<EvalMode>(_mode));
    }

    void set_eval_mode(Key const &_key, EvalMode _mode) {
      m_eval[_key.index()] = _mode;
    }

    void write(ClexParamKey const &_key, std::vector<double> const &_val) override {
      write(*static_cast<Key const *>(_key.ptr()), _val);
    }

    void write(Key const &_key, std::vector<double> const &_val) {
      m_data[_key.index()] = _val;
    }

    void write(ClexParamKey const &_key, size_type _ind, double _val) override {
      write(*static_cast<Key const *>(_key.ptr()), _ind, _val);
    }

    void write(Key const &_key, size_type _ind, double _val) {
      m_data[_key.index()][_ind] = _val;
    }

    Key allocate(std::string const &_keyname, Index _size) {
      auto it = keys().find(_keyname);
      if(it != keys().end())
        throw std::runtime_error("Naming collision in BasicClexParamPack::allocate(), ClexParamPack already managing parameter allocation corresponding to name " + _keyname + ".");

      Key protokey(_keyname, m_data.size());

      m_data.push_back(std::vector<double>(_size));

      m_eval.push_back(EvalMode::DEFAULT);

      m_keys[_keyname] = protokey;
      return protokey;
    }


  private:
    std::vector<std::vector<double> > m_data;
    std::vector<EvalMode> m_eval;
  };


  template<>
  struct traits<BasicClexParamPack::EvalMode> {

    static const std::string name;

    static const std::multimap<BasicClexParamPack::EvalMode, std::vector<std::string> > strval;

  };


  std::ostream &operator<<(std::ostream &sout, const BasicClexParamPack::EvalMode &val) {
    sout << to_string<BasicClexParamPack::EvalMode>(val);
    return sout;
  }

  std::istream &operator>>(std::istream &sin, BasicClexParamPack::EvalMode &val) {
    std::string s;
    sin >> s;
    val = from_string<BasicClexParamPack::EvalMode>(s);
    return sin;
  }


  const std::string traits<BasicClexParamPack::EvalMode>::name = "clex_eval_mode";

  const std::multimap<BasicClexParamPack::EvalMode, std::vector<std::string> > traits<BasicClexParamPack::EvalMode>::strval = {
    {BasicClexParamPack::EvalMode::DEFAULT, {"Default", "DEFAULT", "default"} },
    {BasicClexParamPack::EvalMode::READ, {"Read", "READ", "read"} },
    {BasicClexParamPack::EvalMode::DYNAM, {"Dynamic", "Dynam", "DYNAMIC", "DYNAM", "dynamic", "dynam"} }
  };

  const BasicClexParamPack::EvalMode BasicClexParamPack::DEFAULT = BasicClexParamPack::EvalMode::DEFAULT;
  const BasicClexParamPack::EvalMode BasicClexParamPack::DYNAM = BasicClexParamPack::EvalMode::DYNAM;
  const BasicClexParamPack::EvalMode BasicClexParamPack::READ = BasicClexParamPack::EvalMode::READ;



}
#endif
