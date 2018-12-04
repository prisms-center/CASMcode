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


    BasicClexParamKey(std::string const &_name = "", bool _standalone = false, size_type _ind = -1) :
      ClexParamPack_impl::BaseKey(_name, _standalone),
      m_index(_ind) {}

    ~BasicClexParamKey() {
    }

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


    using Key = BasicClexParamKey;
    using DoubleReference = Eigen::MatrixXd::CoeffReturnType;

    size_type size(ClexParamKey const &_key) const override {
      return size(*static_cast<Key const *>(_key.ptr()));
    }

    size_type size(Key const &_key) const {
      return m_data[_key.index()].cols();
    }

    size_type dim(ClexParamKey const &_key) const override {
      return dim(*static_cast<Key const *>(_key.ptr()));
    }

    size_type dim(Key const &_key) const {
      return m_data[_key.index()].rows();
    }

    std::string eval_mode(ClexParamKey const &_key) const override {
      return to_string(eval_mode(*static_cast<Key const *>(_key.ptr())));
    }

    EvalMode eval_mode(Key const &_key) const {
      return m_eval[_key.index()];
    }

    Eigen::MatrixXd const &read(ClexParamKey const &_key) const override {
      return read(*static_cast<Key const *>(_key.ptr()));
    }

    Eigen::MatrixXd const &read(Key const &_key) const {
      return m_data[_key.index()];
    }

    double const &read(ClexParamKey const &_key, size_type _ind) const override {
      return read(*static_cast<Key const *>(_key.ptr()), _ind);
    }

    double const &read(Key const &_key, size_type _ind) const {
      return m_data[_key.index()](_ind, 0);
    }

    double const &read(ClexParamKey const &_key, size_type _i, size_type _j) const override {
      return read(*static_cast<Key const *>(_key.ptr()), _i, _j);
    }

    double const &read(Key const &_key, size_type _i, size_type _j) const {
      return m_data[_key.index()](_i, _j);
    }

    void set_eval_mode(ClexParamKey const &_key, std::string const &_mode) override {
      set_eval_mode(*static_cast<Key const *>(_key.ptr()), from_string<EvalMode>(_mode));
    }

    void set_eval_mode(Key const &_key, EvalMode _mode) {
      m_eval[_key.index()] = _mode;
    }

    void write(ClexParamKey const &_key, Eigen::Ref<const Eigen::MatrixXd> const &_val) override {
      write(*static_cast<Key const *>(_key.ptr()), _val);
    }

    void write(Key const &_key, Eigen::Ref<const Eigen::MatrixXd> const &_val) {
      m_data[_key.index()] = _val;
    }

    void write(ClexParamKey const &_key, size_type _i, double _val) override {
      write(*static_cast<Key const *>(_key.ptr()), _i, _val);
    }

    void write(Key const &_key, size_type _i, double _val) {
      m_data[_key.index()](_i, 0) = _val;
    }

    void write(ClexParamKey const &_key, size_type _i, size_type _j, double _val) override {
      write(*static_cast<Key const *>(_key.ptr()), _i, _j, _val);
    }

    void write(Key const &_key, size_type _i, size_type _j, double _val) {
      m_data[_key.index()](_i, _j) = _val;
    }

    Key allocate(std::string const &_keyname, Index _rows, Index _cols, bool _independent) {
      auto it = keys().find(_keyname);
      if(it != keys().end())
        throw std::runtime_error("Naming collision in BasicClexParamPack::allocate(), ClexParamPack already managing parameter allocation corresponding to name " + _keyname + ".");

      Key protokey(_keyname, true, m_data.size());

      m_data.push_back(Eigen::MatrixXd::Zero(_rows, _cols));

      m_eval.push_back(EvalMode::DEFAULT);

      m_keys[_keyname] = protokey;
      return protokey;
    }


  private:
    std::vector<Eigen::MatrixXd> m_data;
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
