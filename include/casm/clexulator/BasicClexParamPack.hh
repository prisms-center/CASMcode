#ifndef CASM_clexulator_BasicClexParamPack
#define CASM_clexulator_BasicClexParamPack
#include <cstddef>
#include <iostream>
#include <sstream>
#include <string>

#include "casm/casm_io/enum/json_io.hh"
#include "casm/casm_io/enum/stream_io.hh"
#include "casm/clexulator/ClexParamPack.hh"
#include "casm/global/definitions.hh"

namespace CASM {
namespace clexulator {
/** \ingroup Clexulator
 * @{ */

class BasicClexParamKey : public ClexParamPack_impl::BaseKey {
 public:
  typedef ClexParamPack::size_type size_type;

  BasicClexParamKey(std::string const &_name = "", bool _standalone = false,
                    size_type _ind = -1)
      : ClexParamPack_impl::BaseKey(_name, _standalone), m_index(_ind) {}

  ~BasicClexParamKey() {}

  size_type index() const { return m_index; }

 protected:
  ClexParamPack_impl::BaseKey *_clone() const override {
    return new BasicClexParamKey(*this);
  }

  size_type m_index;
};

/// \brief Abstract base class for reading/writing clexulator parameters
/// Parameters are assume be naturally representable as 1D or 2D arrays
class BasicClexParamPack : public ClexParamPack {
 public:
  enum class EvalMode { DEFAULT, READ, DYNAM };

  static const EvalMode DEFAULT;
  static const EvalMode READ;
  static const EvalMode DYNAM;

  using Key = BasicClexParamKey;
  using DoubleReference = Eigen::MatrixXd::CoeffReturnType;

  template <typename Scalar>
  using Val = ValAccess<Scalar>;

  template <typename Scalar>
  friend struct ValAccess;

  size_type size(ClexParamKey const &_key) const override {
    return size(*static_cast<Key const *>(_key.ptr()));
  }

  size_type size(Key const &_key) const { return m_data[_key.index()].cols(); }

  size_type dim(ClexParamKey const &_key) const override {
    return dim(*static_cast<Key const *>(_key.ptr()));
  }

  size_type dim(Key const &_key) const { return m_data[_key.index()].rows(); }

  std::string eval_mode(ClexParamKey const &_key) const override {
    return to_string(eval_mode(*static_cast<Key const *>(_key.ptr())));
  }

  EvalMode eval_mode(Key const &_key) const { return m_eval[_key.index()]; }

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

  double const &read(ClexParamKey const &_key, size_type _i,
                     size_type _j) const override {
    return read(*static_cast<Key const *>(_key.ptr()), _i, _j);
  }

  double const &read(Key const &_key, size_type _i, size_type _j) const {
    return m_data[_key.index()](_i, _j);
  }

  void set_eval_mode(ClexParamKey const &_key,
                     std::string const &_mode) override {
    set_eval_mode(*static_cast<Key const *>(_key.ptr()),
                  from_string<EvalMode>(_mode));
  }

  void set_eval_mode(Key const &_key, EvalMode _mode) {
    m_eval[_key.index()] = _mode;
  }

  void write(ClexParamKey const &_key,
             Eigen::Ref<const Eigen::MatrixXd> const &_val) override {
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

  void write(ClexParamKey const &_key, size_type _i, size_type _j,
             double _val) override {
    write(*static_cast<Key const *>(_key.ptr()), _i, _j, _val);
  }

  void write(Key const &_key, size_type _i, size_type _j, double _val) {
    m_data[_key.index()](_i, _j) = _val;
  }

  Key allocate(std::string const &_keyname, Index _rows, Index _cols,
               bool _independent) {
    auto it = keys().find(_keyname);
    if (it != keys().end())
      throw std::runtime_error(
          "Naming collision in BasicClexParamPack::allocate(), ClexParamPack "
          "already managing parameter allocation corresponding to name " +
          _keyname + ".");

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

template <>
struct ValAccess<double> {
  using size_type = BasicClexParamPack::size_type;

  static double const &get(BasicClexParamPack const &_pack,
                           BasicClexParamKey const &_key, size_type i) {
    return _pack.m_data[_key.index()](i, 0);
  }

  static double const &get(BasicClexParamPack const &_pack,
                           BasicClexParamKey const &_key, size_type i,
                           size_type j) {
    return _pack.m_data[_key.index()](i, j);
  }

  static void set(BasicClexParamPack &_pack, BasicClexParamKey const &_key,
                  Eigen::Ref<const Eigen::MatrixXd> const &_val) {
    _pack.m_data[_key.index()] = _val;
  }

  template <typename Scalar2>
  static void set(BasicClexParamPack &_pack, BasicClexParamKey const &_key,
                  size_type i, Scalar2 const &_val) {
    _pack.m_data[_key.index()](i, 0) = _val;
  }

  template <typename Scalar2>
  static void set(BasicClexParamPack &_pack, BasicClexParamKey const &_key,
                  size_type i, size_type j, Scalar2 const &_val) {
    _pack.m_data[_key.index()](i, j) = _val;
  }
};

}  // namespace clexulator

template <>
struct traits<clexulator::BasicClexParamPack::EvalMode> {
  static const std::string name;

  static const std::multimap<clexulator::BasicClexParamPack::EvalMode,
                             std::vector<std::string> >
      strval;
};

std::ostream &operator<<(std::ostream &sout,
                         const clexulator::BasicClexParamPack::EvalMode &val) {
  sout << to_string<clexulator::BasicClexParamPack::EvalMode>(val);
  return sout;
}

std::istream &operator>>(std::istream &sin,
                         clexulator::BasicClexParamPack::EvalMode &val) {
  std::string s;
  sin >> s;
  val = from_string<clexulator::BasicClexParamPack::EvalMode>(s);
  return sin;
}

const std::string traits<clexulator::BasicClexParamPack::EvalMode>::name =
    "clex_eval_mode";

const std::multimap<clexulator::BasicClexParamPack::EvalMode,
                    std::vector<std::string> >
    traits<clexulator::BasicClexParamPack::EvalMode>::strval = {
        {clexulator::BasicClexParamPack::EvalMode::DEFAULT,
         {"Default", "DEFAULT", "default"}},
        {clexulator::BasicClexParamPack::EvalMode::READ,
         {"Read", "READ", "read"}},
        {clexulator::BasicClexParamPack::EvalMode::DYNAM,
         {"Dynamic", "Dynam", "DYNAMIC", "DYNAM", "dynamic", "dynam"}}};

namespace clexulator {

const BasicClexParamPack::EvalMode BasicClexParamPack::DEFAULT =
    BasicClexParamPack::EvalMode::DEFAULT;
const BasicClexParamPack::EvalMode BasicClexParamPack::DYNAM =
    BasicClexParamPack::EvalMode::DYNAM;
const BasicClexParamPack::EvalMode BasicClexParamPack::READ =
    BasicClexParamPack::EvalMode::READ;

/** @} */
}  // namespace clexulator
}  // namespace CASM
#endif
