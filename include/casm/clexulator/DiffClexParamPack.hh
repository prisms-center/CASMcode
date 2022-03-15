#ifndef CASM_clexulator_DiffClexParamPack
#define CASM_clexulator_DiffClexParamPack
#include <cstddef>
#include <iostream>
#include <sstream>
#include <string>

#include "casm/casm_io/container/stream_io.hh"
#include "casm/casm_io/enum/json_io.hh"
#include "casm/casm_io/enum/stream_io.hh"
#include "casm/clexulator/ClexParamPack.hh"
#include "casm/external/fadbad/badiff.h"
#include "casm/external/fadbad/fadiff.h"
#include "casm/global/definitions.hh"

namespace CASM {
namespace clexulator {
/** \ingroup Clexulator
 * @{ */

namespace DiffClexParamPack_impl {
enum class EvalMode {
  DEFAULT,
  READ,
  DIFF,
  DYNAM

};

using DiffScalar = fadbad::B<fadbad::F<double> >;

}  // namespace DiffClexParamPack_impl

// DiffClexParamPack provides an interface for reading and writing an automaitic
// differentiation scalar type, which is implemented using the FADBAD++ library.
// By utilizing the autodiff scalar type as the input and output of algebraic
// expressions (as in the clexulator), first- and higher-order derivatives of
// the expressions can also be calculated simply and efficiently.
// DiffClexParamPack stores 'standalone' values, which comprise the independent
// variables utilized by the clexulator and the dependent function returnvalues,
// but it can also be queried for 'derived' values, which are obtained by
// performing additional operations on the standalone values. The standalone and
// derived values are both queried in an identical fashion, by passing a
// 'ClexParamKey' to the DiffClexParamPack corresponding to the quantity of
// interest.
class DiffClexParamPack;

// Container struct for holding a particular standalone value
class DiffScalarContainer;

// Abstract base class for all keys that interact with DiffClexParamPack
class DiffClexParamKey;

// Key implementation for standalone values
class DiffClexParamValKey;

// Key implementation for derived values corresponding to first derivatives
class DiffClexParamGradKey;

// Key implementation for derived values corresponding to second derivatives
class DiffClexParamHessKey;

/// \brief Abstract base class for all keys that interact with DiffClexParamPack
/// DiffClexParamPack values are assumed to be 1D, 2D, or to be naturally
/// accessed via 2D slices Standalone values are assumed to be 1D or 2D (for
/// example, a 3xN displacement field) Derived values typically have more
/// dimensions (for example, the derivative of the Mx1 correlation vector with
/// respect to the displacement field has 3 indices, whild the second derivative
/// wrt the displacement field has 5 indices). Derivatives are always asociated
/// with the first *independent* standalone value in the key name. The key
/// "diff/corr/disp_var" would be associated with "disp_var", which is a 3xN
/// Matrix, and thus the resulting derived value (which is the first derifative
/// of correlations with respect to displacements) will be sliced into 3xN
/// slices. Each slice corresponds to a particular entry in the Mx1 correlation
/// vector (i.e., slice 'j' is the gradient of correlation value j with respect
/// to the entire set of displacement variables).  The 'j' index that specifies
/// the particular slice is called a 'secondary identifier', and it can be
/// specified via ClexParamKey::set_identifiers().
class DiffClexParamKey : public ClexParamPack_impl::BaseKey {
 public:
  typedef ClexParamPack::size_type size_type;

  /// \brief DiffClexParamKey Constructor
  /// \param _name The name of the associated value
  /// \param _ind Index of the standalone value that is the primary identifier
  /// for the associated value \param _offset Constant offset for any secondary
  /// identifiers for the value of interest \param _stride
  DiffClexParamKey(std::string const &_name, bool _standalone, size_type _ind,
                   std::vector<Index> const &_offset = {},
                   std::vector<Index> const &_stride = {})
      : ClexParamPack_impl::BaseKey(_name, _standalone, _offset, _stride),
        m_index(_ind) {
    // std::cout << "Constructed key " << _name << " at index " << m_index << "
    // with offset " << _offset << " and stride " << _stride << std::endl;
  }

  DiffClexParamKey() : DiffClexParamKey("", -1, false) {}

  /// \brief Destructor of abstract base class must be virtual
  virtual ~DiffClexParamKey(){};

  /// \brief Operate on associated DiffScalarContainer to ensure the standalone
  /// or derived value associated with this key is present in its cache matrix.
  virtual Eigen::MatrixXd const &eval(
      DiffScalarContainer const &_data,
      DiffClexParamPack_impl::EvalMode mode) const = 0;

  /// \brief Access the primary identifier of this value (which determines the
  /// associated standalone value and/or the shape of the dataslice
  size_type index() const { return m_index; }

 private:
  size_type m_index;
};

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/// \brief container struct for storing a standalone value as 2D arrays of both
/// autodiff Scalars doubles
struct DiffScalarContainer {
  using DiffScalar = DiffClexParamPack_impl::DiffScalar;

  /// \brief Name of the standalone value stored in this container. Should match
  /// name of associated key
  std::string m_name;

  /// \brief linear offset for start of value index in linear array
  /// FADBAD forward (backward) autodiff types reference independent (dependent)
  /// via a linear index For a particular entry (i,j) into the 2d container for
  /// a specific type (e.g., "occ_site_func"), the linear index is defined as
  /// l=lbegin+i*m_cache.cols()+j
  Index lbegin;

  /// \brief Flag, set to true if this container holds independent values (false
  /// if dependent)
  bool m_independent;

  /// \brief Stored data represented as autodif scalar type
  // In future, migrate to Eigen::Matrix<DiffScalar> ?
  std::vector<std::vector<DiffScalar> > m_data;

  /// \brief Stored data, represented as matrix of doubles
  mutable Eigen::MatrixXd m_data_double;

  /// \brief Derived data cache, represented as matrix of doubles
  mutable Eigen::MatrixXd m_cache;

  /// \brief DiffSclarContainer constructor
  /// \param _name Name of the standalone value stored in this container, should
  /// match name of associated key \param param_dim Dimension of vector
  /// parameter (e.g., 3 for "disp_var"), corresponds to columns of m_data,
  /// m_data_double, etc \param num_params Number of vector parameters (e.g., 1
  /// for global DoF, number of sites for local DoF), corresponds to rows of
  /// m_data, m_data_double, etc
  DiffScalarContainer(std::string const &_name, Index param_dim,
                      Index num_params, Index _lbegin, bool _independent)
      : m_name(_name),
        lbegin(_lbegin),
        m_independent(_independent),
        m_data(param_dim, std::vector<DiffScalar>(num_params)),
        m_data_double(param_dim, num_params),
        m_cache(param_dim, num_params) {
    // std::cout << "Allocated "<< (_independent? "in" : "") <<  "dependent
    // Container '" << _name << "' size " << param_dim << "x" << num_params << "
    // starting at " << _lbegin << std::endl;
  }

  /// \brief Gradient of dependent function index f wrt the independent
  /// variables contained in here \param f secondary identifier denoting linear
  /// index of function F(f) among the set F of all dependent variables \result
  /// matrix G, with G(i,j) = Grad[F(f)](l(i,j)), where 'l' is linear index of
  /// each (i,j) entry in m_data.
  Eigen::MatrixXd const &grad(Index f) const {
    // std::cout << "WORKING ON f = " << f << "\n";
    for (Index i = 0; i < m_data.size(); ++i) {
      for (Index j = 0; j < m_data[i].size(); ++j) {
        m_cache(i, j) = m_data[i][j].d(f).val();
        // std::cout << "FOR (i,j) : (" << i << ", " << j << "), val=" <<
        // m_data[i][j].val().val() << "; grad=" << m_cache(i,j) << std::endl;
      }
    }
    return m_cache;
  }

  /// \brief Hessian components of dependent function index f corresponding to
  /// independent variable index 'a' and all independent variables managed by
  /// this \param f secondary identifier denoting linear index of function F(f)
  /// among the set F of all dependent variables \param a secondary identifier
  /// denoting linear index of independent variable 'a' amont the set of all
  /// independent variables \result matrix H with M(i,j) = H[F(f)](a,l(i,j)),
  /// where 'l' is linear index of each (i,j) entry in m_data.
  Eigen::MatrixXd const &hess(Index f, Index a) const {
    for (Index i = 0; i < m_data.size(); ++i) {
      for (Index j = 0; j < m_data[i].size(); ++j) {
        m_cache(i, j) = m_data[i][j].d(f).d(a);
      }
    }
    return m_cache;
  }

  /// \brief Must be called on each independent parameter before function
  /// evaluation if derivatives are desired. \param Total number of independent
  /// parameters managed by DiffClexParamPack
  void pre_eval(Index N) {
    Index l = lbegin;
    for (auto &dvec : m_data) {
      for (auto &datum : dvec) {
        datum.val().diff(l++, N);
      }
    }
  }

  /// \brief Must be called on each dependent parameter after function
  /// evaluation if derivatives are desired. \param Total number of dependent
  /// parameters managed by DiffClexParamPack
  void post_eval(Index N) {
    Index l = lbegin;
    for (auto &dvec : m_data) {
      for (auto &datum : dvec) {
        datum.diff(l++, N);
      }
    }
  }

  /// \brief Copy values of autodiff scalars from m_data to m_data_double
  /// (matrix of doubles) in preparation for external access
  void eval_double() {
    for (Index i = 0; i < m_cache.rows(); ++i) {
      for (Index j = 0; j < m_cache.cols(); ++j) {
        m_data_double(i, j) = m_data[i][j].val().val();
      }
    }
  }

  /// \brief Returns one past end of max linear index range allowed for this
  /// parameter set
  Index linear_index_end() const { return lbegin + m_cache.size(); }
};

/// \brief Class for managing the set of all dependent and independent paramters
/// used/generated by clexulator
class DiffClexParamPack : public ClexParamPack {
 public:
  using DiffScalar = DiffClexParamPack_impl::DiffScalar;

  using EvalMode = DiffClexParamPack_impl::EvalMode;

  static const EvalMode DEFAULT;
  static const EvalMode READ;
  static const EvalMode DIFF;
  static const EvalMode DYNAM;

  using BaseKey = DiffClexParamKey;
  using Key = DiffClexParamValKey;
  using DoubleReference = Eigen::MatrixXd::CoeffReturnType;

  template <typename Scalar>
  using Val = ValAccess<Scalar>;

  template <typename Scalar>
  friend struct ValAccess;

  /// \brief Default constructor initializes evaluation mode and zeros numbers
  /// of managed parameters
  DiffClexParamPack()
      : m_tot_eval_mode(DEFAULT), m_N_independent(0), m_N_dependent(0) {}

  size_type size(ClexParamKey const &_key) const override {
    return size(*static_cast<BaseKey const *>(_key.ptr()));
  }

  size_type size(BaseKey const &_key) const {
    return m_data[_key.index()].m_cache.cols();
  }

  size_type dim(ClexParamKey const &_key) const override {
    return dim(*static_cast<BaseKey const *>(_key.ptr()));
  }

  size_type dim(BaseKey const &_key) const {
    return m_data[_key.index()].m_cache.rows();
  }

  std::string eval_mode(ClexParamKey const &_key) const override {
    return to_string(eval_mode(*static_cast<BaseKey const *>(_key.ptr())));
  }

  EvalMode eval_mode(BaseKey const &_key) const { return m_eval[_key.index()]; }

  EvalMode eval_mode() const { return m_tot_eval_mode; }

  Eigen::MatrixXd const &read(ClexParamKey const &_key) const override {
    return read(*static_cast<BaseKey const *>(_key.ptr()));
  }

  Eigen::MatrixXd const &read(BaseKey const &_key) const {
    return _key.eval(m_data[_key.index()], eval_mode());
  }

  double const &read(ClexParamKey const &_key, size_type _ind) const override {
    return read(*static_cast<BaseKey const *>(_key.ptr()), _ind);
  }

  double const &read(BaseKey const &_key, size_type _ind) const {
    return _key.eval(m_data[_key.index()], eval_mode())(_ind, 0);
  }

  double const &read(ClexParamKey const &_key, size_type _i,
                     size_type _j) const override {
    return read(*static_cast<BaseKey const *>(_key.ptr()), _i, _j);
  }

  double const &read(BaseKey const &_key, size_type _i, size_type _j) const {
    return _key.eval(m_data[_key.index()], eval_mode())(_i, _j);
  }

  void write(ClexParamKey const &_key,
             Eigen::Ref<const Eigen::MatrixXd> const &_val) override {
    write(*static_cast<BaseKey const *>(_key.ptr()), _val);
  }

  void write(BaseKey const &_key,
             Eigen::Ref<const Eigen::MatrixXd> const &_val);

  void write(ClexParamKey const &_key, size_type _i, double _val) override {
    write(*static_cast<BaseKey const *>(_key.ptr()), _i, _val);
  }

  void write(BaseKey const &_key, size_type _i, double _val);

  void write(ClexParamKey const &_key, size_type _i, size_type _j,
             double _val) override {
    write(*static_cast<BaseKey const *>(_key.ptr()), _i, _j, _val);
  }

  void write(BaseKey const &_key, size_type _i, size_type _j, double _val);

  Key allocate(std::string const &_keyname, Index _rows, Index _cols,
               bool _independent);

  void set_eval_mode(ClexParamKey const &_key,
                     std::string const &_mode) override {
    set_eval_mode(*static_cast<BaseKey const *>(_key.ptr()),
                  from_string<EvalMode>(_mode));
  }

  /// \brief Set evaluation mode of specified parameter type
  /// \param _key Specifies parameter for which evluation mode is being set
  /// \param _mode New value for evaluation mode
  void set_eval_mode(BaseKey const &_key, EvalMode _mode) {
    m_eval[_key.index()] = _mode;
    m_tot_eval_mode = DEFAULT;
    for (auto const &eval : m_eval) {
      if (eval == DIFF) {
        m_tot_eval_mode = DIFF;
        break;
      }
    }
  }

  /// \brief Prepares all independent parameters before function evaluation if
  /// derivatives are desired
  void pre_eval() {
    if (eval_mode() == DIFF) {
      for (DiffScalarContainer &datum : m_data) {
        if (datum.m_independent) {
          datum.pre_eval(m_N_independent);
        }
      }
    }
  }

  /// \brief Processes all dependent parameters after function evaluation if
  /// derivatives are desired
  void post_eval() {
    if (eval_mode() == DIFF) {
      for (DiffScalarContainer &datum : m_data) {
        if (!datum.m_independent) {
          datum.post_eval(m_N_dependent);
        }
      }
    }
  }

 private:
  std::vector<DiffScalarContainer> m_data;
  std::vector<EvalMode> m_eval;
  EvalMode m_tot_eval_mode;
  Index m_N_independent;
  Index m_N_dependent;
};

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class DiffClexParamValKey : public DiffClexParamKey {
 public:
  typedef ClexParamPack::size_type size_type;

  DiffClexParamValKey(std::string const &_name = "", size_type _ind = -1)
      : DiffClexParamKey(_name, true, _ind) {
    // std::cout << "Constructing value key " << _name << "\n";
  }

  Eigen::MatrixXd const &eval(DiffScalarContainer const &_data,
                              DiffClexParamPack::EvalMode mode) const override {
    if (mode == DiffClexParamPack::DIFF) {
      for (Index i = 0; i < _data.m_cache.rows(); ++i) {
        for (Index j = 0; j < _data.m_cache.cols(); ++j) {
          _data.m_cache(i, j) = _data.m_data[i][j].val().val();
        }
      }
    }
    return _data.m_cache;
  }

 protected:
  ClexParamPack_impl::BaseKey *_clone() const override {
    return new DiffClexParamValKey(*this);
  }
};

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class DiffClexParamGradKey : public DiffClexParamKey {
 public:
  typedef ClexParamPack::size_type size_type;

  DiffClexParamGradKey(std::string const &_name = "", size_type _ix = -1,
                       std::vector<Index> const &_offset = {},
                       std::vector<Index> const &_stride = {})
      : DiffClexParamKey(_name, false, _ix, _offset, _stride) {
    // std::cout << "Constructing grad key " << _name << ", ID " << _ix << "\n";
  }

  Eigen::MatrixXd const &eval(DiffScalarContainer const &_data,
                              DiffClexParamPack::EvalMode mode) const override {
    if (mode != DiffClexParamPack::DIFF)
      throw std::runtime_error("Requested Gradient parameter " + name() +
                               " for incompatible evaluation mode " +
                               to_string(mode) + ".");
    _data.grad(_l(0));
    return _data.m_cache;
  }

 protected:
  ClexParamPack_impl::BaseKey *_clone() const override {
    return new DiffClexParamGradKey(*this);
  }
};

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class DiffClexParamHessKey : public DiffClexParamKey {
 public:
  typedef ClexParamPack::size_type size_type;

  DiffClexParamHessKey(std::string const &_name = "", size_type _ix = -1,
                       std::vector<Index> const &_offset = {},
                       std::vector<Index> const &_stride = {})
      : DiffClexParamKey(_name, false, _ix, _offset, _stride) {
    // std::cout << "Constructing hess key " << _name << ", ID " << _ix << "\n";
  }

  Eigen::MatrixXd const &eval(DiffScalarContainer const &_data,
                              DiffClexParamPack::EvalMode mode) const override {
    if (mode != DiffClexParamPack::DIFF)
      throw std::runtime_error("Requested Hessian parameter " + name() +
                               " for incompatible evaluation mode " +
                               to_string(mode) + ".");
    _data.hess(_l(0), _l(1));
    return _data.m_cache;
  }

 protected:
  ClexParamPack_impl::BaseKey *_clone() const override {
    return new DiffClexParamHessKey(*this);
  }
};

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

template <>
struct traits<DiffClexParamPack::EvalMode> {
  static const std::string name;

  static const std::multimap<DiffClexParamPack::EvalMode,
                             std::vector<std::string> >
      strval;
};

template <>
struct ValAccess<double> {
  using size_type = DiffClexParamPack::size_type;

  static double const &get(DiffClexParamPack const &_pack,
                           DiffClexParamKey const &_key, size_type i) {
    return _pack.m_data[_key.index()].m_cache(i, 0);
  }

  static double const &get(DiffClexParamPack const &_pack,
                           DiffClexParamKey const &_key, size_type i,
                           size_type j) {
    return _pack.m_data[_key.index()].m_cache(i, j);
  }

  static void set(DiffClexParamPack &_pack, DiffClexParamKey const &_key,
                  Eigen::Ref<const Eigen::MatrixXd> const &_val) {
    _pack.m_data[_key.index()].m_cache = _val;
  }

  template <typename Scalar2>
  static void set(DiffClexParamPack &_pack, DiffClexParamKey const &_key,
                  size_type i, Scalar2 const &_val) {
    _pack.m_data[_key.index()].m_cache(i, 0) = _val;
  }

  template <typename Scalar2>
  static void set(DiffClexParamPack &_pack, DiffClexParamKey const &_key,
                  size_type i, size_type j, Scalar2 const &_val) {
    _pack.m_data[_key.index()].m_cache(i, j) = _val;
  }
};

template <>
struct ValAccess<typename DiffClexParamPack::DiffScalar> {
  using size_type = DiffClexParamPack::size_type;

  static typename DiffClexParamPack::DiffScalar const &get(
      DiffClexParamPack const &_pack, DiffClexParamKey const &_key,
      size_type i) {
    return _pack.m_data[_key.index()].m_data[i][0];
  }

  static typename DiffClexParamPack::DiffScalar const &get(
      DiffClexParamPack const &_pack, DiffClexParamKey const &_key, size_type i,
      size_type j) {
    return _pack.m_data[_key.index()].m_data[i][j];
  }

  static void set(DiffClexParamPack &_pack, DiffClexParamKey const &_key,
                  Eigen::Ref<const Eigen::MatrixXd> const &_val) {
    assert(_val.rows() == _pack.m_data[_key.index()].m_data.size());
    for (Index i = 0; i < _val.rows(); ++i) {
      assert(_val.cols() == _pack.m_data[_key.index()].m_data[i].size());
      for (Index j = 0; j < _val.cols(); ++j) {
        _pack.m_data[_key.index()].m_data[i][j] = _val(i, j);
      }
    }
  }

  template <typename Scalar2>
  static void set(DiffClexParamPack &_pack, DiffClexParamKey const &_key,
                  size_type i, Scalar2 const &_val) {
    _pack.m_data[_key.index()].m_data[i][0] = _val;
  }

  template <typename Scalar2>
  static void set(DiffClexParamPack &_pack, DiffClexParamKey const &_key,
                  size_type i, size_type j, Scalar2 const &_val) {
    _pack.m_data[_key.index()].m_data[i][j] = _val;
  }
};

inline void DiffClexParamPack::write(
    BaseKey const &_key, Eigen::Ref<const Eigen::MatrixXd> const &_val) {
  if (!_key.standalone()) {
    throw std::runtime_error("Cannot write to dependent parameter " +
                             _key.name() + ".");
  }
  if (eval_mode() == DEFAULT) {
    Val<double>::set(*this, _key, _val);
  } else {
    Val<DiffScalar>::set(*this, _key, _val);
  }
}

inline void DiffClexParamPack::write(BaseKey const &_key, size_type _i,
                                     double _val) {
  if (!_key.standalone()) {
    throw std::runtime_error("Cannot write to dependent parameter " +
                             _key.name() + ".");
  }

  if (eval_mode() == DEFAULT) {
    Val<double>::set(*this, _key, _i, _val);
  } else {
    Val<DiffScalar>::set(*this, _key, _i, _val);
  }
}

inline void DiffClexParamPack::write(BaseKey const &_key, size_type _i,
                                     size_type _j, double _val) {
  if (!_key.standalone()) {
    throw std::runtime_error("Cannot write to dependent parameter " +
                             _key.name() + ".");
  }

  if (eval_mode() == DEFAULT) {
    Val<double>::set(*this, _key, _i, _j, _val);
  } else {
    Val<DiffScalar>::set(*this, _key, _i, _j, _val);
  }
}

inline DiffClexParamPack::Key DiffClexParamPack::allocate(
    std::string const &_keyname, Index _rows, Index _cols, bool _independent) {
  auto it = keys().find(_keyname);
  if (it != keys().end())
    throw std::runtime_error(
        "Naming collision in DiffClexParamPack::allocate(), ClexParamPack "
        "already managing parameter allocation corresponding to name " +
        _keyname + ".");
  Key protokey(_keyname, m_data.size());
  m_keys[_keyname] = protokey;
  if (_independent) {
    m_data.push_back(
        DiffScalarContainer(_keyname, _rows, _cols, m_N_independent, true));
    m_N_independent += _rows * _cols;
    Index i = protokey.index();
    for (Index f = 0; f < m_data.size(); ++f) {
      if (!m_data[f].m_independent) {
        std::string gradname =
            "diff/" + m_data[f].m_name + "/" + m_data[i].m_name;
        m_keys[gradname] = DiffClexParamGradKey(gradname, i, {m_data[f].lbegin},
                                                {m_data[f].m_cache.cols()});

        for (Index j = 0; j < m_data.size(); ++j) {
          if (m_data[j].m_independent) {
            std::string hessname = gradname + "/" + m_data[j].m_name;
            m_keys[hessname] = DiffClexParamHessKey(
                hessname, i, {m_data[f].lbegin, m_data[j].lbegin},
                {m_data[f].m_cache.cols(), m_data[j].m_cache.cols()});
          }
        }
      }
    }
  } else {
    m_data.push_back(
        DiffScalarContainer(_keyname, _rows, _cols, m_N_dependent, false));
    m_N_dependent += _rows * _cols;

    Index f = protokey.index();

    for (Index i = 0; i < m_data.size(); ++i) {
      if (m_data[i].m_independent) {
        std::string gradname =
            "diff/" + m_data[f].m_name + "/" + m_data[i].m_name;
        m_keys[gradname] = DiffClexParamGradKey(gradname, i, {m_data[f].lbegin},
                                                {m_data[f].m_cache.cols()});

        for (Index j = 0; j < m_data.size(); ++j) {
          if (m_data[j].m_independent) {
            std::string hessname = gradname + "/" + m_data[j].m_name;
            m_keys[hessname] = DiffClexParamHessKey(
                hessname, i, {m_data[f].lbegin, m_data[j].lbegin},
                {m_data[f].m_cache.cols(), m_data[j].m_cache.cols()});
          }
        }
      }
    }
  }

  m_eval.push_back(EvalMode::DEFAULT);

  return protokey;
}

}  // namespace clexulator

inline std::ostream &operator<<(
    std::ostream &sout, const clexulator::DiffClexParamPack::EvalMode &val) {
  sout << to_string<clexulator::DiffClexParamPack::EvalMode>(val);
  return sout;
}

inline std::istream &operator>>(std::istream &sin,
                                clexulator::DiffClexParamPack::EvalMode &val) {
  std::string s;
  sin >> s;
  val = from_string<clexulator::DiffClexParamPack::EvalMode>(s);
  return sin;
}

const std::string traits<clexulator::DiffClexParamPack::EvalMode>::name =
    "clex_eval_mode";

const std::multimap<clexulator::DiffClexParamPack::EvalMode,
                    std::vector<std::string> >
    traits<clexulator::DiffClexParamPack::EvalMode>::strval = {
        {clexulator::DiffClexParamPack::EvalMode::DEFAULT,
         {"Default", "DEFAULT", "default"}},
        {clexulator::DiffClexParamPack::EvalMode::DIFF,
         {"Diff", "DIFF", "diff", "Differential", "differential",
          "DIFFERENTIAL"}},
        {clexulator::DiffClexParamPack::EvalMode::READ,
         {"Read", "READ", "read"}},
        {clexulator::DiffClexParamPack::EvalMode::DYNAM,
         {"Dynamic", "Dynam", "DYNAMIC", "DYNAM", "dynamic", "dynam"}}};

namespace clexulator {

const DiffClexParamPack::EvalMode DiffClexParamPack::DEFAULT =
    DiffClexParamPack::EvalMode::DEFAULT;
const DiffClexParamPack::EvalMode DiffClexParamPack::DYNAM =
    DiffClexParamPack::EvalMode::DYNAM;
const DiffClexParamPack::EvalMode DiffClexParamPack::DIFF =
    DiffClexParamPack::EvalMode::DIFF;
const DiffClexParamPack::EvalMode DiffClexParamPack::READ =
    DiffClexParamPack::EvalMode::READ;

/** @} */
}  // namespace clexulator
}  // namespace CASM
#endif
