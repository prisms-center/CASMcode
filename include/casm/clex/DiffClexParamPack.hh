#ifndef DIFFCLEXPARAMPACK_HH
#define DIFFCLEXPARAMPACK_HH
#include <cstddef>
#include <iostream>
#include <string>
#include <sstream>

#include "casm/CASM_global_definitions.hh"
#include "casm/clex/ClexParamPack.hh"
#include "casm/casm_io/EnumIO.hh"
#include "casm/external/fadbad/fadiff.h"
#include "casm/external/fadbad/badiff.h"

namespace CASM {
  namespace DiffClexParamPack_impl {
    enum class EvalMode {
      DEFAULT, READ, DIFF, DYNAM

    };

    using DiffScalar = fadbad::B<fadbad::F<double> >;

  }
  class DiffClexParamPack;
  class DiffScalarContainer;
  class DiffClexParamKey;
  class DiffClexParamValKey;
  class DiffClexParamGradKey;
  class DiffClexParamHessKey;


  class DiffClexParamKey : public ClexParamPack_impl::BaseKey {
  public:
    typedef ClexParamPack::size_type size_type;

    DiffClexParamKey(std::string const &_name = "", bool _standalone = false, size_type _ind = -1) :
      ClexParamPack_impl::BaseKey(_name, _standalone),
      m_index(_ind) {}

    virtual ~DiffClexParamKey() {};

    virtual Eigen::MatrixXd const &eval(DiffScalarContainer const &_data, DiffClexParamPack_impl::EvalMode mode) const = 0;

    size_type index() const {
      return m_index;
    }

  private:
    size_type m_index;
  };

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



  struct DiffScalarContainer {
    using DiffScalar = DiffClexParamPack_impl::DiffScalar;

    std::string m_name;
    Index lbegin;
    bool m_independent;
    std::vector<std::vector<DiffScalar > > m_data;
    mutable Eigen::MatrixXd m_cache;

    DiffScalarContainer(std::string const &_name, Index param_dim, Index num_params, Index _lbegin, bool _independent) :
      m_name(_name),
      lbegin(_lbegin),
      m_independent(_independent),
      m_data(param_dim, std::vector<DiffScalar>(num_params)),
      m_cache(param_dim, num_params) {

    }

    // Gradient of dependent function index f wrt the independent variables contained in here
    // -> result(i,j) = Grad[F(f)](l(i,j)), where 'l' is linear index of each (i,j) entry in m_data.
    Eigen::MatrixXd const &grad(Index f) {
      for(Index i = 0; i < m_data.size(); ++i) {
        for(Index j = 0; j < m_data[i].size(); ++j) {
          m_cache(i, j) = m_data[i][j].d(f).val();
        }
      }
      return m_cache;
    }

    // Hessian components of dependent function index f corresponding to independent variable
    // index 'a' and all independent variables contained in here
    // -> result(i,j) = H[F(f)](a,l), where 'l' is linear index of each (i,j) entry in m_data.
    Eigen::MatrixXd const &hess(Index f, Index a) {
      Index l = lbegin;
      for(Index i = 0; i < m_data.size(); ++i) {
        for(Index j = 0; j < m_data[i].size(); ++j, ++l) {
          m_cache(i, j) = m_data[i][j].d(f).d(l);
        }
      }
      return m_cache;
    }

    void pre_eval(Index N) {
      Index l = lbegin;
      for(auto &dvec : m_data) {
        for(auto &datum : dvec) {
          datum.val().diff(l, N);
        }
      }
    }

    void eval() {
      for(Index i = 0; i < m_cache.rows(); ++i) {
        for(Index j = 0; j < m_cache.cols(); ++j) {
          m_cache(i, j) = m_data[i][j].val().val();
        }
      }
    }

    void post_eval(Index N) {
      Index l = lbegin;
      for(auto &dvec : m_data) {
        for(auto &datum : dvec) {
          datum.diff(l, N);
        }
      }
    }

    Index linear_index_end()const {
      return lbegin + m_cache.size();
    }
  };

  template<typename Scalar>
  struct ValAccess;

  /// \brief Abstract base class for reading/writing clexulator parameters
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

    template<typename Scalar>
    using Val = ValAccess<Scalar>;

    template<typename Scalar>
    friend class ValAccess;

    DiffClexParamPack() :
      m_tot_eval_mode(DEFAULT),
      m_N_independent(0),
      m_N_dependent(0) {}

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

    EvalMode eval_mode(BaseKey const &_key) const {
      return m_eval[_key.index()];
    }

    EvalMode eval_mode() const {
      return m_tot_eval_mode;
    }

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

    double const &read(ClexParamKey const &_key, size_type _i, size_type _j) const override {
      return read(*static_cast<BaseKey const *>(_key.ptr()), _i, _j);
    }

    double const &read(BaseKey const &_key, size_type _i, size_type _j) const {
      return _key.eval(m_data[_key.index()], eval_mode())(_i, _j);
    }

    void write(ClexParamKey const &_key, Eigen::Ref<const Eigen::MatrixXd> const &_val) override {
      write(*static_cast<BaseKey const *>(_key.ptr()), _val);
    }

    void write(BaseKey const &_key, Eigen::Ref<const Eigen::MatrixXd> const &_val);

    void write(ClexParamKey const &_key, size_type _i, double _val) override {
      write(*static_cast<BaseKey const *>(_key.ptr()), _i, _val);
    }

    void write(BaseKey const &_key, size_type _i, double _val);

    void write(ClexParamKey const &_key, size_type _i, size_type _j, double _val) override {
      write(*static_cast<BaseKey const *>(_key.ptr()), _i, _j, _val);
    }

    void write(BaseKey const &_key, size_type _i, size_type _j, double _val);


    Key allocate(std::string const &_keyname, Index _rows, Index _cols, bool _independent);

    void set_eval_mode(ClexParamKey const &_key, std::string const &_mode) override {
      set_eval_mode(*static_cast<BaseKey const *>(_key.ptr()), from_string<EvalMode>(_mode));
    }

    void set_eval_mode(BaseKey const &_key, EvalMode _mode) {
      m_eval[_key.index()] = _mode;
      m_tot_eval_mode = DEFAULT;
      for(auto const &eval : m_eval) {
        if(eval == DIFF) {
          m_tot_eval_mode = DIFF;
          break;
        }
      }
    }

    void pre_eval() {
      if(eval_mode() == DIFF) {
        for(DiffScalarContainer &datum : m_data) {
          if(datum.m_independent) {
            datum.pre_eval(m_N_independent);
          }
        }
      }
    }

    void post_eval() {
      if(eval_mode() == DIFF) {
        for(DiffScalarContainer &datum : m_data) {
          if(!datum.m_independent) {
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

    DiffClexParamValKey(std::string const &_name = "", size_type _ind = -1) :
      DiffClexParamKey(_name, _ind, true) {}

    Eigen::MatrixXd const &eval(DiffScalarContainer const &_data, DiffClexParamPack::EvalMode mode) const override {
      if(mode == DiffClexParamPack::DIFF) {
        for(Index i = 0; i < _data.m_cache.rows(); ++i) {
          for(Index j = 0; j < _data.m_cache.cols(); ++j) {
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

    size_type m_index;
  };

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  class DiffClexParamGradKey : public DiffClexParamKey {
  public:
    typedef ClexParamPack::size_type size_type;

    DiffClexParamGradKey(std::string const &_name = "", size_type _if = -1, size_type _ix = -1) :
      DiffClexParamKey(_name, false, _ix),
      m_if(_if),
      m_ifcompon(-1) {}

    Eigen::MatrixXd const &eval(DiffScalarContainer const &_data, DiffClexParamPack::EvalMode mode) const override {
      if(mode != DiffClexParamPack::DIFF)
        throw std::runtime_error("Requested Gradient parameter " + name() + " for incompatible evaluation mode " + to_string(mode) + ".");
      //Do stuff
      throw std::runtime_error("DiffClexParamGradKey::eval() is not implemented!");
      return _data.m_cache;
    }

  protected:

    ClexParamPack_impl::BaseKey *_clone() const override {
      return new DiffClexParamGradKey(*this);
    }

    const size_type m_if;
    mutable size_type m_ifcompon;
  };

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  class DiffClexParamHessKey : public DiffClexParamKey {
  public:
    typedef ClexParamPack::size_type size_type;


    DiffClexParamHessKey(std::string const &_name = "", size_type _if = -1, size_type _ix = -1, size_type _iy = -1) :
      DiffClexParamKey(_name, false, _ix),
      m_if(_if),
      m_iy(_iy),
      m_ifcompon(-1),
      m_iycompon(-1) {}

    Eigen::MatrixXd const &eval(DiffScalarContainer const &_data, DiffClexParamPack::EvalMode mode) const override {
      if(mode != DiffClexParamPack::DIFF)
        throw std::runtime_error("Requested Hessian parameter " + name() + " for incompatible evaluation mode " + to_string(mode) + ".");
      //Do stuff
      throw std::runtime_error("DiffClexParamGradKey::eval() is not implemented!");
      return _data.m_cache;
    }

  protected:

    ClexParamPack_impl::BaseKey *_clone() const override {
      return new DiffClexParamHessKey(*this);
    }

    const size_type m_if, m_iy;
    mutable size_type m_ifcompon, m_iycompon;
  };

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



  template<>
  struct traits<DiffClexParamPack::EvalMode> {

    static const std::string name;

    static const std::multimap<DiffClexParamPack::EvalMode, std::vector<std::string> > strval;

  };

  template<typename Scalar>
  struct ValAccess { };

  template<>
  struct ValAccess<double> {
    using size_type = DiffClexParamPack::size_type;

    static double const &get(DiffClexParamPack const &_pack, DiffClexParamKey const &_key, size_type i) {
      return _pack.m_data[_key.index()].m_cache(i, 0);
    }

    static double const &get(DiffClexParamPack const &_pack, DiffClexParamKey const &_key, size_type i, size_type j) {
      return _pack.m_data[_key.index()].m_cache(i, j);
    }

    static void set(DiffClexParamPack &_pack, DiffClexParamKey const &_key, Eigen::Ref<const Eigen::MatrixXd> const &_val) {
      _pack.m_data[_key.index()].m_cache = _val;
    }

    template<typename Scalar2>
    static void set(DiffClexParamPack &_pack, DiffClexParamKey const &_key, size_type i, Scalar2 const &_val) {
      _pack.m_data[_key.index()].m_cache(i, 0) = _val;
    }

    template<typename Scalar2>
    static  void set(DiffClexParamPack &_pack, DiffClexParamKey const &_key, size_type i, size_type j, Scalar2 const &_val) {
      _pack.m_data[_key.index()].m_cache(i, j) = _val;
    }
  };

  template<>
  struct ValAccess<typename DiffClexParamPack::DiffScalar> {
    using size_type = DiffClexParamPack::size_type;

    static typename DiffClexParamPack::DiffScalar const &get(DiffClexParamPack const &_pack, DiffClexParamKey const &_key, size_type i) {
      return _pack.m_data[_key.index()].m_data[i][0];
    }

    static typename DiffClexParamPack::DiffScalar const &get(DiffClexParamPack const &_pack, DiffClexParamKey const &_key, size_type i, size_type j) {
      return _pack.m_data[_key.index()].m_data[i][j];
    }

    static void set(DiffClexParamPack &_pack, DiffClexParamKey const &_key, Eigen::Ref<const Eigen::MatrixXd> const &_val) {
      assert(_val.rows() == _pack.m_data[_key.index()].m_data.size());
      for(Index i = 0; i < _val.rows(); ++i) {
        assert(_val.cols() == _pack.m_data[_key.index()].m_data[i].size());
        for(Index j = 0; j < _val.cols(); ++j) {
          _pack.m_data[_key.index()].m_data[i][j] = _val(i, j);
        }
      }
    }

    template<typename Scalar2>
    static void set(DiffClexParamPack &_pack, DiffClexParamKey const &_key, size_type i, Scalar2 const &_val) {
      _pack.m_data[_key.index()].m_data[i][0] = _val;
    }

    template<typename Scalar2>
    static  void set(DiffClexParamPack &_pack, DiffClexParamKey const &_key, size_type i, size_type j, Scalar2 const &_val) {
      _pack.m_data[_key.index()].m_data[i][j] = _val;
    }
  };



  inline
  void DiffClexParamPack::write(BaseKey const &_key, Eigen::Ref<const Eigen::MatrixXd> const &_val) {
    if(!_key.standalone()) {
      throw std::runtime_error("Cannot write to dependent parameter " + _key.name() + ".");
    }
    if(eval_mode() == DEFAULT) {
      Val<double>::set(*this, _key, _val);
    }
    else {
      Val<DiffScalar>::set(*this, _key, _val);
    }
  }

  inline
  void DiffClexParamPack::write(BaseKey const &_key, size_type _i, double _val) {
    if(!_key.standalone()) {
      throw std::runtime_error("Cannot write to dependent parameter " + _key.name() + ".");
    }

    if(eval_mode() == DEFAULT) {
      Val<double>::set(*this, _key, _i, _val);
    }
    else {
      Val<DiffScalar>::set(*this, _key, _i, _val);
    }

  }

  inline
  void DiffClexParamPack::write(BaseKey const &_key, size_type _i, size_type _j, double _val) {
    if(!_key.standalone()) {
      throw std::runtime_error("Cannot write to dependent parameter " + _key.name() + ".");
    }

    if(eval_mode() == DEFAULT) {
      Val<double>::set(*this, _key, _i, _j, _val);
    }
    else {
      Val<DiffScalar>::set(*this, _key, _i, _j, _val);
    }

  }

  inline
  DiffClexParamPack::Key DiffClexParamPack::allocate(std::string const &_keyname, Index _rows, Index _cols, bool _independent) {
    auto it = keys().find(_keyname);
    if(it != keys().end())
      throw std::runtime_error("Naming collision in DiffClexParamPack::allocate(), ClexParamPack already managing parameter allocation corresponding to name " + _keyname + ".");
    Key protokey(_keyname, m_data.size());
    m_keys[_keyname] = protokey;

    if(_independent) {
      m_data.push_back(DiffScalarContainer(_keyname, _rows, _cols, m_N_independent, true));
      m_N_independent += _rows * _cols;
      Index i = protokey.index();
      for(Index f = 0; f < m_data.size(); ++f) {
        if(!m_data[f].m_independent) {
          std::string gradname = "diff/" + m_data[f].m_name + "/" + m_data[i].m_name;
          m_keys[gradname] = DiffClexParamGradKey(gradname, f, i);

          for(Index j = 0; j < m_data.size(); ++j) {
            if(m_data[j].m_independent) {
              std::string hessname = gradname + "/" + m_data[j].m_name;
              m_keys[gradname] = DiffClexParamHessKey(gradname, f, i, j);
            }
          }
        }
      }
    }
    else {
      m_data.push_back(DiffScalarContainer(_keyname, _rows, _cols, m_N_dependent, false));
      m_N_dependent += _rows * _cols;

      Index f = protokey.index();
      for(Index i = 0; i < m_data.size(); ++i) {
        if(m_data[i].m_independent) {
          std::string gradname = "diff/" + m_data[f].m_name + "/" + m_data[i].m_name;
          m_keys[gradname] = DiffClexParamGradKey(gradname, f, i);

          for(Index j = 0; j < m_data.size(); ++j) {
            if(m_data[j].m_independent) {
              std::string hessname = gradname + "/" + m_data[j].m_name;
              m_keys[gradname] = DiffClexParamHessKey(gradname, f, i, j);
            }
          }
        }
      }

    }

    m_eval.push_back(EvalMode::DEFAULT);

    return protokey;
  }

  inline
  std::ostream &operator<<(std::ostream &sout, const DiffClexParamPack::EvalMode &val) {
    sout << to_string<DiffClexParamPack::EvalMode>(val);
    return sout;
  }

  inline
  std::istream &operator>>(std::istream &sin, DiffClexParamPack::EvalMode &val) {
    std::string s;
    sin >> s;
    val = from_string<DiffClexParamPack::EvalMode>(s);
    return sin;
  }


  const std::string traits<DiffClexParamPack::EvalMode>::name = "clex_eval_mode";

  const std::multimap<DiffClexParamPack::EvalMode, std::vector<std::string> > traits<DiffClexParamPack::EvalMode>::strval = {
    {DiffClexParamPack::EvalMode::DEFAULT, {"Default", "DEFAULT", "default"} },
    {DiffClexParamPack::EvalMode::DIFF, {"Diff", "DIFF", "diff", "Differential", "differential", "DIFFERENTIAL"} },
    {DiffClexParamPack::EvalMode::READ, {"Read", "READ", "read"} },
    {DiffClexParamPack::EvalMode::DYNAM, {"Dynamic", "Dynam", "DYNAMIC", "DYNAM", "dynamic", "dynam"} }
  };

  const DiffClexParamPack::EvalMode DiffClexParamPack::DEFAULT = DiffClexParamPack::EvalMode::DEFAULT;
  const DiffClexParamPack::EvalMode DiffClexParamPack::DYNAM = DiffClexParamPack::EvalMode::DYNAM;
  const DiffClexParamPack::EvalMode DiffClexParamPack::DIFF = DiffClexParamPack::EvalMode::DIFF;
  const DiffClexParamPack::EvalMode DiffClexParamPack::READ = DiffClexParamPack::EvalMode::READ;



}
#endif
