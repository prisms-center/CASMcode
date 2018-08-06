#ifndef CASM_ConfigDoFIsEquivalent
#define CASM_ConfigDoFIsEquivalent

#include "casm/casm_io/Log.hh"
#include "casm/clex/ConfigDoF.hh"
#include "casm/clex/Configuration.hh"

namespace CASM {

  /** \ingroup ConfigIsEquivalent
   *  @{
   */

  /// Namespace containing DoF comparison functors
  namespace DoFIsEquivalent {

    /// \brief Base class for functors that compare ConfigDoF
    ///
    /// - Specialized for Integral or Float DoF, and then for type
    ///   (strain, occupation, displacement etc.)
    /// - The derived call operators return the value for equality comparison,
    ///   and if not equivalent, also store the result for less than comparison
    class ConfigDoFIsEquivalentBase {

    public:

      ConfigDoFIsEquivalentBase(const ConfigDoF &_configdof) :
        m_configdof(&_configdof) {}

      const ConfigDoF &configdof() const {
        return *m_configdof;
      }

      Index size() const {
        return configdof().size();
      }

      /// \brief Returns less than comparison
      ///
      /// - Only valid after call operator returns false
      bool is_less() const {
        return m_less;
      }

      /// \brief Return config == other
      bool operator()(const Configuration &other) const {
        return (*this)(other.configdof());
      }

      /// \brief Return config == other
      virtual bool operator()(const ConfigDoF &other) const = 0;

      /// \brief Return config == A*config
      virtual bool operator()(const PermuteIterator &A) const = 0;

      /// \brief Return A*config == B*config
      virtual bool operator()(const PermuteIterator &A, const PermuteIterator &B) const = 0;

      std::unique_ptr<ConfigDoFIsEquivalentBase> clone() const {
        return std::unique_ptr<ConfigDoFIsEquivalentBase>(this->_clone());
      }

    protected:

      mutable bool m_less;

    private:

      virtual ConfigDoFIsEquivalentBase *_clone() const = 0;

      const ConfigDoF *m_configdof;

    };

    /// \brief Abstract base class specialization of ConfigDoFIsEquivalentBase for
    /// integral DoF types
    ///
    /// - The protected '_check' method provides for both checking equality and if
    ///   not equivalent, storing the 'less than' result
    class IntegralIsEquivalent: public ConfigDoFIsEquivalentBase {

    public:
      IntegralIsEquivalent(const ConfigDoF &_configdof) :
        ConfigDoFIsEquivalentBase(_configdof) {}

    protected:
      template<typename T>
      bool _check(const T &A, const T &B) const {
        if(A == B) {
          return true;
        }
        m_less = (A < B);
        return false;
      }
    };

    /// \brief Abstract base class specialization of ConfigDoFIsEquivalentBase for
    /// floating point DoF types
    ///
    /// - The protected '_check' method provides for both checking equality and if
    ///   not equivalent, storing the 'less than' result
    class FloatIsEquivalent: public ConfigDoFIsEquivalentBase {

    public:
      FloatIsEquivalent(const ConfigDoF &_configdof, double _tol) :
        ConfigDoFIsEquivalentBase(_configdof), m_tol(_tol) {}

    protected:
      template<typename T>
      bool _check(const T &A, const T &B) const {
        if(A < B - tol()) {
          m_less = true;
          return false;
        }
        if(A > B + tol()) {
          m_less = false;
          return false;
        }
        return true;
      }

    private:

      double tol() const {
        return m_tol;
      }

      double m_tol;
    };


    /// Compare occupation DoF
    class Occupation : public IntegralIsEquivalent {

    public:

      Occupation(const ConfigDoF &_configdof) :
        IntegralIsEquivalent(_configdof) {}

      Occupation(const Configuration &_config) :
        Occupation(_config.configdof()) {}

      /// \brief Return config == other, store config < other
      bool operator()(const ConfigDoF &other) const override {
        return _for_each(
        [&](Index i) {
          return this->configdof().occ(i);
        },
        [&](Index i) {
          return other.occ(i);
        });
      }

      /// \brief Return config == A*config, store config < A*config
      bool operator()(const PermuteIterator &A) const override {
        return _for_each(
        [&](Index i) {
          return this->configdof().occ(i);
        },
        [&](Index i) {
          return this->configdof().occ(A.permute_ind(i));
        });
      }

      /// \brief Return A*config == B*config, store A*config < B*config
      bool operator()(const PermuteIterator &A, const PermuteIterator &B) const override {
        return _for_each(
        [&](Index i) {
          return this->configdof().occ(A.permute_ind(i));
        },
        [&](Index i) {
          return this->configdof().occ(B.permute_ind(i));
        });
      }

      std::unique_ptr<Occupation> clone() const {
        return std::unique_ptr<Occupation>(this->_clone());
      }

    private:

      template<typename F, typename G>
      bool _for_each(F f, G g) const {
        Index i;
        for(i = 0; i < size(); i++) {
          if(!_check(f(i), g(i))) {
            return false;
          }
        }
        return true;
      }

      Occupation *_clone() const override {
        return new Occupation(*this);
      }
    };


    /// Compare displacement DoF
    class Displacement : public FloatIsEquivalent {

    public:

      Displacement(const ConfigDoF &_configdof, double _tol) :
        FloatIsEquivalent(_configdof, _tol),
        m_fg_index_A(-1),
        m_fg_index_B(-1) {}

      Displacement(const Configuration &_config, double _tol) :
        Displacement(_config.configdof(), _tol) {}

      /// \brief Return config == other, store config < other
      bool operator()(const ConfigDoF &other) const override {
        return _for_each(
        [&](Index i, Index j) {
          return this->configdof().disp(i)[j];
        },
        [&](Index i, Index j) {
          return other.disp(i)[j];
        });
      }

      /// \brief Return config == A*config, store config < A*config
      bool operator()(const PermuteIterator &A) const override {
        _update_A(A);
        return _for_each(
        [&](Index i, Index j) {
          return this->configdof().disp(i)[j];
        },
        [&](Index i, Index j) {
          return this->new_disp_A(j, A.permute_ind(i));
        });
      }

      /// \brief Return A*config == B*config, store A*config < B*config
      bool operator()(const PermuteIterator &A, const PermuteIterator &B) const override {
        _update_A(A);
        _update_B(B);
        return _for_each(
        [&](Index i, Index j) {
          return this->new_disp_A(j, A.permute_ind(i));
        },
        [&](Index i, Index j) {
          return this->new_disp_B(j, B.permute_ind(i));
        });
      }

      std::unique_ptr<Displacement> clone() const {
        return std::unique_ptr<Displacement>(this->_clone());
      }

    private:

      void _update_A(const PermuteIterator &A) const {
        if(A.factor_group_index() != m_fg_index_A) {
          m_fg_index_A = A.factor_group_index();
          m_new_disp_A = A.sym_op().matrix() * configdof().displacement();
        }
      }

      void _update_B(const PermuteIterator &B) const {
        if(B.factor_group_index() != m_fg_index_B) {
          m_fg_index_B = B.factor_group_index();
          m_new_disp_B = B.sym_op().matrix() * configdof().displacement();
        }
      }

      double new_disp_A(Index i, Index j) const {
        return m_new_disp_A(i, j);
      }

      double new_disp_B(Index i, Index j) const {
        return m_new_disp_B(i, j);
      }

      template<typename F, typename G>
      bool _for_each(F f, G g) const {
        Index i, j;
        for(i = 0; i < size(); i++) {
          for(j = 0; j < 3; j++) {
            if(!_check(f(i, j), g(i, j))) {
              return false;
            }
          }
        }
        return true;
      }

      Displacement *_clone() const override {
        return new Displacement(*this);
      }

      mutable Index m_fg_index_A;
      mutable Eigen::MatrixXd m_new_disp_A;

      mutable Index m_fg_index_B;
      mutable Eigen::MatrixXd m_new_disp_B;
    };


    /// Compare strain DoF
    ///
    /// - Compares F.t * F, unrolled, lexicographically
    class Strain : public FloatIsEquivalent {

    public:

      Strain(const ConfigDoF &_configdof, double _tol) :
        FloatIsEquivalent(_configdof, _tol),
        m_def_tensor(_configdof.deformation().transpose() * _configdof.deformation()),
        m_fg_index_A(-1),
        m_fg_index_B(-1) {}

      Strain(const Configuration &_config, double _tol) :
        Strain(_config.configdof(), _tol) {}

      /// \brief Return config == other, store config < other
      bool operator()(const ConfigDoF &other) const override {
        Eigen::MatrixXd other_def_tensor = other.deformation().transpose() * other.deformation();
        return _for_each(
        [&](Index i, Index j) {
          return this->_def_tensor(i, j);
        },
        [&](Index i, Index j) {
          return other_def_tensor(i, j);
        });
      }

      /// \brief Return config == A*config, store config < A*config
      bool operator()(const PermuteIterator &A) const override {
        _update_A(A);
        return _for_each(
        [&](Index i, Index j) {
          return this->_def_tensor(i, j);
        },
        [&](Index i, Index j) {
          return this->_def_tensor_A(i, j);
        });
      }

      /// \brief Return A*config == B*config, store A*config < B*config
      bool operator()(const PermuteIterator &A, const PermuteIterator &B) const override {
        _update_A(A);
        _update_B(B);
        return _for_each(
        [&](Index i, Index j) {
          return this->_def_tensor_A(i, j);
        },
        [&](Index i, Index j) {
          return this->_def_tensor_B(i, j);
        });
      }

      std::unique_ptr<Strain> clone() const {
        return std::unique_ptr<Strain>(this->_clone());
      }

    private:

      Strain *_clone() const override {
        return new Strain(*this);
      }

      void _update_A(const PermuteIterator &A) const {
        if(A.factor_group_index() != m_fg_index_A) {
          m_fg_index_A = A.factor_group_index();
          m_def_tensor_A = A.sym_op().matrix() * m_def_tensor * A.sym_op().matrix().transpose();
        }
      }

      void _update_B(const PermuteIterator &B) const {
        if(B.factor_group_index() != m_fg_index_B) {
          m_fg_index_B = B.factor_group_index();
          m_def_tensor_B = B.sym_op().matrix() * m_def_tensor * B.sym_op().matrix().transpose();
        }
      }

      double _def_tensor(Index i, Index j) const {
        return m_def_tensor(i, j);
      }

      double _def_tensor_A(Index i, Index j) const {
        return m_def_tensor_A(i, j);
      }

      double _def_tensor_B(Index i, Index j) const {
        return m_def_tensor_B(i, j);
      }

      template<typename F, typename G>
      bool _for_each(F f, G g) const {
        Index i, j;
        for(i = 0; i < 3; i++) {
          for(j = 0; j < 3; j++) {
            if(!_check(f(i, j), g(i, j))) {
              return false;
            }
          }
        }
        return true;
      }

      /// config().deformation().transpose()*config().deformation()
      mutable Eigen::MatrixXd m_def_tensor;

      mutable Index m_fg_index_A;
      mutable Eigen::MatrixXd m_def_tensor_A;

      mutable Index m_fg_index_B;
      mutable Eigen::MatrixXd m_def_tensor_B;
    };

  }// end namespace DoF


  /// \brief Wrapper class for generic equality comparison of ConfigDoF
  ///
  /// - Wraps a functor derived from ConfigDoFIsEquivalentBase that is specialized for
  ///   comparison of a particular type of DoF
  class ConfigDoFIsEquivalent {

  public:

    /// \brief Construct a ConfigDoFCompare object for a particular DoF type
    ///
    /// Easiest construction is probably using 'make_dof_compare'.
    ///
    /// Example:
    /// \code
    /// ConfigDoFIsEquivalent strain_is_equivalent = make_dof_is_equivalent<DoFIsEquivalent::Strain>(my_configdof or my_config);
    /// ConfigDoFIsEquivalent occ_is_equivalent = make_dof_is_equivalent<DoFIsEquivalent::Occupation>(my_configdof or my_config);
    /// ConfigDoFIsEquivalent disp_is_equivalent = make_dof_is_equivalent<DoFIsEquivalent::Displacement>(my_configdof or my_config);
    /// \endcode
    template<typename ConfigDoFIsEquivalentType>
    ConfigDoFIsEquivalent(std::unique_ptr<ConfigDoFIsEquivalentType> f) :
      m_f(f) {}

    /// \brief Returns less than comparison
    ///
    /// - Only valid after call operator returns false
    bool is_less() const {
      return m_f->is_less();
    }

    /// \brief Return config == other
    bool operator()(const Configuration &other) const {
      return (*this)(other.configdof());
    }

    /// \brief Return config == other
    bool operator()(const ConfigDoF &other) const {
      return (*m_f)(other);
    }

    /// \brief Return config == A*config
    bool operator()(const PermuteIterator &A) const {
      return (*m_f)(A);
    }

    /// \brief Return A*config == B*config
    bool operator()(const PermuteIterator &A, const PermuteIterator &B) const {
      return (*m_f)(A, B);
    }

  private:
    notstd::cloneable_ptr<DoFIsEquivalent::ConfigDoFIsEquivalentBase> m_f;

  };

  /// Factory function to make ConfigDoFIsEquivalent
  ///
  /// \relates ConfigDoFIsEquivalent
  template<typename ConfigDoFIsEquivalentType, typename ...Args>
  ConfigDoFIsEquivalent make_dof_is_equivalent(Args &&...args) {
    return ConfigDoFIsEquivalent(notstd::make_unique<ConfigDoFIsEquivalentType>(std::forward<Args>(args)...));
  }

  /** @} */
}

#endif
