#include "casm/kinetics/DiffTransConfigIsEqual.hh"
#include "casm/kinetics/DiffTransConfiguration_impl.hh"

namespace {
  using namespace CASM;
  ScelPeriodicDiffTransSymCompare _construct_scel_sym_compare(const Supercell &scel) {
    return ScelPeriodicDiffTransSymCompare(scel.primclex().shared_prim(), scel.prim_grid(), scel.crystallography_tol());
  }
}

namespace CASM {
  namespace Kinetics {

    // --- DiffTransConfigCompareFast implementation helpers ---

    namespace DiffTransConfigIsEqualFastImpl {

      Untransformed make(DiffTransConfigIsEqualFast *_home);

      Untransformed make(DiffTransConfigIsEqualFast *_home, const DiffTransConfiguration &_other);

      Transformed make(DiffTransConfigIsEqualFast *_home, const PermuteIterator &_op);

      Transformed make(DiffTransConfigIsEqualFast *_home, const PermuteIterator &_op, const DiffTransConfiguration &_other);

      /// This implements:
      ///   A*dtconfig ?= B*other
      /// where,
      ///   A=first.op, dtconfig=first.dtconfig,
      ///   B=second.op, other=second.dtconfig
      /// by checking:
      ///   min(A*dtconfig.from_config, A*dtconfig.to_config) ?= min(B*other.from_config, B*other.to_config),
      ///
      template<typename First, typename Second>
      bool _check(const First &first, const Second &second) {
        if(first.is_sorted()) {
          if(second.is_sorted()) {
            return first.check_from_from(second);
          }
          else {
            return first.check_from_to(second);
          }
        }
        else {
          if(second.is_sorted()) {
            return first.check_to_from(second);
          }
          else {
            return first.check_to_to(second);
          }
        }
      }

      /// Base for Untransformed & Transformed
      template<typename Base>
      struct Common : public Base {

        using Base::derived;
        typedef DiffTransConfigIsEqualFast Home;

        Common(const Home *_home) :
          home(_home), dtconfig(home->m_dtconfig) {}

        Common(const Home *_home, const DiffTransConfiguration &other) :
          home(_home), dtconfig(&other) {}

        const Home *home;
        const DiffTransConfiguration *dtconfig;

        const DiffusionTransformation &diff_trans() const {
          return dtconfig->diff_trans();
        }

        const Configuration &from_config() const {
          return dtconfig->from_config();
        }

        const Configuration &to_config() const {
          return dtconfig->to_config();
        }

        /// Checks if this->from_config() == other.from_config(), including ops
        template<typename Other>
        bool check_from_from(const Other &other) const {
          return derived()._check(from_config(), other, other.from_config());
        }

        /// Checks if this->from_config() == other.to_config(), including ops
        template<typename Other>
        bool check_from_to(const Other &other) const {
          return derived()._check(from_config(), other, other.to_config());
        }

        /// Checks if this->to_config() == other.from_config(), including ops
        template<typename Other>
        bool check_to_from(const Other &other) const {
          return derived()._check(to_config(), other, other.from_config());
        }

        /// Checks if this->to_config() == other.to_config(), including ops
        template<typename Other>
        bool check_to_to(const Other &other) const {
          return derived()._check(to_config(), other, other.to_config());
        }
      };

      /// Represents 'op*dtconfig'
      struct Transformed : public Common<CRTPBase<Transformed>> {

        typedef DiffTransConfigIsEqualFast Home;

        Transformed(const Home *_home, const PermuteIterator &_op) :
          Common(_home), op(_op) {}

        Transformed(const Home *_home, const PermuteIterator &_op, const DiffTransConfiguration &_other) :
          Common(_home, _other), op(_op) {}

        const PermuteIterator &op;

        bool is_sorted() const {
          auto f = to_config().less();
          return !f(op, from_config());
        }

        bool _check(const Configuration &_first, const Transformed &other, const Configuration &_second) const {
          // compare diff_trans first
          auto this_diff_trans = home->m_sym_compare.prepare(copy_apply(op, diff_trans()));
          auto other_diff_trans = home->m_sym_compare.prepare(copy_apply(other.op, other.diff_trans()));
          if(this_diff_trans < other_diff_trans) {
            home->set_is_less(true);
            return false;
          }
          else if(other_diff_trans < this_diff_trans) {
            home->set_is_less(false);
            return false;
          }

          // then compare configs
          auto f = _first.equal_to();
          if(!f(op, other.op, _second)) {
            home->set_is_less(f.is_less());
            return false;
          }
          return true;
        }
      };

      /// Represents 'dtconfig'
      struct Untransformed : public Common<CRTPBase<Untransformed>> {

        typedef DiffTransConfigIsEqualFast Home;

        Untransformed(const Home *_home) :
          Common(_home) {}

        Untransformed(const Home *_home, const DiffTransConfiguration &other) :
          Common(_home, other) {}

        bool is_sorted() const {
          return dtconfig->is_sorted();
        }

        bool _check(const Configuration &_first, const Transformed &other, const Configuration &_second) const {
          // compare diff_trans first
          auto other_diff_trans = home->m_sym_compare.prepare(copy_apply(other.op, other.diff_trans()));
          if(diff_trans() < other_diff_trans) {
            home->set_is_less(true);
            return false;
          }
          else if(other_diff_trans < diff_trans()) {
            home->set_is_less(false);
            return false;
          }

          // then compare configs
          auto f = _first.equal_to();
          if(!f(other.op, _second)) {
            home->set_is_less(f.is_less());
            return false;
          }
          return true;
        }

        bool _check(const Configuration &_first, const Untransformed &other, const Configuration &_second) const {
          // compare diff_trans first
          if(diff_trans() < other.diff_trans()) {
            home->set_is_less(true);
            return false;
          }
          else if(other.diff_trans() < diff_trans()) {
            home->set_is_less(false);
            return false;
          }

          // then compare configs
          auto f = _first.equal_to();
          if(!f(_second)) {
            home->set_is_less(f.is_less());
            return false;
          }
          return true;
        }
      };


      Untransformed make(const DiffTransConfigIsEqualFast *_home) {
        return Untransformed(_home);
      }

      Untransformed make(const DiffTransConfigIsEqualFast *_home, const DiffTransConfiguration &_other) {
        return Untransformed(_home, _other);
      }

      Transformed make(const DiffTransConfigIsEqualFast *_home, const PermuteIterator &_op) {
        return Transformed(_home, _op);
      }

      Transformed make(const DiffTransConfigIsEqualFast *_home, const PermuteIterator &_op, const DiffTransConfiguration &_other) {
        return Transformed(_home, _op, _other);
      }

    }


    // --- DiffTransConfigCompareFast ---

    DiffTransConfigIsEqualFast::DiffTransConfigIsEqualFast(const DiffTransConfiguration &_dtconfig) :
      m_dtconfig(&_dtconfig),
      m_sym_compare(::_construct_scel_sym_compare(_dtconfig.supercell())) {}

    /// \brief Check if dtconfig < other (may have different Supercell, may not be sorted)
    bool DiffTransConfigIsEqualFast::operator()(const DiffTransConfiguration &other) const {
      using DiffTransConfigIsEqualFastImpl::make;
      return _check(make(this), make(this, other));
    }

    /// \brief Check if dtconfig < A*dtconfig
    bool DiffTransConfigIsEqualFast::operator()(const PermuteIterator &A) const {
      using DiffTransConfigIsEqualFastImpl::make;
      return _check(make(this), make(this, A, *m_dtconfig));
    }

    /// \brief Check if A*dtconfig < B*dtconfig
    bool DiffTransConfigIsEqualFast::operator()(const PermuteIterator &A, const PermuteIterator &B) const {
      using DiffTransConfigIsEqualFastImpl::make;
      return _check(make(this, A, *m_dtconfig), make(this, B, *m_dtconfig));
    }

    /// \brief Check if dtconfig < A*other
    bool DiffTransConfigIsEqualFast::operator()(const PermuteIterator &A, const DiffTransConfiguration &other) const {
      using DiffTransConfigIsEqualFastImpl::make;
      return _check(make(this), make(this, A, other));
    }

    /// \brief Check if A*dtconfig < B*other
    bool DiffTransConfigIsEqualFast::operator()(const PermuteIterator &A, const PermuteIterator &B, const DiffTransConfiguration &other) const {
      using DiffTransConfigIsEqualFastImpl::make;
      return _check(make(this, A, *m_dtconfig), make(this, B, other));
    }

    bool DiffTransConfigIsEqualFast::is_less() const {
      return m_is_less;
    }

    /// For use by friends
    void DiffTransConfigIsEqualFast::set_is_less(bool value) const {
      m_is_less = value;
    }



    // --- DiffTransConfigIsEqualSimple ---

    DiffTransConfigIsEqualSimple::DiffTransConfigIsEqualSimple(const DiffTransConfiguration &_dtconfig) :
      m_dtconfig(&_dtconfig) {}

    /// \brief Check if dtconfig == other (may have different Supercell, may not be sorted)
    bool DiffTransConfigIsEqualSimple::operator()(const DiffTransConfiguration &other) const {

      // compare diff_trans (diff_trans from DiffTransConfig are always prepared)
      if(m_dtconfig->diff_trans() < other.diff_trans()) {
        m_is_less = true;
        return false;
      }
      else if(other.diff_trans() < m_dtconfig->diff_trans()) {
        m_is_less = false;
        return false;
      }

      auto f = m_dtconfig->sorted().from_config().equal_to();
      if(!f(other.sorted().from_config())) {
        m_is_less = f.is_less();
        return false;
      }
      return true;
    }

    /// \brief Check if dtconfig == A*dtconfig
    bool DiffTransConfigIsEqualSimple::operator()(const PermuteIterator &A) const {
      DiffTransConfiguration tmp {copy_apply(A, *m_dtconfig)};
      return (*this)(tmp);
    }

    /// \brief Check if A*dtconfig == B*dtconfig
    bool DiffTransConfigIsEqualSimple::operator()(const PermuteIterator &A, const PermuteIterator &B) const {
      DiffTransConfiguration tmp {copy_apply(A, *m_dtconfig)};
      DiffTransConfigIsEqualSimple f(tmp);
      if(!f(copy_apply(B, *m_dtconfig))) {
        m_is_less = f.is_less();
        return false;
      }
      return true;
    }

    /// \brief Check if dtconfig == A*other
    bool DiffTransConfigIsEqualSimple::operator()(
      const PermuteIterator &A,
      const DiffTransConfiguration &other) const {
      DiffTransConfiguration tmp {copy_apply(A, other)};
      return (*this)(tmp);
    }

    /// \brief Check if A*dtconfig == B*other
    bool DiffTransConfigIsEqualSimple::operator()(
      const PermuteIterator &A,
      const PermuteIterator &B,
      const DiffTransConfiguration &other) const {
      DiffTransConfiguration tmp {copy_apply(A, *m_dtconfig)};
      DiffTransConfigIsEqualSimple f(tmp);
      if(!f(copy_apply(B, other))) {
        m_is_less = f.is_less();
        return false;
      }
      return true;
    }

    bool DiffTransConfigIsEqualSimple::is_less() const {
      return m_is_less;
    }

  }
}
