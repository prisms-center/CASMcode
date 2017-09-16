#ifndef CASM_DiffTransConfigIsEqual
#define CASM_DiffTransConfigIsEqual

#include "casm/clusterography/ClusterSymCompare.hh"
#include "casm/kinetics/DiffusionTransformationTraits.hh"

namespace CASM {

  class PermuteIterator;

  namespace Kinetics {
    namespace DiffTransConfigIsEqualFastImpl {
      template<typename T> struct Common;
      struct Transformed;
      struct Untransformed;
    }

    class DiffTransConfiguration;

    /// \brief Class for less than comparison of Configurations
    ///
    /// - Possibly more efficient implementation because symops/permutations are
    ///   applied one-by-one as needed until the comparison can be made instead
    ///   of applying to all DoFs. Might reduce #operations, #mem alloc, etc.
    class DiffTransConfigIsEqualFast {

    public:

      DiffTransConfigIsEqualFast(const DiffTransConfiguration &_dtconfig);

      /// Custom tolerance is ignored currently. The controlling tol is primclex.crystallography_tol()
      DiffTransConfigIsEqualFast(const DiffTransConfiguration &_dtconfig, double _tol) :
        DiffTransConfigIsEqualFast(_dtconfig) {}

      const DiffTransConfiguration &config() const {
        return *m_dtconfig;
      }

      /// \brief Check if dtconfig < other (may have different Supercell)
      bool operator()(const DiffTransConfiguration &other) const;

      /// \brief Check if dtconfig < A*dtconfig
      bool operator()(const PermuteIterator &A) const;

      /// \brief Check if A*dtconfig < B*dtconfig
      bool operator()(const PermuteIterator &A, const PermuteIterator &B) const;

      /// \brief Check if dtconfig < A*other
      bool operator()(const PermuteIterator &A, const DiffTransConfiguration &other) const;

      /// \brief Check if A*dtconfig < B*other
      bool operator()(const PermuteIterator &A, const PermuteIterator &B, const DiffTransConfiguration &other) const;

      bool is_less() const;

    private:

      template<typename T> friend struct DiffTransConfigIsEqualFastImpl::Common;
      friend DiffTransConfigIsEqualFastImpl::Untransformed;
      friend DiffTransConfigIsEqualFastImpl::Transformed;

      /// For use by friends
      void set_is_less(bool value) const;


      const DiffTransConfiguration *m_dtconfig;

      ScelPeriodicDiffTransSymCompare m_sym_compare;

      mutable bool m_is_less;
    };

    /// \brief Class for less than comparison of Configurations
    ///
    /// - Easy to understand implementation
    class DiffTransConfigIsEqualSimple {

    public:

      DiffTransConfigIsEqualSimple(const DiffTransConfiguration &_dtconfig);

      /// Custom tolerance is ignored currently. The controlling tol is primclex.crystallography_tol()
      DiffTransConfigIsEqualSimple(const DiffTransConfiguration &_dtconfig, double _tol) :
        DiffTransConfigIsEqualSimple(_dtconfig) {}

      const DiffTransConfiguration &config() const {
        return *m_dtconfig;
      }

      /// \brief Check if dtconfig < other (may have different Supercell)
      bool operator()(const DiffTransConfiguration &other) const;

      /// \brief Check if dtconfig < A*dtconfig
      bool operator()(const PermuteIterator &A) const;

      /// \brief Check if A*dtconfig < B*dtconfig
      bool operator()(const PermuteIterator &A, const PermuteIterator &B) const;

      /// \brief Check if dtconfig < A*other
      bool operator()(const PermuteIterator &A, const DiffTransConfiguration &other) const;

      /// \brief Check if A*dtconfig < B*other
      bool operator()(const PermuteIterator &A, const PermuteIterator &B, const DiffTransConfiguration &other) const;

      bool is_less() const;

    private:

      const DiffTransConfiguration *m_dtconfig;

      mutable bool m_is_less;
    };

  }
}

#endif
