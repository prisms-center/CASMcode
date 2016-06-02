#ifndef CASM_IntegralCluster
#define CASM_IntegralCluster

#include <vector>

#include "casm/symmetry/Orbit.hh"
#include "casm/crystallography/PrimGrid.hh"
#include "casm/crystallography/UnitCellCoord.hh"
#include "casm/clusterography/ClusterInvariants.hh"
#include "casm/clusterography/CoordCluster.hh"

namespace CASM {

  /** \defgroup Clusterography

      \brief Functions and classes related to clusters
  */

  /* -- IntegralCluster Declaration ------------------------------------- */

  /// \brief Cluster of UnitCellCoord
  ///
  /// \ingroup Clusterography
  ///
  typedef CoordCluster<UnitCellCoord> IntegralCluster;


  /* -- BasicUCCCSymCompare Declaration ------------------------------------- */

  /// \brief Abstract base class for IntegralCluster comparisons
  class BasicUCCCSymCompare : public SymCompare<IntegralCluster> {

  public:

    /// \brief Constructor
    ///
    /// \param tol Tolerance for inter_orbit_compare of site-to-site distances
    ///
    BasicUCCCSymCompare(double tol):
      SymCompare(),
      m_tol(tol) {

    }

    /// \brief Orders 'prepared' elements in the same orbit
    ///
    /// - Returns 'true' to indicate A < B
    /// - Equivalence is indicated by \code !compare(A,B) && !compare(B,A) \endcode
    /// - Assumes elements are 'prepared' before being compared
    /// Implementation:
    /// - std::lexicographical_compare of UnitCellCoord in A and B
    bool intra_orbit_compare(const IntegralCluster &A, const IntegralCluster &B) const override {
      return cluster_intra_orbit_compare(A, B);
    }

    /// \brief Orders orbit prototypes in canonical form
    ///
    /// - Returns 'true' to indicate A < B
    /// - Equivalence is indicated by \code !compare(A,B) && !compare(B,A) \endcode
    /// - Assumes elements are in canonical form
    ///
    /// Implementation:
    /// - First compare A.size(), B.size()
    /// - Second compare ClusterInvariants(A), ClusterInvariants(B)
    /// - Finally, std::lexicographical_compare of UnitCellCoord in A and B
    bool inter_orbit_compare(const IntegralCluster &A, const IntegralCluster &B) const override {
      return cluster_inter_orbit_compare(A, B, tol());
    }

    /// \brief Return tolerance
    double tol() const {
      return m_tol;
    }

  protected:

    /// \brief Private virtual apply_sym
    virtual void _apply_sym(const SymOp &op) const override {
      return;
    }

  private:

    double m_tol;

  };


  /* -- LocalSymCompare<IntegralCluster> Declaration ------------------------------------- */

  /// \brief Comparisons of IntegralCluster with aperiodic symmetry
  ///
  /// \relates IntegralCluster
  template<>
  class LocalSymCompare<IntegralCluster> : public BasicUCCCSymCompare {

  public:

    /// \brief Constructor
    ///
    /// \param tol Tolerance for inter_orbit_compare of site-to-site distances
    ///
    LocalSymCompare(double tol):
      BasicUCCCSymCompare(tol) {}

    /// \brief Prepare an element for comparison
    ///
    /// - Sorts UnitCellCoord
    IntegralCluster prepare(IntegralCluster obj) const override {
      std::sort(obj.begin(), obj.end());
      return obj;
    }

    /// \brief Apply symmetry to this
    ///
    /// - Affects no change
    LocalSymCompare &apply_sym(const SymOp &op) {
      this->_apply_sym(op);
      return *this;
    }

    /// \brief Public non-virtual clone
    std::unique_ptr<LocalSymCompare> clone() const {
      return std::unique_ptr<LocalSymCompare>(this->_clone());
    }

  private:

    /// \brief Private virtual clone
    virtual LocalSymCompare *_clone() const override {
      return new LocalSymCompare(*this);
    }

    double m_tol;

  };


  /* -- PrimPeriodicSymCompare<IntegralCluster> Declaration ------------------------------------- */

  /// \brief Comparisons of IntegralCluster with periodic symmetry of the primitive lattice
  ///
  /// \relates IntegralCluster
  template<>
  class PrimPeriodicSymCompare<IntegralCluster> : public BasicUCCCSymCompare {

  public:

    /// \brief Constructor
    ///
    /// \param tol Tolerance for inter_orbit_compare of site-to-site distances
    ///
    PrimPeriodicSymCompare(double tol):
      BasicUCCCSymCompare(tol) {}

    /// \brief Prepare an element for comparison
    ///
    /// - Sorts UnitCellCoord and translates so that obj[0] is in the origin unit cell
    IntegralCluster prepare(IntegralCluster obj) const override {
      std::sort(obj.begin(), obj.end());
      return obj - obj[0].unitcell();
    }

    /// \brief Apply symmetry to this
    ///
    /// - Affects no change
    PrimPeriodicSymCompare &apply_sym(const SymOp &op) {
      this->_apply_sym(op);
      return *this;
    }

    /// \brief Public non-virtual clone
    std::unique_ptr<PrimPeriodicSymCompare> clone() const {
      return std::unique_ptr<PrimPeriodicSymCompare>(this->_clone());
    }


  private:

    /// \brief Private virtual clone
    virtual PrimPeriodicSymCompare *_clone() const override {
      return new PrimPeriodicSymCompare(*this);
    }

  };


  /* -- ScelPeriodicSymCompare<IntegralCluster> Declaration ------------------------------------- */

  /// \brief Comparisons of IntegralCluster with periodic symmetry of a supercell lattice
  ///
  /// \relates IntegralCluster
  template<>
  class ScelPeriodicSymCompare<IntegralCluster> : public BasicUCCCSymCompare {

  public:

    /// \brief Constructor
    ///
    /// \param prim_grid A prim_grid reference
    /// \param tol Tolerance for inter_orbit_compare of site-to-site distances
    ///
    ScelPeriodicSymCompare(const PrimGrid &prim_grid, double tol):
      BasicUCCCSymCompare(tol),
      m_prim_grid(prim_grid) {}

    /// \brief Prepare an element for comparison
    ///
    /// - Sorts UnitCellCoord and translates so that obj[0] is within the supercell
    IntegralCluster prepare(IntegralCluster obj) const override {
      std::sort(obj.begin(), obj.end());
      auto trans = obj[0].unitcell() - m_prim_grid.get_within(obj[0]).unitcell();
      return obj - trans;
    }

    /// \brief Apply symmetry to this
    ///
    /// - Affects no change
    ScelPeriodicSymCompare &apply_sym(const SymOp &op) {
      this->_apply_sym(op);
      return *this;
    }

    /// \brief Public non-virtual clone
    std::unique_ptr<ScelPeriodicSymCompare> clone() const {
      return std::unique_ptr<ScelPeriodicSymCompare>(this->_clone());
    }


  private:

    /// \brief Private virtual clone
    virtual ScelPeriodicSymCompare *_clone() const override {
      return new ScelPeriodicSymCompare(*this);
    }

    const PrimGrid &m_prim_grid;

  };

}

#endif