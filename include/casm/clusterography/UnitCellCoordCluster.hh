#ifndef CASM_UnitCellCoordCluster
#define CASM_UnitCellCoordCluster

#include <vector>
#include <iterator>

#include <boost/function_output_iterator.hpp>

#include "casm/misc/General.hh"
#include "casm/crystallography/PrimMotif.hh"
#include "casm/clusterography/Orbit.hh"
#include "casm/symmetry/UnitCellCoordSymOp.hh"

namespace casm {
  
  /** \defgroup Clusterography
      
      \brief Functions and classes related to clusters
  */
  
  /* -- UnitCellCoordCluster Declaration ------------------------------------- */
  
  typedef GenericCluster<UnitCellCoord> UnitCellCoordCluster;
  
  
  /* -- BasicUCCCSymCompare Declaration ------------------------------------- */
  
  class BasicUCCCSymCompare : public SymCompare<UnitCellCoordCluster> {
  
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
    bool intra_orbit_compare(const UnitCellCoordCluster& A, const UnitCellCoordCluster& B) const override {
      return cluster_intra_orbit_compare(A,B);
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
    bool inter_orbit_compare(const UnitCellCoordCluster& A, const UnitCellCoordCluster& B) const override {
      cluster_inter_orbit_compare(A, B, std::less<UnitCellCoord>(), tol());
    }
    
    /// \brief Return tolerance
    double tol() const {
      return m_tol;
    }
    
    private:
    
    /// \brief Private virtual apply_sym
    virtual void _apply_sym(const SymOp& op) const override {
      return;
    }
    
    double m_tol;
    
  };
  
  
  /* -- LocalSymCompare<UnitCellCoordCluster> Declaration ------------------------------------- */
  
  template<>
  class LocalSymCompare<UnitCellCoordCluster> : public BasicUCCCSymCompare {
  
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
    UnitCellCoordCluster prepare(UnitCellCoordCluster obj) const override {
      std::sort(obj.begin(), obj.end());
      return obj;
    }
    
    /// \brief Apply symmetry to this
    ///
    /// - Affects no change
    LocalSymCompare& apply_sym(const SymOp& op) {
      this->_apply_sym(op);
      return *this;
    }
    
    /// \brief Public non-virtual clone
    std::unique_ptr<LocalSymCompare> clone() const {
      return std::unique_ptr<LocalSymCompare>(this->_clone());
    }
    
    private:
    
    /// \brief Private virtual clone
    virtual LocalSymCompare* _clone() const {
      return new LocalSymCompare(*this);
    }
    
    double m_tol;
    
  };
  
  
  /* -- PrimPeriodicSymCompare<UnitCellCoordCluster> Declaration ------------------------------------- */
  
  template<>
  class PrimPeriodicSymCompare<UnitCellCoordCluster> : public BasicUCCCSymCompare {
  
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
    UnitCellCoordCluster prepare(UnitCellCoordCluster obj) const override {
      std::sort(obj.begin(), obj.end());
      return obj - obj[0].unitcell();
    }
    
    /// \brief Apply symmetry to this
    ///
    /// - Affects no change
    PrimPeriodicSymCompare& apply_sym(const SymOp& op) {
      this->_apply_sym(op);
      return *this;
    }
    
    /// \brief Public non-virtual clone
    std::unique_ptr<PrimPeriodicSymCompare> clone() const {
      return std::unique_ptr<PrimPeriodicSymCompare>(this->_clone());
    }
    
    
    private:
    
    /// \brief Private virtual clone
    virtual PrimPeriodicSymCompare* _clone() const {
      return new PrimPeriodicSymCompare(*this);
    }
    
  };
  
  
  /* -- ScelPeriodicSymCompare<UnitCellCoordCluster> Declaration ------------------------------------- */
  
  template<>
  class ScelPeriodicSymCompare<UnitCellCoordCluster> : public BasicUCCCSymCompare {
  
    public:
    
    /// \brief Constructor
    ///
    /// \param prim_grid A prim_grid reference
    /// \param tol Tolerance for inter_orbit_compare of site-to-site distances
    ///
    ScelPeriodicSymCompare(const PrimGrid& prim_grid, double tol):
      BasicUCCCSymCompare(tol),
      m_prim_grid(prim_grid) {}
    
    /// \brief Prepare an element for comparison
    ///
    /// - Sorts UnitCellCoord and translates so that obj[0] is within the supercell
    UnitCellCoordCluster prepare(UnitCellCoordCluster obj) const override {
      std::sort(obj.begin(), obj.end());
      auto trans = obj[0].unitcell() - m_prim_grid.get_within(obj[0]).unitcell();
      return obj - trans;
    }
    
    /// \brief Apply symmetry to this
    ///
    /// - Affects no change
    PrimPeriodicSymCompare& apply_sym(const SymOp& op) {
      this->_apply_sym(op);
      return *this;
    }
    
    /// \brief Public non-virtual clone
    std::unique_ptr<PrimPeriodicSymCompare> clone() const {
      return std::unique_ptr<PrimPeriodicSymCompare>(this->_clone());
    }
    
    
    private:
    
    /// \brief Private virtual clone
    virtual ScelPeriodicSymCompare* _clone() const {
      return new ScelPeriodicSymCompare(*this);
    }
    
    const PrimGrid& m_prim_grid;
    
  };
  
}

#endif