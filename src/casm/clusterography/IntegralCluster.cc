#include "casm/clusterography/IntegralCluster.hh"
#include "casm/crystallography/Coordinate.hh"
#include "casm/crystallography/Structure.hh"

namespace CASM {

  IntegralCluster::IntegralCluster(PrimType const &prim):
    m_prim_ptr(&prim) {}

  typename IntegralCluster::PrimType const &IntegralCluster::prim() const {
    return *m_prim_ptr;
  }

  /// \brief Access vector of elements
  std::vector<xtal::UnitCellCoord> &IntegralCluster::elements() {
    return m_element;
  }

  /// \brief const Access vector of elements
  const std::vector<xtal::UnitCellCoord> &IntegralCluster::elements() const {
    return m_element;
  }

  /// \brief Return the coordinate corresponding to element(i)
  xtal::Coordinate IntegralCluster::coordinate(size_type i) const {
    return this->element(i).coordinate(prim());
  }

  /// \brief Translate the cluster by a UnitCell translation
  IntegralCluster &IntegralCluster::operator+=(xtal::UnitCell trans) {
    for(auto it = this->begin(); it != this->end(); ++it) {
      *it += trans;
    }
    return *this;
  }
}

#include "casm/symmetry/SymTools.hh"
namespace CASM {

  namespace sym {
    /// Apply SymOp to IntegralCluster
    ///
    /// Specialization of template from casm/symmetry/Symtools.hh:
    ///   template <typename Transform, typename Object, typename... Args>
    ///   Object &apply(const Transform &transformation, Object &obj, const Args &... args);
    template<>
    IntegralCluster &apply(SymOp const &op, IntegralCluster &clust, Structure const &prim) {
      for(auto &e : clust) {
        sym::apply(op, e, prim);
      }
      return clust;
    }
  }
}
