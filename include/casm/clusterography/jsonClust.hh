#ifndef JSONCLUST_HH
#define JSONCLUST_HH

#include "casm/crystallography/Structure.hh"

#include "casm/clusterography/Orbitree.hh"
#include "casm/clusterography/SiteCluster.hh"

namespace CASM {

  template<typename ValueType>
  class ClustJsonHelper {
    ValueType &m_value;
    const Structure &m_struc;
    double m_tol;
  public:
    ClustJsonHelper(ValueType &_value, const Structure &_struc, double _tol = TOL) :
      m_value(_value), m_struc(_struc), m_tol(TOL) {};

    ValueType &value() {
      return m_value;
    };
    const ValueType &value()const {
      return m_value;
    }
    const Structure &struc()const {
      return m_struc;
    }

    double tol() const {
      return m_tol;
    }

    operator ClustJsonHelper<const ValueType>()const {
      return ClustJsonHelper<const ValueType>(m_value, m_struc);
    };
  };

  typedef ClustJsonHelper<SiteCluster> SiteClusterJsonHelper;
  typedef ClustJsonHelper<SiteOrbit> SiteOrbitJsonHelper;
  typedef ClustJsonHelper<SiteOrbitBranch> SiteOrbitBranchJsonHelper;
  typedef ClustJsonHelper<SiteOrbitree> SiteOrbitreeJsonHelper;

  typedef ClustJsonHelper<const SiteCluster> ConstSiteClusterJsonHelper;
  typedef ClustJsonHelper<const SiteOrbit> ConstSiteOrbitJsonHelper;
  typedef ClustJsonHelper<const SiteOrbitBranch> ConstSiteOrbitBranchJsonHelper;
  typedef ClustJsonHelper<const SiteOrbitree> ConstSiteOrbitreeJsonHelper;

  void from_json(SiteClusterJsonHelper clust_helper, const jsonParser &json);
  jsonParser &to_json(const ConstSiteClusterJsonHelper &clust_helper, jsonParser &json);

  /// Lean printing of Orbits (only the prototype gets printed; full orbit can later be reconstructed using struc.factor_group();
  jsonParser &to_json(const ConstSiteOrbitJsonHelper &orbit_helper, jsonParser &json);

  void from_json(SiteOrbitBranchJsonHelper branch_helper, const jsonParser &json);
  jsonParser &to_json(const ConstSiteOrbitBranchJsonHelper &branch_helper, jsonParser &json);

  void from_json(SiteOrbitreeJsonHelper tree_helper, const jsonParser &json);
  jsonParser &to_json(const ConstSiteOrbitreeJsonHelper &tree_helper, jsonParser &json);

  template<typename ValueType>
  ClustJsonHelper<ValueType> jsonHelper(ValueType &_value, const Structure &_struc, double tol = TOL) {
    return ClustJsonHelper<ValueType>(_value, _struc, tol);
  }
}
#endif
