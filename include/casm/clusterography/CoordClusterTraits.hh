#ifndef COORDCLUSTERTRAITS_HH
#define COORDCLUSTERTRAITS_HH

#include <string>
#include "casm/clusterography/ClusterDecl.hh"

namespace CASM {
  namespace sym {
    class CopyApplyWithPrim_f;
    class CopyApplyDefault_f;
  }

  template <typename Base>
  class CopyApplyWithPrim;

  template<typename CoordType>
  struct traits<CoordCluster<CoordType>> {
    typedef CoordType Element;
    typedef ClusterInvariants<CoordCluster<CoordType>> InvariantsType;
    static CoordType position(const CoordCluster<CoordType> &clust);
    typedef unsigned int size_type;
    static const std::string name;

    template <typename Base>
    using CopyApplyType = CopyApplyWithPrim<Base>;
    typedef sym::CopyApplyWithPrim_f copy_apply_f_type;
  };

  template<typename CoordType>
  CoordType traits<CoordCluster<CoordType>>::position(const CoordCluster<CoordType> &clust) {
    return clust[0];
  }

  template<typename CoordType>
  const std::string traits<CoordCluster<CoordType>>::name = "CoordCluster";
}

#endif
