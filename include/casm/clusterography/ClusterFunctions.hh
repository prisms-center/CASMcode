#ifndef CLUSTERFUNCTIONS_HH
#define CLUSTERFUNCTIONS_HH

#include <string>

namespace CASM {

  template<class T> class Array;
  class Structure;
  template<class ClustType> class GenericOrbit;
  template<class ClustType> class GenericOrbitBranch;
  template<class ClustType> class GenericOrbitree;
  class SiteCluster;
  class HopCluster;
  typedef GenericOrbitree<HopCluster> HopOrbitree;
  class jsonParser;


  jsonParser &to_json(const GenericOrbitBranch<SiteCluster> &branch, jsonParser &json);
  /// Assumes the pivot lattice is already set
  void from_json(GenericOrbitBranch<SiteCluster> &branch, const jsonParser &json);

  jsonParser &to_json(const GenericOrbit<SiteCluster> &orbit, jsonParser &json);
  /// Assumes the prototype lattice is already set
  void from_json(GenericOrbit<SiteCluster> &orbit, const jsonParser &json);

};

#endif
