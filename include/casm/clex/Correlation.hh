#ifndef CORRELATION_HH
#define CORRELATION_HH

#include "casm/container/Array.hh"
#include "casm/clusterography/Orbitree.hh"

namespace CASM {
  typedef Array<double> Correlation;
  template<class U>
  void match_shape(Correlation &corr, const Array<Array<Array<U> > > &thing) {

    //resize(thing.size());
    Index num_corr = 0;
    for(Index i = 0; i < thing.size(); i++) {
      //at(i).resize(thing[i].size());
      for(Index j = 0; j < thing[i].size(); j++) {
        num_corr += thing[i][j].size();
      }
    }
    corr.resize(num_corr, 0);
  }

  template<class ClustType>
  void match_shape(Correlation &corr, const GenericOrbitree<ClustType> &tree) {
    Index num_corr = 0;
    for(Index i = 0; i < tree.size(); i++) {
      for(Index j = 0; j < tree[i].size(); j++) {
        num_corr += tree.prototype(i, j).clust_basis.size();
      }
    }
    corr.resize(num_corr, 0);
  }
}

#endif
