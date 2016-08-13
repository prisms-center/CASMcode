#ifndef CASM_CORRELATION_HH
#define CASM_CORRELATION_HH

#include "casm/external/Eigen/Dense"
#include "casm/clusterography/Orbitree.hh"

namespace CASM {

  typedef Eigen::VectorXd Correlation;

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
    corr = Eigen::VectorXd::Zero(num_corr);
  }

  template<class ClustType>
  void match_shape(Correlation &corr, const GenericOrbitree<ClustType> &tree) {
    Index num_corr = 0;
    for(Index i = 0; i < tree.size(); i++) {
      for(Index j = 0; j < tree[i].size(); j++) {
        num_corr += tree.prototype(i, j).clust_basis.size();
      }
    }
    corr = Eigen::VectorXd::Zero(num_corr);
  }
}

#endif
