#include "casm/clusterography/ClusterInvariants.hh"

namespace casm {
  
  /// \brief Check if ClusterInvariants are equal
  bool almost_equal(const ClusterInvariants &A, const ClusterInvariants &B, double tol) {
    return A.size() == B.size() &&
           std::equal(A.displacement().cbegin(),
                      A.displacement().cend(),
                      B.displacement().cbegin(),
                      [&](double a, double b) {
                        return almost_equal(a,b,tol);
                      });
  }
  
  /// \brief Compare ClusterInvariants
  ///
  /// \returns True if A < B, to specified tolerance
  ///
  /// - First compares by number of sites in cluster
  /// - Then compare all displacements, from longest to shortest
  ///
  bool compare(const ClusterInvariants &A, const ClusterInvariants &B, double tol) {
    
    // first sort by number of sites in cluster
    if(A.size() < B.size()) {
      return true;
    }
    if(A.size() > B.size()) {
      return false;
    }
    
    // all displacements
    for(int i=A.displacement().size()-1; i>=0; i--) {
      if(almost_equal(A.displacement()[i], B.displacement()[i], tol)) {
        continue;
      }
      if( A.displacement()[i] < B.displacement()[i]) {
        return true;
      }
      if( A.displacement()[i] > B.displacement()[i]) {
        return false;
      }
    }
    return false;
    
  }
  
  /// \brief Print ClusterInvariants
  std::ostream& operator<<(std::ostream &sout, const ClusterInvariants &invariants) {
    
    if(invariants.size() <= 1) {
      sout << "  #Points: " << invariants.size();
    }
    else {
      sout << "  #Points: " << invariants.size();
      sout << "  Site Distances: {";
      for(int i=0; i<invariants.displacement().size(); i++) {
        if(i != 0) {
          sout << ", ";
        }
        sout << std::setprecision(5) << invariants.displacement()[i];
      }
      sout << "}";
    }
    return sout;
  }
  
}
