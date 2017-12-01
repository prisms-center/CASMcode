#ifndef CASM_DiffTransOrbitDatabase
#define CASM_DiffTransOrbitDatabase

#include "casm/database/Database.hh"
#include "casm/kinetics/DiffusionTransformationTraits.hh"
#include "casm/symmetry/Orbit.hh"
#include "casm/clusterography/ClusterSymCompare.hh"

namespace CASM {

  namespace DB {

    /// Derived DiffTransOrbitDatabase must implement public methods:
    /// - std::pair<iterator,bool> insert(const PrimPeriodicDiffTransOrbit &orbit)
    ///
    ///
    /// Orbits are always generated in some sort of canonical form?
    ///
    template<>
    class Database<PrimPeriodicDiffTransOrbit> : public ValDatabase<PrimPeriodicDiffTransOrbit> {

    public:

      Database(const PrimClex &_primclex) :
        ValDatabase<PrimPeriodicDiffTransOrbit>(_primclex) {}

      virtual ~Database() {}

      /// Find PrimPeriodicDiffTransOrbit in database by comparing prototype
      virtual iterator search(const PrimPeriodicDiffTransOrbit &orbit) const;

      /// Find DiffusionTransformation in database by comparing to orbit prototypes
      virtual iterator search(const Kinetics::DiffusionTransformation &diff_trans) const;

      /// Range of PrimPeriodicDiffTransOrbit that were created from a given IntegralCluster
      ///
      /// - Should return range {end(), end()} if no PrimPeriodicDiffTransOrbit on specified IntegralCluster
      /// - Note: boost::iterator_range<iterator>::size is not valid for
      ///   DatabaseIterator.  Use boost::distance instead.
      /// virtual boost::iterator_range<iterator> cluster_range(const IntegralCluster &cluster) const = 0;

      /// Number of PrimPeriodicDiffTransOrbit in a particular cluster size
      /// Index cluster_range_size(const IntegralCluster &cluster) const;

      /// Number of PrimPeriodicDiffTransOrbit in a particular cluster size
      /// Index orbit_branch_size(const Index) const;

      /// Range of PrimPeriodicDiffTransOrbit in a particular cluster size
      /// virtual boost::iterator_range<iterator> orbit_branch_range(const Index) const;
    };

  }
}

#endif
