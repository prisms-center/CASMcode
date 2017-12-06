#ifndef CASM_DiffusionTransformation_impl
#define CASM_DiffusionTransformation_impl

#include "casm/kinetics/DiffusionTransformation.hh"
#include "casm/clex/HasCanonicalForm_impl.hh"
#include "casm/symmetry/Orbit_impl.hh"
#include "casm/clusterography/ClusterSymCompare_impl.hh"
#include "casm/symmetry/OrbitGeneration_impl.hh"

namespace CASM {

  namespace {
    struct IsIncluded {
      IsIncluded(const std::map<AtomSpecie, Index> &_specie_count) :
        specie_count(_specie_count) {}
      const std::map<AtomSpecie, Index> &specie_count;

      bool operator()(std::string name) {
        auto lambda = [&](const std::pair<AtomSpecie, Index> &val) {
          return (val.first.name() == name) && val.second != 0;
        };
        return std::any_of(specie_count.begin(), specie_count.end(), lambda);
      }
    };
  }

  template<typename NameIterator>
  bool includes_all(const std::map<AtomSpecie, Index> specie_count, NameIterator begin, NameIterator end) {
    return std::all_of(begin, end, IsIncluded(specie_count));
  }

  template<typename NameIterator>
  bool excludes_all(const std::map<AtomSpecie, Index> specie_count, NameIterator begin, NameIterator end) {
    return !std::any_of(begin, end, IsIncluded(specie_count));
  }
}

#endif
