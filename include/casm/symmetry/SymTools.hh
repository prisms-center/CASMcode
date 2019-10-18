#ifndef SYMTOOLS_HH
#define SYMTOOLS_HH

#include<vector>

namespace CASM {
  class SymGroup;
  class SymOp;
  namespace xtal {
    class Lattice;
  }
  namespace sym {
    ///Returns the subgroup of the given group that keeps the lattice invariant
    SymGroup invariant_subgroup(const SymGroup &super_group, const xtal::Lattice &lat);

    std::vector<SymOp> invariant_subgroup(std::vector<SymOp>::const_iterator begin, std::vector<SymOp>::const_iterator end, const xtal::Lattice &lat);

    template<typename OutputIt>
    OutputIt invariant_subgroup(std::vector<SymOp>::const_iterator begin, std::vector<SymOp>::const_iterator end, const xtal::Lattice &lat, OutputIt result);
  }
}

#endif
