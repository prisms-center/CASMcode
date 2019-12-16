#ifndef SYMTOOLS_HH
#define SYMTOOLS_HH

#include <vector>

namespace CASM {
  class SymGroup;
  class SymOp;
  namespace xtal {
    class Lattice;
    class UnitCellCoord;
    class Structure;
  } // namespace xtal

  namespace Kinetics {
    class DiffusionTransformation;
  }

  namespace sym {
    /// Returns the subgroup of the given group that keeps the lattice invariant
    SymGroup invariant_subgroup(const SymGroup &super_group, const xtal::Lattice &lat);

    template <typename OutputIt>
    OutputIt invariant_subgroup(const std::vector<SymOp> &super_group, const xtal::Lattice &lat, OutputIt result);

    // TODO: Do we keep passing by reference or do we want to change our ways here
    // and start passing by pointer?
    // it could just be:
    // sym::apply(op, my_ucc, args)  // returns new value
    // sym::apply(op, &my_ucc, args) // modifies value
    /// Apply a transformation, in place, return reference to the provided object.
    template <typename Transform, typename Object, typename... Args>
    Object &apply(const Transform &transformation, Object &obj, const Args &... args);

    /// Copy and apply a transformation, retun a new transformed copy.
    template <typename Transform, typename Object, typename... Args>
    Object copy_apply(const Transform &transformation, Object obj_copy, const Args &... args) {
      apply(transformation, obj_copy, args...);
      return obj_copy;
    }
  } // namespace sym
} // namespace CASM

#endif
