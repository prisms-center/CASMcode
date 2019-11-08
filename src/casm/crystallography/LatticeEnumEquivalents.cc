#include "casm/crystallography/LatticeEnumEquivalents.hh"
#include "casm/crystallography/Lattice_impl.hh"
#include "casm/symmetry/SymTools.hh"
#include "casm/symmetry/SymTools_impl.hh"

namespace CASM {
  namespace xtal {

    namespace {

      struct MakeInvariantSubgroup {

        MakeInvariantSubgroup() {}

        template<typename SymOpIterator, typename SymOpOutputIterator>
        SymOpOutputIterator operator()(
          const Lattice &lat,
          SymOpIterator begin,
          SymOpIterator end,
          SymOpOutputIterator result) {
          return sym::invariant_subgroup(std::vector<CASM::SymOp>(begin, end), lat, result);
        }

      };

    }

    const std::string LatticeEnumEquivalents::enumerator_name = "LatticeEnumEquivalents";

    /// \brief Constructor
    ///
    /// \param lat Lattice to generate equivalents of
    /// \param super_g The super group to use to generate equivalents. Must have a valid MasterSymGroup.
    /// \param tol Tolerance used for checking equivalence
    ///
    /// \throws std::runtime_error if super_g does not have a MasterSymGroup
    ///
    LatticeEnumEquivalents::LatticeEnumEquivalents(const Lattice &lat, const SymGroup &super_g) :
      EnumEquivalents<Lattice, std::vector<CASM::SymOp>::const_iterator, CASM::SymOp, SymRepIndexCompare>(
        canonical::equivalent(lat, super_g), super_g.begin(), super_g.end(), MakeInvariantSubgroup()) {

      if(!super_g.begin()->has_valid_master()) {
        throw std::runtime_error("Error constructing LatticeEnumEquivalents: SymGroup has no MasterSymGroup");
      }

    }

  }
}
