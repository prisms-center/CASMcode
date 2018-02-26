#ifndef CASM_OccupationTransformation
#define CASM_OccupationTransformation

#include "casm/kinetics/DoFTransformation.hh"
#include "casm/kinetics/OccupationTransformationTraits.hh"
#include "casm/misc/Comparisons.hh"
#include "casm/crystallography/UnitCellCoord.hh"
#include "casm/app/AppIO.hh"

namespace CASM {

  class Configuration;
  class SymOp;
  class AtomSpecie;
  class Molecule;

  namespace Kinetics {

    /// \brief Describes how occupation values transform
    /// An OccupationTransformation object is a representation of a occupant on a given
    /// site. This is both a sub component of DiffusionTransformations (shuffling occupants)
    /// and a sub component of a perturbation (changing one configuration to another)
    class OccupationTransformation :
      public Comparisons<Translatable<DoFTransformation<CRTPBase<OccupationTransformation>>>> {

    public:
      /// Constructor Create an OccupationTransformation on a site indicating what the current
      /// occupant is on that site and the one it will change to
      OccupationTransformation(const UnitCellCoord &_uccoord,
                               Index _from_value,
                               Index _to_value);

      Index from_value;
      Index to_value;
      UnitCellCoord uccoord;

      /// The tiling unit of the infinite crystal that the site of this
      /// OccupationTransformation lives in
      const Structure &prim() const;

      /// The current occupant of the site
      const Molecule &from_mol() const;

      /// The future occupant of the site
      const Molecule &to_mol() const;

      /// Lexicographical comparison of OccupationTransformation for sorting purposes
      bool operator<(const OccupationTransformation &B) const;

      /// Rigidly shifts the coordinate of this transformation by a lattice shift
      OccupationTransformation &operator+=(UnitCell Frac);

      /// Applies symmetry to the coordinate of this transformation
      OccupationTransformation &apply_sym(const SymOp &op);

      /// Transform the occupation of config and return a configuration
      /// that is perturbed according to this transformation
      Configuration &apply_to(Configuration &config) const;

      /// Flip the initial and final state of this transformation
      void reverse();

    private:

      friend DoFTransformation<CRTPBase<OccupationTransformation>>;

      Configuration &apply_reverse_to_impl(Configuration &config) const;

      std::tuple<UnitCellCoord, Index, Index> _tuple() const;

    };

    /// \brief Print OccupationTransformation to stream, using default Printer<Kinetics::OccupationTransformation>
    std::ostream &operator<<(std::ostream &sout, const OccupationTransformation &trans);

  }
  /// Returns an empty map for a given infinite crystal
  std::map<AtomSpecie, Index> empty_specie_count(const Structure &prim);

  /// Iterates over a vector of Occupation Transformation to give the count map of initial species
  template<typename OccTransfIt>
  std::map<AtomSpecie, Index> from_specie_count(OccTransfIt begin, OccTransfIt end);

  /// Iterates over a vector of OccupationTransformation to give the count map of final species
  template<typename OccTransfIt>
  std::map<AtomSpecie, Index> to_specie_count(OccTransfIt begin, OccTransfIt end);


  /// \brief Write OccupationTransformation to JSON object
  jsonParser &to_json(const Kinetics::OccupationTransformation &trans, jsonParser &json);

  template<>
  struct jsonConstructor<Kinetics::OccupationTransformation> {

    static Kinetics::OccupationTransformation from_json(const jsonParser &json, const Structure &prim);
  };

  void from_json(Kinetics::OccupationTransformation &fill_value, const jsonParser &read_json);

  template<>
  struct Printer<Kinetics::OccupationTransformation> : public PrinterBase {

    typedef Kinetics::OccupationTransformation Element;
    static const std::string element_name;

    Printer(int _indent_space = 6, char _delim = '\n', COORD_TYPE _mode = INTEGRAL) :
      PrinterBase(_indent_space, _delim, _mode) {}

    void print(const Element &element, std::ostream &out);
  };
}

#endif
