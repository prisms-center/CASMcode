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
    ///
    /// - Used more generally as an Element in OccPerturbation
    class OccupationTransformation :
      public Comparisons<Translatable<DoFTransformation<CRTPBase<OccupationTransformation>>>> {

    public:

      OccupationTransformation(const UnitCellCoord &_uccoord,
                               Index _from_value,
                               Index _to_value);

      Index from_value;
      Index to_value;
      UnitCellCoord uccoord;

      const Structure &prim() const;

      const Molecule &from_mol() const;

      const Molecule &to_mol() const;

      bool operator<(const OccupationTransformation &B) const;

      OccupationTransformation &operator+=(UnitCell Frac);

      OccupationTransformation &apply_sym(const SymOp &op);

      Configuration &apply_to(Configuration &config) const;

      void reverse();

    private:

      friend DoFTransformation<CRTPBase<OccupationTransformation>>;

      Configuration &apply_reverse_to_impl(Configuration &config) const;

      std::tuple<UnitCellCoord, Index, Index> _tuple() const;

    };

    /// \brief Print OccupationTransformation to stream, using default Printer<Kinetics::OccupationTransformation>
    std::ostream &operator<<(std::ostream &sout, const OccupationTransformation &trans);

  }

  std::map<AtomSpecie, Index> empty_specie_count(const Structure &prim);

  template<typename OccTransfIt>
  std::map<AtomSpecie, Index> from_specie_count(OccTransfIt begin, OccTransfIt end);

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

    void print(const Element &element, Log &out);
  };
}

#endif
