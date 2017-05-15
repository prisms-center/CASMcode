#ifndef CASM_OccupationTransformation
#define CASM_OccupationTransformation

#include "casm/kinetics/DoFTransformation.hh"
#include "casm/app/AppIO.hh"

namespace CASM {

  class Configuration;
  class SymOp;
  class Molecule;
  class Specie;

  namespace Kinetics {

    /// \brief Describes how occupation values transform
    class OccupationTransformation :
      public SiteDoFTransformation,
      public Comparisons<OccupationTransformation> {

    public:

      OccupationTransformation(const UnitCellCoord &uccoord,
                               Index from_value,
                               Index to_value);

      OccupationTransformation &apply_sym(const SymOp &op);

      std::unique_ptr<OccupationTransformation> clone() const;

      Index from_value;
      Index to_value;

      const Molecule &from_mol() const;

      const Molecule &to_mol() const;

      bool operator<(const OccupationTransformation &B) const;

    private:

      Configuration &apply_to_impl(Configuration &config) const override;

      Configuration &apply_reverse_to_impl(Configuration &config) const override;

      void reverse_impl() override;

      OccupationTransformation *_clone() const override;

      std::tuple<UnitCellCoord, Index, Index> _tuple() const;

    };

    /// \brief Print OccupationTransformation to stream, using default Printer<Kinetics::OccupationTransformation>
    std::ostream &operator<<(std::ostream &sout, const OccupationTransformation &trans);

  }

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
