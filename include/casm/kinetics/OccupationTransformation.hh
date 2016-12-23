#ifndef CASM_OccupationTransformation
#define CASM_OccupationTransformation

#include "casm/kinetics/DoFTransformation.hh"

namespace CASM {

  class Configuration;
  class SymOp;

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

      bool operator<(const OccupationTransformation &B) const {
        return _tuple() < B._tuple();
      }

    private:

      Configuration &apply_to_impl(Configuration &config) const override;

      Configuration &apply_reverse_to_impl(Configuration &config) const override;

      void apply_sym_impl(const SymOp &op) override;

      void reverse_impl() override;

      OccupationTransformation *_clone() const override;

      std::tuple<UnitCellCoord, Index, Index> _tuple() const {
        return std::make_tuple(uccoord, from_value, to_value);
      }

    };
  }
}

#endif
