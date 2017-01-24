#ifndef CASM_DoFTransformation
#define CASM_DoFTransformation

#include "casm/CASM_global_definitions.hh"
#include "casm/crystallography/UnitCellCoord.hh"

namespace CASM {

  class Structure;
  class Configuration;
  class SymOp;

  namespace Kinetics {

    /// \brief Abstract base class for kinetics
    class DoFTransformation {

    public:

      typedef Structure PrimType;

      DoFTransformation(const PrimType &prim);

      virtual ~DoFTransformation() {};

      const PrimType &prim() const {
        return *m_prim;
      }

      Configuration &apply_to(Configuration &config) const;

      Configuration &apply_reverse_to(Configuration &config) const;

      DoFTransformation &apply_sym(const SymOp &op);

      void reverse();

      std::unique_ptr<DoFTransformation> clone() const;

    private:

      virtual Configuration &apply_to_impl(Configuration &config) const = 0;

      virtual Configuration &apply_reverse_to_impl(Configuration &config) const = 0;

      virtual void apply_sym_impl(const SymOp &op) = 0;

      virtual void reverse_impl() = 0;

      virtual DoFTransformation *_clone() const = 0;


      const PrimType *m_prim;
    };


    class SiteDoFTransformation : public DoFTransformation {

    public:

      SiteDoFTransformation(const UnitCellCoord &_uccoord);

      virtual ~SiteDoFTransformation();

      SiteDoFTransformation &operator+=(UnitCell frac);

      SiteDoFTransformation &operator-=(UnitCell frac);

      SiteDoFTransformation &apply_sym(const SymOp &op);

      std::unique_ptr<SiteDoFTransformation> clone() const;

      UnitCellCoord uccoord;

    private:

      // inherited virtual functions:
      //virtual Configuration& apply_to_impl(Configuration& config) const = 0;
      //virtual Configuration& apply_reverse_to_impl(Configuration& config) const = 0;
      virtual SiteDoFTransformation *_clone() const override = 0;
      //virtual void reverse_impl() = 0;

      virtual void apply_sym_impl(const SymOp &op) override;

    };
  }
}

#endif
