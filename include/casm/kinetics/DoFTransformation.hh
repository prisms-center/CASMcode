#ifndef CASM_DoFTransformation
#define CASM_DoFTransformation

#include "casm/CASM_global_definitions.hh"
#include "casm/crystallography/UnitCellCoord.hh"

namespace CASM {

  class Structure;
  class Configuration;
  class SymOp;

  namespace Kinetics {

    class DoFTransformation {

    public:

      typedef Structure PrimType;

      DoFTransformation(const PrimType &prim);

      virtual ~DoFTransformation() {};

      Configuration &apply_to(Configuration &config) const;

      Configuration &apply_reverse_to(Configuration &config) const;

      DoFTransformation &apply_sym(const SymOp &op);

      void reverse();

      std::unique_ptr<DoFTransformation> clone() const;

    private:

      virtual Configuration &apply_to_impl(Configuration &config) const = 0;

      virtual Configuration &apply_reverse_to_impl(Configuration &config) const = 0;

      virtual void apply_sym_impl(const SymOp &op) = 0;

      virtual reverse_impl() = 0;

      virtual DoFTransformation *_clone() const = 0;


      const PrimType *m_prim;
    };


    class SiteDoFTransformation : public DoFTransformation {

    public:

      SiteDoFTransformation(const UnitCellCoord &uccoord);

      virtual ~SiteDoFTransformation();

      UnitCellCoord &uccoord();

      const UnitCellCoord &uccoord() const;

      SiteDoFTransformation &operator+=(UnitCell frac);

      SiteDoFTransformation &operator-=(UnitCell frac);

      SiteDoFTransformation &apply_sym(const SymOp &op);

      std::unique_ptr<SiteDoFTransformation> clone() const;


    private:

      // inherited virtual functions:
      //virtual Configuration& apply_to_impl(Configuration& config) const = 0;
      //virtual Configuration& apply_reverse_to_impl(Configuration& config) const = 0;
      //virtual DoFTransformation* _clone() const = 0;
      //virtual reverse_impl() = 0;

      virtual void apply_sym_impl(const SymOp &op) override;

      UnitCellCoord m_uccoord;

    };


    /// \brief Describes how occupantation values transform
    class OccupationTransformation : public SiteDoFTransformation {

    public:

      OccupationTransformation(const UnitCellCoord &uccoord,
                               Index from_value,
                               Index to_value);

      Index &from_value();

      const Index &from_value() const;

      Index &to_value();

      const Index &to_value() const;

      OccupationTransformation &apply_sym(const SymOp &op);

      std::unique_ptr<OccupationTransformation> clone() const;


    private:

      Configuration &apply_to_impl(Configuration &config) const override;

      Configuration &apply_reverse_to_impl(Configuration &config) const override;

      void apply_sym_impl(const SymOp &op) override;

      void reverse_impl() override;

      OccupationTransformation *_clone() const override;

      Index m_from_value;
      Index m_to_value;

    };


    /// \brief Describes how atoms move
    class AtomTrajectory : public DoFTransformation {

    public:

      AtomTrajectory(const UnitCellCoord &from_uccoord,
                     Index from_value,
                     Index from_atom_index,
                     const UnitCellCoord &to_uccoord,
                     Index to_value,
                     Index to_atom_index);


      UnitCellCoord &from_uccoord();

      const UnitCellCoord &from_uccoord() const;

      UnitCellCoord &to_uccoord();

      const UnitCellCoord &to_uccoord() const;


      Index &from_value();

      const Index &from_value() const;

      Index &to_value();

      const Index &to_value() const;


      Index &from_atom_index();

      const Index &from_atom_index() const;

      Index &to_atom_index();

      const Index &to_atom_index() const;

      AtomTrajectory &operator+=(UnitCell frac);

      AtomTrajectory &operator-=(UnitCell frac);

      std::unique_ptr<AtomTrajectory> clone() const;

    private:

      Configuration &apply_to_impl(Configuration &config) const override;

      Configuration &apply_reverse_to_impl(Configuration &config) const override;

      void apply_sym_impl(const SymOp &op) override;

      void reverse_impl() override;

      AtomTrajectory *_clone() const override;

      UnitCellCoord m_from_uccoord;
      Index m_from_value;
      Index m_from_atom_index;

      UnitCellCoord m_to_uccoord;
      Index m_to_value;
      Index m_to_atom_index;

    };

    /// \brief Describes how atoms move
    class DiffusionTransformation : public DoFTransformation {

    public:

      DiffusionTransformation(const PrimType &prim);


      std::vector<OccupationTransformation> &occ_transform();

      const std::vector<OccupationTransformation> &occ_transform() const;


      std::vector<AtomTrajectory> &atom_traj();

      const std::vector<AtomTrajectory> &atom_traj() const;


      AtomTrajectory &operator+=(UnitCell frac);

      AtomTrajectory &operator-=(UnitCell frac);

      std::unique_ptr<DiffusionTransformation> clone() const;

    private:

      Configuration &apply_to_impl(Configuration &config) const override;

      Configuration &apply_reverse_to_impl(Configuration &config) const override;

      DoFTransformation &apply_sym_impl(const SymOp &op) override;

      void reverse_impl() override;

      DiffusionTransformation *_clone() const override;


      std::vector<OccupationTransformation> m_occ_transform;
      std::vector<AtomTrajectory> m_atom_traj;

    };
  }
}

#endif CASM_DoFTransformation
