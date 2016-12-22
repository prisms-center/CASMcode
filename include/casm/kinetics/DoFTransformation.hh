#ifndef CASM_DoFTransformation
#define CASM_DoFTransformation

#include "casm/CASM_global_definitions.hh"
#include "casm/misc/CASM_TMP.hh"
#include "casm/misc/cloneable_ptr.hh"
#include "casm/crystallography/UnitCellCoord.hh"
#include "casm/clusterography/ClusterInvariants.hh"
#include "casm/clusterography/IntegralCluster.hh"

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

    /// \brief Specifies a particular specie
    struct SpecieLocation : public Comparisons<SpecieLocation> {

      SpecieLocation(const UnitCellCoord &_uccoord, Index _occ, Index _pos);

      UnitCellCoord uccoord;

      /// Occupant index
      Index occ;

      /// Position of specie in Molecule
      Index pos;

      bool operator<(const SpecieLocation &B) const {
        return _tuple() < B._tuple();
      }

      const Molecule &mol() const {
        return uccoord.site().site_occupant()[occ];
      }

      const Specie &specie() const {
        return mol()[pos].specie;
      }

    private:

      std::tuple<UnitCellCoord, Index, Index> _tuple() const {
        return std::make_tuple(uccoord, occ, pos);
      }
    };

    /// \brief Describes how one specie moves
    class SpecieTrajectory : public Comparisons<SpecieTrajectory> {

    public:

      SpecieTrajectory(const SpecieLocation &_from, const SpecieLocation &_to);

      SpecieTrajectory &operator+=(UnitCell frac);

      SpecieTrajectory &operator-=(UnitCell frac);

      bool specie_types_map() const {
        return from.specie() == to.specie();
      }

      bool is_no_change() const {
        return from == to;
      }

      SpecieLocation from;
      SpecieLocation to;

      bool operator<(const SpecieTrajectory &B) const {
        return _tuple() < B._tuple();
      }

      void apply_sym(const SymOp &op);

      void reverse();

    private:

      std::tuple<SpecieLocation, SpecieLocation> _tuple() const {
        return std::make_tuple(from, to);
      }

    };

    class DiffusionTransformation;

    /// \brief Invariants of a DiffusionTransformation, used to sort orbits
    class DiffusionTransformationInvariants {

    public:

      DiffusionTransformationInvariants(const DiffusionTransformation &trans);

      ClusterInvariants<IntegralCluster> cluster_invariants;
      std::map<Specie, Index> specie_count;

    };
  }


  /// \brief Check if DiffusionTransformationInvariants are equal
  bool almost_equal(const Kinetics::DiffusionTransformationInvariants &A,
                    const Kinetics::DiffusionTransformationInvariants &B,
                    double tol);

  /// \brief Compare DiffusionTransformationInvariants
  bool compare(const Kinetics::DiffusionTransformationInvariants &A,
               const Kinetics::DiffusionTransformationInvariants &B,
               double tol);

  /// \brief Print DiffusionTransformationInvariants
  std::ostream &operator<<(std::ostream &sout,
                           const Kinetics::DiffusionTransformationInvariants &obj);

  namespace Kinetics {
    class PrimPeriodicDiffTransSymCompare;
  }

  namespace CASM_TMP {

    /// \brief Traits class for any ClusterSymCompare derived class
    ///
    template<>
    struct traits<Kinetics::PrimPeriodicDiffTransSymCompare> {

      typedef typename Kinetics::PrimPeriodicDiffTransSymCompare MostDerived;
      typedef typename Kinetics::DiffusionTransformation Element;
      typedef typename Kinetics::DiffusionTransformationInvariants InvariantsType;

    };
  }

  namespace Kinetics {

    /// \brief Used to sort orbits
    class PrimPeriodicDiffTransSymCompare : public SymCompare<PrimPeriodicDiffTransSymCompare> {

    public:

      typedef CASM_TMP::traits<PrimPeriodicDiffTransSymCompare>::MostDerived MostDerived;
      typedef CASM_TMP::traits<PrimPeriodicDiffTransSymCompare>::Element Element;
      typedef CASM_TMP::traits<PrimPeriodicDiffTransSymCompare>::InvariantsType InvariantsType;

      PrimPeriodicDiffTransSymCompare(double tol);

      double tol() const {
        return m_tol;
      }

    private:

      Element prepare_impl(const Element &A) const;

      bool compare_impl(const Element &A, const Element &B) const;

      bool invariants_compare_impl(const Element &A, const Element &B) const;

      double m_tol;

    };
  }

  namespace CASM_TMP {

    /// \brief Traits class for DiffusionTransformation
    ///
    template<>
    struct traits<Kinetics::DiffusionTransformation> {

      typedef typename Kinetics::DiffusionTransformation MostDerived;
      typedef typename Kinetics::DiffusionTransformation Element;
      typedef typename Kinetics::DiffusionTransformationInvariants InvariantsType;

    };
  }

  namespace Kinetics {

    /// \brief Describes how species move
    class DiffusionTransformation :
      public DoFTransformation,
      public Comparisons<DiffusionTransformation>,
      public SymComparable<DiffusionTransformation> {

    public:

      DiffusionTransformation(const PrimType &prim);


      DiffusionTransformation &operator+=(UnitCell frac);

      DiffusionTransformation &operator-=(UnitCell frac);

      std::unique_ptr<DiffusionTransformation> clone() const;

      bool is_valid_occ_transform() const;

      bool is_valid_specie_traj() const;

      bool is_valid() const;

      bool is_canonical() const;
      bool is_canonical(const SymGroup &g) const;

      SymOp to_canonical() const;
      SymOp to_canonical(const SymGroup &g) const;

      SymOp from_canonical() const;
      SymOp from_canonical(const SymGroup &g) const;

      DiffusionTransformation canonical_form() const;
      DiffusionTransformation canonical_form(const SymGroup &g) const;

      std::vector<OccupationTransformation> &occ_transform();
      const std::vector<OccupationTransformation> &occ_transform() const;

      std::vector<SpecieTrajectory> &specie_traj();
      const std::vector<SpecieTrajectory> &specie_traj() const;

      const IntegralCluster &cluster() const;
      const std::map<Specie, Index> &specie_count() const;

      bool operator<(const DiffusionTransformation &B) const;

      DiffusionTransformation &prepare();

      DiffusionTransformation prepared() const;

      DiffusionTransformation copy_apply_sym(const SymOp &op) const;

    private:

      Configuration &apply_to_impl(Configuration &config) const override;

      Configuration &apply_reverse_to_impl(Configuration &config) const override;

      void apply_sym_impl(const SymOp &op) override;

      void reverse_impl() override;

      DiffusionTransformation *_clone() const override;

      bool _specie_types_map() const;

      bool _breaks_indivisible_mol() const;

      bool _is_subcluster_transformation() const;

      void _forward_sort();

      bool _lt(const DiffusionTransformation &B) const;

      /// \brief Reset mutable members, cluster and invariants, when necessary
      void _reset();

      std::map<Specie, Index> _from_specie_count() const;
      std::map<Specie, Index> _to_specie_count() const;
      std::map<Specie, Index> _empty_specie_count() const;


      std::vector<OccupationTransformation> m_occ_transform;
      std::vector<SpecieTrajectory> m_specie_traj;

      // stores IntegralCluster, based on occ_transform uccoord
      mutable notstd::cloneable_ptr<IntegralCluster> m_cluster;

      // stores Specie -> count, using 'from' specie
      // - is equal to 'to' specie count if is_valid_occ_transform() == true
      mutable notstd::cloneable_ptr<std::map<Specie, Index> > m_specie_count;

    };
  }
}

#endif
