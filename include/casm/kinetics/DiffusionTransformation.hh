#ifndef CASM_DiffusionTransformation
#define CASM_DiffusionTransformation

#include "casm/kinetics/DoFTransformation.hh"
#include "casm/kinetics/OccupationTransformation.hh"
#include "casm/misc/CASM_TMP.hh"
#include "casm/misc/cloneable_ptr.hh"
#include "casm/crystallography/UnitCellCoord.hh"
#include "casm/clusterography/ClusterInvariants.hh"
#include "casm/clusterography/IntegralCluster.hh"
#include "casm/casm_io/jsonParser.hh"
#include "casm/app/AppIO.hh"

namespace CASM {

  class Structure;
  class Configuration;
  class SymOp;

  namespace Kinetics {

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
        return uccoord.sublat_site().site_occupant()[occ];
      }

      const Specie &specie() const {
        return mol()[pos].specie;
      }

    private:

      std::tuple<UnitCellCoord, Index, Index> _tuple() const {
        return std::make_tuple(uccoord, occ, pos);
      }
    };

    /// \brief Print DiffusionTransformationInvariants
    std::ostream &operator<<(std::ostream &sout, const SpecieLocation &obj);

  }

  jsonParser &to_json(const Kinetics::SpecieLocation &obj, jsonParser &json);

  template<>
  struct jsonConstructor<Kinetics::SpecieLocation> {

    static Kinetics::SpecieLocation from_json(const jsonParser &json, const Structure &prim);
  };

  void from_json(Kinetics::SpecieLocation &obj, const jsonParser &json);


  namespace Kinetics {

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
  }

  jsonParser &to_json(const Kinetics::SpecieTrajectory &traj, jsonParser &json);

  template<>
  struct jsonConstructor<Kinetics::SpecieTrajectory> {

    static Kinetics::SpecieTrajectory from_json(const jsonParser &json, const Structure &prim);
  };

  void from_json(Kinetics::SpecieTrajectory &traj, const jsonParser &json);


  namespace Kinetics {

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

      friend class SymCompare<PrimPeriodicDiffTransSymCompare>;

      Element prepare_impl(const Element &A) const;

      bool compare_impl(const Element &A, const Element &B) const;

      bool invariants_compare_impl(const InvariantsType &A, const InvariantsType &B) const;

      /// \brief Apply symmetry to this
      ///
      /// - Affects no change
      void apply_sym_impl(const SymOp &op) {
        return;
      }

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

      /// \brief Check specie_types_map() && !breaks_indivisible_mol() && !is_subcluster_transformation()
      bool is_valid_specie_traj() const;

      bool specie_types_map() const;

      bool breaks_indivisible_mol() const;

      bool is_subcluster_transformation() const;

      bool is_valid() const;

      std::vector<OccupationTransformation> &occ_transform();
      const std::vector<OccupationTransformation> &occ_transform() const;

      std::vector<SpecieTrajectory> &specie_traj();
      const std::vector<SpecieTrajectory> &specie_traj() const;

      const IntegralCluster &cluster() const;
      const std::map<Specie, Index> &specie_count() const;

      /// \brief Compare DiffusionTransformation
      ///
      /// - Comparison is made using the sorted forms
      bool operator<(const DiffusionTransformation &B) const;

      DiffusionTransformation &sort();

      DiffusionTransformation sorted() const;

      bool is_sorted() const;

      /// \brief Return the cluster size
      Index size() const {
        return cluster().size();
      }

      /// \brief Return the min pair distance, or 0.0 if size() <= 1
      double min_length() const {
        return cluster().min_length();
      }

      /// \brief Return the max pair distance, or 0.0 if size() <= 1
      double max_length() const {
        return cluster().max_length();
      }

      /// \brief Return a standardized name for this diffusion transformation
      std::string name() const;

    private:

      Configuration &apply_to_impl(Configuration &config) const override;

      Configuration &apply_reverse_to_impl(Configuration &config) const override;

      void apply_sym_impl(const SymOp &op) override;

      void reverse_impl() override;

      DiffusionTransformation *_clone() const override;

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

    /// \brief Print DiffusionTransformation to stream, using default Printer<Kinetics::DiffusionTransformation>
    std::ostream &operator<<(std::ostream &sout, const DiffusionTransformation &trans);

    typedef Orbit<DiffusionTransformation, PrimPeriodicDiffTransSymCompare> PrimPeriodicDiffTransOrbit;

  }

  /// \brief Write DiffusionTransformation to JSON object
  jsonParser &to_json(const Kinetics::DiffusionTransformation &trans, jsonParser &json);

  template<>
  struct jsonConstructor<Kinetics::DiffusionTransformation> {

    static Kinetics::DiffusionTransformation from_json(const jsonParser &json, const Structure &prim);
  };

  /// \brief Read from JSON
  void from_json(Kinetics::DiffusionTransformation &trans, const jsonParser &json, const Structure &prim);

  template<>
  struct Printer<Kinetics::DiffusionTransformation> : public PrinterBase {

    typedef Kinetics::DiffusionTransformation Element;
    static const std::string element_name;

    Printer(int _indent_space = 6, char _delim = '\n', COORD_TYPE _mode = INTEGRAL) :
      PrinterBase(_indent_space, _delim, _mode) {}

    void print(const Element &element, std::ostream &out);
  };

}

#endif
