#ifndef CASM_DiffusionTransformation
#define CASM_DiffusionTransformation

#include "casm/kinetics/DoFTransformation.hh"
#include "casm/kinetics/OccupationTransformation.hh"
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

      bool operator<(const SpecieLocation &B) const;

      const Molecule &mol() const;

      const Specie &specie() const;

    private:

      std::tuple<UnitCellCoord, Index, Index> _tuple() const;
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

      bool specie_types_map() const;

      bool is_no_change() const;

      SpecieLocation from;
      SpecieLocation to;

      bool operator<(const SpecieTrajectory &B) const;

      void apply_sym(const SymOp &op);

      void reverse();

    private:

      std::tuple<SpecieLocation, SpecieLocation> _tuple() const;

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
    class LocalDiffTransSymCompare;
    class ScelPeriodicDiffTransSymCompare;
  }

  /// \brief Traits class for any ClusterSymCompare derived class
  ///
  template<>
  struct traits<Kinetics::PrimPeriodicDiffTransSymCompare> {

    typedef typename Kinetics::PrimPeriodicDiffTransSymCompare MostDerived;
    typedef typename Kinetics::DiffusionTransformation Element;
    typedef typename Kinetics::DiffusionTransformationInvariants InvariantsType;

  };


  /// \brief Traits class for any ClusterSymCompare derived class
  ///
  template<>
  struct traits<Kinetics::LocalDiffTransSymCompare> {

    typedef typename Kinetics::LocalDiffTransSymCompare MostDerived;
    typedef typename Kinetics::DiffusionTransformation Element;
    typedef typename Kinetics::DiffusionTransformationInvariants InvariantsType;
  };

  /// \brief Traits class for any ClusterSymCompare derived class
  ///
  template<>
  struct traits<Kinetics::ScelPeriodicDiffTransSymCompare> {

    typedef typename Kinetics::ScelPeriodicDiffTransSymCompare MostDerived;
    typedef typename Kinetics::DiffusionTransformation Element;
    typedef typename Kinetics::DiffusionTransformationInvariants InvariantsType;
  };


  namespace Kinetics {

    /// \brief Used to sort orbits
    class PrimPeriodicDiffTransSymCompare : public SymCompare<PrimPeriodicDiffTransSymCompare> {

    public:

      typedef traits<PrimPeriodicDiffTransSymCompare>::MostDerived MostDerived;
      typedef traits<PrimPeriodicDiffTransSymCompare>::Element Element;
      typedef traits<PrimPeriodicDiffTransSymCompare>::InvariantsType InvariantsType;

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

    /// \brief Used to sort orbits
    class LocalDiffTransSymCompare : public SymCompare<LocalDiffTransSymCompare> {

    public:

      typedef CASM_TMP::traits<LocalDiffTransSymCompare>::MostDerived MostDerived;
      typedef CASM_TMP::traits<LocalDiffTransSymCompare>::Element Element;
      typedef CASM_TMP::traits<LocalDiffTransSymCompare>::InvariantsType InvariantsType;

      LocalDiffTransSymCompare(double tol);

      double tol() const {
        return m_tol;
      }

    private:

      friend class SymCompare<LocalDiffTransSymCompare>;

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

    /// \brief Used to canonicalize DiffusionTransformations
    class ScelPeriodicDiffTransSymCompare : public SymCompare<ScelPeriodicDiffTransSymCompare> {

    public:

      typedef CASM_TMP::traits<ScelPeriodicDiffTransSymCompare>::MostDerived MostDerived;
      typedef CASM_TMP::traits<ScelPeriodicDiffTransSymCompare>::Element Element;
      typedef CASM_TMP::traits<ScelPeriodicDiffTransSymCompare>::InvariantsType InvariantsType;

      ScelPeriodicDiffTransSymCompare(const PrimGrid &prim_grid, double tol);

      double tol() const {
        return m_tol;
      }

    private:

      friend class SymCompare<ScelPeriodicDiffTransSymCompare>;

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

      const PrimGrid &m_prim_grid;

    };



  }

  /// \brief Traits class for DiffusionTransformation
  ///
  template<>
  struct traits<Kinetics::DiffusionTransformation> {

    typedef typename Kinetics::DiffusionTransformation MostDerived;
    typedef typename Kinetics::DiffusionTransformation Element;
    typedef typename Kinetics::DiffusionTransformationInvariants InvariantsType;

  };

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
      Index size() const;

      /// \brief Return the min pair distance, or 0.0 if size() <= 1
      double min_length() const;

      /// \brief Return the max pair distance, or 0.0 if size() <= 1
      double max_length() const;

      void apply_sym(const PermuteIterator &it) {
        apply_sym_impl(it.sym_op());
      }

      void apply_sym(const SymOp &op) {
        apply_sym_impl(op);
      }

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

    /// \brief Return a standardized name for this diffusion transformation orbit
    std::string orbit_name(const PrimPeriodicDiffTransOrbit &orbit);

    // \brief Returns the distance from uccoord to the closest point on a linearly
    /// interpolated diffusion path. (Could be an end point)
    double dist_to_path(const DiffusionTransformation &diff_trans, const UnitCellCoord &uccoord);

    /// \brief Determines which site is closest to the diffusion transformation and the distance
    std::pair<UnitCellCoord, double> _path_nearest_neighbor(const DiffusionTransformation &diff_trans) ;

    /// \brief Determines which site is closest to the diffusion transformation
    UnitCellCoord path_nearest_neighbor(const DiffusionTransformation &diff_trans);

    /// \brief Determines the nearest site distance to the diffusion path
    double min_dist_to_path(const DiffusionTransformation &diff_trans);

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
