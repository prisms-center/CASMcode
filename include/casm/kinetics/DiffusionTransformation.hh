#ifndef CASM_DiffusionTransformation
#define CASM_DiffusionTransformation

#include "casm/kinetics/DoFTransformation.hh"
#include "casm/kinetics/OccupationTransformation.hh"
#include "casm/misc/cloneable_ptr.hh"
#include "casm/symmetry/PermuteIterator.hh"
#include "casm/crystallography/UnitCellCoord.hh"
#include "casm/clusterography/ClusterInvariants.hh"
#include "casm/clusterography/IntegralCluster.hh"

namespace CASM {

  class AtomSpecie;
  class Structure;
  class Configuration;
  class SymOp;
  class jsonParser;
  template<typename T> struct jsonConstructor;

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

      const AtomSpecie &specie() const;

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
      std::map<AtomSpecie, Index> specie_count;

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
    template<typename Derived>
    class DiffTransSymCompare : public SymCompare<DiffTransSymCompare<Derived> > {
    public:
      DiffTransSymCompare(double _tol) :
        m_tol(_tol) {}

      double tol() const {
        return m_tol;
      }

    private:
      double m_tol;
    };
  }

  /// \brief Traits class for any ClusterSymCompare derived class
  ///
  template<typename Derived>
  struct traits<Kinetics::DiffTransSymCompare<Derived> > {

    typedef Derived MostDerived;
    typedef Kinetics::DiffusionTransformation Element;
    typedef Kinetics::DiffusionTransformationInvariants InvariantsType;

  };

  /// \brief Used to sort orbits
  template<>
  class PrimPeriodicSymCompare<Kinetics::DiffusionTransformation> :
    public Kinetics::DiffTransSymCompare<PrimPeriodicSymCompare<Kinetics::DiffusionTransformation> > {

  public:

    typedef Kinetics::DiffusionTransformation Element;
    typedef Kinetics::DiffusionTransformationInvariants InvariantsType;

    PrimPeriodicSymCompare(double tol);

  private:

    friend class SymCompare<Kinetics::DiffTransSymCompare<PrimPeriodicSymCompare<Kinetics::DiffusionTransformation> > >;

    Element prepare_impl(const Element &A) const;

    bool compare_impl(const Element &A, const Element &B) const;

    bool invariants_compare_impl(const InvariantsType &A, const InvariantsType &B) const;

    /// \brief Apply symmetry to this
    ///
    /// - Affects no change
    void apply_sym_impl(const SymOp &op) {
      return;
    }

  };


  /// \brief Used to sort orbits
  template<>
  class LocalSymCompare<Kinetics::DiffusionTransformation> :
    public Kinetics::DiffTransSymCompare<LocalSymCompare<Kinetics::DiffusionTransformation> > {

  public:

    typedef Kinetics::DiffusionTransformation Element;
    typedef Kinetics::DiffusionTransformationInvariants InvariantsType;

    LocalSymCompare(double tol);

  private:

    friend class SymCompare<Kinetics::DiffTransSymCompare<LocalSymCompare<Kinetics::DiffusionTransformation> > >;

    Element prepare_impl(const Element &A) const;

    bool compare_impl(const Element &A, const Element &B) const;

    bool invariants_compare_impl(const InvariantsType &A, const InvariantsType &B) const;

    /// \brief Apply symmetry to this
    ///
    /// - Affects no change
    void apply_sym_impl(const SymOp &op) {
      return;
    }

  };


  /// \brief Used to canonicalize DiffusionTransformations
  template<>
  class ScelPeriodicSymCompare<Kinetics::DiffusionTransformation> :
    public Kinetics::DiffTransSymCompare<ScelPeriodicSymCompare<Kinetics::DiffusionTransformation>> {

  public:

    typedef Kinetics::DiffusionTransformation Element;
    typedef Kinetics::DiffusionTransformationInvariants InvariantsType;

    ScelPeriodicSymCompare(const PrimGrid &prim_grid, double tol);

  private:

    friend class SymCompare<Kinetics::DiffTransSymCompare<ScelPeriodicSymCompare<Kinetics::DiffusionTransformation> > >;

    Element prepare_impl(const Element &A) const;

    bool compare_impl(const Element &A, const Element &B) const;

    bool invariants_compare_impl(const InvariantsType &A, const InvariantsType &B) const;

    /// \brief Apply symmetry to this
    ///
    /// - Affects no change
    void apply_sym_impl(const SymOp &op) {
      return;
    }

    const PrimGrid &m_prim_grid;

  };

  namespace Kinetics {
    typedef PrimPeriodicSymCompare<Kinetics::DiffusionTransformation> PrimPeriodicDiffTransSymCompare;
    typedef ScelPeriodicSymCompare<Kinetics::DiffusionTransformation> ScelPeriodicDiffTransSymCompare;
    typedef LocalSymCompare<Kinetics::DiffusionTransformation> LocalDiffTransSymCompare;
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
      const std::map<AtomSpecie, Index> &specie_count() const;

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

      std::map<AtomSpecie, Index> _from_specie_count() const;
      std::map<AtomSpecie, Index> _to_specie_count() const;
      std::map<AtomSpecie, Index> _empty_specie_count() const;


      std::vector<OccupationTransformation> m_occ_transform;
      std::vector<SpecieTrajectory> m_specie_traj;

      // stores IntegralCluster, based on occ_transform uccoord
      mutable notstd::cloneable_ptr<IntegralCluster> m_cluster;

      // stores Specie -> count, using 'from' specie
      // - is equal to 'to' specie count if is_valid_occ_transform() == true
      mutable notstd::cloneable_ptr<std::map<AtomSpecie, Index> > m_specie_count;

    };

    /// \brief Print DiffusionTransformation to stream, using default Printer<Kinetics::DiffusionTransformation>
    std::ostream &operator<<(std::ostream &sout, const DiffusionTransformation &trans);

    /// \brief Return a standardized name for this diffusion transformation orbit
    //std::string orbit_name(const PrimPeriodicDiffTransOrbit &orbit);

    // \brief Returns the distance from uccoord to the closest point on a linearly
    /// interpolated diffusion path. (Could be an end point)
    double dist_to_path(const DiffusionTransformation &diff_trans, const UnitCellCoord &uccoord);

    // \brief Returns the vector from uccoord to the closest point on a linearly
    /// interpolated diffusion path. (Could be an end point)
    Eigen::Vector3d vector_to_path(const DiffusionTransformation &diff_trans, const UnitCellCoord &uccoord);

    /// \brief Determines which site is closest to the diffusion transformation and the vector to take it to the path
    std::pair<UnitCellCoord, Eigen::Vector3d> _path_nearest_neighbor(const DiffusionTransformation &diff_trans) ;

    /// \brief Determines which site is closest to the diffusion transformation
    UnitCellCoord path_nearest_neighbor(const DiffusionTransformation &diff_trans);

    /// \brief Determines the nearest site distance to the diffusion path
    double min_dist_to_path(const DiffusionTransformation &diff_trans);

    /// \brief Determines the vector from the nearest site to the diffusion path in cartesian coordinates
    Eigen::Vector3d min_vector_to_path(const DiffusionTransformation &diff_trans);

    /// \brief Determines whether the atoms moving in the diffusion transformation will collide on a linearly interpolated path
    bool path_collision(const DiffusionTransformation &diff_trans);

  }

  typedef Orbit<Kinetics::DiffusionTransformation, Kinetics::PrimPeriodicDiffTransSymCompare> PrimPeriodicDiffTransOrbit;

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
