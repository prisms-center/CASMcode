#ifndef CASM_DiffusionTransformation
#define CASM_DiffusionTransformation

#include "casm/misc/cloneable_ptr.hh"
#include "casm/misc/Comparisons.hh"
#include "casm/symmetry/SymCompare.hh"
#include "casm/crystallography/UnitCellCoord.hh"
#include "casm/clusterography/ClusterInvariants.hh"
#include "casm/clusterography/IntegralCluster.hh"
#include "casm/clex/HasPrimClex.hh"
#include "casm/clex/HasCanonicalForm.hh"
#include "casm/kinetics/DiffusionTransformationTraits.hh"
#include "casm/kinetics/DoFTransformation.hh"
#include "casm/kinetics/OccupationTransformation.hh"

namespace CASM {

  class AtomSpecie;
  class Structure;
  class Configuration;
  class SymOp;
  class jsonParser;
  class PermuteIterator;
  template<typename T> struct jsonConstructor;

  namespace Kinetics {

    /// \brief Specifies a particular specie
    struct SpecieLocation : public Comparisons<CRTPBase<SpecieLocation>> {

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

    /// \brief Print DiffTransInvariants
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
    class SpecieTrajectory : public Comparisons<CRTPBase<SpecieTrajectory>> {

    public:

      SpecieTrajectory(const SpecieLocation &_from, const SpecieLocation &_to);

      SpecieTrajectory &operator+=(UnitCell frac);

      SpecieTrajectory &operator-=(UnitCell frac);

      bool specie_types_map() const;

      bool is_no_change() const;

      SpecieLocation from;
      SpecieLocation to;

      bool operator<(const SpecieTrajectory &B) const;

      SpecieTrajectory &apply_sym(const SymOp &op);

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
    class DiffTransInvariants {

    public:

      DiffTransInvariants(const DiffusionTransformation &trans);

      ClusterInvariants<IntegralCluster> cluster_invariants;
      std::map<AtomSpecie, Index> specie_count;

    };
  }


  /// \brief Check if DiffTransInvariants are equal
  bool almost_equal(const Kinetics::DiffTransInvariants &A,
                    const Kinetics::DiffTransInvariants &B,
                    double tol);

  /// \brief Compare DiffTransInvariants
  bool compare(const Kinetics::DiffTransInvariants &A,
               const Kinetics::DiffTransInvariants &B,
               double tol);

  /// \brief Print DiffTransInvariants
  std::ostream &operator<<(std::ostream &sout,
                           const Kinetics::DiffTransInvariants &obj);


  namespace Kinetics {

    typedef DoFTransformation <
    CanonicalForm <
    Comparisons <
    Translatable <
    SymComparable <
    CRTPBase<DiffusionTransformation >>> >>> DiffTransBase;

    /// \brief Describes how species move
    class DiffusionTransformation : public DiffTransBase {

    public:

      DiffusionTransformation(const Structure &_prim);


      const Structure &prim() const;

      DiffusionTransformation &operator+=(UnitCell frac);

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

      DiffusionTransformation &apply_sym(const SymOp &op);

      DiffusionTransformation &apply_sym(const PermuteIterator &it);

      Configuration &apply_to(Configuration &config) const;

      void reverse();


    private:

      friend DiffTransBase;

      Configuration &apply_reverse_to_impl(Configuration &config) const;

      void _forward_sort();

      bool _lt(const DiffusionTransformation &B) const;

      /// \brief Reset mutable members, cluster and invariants, when necessary
      void _reset();

      std::map<AtomSpecie, Index> _from_specie_count() const;
      std::map<AtomSpecie, Index> _to_specie_count() const;

      const Structure *m_prim_ptr;

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
