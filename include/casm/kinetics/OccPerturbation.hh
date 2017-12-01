#ifndef CASM_OccPerturbation
#define CASM_OccPerturbation

#include "casm/clusterography/ClusterInvariants.hh"
#include "casm/clusterography/GenericCluster.hh"
#include "casm/clusterography/CoordCluster.hh"
#include "casm/clex/HasCanonicalForm.hh"
#include "casm/kinetics/DoFTransformation.hh"
#include "casm/kinetics/OccupationTransformation.hh"
#include "casm/kinetics/OccupationTransformationTraits.hh"

namespace CASM {

  class Configuration;

  /// \brief Invariants of a OccPerturbation
  class OccPerturbationInvariants {

  public:

    OccPerturbationInvariants(const OccPerturbation &perturb);

    ClusterInvariants<IntegralCluster> cluster_invariants;
    std::map<AtomSpecie, Index> from_specie_count;
    std::map<AtomSpecie, Index> to_specie_count;
  };

  /// \brief Check if DiffTransInvariants are equal
  bool almost_equal(const OccPerturbationInvariants &A,
                    const OccPerturbationInvariants &B,
                    double tol);

  /// \brief Compare DiffTransInvariants
  bool compare(const OccPerturbationInvariants &A,
               const OccPerturbationInvariants &B,
               double tol);

  /// \brief Print DiffTransInvariants
  std::ostream &operator<<(std::ostream &sout, const OccPerturbationInvariants &obj);


  class OccPerturbation : public
    CanonicalForm <
    ElementWiseSymApply <
    DoFTransformation <
    GenericCoordCluster <
    CRTPBase<OccPerturbation >>> >> {

  public:

    typedef Structure PrimType;
    typedef traits<OccPerturbation>::Element Element;
    typedef traits<OccPerturbation>::InvariantsType InvariantsType;
    typedef traits<OccPerturbation>::size_type size_type;

    /// \brief Construct an empty UnitCellCoordCluster
    explicit OccPerturbation(const PrimType &_prim);

    /// \brief Construct a CoordCluster with a range of CoordType
    template<typename InputIterator>
    OccPerturbation(const PrimType &_prim,
                    InputIterator _begin,
                    InputIterator _end) :
      m_prim_ptr(&_prim),
      m_element(_begin, _end) {}

    const PrimType &prim() const;

    /// \brief Access vector of elements
    std::vector<Element> &elements();

    /// \brief const Access vector of elements
    const std::vector<Element> &elements() const;

    const IntegralCluster &cluster() const;

    Configuration &apply_to(Configuration &config) const;

    void reverse();

  protected:

    friend GenericCoordCluster<CRTPBase<OccPerturbation>>;
    friend DoFTransformation<GenericCoordCluster<CRTPBase<OccPerturbation>>>;

    Configuration &apply_reverse_to_impl(Configuration &config) const;

    Coordinate coordinate_impl(size_type i) const;

  private:

    std::vector<OccupationTransformation> m_element;
    const PrimType *m_prim_ptr;

    // stores IntegralCluster, based on occ_transform uccoord
    mutable notstd::cloneable_ptr<IntegralCluster> m_cluster;

  };

  /// \brief Write OccPerturbation to JSON object
  jsonParser &to_json(const OccPerturbation &trans, jsonParser &json);

  template<>
  struct jsonConstructor<OccPerturbation> {

    static OccPerturbation from_json(const jsonParser &json, const Structure &prim);
  };

  /// \brief Print OccPerturbation to stream, using default Printer<OccPerturbation>
  std::ostream &operator<<(std::ostream &sout, const OccPerturbation &perturb);

  template<>
  struct Printer<OccPerturbation> : public PrinterBase {

    typedef OccPerturbation Element;
    static const std::string element_name;

    Printer(int _indent_space = 6, char _delim = '\n', COORD_TYPE _mode = INTEGRAL) :
      PrinterBase(_indent_space, _delim, _mode) {}

    void print(const Element &element, std::ostream &out);
  };
}

#endif
