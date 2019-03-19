#include "casm/kinetics/OccPerturbation.hh"

#include "casm/symmetry/ScelOrbitGeneration_impl.hh"
#include "casm/clusterography/ClusterSymCompare_impl.hh"
#include "casm/clex/HasCanonicalForm_impl.hh"

namespace CASM {

  template struct ScelIsCanonical<OccPerturbation>;
  template bool CanonicalForm<ElementWiseSymApply<Kinetics::DoFTransformation<GenericCoordCluster<CRTPBase<CASM::OccPerturbation> > > > >::
  is_canonical<std::vector<PermuteIterator>::iterator>(
    Supercell const &,
    std::vector<PermuteIterator>::iterator,
    std::vector<PermuteIterator>::iterator) const;
  template bool CanonicalForm<ElementWiseSymApply<Kinetics::DoFTransformation<GenericCoordCluster<CRTPBase<CASM::OccPerturbation> > > > >::
  is_canonical<std::vector<PermuteIterator>::const_iterator>(
    Supercell const &,
    std::vector<PermuteIterator>::const_iterator,
    std::vector<PermuteIterator>::const_iterator) const;


  OccPerturbationInvariants::OccPerturbationInvariants(
    const OccPerturbation &perturb) :
    cluster_invariants(perturb.cluster().invariants()),
    from_species_count(CASM::from_species_count(perturb.begin(), perturb.end())),
    to_species_count(CASM::to_species_count(perturb.begin(), perturb.end())) {}

  /// \brief Check if DiffTransInvariants are equal
  bool almost_equal(const OccPerturbationInvariants &A, const OccPerturbationInvariants &B, double tol) {
    return almost_equal(A.cluster_invariants, B.cluster_invariants, tol)
           && A.from_species_count == B.from_species_count
           && A.to_species_count == B.to_species_count;
  }

  /// \brief Compare DiffTransInvariants
  bool compare(const OccPerturbationInvariants &A, const OccPerturbationInvariants &B, double tol) {
    if(compare(A.cluster_invariants, B.cluster_invariants, tol)) {
      return true;
    }
    if(compare(B.cluster_invariants, A.cluster_invariants, tol)) {
      return false;
    }
    if(A.from_species_count < B.from_species_count) {
      return true;
    }
    if(B.from_species_count < A.from_species_count) {
      return false;
    }
    return A.to_species_count < B.to_species_count;
  }

  /// \brief Print DiffTransInvariants
  std::ostream &operator<<(std::ostream &sout, const OccPerturbationInvariants &obj) {
    sout << obj.cluster_invariants;
    sout << " from(";
    if(obj.from_species_count.size() > 0) {
      for(const auto &t : obj.from_species_count) {
        sout << " " << t.first << ":" << t.second;
      }
    }
    sout << ") to(";
    if(obj.to_species_count.size() > 0) {
      for(const auto &t : obj.to_species_count) {
        sout << " " << t.first << ":" << t.second;
      }
    }
    sout << ")";
    return sout;
  }

  UnitCellCoord traits<OccPerturbation>::position(
    const OccPerturbation &perturb) {
    return perturb[0].uccoord;
  }


  /// \brief Construct an empty UnitCellCoordCluster
  OccPerturbation::OccPerturbation(const PrimType &_prim) :
    m_prim_ptr(&_prim) {}

  const OccPerturbation::PrimType &OccPerturbation::prim() const {
    return *m_prim_ptr;
  }

  /// \brief Access vector of elements
  std::vector<OccPerturbation::Element> &OccPerturbation::elements() {
    this->reset_invariants();
    return m_element;
  }

  /// \brief const Access vector of elements
  const std::vector<OccPerturbation::Element> &OccPerturbation::elements() const {
    return m_element;
  }

  /// \brief IntegralCluster as determined from sites in occ_transform()
  const IntegralCluster &OccPerturbation::cluster() const {
    if(!m_cluster) {
      m_cluster = notstd::make_cloneable<IntegralCluster>(prim());
      for(const auto &e : elements()) {
        m_cluster->elements().push_back(e.uccoord);
      }
    }
    return *m_cluster;
  }

  Configuration &OccPerturbation::apply_to(Configuration &config) const {
    for(auto it = this->begin(); it != this->end(); ++it) {
      it->apply_to(config);
    }
    return config;
  }

  Configuration &OccPerturbation::apply_reverse_to_impl(Configuration &config) const {
    for(auto it = this->begin(); it != this->end(); ++it) {
      it->apply_reverse_to(config);
    }
    return config;
  }

  void OccPerturbation::reverse() {
    OccPerturbation rev {*this};
    for(auto &occ_transform : rev) {
      occ_transform.reverse();
    }
  }

  Coordinate OccPerturbation::coordinate_impl(size_type i) const {
    return static_cast<Coordinate>(element(i).uccoord);
  }


  /// \brief Write OccupationTransformation to JSON object
  jsonParser &to_json(const OccPerturbation &perturb, jsonParser &json) {
    json.put_array(perturb.size());
    for(Index i = 0; i != perturb.size(); ++i) {
      to_json(perturb[i], json[i]);
    }
    return json;
  }

  OccPerturbation jsonConstructor<OccPerturbation>::from_json(const jsonParser &json, const Structure &prim) {
    OccPerturbation perturb(prim);
    for(Index i = 0; i < json.size(); ++i) {
      perturb.elements().push_back(
        jsonConstructor<Kinetics::OccupationTransformation>::from_json(json[i], prim));
    }
    return perturb;
  }

  const std::string Printer<OccPerturbation>::element_name = "OccPerturbation";

  /// \brief Print OccupationTransformation to stream, using default Printer<Kinetics::OccupationTransformation>
  std::ostream &operator<<(std::ostream &sout, const OccPerturbation &perturb) {
    Printer<OccPerturbation> printer;
    Log out(sout);
    printer.print(perturb, out);
    return sout;
  }

  void Printer<OccPerturbation>::print(const OccPerturbation &perturb, Log &out) {
    if(!out.print()) {
      return;
    }
    COORD_MODE printer_mode(this->opt.coord_type);

    Printer<Kinetics::OccupationTransformation> printer;
    for(const auto &trans : perturb) {
      printer.print(trans, out);
    }
  }
}
