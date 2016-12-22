#include "casm/kinetics/DoFTransformation.hh"

#include "casm/crystallography/Structure.hh"
#include "casm/clex/Configuration.hh"
#include "casm/symmetry/SymOp.hh"

namespace CASM {

  namespace Kinetics {

    // DoFTransformation

    DoFTransformation::DoFTransformation(const PrimType &prim) :
      m_prim(&prim) {};

    Configuration &DoFTransformation::apply_to(Configuration &config) const {
      return this->apply_to_impl(config);
    }

    Configuration &DoFTransformation::apply_reverse_to(Configuration &config) const {
      return this->apply_reverse_to_impl(config);
    }

    DoFTransformation &DoFTransformation::apply_sym(const SymOp &op) {
      this->apply_sym_impl(op);
      return *this;
    }

    void DoFTransformation::reverse() {
      this->reverse_impl();
    }

    std::unique_ptr<DoFTransformation> DoFTransformation::clone() const {
      return std::unique_ptr<DoFTransformation>(this->_clone());
    }


    // SiteDoFTransformation

    SiteDoFTransformation::SiteDoFTransformation(const UnitCellCoord &_uccoord) :
      DoFTransformation(_uccoord.unit()),
      uccoord(_uccoord) {}

    SiteDoFTransformation::~SiteDoFTransformation() {}

    SiteDoFTransformation &SiteDoFTransformation::operator+=(UnitCell frac) {
      uccoord += frac;
      return *this;
    }

    SiteDoFTransformation &SiteDoFTransformation::operator-=(UnitCell frac) {
      uccoord -= frac;
      return *this;
    }

    SiteDoFTransformation &SiteDoFTransformation::apply_sym(const SymOp &op) {
      this->apply_sym_impl(op);
      return *this;
    }

    std::unique_ptr<SiteDoFTransformation> SiteDoFTransformation::clone() const {
      return std::unique_ptr<SiteDoFTransformation>(this->_clone());
    }

    void SiteDoFTransformation::apply_sym_impl(const SymOp &op) {
      uccoord.apply_sym(op);
    }


    // OccupationTransformation

    OccupationTransformation::OccupationTransformation(
      const UnitCellCoord &uccoord,
      Index from_value,
      Index to_value) :
      SiteDoFTransformation(uccoord),
      from_value(from_value),
      to_value(to_value) {}

    OccupationTransformation &OccupationTransformation::apply_sym(const SymOp &op) {
      this->apply_sym_impl(op);
      return *this;
    }

    std::unique_ptr<OccupationTransformation> OccupationTransformation::clone() const {
      return std::unique_ptr<OccupationTransformation>(this->_clone());
    }

    Configuration &OccupationTransformation::apply_to_impl(Configuration &config) const {
      config.set_occ(config.linear_index(uccoord), to_value);
      return config;
    }

    Configuration &OccupationTransformation::apply_reverse_to_impl(Configuration &config) const {
      config.set_occ(config.linear_index(uccoord), from_value);
      return config;
    }

    void OccupationTransformation::apply_sym_impl(const SymOp &op) {
      static_cast<SiteDoFTransformation &>(*this).apply_sym(op);

      //MOLECULE_SUPPORT: apply permutation to from_value & to_value
    }

    void OccupationTransformation::reverse_impl() {
      using std::swap;
      swap(from_value, to_value);
    }

    OccupationTransformation *OccupationTransformation::_clone() const {
      return new OccupationTransformation(*this);
    }


    // SpecieLocation

    SpecieLocation::SpecieLocation(const UnitCellCoord &_uccoord, Index _occ, Index _pos) :
      uccoord(_uccoord),
      occ(_occ),
      pos(_pos) {}


    // SpecieTrajectory

    SpecieTrajectory::SpecieTrajectory(const SpecieLocation &_from,
                                       const SpecieLocation &_to) :
      from(_from),
      to(_to) {}

    SpecieTrajectory &SpecieTrajectory::operator+=(UnitCell frac) {
      from.uccoord += frac;
      to.uccoord += frac;
      return *this;
    }

    SpecieTrajectory &SpecieTrajectory::operator-=(UnitCell frac) {
      from.uccoord -= frac;
      to.uccoord -= frac;
      return *this;
    }

    void SpecieTrajectory::apply_sym(const SymOp &op) {
      from.uccoord.apply_sym(op);
      to.uccoord.apply_sym(op);

      //MOLECULE_SUPPORT: apply permutation to to/from_value & to/from_specie_index
    }

    void SpecieTrajectory::reverse() {
      using std::swap;
      swap(from, to);
    }


    // DiffusionTransformationInvariants

    DiffusionTransformationInvariants::DiffusionTransformationInvariants(
      const DiffusionTransformation &trans) :
      cluster_invariants(trans.cluster().invariants()),
      specie_count(trans.specie_count()) {}
  }

  /// \brief Check if DiffusionTransformationInvariants are equal
  bool almost_equal(const Kinetics::DiffusionTransformationInvariants &A, const Kinetics::DiffusionTransformationInvariants &B, double tol) {
    return almost_equal(A.cluster_invariants, B.cluster_invariants, tol) &&
           A.specie_count == B.specie_count;
  }

  /// \brief Compare DiffusionTransformationInvariants
  bool compare(const Kinetics::DiffusionTransformationInvariants &A, const Kinetics::DiffusionTransformationInvariants &B, double tol) {
    if(compare(A.cluster_invariants, B.cluster_invariants, tol)) {
      return true;
    }
    if(compare(B.cluster_invariants, A.cluster_invariants, tol)) {
      return false;
    }
    return A.specie_count < B.specie_count;
  }

  /// \brief Print DiffusionTransformationInvariants
  std::ostream &operator<<(std::ostream &sout, const Kinetics::DiffusionTransformationInvariants &obj) {
    sout << obj.cluster_invariants;
    if(obj.specie_count.size() > 0) {
      for(const auto &t : obj.specie_count) {
        sout << " " << t.first.name << ":" << t.second;
      }
    }
    return sout;
  }

  namespace Kinetics {

    // PrimPeriodicDiffTransSymCompare

    PrimPeriodicDiffTransSymCompare::PrimPeriodicDiffTransSymCompare(double tol) :
      SymCompare<PrimPeriodicDiffTransSymCompare>(),
      m_tol(tol) {}

    PrimPeriodicDiffTransSymCompare::Element PrimPeriodicDiffTransSymCompare::prepare_impl(const Element &A) const {
      Element tmp = A.prepared();
      return tmp -= tmp.occ_transform()[0].uccoord.unitcell();
    }

    bool PrimPeriodicDiffTransSymCompare::compare_impl(const Element &A, const Element &B) const {
      return A < B;
    }

    bool PrimPeriodicDiffTransSymCompare::invariants_compare_impl(const Element &A, const Element &B) const {
      return CASM::compare(A.invariants(), B.invariants(), tol());
    }


    // DiffusionTransformation

    DiffusionTransformation::DiffusionTransformation(const PrimType &prim) :
      DoFTransformation(prim) {
    }


    DiffusionTransformation &DiffusionTransformation::operator+=(UnitCell frac) {
      for(auto &t : m_occ_transform) {
        t += frac;
      }

      for(auto &t : m_specie_traj) {
        t += frac;
      }
      return *this;
    }

    DiffusionTransformation &DiffusionTransformation::operator-=(UnitCell frac) {
      for(auto &t : m_occ_transform) {
        t -= frac;
      }

      for(auto &t : m_specie_traj) {
        t -= frac;
      }
      return *this;
    }

    std::unique_ptr<DiffusionTransformation> DiffusionTransformation::clone() const {
      return std::unique_ptr<DiffusionTransformation>(this->_clone());
    }

    /// \brief Check if valid occupation transform
    ///
    /// Returns true if:
    /// - number of species of each type remains constant
    bool DiffusionTransformation::is_valid_occ_transform() const {
      return _from_specie_count() == _to_specie_count();
    }

    /// \brief Check if valid specie trajectories
    ///
    /// Returns true if:
    /// - species map to the correct specie type for occupation values,
    /// - no indivisble molecules are broken up,
    /// - some change occurs on every unitcell site (not a sub-hopcluster)
    bool DiffusionTransformation::is_valid_specie_traj() const {
      if(!_specie_types_map()) {
        return false;
      }
      if(_breaks_indivisible_mol()) {
        return false;
      }
      if(_is_subcluster_transformation()) {
        return false;
      }
      return true;
    }

    /// \brief Check if any SpecieTrajectory maps a Specie onto the wrong type
    bool DiffusionTransformation::_specie_types_map() const {
      return std::all_of(
               specie_traj().begin(),
               specie_traj().end(),
      [ = ](const SpecieTrajectory & t) {
        return t.specie_types_map();
      });
    }

    /// \brief Check if any indivisible Molecules are broken
    bool DiffusionTransformation::_breaks_indivisible_mol() const {

      // sort by 'from' uccoord (this may not typically be necessary, but let's be safe)
      auto tmp = specie_traj();
      std::sort(tmp.begin(), tmp.end());

      // if species from the same 'from' molecule end up on different 'to' molecule (checked via uccoord),
      //   and the 'from' molecule is indivisible -> true

      auto f = [ = ](const SpecieTrajectory & A, const SpecieTrajectory & B) {
        return A.from.uccoord == B.from.uccoord && A.to.uccoord != B.to.uccoord && A.from.mol().is_indivisible();
      };

      return std::adjacent_find(tmp.begin(), tmp.end(), f) != tmp.end();
    }

    /// \brief Check if removing any site would leave the DiffusionTransformation unchanged
    ///
    /// - Is a subcluster if there is any unitcell site on which no species change positions
    bool DiffusionTransformation::_is_subcluster_transformation() const {

      // sort SpecieTrajectory by 'from' molecule
      typedef std::map<UnitCellCoord, std::vector<SpecieTrajectory> > map_type;
      map_type m;
      for(const auto &t : m_specie_traj) {
        m[t.from.uccoord].push_back(t);
      }

      // lambda checks if all SpecieTrajectory from a molecule are 'no_change' trajectories
      auto is_no_change_mol = [ = ](const map_type::value_type & v) {
        return std::all_of(v.second.begin(), v.second.end(), [ = ](const SpecieTrajectory & t) {
          return t.is_no_change();
        });
      };

      // check if any Molecule is unchanged
      return std::any_of(m.begin(), m.end(), is_no_change_mol);
    }

    /// \brief Check if occ_transform and specie_traj are valid
    bool DiffusionTransformation::is_valid() const {
      return is_valid_occ_transform() && is_valid_specie_traj();
    }

    /// \brief Check if in canonical form
    ///
    /// - Uses prim factor group
    /// - Being in canonical form does not imply that this is sorted
    bool DiffusionTransformation::is_canonical() const {
      return is_canonical(prim().factor_group());
    }

    /// \brief Check if in canonical form
    ///
    /// - Being in canonical form does not imply that this is sorted
    bool DiffusionTransformation::is_canonical(const SymGroup &g) const {
      return std::none_of(
               g.begin(),
               g.end(),
      [&](const SymOp & op) {
        return *this < this->copy_apply_sym(op);
      });
    }

    /// \brief Returns the operation that applied to *this returns the canonical form
    ///
    /// - Uses prim factor group
    SymOp DiffusionTransformation::to_canonical() const {
      return to_canonical(prim().factor_group());
    }

    /// \brief Returns the operation that applied to *this returns the canonical form
    SymOp DiffusionTransformation::to_canonical(const SymGroup &g) const {
      auto it = std::max_element(
                  g.begin(),
                  g.end(),
      [&](const SymOp & op_A, const SymOp & op_B) {
        return this->copy_apply_sym(op_A) < this->copy_apply_sym(op_B);
      });
      return *it;
    }

    /// \brief Returns the operation that applied to the the canonical form returns *this
    ///
    /// - Uses prim factor group
    SymOp DiffusionTransformation::from_canonical() const {
      return from_canonical(prim().factor_group());
    }

    /// \brief Returns the operation that applied to the the canonical form returns *this
    SymOp DiffusionTransformation::from_canonical(const SymGroup &g) const {
      return to_canonical(g).inverse();
    }

    /// \brief Return canonical form
    ///
    /// - Uses prim factor group
    DiffusionTransformation DiffusionTransformation::canonical_form() const {
      return canonical_form(prim().factor_group());
    }

    /// \brief Return canonical form
    DiffusionTransformation DiffusionTransformation::canonical_form(const SymGroup &g) const {
      return this->copy_apply_sym(to_canonical(g));
    }

    std::vector<OccupationTransformation> &DiffusionTransformation::occ_transform() {
      _reset_invariants();
      return m_occ_transform;
    }

    const std::vector<OccupationTransformation> &DiffusionTransformation::occ_transform() const {
      return m_occ_transform;
    }

    std::vector<SpecieTrajectory> &DiffusionTransformation::specie_traj() {
      _reset_invariants();
      return m_specie_traj;
    }

    const std::vector<SpecieTrajectory> &DiffusionTransformation::specie_traj() const {
      return m_specie_traj;
    }

    /// \brief IntegralCluster as determined from sites in occ_transform()
    const IntegralCluster &DiffusionTransformation::cluster() const {
      if(!m_cluster) {
        m_cluster = notstd::make_cloneable<IntegralCluster>(prim());
        for(const auto &t : m_occ_transform) {
          m_cluster->elements().push_back(t.uccoord);
        }
      }
      return *m_cluster;
    }

    /// \brief Specie count as determined from sites in occ_transform()
    ///
    /// - Uses occ_transform() 'from' specie
    /// - Is equal to 'to' specie count if is_valid_occ_transform() == true
    const std::map<Specie, Index> &DiffusionTransformation::specie_count() const {
      if(!m_specie_count) {
        m_specie_count = notstd::make_cloneable<std::map<Specie, Index> >(_from_specie_count());
      }
      return *m_specie_count;
    }

    /// \brief Compare DiffusionTransformation
    ///
    /// - lexicographic comparison of [size, occ_transform, specie_traj], for the sorted
    ///   versions of this and B.
    bool DiffusionTransformation::operator<(const DiffusionTransformation &B) const {
      return this->prepared()._lt(B.prepared());
    }

    /// \brief Puts this in a prepared form, to prepare for comparisons
    ///
    /// - the forward and reverse occ_transform and specie_traj are sorted in
    ///   ascending order
    /// - this becomes the minimum of the forward and reverse
    ///
    DiffusionTransformation &DiffusionTransformation::prepare() {
      _forward_sort();
      DiffusionTransformation rev {*this};
      rev.reverse();
      rev._forward_sort();
      if(rev._lt(*this)) {
        *this = rev;
      }
      return *this;
    }

    /// \brief Returns the sorted version of this
    ///
    DiffusionTransformation DiffusionTransformation::prepared() const {
      DiffusionTransformation tmp {*this};
      return tmp.prepare();
    }

    DiffusionTransformation DiffusionTransformation::copy_apply_sym(const SymOp &op) const {
      DiffusionTransformation t {*this};
      t.apply_sym(op);
      return t;
    }

    Configuration &DiffusionTransformation::apply_to_impl(Configuration &config) const {

      // create the final specie id vectors in a temporary map

      // map of 'to' linear index -> 'to' specie_id
      std::map<Index, std::vector<Index> > _specie_id;

      for(const auto &t : m_specie_traj) {

        // linear indices of 'from' and 'to' sites
        Index from_l = config.linear_index(t.from.uccoord);
        Index to_l = config.linear_index(t.to.uccoord);

        // if 'to' linear index not yet in _specie_id, construct with correct length
        auto it = _specie_id.find(to_l);
        if(it == _specie_id.end()) {
          it = _specie_id.insert(std::make_pair(to_l, std::vector<Index>(t.to.mol().size()))).first;
        }

        // copy the specie id
        it->second[t.to.pos] = config.specie_id(from_l)[t.from.pos];
      }

      // transform the occupation variables
      for(const auto &t : m_occ_transform) {
        t.apply_to(config);
      }

      // copy the temporary specie_id
      for(const auto &t : _specie_id) {
        config.specie_id(t.first) = t.second;
      }
      return config;
    }

    Configuration &DiffusionTransformation::apply_reverse_to_impl(Configuration &config) const {
      // create the final specie id vectors in a temporary map

      // map of 'from' linear index -> 'from' specie_id
      std::map<Index, std::vector<Index> > _specie_id;

      for(const auto &t : m_specie_traj) {

        // linear indices of 'from' and 'to' sites
        Index from_l = config.linear_index(t.from.uccoord);
        Index to_l = config.linear_index(t.to.uccoord);

        // if 'from' linear index not yet in _specie_id, construct with correct length
        auto it = _specie_id.find(from_l);
        if(it == _specie_id.end()) {
          it = _specie_id.insert(std::make_pair(from_l, std::vector<Index>(t.from.mol().size()))).first;
        }

        // copy the specie id
        it->second[t.from.pos] = config.specie_id(to_l)[t.to.pos];
      }

      // transform the occupation variables
      for(const auto &t : m_occ_transform) {
        t.apply_reverse_to(config);
      }

      // copy the temporary specie_id
      for(const auto &t : _specie_id) {
        config.specie_id(t.first) = t.second;
      }
      return config;
    }

    void DiffusionTransformation::apply_sym_impl(const SymOp &op) {
      for(auto &t : m_occ_transform) {
        t.apply_sym(op);
      }

      for(auto &t : m_specie_traj) {
        t.apply_sym(op);
      }
    }

    void DiffusionTransformation::reverse_impl() {
      for(auto &t : m_occ_transform) {
        t.reverse();
      }

      for(auto &t : m_specie_traj) {
        t.reverse();
      }
    }

    DiffusionTransformation *DiffusionTransformation::_clone() const {
      return new DiffusionTransformation(*this);
    }

    /// \brief Puts this in a sorted form, without considering the reverse
    void DiffusionTransformation::_forward_sort() {
      std::sort(occ_transform().begin(), occ_transform().end());
      std::sort(specie_traj().begin(), specie_traj().end());
    }

    /// \brief Comparison of this and B, without sorting or considering reverse
    bool DiffusionTransformation::_lt(const DiffusionTransformation &B) const {
      if(occ_transform().size() < B.occ_transform().size()) {
        return true;
      }
      if(occ_transform().size() > B.occ_transform().size()) {
        return false;
      }

      {
        auto it = occ_transform().begin();
        auto B_it = B.occ_transform().begin();
        for(; it != occ_transform().end(); ++it, ++B_it) {
          if(*it < *B_it) {
            return true;
          }
          if(*it > *B_it) {
            return false;
          }
        }
      }

      {
        auto it = specie_traj().begin();
        auto B_it = B.specie_traj().begin();
        for(; it != specie_traj().end(); ++it, ++B_it) {
          if(*it < *B_it) {
            return true;
          }
          if(*it > *B_it) {
            return false;
          }
        }
      }
      return false;
    }

    /// \brief Reset mutable members, cluster and invariants, when necessary
    void DiffusionTransformation::_reset() {
      m_cluster.reset();
      _reset_invariants();
      m_specie_count.reset();
    }

    std::map<Specie, Index> DiffusionTransformation::_from_specie_count() const {
      std::map<Specie, Index> _specie_count = _empty_specie_count();
      for(const auto &t : m_occ_transform) {
        const Molecule &mol = t.uccoord.site().site_occupant()[t.from_value];
        for(const AtomPosition &specie_pos : mol) {
          _specie_count[specie_pos.specie]++;
        }
      }
      return _specie_count;
    }

    std::map<Specie, Index> DiffusionTransformation::_to_specie_count() const {
      std::map<Specie, Index> _specie_count = _empty_specie_count();
      for(const auto &t : m_occ_transform) {
        const Molecule &mol = t.uccoord.site().site_occupant()[t.to_value];
        for(const AtomPosition &specie_pos : mol) {
          _specie_count[specie_pos.specie]++;
        }
      }
      return _specie_count;
    }

    std::map<Specie, Index> DiffusionTransformation::_empty_specie_count() const {
      auto struc_specie = prim().struc_specie();
      std::map<Specie, Index> _specie_count;
      for(const Specie &s : struc_specie) {
        _specie_count[s] = 0;
      }
      return _specie_count;
    }

  }
}

