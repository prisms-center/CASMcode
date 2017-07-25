#include "casm/kinetics/DiffusionTransformation.hh"

#include "casm/casm_io/jsonParser.hh"
#include "casm/casm_io/json_io/container.hh"
#include "casm/symmetry/Orbit_impl.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/crystallography/Molecule.hh"
#include "casm/basis_set/DoF.hh"
#include "casm/clex/Configuration.hh"
#include "casm/clex/NeighborList.hh"


namespace CASM {

  template class Orbit<Kinetics::DiffusionTransformation, Kinetics::PrimPeriodicDiffTransSymCompare>;

  namespace Kinetics {

    namespace {
      std::ostream &operator<<(std::ostream &sout, const std::map<AtomSpecie, Index> &count) {
        for(const auto &t : count) {
          sout << "  " << t.first.name() << ": " << t.second << std::endl;
        }
        return sout;
      }
    }

    // SpecieLocation

    SpecieLocation::SpecieLocation(const UnitCellCoord &_uccoord, Index _occ, Index _pos) :
      uccoord(_uccoord),
      occ(_occ),
      pos(_pos) {}

    bool SpecieLocation::operator<(const SpecieLocation &B) const {
      return _tuple() < B._tuple();
    }

    const Molecule &SpecieLocation::mol() const {
      return uccoord.sublat_site().site_occupant()[occ];
    }

    const AtomSpecie &SpecieLocation::specie() const {
      return mol().atom(pos).specie();
    }

    std::tuple<UnitCellCoord, Index, Index> SpecieLocation::_tuple() const {
      return std::make_tuple(uccoord, occ, pos);
    }

    /// \brief Print DiffusionTransformationInvariants
    std::ostream &operator<<(std::ostream &sout, const SpecieLocation &obj) {
      sout << obj.uccoord << " : " << obj.occ << " " << obj.pos;
      return sout;
    }

  }

  jsonParser &to_json(const Kinetics::SpecieLocation &obj, jsonParser &json) {
    json.put_obj();
    json["uccoord"] = obj.uccoord;
    json["occ"] = obj.occ;
    json["pos"] = obj.pos;
    return json;
  }

  Kinetics::SpecieLocation jsonConstructor<Kinetics::SpecieLocation>::from_json(const jsonParser &json, const Structure &prim) {
    return Kinetics::SpecieLocation {
      jsonConstructor<UnitCellCoord>::from_json(json["uccoord"], prim),
      json["occ"].get<Index>(),
      json["pos"].get<Index>()
    };
  }

  void from_json(Kinetics::SpecieLocation &obj, const jsonParser &json) {
    from_json(obj.uccoord, json["uccoord"]);
    from_json(obj.occ, json["occ"]);
    from_json(obj.pos, json["pos"]);
  }


  namespace Kinetics {


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

    bool SpecieTrajectory::specie_types_map() const {
      return from.specie() == to.specie();
    }

    bool SpecieTrajectory::is_no_change() const {
      return from == to;
    }

    bool SpecieTrajectory::operator<(const SpecieTrajectory &B) const {
      return _tuple() < B._tuple();
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

    std::tuple<SpecieLocation, SpecieLocation> SpecieTrajectory::_tuple() const {
      return std::make_tuple(from, to);
    }
  }

  jsonParser &to_json(const Kinetics::SpecieTrajectory &traj, jsonParser &json) {
    json.put_obj();
    to_json(traj.to, json["to"]);
    to_json(traj.from, json["from"]);
    return json;
  }

  Kinetics::SpecieTrajectory jsonConstructor<Kinetics::SpecieTrajectory>::from_json(const jsonParser &json, const Structure &prim) {
    return Kinetics::SpecieTrajectory {
      jsonConstructor<Kinetics::SpecieLocation>::from_json(json["from"], prim),
      jsonConstructor<Kinetics::SpecieLocation>::from_json(json["to"], prim)
    };
  }

  void from_json(Kinetics::SpecieTrajectory &traj, const jsonParser &json) {
    from_json(traj.from, json["from"]);
    from_json(traj.to, json["to"]);
  }


  namespace Kinetics {


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
        sout << " " << t.first.name() << ":" << t.second;
      }
    }
    return sout;
  }

  PrimPeriodicSymCompare<Kinetics::DiffusionTransformation>::PrimPeriodicSymCompare(double _tol) :
    Kinetics::DiffTransSymCompare<PrimPeriodicSymCompare<Kinetics::DiffusionTransformation> >(_tol) {}

  PrimPeriodicSymCompare<Kinetics::DiffusionTransformation>::Element
  PrimPeriodicSymCompare<Kinetics::DiffusionTransformation>::prepare_impl(const Element &A) const {
    if(A.occ_transform().size()) {
      Element tmp = A.sorted();
      m_integral_tau = -(tmp.occ_transform()[0].uccoord.unitcell());
      tmp -= tmp.occ_transform()[0].uccoord.unitcell();
      return tmp;
    }
    else {
      return A;
    }
  }

  bool PrimPeriodicSymCompare<Kinetics::DiffusionTransformation>::compare_impl(
    const Element &A,
    const Element &B) const {
    return A < B;
  }

  bool PrimPeriodicSymCompare<Kinetics::DiffusionTransformation>::invariants_compare_impl(
    const InvariantsType &A,
    const InvariantsType &B) const {
    return CASM::compare(A, B, tol());
  }

  // LocalDiffTransSymCompare

  LocalSymCompare<Kinetics::DiffusionTransformation>::LocalSymCompare(double _tol) :
    DiffTransSymCompare<LocalSymCompare<Kinetics::DiffusionTransformation> >(_tol) {}

  LocalSymCompare<Kinetics::DiffusionTransformation>::Element
  LocalSymCompare<Kinetics::DiffusionTransformation>::prepare_impl(const Element &A) const {
    if(A.occ_transform().size()) {
      Element tmp = A.sorted();
      return tmp;
    }
    else {
      return A;
    }
  }

  bool LocalSymCompare<Kinetics::DiffusionTransformation>::compare_impl(
    const Element &A,
    const Element &B) const {
    return A < B;
  }

  bool LocalSymCompare<Kinetics::DiffusionTransformation>::invariants_compare_impl(
    const InvariantsType &A,
    const InvariantsType &B) const {
    return CASM::compare(A, B, tol());
  }


  // ScelPeriodicDiffTransSymCompare

  ScelPeriodicSymCompare<Kinetics::DiffusionTransformation>::ScelPeriodicSymCompare(
    const PrimGrid &prim_grid,
    double _tol) :
    DiffTransSymCompare<ScelPeriodicSymCompare<Kinetics::DiffusionTransformation> >(_tol),
    m_prim_grid(prim_grid) {}

  ScelPeriodicSymCompare<Kinetics::DiffusionTransformation>::Element
  ScelPeriodicSymCompare<Kinetics::DiffusionTransformation>::prepare_impl(const Element &A) const {
    if(A.occ_transform().size()) {
      Element tmp = A.sorted();
      m_integral_tau = (m_prim_grid.within(tmp.occ_transform()[0].uccoord).unitcell())
                       - tmp.occ_transform()[0].uccoord.unitcell();
      tmp += m_integral_tau;
      return tmp;
    }
    else {
      return A;
    }
  }

  bool ScelPeriodicSymCompare<Kinetics::DiffusionTransformation>::compare_impl(
    const Element &A, const Element &B) const {
    return A < B;
  }

  bool ScelPeriodicSymCompare<Kinetics::DiffusionTransformation>::invariants_compare_impl(
    const InvariantsType &A,
    const InvariantsType &B) const {
    return CASM::compare(A, B, tol());
  }

  namespace Kinetics {

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
      return specie_types_map() && !breaks_indivisible_mol() && !is_subcluster_transformation();
    }

    /// \brief Check if any SpecieTrajectory maps a Specie onto the wrong type
    bool DiffusionTransformation::specie_types_map() const {
      return std::all_of(
               specie_traj().begin(),
               specie_traj().end(),
      [ = ](const SpecieTrajectory & t) {
        return t.specie_types_map();
      });
    }

    /// \brief Check if any indivisible Molecules are broken
    bool DiffusionTransformation::breaks_indivisible_mol() const {

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
    /// - Is a subcluster if there is any unitcell site on which no species change positions,
    ///   of if any site begins and ends as a vacancy
    bool DiffusionTransformation::is_subcluster_transformation() const {

      // check if any vacancy transforms into a "different" vacancy
      auto is_va_to_va = [&](const OccupationTransformation & t) {
        return t.from_mol().is_vacancy() && t.to_mol().is_vacancy();
      };

      if(std::any_of(m_occ_transform.begin(), m_occ_transform.end(), is_va_to_va)) {
        return true;
      }

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

      // check if any Molecule is completely unchanged
      return std::any_of(m.begin(), m.end(), is_no_change_mol);
    }

    /// \brief Check if occ_transform and specie_traj are valid
    bool DiffusionTransformation::is_valid() const {
      return is_valid_occ_transform() && is_valid_specie_traj();
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
    const std::map<AtomSpecie, Index> &DiffusionTransformation::specie_count() const {
      if(!m_specie_count) {
        m_specie_count = notstd::make_cloneable<std::map<AtomSpecie, Index> >(_from_specie_count());
      }
      return *m_specie_count;
    }

    /// \brief Compare DiffusionTransformation
    ///
    /// - lexicographic comparison of [size, occ_transform, specie_traj], for the sorted
    ///   versions of this and B.
    bool DiffusionTransformation::operator<(const DiffusionTransformation &B) const {
      return this->sorted()._lt(B.sorted());
    }

    /// \brief Puts this in a sorted form, to enable comparisons
    ///
    /// - the forward and reverse occ_transform and specie_traj are sorted in
    ///   ascending order
    /// - this becomes the minimum of the forward and reverse
    ///
    DiffusionTransformation &DiffusionTransformation::sort() {
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
    DiffusionTransformation DiffusionTransformation::sorted() const {
      DiffusionTransformation tmp {*this};
      return tmp.sort();
    }

    /// \brief Check if sorted
    bool DiffusionTransformation::is_sorted() const {
      DiffusionTransformation _tmp = sorted();
      return !this->_lt(_tmp) && !_tmp._lt(*this);
    }

    /// \brief Return the cluster size
    Index DiffusionTransformation::size() const {
      return cluster().size();
    }

    /// \brief Return the min pair distance, or 0.0 if size() <= 1
    double DiffusionTransformation::min_length() const {
      return cluster().min_length();
    }

    /// \brief Return the max pair distance, or 0.0 if size() <= 1
    double DiffusionTransformation::max_length() const {
      return cluster().max_length();
    }

    /*
    DiffusionTransformation DiffusionTransformation::copy_apply_sym(const SymOp &op) const {
      DiffusionTransformation t {*this};
      t.apply_sym(op);
      return t;
    }
    */
    /*
    std::string orbit_name(const PrimPeriodicDiffTransOrbit &orbit) {
      Structure prim(orbit.prototype().specie_traj().begin()->from.uccoord.unit());
      std::set<int> sublat_indices;
      for(int i = 0; i < prim.basis.size(); i++) {
        sublat_indices.insert(i);
      }

      // construct
      PrimNeighborList nlist(
        PrimNeighborList::make_weight_matrix(prim.lattice().lat_column_mat(), 10, TOL),
        sublat_indices.begin(),
        sublat_indices.end()
      );
      int prev_size;
      int post_size;
      std::map<int, std::map<int, std::map<int, int>>> unique_inds;
      //May need to sort first?
      DiffusionTransformation s_this = orbit.prototype().sorted();
      std::string result = "DT" + std::to_string(orbit.prototype().size());
      std::vector<int> totals = {0, 0, 0, 0};
      for(auto it = s_this.specie_traj().begin(); it != s_this.specie_traj().end(); it++) {
        prev_size = nlist.size();
        nlist.expand(it->from.uccoord);
        nlist.expand(it->to.uccoord);
        post_size = nlist.size();
        if(prev_size != post_size) {
          int idx = 0;
          for(auto n_it = nlist.begin(); n_it != nlist.end(); n_it++) {
            unique_inds[(*n_it)(0)][(*n_it)(1)][(*n_it)(2)] = idx;
            idx++;
          }
        }
        result += "_" + it->from.specie().name + "-" + std::to_string(unique_inds[it->from.uccoord.unitcell(0)]
                                                                      [it->from.uccoord.unitcell(1)]
                                                                      [it->from.uccoord.unitcell(2)] + it->from.uccoord.sublat() * nlist.size()) ;

        result += ",";
        result += std::to_string(unique_inds[it->to.uccoord.unitcell(0)]
                                 [it->to.uccoord.unitcell(1)]
                                 [it->to.uccoord.unitcell(2)] + it->to.uccoord.sublat() * nlist.size()) + "-";
      }
      return result;
    }*/

    /// \brief Returns the distance from uccoord to the closest point on a linearly
    /// interpolated diffusion path. (Could be an end point)
    double dist_to_path(const DiffusionTransformation &diff_trans, const UnitCellCoord &uccoord) {
      return vector_to_path(diff_trans, uccoord).norm();
    }

    /// \brief Returns the vector from uccoord to the closest point on a linearly
    /// interpolated diffusion path. (Could be an end point)
    Eigen::Vector3d vector_to_path(const DiffusionTransformation &diff_trans, const UnitCellCoord &uccoord) {
      double dist = std::numeric_limits<double>::max();
      Eigen::Vector3d result;
      for(auto it = diff_trans.specie_traj().begin(); it != diff_trans.specie_traj().end(); it++) {
        //vector from -> input
        Coordinate v1 = (uccoord.coordinate() - it->from.uccoord.coordinate());
        //vector from -> to
        Coordinate v2 = (it->to.uccoord.coordinate() - it->from.uccoord.coordinate());
        // projection of v1 onto v2
        Eigen::Vector3d v3 = v1.const_cart().dot(v2.const_cart()) / (v1.const_cart().norm()) / (v2.const_cart().norm()) * v2.const_cart();
        double curr_dist;
        Eigen::Vector3d curr_vec;
        //if v3 length is greater than v2 then input is closer to "to" than the path
        if(v3.norm() > v2.const_cart().norm()) {
          curr_vec = it->to.uccoord.coordinate().const_cart() - uccoord.coordinate().const_cart();
          curr_dist = curr_vec.norm();
        }
        //if v3 is in opposite direction of v2 then input is closer to "from" than the path
        else if(v3.dot(v2.const_cart()) < 0) {
          curr_vec = -v1.const_cart();
          curr_dist = curr_vec.norm();
        }
        else {
          // find magnitude of v1-v3 and set to current distance
          curr_vec = v3 - v1.const_cart();
          curr_dist = curr_vec.norm();
        }
        if(curr_dist < dist) {
          dist = curr_dist;
          result = curr_vec;
        }
      }
      return result;
    }

    /// \brief Determines the nearest site distance to the diffusion path
    std::pair<UnitCellCoord, Eigen::Vector3d> _path_nearest_neighbor(const DiffusionTransformation &diff_trans) {
      double dist = std::numeric_limits<double>::max();
      Eigen::Vector3d ret_vec;
      Structure prim(diff_trans.specie_traj().begin()->from.uccoord.unit());
      std::set<int> sublat_indices;
      for(int i = 0; i < prim.basis.size(); i++) {
        sublat_indices.insert(i);
      }
      UnitCellCoord ret_coord(prim);
      // construct
      PrimNeighborList nlist(
        PrimNeighborList::make_weight_matrix(prim.lattice().lat_column_mat(), 10, TOL),
        sublat_indices.begin(),
        sublat_indices.end()
      );
      UnitCell pos(1, 1, 1);
      for(auto it = diff_trans.specie_traj().begin(); it != diff_trans.specie_traj().end(); it++) {
        UnitCellCoord fromcoord = it->from.uccoord;
        UnitCellCoord tocoord = it->to.uccoord;

        nlist.expand(fromcoord);
        fromcoord += pos;
        nlist.expand(fromcoord);
        fromcoord -= pos;
        fromcoord -= pos;
        nlist.expand(fromcoord);
        nlist.expand(tocoord);
        tocoord += pos;
        nlist.expand(tocoord);
        tocoord -= pos;
        tocoord -= pos;
        nlist.expand(tocoord);
      }
      for(auto n_it = nlist.begin(); n_it != nlist.end(); n_it++) {
        for(int b = 0; b < prim.basis.size(); b++) {
          UnitCellCoord uccoord(prim, b, *n_it);
          bool in_diff_trans = false;
          for(auto it = diff_trans.specie_traj().begin(); it != diff_trans.specie_traj().end(); it++) {
            if(uccoord == it->from.uccoord || uccoord == it->to.uccoord) {
              in_diff_trans = true;
            }
          }

          if(!in_diff_trans) {
            double curr_dist = dist_to_path(diff_trans, uccoord);
            Eigen::Vector3d curr_vec = vector_to_path(diff_trans, uccoord);
            if(curr_dist < dist) {
              dist = curr_dist;
              ret_vec = curr_vec;
              ret_coord = uccoord;
            }
          }
        }
      }

      std::pair<UnitCellCoord, Eigen::Vector3d> pair(ret_coord, ret_vec);
      return pair;
    }

    /// \brief Determines which site is closest to the diffusion transformation
    UnitCellCoord path_nearest_neighbor(const DiffusionTransformation &diff_trans) {
      return _path_nearest_neighbor(diff_trans).first;
    }

    /// \brief Determines the nearest site distance to the diffusion path
    double min_dist_to_path(const DiffusionTransformation &diff_trans) {
      return _path_nearest_neighbor(diff_trans).second.norm();
    }

    /// \brief Determines the vector from the nearest site to the diffusion path
    Eigen::Vector3d min_vector_to_path(const DiffusionTransformation &diff_trans) {
      return _path_nearest_neighbor(diff_trans).second;
    }

    Configuration &DiffusionTransformation::apply_to_impl(Configuration &config) const {

      if(config.has_specie_id()) {
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

        // copy the temporary specie_id
        for(const auto &t : _specie_id) {
          config.specie_id(t.first) = t.second;
        }
      }

      // transform the occupation variables
      for(const auto &t : m_occ_transform) {
        t.apply_to(config);
      }


      return config;
    }

    Configuration &DiffusionTransformation::apply_reverse_to_impl(Configuration &config) const {

      if(config.has_specie_id()) {
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
        // copy the temporary specie_id
        for(const auto &t : _specie_id) {
          config.specie_id(t.first) = t.second;
        }
      }
      // transform the occupation variables
      for(const auto &t : m_occ_transform) {
        t.apply_reverse_to(config);
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

    std::map<AtomSpecie, Index> DiffusionTransformation::_from_specie_count() const {
      std::map<AtomSpecie, Index> _specie_count = _empty_specie_count();
      for(const auto &t : m_occ_transform) {
        const Molecule &mol = t.uccoord.sublat_site().site_occupant()[t.from_value];
        for(const AtomPosition &specie_pos : mol.atoms()) {
          _specie_count[specie_pos.specie()]++;
        }
      }
      return _specie_count;
    }

    std::map<AtomSpecie, Index> DiffusionTransformation::_to_specie_count() const {
      std::map<AtomSpecie, Index> _specie_count = _empty_specie_count();
      for(const auto &t : m_occ_transform) {
        const Molecule &mol = t.uccoord.sublat_site().site_occupant()[t.to_value];
        for(const AtomPosition &specie_pos : mol.atoms()) {
          _specie_count[specie_pos.specie()]++;
        }
      }
      return _specie_count;
    }

    std::map<AtomSpecie, Index> DiffusionTransformation::_empty_specie_count() const {
      auto struc_specie = prim().struc_specie();
      std::map<AtomSpecie, Index> _specie_count;
      for(const AtomSpecie &s : struc_specie) {
        _specie_count[s] = 0;
      }
      return _specie_count;
    }

    /// \brief Print DiffusionTransformation to stream, using default Printer<Kinetics::DiffusionTransformation>
    std::ostream &operator<<(std::ostream &sout, const DiffusionTransformation &trans) {
      Printer<Kinetics::DiffusionTransformation> printer;
      printer.print(trans, sout);
      return sout;
    }
  }

  /// \brief Write DiffusionTransformation to JSON object
  jsonParser &to_json(const Kinetics::DiffusionTransformation &trans, jsonParser &json) {
    json.put_obj();
    json["occ_transform"].put_array(trans.occ_transform().begin(), trans.occ_transform().end());
    json["specie_traj"].put_array(trans.specie_traj().begin(), trans.specie_traj().end());
    return json;
  }

  Kinetics::DiffusionTransformation jsonConstructor<Kinetics::DiffusionTransformation>::from_json(const jsonParser &json, const Structure &prim) {
    Kinetics::DiffusionTransformation trans {prim};
    CASM::from_json(trans, json, prim);
    return trans;
  }

  /// \brief Read from JSON
  void from_json(Kinetics::DiffusionTransformation &trans, const jsonParser &json, const Structure &prim) {
    if(json["occ_transform"].size() > 0) {
      trans.occ_transform().clear();
      for(auto it = json["occ_transform"].begin(); it != json["occ_transform"].end(); ++it) {
        trans.occ_transform().push_back(jsonConstructor<Kinetics::OccupationTransformation>::from_json(*it, prim));
      }
      /*from_json(
        trans.occ_transform(),
        json["occ_transform"],
        json["occ_transform"][0].get<Kinetics::OccupationTransformation>(prim));*/
    }
    if(json["specie_traj"].size() > 0) {
      trans.specie_traj().clear();
      for(auto it = json["specie_traj"].begin(); it != json["specie_traj"].end(); ++it) {
        trans.specie_traj().push_back(jsonConstructor<Kinetics::SpecieTrajectory>::from_json(*it, prim));
      }
      /*from_json(
        trans.specie_traj(),
        json["specie_traj"],
        json["specie_traj"][0].get<Kinetics::SpecieTrajectory>(prim));*/
    }
  }

  const std::string Printer<Kinetics::DiffusionTransformation>::element_name = "DiffusionTransformation";

  void Printer<Kinetics::DiffusionTransformation>::print(const Kinetics::DiffusionTransformation &trans, std::ostream &out) {
    COORD_MODE printer_mode(mode);

    if(trans.is_valid()) {
      for(const auto &traj : trans.specie_traj()) {
        out << indent() << indent() << indent();
        out << traj.from.specie().name() + ": " << traj.from << "  ->  " << traj.to;

        if(delim)
          out << delim;
        out << std::flush;
      }
    }
    else {
      out << indent() << indent() << indent() << "occupation transformation:" << delim;
      for(const auto &t : trans.occ_transform()) {
        out << t;
      }
      out << indent() << indent() << indent() << "specie trajectory:" << delim;
      for(const auto &traj : trans.specie_traj()) {
        out << indent() << indent() << indent();
        out << traj.from << " (" << traj.from.specie().name() << ")";
        out << "  ->  ";
        out << traj.to << " (" << traj.to.specie().name() << ")";

        if(delim)
          out << delim;
        out << std::flush;
      }
    }
  }

}

