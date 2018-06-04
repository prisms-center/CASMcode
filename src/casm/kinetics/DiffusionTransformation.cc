#include "casm/kinetics/DiffusionTransformation_impl.hh"
#include "casm/clex/NeighborList.hh"
#include "casm/database/Named_impl.hh"
#include "casm/database/DiffTransOrbitDatabase.hh"

namespace CASM {

  namespace DB {
    template class Named<CRTPBase<PrimPeriodicDiffTransOrbit> >;
    template class Indexed<CRTPBase<PrimPeriodicDiffTransOrbit> >;
  }

  template class ClusterSymCompare<SymCompare<CRTPBase<AperiodicSymCompare<Kinetics::DiffusionTransformation> > > >;
  template class AperiodicSymCompare<Kinetics::DiffusionTransformation>;

  template class ClusterSymCompare<SymCompare<CRTPBase<PrimPeriodicSymCompare<Kinetics::DiffusionTransformation> > > >;
  template class PrimPeriodicSymCompare<Kinetics::DiffusionTransformation>;

  template class ClusterSymCompare<SymCompare<CRTPBase<ScelPeriodicSymCompare<Kinetics::DiffusionTransformation> > > >;
  template class ScelPeriodicSymCompare<Kinetics::DiffusionTransformation>;

  namespace Kinetics {

    namespace {
      std::ostream &operator<<(std::ostream &sout, const std::map<AtomSpecies, Index> &count) {
        for(const auto &t : count) {
          sout << "  " << t.first.name() << ": " << t.second << std::endl;
        }
        return sout;
      }
    }

    // SpecieLocation

    SpeciesLocation::SpeciesLocation(const UnitCellCoord &_uccoord, Index _occ, Index _pos) :
      uccoord(_uccoord),
      occ(_occ),
      pos(_pos) {}

    bool SpeciesLocation::operator<(const SpeciesLocation &B) const {
      return _tuple() < B._tuple();
    }

    const Molecule &SpeciesLocation::mol() const {
      return uccoord.sublat_site().site_occupant()[occ];
    }

    const AtomSpecies &SpeciesLocation::species() const {
      return mol().atom(pos).species();
    }

    std::tuple<UnitCellCoord, Index, Index> SpeciesLocation::_tuple() const {
      return std::make_tuple(uccoord, occ, pos);
    }

    /// \brief Print DiffTransInvariants
    std::ostream &operator<<(std::ostream &sout, const SpeciesLocation &obj) {
      sout << obj.uccoord << " : " << obj.occ << " " << obj.pos;
      return sout;
    }

  }

  jsonParser &to_json(const Kinetics::SpeciesLocation &obj, jsonParser &json) {
    json.put_obj();
    json["uccoord"] = obj.uccoord;
    json["occ"] = obj.occ;
    json["pos"] = obj.pos;
    return json;
  }

  Kinetics::SpeciesLocation jsonConstructor<Kinetics::SpeciesLocation>::from_json(const jsonParser &json, const Structure &prim) {
    return Kinetics::SpeciesLocation {
      jsonConstructor<UnitCellCoord>::from_json(json["uccoord"], prim),
      json["occ"].get<Index>(),
      json["pos"].get<Index>()
    };
  }

  void from_json(Kinetics::SpeciesLocation &obj, const jsonParser &json) {
    from_json(obj.uccoord, json["uccoord"]);
    from_json(obj.occ, json["occ"]);
    from_json(obj.pos, json["pos"]);
  }


  namespace Kinetics {


    // SpecieTrajectory

    SpecieTrajectory::SpecieTrajectory(const SpeciesLocation &_from,
                                       const SpeciesLocation &_to) :
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

    bool SpecieTrajectory::species_types_map() const {
      return from.species() == to.species();
    }

    bool SpecieTrajectory::is_no_change() const {
      return from == to;
    }

    bool SpecieTrajectory::operator<(const SpecieTrajectory &B) const {
      return _tuple() < B._tuple();
    }

    SpecieTrajectory &SpecieTrajectory::apply_sym(const SymOp &op) {
      from.uccoord.apply_sym(op);
      to.uccoord.apply_sym(op);

      //MOLECULE_SUPPORT: apply permutation to to/from_value & to/from_species_index
      return *this;
    }

    void SpecieTrajectory::reverse() {
      using std::swap;
      swap(from, to);
    }

    std::tuple<SpeciesLocation, SpeciesLocation> SpecieTrajectory::_tuple() const {
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
      jsonConstructor<Kinetics::SpeciesLocation>::from_json(json["from"], prim),
      jsonConstructor<Kinetics::SpeciesLocation>::from_json(json["to"], prim)
    };
  }

  void from_json(Kinetics::SpecieTrajectory &traj, const jsonParser &json) {
    from_json(traj.from, json["from"]);
    from_json(traj.to, json["to"]);
  }


  namespace Kinetics {


    // DiffTransInvariants

    DiffTransInvariants::DiffTransInvariants(
      const DiffusionTransformation &trans) :
      cluster_invariants(trans.cluster().invariants()),
      species_count(trans.species_count()) {}
  }

  /// \brief Check if DiffTransInvariants are equal
  bool almost_equal(const Kinetics::DiffTransInvariants &A, const Kinetics::DiffTransInvariants &B, double tol) {
    return almost_equal(A.cluster_invariants, B.cluster_invariants, tol) &&
           A.species_count == B.species_count;
  }

  /// \brief Compare DiffTransInvariants
  bool compare(const Kinetics::DiffTransInvariants &A, const Kinetics::DiffTransInvariants &B, double tol) {
    if(compare(A.cluster_invariants, B.cluster_invariants, tol)) {
      return true;
    }
    if(compare(B.cluster_invariants, A.cluster_invariants, tol)) {
      return false;
    }
    return A.species_count < B.species_count;
  }

  /// \brief Print DiffTransInvariants
  std::ostream &operator<<(std::ostream &sout, const Kinetics::DiffTransInvariants &obj) {
    sout << obj.cluster_invariants;
    if(obj.species_count.size() > 0) {
      for(const auto &t : obj.species_count) {
        sout << " " << t.first.name() << ":" << t.second;
      }
    }
    return sout;
  }

  namespace Kinetics {

    // DiffusionTransformation

    DiffusionTransformation::DiffusionTransformation(const Structure &_prim) :
      m_prim_ptr(&_prim) {
    }

    const Structure &DiffusionTransformation::prim() const {
      return *m_prim_ptr;
    }

    DiffusionTransformation &DiffusionTransformation::operator+=(UnitCell frac) {
      m_cluster.reset();
      for(auto &t : m_occ_transform) {
        t += frac;
      }

      for(auto &t : m_species_traj) {
        t += frac;
      }
      return *this;
    }

    /// \brief Check if valid occupation transform
    ///
    /// Returns true if:
    /// - number of species of each type remains constant
    bool DiffusionTransformation::is_valid_occ_transform() const {
      return _from_species_count() == _to_species_count();
    }

    /// \brief Check if valid specie trajectories
    ///
    /// Returns true if:
    /// - species map to the correct specie type for occupation values,
    /// - no indivisble molecules are broken up,
    /// - some change occurs on every unitcell site (not a sub-hopcluster)
    bool DiffusionTransformation::is_valid_species_traj() const {
      return species_types_map() && !breaks_indivisible_mol() && !is_subcluster_transformation();
    }

    /// \brief Check if any SpecieTrajectory maps a Specie onto the wrong type
    bool DiffusionTransformation::species_types_map() const {
      return std::all_of(
               species_traj().begin(),
               species_traj().end(),
      [ = ](const SpecieTrajectory & t) {
        return t.species_types_map();
      });
    }

    /// \brief Check if any indivisible Molecules are broken
    bool DiffusionTransformation::breaks_indivisible_mol() const {

      // sort by 'from' uccoord (this may not typically be necessary, but let's be safe)
      auto tmp = species_traj();
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
      for(const auto &t : m_species_traj) {
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

    /// \brief Check if species_traj() and occ_transform() are consistent
    ///
    /// - Checks if species_traj occ indices match occ_transform indices
    ///   and that there as many traj as AtomSpecie in a Molecule
    ///
    bool DiffusionTransformation::is_self_consistent() const {
      for(const auto &trans : occ_transform()) {
        auto is_from_match = [&](const SpecieTrajectory & traj) {
          return trans.uccoord == traj.from.uccoord && trans.from_mol() == traj.from.mol();
        };
        auto from_match_count = std::count_if(species_traj().begin(), species_traj().end(), is_from_match);
        if(from_match_count != trans.from_mol().size()) {
          return false;
        }

        auto is_to_match = [&](const SpecieTrajectory & traj) {
          return trans.uccoord == traj.to.uccoord && trans.to_mol() == traj.to.mol();
        };
        auto to_match_count = std::count_if(species_traj().begin(), species_traj().end(), is_to_match);
        if(from_match_count != trans.to_mol().size()) {
          return false;
        }
      }
      return true;
    }

    /// \brief Check if occ_transform and species_traj are valid and self consistent
    bool DiffusionTransformation::is_valid() const {
      return is_valid_occ_transform() && is_valid_species_traj() && is_self_consistent();
    }

    std::vector<OccupationTransformation> &DiffusionTransformation::occ_transform() {
      _reset();
      return m_occ_transform;
    }

    const std::vector<OccupationTransformation> &DiffusionTransformation::occ_transform() const {
      return m_occ_transform;
    }

    std::vector<SpecieTrajectory> &DiffusionTransformation::species_traj() {
      _reset();
      return m_species_traj;
    }

    const std::vector<SpecieTrajectory> &DiffusionTransformation::species_traj() const {
      return m_species_traj;
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
    const std::map<AtomSpecies, Index> &DiffusionTransformation::species_count() const {
      if(!m_species_count) {
        m_species_count = notstd::make_cloneable<std::map<AtomSpecies, Index> >(_from_species_count());
      }
      return *m_species_count;
    }

    /// \brief Compare DiffusionTransformation
    ///
    /// - lexicographic comparison of [size, occ_transform, species_traj], for the sorted
    ///   versions of this and B.
    bool DiffusionTransformation::operator<(const DiffusionTransformation &B) const {
      return this->sorted()._lt(B.sorted());
    }

    Permutation DiffusionTransformation::sort_permutation() const {
      throw std::runtime_error("DiffusionTransformation::sort_permutation() not yet implemented.");
      return Permutation(0);

    }

    /// \brief Puts this in a sorted form, to enable comparisons
    ///
    /// - the forward and reverse occ_transform and species_traj are sorted in
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

    Configuration &DiffusionTransformation::apply_to(Configuration &config) const {
      // transform the occupation variables
      for(const auto &t : m_occ_transform) {
        t.apply_to(config);
      }
      return config;
    }

    DiffusionTransformation &DiffusionTransformation::apply_sym(const SymOp &op) {
      m_cluster.reset();
      for(auto &t : m_occ_transform) {
        t.apply_sym(op);
      }

      for(auto &t : m_species_traj) {
        t.apply_sym(op);
      }
      return *this;
    }

    DiffusionTransformation &DiffusionTransformation::apply_sym(const PermuteIterator &it) {
      apply_sym(it.sym_op());
      return *this;
    }

    void DiffusionTransformation::reverse() {
      for(auto &t : m_occ_transform) {
        t.reverse();
      }

      for(auto &t : m_species_traj) {
        t.reverse();
      }
    }

    Configuration &DiffusionTransformation::apply_reverse_to_impl(Configuration &config) const {
      // transform the occupation variables
      for(const auto &t : m_occ_transform) {
        t.apply_reverse_to(config);
      }
      return config;
    }

    /// \brief Puts this in a sorted form, without considering the reverse
    void DiffusionTransformation::_forward_sort() {
      std::sort(occ_transform().begin(), occ_transform().end());
      std::sort(species_traj().begin(), species_traj().end());
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
        auto it = species_traj().begin();
        auto B_it = B.species_traj().begin();
        for(; it != species_traj().end(); ++it, ++B_it) {
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
      reset_invariants();
      m_species_count.reset();
    }

    std::map<AtomSpecies, Index> DiffusionTransformation::_from_species_count() const {
      return from_species_count(m_occ_transform.begin(), m_occ_transform.end());
    }

    std::map<AtomSpecies, Index> DiffusionTransformation::_to_species_count() const {
      return to_species_count(m_occ_transform.begin(), m_occ_transform.end());
    }

    /// \brief Print DiffusionTransformation to stream, using default Printer<Kinetics::DiffusionTransformation>
    std::ostream &operator<<(std::ostream &sout, const DiffusionTransformation &trans) {
      Printer<Kinetics::DiffusionTransformation> printer;
      Log out(sout);
      printer.print(trans, out);
      return sout;
    }

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
      for(auto it = diff_trans.species_traj().begin(); it != diff_trans.species_traj().end(); it++) {
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
      Structure prim(diff_trans.species_traj().begin()->from.uccoord.unit());
      std::set<int> sublat_indices;
      for(int i = 0; i < prim.basis().size(); i++) {
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
      for(auto it = diff_trans.species_traj().begin(); it != diff_trans.species_traj().end(); ++it) {
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
        for(int b = 0; b < prim.basis().size(); b++) {
          UnitCellCoord uccoord(prim, b, *n_it);
          bool in_diff_trans = false;
          for(auto it = diff_trans.species_traj().begin(); it != diff_trans.species_traj().end(); it++) {
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

    /// \brief Determines whether the atoms moving in the diffusion transformation will collide on a linearly interpolated path
    bool path_collision(const DiffusionTransformation &diff_trans) {
      std::vector<SpecieTrajectory> paths_to_check;
      for(auto it = diff_trans.species_traj().begin(); it != diff_trans.species_traj().end(); ++it) {
        if(!is_vacancy(it->from.species().name())) {
          paths_to_check.push_back(*it);
        }
      }
      for(int i = 0; i < paths_to_check.size(); ++i) {
        for(int j = i + 1; j < paths_to_check.size(); ++j) {
          //vector from -> to path 1
          Eigen::Vector3d v1 = (paths_to_check[i].to.uccoord.coordinate() - paths_to_check[i].from.uccoord.coordinate()).const_cart();
          //vector from -> to path 2
          Eigen::Vector3d v2 = (paths_to_check[j].to.uccoord.coordinate() - paths_to_check[j].from.uccoord.coordinate()).const_cart();
          // simplification of the following problem
          // parametric representation of path 1 = parametric representation of path 2
          // paths_to_check[i].from.uccoord.coordinate() + t*v1 = paths_to_check[j].from.uccoord.coordinate() + s*v2
          // v3 =  paths_to_check[j].from.uccoord.coordinate() - paths_to_check[i].from.uccoord.coordinate()
          Eigen::Vector3d v3 = (paths_to_check[j].from.uccoord.coordinate() - paths_to_check[i].from.uccoord.coordinate()).const_cart();
          Eigen::MatrixXd soln(2, 1);
          Eigen::Matrix2d m;
          Eigen::MatrixXd b(2, 1);
          m << v1[0], v2[0],
          v1[1], v2[1];
          b << v3[0],
          v3[1];
          int excluded_index = 2;
          try {
            Eigen::FullPivLU<Eigen::Matrix2d> lu(m);
            if(lu.rank() < 2) {
              m << v1[1], v2[1],
              v1[2], v2[2];
              b << v3[1],
              v3[2];
              excluded_index = 0;
              Eigen::FullPivLU<Eigen::Matrix2d> lu2(m);
              if(lu2.rank() < 2) {
                return true;
              }
            }
            soln = m.inverse() * b;
            if(soln(0, 0)*v1[excluded_index] + soln(1, 0)*v2[excluded_index] != v3[excluded_index]) {
              continue; //if solution for x and y coordinates but not z solution is not viable
            }
            if(soln(0, 0) < 1 && soln(0, 0) > 0 && soln(1, 0) < 0 && soln(1, 0) > -1) {
              return true;  // there is an intersection of paths not at the end points
            }

            if(soln(0, 0) == 1 || soln(0, 0) == 0 || soln(1, 0) == 0 || soln(1, 0) == -1) {
              if(v1 == -v2) {
                return true; // the paths are overlapping in opposition directions perfectly
              }
            }
          }
          catch(...) {
            throw;
          }
        }
      }
      return false;
    }
  }

  /// \brief Write DiffusionTransformation to JSON object
  jsonParser &to_json(const Kinetics::DiffusionTransformation &trans, jsonParser &json) {
    json.put_obj();
    json["occ_transform"].put_array(trans.occ_transform().begin(), trans.occ_transform().end());
    json["species_traj"].put_array(trans.species_traj().begin(), trans.species_traj().end());
    return json;
  }

  Kinetics::DiffusionTransformation jsonConstructor<Kinetics::DiffusionTransformation>::from_json(const jsonParser &json, const Structure &prim) {
    Kinetics::DiffusionTransformation trans {prim};
    CASM::from_json(trans, json, prim);
    return trans;
  }

  Kinetics::DiffusionTransformation jsonConstructor<Kinetics::DiffusionTransformation>::from_json(
    const jsonParser &json, const PrimClex &primclex) {
    Kinetics::DiffusionTransformation trans {primclex.prim()};
    CASM::from_json(trans, json, primclex.prim());
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
    if(json["species_traj"].size() > 0) {
      trans.species_traj().clear();
      for(auto it = json["species_traj"].begin(); it != json["species_traj"].end(); ++it) {
        trans.species_traj().push_back(jsonConstructor<Kinetics::SpecieTrajectory>::from_json(*it, prim));
      }
      /*from_json(
        trans.species_traj(),
        json["species_traj"],
        json["species_traj"][0].get<Kinetics::SpecieTrajectory>(prim));*/
    }
  }

  const std::string Printer<Kinetics::DiffusionTransformation>::element_name = "DiffusionTransformation";

  void Printer<Kinetics::DiffusionTransformation>::print(const Kinetics::DiffusionTransformation &trans, Log &out) const {
    if(!out.print()) {
      return;
    }

    COORD_MODE printer_mode(mode);

    if(trans.is_valid()) {
      if(mode != INTEGRAL) {
        // calculate nice widths
        int name_width = 0;
        int prec = 7;
        int width = prec;
        out.ostream().precision(prec);
        out.ostream().flags(std::ios::showpoint | std::ios::fixed | std::ios::right);
        for(const auto &traj : trans.species_traj()) {
          if(traj.from.species().name().length() > name_width) name_width = traj.from.species().name().length();
          Eigen::Vector3d vec_from, vec_to;
          if(mode == CART) {
            vec_from = traj.from.uccoord.coordinate().cart();
            vec_to = traj.to.uccoord.coordinate().cart();
          }
          else {
            vec_from = traj.from.uccoord.coordinate().frac();
            vec_to = traj.to.uccoord.coordinate().frac();
          }
          width = print_matrix_width(out, vec_from, width);
          width = print_matrix_width(out, vec_to, width);
        }

        // print
        Eigen::IOFormat format(prec, width + 1);
        for(const auto &traj : trans.species_traj()) {
          out << out.indent_str() << indent();
          out << std::setw(name_width) << traj.from.species().name() << ": ";
          {
            const auto &obj = traj.from;
            obj.uccoord.coordinate().print(out, 0, format);
            out << " : " << obj.occ << " " << obj.pos;
          }
          out << "  ->  ";
          {
            const auto &obj = traj.to;
            obj.uccoord.coordinate().print(out, 0, format);
            out << " : " << obj.occ << " " << obj.pos;
          }
          if(delim)
            out << delim;
          out << std::flush;
        }
      }
      else {
        // calculate nice widths
        int name_width = 0;
        int prec = 1;
        int width = prec;
        out.ostream().precision(prec);
        out.ostream().flags(std::ios::showpoint | std::ios::fixed | std::ios::right);
        for(const auto &traj : trans.species_traj()) {
          if(traj.from.species().name().length() > name_width) name_width = traj.from.species().name().length();
          width = print_matrix_width(out, traj.from.uccoord.unitcell(), width);
          width = print_matrix_width(out, traj.to.uccoord.unitcell(), width);
        }

        // print
        Eigen::IOFormat format(prec, width);
        for(const auto &traj : trans.species_traj()) {
          out << out.indent_str() << indent();
          out << std::setw(name_width) << traj.from.species().name() << ": ";
          {
            const auto &obj = traj.from;
            out << obj.uccoord.sublat() << ", " << obj.uccoord.unitcell().transpose().format(format) << " : " << obj.occ << " " << obj.pos;
          }
          out << "  ->  ";
          {
            const auto &obj = traj.to;
            out << obj.uccoord.sublat() << ", " << obj.uccoord.unitcell().transpose().format(format) << " : " << obj.occ << " " << obj.pos;
          }
          if(delim)
            out << delim;
          out << std::flush;
        }
      }
    }
    else {
      out << out.indent_str() << indent() << "occupation transformation:" << delim;
      for(const auto &t : trans.occ_transform()) {
        out << t;
      }
      out << out.indent_str() << indent() << "specie trajectory:" << delim;
      for(const auto &traj : trans.species_traj()) {
        out << out.indent_str() << indent();
        out << traj.from << " (" << traj.from.species().name() << ")";
        out << "  ->  ";
        out << traj.to << " (" << traj.to.species().name() << ")";

        if(delim)
          out << delim;
        out << std::flush;
      }
    }
  }

}

