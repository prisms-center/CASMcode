#include "casm/crystallography/StrucMapping.hh"

#include "casm/crystallography/SuperlatticeEnumerator.hh"
#include "casm/crystallography/StrucMapCalculatorInterface.hh"
#include "casm/container/algorithm.hh"
#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/crystallography/Strain.hh"
#include "casm/crystallography/Lattice.hh"
#include "casm/crystallography/Lattice_impl.hh"
#include "casm/crystallography/SymTools.hh"
#include "casm/crystallography/Coordinate.hh"
#include "casm/crystallography/CanonicalForm.hh"
#include "casm/crystallography/LatticeMap.hh"

namespace CASM {
  namespace xtal {

    namespace Local {

      static bool lex_lt(Eigen::Matrix<long, 3, 3> const &A, Eigen::Matrix<long, 3, 3> const &B) {
        return std::lexicographical_compare(A.data(), A.data() + 9, B.data(), B.data() + 9);
      }
    }
    //*******************************************************************************************

    namespace StrucMapping {
      double strain_cost(double relaxed_lat_vol, const MappingNode &mapped_result, const Index Nsites) {
        return StrainCostCalculator::iso_strain_cost(mapped_result.lat_node.stretch, relaxed_lat_vol / double(max(Nsites, Index(1))));
      }

      //*******************************************************************************************

      double basis_cost(const MappingNode &mapped_result, Index Nsites) {
        // mean square displacement distance in deformed coordinate system
        return (mapped_result.lat_node.stretch.inverse() * mapped_result.atom_displacement).squaredNorm() / double(max(Nsites, Index(1)));
      }
    }

    //*******************************************************************************************
    namespace Local {
      // Local helper function for StrucMapper::k_best_maps_better_than
      template<typename OutputIterator>
      static bool initial_basis_maps(SimpleStructure child_struc,
                                     MappingNode const &seed,
                                     StrucMapCalculatorInterface const &calculator,
                                     double max_cost,
                                     OutputIterator it) {

        //derotate first
        child_struc.rotate_coords(seed.isometry());

        //Then undeform by inverse of right stretch
        child_struc.deform_coords(seed.stretch());

        // We want to get rid of translations.
        // define translation such that:
        //    IDEAL = RELAXED + translation
        // and use it when calculating cost matrix

        for(Eigen::Vector3d const &translation : calculator.translations(seed, child_struc)) {

          MappingNode node = seed;
          node.basis_node.translation = translation;
          if(!calculator.populate_cost_mat(node,
                                           child_struc)) {
            // Indicates that structure is incompatible with supercell, regardless of translation so return false
            return false;
          }

          // The mapping routine is called here
          node.calc();

          // if assignment is smaller than child_struc.basis().size(), then child_struc is incompattible with supercell
          // (assignment.size()==0 if the hungarian routine detects an incompatibility, regardless of translation)
          if(!node.is_viable) {
            return false;
          }
          // Now we are filling up displacements
          calculator.finalize(node, child_struc);

          // add small penalty (~_tol) for larger translation distances, so that shortest equivalent translation is used
          if(node.cost < max_cost) {
            *it = node;
          }

        }

        return true;
      }

      //*******************************************************************************************
      // Local helper function for StrucMapper::k_best_maps_better_than
      template<typename OutputIterator>
      static void partition_node(MappingNode const &_node,
                                 StrucMapCalculatorInterface const &_calculator,
                                 SimpleStructure child_struc,
                                 OutputIterator it) {
        //derotate first
        child_struc.rotate_coords(_node.isometry());

        //Then undeform by inverse of right stretch
        child_struc.deform_coords(_node.stretch());

        Index cN = _node.lat_node.child.size() * _calculator.struc_info(child_struc).size();

        Index j, jj, currj, realj;
        Index n = _node.basis_node.assignment.size();
        MappingNode t1(_node),
                    t2(_node);

        _node.is_partitioned = true;
        t1.is_partitioned = false;
        t2.is_partitioned = false;

        MappingNode *p1 = &t1;
        MappingNode *p2 = &t2;
        for(Index m = n; 1 < m && (p1->is_viable || p2->is_viable); --m) {
          AssignmentNode &n1(p1->basis_node);
          AssignmentNode &n2(p2->basis_node);
          //clear assignment and cost_mat
          n2.assignment.clear();
          n2.irow.clear();
          n2.icol.clear();
          n2.cost_mat.resize(m - 1, m - 1);
          n2.forced_on = n1.forced_on;

          // We are forcing on the first site assignment in t1: i.e.,  [0,currj]
          // this involves striking row 0 and col currj from t2's cost_mat
          // and augmenting the cost_offset of t2 by the cost of [0,currj]
          currj = n1.assignment[0];
          realj = n1.icol[currj];
          n2.cost_offset = n1.cost_offset + n1.cost_mat(0, currj);



          // [0,currj] have local context only. We store the forced assignment in forced_on
          // using the original indexing from n1
          n2.forced_on.emplace(n1.irow[0], n1.icol[currj]);

          // Strike row 0 and col currj to form new AssignmentNode for t2
          // (i,j) indexes the starting cost_mat, (i-1, jj) indexes the resulting cost_mat
          n2.irow = std::vector<Index>(++n1.irow.begin(), n1.irow.end());
          // We will also store an updated assignment vector in t2, which will be
          // used to construct next node of partition
          n2.assignment = std::vector<Index>(++n1.assignment.begin(), n1.assignment.end());
          for(j = 0, jj = 0; j < m; ++j, ++jj) {
            if(j == currj) {
              --jj;
              continue;
            }
            n2.icol.push_back(n1.icol[j]);

            // We will also store an updated assignment vector in t2, which will be
            // used to construct next node of partition
            if(n2.assignment[jj] > currj)
              n2.assignment[jj]--;

            // Fill col jj of t2's cost mat
            for(Index i = 1; i < m; ++i)
              n2.cost_mat(i - 1, jj) = n1.cost_mat(i, j);
          }
          //t2 properly initialized; we can now force OFF [0,currj] in t1, and add it to node list
          n1.cost_mat(0, currj) = StrucMapping::big_inf();
          // IMPORTANT: If realj>=cN, it is a virtual vacancy.
          // If we exclude a single virtual vacancy from occupying this parent site, the marginal cost
          // of assigning a different virtual vacancy to the same site is zero, and the mapping is equivalent.
          // So, we need to exclude ALL virtual vacancies from occupying this parent site:
          if(realj >= cN) {
            for(j = n1.icol.size() - 1; j >= 0 && n1.icol[j] >= cN; --j) {
              n1.cost_mat(0, j) = StrucMapping::big_inf();
            }
          }
          n1.assignment.clear();
          p1->is_viable = true;
          p1->calc();
          if(p1->is_viable) {
            //even if p1 is unviable, p2 may still be viable, so we continue
            _calculator.finalize(*p1, child_struc);
            it = *p1;
          }
          std::swap(p1, p2);
        }
      }
    }

    //*******************************************************************************************

    LatticeNode::LatticeNode(Lattice const &parent_prim,
                             Lattice const &parent_scel,
                             Lattice const &child_prim,
                             Lattice const &child_scel,
                             Index child_N_atom,
                             double _cost /*=StrucMapping::big_inf()*/) :
      parent(parent_prim, parent_scel),
      // Transform child_prim lattice to its idealized state using same F.inverse
      // as below, but inline:
      child(Lattice((parent_scel.lat_column_mat() * child_scel.inv_lat_column_mat())
                    * child_prim.lat_column_mat(), parent_prim.tol()),
            parent_scel),
      cost(_cost) {

      // F is from ideal parent to child ( child_scel = F * parent_scel )
      Eigen::Matrix3d F = child_scel.lat_column_mat() * parent_scel.inv_lat_column_mat();

      // stretch is from (de-rotated, strained) child to ideal parent
      // child_scel = F * parent_scel = isometry.transpose() * stretch.inverse() * parent_scel
      // OR: parent_scel = stretch * isometry * child_scel
      stretch = strain::right_stretch_tensor(F).inverse();

      // isometry is from child to strained parent
      isometry = (F * stretch).transpose();

      if(StrucMapping::is_inf(cost))
        cost = StrainCostCalculator::iso_strain_cost(stretch, child_prim.volume() / double(max(child_N_atom, Index(1))));
    }

    //*******************************************************************************************

    LatticeNode::LatticeNode(LatticeMap const &_lat_map,
                             Lattice const &parent_prim,
                             Lattice const &child_prim) :
      // stretch is from (de-rotated, strained) child to ideal parent
      stretch(polar_decomposition(_lat_map.matrixF()).inverse()),
      // isometry is from child to strained parent
      isometry((_lat_map.matrixF() * stretch).transpose()),
      parent(parent_prim, Lattice(_lat_map.parent_matrix(), parent_prim.tol())),
      child(Lattice(_lat_map.matrixF().inverse() * child_prim.lat_column_mat(), parent_prim.tol()), Lattice(_lat_map.parent_matrix(), parent_prim.tol())),
      cost(_lat_map.strain_cost()) {
    }

    //*******************************************************************************************

    /// \brief Compare two LatticeMap objects, based on their mapping cost
    bool less(LatticeNode const &A, LatticeNode const &B, double cost_tol) {
      if(!almost_equal(A.cost, B.cost, cost_tol))
        return A.cost < B.cost;
      if(A.child.transformation_matrix() != B.child.transformation_matrix())
        return Local::lex_lt(A.child.transformation_matrix(), B.child.transformation_matrix());
      if(A.parent.transformation_matrix() != B.parent.transformation_matrix())
        return Local::lex_lt(A.parent.transformation_matrix(), B.parent.transformation_matrix());
      return false;
    }

    //*******************************************************************************************

    bool identical(LatticeNode const &A, LatticeNode const &B, double cost_tol) {
      if(!almost_equal(A.cost, B.cost, cost_tol))
        return false;
      if(A.parent.transformation_matrix() != B.parent.transformation_matrix())
        return false;
      if(A.child.transformation_matrix() != B.child.transformation_matrix())
        return false;
      return true;
    }

    //*******************************************************************************************

    bool AssignmentNode::operator <(AssignmentNode const &B) const {
      if(empty() != B.empty())
        return empty();
      if(time_reversal != B.time_reversal)
        return B.time_reversal;
      if(!almost_equal(translation, B.translation, 1e-6))
        return float_lexicographical_compare(translation, B.translation, 1e-6);
      return false;
    }

    //*******************************************************************************************

    bool identical(AssignmentNode const &A, AssignmentNode const &B) {
      if(A.empty() != B.empty())
        return false;
      if(A.time_reversal != B.time_reversal)
        return false;
      if(!almost_equal(A.translation, B.translation, 1e-6))
        return false;
      return true;
    }
    //*******************************************************************************************

    MappingNode MappingNode::invalid() {
      static MappingNode result(LatticeNode(Lattice::cubic(),
                                            Lattice::cubic(),
                                            Lattice::cubic(),
                                            Lattice::cubic(),
                                            1),
                                0.5);
      result.is_viable = false;
      result.is_valid = false;
      result.is_partitioned = false;
      return result;
    }

    //*******************************************************************************************

    void MappingNode::calc() {

      if(is_viable) {
        if(basis_node.irow.empty())
          basis_node.irow = sequence<Index>(0, basis_node.cost_mat.rows() - 1);
        if(basis_node.icol.empty())
          basis_node.icol = sequence<Index>(0, basis_node.cost_mat.cols() - 1);
        basis_node.cost = hungarian_method(basis_node.cost_mat, basis_node.assignment, cost_tol()) + basis_node.cost_offset;
        if(StrucMapping::is_inf(basis_node.cost)) {
          is_viable = false;
          cost = StrucMapping::big_inf();
        }
        else {
          cost = strain_weight * lat_node.cost + basis_weight * basis_node.cost;
        }
      }
      else
        cost = StrucMapping::big_inf();

    }

    bool MappingNode::operator<(MappingNode const &B) const {
      if(!almost_equal(cost, B.cost)) {
        return cost < B.cost;
      }
      if(!almost_equal(lat_node.cost, B.lat_node.cost)) {
        return lat_node.cost < B.lat_node.cost;
      }
      if(basis_node.empty() != B.basis_node.empty()) {
        return basis_node.empty();
      }
      if(!identical(lat_node, B.lat_node)) {
        return lat_node < B.lat_node;
      }
      if(!identical(basis_node, B.basis_node)) {
        return basis_node < B.basis_node;
      }
      if(atom_permutation != B.atom_permutation)
        return std::lexicographical_compare(atom_permutation.begin(),
                                            atom_permutation.end(),
                                            B.atom_permutation.begin(),
                                            B.atom_permutation.end());

      return false;
    }

    //*******************************************************************************************

    StrucMapper::StrucMapper(StrucMapCalculatorInterface const &calculator,
                             double _strain_weight /*= 0.5*/,
                             double _max_volume_change /*= 0.5*/,
                             int _options /*= robust*/, // this should actually be a bitwise-OR of StrucMapper::Options
                             double _cost_tol /*= TOL*/,
                             double _min_va_frac /*= 0.*/,
                             double _max_va_frac /*= 1.*/) :
      m_calc_ptr(calculator.clone()),
      m_max_volume_change(_max_volume_change),
      m_options(_options),
      m_cost_tol(max(1e-10, _cost_tol)),
      m_xtal_tol(TOL),
      m_filtered(false) {

      set_min_va_frac(_min_va_frac);
      set_max_va_frac(_max_va_frac);

      //squeeze strain_weight into (0,1] if necessary
      set_strain_weight(_strain_weight);

      //Make sure that max_volume_change is positive
      m_max_volume_change = max(3 * xtal_tol(), _max_volume_change);
    }

    //*******************************************************************************************

    Index StrucMapper::_n_species(SimpleStructure const &sstruc) const {
      return calculator().struc_info(sstruc).size();
    }

    //*******************************************************************************************

    std::pair<Index, Index> StrucMapper::_vol_range(const SimpleStructure &child_struc) const {
      Index min_vol(0), max_vol(0);
      //mapped_result.clear();
      SimpleStructure::Info const &c_info(calculator().struc_info(child_struc));

      if(calculator().fixed_species().size() > 0) {
        std::string tcompon = calculator().fixed_species().begin()->first;
        int ncompon(0);
        for(std::string const &sp : c_info.names) {
          if(sp == tcompon)
            ncompon++;
        }
        min_vol = ncompon / int(calculator().fixed_species().begin()->second);
        max_vol = min_vol;
      }
      else {

        // Try to narrow the range of supercell volumes -- the best bounds are obtained from
        // the convex hull of the end-members, but we need to wait for improvements to convex hull
        // routines

        int max_n_va = calculator().max_n_va();

        // Using _n_species() instead of parent().n_mol() here.
        // For mapping molecules, may need to push this entire routine into the StrucMapCalculator?
        double N_sites = double(_n_species(parent()));

        // Absolute largest Va fraction
        double max_va_frac_limit = double(max_n_va) / N_sites;

        double t_min_va_frac = min(min_va_frac(), max_va_frac_limit);
        double t_max_va_frac = min(max_va_frac(), max_va_frac_limit);

        // Use vacancy range to narrow volume range. Ceiling function is used because the floating-point
        // volume is a lower-bound value in both cases.

        // the number of sites implied in order to realize va_frac >= min_va_frac (TOL is not super important):
        double min_n_sites = floor(_n_species(child_struc) / (1. - t_min_va_frac) - TOL);

        // the number of sites implied to realize va_frac <= max_va_frac  (TOL is not super important):
        double max_n_sites = ceil(_n_species(child_struc) / (1. - t_max_va_frac) + TOL);

        // min_vol assumes min number vacancies -- best case scenario (buffer by half an atom)
        min_vol = ceil((min_n_sites - 0.5) / N_sites);

        // max_vol assumes min number vacancies -- best case scenario (buffer by half an atom)
        max_vol = floor((max_n_sites + 0.5) / N_sites);

        // If volume is not fully determined, use volume ratio of child to parent to narrow search
        if(min_vol != max_vol) {
          //Nvol is approximate integer volume-- assume that actual integer is within m_max_volume_change*100% of this volume, and use it to tighten our bounds
          double Nvol = (std::abs(child_struc.lat_column_mat.determinant() / parent().lat_column_mat.determinant()));
          min_vol = min(max_vol, max<Index>(round((1.0 - m_max_volume_change) * double(Nvol)), min_vol));
          max_vol = max(min_vol, min<Index>(round((1.0 + m_max_volume_change) * double(Nvol)), max_vol));
        }
      }

      // make sure volume range is above zero
      min_vol = max<Index>(min_vol, 1);
      max_vol = max<Index>(max_vol, 1);

      return std::pair<Index, Index>(min_vol, max_vol);

    }

    //*******************************************************************************************

    SimpleStructure const &StrucMapper::parent() const {
      return calculator().parent();
    }


    //*******************************************************************************************

    std::set<MappingNode> StrucMapper::_seed_from_vol_range(SimpleStructure const &child_struc,
                                                            Index k,
                                                            Index min_vol,
                                                            Index max_vol) const {

      if(!valid_index(min_vol) || !valid_index(min_vol) || max_vol < min_vol) {
        auto vol_range = _vol_range(child_struc);
        min_vol = vol_range.first;
        max_vol = vol_range.second;
      }
      Lattice child_lat(child_struc.lat_column_mat, xtal_tol());
      std::set<MappingNode> mapping_seed;
      for(Index i_vol = min_vol; i_vol <= max_vol; i_vol++) {
        std::vector<Lattice> lat_vec;
        for(Lattice const &lat : _lattices_of_vol(i_vol)) {
          if(m_filtered && !_filter_lat(lat, child_lat)) {
            continue;
          }
          lat_vec.push_back(lat);
        }

        std::set<MappingNode> t_seed = _seed_k_best_from_super_lats(child_struc,
                                                                    lat_vec,
        {Lattice(child_struc.lat_column_mat, xtal_tol())},
        k);

        mapping_seed.insert(std::make_move_iterator(t_seed.begin()), std::make_move_iterator(t_seed.end()));

      }
      return mapping_seed;
    }

    //*******************************************************************************************
    /*
     * Given a structure and a mapping node, find a perfect supercell of the prim that is equivalent to structure's lattice
     * and then try to map the structure's basis onto that supercell
     *
     * Returns false if no mapping is possible, or if the lattice is not ideal
     *
     * What this does NOT do:
     *    -Check if the imported Structure is the same as one in a smaller Supercell
     *
     */
    //*******************************************************************************************

    std::set<MappingNode> StrucMapper::map_ideal_struc(const SimpleStructure &child_struc,
                                                       Index k,
                                                       double max_cost /*=StrucMapping::big_inf()*/,
                                                       double min_cost /*=-TOL*/,
                                                       bool keep_invalid /*=false*/) const {

      // xtal::is_superlattice isn't very smart right now, and will return
      // false if the two lattices differ by a rigid rotation
      // In the future this may not be the case, so we will assume that child_struc may
      // be rigidly rotated relative to prim
      Eigen::Matrix3d trans_mat;

      // c_lat must be an ideal supercell of the parent lattice, but it need not be canonical
      // We will account for the difference in orientation between c_lat and the canonical supercell,
      // which must be related by a point group operation
      Lattice c_lat(child_struc.lat_column_mat, xtal_tol());


      bool c_lat_is_supercell_of_parent;
      std::tie(c_lat_is_supercell_of_parent, trans_mat) = xtal::is_superlattice(c_lat, Lattice(parent().lat_column_mat), xtal_tol());
      if(!c_lat_is_supercell_of_parent) {
        return {};
      }

      // We know child_struc.lattice() is a supercell of the prim, now we have to
      // reorient 'child_struc' by a point-group operation of the parent to match canonical lattice vectors
      // This may not be a rotation in the child structure's point group
      Lattice derot_c_lat(canonical::equivalent(Lattice(parent().lat_column_mat * trans_mat, xtal_tol()), calculator().point_group()));

      // We now find a transformation matrix of c_lat so that, after transformation, it is related
      // to derot_c_lat by rigid rotation only. Following line finds R and T such that derot_c_lat = R*c_lat*T
      auto res = xtal::is_equivalent_superlattice(derot_c_lat, c_lat, calculator().point_group().begin(), calculator().point_group().end(), xtal_tol());
      LatticeNode lat_node(Lattice(parent().lat_column_mat, xtal_tol()),
                           derot_c_lat,
                           c_lat,
                           Lattice(child_struc.lat_column_mat * res.second.cast<double>(), xtal_tol()),
                           _n_species(child_struc),
                           0. /*strain_cost is zero in ideal case*/);

      return map_deformed_struc_impose_lattice_node(child_struc, lat_node, k, max_cost, min_cost, keep_invalid);

    }

    //*******************************************************************************************

    std::set<MappingNode> StrucMapper::map_deformed_struc(const SimpleStructure &child_struc,
                                                          Index k /*=1*/,
                                                          double max_cost /*=StrucMapping::big_inf()*/,
                                                          double min_cost /*=-TOL*/,
                                                          bool keep_invalid /*=false*/) const {
      auto vols = _vol_range(child_struc);
      return map_deformed_struc_impose_lattice_vols(child_struc, vols.first, vols.second, k, max_cost, min_cost, keep_invalid);
    }

    //*******************************************************************************************

    std::set<MappingNode> StrucMapper::map_deformed_struc_impose_lattice_vols(const SimpleStructure &child_struc,
                                                                              Index min_vol,
                                                                              Index max_vol,
                                                                              Index k /*=1*/,
                                                                              double max_cost /*=StrucMapping::big_inf()*/,
                                                                              double min_cost /*=-TOL*/,
                                                                              bool keep_invalid /*=false*/) const {

      int seed_k = 10 + 5 * k;
      std::set<MappingNode> mapping_seed = _seed_from_vol_range(child_struc, seed_k, min_vol, max_vol);

      std::cout << "LatticeMap seed size: " << mapping_seed.size() << std::endl;

      bool no_partition = !(robust & options()) && k <= 1;
      k_best_maps_better_than(child_struc, mapping_seed, k, max_cost, min_cost, keep_invalid, false, no_partition);

      return mapping_seed;
    }

    //*******************************************************************************************

    std::set<MappingNode> StrucMapper::map_deformed_struc_impose_lattice(const SimpleStructure &child_struc,
                                                                         const Lattice &imposed_lat,
                                                                         Index k /*=1*/,
                                                                         double max_cost /*=StrucMapping::big_inf()*/,
                                                                         double min_cost /*=-TOL*/,
                                                                         bool keep_invalid /*=false*/) const {

      std::set<MappingNode> mapping_seed = _seed_k_best_from_super_lats(child_struc,
      {imposed_lat},
      {Lattice(child_struc.lat_column_mat, xtal_tol())},
      k,
      max(min_cost, 1e-5));
      bool no_partition = !(robust & options()) && k <= 1;
      k_best_maps_better_than(child_struc, mapping_seed, k, max_cost, min_cost, keep_invalid, false, no_partition);
      return mapping_seed;
    }

    //*******************************************************************************************

    std::set<MappingNode> StrucMapper::map_deformed_struc_impose_lattice_node(const SimpleStructure &child_struc,
                                                                              const LatticeNode &imposed_node,
                                                                              Index k /*=1*/,
                                                                              double max_cost /*=StrucMapping::big_inf()*/,
                                                                              double min_cost /*=-TOL*/,
                                                                              bool keep_invalid /*=false*/) const {

      std::set<MappingNode> mapping_seed;
      mapping_seed.emplace(imposed_node, m_strain_weight);
      bool no_partition = !(robust & options()) && k <= 1;
      k_best_maps_better_than(child_struc, mapping_seed, k, max_cost, min_cost, keep_invalid, false, no_partition);
      return mapping_seed;
    }

    //*******************************************************************************************

    std::vector<Lattice> StrucMapper::_lattices_of_vol(Index prim_vol) const {

      if(!valid_index(prim_vol)) {
        throw std::runtime_error("Cannot enumerate lattice of volume " + std::to_string(prim_vol) + ", which is out of bounds.\n");
      }

      //If you specified that you wanted certain lattices, return those, otherwise do the
      //usual enumeration
      if(this->lattices_constrained()) {
        //This may very well return an empty vector, saving painful time enumerating things
        return m_allowed_superlat_map[prim_vol];
      }

      //If we already have candidate lattices for the given volume, return those
      auto it = m_superlat_map.find(prim_vol);
      if(it != m_superlat_map.end())
        return it->second;

      //We don't have any lattices for the provided volume, enumerate them all!!!
      std::vector<Lattice> &lat_vec = m_superlat_map[prim_vol];

      auto pg = calculator().point_group();
      SuperlatticeEnumerator enumerator(pg.begin(), pg.end(), Lattice(parent().lat_column_mat, xtal_tol()),
                                        ScelEnumProps(prim_vol, prim_vol + 1));

      for(auto it = enumerator.begin(); it != enumerator.end(); ++it) {
        Lattice canon_lat = *it;
        if(canonical::check(canon_lat, calculator().point_group())) {
          canon_lat = canonical::equivalent(canon_lat, calculator().point_group());
        }
        lat_vec.push_back(canon_lat);
      }

      return lat_vec;
    }

    //*******************************************************************************************
    /*
     * Starting from a child structure and a queue of (possibly unfinalized) mapping nodes,
     * visit each node in the queue, test its fitness, and add any of its viable derivative nodes to the queue
     *
     * The queue is ordered by the lower bound of a node's mapping cost. This lower cost is less than or equal
     * to its finalized mapping cost or the finalized mapping costs of any of its derivatives node. As such,
     * the queue only needs to be iterated over once.
     *
     * At the end of the process, the queue is pruned to have at most 'k' valid nodes that are the
     * best mappings of child_struc onto the parent structure.
     *
     * Returns the number of mappings found.
     *
     *   @param max_cost: maximum allowed mapping cost. Will ignore any mapping greater than this value
     *   @param min_cost: minimum mapping cost. For each node with a cost less than min_cost, 'k' will be increased by 1
     *   @param keep_invalid: Invalid nodes are retained in queue if true
     *   @param erase_tail: Tail end of queue is retained if true
     *
     * What this does NOT do:
     *    -Check if the imported Structure is the same as one in a smaller Supercell
     *
     */
    //*******************************************************************************************
    // This is where the magic happens, part 1
    Index StrucMapper::k_best_maps_better_than(SimpleStructure const &child_struc,
                                               std::set<MappingNode> &queue,
                                               Index k,
                                               double max_cost /*=StrucMapping::big_inf()*/,
                                               double min_cost /*=-TOL*/,
                                               bool keep_invalid /*=false*/,
                                               bool keep_tail /*= false*/,
                                               bool no_partition /*= false*/) const {
      int nfound = 0;
      // Track pairs of supercell volumes that are chemically incompatible
      std::set<std::pair<Index, Index> > vol_mismatch;

      if(k == 0) {
        // If k==0, then we only keep values less than min_cost
        // However, max_cost controls search loop, so we set max_cost to min_cost
        max_cost = min_cost;
      }

      auto it = queue.begin();
      while(it != queue.end() && it->cost <= (max_cost + cost_tol())) {

        auto current = it;
        bool erase = true;
        // If supercell volumes have already been determined incompatible, we do nothing;
        // current node is deleted
        if(!vol_mismatch.count(current->vol_pair())) {
          std::cout << "Before -- nfound: " << nfound << "; k: " << k <<  "; current cost: " << current->cost << std::endl;
          std::cout << "min_cost: " << min_cost << ";  max_cost: " << max_cost << "\n";
          // Consider two exlusive cases (and base case, in which current node isn't even viable)
          if(current->basis_node.empty()) {
            std::cout << "IS LATTICE MAP! Size before: " << queue.size() << "\n";
            // Case 1: Current node only describes a lattice mapping with no basis mapping
            //         Perform basis mapping, insert result into queue, and erase current
            //         (new node must have cost greather than the current node, so will
            //          appear later in the queue)
            if(!Local::initial_basis_maps(child_struc,
                                          *current,
                                          calculator(),
                                          max_cost,
                                          std::inserter(queue, current))) {
              // If no basis maps are viable, it indicates volume mismatch; add to vol_mismatch
              vol_mismatch.insert(current->vol_pair());
            }
            std::cout << "Size after: " << queue.size() << "\n";
          }
          else if(current->is_viable) {

            std::cout << "COMPLETE MAPPING! Size before: " << queue.size() << std::endl;
            // Case 2: Current node is a complete mapping and is viable
            //         Either it is a valid node, and thus part of the solution set,
            //         or it is invalid, but we must add its partition to the queue
            //         because it may have a suboptimal derivative mapping that is valid
            //         and is part of the solution set

            if(current->is_valid && current->cost > min_cost && nfound < k) {
              // current node is part of solution set
              ++nfound;

              // Need to account for case where mapping k+1 is as good as mapping k
              // A few ways to make this happen, but we will do it by increasing k when nfound==k,
              // and shrink max_cost to be *just* higher than current cost
              if(nfound == k) {
                ++k;
                max_cost = current->cost;// + tol();
              }
            }


            // Regardless of validity, we partition current node and add the derivative nodes to queue
            // (but only if we haven't reached stopping condition)
            // Skip partitioning if node is already partitioned, or if caller has asked not to
            if(nfound < k || current->cost <= min_cost) {
              if(!(no_partition || current->is_partitioned)) {
                Local::partition_node(*current,
                                      calculator(),
                                      child_struc,
                                      std::inserter(queue, current));
              }


              // Keep current node if it is in the solution set if we have been asked to keep invalids
              if(current->is_valid || keep_invalid)
                erase = false;
            }
          }
          //Never keep unviable nodes or incomplete nodes
        }

        // Safe to increment here:
        //  1) No continue/break statements
        //  2) Nothing has been deleted yet
        //  3) Everything that needs to be inserted has been inserted
        ++it;

        std::cout << "After -- nfound: " << nfound << "; k: " << k <<  "; next cost: " << it->cost << std::endl;
        std::cout << "min_cost: " << min_cost << ";  max_cost: " << max_cost << "\n";

        // Erase current if no longer needed
        if(erase)
          queue.erase(current);
      }
      // Erase everything worse than the solution set, unless asked not to
      if(!keep_tail && it != queue.end()) {
        queue.erase(it, queue.end());
      }
      return nfound;
    }

    //****************************************************************************************************************

    // Find all Lattice mappings better than min_cost and at most the k best mappings in range [min_cost,max_cost]
    std::set<MappingNode> StrucMapper::_seed_k_best_from_super_lats(SimpleStructure const &child_struc,
                                                                    std::vector<Lattice> const &_parent_scels,
                                                                    std::vector<Lattice> const &_child_scels,
                                                                    Index k,
                                                                    double min_cost /*=1e-6*/,
                                                                    double max_cost /*=StrucMapping::small_inf()*/) const {
      Lattice p_prim_lat(parent().lat_column_mat, xtal_tol());
      Lattice c_prim_lat(child_struc.lat_column_mat, xtal_tol());
      std::set<MappingNode> result;

      if(k == 0)
        max_cost = min_cost;

      for(Lattice const &c_lat : _child_scels) {
        for(Lattice const &p_lat : _parent_scels) {
          int n_child_atom = round(std::abs(volume(c_lat) / volume(c_prim_lat))) * _n_species(child_struc);
          LatticeMap strain_map(p_lat, c_lat, n_child_atom, 1, calculator().point_group(), m_strain_gram_mat, max_cost);
          // strain_map is initialized to first mapping better than 'max_cost', if such a mapping exists
          // We will continue checking possibilities until all such mappings are exhausted
          while(strain_map.strain_cost() < max_cost) {

            // Make k bigger if we find really exceptional mappings, with cost better than min-cost
            if(strain_map.strain_cost() < min_cost
               || (!result.empty() && (result.size() == k && almost_equal(result.rbegin()->cost, strain_map.strain_cost(), 1e-6)))) {
              ++k;
            }

            //auto res =
            result.emplace(LatticeNode(strain_map, p_prim_lat, c_prim_lat), strain_weight());
            if(result.size() > k) {
              result.erase(std::next(result.rbegin()).base());
              if(!result.empty())
                max_cost = max(min_cost, (result.rbegin())->cost);
            }
            strain_map.next_mapping_better_than(max_cost);
          }
        }
      }
      return result;
    }

  }
}
