#include "casm/crystallography/StrucMapping.hh"

#include "casm/crystallography/SuperlatticeEnumerator.hh"
#include "casm/crystallography/StrucMapCalculatorInterface.hh"
#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/strain/StrainConverter.hh"
#include "casm/crystallography/Lattice.hh"
#include "casm/crystallography/Lattice_impl.hh"
#include "casm/crystallography/Coordinate.hh"
#include "casm/crystallography/CanonicalForm.hh"
#include "casm/crystallography/LatticeMap.hh"

namespace CASM {
  namespace xtal {
    //*******************************************************************************************

    namespace StrucMapping {
      double strain_cost(double relaxed_lat_vol, const MappingNode &mapped_result, const Index Nsites) {
        return StrainCostCalculator::iso_strain_cost(mapped_result.lat_node.stretch, relaxed_lat_vol / double(max(Nsites, Index(1))));
      }

      //*******************************************************************************************

      double basis_cost(const MappingNode &mapped_result, Index Nsites) {
        // mean square displacement distance in deformed coordinate system
        return (mapped_result.lat_node.stretch.inverse() * mapped_result.displacement).squaredNorm() / double(max(Nsites, Index(1)));
      }
    }

    //*******************************************************************************************
    namespace Local {
      // Local helper function for StrucMapper::k_best_maps_better_than
      template<typename OutputIterator>
      static bool insert_at_most_k_maps(SimpleStructure child_struc,
                                        MappingNode const &seed,
                                        StrucMapCalculatorInterface const &calculator,
                                        double max_cost,
                                        Index k,
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
          //node.cost = strain_weight() * node.lat_node.cost + basis_weight() * (StrucMapping::basis_cost(node, c_info.size() * node.lat_node.child.size()) + m_tol * node.basis_node.translation.norm() / 10.0);
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
                                 SimpleStructure const &_child_struc,
                                 OutputIterator it) {
        Index j, jj, currj, tj;
        Index n = _node.basis_node.assignment.size();
        MappingNode t1(_node),
                    t2(_node.lat_node, _node.strain_weight);

        MappingNode *p1 = &t1;
        MappingNode *p2 = &t2;
        for(Index m = 0; m < (n - 1); ++m) {
          HungarianNode &n1(p1->basis_node);
          HungarianNode &n2(p2->basis_node);
          //clear assignment and cost_mat
          n2.assignment.clear();
          n2.cost_mat.resize(n - m - 1, n - m - 1);

          // We are forcing on the first site assignment in t1: i.e.,  [0,currj]
          // this involves striking row 0 and col currj from t2's cost_mat
          // and augmenting the cost_offset of t2 by the cost of [0,currj]
          currj = n1.assignment[0];
          n2.cost_offset = n1.cost_offset + n1.cost_mat(0, currj);

          // [0,currj] have local context only. We store the forced assignment in forced_on
          // using the original indexing
          n2.forced_on.emplace(n2.irow[0], n2.icol[currj]);

          // Strike row 0 and col currj to form new HungarianNode for t2
          // (i,j) indexes the starting cost_mat, (i-1, jj) indexes the resulting cost_mat
          n2.irow = std::vector<Index>(++n1.irow.begin(), n1.irow.end());
          for(j = 0, jj = 0; j < (n - m); ++j, ++jj) {
            if(j == currj) {
              --jj;
              continue;
            }
            n2.icol.push_back(n1.icol[j]);

            // We will also store an updated assignment vector in t2, which will be
            // used to construct next node of partition
            tj = n1.assignment[j];
            if(tj > currj)
              tj--;
            n2.assignment.push_back(tj);

            // Fill col jj of t2's cost mat
            for(Index i = 1; i < n; ++i)
              n2.cost_mat(i - 1, jj) = n1.cost_mat(i, j);
          }
          //t2 properly initialized; we can now force OFF [0,currj] in t1, and add it to node list
          n1.cost_mat(0, currj) = StrucMapping::big_inf();
          n1.assignment.clear();
          p1->calc();
          _calculator.finalize(*p1, _child_struc);
          it = *p1;
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
      child(Lattice((parent_scel.lat_column_mat() * child_scel.inv_lat_column_mat())
                    * child_prim.lat_column_mat()),
            parent_scel) {

      // F is from ideal parent to child
      Eigen::Matrix3d F = child_scel.lat_column_mat() * parent_scel.inv_lat_column_mat();

      // stretch is from (de-rotated, strained) child to ideal parent
      stretch = StrainConverter::right_stretch_tensor(F).inverse();

      // isometry is from child to strained parent
      isometry = (F * stretch).transpose();

      if(StrucMapping::is_inf(_cost))
        _cost = StrainCostCalculator::iso_strain_cost(stretch, child_prim.vol() / double(max(child_N_atom, Index(1))));
      cost = _cost;
    }

    //*******************************************************************************************

    LatticeNode::LatticeNode(LatticeMap const &_lat_map,
                             Lattice const &parent_prim,
                             Lattice const &child_prim,
                             Index child_N_atom) :
      // stretch is from (de-rotated, strained) child to ideal parent
      stretch(StrainConverter::right_stretch_tensor(_lat_map.matrixF()).inverse()),

      // isometry is from child to strained parent
      isometry((_lat_map.matrixF() * stretch).transpose()),
      parent(parent_prim, Lattice(_lat_map.parent_matrix(), parent_prim.tol())),
      child(Lattice(_lat_map.matrixF().inverse() * child_prim.lat_column_mat()), Lattice(_lat_map.parent_matrix(), parent_prim.tol())),
      cost(_lat_map.strain_cost()) {

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
        basis_node.cost = hungarian_method(basis_node.cost_mat, basis_node.assignment, tol()) + basis_node.cost_offset;
        if(StrucMapping::is_inf(basis_node.cost)) {
          is_viable = false;
          cost = StrucMapping::big_inf();
        }
        else {
          cost += basis_weight + basis_node.cost;
        }
      }
      else
        cost = StrucMapping::big_inf();
    }

    //*******************************************************************************************

    StrucMapper::StrucMapper(StrucMapCalculatorInterface const &calculator,
                             double _strain_weight /*= 0.5*/,
                             Index _Nbest/*= 1*/,
                             double _max_volume_change /*= 0.5*/,
                             int _options /*= robust*/, // this should actually be a bitwise-OR of StrucMapper::Options
                             double _tol /*= TOL*/,
                             double _min_va_frac /*= 0.*/,
                             double _max_va_frac /*= 1.*/) :
      m_calc_ptr(calculator.clone()),
      //squeeze strain_weight into (0,1] if necessary
      m_strain_weight(max(min(_strain_weight, 1.0), 1e-9)),
      m_Nbest(_Nbest),
      m_max_volume_change(_max_volume_change),
      m_options(_options),
      m_tol(max(1e-9, _tol)),
      m_min_va_frac(0.),
      m_max_va_frac(1.) {

      //ParamComposition param_comp(_pclex.prim());
      m_max_volume_change = max(m_tol, _max_volume_change);
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

    MappingNode StrucMapper::map_ideal_struc(const SimpleStructure &child_struc) const {

      // Lattice::is_supercell_of() isn't very smart right now, and will return
      // false if the two lattices differ by a rigid rotation
      // In the future this may not be the case, so we will assume that child_struc may
      // be rigidly rotated relative to prim
      Eigen::Matrix3d trans_mat;

      // c_lat must be an ideal supercell of the parent lattice, but it need not be canonical
      // We will account for the difference in orientation between c_lat and the canonical supercell,
      // which must be related by a point group operation
      Lattice c_lat(child_struc.lat_column_mat, m_tol);


      if(!c_lat.is_supercell_of(Lattice(parent().lat_column_mat), trans_mat)) {
        /*std::cerr << "CRITICAL ERROR: In map_ideal_struc(), primitive structure does not tile the provided\n"
          << "                superstructure. Please use map_deformed_struc() instead.\n"
          << "                Exiting...\n";
        */
        return MappingNode::invalid();
      }

      // tstruc becomes idealized structure
      //SimpleStructure tstruc(child_struc);


      // We know child_struc.lattice() is a supercell of the prim, now we have to
      // reorient 'child_struc' by a point-group operation of the parent to match canonical lattice vectors
      // This may not be a rotation in the child structure's point group
      Lattice derot_c_lat(canonical::equivalent(Lattice(parent().lat_column_mat * trans_mat, m_tol), calculator().point_group()));

      // We now find a transformation matrix of c_lat so that, after transformation, it is related
      // to derot_c_lat by rigid rotation only. Following line finds R and T such that derot_c_lat = R*c_lat*T
      auto res = is_supercell(derot_c_lat, c_lat, calculator().point_group().begin(), calculator().point_group().end(), m_tol);

      std::set<MappingNode> mapping_seed({MappingNode(LatticeNode(Lattice(parent().lat_column_mat, m_tol),
                                                                  derot_c_lat,
                                                                  c_lat,
                                                                  Lattice(child_struc.lat_column_mat * res.second.cast<double>(), m_tol),
                                                                  child_struc.n_atom(),
                                                                  0. /*strain_cost is zero in ideal case*/),
                                                      m_strain_weight)
                                         });


      Index k = k_best_maps_better_than(child_struc, mapping_seed, m_tol, false);
      if(k == 0) {
        return MappingNode::invalid();
      }

      return *(mapping_seed.begin());
    }

    //*******************************************************************************************

    std::pair<Index, Index> StrucMapper::_vol_range(const SimpleStructure &child_struc) const {
      Index min_vol(0), max_vol(0);
      //mapped_result.clear();

      if(calculator().fixed_species().size() > 0) {
        std::string tcompon = calculator().fixed_species().begin()->first;
        int ncompon(0);
        for(Index i = 0; i < child_struc.n_mol(); i++) {
          if(child_struc.mol_info.names[i] == tcompon)
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
        double N_sites = double(parent().n_mol());
        double max_va_frac_limit = double(max_n_va) / N_sites;
        double t_min_va_frac = min(min_va_frac(), max_va_frac_limit);
        double t_max_va_frac = min(max_va_frac(), max_va_frac_limit);
        // min_vol assumes min number vacancies -- best case scenario
        min_vol = ceil(child_struc.n_atom() / (N_sites * (1. - t_min_va_frac)) - m_tol);

        // This is for the worst case scenario -- lots of vacancies
        max_vol = ceil(child_struc.n_atom() / (N_sites * (1. - t_max_va_frac)) - m_tol);

        if(t_max_va_frac > TOL) {
          //Nvol is rounded integer volume-- assume that answer is within 30% of this volume, and use it to tighten our bounds
          int Nvol = round(std::abs(child_struc.lat_column_mat.determinant() / parent().lat_column_mat.determinant()));
          min_vol = min(max_vol, max<Index>(round((1.0 - m_max_volume_change) * double(Nvol)), min_vol));
          max_vol = max(min_vol, min<Index>(round((1.0 + m_max_volume_change) * double(Nvol)), max_vol));
        }
      }

      min_vol = max<Index>(min_vol, 1);
      max_vol = max<Index>(max_vol, 1);
      return std::pair<Index, Index>(min_vol, max_vol);

    }

    //*******************************************************************************************

    SimpleStructure const &StrucMapper::parent() const {
      return calculator().parent();
    }


    //*******************************************************************************************

    std::set<MappingNode> StrucMapper::seed_from_vol_range(SimpleStructure const &child_struc,
                                                           Index min_vol,
                                                           Index max_vol) const {
      int Nkeep = 10 + 5 * m_Nbest;
      if(!valid_index(min_vol) || !valid_index(min_vol) || max_vol < min_vol) {
        auto vol_range = _vol_range(child_struc);
        min_vol = vol_range.first;
        max_vol = vol_range.second;
      }

      std::set<MappingNode> mapping_seed;
      for(Index i_vol = min_vol; i_vol <= max_vol; i_vol++) {
        std::vector<Lattice> lat_vec;
        lat_vec = _lattices_of_vol(i_vol);

        std::set<MappingNode> t_seed = seed_k_best_from_super_lats(child_struc,
                                                                   lat_vec,
        {Lattice(child_struc.lat_column_mat)},
        Nkeep);

        mapping_seed.insert(std::make_move_iterator(t_seed.begin()), std::make_move_iterator(t_seed.end()));

      }
      return mapping_seed;
    }

    //*******************************************************************************************

    std::set<MappingNode> StrucMapper::map_deformed_struc(const SimpleStructure &child_struc,
                                                          double best_cost /*=1e20*/,
                                                          bool keep_invalid) const {
      auto vols = _vol_range(child_struc);
      std::set<MappingNode> mapping_seed = seed_from_vol_range(child_struc, vols.first, vols.second);
      k_best_maps_better_than(child_struc, mapping_seed, best_cost, keep_invalid);
      return mapping_seed;
    }

    //*******************************************************************************************

    std::set<MappingNode> StrucMapper::map_deformed_struc_impose_lattice(const SimpleStructure &child_struc,
                                                                         const Lattice &imposed_lat,
                                                                         double best_cost,
                                                                         bool keep_invalid) const {

      std::set<MappingNode> mapping_seed = seed_k_best_from_super_lats(child_struc,
      {imposed_lat},
      {Lattice(child_struc.lat_column_mat)},
      m_Nbest);

      k_best_maps_better_than(child_struc, mapping_seed, best_cost, keep_invalid);
      return mapping_seed;
    }

    //*******************************************************************************************

    std::set<MappingNode> StrucMapper::map_deformed_struc_impose_lattice_node(const SimpleStructure &child_struc,
                                                                              const LatticeNode &imposed_node,
                                                                              double best_cost,
                                                                              bool keep_invalid) const {

      std::set<MappingNode> mapping_seed;
      mapping_seed.emplace(imposed_node, m_strain_weight);
      k_best_maps_better_than(child_struc, mapping_seed, best_cost, keep_invalid);
      return mapping_seed;
    }

    //*******************************************************************************************

    std::vector<Lattice> StrucMapper::_lattices_of_vol(Index prim_vol) const {
      if(!m_restricted) {
        //If you specified that you wanted certain lattices, return those, otherwise do the
        //usual enumeration
        if(this->lattices_constrained()) {
          //This may very well return an empty vector, saving painful time enumerating things
          return m_allowed_superlat_map[prim_vol];
        }

        if(!valid_index(prim_vol)) {
          throw std::runtime_error("Cannot enumerate lattice of volume " + std::to_string(prim_vol) + ", which is out of bounds.\n");
        }

        //If we already have candidate lattices for the given volume, return those
        auto it = m_superlat_map.find(prim_vol);
        if(it != m_superlat_map.end())
          return it->second;
      }

      //We don't have any lattices for the provided volume, enumerate them all!!!
      std::vector<Lattice> tlat_vec;
      std::vector<Lattice> &lat_vec = (m_restricted ? tlat_vec : m_superlat_map[prim_vol]);

      auto pg = calculator().point_group();
      SuperlatticeEnumerator enumerator(pg.begin(), pg.end(), Lattice(parent().lat_column_mat),
                                        ScelEnumProps(prim_vol, prim_vol + 1));

      for(auto it = enumerator.begin(); it != enumerator.end(); ++it) {
        if(m_restricted && !_filter_lat(*it)) {
          continue;
        }
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
     * Given a structure and a set of mapping nodes, iterate over many supercells of the prim and for ideal supercells that
     * nearly match the structure's lattice, try to map the structure's basis onto that supercell
     *
     * This process iteratively improves the mapping cost function -> total_cost = w*strain_cost + (1-w)*basis_cost
     * where 'w' is the lattice-cost weight parameter. It's unlikely that there is an objective way to choose 'w'
     * without information regarding the physics of the system (e.g., where is the extremum of the energy barrier
     * when going from a particular ideal configuration to the specified deformed structure)
     *
     * Returns false if no mapping is possible (which should only happen if the passed structure has an incompatible
     * chemical composition).
     *
     * What this does NOT do:
     *    -Check if the imported Structure is the same as one in a smaller Supercell
     *
     */
    //*******************************************************************************************
    // This is where the magic happens, part 1
    Index StrucMapper::k_best_maps_better_than(SimpleStructure const &child_struc,
                                               std::set<MappingNode> &queue,
                                               double max_cost,
                                               bool keep_invalid,
                                               bool erase_tail /*= true*/) const {
      int k = 0;

      std::set<std::pair<Index, Index> > vol_mismatch;

      auto it = queue.begin();
      while(it != queue.end() && k < m_Nbest) {
        auto current = it;
        ++it;

        if(vol_mismatch.find(current->vol_pair()) == vol_mismatch.end()) {
          queue.erase(current);
        }
        else if(current->basis_node.empty()) {
          if(!Local::insert_at_most_k_maps(child_struc,
                                           *current,
                                           calculator(),
                                           max_cost,
                                           m_Nbest,
                                           std::inserter(queue, current))) {
            vol_mismatch.insert(current->vol_pair());
          }
          queue.erase(current);
        }
        else if(current->cost > max_cost) {
          return k;
        }
        else if(current->is_viable) {
          if(current->is_valid)
            ++k;

          if(k < m_Nbest && !current->is_partitioned)
            Local::partition_node(*current,
                                  calculator(),
                                  child_struc,
                                  std::inserter(queue, current));

          if(!current->is_valid && ! keep_invalid)
            queue.erase(current);
        }
      }

      if(erase_tail && it != queue.end()) {
        queue.erase(it, queue.end());
      }
      return k;
    }

    //****************************************************************************************************************

    // Find all Lattice mappings better than min_cost and at most the k best mappings in range [min_cost,max_cost]
    std::set<MappingNode> StrucMapper::seed_k_best_from_super_lats(SimpleStructure const &child_struc,
                                                                   std::vector<Lattice> const &_parent_scels,
                                                                   std::vector<Lattice> const &_child_scels,
                                                                   Index k,
                                                                   double min_cost /*=1e-6*/,
                                                                   double max_cost /*=StrucMapping::small_inf()*/) const {
      Lattice p_lat(parent().lat_column_mat);
      Lattice c_lat(child_struc.lat_column_mat);
      std::set<MappingNode> result;

      for(Lattice const &c_lat : _child_scels) {
        for(Lattice const &p_lat : _parent_scels) {
          LatticeMap strain_map(p_lat, c_lat, child_struc.n_atom(), tol(), 1, calculator().point_group(), m_strain_gram_mat, max_cost);

          // strain_map is initialized to first mapping better than 'max_cost', if such a mapping exists
          // We will continue checking possibilities until all such mappings are exhausted
          while(strain_map.strain_cost() < max_cost) {

            // Make k bigger if we find really exception mappings
            if(strain_map.strain_cost() < min_cost)
              ++k;

            result.emplace(LatticeNode(strain_map, p_lat, c_lat, child_struc.n_atom()), strain_weight());
            if(result.size() > k) {
              result.erase(std::next(result.rbegin()).base());
              max_cost = (result.rbegin())->cost;
            }
            strain_map.next_mapping_better_than(max_cost);
          }
        }
      }
      return result;
    }

  }
}
