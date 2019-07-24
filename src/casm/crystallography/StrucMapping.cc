#include "casm/crystallography/StrucMapping.hh"

#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/strain/StrainConverter.hh"
#include "casm/crystallography/Lattice.hh"
#include "casm/crystallography/Coordinate.hh"
#include "casm/crystallography/Niggli.hh"
#include "casm/crystallography/LatticeMap.hh"

namespace CASM {
  //*******************************************************************************************

  namespace StrucMapping {
    double strain_cost(double relaxed_lat_vol, const MappingNode &mapped_result, const Index Nsites) {
      return StrainCostCalculator::iso_strain_cost(mapped_result.lat_node.rstretch, relaxed_lat_vol / double(max(Nsites, Index(1))));
    }

    //*******************************************************************************************

    double basis_cost(const MappingNode &mapped_result, Index Nsites) {
      // mean square displacement distance in deformed coordinate system
      return (mapped_result.lat_node.rstretch * mapped_result.displacement * mapped_result.displacement.transpose() * mapped_result.lat_node.rstretch.transpose()).trace() / double(max(Nsites, Index(1)));
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


    Eigen::Matrix3d F = child_scel.lat_column_mat() * parent_scel.inv_lat_column_mat();
    rstretch = StrainConverter::right_stretch_tensor(F);
    isometry = F * rstretch.inverse();

    if(StrucMapping::is_inf(_cost))
      _cost = StrainCostCalculator::iso_strain_cost(rstretch, child_prim.vol() / double(max(child_N_atom, Index(1))));
    cost = _cost;
  }

  //*******************************************************************************************

  LatticeNode::LatticeNode(LatticeMap const &_lat_map,
                           Lattice const &parent_prim,
                           Lattice const &child_prim,
                           Index child_N_atom) :
    rstretch(StrainConverter::right_stretch_tensor(_lat_map.matrixF())),
    isometry(_lat_map.matrixF() * rstretch.inverse()),
    parent(parent_prim, Lattice(_lat_map.parent_matrix(), parent_prim.tol())),
    child(Lattice(_lat_map.matrixF().inverse() * child_prim.lat_column_mat()), Lattice(_lat_map.parent_matrix(), parent_prim.tol())),
    cost(_lat_map.strain_cost()) {

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
      return MappingNode();
    }

    // tstruc becomes idealized structure
    //SimpleStructure tstruc(child_struc);


    // We know child_struc.lattice() is a supercell of the prim, now we have to
    // reorient 'child_struc' by a point-group operation of the parent to match canonical lattice vectors
    // This may not be a rotation in the child structure's point group
    Lattice derot_c_lat(Lattice(parent().lat_column_mat * trans_mat, m_tol).canonical_form(_calculator().point_group()));

    // We now find a transformation matrix of c_lat so that, after transformation, it is related
    // to derot_c_lat by rigid rotation only. Following line finds R and T such that derot_c_lat = R*c_lat*T
    auto res = is_supercell(derot_c_lat, c_lat, _calculator().point_group().begin(), _calculator().point_group().end(), m_tol);

    std::set<MappingNode> mapping_seed({MappingNode(LatticeNode(Lattice(parent().lat_column_mat, m_tol),
                                                                derot_c_lat,
                                                                c_lat,
                                                                Lattice(child_struc.lat_column_mat * res.second.cast<double>(), m_tol),
                                                                child_struc.n_atom(),
                                                                0. /*strain_cost is zero in ideal case*/),
                                                    m_strain_weight)});


    Index k = k_best_maps_better_than(child_struc, mapping_seed, m_tol, false);
    if(k == 0) {
      return MappingNode();
    }

    return *(mapping_seed.begin());
  }

  //*******************************************************************************************

  std::pair<Index, Index> StrucMapper::_vol_range(const SimpleStructure &child_struc) const {
    Index min_vol(0), max_vol(0);
    //mapped_result.clear();

    if(_calculator().fixed_species().size() > 0) {
      std::string tcompon = _calculator().fixed_species().begin()->first;
      int ncompon(0);
      for(Index i = 0; i < child_struc.n_mol(); i++) {
        if(child_struc.mol_info.names[i] == tcompon)
          ncompon++;
      }
      min_vol = ncompon / int(_calculator().fixed_species().begin()->second);
      max_vol = min_vol;
    }
    else {

      // Try to narrow the range of supercell volumes -- the best bounds are obtained from
      // the convex hull of the end-members, but we need to wait for improvements to convex hull
      // routines

      int max_n_va = _calculator().max_n_va();
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

    SupercellEnumerator<Lattice> enumerator(Lattice(parent().lat_column_mat),
                                            _calculator().point_group(),
                                            ScelEnumProps(prim_vol, prim_vol + 1));

    for(auto it = enumerator.begin(); it != enumerator.end(); ++it) {
      if(m_restricted && !_filter_lat(*it)) {
        continue;
      }
      Lattice canon_lat = *it;
      if(!canon_lat.is_canonical(_calculator().point_group())) {
        canon_lat = canon_lat.canonical_form(_calculator().point_group());
      }
      lat_vec.push_back(canon_lat);
    }

    return lat_vec;
  }

  //***************************************************************************************************

  std::vector<Eigen::Vector3d> SimpleStrucMapCalculator::translations(SimpleStructure const &ichild_struc,
                                                                      LatticeNode const &lat_node) const {

    SimpleStructure::Info const &p_info(info(parent()));
    SimpleStructure::Info const &c_info(info(ichild_struc));

    Index i_trans = 0;

    std::map<std::string, std::vector<Index> > s_map;
    for(std::string const &sp : c_info.names) {
      s_map[sp].push_back(i_trans++);
    }

    auto b_it = s_map.begin();

    i_trans = p_info.names.size();
    if(m_va_allowed.size() < p_info.names.size()) {
      for(Index n = 0; n < p_info.names.size(); ++n) {
        if(!m_va_allowed.count(n)) {
          auto it = s_map.find(p_info.names[n]);

          //It's not possible to map these structures
          if(it == s_map.end())
            return {};

          if(b_it == s_map.end() || (it->second).size() < (b_it->second).size()) {
            b_it = it;
            i_trans = n;
          }
        }
      }
    }
    else {
      throw std::runtime_error(std::string("Mapping child structure onto a parent structure that has no fixed sublattice. Robust handling of this situation has not been implement."));
    }

    std::vector<Eigen::Vector3d> result;
    if(i_trans == p_info.names.size())
      return result;

    result.reserve((b_it->second).size());

    Coordinate translation(lat_node.parent.scel_lattice());
    for(Index n : b_it->second) {
      // Try translating first basis atom onto each chemically compatible site of PRIM
      translation.cart() = p_info.coords.col(i_trans) - c_info.coords.col(n);
      translation.voronoi_within();
      result.push_back(translation.const_cart());
    }

    return result;
  }

  void SimpleStrucMapCalculator::finalize(MappingNode &_node,
                                          SimpleStructure const &child_struc) const {

    populate_displacement(_node, child_struc);
    _node.cost = _node.basis_weight * StrucMapping::basis_cost(_node, info(child_struc).size()) + _node.strain_weight * _node.lat_node.cost;
    return;
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
        if(!insert_at_most_k_maps(child_struc,
                                  *current,
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
        if(_calculator().validate(*current))
          ++k;

        if(k < m_Nbest && !current->is_partitioned)
          partition_node(*current,
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
  //            Assignment Problem methods
  //****************************************************************************************************************

  void SimpleStrucMapCalculator::populate_displacement(MappingNode &_node,
                                                       SimpleStructure const &child_struc) const {

    PrimGrid const &pgrid(_node.lat_node.parent);
    PrimGrid const &cgrid(_node.lat_node.child);
    SimpleStructure::Info const &p_info(info(parent()));
    SimpleStructure::Info const &c_info(info(child_struc));


    _node.permutation = _node.basis_node.permutation();

    // initialize displacement matrix with all zeros
    _node.displacement.setZero(3, pgrid.size()*p_info.size());

    Eigen::Vector3d avg_disp(0, 0, 0);

    Coordinate disp_coord(pgrid.scel_lattice());

    // Populate displacements given as the difference in the Coordinates
    // as described by node.permutation.
    for(Index i = 0; i < _node.permutation.size(); i++) {

      // If we are dealing with a vacancy, its displacment must be zero.
      //if(node.permutation(i) >= child_struc.n_mol()) {
      //  --DO NOTHING--
      //}


      // Using min_dist routine to calculate the displacement vector that corresponds
      // to the distance used in the Cost Matrix and Hungarian Algorithm
      // The method returns the displacement vector pointing from the
      // IDEAL coordinate to the RELAXED coordinate
      if(_node.permutation[i] < c_info.size()*cgrid.size()) {

        Coordinate relaxed_coord(c_info.coords.col(_node.permutation[i] % cgrid.size())
                                 + cgrid.scel_coord(_node.permutation[i] / cgrid.size()).const_cart()
                                 + _node.basis_node.translation,
                                 pgrid.scel_lattice(), CART);

        Coordinate ideal_coord = pgrid.scel_coord(i / pgrid.size());
        ideal_coord.cart() += p_info.coords.col(i % cgrid.size());

        relaxed_coord.min_dist(ideal_coord, disp_coord);
        _node.disp(i) = disp_coord.const_cart();

        avg_disp += _node.disp(i);
      }
    }

    avg_disp /= max<double>(double(c_info.size() * cgrid.size()), 1.);

    _node.displacement.colwise() -= avg_disp;
    // End of filling displacements
  }

  /*
   * Finding the cost_matrix given the relaxed structure
   * This will always return a square matrix with the extra elements
   * reflecting the vacancies specified in the ideal supercell.
   * Costs are calculated in context of the lattice.
   * cost_matrix(i,j) is cost of mapping child site 'j' onto parent site 'i'
   */
  //****************************************************************************************************************
  bool SimpleStrucMapCalculator::populate_cost_mat(MappingNode &_node,
                                                   SimpleStructure const &child_struc) const {
    PrimGrid const &pgrid(_node.lat_node.parent);
    PrimGrid const &cgrid(_node.lat_node.child);
    Eigen::Vector3d const &translation(_node.basis_node.translation);
    Eigen::MatrixXd &cost_matrix(_node.basis_node.cost_mat);
    Eigen::Matrix3d metric = _node.lat_node.rstretch.transpose() * _node.lat_node.rstretch;

    SimpleStructure::Info const &p_info(info(parent()));
    SimpleStructure::Info const &c_info(info(child_struc));

    Index pN = p_info.size() * pgrid.size();
    Index cN = c_info.size() * cgrid.size();
    if(pN != cN)
      return false;

    //if(cost_matrix.rows()!=scel.num_sites() || cost_matrix.cols()!=scel.num_sites())
    cost_matrix = Eigen::MatrixXd::Constant(pN, pN, StrucMapping::small_inf());
    Index inf_counter;
    // loop through all the sites of the structure

    Index j = 0;
    Index l = 0;
    for(; j < c_info.size(); j++) {
      for(Index n = 0; n < cgrid.size(); ++n, ++l) {
        Coordinate current_relaxed_coord(c_info.coords.col(j) + cgrid.scel_coord(n).const_cart() + translation, pgrid.scel_lattice(), CART);
        // loop through all the sites in the supercell
        inf_counter = 0;
        Index k = 0;
        for(Index i = 0; i < p_info.size(); ++i) {
          if(!_allowed_species()[i].count(c_info.names[j])) {
            k += pgrid.size();
            ++inf_counter;
            continue;
          }
          Coordinate curr_parent_coord(p_info.coords.col(i), pgrid.scel_lattice(), CART);

          for(Index m = 0; m < pgrid.size(); ++m, ++k) {
            // Check if relaxed atom j is allowed on site i
            // If so, populate cost_matrix normally
            cost_matrix(k, l) = (curr_parent_coord + pgrid.scel_coord(m)).min_dist2(current_relaxed_coord, metric);
          }
        }
        if(inf_counter == p_info.size()) {
          //std:: cerr << "Bail at 1\n";
          return false;
        }
      }
    }

    // If there are unvisited columns of cost_mat (because fewer sites in the child structure than in parent),
    // we will treat them as a vacant species and set them to zero cost for mapping onto parent sites that can
    // host vacancies
    for(; j < p_info.size(); j++) {
      if(m_va_allowed.empty())
        return false;

      for(Index n = 0; n < cgrid.size(); ++n, ++l) {
        inf_counter = 0;

        Index k = 0;
        for(Index i = 0; i < p_info.size(); ++i) {
          if(!m_va_allowed.count(i)) {
            k += pgrid.size();
            ++inf_counter;
            continue;
          }

          for(Index m = 0; m < pgrid.size(); ++m, ++k) {
            // Check if relaxed atom j is allowed on site i
            // If so, populate cost_matrix normally
            cost_matrix(k, l) = 0.;
          }
        }
      }
      if(inf_counter == p_info.size()) {
        //std:: cerr << "Bail at 2\n";
        return false;
      }
    }

    // JCT: I'm not sure if there's an easy way to check if the cost matrix is viable in all cases
    //      Some of the simpler checks I could think of failed for edge cases with vacancies.
    //      If we return an invalid cost matrix, the Hungarian routines will detect that it is invalid,
    //      so maybe there's no point in doing additional checks here.
    return true;
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
