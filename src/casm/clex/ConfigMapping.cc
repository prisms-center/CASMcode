#include "casm/clex/PrimClex.hh"
#include "casm/clex/ConfigMapping.hh"
#include "casm/strain/StrainConverter.hh"
#include "casm/crystallography/Lattice.hh"
#include "casm/crystallography/Niggli.hh"
#include "casm/crystallography/LatticeMap.hh"
#include "casm/crystallography/SupercellEnumerator.hh"

namespace CASM {
  namespace ConfigMapping {
    double strain_cost(const Lattice &relaxed_lat, const ConfigDoF &_dof, const Index Nsites) {
      return LatticeMap::calc_strain_cost(_dof.deformation(), relaxed_lat.vol() / double(max(Nsites, Index(1))));
    }

    //*******************************************************************************************

    double basis_cost(const ConfigDoF &_dof, Index Nsites) {
      // mean square displacement distance in deformed coordinate system
      return (_dof.deformation() * _dof.displacement() * _dof.displacement().transpose() * _dof.deformation().transpose()).trace() / double(max(Nsites, Index(1)));
    }

    //*******************************************************************************************
    Lattice find_nearest_super_lattice(const Lattice &prim_lat,
                                       const Lattice &relaxed_lat,
                                       const SymGroup &sym_group,
                                       Eigen::Matrix3d &deformation,
                                       Eigen::Matrix3d &trans_mat,
                                       const std::vector<Lattice> &from_range,
                                       double _tol) {
      Lattice best_lat;

      double best_cost = 10e10;

      //std::cout << "min_vol is " << min_vol << "max_vol is " << max_vol << "\n";
      for(auto it = from_range.cbegin(); it != from_range.cend(); ++it) {
        //Supercell matches size, now check to see see if lattices are related to each other
        LatticeMap strainmap(*it, relaxed_lat, 1, _tol, 1);
        if(strainmap.best_strain_mapping().strain_cost() < best_cost) {
          trans_mat = strainmap.matrixN();
          deformation = strainmap.matrixF();
          best_cost = strainmap.strain_cost();
          best_lat = *it;
        }
      }
      //std::cout << "Best strain mapping is " << best_cost << " with F = \n" << deformation << "\n";

      return best_lat;
    }

    //*******************************************************************************************
    Lattice find_nearest_super_lattice(const Lattice &prim_lat,
                                       const Lattice &relaxed_lat,
                                       const SymGroup &sym_group,
                                       Eigen::Matrix3d &deformation,
                                       Eigen::Matrix3d &trans_mat,
                                       Index min_vol,
                                       Index max_vol,
                                       double _tol) {
      Lattice best_lat;
      if(!valid_index(max_vol))
        max_vol = ceil(std::abs(relaxed_lat.vol()) / std::abs(prim_lat.vol()));

      double best_cost = 10e10;
      //We only bother checking pre-existing supercells of min_vol <= volume <=max_vol;
      SupercellEnumerator<Lattice> enumerator(prim_lat, sym_group, ScelEnumProps(min_vol, max_vol + 1));

      //std::cout << "min_vol is " << min_vol << "max_vol is " << max_vol << "\n";
      for(auto it = enumerator.begin(); it != enumerator.end(); ++it) {
        //std::cout << "Enumeration step " << l++ << " best cost is " << best_cost << "\n";
        Lattice tlat = canonical_equivalent_lattice(*it, sym_group, _tol);

        //Supercell matches size, now check to see see if lattices are related to each other
        LatticeMap strainmap(tlat, relaxed_lat, 1, _tol, 1);
        if(strainmap.best_strain_mapping().strain_cost() < best_cost) {
          trans_mat = strainmap.matrixN();
          deformation = strainmap.matrixF();
          best_cost = strainmap.strain_cost();
          best_lat = tlat;
        }
      }
      //std::cout << "Best strain mapping is " << best_cost << " with F = \n" << deformation << "\n";

      return best_lat;
    }


  }

  //*******************************************************************************************

  ConfigMapper::ConfigMapper(PrimClex &_pclex,
                             double _lattice_weight,
                             double _max_volume_change/*=0.5*/,
                             int options/*=robust*/,
                             double _tol/*=TOL*/) :
    m_pclex(&_pclex),
    m_lattice_weight(_lattice_weight),
    m_max_volume_change(_max_volume_change),
    m_min_va_frac(0.),
    m_max_va_frac(1.),
    m_robust_flag(options & robust),
    m_strict_flag(options & strict),
    m_rotate_flag(options & rotate),
    m_tol(max(1e-9, _tol)) {
    //squeeze lattice_weight into (0,1] if necessary
    m_lattice_weight = max(min(_lattice_weight, 1.0), 1e-9);
    ParamComposition param_comp(_pclex.get_prim());
    m_fixed_components = param_comp.fixed_components();
    m_max_volume_change = max(m_tol, _max_volume_change);
  }

  //*******************************************************************************************

  bool ConfigMapper::import_structure_occupation(const fs::path &pos_path,
                                                 std::string &imported_name,
                                                 jsonParser &relaxation_properties,
                                                 std::vector<Index> &best_assignment,
                                                 Eigen::Matrix3d &cart_op) const {

    try {
      BasicStructure<Site> tstruc(pos_path);
      return import_structure_occupation(tstruc,
                                         imported_name,
                                         relaxation_properties,
                                         best_assignment,
                                         cart_op);
    }
    catch(const std::exception &ex) {
      throw std::runtime_error(std::string("Could not successfully import structure ") + pos_path.string() + ":\n" + ex.what());
    }


  }

  //*******************************************************************************************

  bool ConfigMapper::import_structure_occupation(BasicStructure<Site> &_struc,
                                                 std::string &imported_name,
                                                 jsonParser &relaxation_properties,
                                                 std::vector<Index> &best_assignment,
                                                 Eigen::Matrix3d &cart_op,
                                                 bool update_struc) const {
    return import_structure_occupation(_struc,
                                       nullptr,
                                       imported_name,
                                       relaxation_properties,
                                       best_assignment,
                                       cart_op,
                                       update_struc);
  }

  //*******************************************************************************************

  bool ConfigMapper::import_structure_occupation(BasicStructure<Site> &_struc,
                                                 const Configuration *hint_ptr,
                                                 std::string &imported_name,
                                                 jsonParser &relaxation_properties,
                                                 std::vector<Index> &best_assignment,
                                                 Eigen::Matrix3d &cart_op,
                                                 bool update_struc) const {

    //Indices for Configuration index and permutation operation index
    ConfigDoF tconfigdof, suggested_configdof;
    Lattice mapped_lat;
    bool is_new_config(true);
    double bc(1e20), sc(1e20), hint_cost = 1e20;

    relaxation_properties.put_obj();
    //std::vector<Index> best_assignment;

    if(hint_ptr != nullptr) {
      if(ConfigMap_impl::struc_to_configdof(*hint_ptr,
                                            _struc,
                                            suggested_configdof,
                                            best_assignment,
                                            m_robust_flag, // translate_flag -- not sure what to use for this
                                            m_tol)) {
        mapped_lat = (hint_ptr->get_supercell()).get_real_super_lattice();
        bc = ConfigMapping::basis_cost(suggested_configdof, _struc.basis.size());
        sc = ConfigMapping::strain_cost(_struc.lattice(), suggested_configdof, _struc.basis.size());
        relaxation_properties["suggested_mapping"]["basis_deformation"] = bc;
        relaxation_properties["suggested_mapping"]["lattice_deformation"] = sc;
        relaxation_properties["suggested_mapping"]["volume_relaxation"] = suggested_configdof.deformation().determinant();
        hint_cost = m_lattice_weight * sc + (1.0 - m_lattice_weight) * bc - m_tol;
      }
      else {
        //relaxation_properties["suggested_mapping"] = "unknown";
      }
    }
    if(struc_to_configdof(_struc,
                          tconfigdof,
                          mapped_lat,
                          best_assignment,
                          cart_op,
                          hint_cost)) {

      bc = ConfigMapping::basis_cost(tconfigdof, _struc.basis.size());
      sc = ConfigMapping::strain_cost(_struc.lattice(), tconfigdof, _struc.basis.size());
      // robust_cost = m_lattice_weight * sc + (1.0 - m_lattice_weight) * bc - m_tol;
      relaxation_properties["best_mapping"]["basis_deformation"] = bc;
      relaxation_properties["best_mapping"]["lattice_deformation"] = sc;
      relaxation_properties["best_mapping"]["volume_relaxation"] = tconfigdof.deformation().determinant();

    }
    else {
      if(hint_cost > 1e10)
        throw std::runtime_error("Structure is incompatible with PRIM.");

      swap(tconfigdof, suggested_configdof);
      relaxation_properties["best_mapping"] = relaxation_properties["suggested_mapping"];
    }


    ConfigDoF relaxed_occ;

    relaxed_occ.set_occupation(tconfigdof.occupation());
    Supercell::permute_const_iterator it_canon;

    if(hint_ptr != nullptr) {
      Supercell &scel(hint_ptr->get_supercell());
      if(mapped_lat.is_equivalent(scel.get_real_super_lattice(), m_tol)) {
        if(m_strict_flag && relaxed_occ.occupation() == (hint_ptr->configdof()).occupation()) {
          // config is unchanged
          imported_name = hint_ptr->name();
          is_new_config = false;
          it_canon = hint_ptr->get_supercell().permute_begin();
        }
        else {

          Supercell::permute_const_iterator relaxed_it_canon = Configuration(scel, jsonParser(), relaxed_occ).to_canonical();
          Supercell::permute_const_iterator ideal_rev_it_canon = hint_ptr->from_canonical();
          it_canon = ideal_rev_it_canon * relaxed_it_canon;

          if(relaxed_occ.occupation() == copy_apply(it_canon.inverse(), *hint_ptr).occupation()) {
            // config is unchanged
            imported_name = hint_ptr->name();
            is_new_config = false;
          }
        }
      }
    }

    if(is_new_config) {
      Index import_scel_index = primclex().add_supercell(mapped_lat), import_config_index;

      Configuration import_config(primclex().get_supercell(import_scel_index), jsonParser(), relaxed_occ);

      if(m_strict_flag) {
        it_canon = primclex().get_supercell(import_scel_index).permute_begin();
        is_new_config = primclex().get_supercell(import_scel_index).add_canon_config(import_config, import_config_index);
        imported_name = primclex().get_supercell(import_scel_index).get_config(import_config_index).name();
      }
      else {
        is_new_config = primclex().get_supercell(import_scel_index).add_config(import_config, import_config_index, it_canon);
        imported_name = primclex().get_supercell(import_scel_index).get_config(import_config_index).name();
      }
    }
    else {

    }
    // transform deformation tensor to match canonical form and apply operation to cart_op
    ConfigDoF trans_configdof = copy_apply(it_canon, tconfigdof);
    relaxation_properties["best_mapping"]["relaxation_deformation"] = trans_configdof.deformation();
    relaxation_properties["best_mapping"]["relaxation_displacement"] = trans_configdof.displacement().transpose();

    cart_op = it_canon.sym_op().matrix() * cart_op;

    // compose permutations
    std::vector<Index>tperm = it_canon.combined_permute().permute(best_assignment);

    //copy non-vacancy part of permutation into best_assignment
    best_assignment.resize(_struc.basis.size());
    Index num_atoms = _struc.basis.size();
    std::copy_if(tperm.cbegin(), tperm.cend(),
                 best_assignment.begin(),
    [num_atoms](Index i) {
      return i < num_atoms;
    });

    if(update_struc) {
      _struc.set_lattice(Lattice(cart_op.transpose()*tconfigdof.deformation()*mapped_lat.lat_column_mat()), CART);
      _struc.set_lattice(Lattice(tconfigdof.deformation()*mapped_lat.lat_column_mat()), FRAC);
    }
    return is_new_config;
  }

  //*******************************************************************************************

  bool ConfigMapper::import_structure(const fs::path &pos_path,
                                      std::string &imported_name,
                                      jsonParser &relaxation_properties,
                                      std::vector<Index> &best_assignment,
                                      Eigen::Matrix3d &cart_op) const {

    try {
      return import_structure(BasicStructure<Site>(pos_path),
                              imported_name,
                              relaxation_properties,
                              best_assignment,
                              cart_op);
    }
    catch(const std::exception &ex) {
      throw std::runtime_error(std::string("Could not successfully import structure ") + pos_path.string() + ":\n" + ex.what());
    }


  }

  //*******************************************************************************************
  bool ConfigMapper::import_structure(const BasicStructure<Site> &_struc,
                                      std::string &imported_name,
                                      jsonParser &relaxation_properties,
                                      std::vector<Index> &best_assignment,
                                      Eigen::Matrix3d &cart_op) const {

    //Indices for Configuration index and permutation operation index
    Supercell::permute_const_iterator it_canon;

    ConfigDoF tconfigdof;
    Lattice mapped_lat;
    bool new_config_flag;
    //std::vector<Index> best_assignment;
    if(!struc_to_configdof(_struc,
                           tconfigdof,
                           mapped_lat,
                           best_assignment,
                           cart_op))
      throw std::runtime_error("Structure is incompatible with PRIM.");

    relaxation_properties["best_mapping"]["basis_deformation"] = ConfigMapping::basis_cost(tconfigdof, _struc.basis.size());
    relaxation_properties["best_mapping"]["lattice_deformation"] = ConfigMapping::strain_cost(_struc.lattice(), tconfigdof, _struc.basis.size());
    relaxation_properties["best_mapping"]["volume_change"] = tconfigdof.deformation().determinant();

    Index import_scel_index = primclex().add_supercell(mapped_lat), import_config_index;

    Configuration import_config(primclex().get_supercell(import_scel_index), jsonParser(), tconfigdof);

    if(m_strict_flag) {
      it_canon = primclex().get_supercell(import_scel_index).permute_begin();
      new_config_flag = primclex().get_supercell(import_scel_index).add_canon_config(import_config, import_config_index);
      imported_name = primclex().get_supercell(import_scel_index).get_config(import_config_index).name();
    }
    else {
      it_canon = primclex().get_supercell(import_scel_index).permute_begin();
      new_config_flag = primclex().get_supercell(import_scel_index).add_config(import_config, import_config_index, it_canon);
      imported_name = primclex().get_supercell(import_scel_index).get_config(import_config_index).name();
    }

    relaxation_properties["best_mapping"]["relaxation_deformation"] = it_canon.sym_op().matrix() * tconfigdof.deformation() * it_canon.sym_op().matrix().transpose();

    cart_op = it_canon.sym_op().matrix() * cart_op;

    // compose permutations
    std::vector<Index>tperm = it_canon.combined_permute().permute(best_assignment);

    //copy non-vacancy part of permutation into best_assignment
    best_assignment.resize(_struc.basis.size());
    Index num_atoms = _struc.basis.size();
    std::copy_if(tperm.cbegin(), tperm.cend(),
                 best_assignment.begin(),
    [num_atoms](Index i) {
      return i < num_atoms;
    });
    return new_config_flag;

  }

  //*******************************************************************************************
  bool ConfigMapper::struc_to_configdof(const BasicStructure<Site> &struc,
                                        ConfigDoF &mapped_configdof,
                                        Lattice &mapped_lat) const {
    std::vector<Index> t_assign;
    Eigen::Matrix3d t_op;
    return struc_to_configdof(struc, mapped_configdof, mapped_lat, t_assign, t_op);
  }
  //*******************************************************************************************
  bool ConfigMapper::struc_to_configdof(const BasicStructure<Site> &struc,
                                        ConfigDoF &mapped_configdof,
                                        Lattice &mapped_lat,
                                        std::vector<Index> &best_assignment,
                                        Eigen::Matrix3d &cart_op,
                                        double best_cost /*=1e20*/) const {

    bool valid_mapping(false);
    // If structure's lattice is a supercell of the primitive lattice, then import as ideal_structure
    if(!m_robust_flag && struc.lattice().is_supercell_of(primclex().get_prim().lattice(), m_tol)) {
      valid_mapping = ideal_struc_to_configdof(struc,
                                               mapped_configdof,
                                               mapped_lat,
                                               best_assignment,
                                               cart_op);
      //std::cout << "valid_mapping is " << valid_mapping << "\n";
      valid_mapping = valid_mapping && ConfigMapping::basis_cost(mapped_configdof, struc.basis.size()) < (10 * m_tol);
    }

    // If structure's lattice is not a supercell of the primitive lattice, then import as deformed_structure
    if(!valid_mapping) { // if not a supercell or m_robust_flag=true, treat as deformed
      valid_mapping = deformed_struc_to_configdof(struc,
                                                  mapped_configdof,
                                                  mapped_lat,
                                                  best_assignment,
                                                  cart_op,
                                                  best_cost);
    }
    return valid_mapping;
  }

  //*******************************************************************************************
  /*
   * Given a structure and a primclex, find a perfect supercell of the prim that is equivalent to structure's lattice
   * and then try to map the structure's basis onto that supercell
   *
   * Returns false if no mapping is possible, or if the lattice is not ideal
   *
   * What this does NOT do:
   *    -Check if the imported Structure is the same as one in a smaller Supercell
   *
   */
  //*******************************************************************************************
  bool ConfigMapper::ideal_struc_to_configdof(const BasicStructure<Site> &struc,
                                              ConfigDoF &mapped_configdof,
                                              Lattice &mapped_lat,
                                              std::vector<Index> &best_assignment,
                                              Eigen::Matrix3d &cart_op) const {
    // Lattice::is_supercell_of() isn't very smart right now, and will return false if the two lattices differ by a rigid rotation
    // In the future this may not be the case, so we will assume that struc may be rigidly rotated relative to prim
    Eigen::Matrix3d trans_mat;
    if(!struc.lattice().is_supercell_of(primclex().get_prim().lattice(), trans_mat,  m_tol)) {
      /*std::cerr << "CRITICAL ERROR: In ideal_struc_to_configdof(), primitive structure does not tile the provided\n"
        << "                superstructure. Please use deformed_struc_to_configdof() instead.\n"
        << "                Exiting...\n";
      */
      return false;
    }
    BasicStructure<Site> tstruc(struc);

    // We know struc.lattice() is a supercell of the prim, now we have to reorient 'struc' to match canonical lattice vectors
    mapped_lat = canonical_equivalent_lattice(Lattice(primclex().get_prim().lattice().lat_column_mat() * trans_mat), primclex().get_prim().point_group(), m_tol);
    Supercell scel(&primclex(), mapped_lat);

    // note: trans_mat gets recycled here
    mapped_lat.is_supercell_of(struc.lattice(), primclex().get_prim().point_group(), trans_mat,  m_tol);
    tstruc.set_lattice(Lattice(tstruc.lattice().lat_column_mat()*trans_mat), CART);

    //cart_op goes from imported coordinate system to ideal (PRIM) coordinate system
    cart_op = mapped_lat.lat_column_mat() * tstruc.lattice().inv_lat_column_mat();
    tstruc.set_lattice(mapped_lat, FRAC);
    if(!m_rotate_flag) {
      cart_op = Eigen::Matrix3d::Identity(3, 3);
    }
    return ConfigMap_impl::preconditioned_struc_to_configdof(scel,
                                                             tstruc,
                                                             cart_op.transpose(),//cart_op.transpose() is deformation in this case
                                                             mapped_configdof,
                                                             best_assignment,
                                                             true,
                                                             m_tol);
  }

  //*******************************************************************************************
  /*
   * Given a structure and a primclex, iterate over many supercells of the prim and for ideal supercells that
   * nearly match the structure's lattice, try to map the structure's basis onto that supercell
   *
   * This process iteratively improves the mapping cost function -> total_cost = w*lattice_cost + (1-w)*basis_cost
   * where 'w' is the lattice-cost weight parameter. It's unlikely that there is an objective way to choose 'w'
   * without information regarding the physics of the system (e.g., where is the extremum of the energy barrier
   * when going from a particular ideal configuration to the specified deformed structure)
   *
   * Returns false if no mapping is possible (which should only happen if the passed structure has a composition
   * that cannot be realized in the system described by the primclex).
   *
   * What this does NOT do:
   *    -Check if the imported Structure is the same as one in a smaller Supercell
   *
   */
  //*******************************************************************************************
  bool ConfigMapper::deformed_struc_to_configdof(const BasicStructure<Site> &struc,
                                                 ConfigDoF &mapped_configdof,
                                                 Lattice &mapped_lat,
                                                 std::vector<Index> &best_assignment,
                                                 Eigen::Matrix3d &cart_op,
                                                 double best_cost /*=1e20*/) const {

    //squeeze lattice_weight into [0,1] if necessary
    double lw = m_lattice_weight;
    double bw = 1.0 - lw;


    std::vector<Index> assignment;
    //Add new Supercell if it doesn't exist already. Use primitive point group to check for equivalence and
    //store transformation matrix
    Eigen::Matrix3d deformation;
    double num_atoms = double(struc.basis.size());
    int min_vol, max_vol;

    mapped_configdof.clear();
    if(m_fixed_components.size() > 0) {
      std::string tcompon = m_fixed_components[0].first;
      int ncompon(0);
      for(Index i = 0; i < struc.basis.size(); i++) {
        if(struc.basis[i].occ_name() == tcompon)
          ncompon++;
      }
      min_vol = ncompon / int(m_fixed_components[0].second);
      max_vol = min_vol;
    }
    else {
      // Try to narrow the range of supercell volumes -- the best bounds are obtained from
      // the convex hull of the end-members, but we need to wait for improvements to convex hull
      // routines

      int max_n_va = primclex().get_prim().max_possible_vacancies();
      double max_va_frac_limit = double(max_n_va) / double(primclex().get_prim().basis.size());
      double t_min_va_frac = min(min_va_frac(), max_va_frac_limit);
      double t_max_va_frac = min(max_va_frac(), max_va_frac_limit);

      // min_vol assumes min number vacancies -- best case scenario
      min_vol = ceil((num_atoms / (double(primclex().get_prim().basis.size())) * 1. - t_min_va_frac) - m_tol);

      // This is for the worst case scenario -- lots of vacancies
      max_vol = ceil(num_atoms / (double(primclex().get_prim().basis.size()) * (1.0 - t_max_va_frac)) - m_tol);

      if(t_max_va_frac > TOL) {
        //Nvol is rounded integer volume-- assume that answer is within 30% of this volume, and use it to tighten our bounds
        int Nvol = round(std::abs(struc.lattice().vol() / primclex().get_prim().lattice().vol()));
        int new_min_vol = min(max_vol, max(round((1.0 - m_max_volume_change) * double(Nvol)), min_vol));
        int new_max_vol = max(min_vol, min(round((1.0 + m_max_volume_change) * double(Nvol)), max_vol));
        max_vol = new_max_vol;
        min_vol = new_min_vol;
      }
    }

    min_vol = max(min_vol, 1);
    max_vol = max(max_vol, 1);


    Eigen::Matrix3d ttrans_mat, tF, rotF;

    double strain_cost(1e10), basis_cost(1e10), tot_cost;//, best_strain_cost, best_basis_cost;
    ConfigDoF tdof;
    BasicStructure<Site> tstruc(struc);
    Lattice tlat;

    // First pass:  Find a reasonable upper bound
    for(Index i_vol = min_vol; i_vol <= max_vol; i_vol++) {
      tlat = ConfigMapping::find_nearest_super_lattice(primclex().get_prim().lattice(),
                                                       struc.lattice(),
                                                       primclex().get_prim().point_group(),
                                                       tF,
                                                       ttrans_mat,
                                                       _lattices_of_vol(i_vol),
                                                       m_tol);
      strain_cost = lw * LatticeMap::calc_strain_cost(tF, struc.lattice().vol() / max(num_atoms, 1.));

      if(best_cost < strain_cost) {
        continue;
      }

      tstruc = struc;

      //make tstruc an un-rotated, un-strained version of struc
      tstruc.set_lattice(Lattice(tlat.lat_column_mat()*ttrans_mat), FRAC);
      tstruc.set_lattice(tlat, CART);
      rotF = tF;
      if(m_rotate_flag) {
        rotF = StrainConverter::right_stretch_tensor(tF);
      }

      Supercell scel(&primclex(), tlat);
      if(!ConfigMap_impl::preconditioned_struc_to_configdof(scel,
                                                            tstruc,
                                                            rotF,
                                                            tdof,
                                                            assignment,
                                                            true,
                                                            m_tol))
        continue;
      basis_cost = bw * ConfigMapping::basis_cost(tdof, struc.basis.size());
      tot_cost = strain_cost + basis_cost;

      if(tot_cost < best_cost) {
        best_cost = tot_cost - m_tol;
        swap(best_assignment, assignment);

        cart_op = rotF * tF.inverse();

        swap(mapped_configdof, tdof);

        mapped_lat = tlat;
      }
    }


    //Second pass: Find the absolute best mapping
    for(Index i_vol = min_vol; i_vol <= max_vol; i_vol++) {
      //std::cout << "  vol = " << i_vol << "\n";
      const std::vector<Lattice> &lat_vec = _lattices_of_vol(i_vol);

      for(auto it = lat_vec.cbegin(); it != lat_vec.cend(); ++it) {
        if(!deformed_struc_to_configdof_of_lattice(struc,
                                                   *it,
                                                   best_cost,
                                                   mapped_configdof,
                                                   mapped_lat,
                                                   best_assignment,
                                                   cart_op))
          break;
      }
    }

    //std::cout << "best_cost = " << best_cost << "    best_strain_cost = " << best_strain_cost << "    best_basis_cost = " << best_basis_cost << "\n";
    if(best_cost > 1e9) {
      //std::cerr << "WARNING: In Supercell::import_deformed_structure(), no successful mapping was found for Structure " << src << "\n"
      //        << "         This Structure may be incompatible with the ideal crystal specified by the PRIM file. \n";
      return false;
    }

    // If mapped_configdof is empty, it means that nothing better than best_cost was found
    return mapped_configdof.size() > 0;
  }

  //*******************************************************************************************
  bool ConfigMapper::deformed_struc_to_configdof_of_lattice(const BasicStructure<Site> &struc,
                                                            const Lattice &imposed_lat,
                                                            double &best_cost,
                                                            ConfigDoF &mapped_configdof,
                                                            Lattice &mapped_lat,
                                                            std::vector<Index> &best_assignment,
                                                            Eigen::Matrix3d &cart_op) const {
    double strain_cost, basis_cost, tot_cost;
    ConfigDoF tdof;
    BasicStructure<Site> tstruc(struc);
    Eigen::Matrix3d tF, rotF;
    std::vector<Index> assignment;
    Supercell scel(&primclex(), imposed_lat);

    double lw = m_lattice_weight;
    double bw = 1.0 - lw;
    double num_atoms = double(struc.basis.size());
    //Determine best mapping for this supercell

    //Initialize with simplest mapping onto supercell 'i', so that we don't change the crystal setting unnecessarily
    tF = struc.lattice().lat_column_mat() * imposed_lat.inv_lat_column_mat();

    strain_cost = lw * LatticeMap::calc_strain_cost(tF, struc.lattice().vol() / max(num_atoms, 1.));

    // If simplest mapping seems viable, check it further
    if(strain_cost < best_cost) {
      tstruc = struc;
      tstruc.set_lattice(imposed_lat, FRAC);
      rotF = tF;
      if(m_rotate_flag) {
        rotF = StrainConverter::right_stretch_tensor(tF);
      }

      if(!ConfigMap_impl::preconditioned_struc_to_configdof(scel,
                                                            tstruc,
                                                            rotF,
                                                            tdof,
                                                            assignment,
                                                            true,
                                                            m_tol))
        return false;

      basis_cost = bw * ConfigMapping::basis_cost(tdof, struc.basis.size());
      //std::cout << "\n**Starting strain_cost = " << strain_cost << ";   and basis_cost = " << basis_cost << "  TOTAL: " << strain_cost + basis_cost << "\n";
      tot_cost = strain_cost + basis_cost;

      //std::cout << "    scel.name = " << scel.get_name()  << "  simple map: strain_cost " << strain_cost << "   best_cost " << best_cost << "    tot_cost " << tot_cost << "\n";
      if(tot_cost < best_cost) {
        best_cost = tot_cost - m_tol;
        //best_strain_cost = strain_cost;
        //best_basis_cost = basis_cost;
        swap(best_assignment, assignment);
        cart_op = rotF * tF.inverse();
        //best_trans = Matrix3<double>::identity();
        //std::cout << "tF is:\n" << tF << "\n and N is:\n" << best_trans << "\n";
        swap(mapped_configdof, tdof);
        mapped_lat = imposed_lat;
      }
    } // Done checking simplest mapping

    // If the simplest mapping is best, we have avoided a lot of extra work, but we still need to check for
    // non-trivial mappings that are better than both the simplest mapping and the best found mapping
    LatticeMap strainmap(imposed_lat, struc.lattice(), round(num_atoms), m_tol, 1);
    strain_cost = lw * strainmap.strain_cost();
    if(best_cost < strain_cost)
      strain_cost = lw * strainmap.next_mapping_better_than(best_cost).strain_cost();

    while(strain_cost < best_cost) {  // only enter loop if there's a chance of improving on current best

      tstruc = struc;
      // We modify the deformed structure so that its lattice is a deformed version of the nearest ideal lattice
      // Don't need matrixN if we use set_lattice(CART), because matrixF depends on matrixN implicitly
      tstruc.set_lattice(Lattice(imposed_lat.lat_column_mat()*strainmap.matrixN()), FRAC);
      tstruc.set_lattice(imposed_lat, CART);
      tF = strainmap.matrixF();
      rotF = tF;
      if(m_rotate_flag) {
        rotF = StrainConverter::right_stretch_tensor(tF);
      }

      if(!ConfigMap_impl::preconditioned_struc_to_configdof(scel,
                                                            tstruc,
                                                            rotF,
                                                            tdof,
                                                            assignment,
                                                            true,
                                                            m_tol)) {
        return false;
        //no longer unexpected
        //throw std::runtime_error("Unexpected error in deformed_struc_to_config_dof(). This should never happen!\n");
      }
      basis_cost = bw * ConfigMapping::basis_cost(tdof, struc.basis.size());
      //std::cout << "New strain_cost = " << strain_cost << ";   and basis_cost = " << basis_cost << "  TOTAL: " << strain_cost + basis_cost << "\n";
      //std::cout << "  Compare -> best_cost = " << best_cost << "\n";

      tot_cost = strain_cost + basis_cost;
      //std::cout << "      complex map: best_cost " << best_cost << "    strain_cost " << strain_cost << "    tot_cost " << tot_cost << "\n";

      if(tot_cost < best_cost) {
        //std::cout << "Old best_trans, with cost " << best_cost << ":\n" << best_trans << "\n";
        best_cost = tot_cost - m_tol;
        swap(best_assignment, assignment);
        //best_strain_cost = strain_cost;
        //best_basis_cost = basis_cost;
        cart_op = rotF * tF.inverse();
        //best_trans = strainmap.matrixN();
        //std::cout << "New best_trans, with cost " << best_cost << ":\n" << best_trans << "\n";
        swap(mapped_configdof, tdof);
        mapped_lat = imposed_lat;
      }
      // This finds first decomposition:
      //            struc.lattice() =  deformation*(supercell_list[import_scel_index].get_real_super_lattice())*equiv_mat
      //   that has cost function less than best_cost
      strain_cost = lw * strainmap.next_mapping_better_than(best_cost).strain_cost();
    }
    return true;
  }

  //*******************************************************************************************
  const std::vector<Lattice> &ConfigMapper::_lattices_of_vol(Index prim_vol) const {
    if(!valid_index(prim_vol)) {
      throw std::runtime_error("Cannot enumerate lattice of volume " + std::to_string(prim_vol) + ", which is out of bounds.\n");
    }
    auto it = m_superlat_map.find(prim_vol);
    if(it != m_superlat_map.end())
      return it->second;

    std::vector<Lattice> lat_vec;
    SupercellEnumerator<Lattice> enumerator(
      primclex().get_prim().lattice(),
      primclex().get_prim().point_group(),
      ScelEnumProps(prim_vol, prim_vol + 1));

    //std::cout << "min_vol is " << min_vol << "max_vol is " << max_vol << "\n";
    for(auto it = enumerator.begin(); it != enumerator.end(); ++it) {
      //std::cout << "Enumeration step " << l++ << " best cost is " << best_cost << "\n";
      lat_vec.push_back(canonical_equivalent_lattice(*it, primclex().get_prim().point_group(), m_tol));
    }

    return m_superlat_map[prim_vol] = std::move(lat_vec);

  }

  //****************************************************************************************************************

  namespace ConfigMap_impl {

    //****************************************************************************************************************
    //            Assignment Problem methods
    //****************************************************************************************************************
    /*
     * Finding the cost_matrix given the relaxed structure
     * This will always return a square matrix with the extra elements
     * reflecting the vacancies specified in the ideal supercell.
     * Costs are calculated in context of the real_super_lattice.
     */
    //****************************************************************************************************************

    bool calc_cost_matrix(const Supercell &scel,
                          const BasicStructure<Site> &rstruc,
                          const Coordinate &trans,
                          const Eigen::Matrix3d &metric,
                          Eigen::MatrixXd &cost_matrix) {

      if(rstruc.basis.size() > scel.num_sites())
        return false;
      double inf = 10E10;
      //if(cost_matrix.rows()!=scel.num_sites() || cost_matrix.cols()!=scel.num_sites())
      cost_matrix = Eigen::MatrixXd::Constant(scel.num_sites(), scel.num_sites(), inf);
      Index inf_counter;
      // loop through all the sites of the structure
      Index j = 0;
      for(; j < rstruc.basis.size(); j++) {
        Coordinate current_relaxed_coord(rstruc.basis[j].frac(), scel.get_real_super_lattice(), FRAC);
        current_relaxed_coord.cart() += trans.cart();
        // loop through all the sites in the supercell
        inf_counter = 0;
        for(Index i = 0; i < scel.num_sites(); i++) {

          // Check if relaxed atom j is allowed on site i
          // If so, populate cost_matrix normally
          if(scel.get_prim().basis[scel.get_b(i)].contains(rstruc.basis[j].occ_name())) {
            cost_matrix(i, j) = scel.coord(i).min_dist2(current_relaxed_coord, metric);
          }
          // If not, set cost_matrix (i,j) = inf
          else {
            cost_matrix(i, j) = inf;
            inf_counter++;
          }
        }
        if(inf_counter == scel.num_sites()) {
          //std:: cerr << "Bail at 1\n";
          return false;
        }
      }


      for(; j < scel.num_sites(); j++) {
        inf_counter = 0;
        for(Index i = 0; i < scel.num_sites(); i++) {

          // Check if vacancies are allowed at each position in the supercell
          if(scel.get_prim().basis[scel.get_b(i)].contains("Va")) {
            cost_matrix(i, j) = 0;
          }
          else {
            cost_matrix(i, j) = inf;
            inf_counter++;
          }
        }
        if(inf_counter == scel.num_sites()) {
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
    bool calc_cost_matrix(const Configuration &config,
                          const BasicStructure<Site> &rstruc,
                          const Coordinate &trans,
                          const Eigen::Matrix3d &metric,
                          Eigen::MatrixXd &cost_matrix) {


      double inf = 10E10;
      const Supercell &scel(config.get_supercell());
      //if(cost_matrix.rows()!=scel.num_sites() || cost_matrix.cols()!=scel.num_sites())
      cost_matrix = Eigen::MatrixXd::Constant(scel.num_sites(), scel.num_sites(), inf);
      Index inf_counter;
      // loop through all the sites of the structure
      Index j;
      for(j = 0; j < rstruc.basis.size(); j++) {
        Coordinate current_relaxed_coord(rstruc.basis[j].frac(), scel.get_real_super_lattice(), FRAC);
        current_relaxed_coord.cart() += trans.cart();
        // loop through all the sites in the supercell
        inf_counter = 0;
        for(Index i = 0; i < scel.num_sites(); i++) {

          // Check if relaxed atom j is allowed on site i
          // If so, populate cost_matrix normally
          if(config.get_mol(i).name == rstruc.basis[j].occ_name()) {
            cost_matrix(i, j) = scel.coord(i).min_dist2(current_relaxed_coord, metric);
          }
          // If not, set cost_matrix (i,j) = inf
          else {
            cost_matrix(i, j) = inf;
            inf_counter++;
          }
        }
        if(inf_counter == scel.num_sites()) {
          //std:: cerr << "Bail at 1\n";
          return false;
        }
      }


      for(; j < scel.num_sites(); j++) {
        inf_counter = 0;
        for(Index i = 0; i < scel.num_sites(); i++) {

          // Check if vacancies are allowed at each position in the supercell
          if(config.get_mol(i).name == "Va") {
            cost_matrix(i, j) = 0;
          }
          else {
            cost_matrix(i, j) = inf;
            inf_counter++;
          }
        }
        if(inf_counter == scel.num_sites()) {
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

    //***************************************************************************************************
    // Options:
    //
    // translate_flag = true means that the rigid-translations are removed. (typically this option should be used)
    //
    // translate_flag = false means that rigid translations are not considered. (probably don't want to use this since structures should be considered equal if they are related by a rigid translation).
    //
    bool struc_to_configdof(const Supercell &scel,
                            BasicStructure<Site> rstruc,
                            ConfigDoF &config_dof,
                            std::vector<Index> &best_assignments,
                            const bool translate_flag,
                            const double _tol) {
      Eigen::Matrix3d deformation = rstruc.lattice().lat_column_mat() * scel.get_real_super_lattice().inv_lat_column_mat();
      // un-deform rstruc
      rstruc.set_lattice(scel.get_real_super_lattice(), FRAC);
      return preconditioned_struc_to_configdof(scel, rstruc, deformation, config_dof, best_assignments, translate_flag, _tol);
    }

    //***************************************************************************************************
    // Options:
    //
    // translate_flag = true means that the rigid-translations are removed. (typically this option should be used)
    //
    // translate_flag = false means that rigid translations are not considered. (probably don't want to use this since structures should be considered equal if they are related by a rigid translation).
    //
    bool struc_to_configdof(const Configuration &config,
                            BasicStructure<Site> rstruc,
                            ConfigDoF &config_dof,
                            std::vector<Index> &best_assignments,
                            const bool translate_flag,
                            const double _tol) {
      const Lattice &mapped_lat(config.get_supercell().get_real_super_lattice());
      Eigen::Matrix3d deformation = rstruc.lattice().lat_column_mat() * mapped_lat.inv_lat_column_mat();
      // un-deform rstruc
      rstruc.set_lattice(mapped_lat, FRAC);
      return preconditioned_struc_to_configdof(config, rstruc, deformation, config_dof, best_assignments, translate_flag, _tol);
    }

    //***************************************************************************************************
    // Options:
    //
    // translate_flag = true means that the rigid-translations are removed. (typically this option should be used)
    //
    // translate_flag = false means that rigid translations are not considered. (probably don't want to use this since structures should be considered equal if they are related by a rigid translation).
    //
    bool preconditioned_struc_to_configdof(const Supercell &scel,
                                           const BasicStructure<Site> &rstruc,
                                           const Eigen::Matrix3d &deformation,
                                           ConfigDoF &config_dof,
                                           std::vector<Index> &best_assignments,
                                           const bool translate_flag,
                                           const double _tol) {
      //std::cout << "CONVERTING TO CONFIGDOF:\n";
      //rstruc.print(std::cout);
      //std::cout << std::endl;
      // clear config_dof and set its deformation
      config_dof.clear();

      config_dof.set_deformation(deformation);
      Eigen::Matrix3d metric(deformation.transpose()*deformation);
      //Initialize everything

      Eigen::MatrixXd cost_matrix;
      std::vector<Index> optimal_assignments;
      //BasicStructure<Site> best_ideal_struc(rstruc);
      Coordinate best_trans(rstruc.lattice());
      double min_mean = 10E10;

      // We want to get rid of translations.
      // define translation such that:
      //    IDEAL = RELAXED + translation
      // and use it when calculating cost matrix

      Index num_translations(1);

      if(rstruc.basis.size())
        num_translations += scel.get_prim().basis.size();

      //std::cout << "num_translations is " << num_translations << "\n";
      for(Index n = 0; n < num_translations; n++) {
        double mean;

        //shift_struc has **ideal lattice**
        //BasicStructure<Site> shift_struc(rstruc);


        if(n > 0 && !scel.get_prim().basis[n - 1].contains(rstruc.basis[0].occ_name()))
          continue;

        Coordinate translation(scel.get_prim().lattice());

        // Always try the non-translated case (n==0), in case it gives best result
        // Also try translating first basis atom onto each chemically compatible site of PRIM (n>0)
        if(n > 0) {
          translation.cart() = scel.coord((n - 1) * scel.volume()).const_cart() - rstruc.basis[0].const_cart();
          translation.voronoi_within();
        }

        if(!ConfigMap_impl::calc_cost_matrix(scel, rstruc, translation, metric, cost_matrix)) {
          /// Indicates that structure is incompatible with supercell, so return false
          return false;
        }

        // The mapping routine is called here
        mean = hungarian_method(cost_matrix, optimal_assignments, _tol);

        // if optimal_assignments is smaller than rstruc.basis.size(), then rstruc is incompattible with supercell
        // (optimal_assignments.size()==0 if the hungarian routine detects an incompatibility)
        if(optimal_assignments.size() < rstruc.basis.size()) {
          return false;
        }


        // add small penalty (~_tol) for larger translation distances, so that shortest equivalent translation is used
        mean += _tol * translation.const_cart().norm() / 10.0;

        if(mean < min_mean) {
          //std::cout << "mean " << mean << " is better than min_mean " << min_mean <<"\n";

          // new best
          swap(best_assignments, optimal_assignments);

          // best shifted structure
          best_trans.cart() = translation.cart();

          // update the minimum mean costs
          min_mean = mean;
        }
      }

      // Now we are filling up displacements
      //
      // Make zero_vector for special vacancy cases.
      Eigen::Vector3d zero_vector(0.0, 0.0, 0.0);

      // initialize displacement matrix with all zeros
      config_dof.set_displacement(ConfigDoF::displacement_matrix_t::Zero(3, scel.num_sites()));

      Eigen::Vector3d avg_disp(0, 0, 0);

      Coordinate disp_coord(rstruc.lattice());

      // Populate displacements given as the difference in the Coordinates
      // as described by best_assignments.
      for(Index i = 0; i < best_assignments.size(); i++) {

        // If we are dealing with a vacancy, its displacment must be zero.
        //if(best_assignments(i) >= rstruc.basis.size()) {
        //  --DO NOTHING--
        //}


        // Using min_dist routine to calculate the displacement vector that corresponds
        // to the distance used in the Cost Matrix and Hungarian Algorithm
        // The method returns the displacement vector pointing from the
        // IDEAL coordinate to the RELAXED coordinate
        if(best_assignments[i] < rstruc.basis.size()) {

          Coordinate ideal_coord(scel.coord(i).frac(), rstruc.lattice(), FRAC);

          (rstruc.basis[best_assignments[i]] + best_trans).min_dist(ideal_coord, disp_coord);
          //std::cout << "min_dist" << std::endl;
          //std::cout << rstruc.basis[optimal_assignments(i)].min_dist(pos_coord, disp_coord);
          config_dof.disp(i) = disp_coord.const_cart();

          avg_disp += config_dof.disp(i);
        }
      }

      avg_disp /= max(double(rstruc.basis.size()), 1.);


      // End of filling displacements


      // Make the assignment bitstring
      //
      // Loop through all supercell sites
      config_dof.set_occupation(Array<int>(scel.num_sites()));
      std::string rel_basis_atom;
      for(Index i = 0; i < best_assignments.size(); i++) {
        // subtract off average displacement
        if(best_assignments[i] < rstruc.basis.size()) {
          config_dof.disp(i) -= avg_disp;
          // suppress small ugly numbers.
          for(Index j = 0; j < 3; j++) {
            if(almost_zero(config_dof.disp(i)[j], 1e-8))
              config_dof.disp(i)[j] = 0;
          }

          //Record basis atom
          rel_basis_atom = rstruc.basis[best_assignments[i]].occ_name();
        }
        else {
          // Any value of the assignment vector larger than the number
          // of sites in the relaxed structure is by construction
          // specified as a vacancy.
          rel_basis_atom = "Va";
        }

        // set occupant and check for errors
        if(!scel.get_prim().basis[scel.get_b(i)].contains(rel_basis_atom, config_dof.occ(i))) {
          //std::cout << "best_assignments is " << best_assignments << "\n";
          //std::cout << "at site " << i << " corresponding to basis " << scel.get_b(i) << "  attempting to assign type " << rel_basis_atom << "\n";
          //std::cout << "Cost Matrix is \n" << cost_matrix << "\n";
          //std::cerr << "CRITICAL ERROR: In Supercell::struc_to_configdof atoms of relaxed/custom structure are incompatible\n"
          //        << "                with the number or type of atomic species allowed in PRIM. Exiting...\n";
          //exit(1);

          return false;
        }
      }

      return true;
    }



    //***************************************************************************************************
    // Options:
    //
    // translate_flag = true means that the rigid-translations are removed. (typically this option should be used)
    //
    // translate_flag = false means that rigid translations are not considered. (probably don't want to use this since structures should be considered equal if they are related by a rigid translation).
    //
    bool preconditioned_struc_to_configdof(const Configuration &config,
                                           const BasicStructure<Site> &rstruc,
                                           const Eigen::Matrix3d &deformation,
                                           ConfigDoF &config_dof,
                                           std::vector<Index> &best_assignments,
                                           const bool translate_flag,
                                           const double _tol) {

      const Supercell &scel(config.get_supercell());

      // clear config_dof and set its deformation
      config_dof.clear();

      config_dof.set_deformation(deformation);
      Eigen::Matrix3d metric(deformation.transpose()*deformation);
      //Initialize everything

      Eigen::MatrixXd cost_matrix;
      std::vector<Index> optimal_assignments;
      //BasicStructure<Site> best_ideal_struc(rstruc);
      Coordinate best_trans(rstruc.lattice());

      double min_mean = 10E10;

      // We want to get rid of translations.
      // define translation such that:
      //    IDEAL = RELAXED + translation
      // and use it when calculating cost matrix

      Index num_translations(1);

      num_translations += rstruc.basis.size();

      //num_translations = rstruc.basis.size();
      //std::cout << "num_translations is " << num_translations << "\n";
      for(Index n = 0; n < num_translations; n++) {
        double mean;

        if(n > 0 && config.get_mol(0).name != rstruc.basis[n - 1].occ_name())
          continue;

        Coordinate translation(scel.get_real_super_lattice());

        // Always try the non-translated case (n==0), in case it gives best result
        // Also try translating first basis atom onto each chemically compatible site of PRIM (n>0)
        if(n > 0) {
          translation.cart() = scel.coord(0).const_cart() - rstruc.basis[n - 1].const_cart();
          translation.voronoi_within();
        }

        if(!ConfigMap_impl::calc_cost_matrix(config, rstruc, translation, metric, cost_matrix)) {
          //std::cerr << "In Supercell::struc_to_config. Cannot construct cost matrix." << std::endl;
          //std::cerr << "This message is probably OK, if you are using translate_flag == true." << std::endl;
          //continue;
          return false;
        }

        // The mapping routine is called here
        mean = hungarian_method(cost_matrix, optimal_assignments, _tol);

        // if optimal_assignments is smaller than rstruc.basis.size(), then rstruc is incompattible
        // with the supercell (optimal_assignments.size()==0 if the hungarian routine detects an incompatibility)
        if(optimal_assignments.size() < rstruc.basis.size()) {
          return false;
        }

        //std::cout << "mean is " << mean << " and stddev is " << stddev << "\n";
        // add small penalty (~_tol) for larger translation distances, so that shortest equivalent translation is used
        mean += _tol * translation.const_cart().norm() / 10.0;
        if(mean < min_mean) {
          //std::cout << "mean " << mean << " is better than min_mean " << min_mean <<"\n";

          // new best
          swap(best_assignments, optimal_assignments);

          // best shifted structure
          best_trans = translation;

          // update the minimum mean costs
          min_mean = mean;
        }
      }

      // Now we are filling up displacements
      //
      // Make zero_vector for special vacancy cases.
      Eigen::Vector3d zero_vector(0.0, 0.0, 0.0);

      // initialize displacement matrix with all zeros
      config_dof.set_displacement(ConfigDoF::displacement_matrix_t::Zero(3, scel.num_sites()));

      Eigen::Vector3d avg_disp(0, 0, 0);

      Coordinate disp_coord(rstruc.lattice());

      // Populate displacements given as the difference in the Coordinates
      // as described by best_assignments.
      for(Index i = 0; i < best_assignments.size(); i++) {

        // If we are dealing with a vacancy, its displacment must be zero.
        //if(best_assignments(i) >= rstruc.basis.size()) {
        //  --DO NOTHING--
        //}


        // Using min_dist routine to calculate the displacement vector that corresponds
        // to the distance used in the Cost Matrix and Hungarian Algorithm
        // The method returns the displacement vector pointing from the
        // IDEAL coordinate to the RELAXED coordinate
        if(best_assignments[i] < rstruc.basis.size()) {

          Coordinate ideal_coord(scel.coord(i).frac(), rstruc.lattice(), FRAC);

          (rstruc.basis[best_assignments[i]] + best_trans).min_dist(ideal_coord, disp_coord);
          //std::cout << "min_dist" << std::endl;
          //std::cout << rstruc.basis[optimal_assignments(i)].min_dist(pos_coord, disp_coord);
          config_dof.disp(i) = disp_coord.const_cart();

          avg_disp += config_dof.disp(i);
        }
      }

      avg_disp /= max(double(rstruc.basis.size()), 1.);


      // End of filling displacements


      // Make the assignment bitstring
      //
      // Loop through all supercell sites
      config_dof.set_occupation(Array<int>(scel.num_sites()));
      std::string rel_basis_atom;
      for(Index i = 0; i < best_assignments.size(); i++) {
        // subtract off average displacement (non-vacant sites only)
        // Any value of the assignment vector larger than the number
        // of sites in the relaxed structure is by construction
        // specified as a vacancy.
        if(best_assignments[i] < rstruc.basis.size()) {
          config_dof.disp(i) -= avg_disp;
          // suppress small ugly numbers.
          for(Index j = 0; j < 3; j++) {
            if(almost_zero(config_dof.disp(i)[j], 1e-8))
              config_dof.disp(i)[j] = 0;
          }
        }
        config_dof.occ(i) = config.occ(i);
      }

      return true;
    }

  }
}
