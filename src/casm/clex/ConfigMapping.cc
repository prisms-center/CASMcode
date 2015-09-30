#include "casm/clex/PrimClex.hh"
#include "casm/clex/ConfigMapping.hh"
#include "casm/strain/StrainConverter.hh"
#include "casm/crystallography/Lattice.hh"
#include "casm/crystallography/LatticeMap.hh"
#include "casm/crystallography/SupercellEnumerator.hh"

namespace CASM {
  //*******************************************************************************************
  namespace ConfigMapping {
    double strain_cost(const Lattice &relaxed_lat, const ConfigDoF &_dof) {
      return LatticeMap::calc_strain_cost(_dof.deformation(), relaxed_lat.vol() / double(_dof.size()));
    }

    //*******************************************************************************************

    double basis_cost(const ConfigDoF &_dof) {
      // mean square displacement distance in deformed coordinate system
      return (_dof.deformation() * _dof.displacement() * _dof.displacement().transpose() * _dof.deformation().transpose()).trace() / double(_dof.size());
    }
  }

  ConfigMapper::ConfigMapper(PrimClex &_pclex, double _lattice_weight, double _max_volume_change/*=0.25*/, double _tol/*=TOL*/) :
    m_pclex(&_pclex), m_lattice_weight(_lattice_weight), m_max_volume_change(_max_volume_change),m_tol(max(1e-9,_tol)){
    //squeeze lattice_weight into [0,1] if necessary
    m_lattice_weight = max(min(_lattice_weight, 1.0), 1e-9);
    ParamComposition param_comp(_pclex.get_prim());
    m_fixed_components=param_comp.fixed_components();

  }
  
  //*******************************************************************************************
  const std::vector<Lattice>& ConfigMapper::_lattices_of_vol(Index prim_vol) const{
    if(!valid_index(prim_vol)){
      throw std::runtime_error("Cannot enumerate lattice of volume " + std::to_string(prim_vol) + ", which is out of bounds.\n");
    }
    auto it = m_superlat_map.find(prim_vol);
    if(it !=m_superlat_map.end())
      return it->second;

    std::vector<Lattice> lat_vec;
    SupercellEnumerator<Lattice> enumerator(primclex().get_prim().lattice(), primclex().get_prim().point_group(), prim_vol, prim_vol + 1);

    Index l = 0;
    //std::cout << "min_vol is " << min_vol << "max_vol is " << max_vol << "\n";
    for(auto it = enumerator.begin(); it != enumerator.end(); ++it) {
      //std::cout << "Enumeration step " << l++ << " best cost is " << best_cost << "\n";
      lat_vec.push_back(niggli(*it, primclex().get_prim().point_group(), m_tol));
    }
    
    return m_superlat_map[prim_vol]=std::move(lat_vec);

  }
  //*******************************************************************************************

  bool ConfigMapper::import_structure_occupation(const fs::path &pos_path,
                                                 std::string &imported_name,
                                                 jsonParser &relaxation_properties,
                                                 bool robust_flag,
                                                 bool rotate_flag,
                                                 bool strict_flag,
                                                 double vol_tol) const {

    try {
      return import_structure_occupation(BasicStructure<Site>(pos_path),
                                         imported_name,
                                         relaxation_properties,
                                         robust_flag,
                                         rotate_flag,
                                         strict_flag,
                                         vol_tol);
    }
    catch(const std::exception &ex) {
      throw std::runtime_error(std::string("Could not successfully import structure ") + pos_path.string() + ":\n" + ex.what());
    }


  }

  //*******************************************************************************************

  bool ConfigMapper::import_structure_occupation(const BasicStructure<Site> &_struc,
                                                 std::string &imported_name,
                                                 jsonParser &relaxation_properties,
                                                 bool robust_flag,
                                                 bool rotate_flag,
                                                 bool strict_flag,
                                                 double vol_tol) const{
    return import_structure_occupation(_struc,
                                       nullptr,
                                       imported_name,
                                       relaxation_properties,
                                       robust_flag,
                                       rotate_flag,
                                       strict_flag,
                                       vol_tol);
  }
  //*******************************************************************************************

  bool ConfigMapper::import_structure_occupation(const BasicStructure<Site> &_struc,
                                                 const Configuration *hint_ptr,
                                                 std::string &imported_name,
                                                 jsonParser &relaxation_properties,
                                                 bool robust_flag,
                                                 bool rotate_flag,
                                                 bool strict_flag,
                                                 double vol_tol) const{

    //Indices for Configuration index and permutation operation index
    ConfigDoF tconfigdof;
    Lattice mapped_lat;
    bool new_config_flag;

    relaxation_properties.put_obj();

    if(!struc_to_configdof(_struc,
                           tconfigdof,
                           mapped_lat,
                           robust_flag,
                           rotate_flag,
                           vol_tol))
      throw std::runtime_error("Structure is incompatible with PRIM.");

    relaxation_properties["basis_deformation"] = ConfigMapping::basis_cost(tconfigdof);
    relaxation_properties["lattice_deformation"] = ConfigMapping::strain_cost(_struc.lattice(), tconfigdof);
    relaxation_properties["volume_relaxation"] = tconfigdof.deformation().determinant();
    Eigen::Matrix3d E = StrainConverter::green_lagrange(tconfigdof.deformation());
    std::vector<double> Evec(6);
    Evec[0] = (E(0, 0));
    Evec[1] = (E(1, 1));
    Evec[2] = (E(2, 2));
    Evec[3] = (sqrt(2.0) * E(1, 2));
    Evec[4] = (sqrt(2.0) * E(0, 2));
    Evec[5] = (sqrt(2.0) * E(0, 1));
    relaxation_properties["relaxation_strain"] = Evec;

    ConfigDoF relaxed_occ;

    relaxed_occ.set_occupation(tconfigdof.occupation());
    if(hint_ptr != nullptr) {
      ConfigDoF canon_relaxed_occ, canon_ideal_occ;
      Supercell const &scel(hint_ptr->get_supercell());
      if(mapped_lat.is_equivalent(scel.get_real_super_lattice())) {
        if(strict_flag && relaxed_occ.occupation() == (hint_ptr->configdof()).occupation()) {
          // config is unchanged
          imported_name = hint_ptr->name();
          return false;
        }
        canon_relaxed_occ = relaxed_occ.canonical_form(scel.permute_begin(), scel.permute_end(), m_tol);

        canon_ideal_occ = (hint_ptr->configdof()).canonical_form(scel.permute_begin(), scel.permute_end(), m_tol);
        //std::cout << "canon_relaxed_occ.occupation() is " << canon_relaxed_occ.occupation() << "\n";
        //std::cout << "canon_ideal_occ.occupation() is " << canon_ideal_occ.occupation() << "\n";

        if(canon_relaxed_occ.occupation() == canon_ideal_occ.occupation()) {
          // config is unchanged
          imported_name = hint_ptr->name();
          return false;
        }
      }
    }
    Supercell::permute_const_iterator permute_it;

    Index import_scel_index = primclex().add_supercell(mapped_lat), import_config_index;

    Configuration import_config(primclex().get_supercell(import_scel_index), jsonParser(), relaxed_occ);

    if(strict_flag) {
      permute_it = primclex().get_supercell(import_scel_index).permute_begin();
      new_config_flag = primclex().get_supercell(import_scel_index).add_canon_config(import_config, import_config_index);
      imported_name = primclex().get_supercell(import_scel_index).get_config(import_config_index).name();
    }
    else {
      permute_it = primclex().get_supercell(import_scel_index).permute_begin();
      new_config_flag = primclex().get_supercell(import_scel_index).add_config(import_config, import_config_index, permute_it);
      imported_name = primclex().get_supercell(import_scel_index).get_config(import_config_index).name();
    }

    return new_config_flag;
  }

  //*******************************************************************************************

  bool ConfigMapper::import_structure(const fs::path &pos_path,
                                      std::string &imported_name,
                                      jsonParser &relaxation_properties,
                                      bool robust_flag,
                                      bool rotate_flag,
                                      bool strict_flag,
                                      double vol_tol) const{

    try {
      return import_structure(BasicStructure<Site>(pos_path),
                              imported_name,
                              relaxation_properties,
                              robust_flag,
                              rotate_flag,
                              strict_flag,
                              vol_tol);
    }
    catch(const std::exception &ex) {
      throw std::runtime_error(std::string("Could not successfully import structure ") + pos_path.string() + ":\n" + ex.what());
    }


  }

  //*******************************************************************************************
  bool ConfigMapper::import_structure(const BasicStructure<Site> &_struc,
                                      std::string &imported_name,
                                      jsonParser &relaxation_properties,
                                      bool robust_flag,
                                      bool rotate_flag,
                                      bool strict_flag,
                                      double vol_tol) const{

    //Indices for Configuration index and permutation operation index
    Supercell::permute_const_iterator permute_it;

    ConfigDoF tconfigdof;
    Lattice mapped_lat;
    bool new_config_flag;

    if(!struc_to_configdof(_struc,
                           tconfigdof,
                           mapped_lat,
                           robust_flag,
                           rotate_flag,
                           vol_tol))
      throw std::runtime_error("Structure is incompatible with PRIM.");

    relaxation_properties["basis_deformation"] = ConfigMapping::basis_cost(tconfigdof);
    relaxation_properties["lattice_deformation"] = ConfigMapping::strain_cost(_struc.lattice(), tconfigdof);
    relaxation_properties["volume_change"] = tconfigdof.deformation().determinant();

    Index import_scel_index = primclex().add_supercell(mapped_lat), import_config_index;

    Configuration import_config(primclex().get_supercell(import_scel_index), jsonParser(), tconfigdof);

    if(strict_flag) {
      permute_it = primclex().get_supercell(import_scel_index).permute_begin();
      new_config_flag = primclex().get_supercell(import_scel_index).add_canon_config(import_config, import_config_index);
      imported_name = primclex().get_supercell(import_scel_index).get_config(import_config_index).name();
    }
    else {
      permute_it = primclex().get_supercell(import_scel_index).permute_begin();
      new_config_flag = primclex().get_supercell(import_scel_index).add_config(import_config, import_config_index, permute_it);
      imported_name = primclex().get_supercell(import_scel_index).get_config(import_config_index).name();
    }

    return new_config_flag;

  }

  //*******************************************************************************************
  Lattice find_nearest_super_lattice(const Lattice &prim_lat,
                                     const Lattice &relaxed_lat,
                                     const SymGroup &sym_group,
                                     Eigen::MatrixXd &deformation,
                                     Eigen::MatrixXd &trans_mat,
                                     Index min_vol,
                                     Index max_vol,
                                     double _tol) {
    Lattice best_lat;
    if(!valid_index(max_vol))
      max_vol = ceil(std::abs(relaxed_lat.vol()) / std::abs(prim_lat.vol()));

    double best_cost = 10e10;
    //We only bother checking pre-existing supercells of min_vol <= volume <=max_vol;
    SupercellEnumerator<Lattice> enumerator(prim_lat, sym_group, min_vol, max_vol + 1);

    Index l = 0;
    //std::cout << "min_vol is " << min_vol << "max_vol is " << max_vol << "\n";
    for(auto it = enumerator.begin(); it != enumerator.end(); ++it) {
      //std::cout << "Enumeration step " << l++ << " best cost is " << best_cost << "\n";
      Lattice tlat = niggli(*it, sym_group, _tol);

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

  //*******************************************************************************************
  Lattice find_nearest_super_lattice(const Lattice &prim_lat,
                                     const Lattice &relaxed_lat,
                                     const SymGroup &sym_group,
                                     Eigen::MatrixXd &deformation,
                                     Eigen::MatrixXd &trans_mat,
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
  bool ConfigMapper::struc_to_configdof(const BasicStructure<Site> &struc,
                                        ConfigDoF &mapped_configdof,
                                        Lattice &mapped_lat,
                                        bool robust_flag,
                                        bool rotate_flag,
                                        double vol_tol) const{
    
      bool valid_mapping(false);
      // If structure's lattice is a supercell of the primitive lattice, then import as ideal_structure
      if(!robust_flag && struc.lattice().is_supercell_of(primclex().get_prim().lattice(), m_tol)) {
        valid_mapping = ideal_struc_to_configdof(struc, mapped_configdof, mapped_lat);
        //std::cout << "valid_mapping is " << valid_mapping << "\n";
        valid_mapping = valid_mapping && ConfigMapping::basis_cost(mapped_configdof) < (10 * m_tol);
      }

      // If structure's lattice is not a supercell of the primitive lattice, then import as deformed_structure
      if(!valid_mapping) // if not a supercell or robust_flag=true, treat as deformed
        valid_mapping = deformed_struc_to_configdof(struc,
                                                    mapped_configdof,
                                                    mapped_lat,
                                                    rotate_flag,
                                                    vol_tol);

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
    bool ConfigMapper::ideal_struc_to_configdof(BasicStructure<Site> struc,
                                                ConfigDoF &mapped_configdof,
                                                Lattice &mapped_lat) const{
      if(!struc.lattice().is_supercell_of(primclex().get_prim().lattice(),  m_tol)) {
        std::cerr << "CRITICAL ERROR: In ideal_struc_to_configdof(), primitive structure does not tile the provided\n"
                  << "                superstructure. Please use deformed_struc_to_configdof() instead.\n"
                  << "                Exiting...\n";
        exit(1);
      }

      mapped_lat = niggli(struc.lattice(), primclex().get_prim().point_group(), m_tol);
      Supercell scel(&primclex(), mapped_lat);

      // We modify the idealized structure so that the lattice matches the one in the list
      Matrix3<double> trans_mat;
      mapped_lat.is_supercell_of(struc.lattice(), primclex().get_prim().point_group(), trans_mat,  m_tol);
      struc.set_lattice(Lattice(struc.lattice().lat_column_mat()*trans_mat), CART);
      struc.set_lattice(mapped_lat, FRAC);

      return ConfigMap_impl::struc_to_configdof(scel, struc, mapped_configdof, true, m_tol);
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
     * that cannot be realized in the system descript by the primclex).
     *
     * What this does NOT do:
     *    -Check if the imported Structure is the same as one in a smaller Supercell
     *
     */
    //*******************************************************************************************
    bool ConfigMapper::deformed_struc_to_configdof(const BasicStructure<Site> &struc,
                                     ConfigDoF &mapped_configdof,
                                     Lattice &mapped_lat,
                                     bool rotate_flag,
                                     double vol_tol) const{
      //squeeze lattice_weight into [0,1] if necessary
      double lw=m_lattice_weight;
      double bw = 1.0 - lw;
      // make vol_tol >= m_tol
      vol_tol = max(m_tol, vol_tol);

      //Add new Supercell if it doesn't exist already. Use primitive point group to check for equivalence and
      //store transformation matrix
      Eigen::MatrixXd deformation;
      double num_atoms = double(struc.basis.size());
      int min_vol, max_vol;

      if(m_fixed_components.size()>0){
        std::string tcompon=m_fixed_components[0].first;
        int ncompon(0);
        for(Index i=0; i<struc.basis.size(); i++){
          if(struc.basis[i].occ_name()==tcompon)
            ncompon++;
        }
        min_vol=ncompon/int(m_fixed_components[0].second);
        max_vol=min_vol;
      }
      else{
        // Try to narrow the range of supercell volumes -- the best bounds are obtained from
        // the convex hull of the end-members, but we need to wait for improvements to convex hull
        // routines
        
        // min_vol assumes no vacancies -- best case scenario
        min_vol=ceil((num_atoms / double(primclex().get_prim().basis.size())) - m_tol);
        //For practical purposes, set maximum Va fraction to 0.75
        double max_va_fraction = min(0.75, double(primclex().get_prim().max_possible_vacancies()) / double(primclex().get_prim().basis.size()));
        // This is for the worst case scenario -- lots of vacancies
        max_vol = ceil(num_atoms / (double(primclex().get_prim().basis.size()) * (1.0 - max_va_fraction)) - m_tol);

        if(max_va_fraction > TOL) {
          //Nvol is rounded integer volume-- assume that answer is within 30% of this volume, and use it to tighten our bounds
          int Nvol = round(std::abs(struc.lattice().vol() / primclex().get_prim().lattice().vol()));
          int new_min_vol = min(max_vol, max(round((1.0 - vol_tol) * double(Nvol)), min_vol));
          int new_max_vol = max(min_vol, min(round((1.0 + vol_tol) * double(Nvol)), max_vol));
          max_vol = new_max_vol;
          min_vol = new_min_vol;
        }
      }

      //std::cout << "max_va_fraction: " << max_va_fraction << "   Volume range: " << min_vol << " to " << max_vol << "\n";
      Eigen::MatrixXd ttrans_mat, tF, best_trans;
      double strain_cost(1e10), basis_cost(1e10), best_cost(1e20), tot_cost, best_strain_cost, best_basis_cost;
      ConfigDoF tdof;
      BasicStructure<Site> tstruc(struc);
      Lattice tlat;
      //std::cout << "First pass: ";
      for(Index i_vol = min_vol; i_vol <= max_vol; i_vol++) {
        //std::cout << "v=" << i_vol << "   ";
        tlat = find_nearest_super_lattice(primclex().get_prim().lattice(), struc.lattice(), primclex().get_prim().point_group(), tF, ttrans_mat, _lattices_of_vol(i_vol), m_tol);
        strain_cost = lw * LatticeMap::calc_strain_cost(tF, struc.lattice().vol() / num_atoms);

        if(best_cost < strain_cost)
          continue;
        //std::cout << "best_cost " << best_cost << "   strain_cost " << strain_cost << "\n";
        tstruc = struc;

        tstruc.set_lattice(Lattice(tF * Eigen::MatrixXd(tlat.lat_column_mat())), CART);

        if(rotate_flag) {
          tF = StrainConverter::right_stretch_tensor(tF);
          tstruc.set_lattice(Lattice(tF * Eigen::MatrixXd(tlat.lat_column_mat())), FRAC);
        }

        Supercell scel(&primclex(), tlat);

        if(!ConfigMap_impl::struc_to_configdof(scel, tstruc, tdof, true, m_tol))
          continue;
        basis_cost = bw * ConfigMapping::basis_cost(tdof);
        tot_cost = strain_cost + basis_cost;
        //std::cout << "\n**Starting strain_cost = " << strain_cost << ";   and basis_cost = " << basis_cost << "  TOTAL: " << strain_cost + basis_cost << "\n";
        if(tot_cost < best_cost) {
          best_cost = tot_cost - m_tol;
          best_strain_cost = strain_cost;
          best_basis_cost = basis_cost;
          best_trans = Eigen::MatrixXd(tlat.inv_lat_column_mat()) * tF.inverse() * Eigen::MatrixXd(struc.lattice().lat_column_mat());
          //std::cout << "tF is:\n" << tF << "\n and N is:\n" << best_trans << "\n";
          mapped_configdof = tdof;

          mapped_lat = tlat;
        }
      }
      //std::cout << "\n\n";

      //std::cout << "Second pass:\n";
      for(Index i_vol = min_vol; i_vol <= max_vol; i_vol++) {
        //std::cout << "  vol = " << i_vol << "\n";
        const std::vector<Lattice>& lat_vec=_lattices_of_vol(i_vol);
        bool break_early(false);
        for(auto it = lat_vec.cbegin(); it != lat_vec.cend() && !break_early; ++it) {
          Supercell scel(&primclex(), *it);

          //Determine best mapping for this supercell

          //Initialize with simplest mapping onto supercell 'i', so that we don't change the crystal setting unnecessarily
          tF = struc.lattice().lat_column_mat() * it->inv_lat_column_mat();

          strain_cost = lw * LatticeMap::calc_strain_cost(tF, struc.lattice().vol() / num_atoms);

          // If simplest mapping seems viable, check it further
          if(strain_cost < best_cost) {
            tstruc = struc;
            tstruc.set_lattice(Lattice(tF * Eigen::MatrixXd(it->lat_column_mat())), CART);
            if(rotate_flag) {
              tF = StrainConverter::right_stretch_tensor(tF);
              tstruc.set_lattice(Lattice(tF * Eigen::MatrixXd(it->lat_column_mat())), FRAC);
            }

            if(!ConfigMap_impl::struc_to_configdof(scel, tstruc, tdof, true, m_tol))
              break;
            basis_cost = bw * ConfigMapping::basis_cost(tdof);
            //std::cout << "\n**Starting strain_cost = " << strain_cost << ";   and basis_cost = " << basis_cost << "  TOTAL: " << strain_cost + basis_cost << "\n";
            tot_cost = strain_cost + basis_cost;
            //std::cout << "    scel.name = " << scel.get_name()  << "  simple map: strain_cost " << strain_cost << "   best_cost " << best_cost << "    tot_cost " << tot_cost << "\n";
            if(tot_cost < best_cost) {
              best_cost = tot_cost - m_tol;
              best_strain_cost = strain_cost;
              best_basis_cost = basis_cost;
              best_trans = Matrix3<double>::identity();
              //std::cout << "tF is:\n" << tF << "\n and N is:\n" << best_trans << "\n";
              mapped_configdof = tdof;

              mapped_lat = *it;
            }
          } // Done checking simplest mapping

          // If the simplest mapping is best, we have avoided a lot of extra work, but we still need to check for
          // non-trivial mappings that are better than both the simplest mapping and the best found mapping
          LatticeMap strainmap(*it, struc.lattice(), round(num_atoms), m_tol, 1);
          strain_cost = lw * strainmap.strain_cost();
          if(best_cost < strain_cost)
            strain_cost = lw * strainmap.next_mapping_better_than(best_cost).strain_cost();

          while(strain_cost < best_cost) {  // only enter loop if there's a chance of improving on current best
            tstruc = struc;
            // We modify the deformed structure so that its lattice is a deformed version of the nearest ideal lattice
            // Don't need matrixN if we use set_lattice(CART), because matrixF depends on matrixN implicitly
            tstruc.set_lattice(Lattice(strainmap.matrixF() * Eigen::MatrixXd(it->lat_column_mat())), CART);
            if(rotate_flag) {
              tF = StrainConverter::right_stretch_tensor(strainmap.matrixF());
              tstruc.set_lattice(Lattice(tF * Eigen::MatrixXd(it->lat_column_mat())), FRAC);
            }

            if(!ConfigMap_impl::struc_to_configdof(scel, tstruc, tdof, true, m_tol)) {
              break_early = true;
              break;
              //no longer unexpected
              //throw std::runtime_error("Unexpected error in deformed_struc_to_config_dof(). This should never happen!\n");
            }
            //std::cout << "New strain_cost = " << strain_cost << ";   and basis_cost = " << basis_cost << "  TOTAL: " << strain_cost + basis_cost << "\n";
            //std::cout << "  Compare -> best_cost = " << best_cost << "\n";
            basis_cost = bw * ConfigMapping::basis_cost(tdof);
            tot_cost = strain_cost + basis_cost;
            //std::cout << "      complex map: best_cost " << best_cost << "    strain_cost " << strain_cost << "    tot_cost " << tot_cost << "\n";
            if(tot_cost < best_cost) {
              //std::cout << "Old best_trans, with cost " << best_cost << ":\n" << best_trans << "\n";
              best_cost = tot_cost - m_tol;
              best_strain_cost = strain_cost;
              best_basis_cost = basis_cost;
              best_trans = strainmap.matrixN();
              //std::cout << "New best_trans, with cost " << best_cost << ":\n" << best_trans << "\n";
              mapped_configdof = tdof;

              mapped_lat = *it;
            }
            // This finds first decomposition:
            //            struc.lattice() =  deformation*(supercell_list[import_scel_index].get_real_super_lattice())*equiv_mat
            //   that has cost function less than best_cost
            strain_cost = lw * strainmap.next_mapping_better_than(best_cost).strain_cost();
          }
        }
      }
      //std::cout << "best_cost = " << best_cost << "    best_strain_cost = " << best_strain_cost << "    best_basis_cost = " << best_basis_cost << "\n";
      if(best_cost > 1e9) {
        //std::cerr << "WARNING: In Supercell::import_deformed_structure(), no successful mapping was found for Structure " << src << "\n"
        //        << "         This Structure may be incompatible with the ideal crystal specified by the PRIM file. \n";
        return false;
      }
      //std::cout << "FINAL COST IS: " << best_cost << "\nFINAL TRANS MAT IS:\n" << best_trans << "\n\n";
      return true;

    }

  namespace ConfigMap_impl{
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
                          Eigen::MatrixXd &cost_matrix) {


      double inf = 10E10;
      //if(cost_matrix.rows()!=scel.num_sites() || cost_matrix.cols()!=scel.num_sites())
      cost_matrix = Eigen::MatrixXd::Constant(scel.num_sites(), scel.num_sites(), inf);
      Index inf_counter;
      double dist;
      // loop through all the sites of the structure
      Index j = 0;
      for(; j < rstruc.basis.size(); j++) {
        Coordinate current_relaxed_coord(rstruc.basis[j](FRAC), scel.get_real_super_lattice(), FRAC);
        current_relaxed_coord(CART) += trans(CART);
        // loop through all the sites in the supercell
        inf_counter = 0;
        for(Index i = 0; i < scel.num_sites(); i++) {

          // Check if relaxed atom j is allowed on site i
          // If so, populate cost_matrix normally
          if(scel.get_prim().basis[scel.get_b(i)].contains(rstruc.basis[j].occ_name())) {
            dist = scel.coord(i).min_dist(current_relaxed_coord);
            cost_matrix(i, j) = dist * dist;
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

      // See if a very simple (non-optimal) assignment is possible using the cost matrix.
      // if not, the cost matrix is not viable -- only an issue in the presence of vacancies
      //JCT: I'm not sure if this could fail in more complicated cases
      //     (e.g., A-B alloying on one sublattice and A-Va disorder on another sublattice)
      //  -- Update: it can fail.  For now, just rely on hungarian method to decide if assignment is possible
      /*
        std::vector<bool> assignable(scel.num_sites(), true);
        for(Index i = 0; i < scel.num_sites(); i++) {
        Index ii = 0;
        for(; ii < scel.num_sites(); ii++) {
        if(assignable[ii] && cost_matrix(i, ii) < (inf - TOL)) {
        assignable[ii] = false;
        break;
        }
        }
        if(ii == scel.num_sites()) {
        std::cerr << "Bail at 3\n"
        << " cost_matrix is \n"<< cost_matrix << "\n";
        return false;
        }
        }
      */
      return true;
    }

    //***************************************************************************************************
    // New mapping routine. Return an ideal configuration corresponding to a relaxed structure.
    // Options:
    //
    // translate_flag = true means that the rigid-translations are removed. (typically this option should be used)
    //
    // translate_flag = false means that rigid translations are not considered. (probably don't want to use this since structures should be considered equal if they are related by a rigid translation).
    //
    bool struc_to_configdof(const Supercell &scel,
                            BasicStructure<Site> rstruc,
                            ConfigDoF &config_dof,
                            const bool translate_flag,
                            const double _tol) {
      //std::cout << "CONVERTING TO CONFIGDOF:\n";
      //rstruc.print(std::cout);
      //std::cout << std::endl;
      // clear config_dof and set its deformation
      config_dof.clear();

      Eigen::Matrix3d deformation = Eigen::Matrix3d(rstruc.lattice().coord_trans(FRAC) * scel.get_real_super_lattice().coord_trans(CART));

      config_dof.set_deformation(deformation);

      // un-deform rstruc
      rstruc.set_lattice(scel.get_real_super_lattice(), FRAC);


      //Initialize everything

      Eigen::MatrixXd cost_matrix(scel.num_sites(), scel.num_sites());
      std::vector<Index> optimal_assignments(scel.num_sites()), best_assignments(scel.num_sites());
      //BasicStructure<Site> best_ideal_struc(rstruc);
      Coordinate ttrans(Vector3<double>(0, 0, 0), rstruc.lattice(), FRAC), best_trans(Vector3<double>(0, 0, 0), rstruc.lattice(), FRAC);
      Array<int> assignment_bitstring(scel.num_sites());
      double min_mean = 10E10;
      double trans_dist;

      // We want to get rid of translations.
      // trans_coord is a vector from IDEAL to RELAXED
      // Subtract this from every rstruc coordinate
      Index num_translations(1);
      //if(translate_flag == true)
      num_translations += scel.get_prim().basis.size();
      //num_translations = rstruc.basis.size();
      //std::cout << "num_translations is " << num_translations << "\n";
      for(Index n = 0; n < num_translations; n++) {
        double mean;

        //shift_struc has **ideal lattice**
        //BasicStructure<Site> shift_struc(rstruc);


        if(n > 0 && !scel.get_prim().basis[n - 1].contains(rstruc.basis[0].occ_name()))
          continue;

        Coordinate ref_coord(rstruc.basis[0]);

        if(n > 0)
          ref_coord(FRAC) = scel.coord((n - 1) * scel.volume())(FRAC);


        trans_dist = rstruc.basis[0].min_dist(ref_coord, ttrans);

        ttrans.set_lattice(scel.get_prim().lattice(), CART);
        //std::cout << "Before:  ttrans " << ttrans(FRAC) << "; V_number: " << ttrans.voronoi_number() << "\n";
        ttrans.within();// <-- should be voronoi_within()?
        //std::cout << "After:  ttrans " << ttrans(FRAC) << "; V_number: " << ttrans.voronoi_number() << "\n\n\n";
        ttrans.set_lattice(rstruc.lattice(), CART);
        trans_dist = ttrans(CART).length();
        //shift_struc -= ttrans;
        if(!ConfigMap_impl::calc_cost_matrix(scel, rstruc, ttrans, cost_matrix)) {
          //std::cerr << "In Supercell::struc_to_config. Cannot construct cost matrix." << std::endl;
          //std::cerr << "This message is probably OK, if you are using translate_flag == true." << std::endl;
          //continue;
          return false;
        }

        //std::cout << "cost_matrix is\n" << cost_matrix <<  "\n\n";
        // The mapping routine is called here
        mean = hungarian_method(cost_matrix, optimal_assignments, _tol);
        //std::cout << "mean is " << mean << " and stddev is " << stddev << "\n";
        // add small penalty (~_tol) for larger translation distances, so that shortest equivalent translation is used
        mean += _tol * trans_dist / 10.0;
        if(mean < min_mean) {
          //std::cout << "mean " << mean << " is better than min_mean " << min_mean <<"\n";

          // new best
          best_assignments = optimal_assignments;

          // best shifted structure
          best_trans = ttrans;

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
      double avg2_disp(0);
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

          Coordinate ideal_coord(scel.coord(i)(FRAC), rstruc.lattice(), FRAC);

          (rstruc.basis[best_assignments[i]] + best_trans).min_dist(ideal_coord, disp_coord);
          //std::cout << "min_dist" << std::endl;
          //std::cout << rstruc.basis[optimal_assignments(i)].min_dist(pos_coord, disp_coord);
          for(Index j = 0; j < 3; j++) {
            config_dof.disp(i)[j] = disp_coord(CART)[j];
          }
          avg_disp += config_dof.disp(i);
          avg2_disp += config_dof.disp(i).squaredNorm();
        }
      }

      avg_disp /= double(rstruc.basis.size());

      //std::cout << "rstruc is:\n";
      //rstruc.print(std::cout);
      //std::cout << "\n\nand avg_disp is" << avg_disp.transpose() << "\n\n";
      //std::cout << "\n\nand displacement field (Nx3) is\n" << config_dof.displacement().transpose() << "\n\n";

      // subtract off average displacement
      for(Index i = 0; i < best_assignments.size(); i++) {
        if(best_assignments[i] < rstruc.basis.size()) {
          config_dof.disp(i) -= avg_disp;
          // suppress small ugly numbers.
          for(Index j = 0; j < 3; j++) {
            if(almost_zero(config_dof.disp(i)[j], 1e-8))
              config_dof.disp(i)[j] = 0;
          }
        }
      }
      // End of filling displacements


      // Make the assignment bitstring
      //
      // Loop through all supercell sites
      for(Index i = 0; i < scel.num_sites(); i++) {

        // Retrieve the identity of the atom assigned to site i
        // Any value of the assignment vector larger than the number
        // of sites in the relaxed structure is by construction
        // specified as a vacancy.
        std::string rel_basis_atom;
        if(best_assignments[i] >= rstruc.basis.size()) {
          rel_basis_atom = "Va";
        }
        else {
          rel_basis_atom = rstruc.basis[best_assignments[i]].occ_name();
        }

        // How many possible atoms are allowed at the basis site?


        if(!scel.get_prim().basis[scel.get_b(i)].contains(rel_basis_atom, assignment_bitstring[i])) {
          //std::cout << "best_assignments is " << best_assignments << "\n";
          //std::cout << "at site " << i << " corresponding to basis " << scel.get_b(i) << "  attempting to assign type " << rel_basis_atom << "\n";
          //std::cout << "Cost Matrix is \n" << cost_matrix << "\n";
          //std::cerr << "CRITICAL ERROR: In Supercell::struc_to_configdof atoms of relaxed/custom structure are incompatible\n"
          //        << "                with the number or type of atomic species allowed in PRIM. Exiting...\n";
          //exit(1);
          return false;
        }
      }

      config_dof.set_occupation(assignment_bitstring);
      return true;
    }
  }
}
