#include "casm/clex/ConfigMapping.hh"

#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/casm_io/json_io/container.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/Configuration.hh"
#include "casm/clex/Supercell.hh"
#include "casm/clex/ParamComposition.hh"
#include "casm/strain/StrainConverter.hh"
#include "casm/crystallography/Lattice.hh"
#include "casm/crystallography/Niggli.hh"
#include "casm/crystallography/LatticeMap.hh"
#include "casm/crystallography/SupercellEnumerator.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/symmetry/PermuteIterator.hh"
#include "casm/completer/Handlers.hh"

namespace CASM {
  ConfigDoF MappedConfig::to_configdof() const {
    throw std::runtime_error("MappedConfig::to_configdof() not yet implemented");
    return ConfigDoF(0, 0, std::map<DoFKey, DoFSetInfo>(), std::map<DoFKey, std::vector<DoFSetInfo> >(), 0.);
  }

  MappedConfig &MappedConfig::apply_sym(PermuteIterator const &it) {
    throw std::runtime_error("MappedConfig::apply_sym() is not implemented!");
    return *this;
  }


  namespace ConfigMapping {
    double strain_cost(double relaxed_lat_vol, const MappedConfig &_dof, const Index Nsites) {
      return LatticeMap::calc_strain_cost(_dof.deformation, relaxed_lat_vol / double(max(Nsites, Index(1))));
    }

    //*******************************************************************************************

    double basis_cost(const MappedConfig &_dof, Index Nsites) {
      // mean square displacement distance in deformed coordinate system
      return (_dof.deformation * _dof.displacement * _dof.displacement.transpose() * _dof.deformation.transpose()).trace() / double(max(Nsites, Index(1)));
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

      for(auto it = enumerator.begin(); it != enumerator.end(); ++it) {

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

      return best_lat;
    }

    ConfigMapperResult structure_mapping(Structure &host, Structure &other, double lattice_weight) {
      const PrimClex &pclex = PrimClex(host, null_log());
      ConfigMapper tmp_mapper(pclex, lattice_weight, 0.0);
      tmp_mapper.set_max_va_frac(0.0);
      return tmp_mapper.import_structure_occupation(SimpleStructure(other));
    }
  }

  //*******************************************************************************************

  ConfigMapper::ConfigMapper(const PrimClex &_pclex,
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
    ParamComposition param_comp(_pclex.prim());
    m_fixed_components = param_comp.fixed_components();
    m_max_volume_change = max(m_tol, _max_volume_change);
  }

  //*******************************************************************************************

  void ConfigMapper::force_lattices(const std::vector<std::string> &lattice_names) const {
    //Turn the lattice names into actual lattices
    for(auto &name : lattice_names) {
      //Somehow adding new supercells to the database is considered const (?)
      const Supercell &force_scel = CASM::make_supercell(*m_pclex, name);
      const Lattice &force_lat = force_scel.lattice();
      Index force_vol = force_scel.volume();
      m_forced_superlat_map[force_vol].push_back(force_lat);
    }
    return;
  }

  //*******************************************************************************************

  void ConfigMapper::unforce_lattices() const {
    m_forced_superlat_map.clear();
    return;
  }

  //*******************************************************************************************

  ConfigMapperResult ConfigMapper::import_structure_occupation(const fs::path &pos_path) const {
    BasicStructure<Site> tstruc(pos_path);
    return import_structure_occupation(SimpleStructure(tstruc));
  }

  //*******************************************************************************************

  ConfigMapperResult ConfigMapper::import_structure_occupation(const SimpleStructure &_struc) const {
    return import_structure_occupation(_struc, nullptr);
  }

  //*******************************************************************************************

  ConfigMapperResult ConfigMapper::import_structure_occupation(
    const SimpleStructure &_struc,
    const Configuration *hint_ptr) const {

    ConfigMapperResult result;
    result.structure = _struc;
    result.relaxation_properties.put_obj();

    //Indices for Configuration index and permutation operation index
    MappedConfig best_configdof, suggested_configdof;
    Lattice mapped_lat;
    double bc(1e20), sc(1e20), best_cost = 1e20, robust_cost = 1e20;
    bool valid_mapping;

    // -------------
    // if given a hint,
    // - try to map to hint configuration,
    // - get "best_cost" which will be used to limit the number of other
    //   possibly mappings that need to be considered
    // - get mapped_lat
    // - store "suggested_mapping" in 'relaxation_properties, storing mapping
    //   score for the "hinted" configuration
    if(hint_ptr != nullptr) {

      valid_mapping = ConfigMap_impl::struc_to_configdof(
                        *hint_ptr,
                        result.structure,
                        suggested_configdof,
                        result.best_assignment,
                        m_robust_flag, // translate_flag -- not sure what to use for this
                        m_tol);

      if(valid_mapping) {

        mapped_lat = (hint_ptr->supercell()).lattice();
        bc = ConfigMapping::basis_cost(suggested_configdof, result.structure.n_mol());
        sc = ConfigMapping::strain_cost(
               result.structure.lat_column_mat.determinant(),
               suggested_configdof,
               result.structure.n_mol());

        result.relaxation_properties["suggested_mapping"]["basis_deformation"] = bc;
        result.relaxation_properties["suggested_mapping"]["lattice_deformation"] = sc;
        result.relaxation_properties["suggested_mapping"]["volume_relaxation"] =
          suggested_configdof.deformation.determinant();

        best_cost = m_lattice_weight * sc + (1.0 - m_lattice_weight) * bc - m_tol;

      }
    }

    // -------------
    // do a general mapping, but with limit to configurations that will improve
    // on the "hint_cost"
    // - try to map to best configuration,
    // - get mapped_lat (which is the canonical_equivalent)
    // - store "suggested_mapping" in 'relaxation_properties, storing mapping
    //   score for the "hinted" configuration
    valid_mapping = struc_to_configdof(
                      result.structure,
                      best_configdof,
                      mapped_lat,
                      result.best_assignment,
                      result.cart_op,
                      best_cost);

    // if a valid mapping was found, it had a better cost than "hint"
    // - store relaxation properties in "best_mapping"
    if(valid_mapping) {

      bc = ConfigMapping::basis_cost(best_configdof, result.structure.n_mol());
      sc = ConfigMapping::strain_cost(
             result.structure.lat_column_mat.determinant(),
             best_configdof,
             result.structure.n_mol());
      robust_cost = m_lattice_weight * sc + (1.0 - m_lattice_weight) * bc - m_tol;
      result.relaxation_properties["best_mapping"]["basis_deformation"] = bc;
      result.relaxation_properties["best_mapping"]["lattice_deformation"] = sc;
      result.relaxation_properties["best_mapping"]["volume_relaxation"] = best_configdof.deformation.determinant();

    }
    // if no valid mapping was found, "suggested" (the mapping to hint config) is also "best"
    else {
      if(best_cost > 1e10) {
        result.success = false;
        result.fail_msg = "Structure is incompatible with PRIM.";
      }

      std::swap(best_configdof, suggested_configdof);
      result.relaxation_properties["best_mapping"] = result.relaxation_properties["suggested_mapping"];
    }

    // -------------
    // construct the 'relaxed' configuration from the mapping results
    // and get the 'it_canon'
    MappedConfig relaxed_occ;
    relaxed_occ.occupation = best_configdof.occupation;

    PermuteIterator it_canon;
    bool is_new_config(true);

    if(hint_ptr != nullptr) {
      const Supercell &scel(hint_ptr->supercell());
      if(mapped_lat.is_equivalent(scel.lattice())) {
        if(m_strict_flag && relaxed_occ.occupation == (hint_ptr->configdof()).occupation()) {
          // mapped config is "hint" config
          is_new_config = false;
          result.config = notstd::make_unique<Configuration>(*hint_ptr);

          it_canon = hint_ptr->supercell().sym_info().permute_begin();
        }
        else {

          PermuteIterator relaxed_it_canon = Configuration(scel, jsonParser(), relaxed_occ.to_configdof()).to_canonical();
          PermuteIterator ideal_rev_it_canon = hint_ptr->from_canonical();
          it_canon = ideal_rev_it_canon * relaxed_it_canon;

          if(relaxed_occ.occupation == copy_apply(it_canon.inverse(), *hint_ptr).occupation()) {
            // mapped config is "hint" config
            is_new_config = false;
            result.config = notstd::make_unique<Configuration>(*hint_ptr);

          }
        }
      }
    }

    if(is_new_config) {
      std::shared_ptr<Supercell> shared_scel = std::make_shared<Supercell>(&primclex(), mapped_lat);
      result.config = notstd::make_unique<Configuration>(shared_scel, jsonParser(), relaxed_occ.to_configdof());

      if(m_strict_flag) {
        it_canon = shared_scel->sym_info().permute_begin();
      }
      else {
        it_canon = result.config->to_canonical();
      }
    }


    // calculate and store:
    // - 'relaxation_deformation'
    // - 'relaxation_displacement'
    // - cart_op
    // - best_assignement (excluding vacancies)

    // transform deformation tensor to match canonical form and apply operation to cart_op

    MappedConfig trans_configdof = copy_apply(it_canon, best_configdof);
    result.relaxation_properties["best_mapping"]["relaxation_deformation"] = trans_configdof.deformation;
    result.relaxation_properties["best_mapping"]["relaxation_displacement"] = trans_configdof.displacement.transpose();

    result.cart_op = it_canon.sym_op().matrix() * result.cart_op;

    // compose permutations
    std::vector<Index> tperm = it_canon.combined_permute().permute(result.best_assignment);

    //copy non-vacancy part of permutation into best_assignment
    result.best_assignment.resize(result.structure.n_mol());
    Index num_atoms = result.structure.n_mol();
    std::copy_if(tperm.cbegin(), tperm.cend(),
                 result.best_assignment.begin(),
    [num_atoms](Index i) {
      return i < num_atoms;
    });


    //result.structure.set_lattice(Lattice(result.cart_op.transpose()*best_configdof.deformation * mapped_lat.lat_column_mat()), CART);
    //result.structure.set_lattice(Lattice(best_configdof.deformation * mapped_lat.lat_column_mat()), FRAC);
    result.structure.lat_column_mat = result.cart_op.transpose() * best_configdof.deformation * mapped_lat.lat_column_mat();
    result.structure.mol_info.coords = result.cart_op * result.structure.mol_info.coords;
    result.structure.atom_info.coords = result.cart_op * result.structure.atom_info.coords;

    result.success = true;

    return result;
  }

  //*******************************************************************************************

  ConfigMapperResult ConfigMapper::import_structure(const fs::path &pos_path) const {
    return import_structure(SimpleStructure(pos_path));
  }

  //*******************************************************************************************
  ConfigMapperResult ConfigMapper::import_structure(const SimpleStructure &_struc) const {

    ConfigMapperResult result;
    result.structure = _struc;

    //Indices for Configuration index and permutation operation index
    PermuteIterator it_canon;

    MappedConfig best_configdof;
    Lattice mapped_lat;

    // -------------
    // do a general mapping
    // - try to map to best configuration,
    // - get mapped_lat (which is the canonical_equivalent)
    bool valid_mapping = struc_to_configdof(
                           result.structure,
                           best_configdof,
                           mapped_lat, // mappe
                           result.best_assignment,
                           result.cart_op);
    if(!valid_mapping) {
      result.success = false;
      result.fail_msg = "Structure is incompatible with PRIM.";
      return result;
    }

    // store "best_mapping" in 'relaxation_properties
    result.relaxation_properties["best_mapping"]["basis_deformation"] =
      ConfigMapping::basis_cost(best_configdof, result.structure.n_mol());
    result.relaxation_properties["best_mapping"]["lattice_deformation"] =
      ConfigMapping::strain_cost(result.structure.lat_column_mat.determinant(), best_configdof, result.structure.n_mol());
    result.relaxation_properties["best_mapping"]["volume_change"] =
      best_configdof.deformation.determinant();

    // store mapped Configuration
    std::shared_ptr<Supercell> shared_scel = std::make_shared<Supercell>(&primclex(), mapped_lat);

    Configuration import_config(shared_scel, jsonParser(), best_configdof.to_configdof());
    it_canon = import_config.to_canonical();
    result.config = notstd::make_unique<Configuration>(copy_apply(it_canon, import_config));

    // store relaxation_properties
    result.relaxation_properties["best_mapping"]["relaxation_deformation"] =
      it_canon.sym_op().matrix() * best_configdof.deformation * it_canon.sym_op().matrix().transpose();

    // store cart op
    result.cart_op = it_canon.sym_op().matrix() * result.cart_op;
    //store trans
    result.trans = it_canon.sym_op().tau();
    // compose permutations
    std::vector<Index> tperm = it_canon.combined_permute().permute(result.best_assignment);

    //copy non-vacancy part of permutation into best_assignment
    result.best_assignment.resize(result.structure.n_mol());
    Index num_atoms = result.structure.n_mol();
    std::copy_if(tperm.cbegin(), tperm.cend(),
                 result.best_assignment.begin(),
    [num_atoms](Index i) {
      return i < num_atoms;
    });

    result.success = true;
    return result;
  }

  //*******************************************************************************************
  bool ConfigMapper::struc_to_configdof(const SimpleStructure &struc,
                                        MappedConfig &mapped_configdof,
                                        Lattice &mapped_lat) const {
    std::vector<Index> t_assign;
    Eigen::Matrix3d t_op;
    return struc_to_configdof(struc, mapped_configdof, mapped_lat, t_assign, t_op);
  }
  //*******************************************************************************************
  bool ConfigMapper::struc_to_configdof(const SimpleStructure &struc,
                                        MappedConfig &mapped_configdof,
                                        Lattice &mapped_lat,
                                        std::vector<Index> &best_assignment,
                                        Eigen::Matrix3d &cart_op,
                                        double best_cost /*=1e20*/) const {

    bool valid_mapping(false);
    // If structure's lattice is a supercell of the primitive lattice, then import as ideal_structure
    if(!m_robust_flag && Lattice(struc.lat_column_mat).is_supercell_of(primclex().prim().lattice())) {
      valid_mapping = ideal_struc_to_configdof(struc,
                                               mapped_configdof,
                                               mapped_lat,
                                               best_assignment,
                                               cart_op);
      valid_mapping = valid_mapping && ConfigMapping::basis_cost(mapped_configdof, struc.n_mol()) < (10 * m_tol);
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
  bool ConfigMapper::ideal_struc_to_configdof(const SimpleStructure &struc,
                                              MappedConfig &mapped_configdof,
                                              Lattice &mapped_lat,
                                              std::vector<Index> &best_assignment,
                                              Eigen::Matrix3d &cart_op) const {
    // Lattice::is_supercell_of() isn't very smart right now, and will return
    // false if the two lattices differ by a rigid rotation
    // In the future this may not be the case, so we will assume that struc may
    // be rigidly rotated relative to prim
    Eigen::Matrix3d trans_mat;
    if(!Lattice(struc.lat_column_mat).is_supercell_of(primclex().prim().lattice(), trans_mat)) {
      /*std::cerr << "CRITICAL ERROR: In ideal_struc_to_configdof(), primitive structure does not tile the provided\n"
        << "                superstructure. Please use deformed_struc_to_configdof() instead.\n"
        << "                Exiting...\n";
      */
      return false;
    }
    SimpleStructure tstruc(struc);

    // We know struc.lattice() is a supercell of the prim, now we have to
    // reorient 'struc' to match canonical lattice vectors
    mapped_lat = Lattice(primclex().prim().lattice().lat_column_mat() * trans_mat, m_tol).canonical_form(primclex().prim().point_group());
    Supercell scel(&primclex(), mapped_lat);

    // note: trans_mat gets recycled here
    mapped_lat.is_supercell_of(Lattice(struc.lat_column_mat), primclex().prim().point_group(), trans_mat);
    tstruc.lat_column_mat = tstruc.lat_column_mat * trans_mat;

    //cart_op goes from imported coordinate system to ideal (PRIM) coordinate system
    cart_op = mapped_lat.lat_column_mat() * tstruc.lat_column_mat.inverse();
    tstruc.deform(cart_op);
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
  bool ConfigMapper::deformed_struc_to_configdof(const SimpleStructure &struc,
                                                 MappedConfig &mapped_configdof,
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
    double num_atoms = double(struc.n_mol());
    int min_vol, max_vol;
    std::unordered_set<Eigen::Matrix3i, HermiteHash> ref_hermites;

    Eigen::Matrix3i ref_hermite = hermite_normal_form(iround(primclex().prim().lattice().inv_lat_column_mat() *
                                                             niggli(Lattice(struc.lat_column_mat), primclex().crystallography_tol()).lat_column_mat())).first;
    ref_hermites.insert(ref_hermite);
    if(m_restricted) {
      for(auto &g : primclex().prim().point_group()) {
        Eigen::Matrix3i transformed = iround(primclex().prim().lattice().lat_column_mat().inverse() * g.matrix() * primclex().prim().lattice().lat_column_mat()) * ref_hermite;
        Eigen::Matrix3i H_transformed = hermite_normal_form(transformed).first;
        ref_hermites.insert(H_transformed);
      }
    }
    //mapped_configdof.clear();
    if(m_fixed_components.size() > 0) {
      std::string tcompon = m_fixed_components[0].first;
      int ncompon(0);
      for(Index i = 0; i < struc.n_mol(); i++) {
        if(struc.mol_info.names[i] == tcompon)
          ncompon++;
      }
      min_vol = ncompon / int(m_fixed_components[0].second);
      max_vol = min_vol;
    }
    else {
      // Try to narrow the range of supercell volumes -- the best bounds are obtained from
      // the convex hull of the end-members, but we need to wait for improvements to convex hull
      // routines

      int max_n_va = primclex().prim().max_possible_vacancies();
      double max_va_frac_limit = double(max_n_va) / double(primclex().prim().basis().size());
      double t_min_va_frac = min(min_va_frac(), max_va_frac_limit);
      double t_max_va_frac = min(max_va_frac(), max_va_frac_limit);
      // min_vol assumes min number vacancies -- best case scenario
      min_vol = ceil((num_atoms / (double(primclex().prim().basis().size())) * 1. - t_min_va_frac) - m_tol);

      // This is for the worst case scenario -- lots of vacancies
      max_vol = ceil(num_atoms / (double(primclex().prim().basis().size()) * (1.0 - t_max_va_frac)) - m_tol);

      if(t_max_va_frac > TOL) {
        //Nvol is rounded integer volume-- assume that answer is within 30% of this volume, and use it to tighten our bounds
        int Nvol = round(std::abs(struc.lat_column_mat.determinant() / primclex().prim().lattice().vol()));
        int new_min_vol = min(max_vol, max(round((1.0 - m_max_volume_change) * double(Nvol)), min_vol));
        int new_max_vol = max(min_vol, min(round((1.0 + m_max_volume_change) * double(Nvol)), max_vol));
        max_vol = new_max_vol;
        min_vol = new_min_vol;
      }
    }

    min_vol = max(min_vol, 1);
    max_vol = max(max_vol, 1);


    Eigen::Matrix3d ttrans_mat, tF, rotF, U;

    double strain_cost(1e10), basis_cost(1e10), tot_cost;//, best_strain_cost, best_basis_cost;
    MappedConfig tdof;
    SimpleStructure tstruc(struc);
    Lattice tlat;
    //
    // First pass:  Find a reasonable upper bound
    for(Index i_vol = min_vol; i_vol <= max_vol; i_vol++) {

      std::vector<Lattice> lattice_candidates;
      if(m_restricted) {
        lattice_candidates = _lattices_of_vol_restricted(i_vol, ref_hermites);
      }
      else {
        lattice_candidates = _lattices_of_vol(i_vol);
      }
      if(lattice_candidates.size() == 0) {
        //This means you forced lattices on
        continue;
      }
      tlat = ConfigMapping::find_nearest_super_lattice(primclex().prim().lattice(),
                                                       Lattice(struc.lat_column_mat),
                                                       primclex().prim().point_group(),
                                                       tF,
                                                       ttrans_mat,
                                                       lattice_candidates,
                                                       m_tol);
      strain_cost = lw * LatticeMap::calc_strain_cost(tF, struc.lat_column_mat.determinant() / max(num_atoms, 1.));

      if(best_cost < strain_cost) {
        continue;
      }

      tstruc = struc;

      rotF = tF;
      U = StrainConverter::right_stretch_tensor(tF);
      if(m_rotate_flag) {
        rotF = U;
      }

      //make tstruc an un-rotated, un-strained version of struc
      // Math: tF = R*U -> R=tF*U.inv() -> R.inv = U*tF.inv()
      tstruc.deform(U * tF.inverse());
      tstruc.lat_column_mat = tlat.lat_column_mat();

      Supercell scel(&primclex(), tlat);
      if(!ConfigMap_impl::preconditioned_struc_to_configdof(scel,
                                                            tstruc,
                                                            rotF,
                                                            tdof,
                                                            assignment,
                                                            true,
                                                            m_tol))
        continue;
      basis_cost = bw * ConfigMapping::basis_cost(tdof, struc.n_mol());
      tot_cost = strain_cost + basis_cost;
      if(tot_cost < best_cost) {
        best_cost = tot_cost - m_tol;
        swap(best_assignment, assignment);

        cart_op = rotF * tF.inverse();

        std::swap(mapped_configdof, tdof);

        mapped_lat = tlat;
      }
    }

    //Second pass: Find the absolute best mapping
    for(Index i_vol = min_vol; i_vol <= max_vol; i_vol++) {
      std::vector<Lattice> lat_vec;
      if(m_restricted) {
        lat_vec = _lattices_of_vol_restricted(i_vol, ref_hermites);
      }
      else {
        lat_vec = _lattices_of_vol(i_vol);
      }
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

    if(best_cost > 1e9) {
      return false;
    }

    // If mapped_configdof is empty, it means that nothing better than best_cost was found
    return mapped_configdof.occupation.size() > 0;
  }

  //*******************************************************************************************
  bool ConfigMapper::deformed_struc_to_configdof_of_lattice(const SimpleStructure &struc,
                                                            const Lattice &imposed_lat,
                                                            double &best_cost,
                                                            MappedConfig &mapped_configdof,
                                                            Lattice &mapped_lat,
                                                            std::vector<Index> &best_assignment,
                                                            Eigen::Matrix3d &cart_op) const {
    double strain_cost, basis_cost, tot_cost;
    MappedConfig tdof;
    SimpleStructure tstruc(struc);
    Eigen::Matrix3d tF, rotF, U;
    std::vector<Index> assignment;
    Supercell scel(&primclex(), imposed_lat);
    double lw = m_lattice_weight;
    double bw = 1.0 - lw;
    double num_atoms = double(struc.n_mol());
    //Determine best mapping for this supercell

    //Initialize with simplest mapping onto supercell 'i', so that we don't change the crystal setting unnecessarily
    tF = struc.lat_column_mat * imposed_lat.inv_lat_column_mat();

    strain_cost = lw * LatticeMap::calc_strain_cost(tF, struc.lat_column_mat.determinant() / max(num_atoms, 1.));

    // If simplest mapping seems viable, check it further
    if(strain_cost < best_cost) {
      tstruc = struc;
      tstruc.deform(tF.inverse());
      //tstruc.set_lattice(imposed_lat, FRAC);
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

      basis_cost = bw * ConfigMapping::basis_cost(tdof, struc.n_mol());
      tot_cost = strain_cost + basis_cost;

      if(tot_cost < best_cost) {
        best_cost = tot_cost - m_tol;
        //best_strain_cost = strain_cost;
        //best_basis_cost = basis_cost;
        swap(best_assignment, assignment);
        cart_op = rotF * tF.inverse();
        //best_trans = Matrix3<double>::identity();
        std::swap(mapped_configdof, tdof);
        mapped_lat = imposed_lat;
      }
    } // Done checking simplest mapping

    // If the simplest mapping is best, we have avoided a lot of extra work, but we still need to check for
    // non-trivial mappings that are better than both the simplest mapping and the best found mapping
    LatticeMap strainmap(imposed_lat, Lattice(struc.lat_column_mat), round(num_atoms), m_tol, 1);
    strain_cost = lw * strainmap.strain_cost();
    if(best_cost < strain_cost)
      strain_cost = lw * strainmap.next_mapping_better_than(best_cost).strain_cost();

    while(strain_cost < best_cost) {  // only enter loop if there's a chance of improving on current best

      tstruc = struc;
      tF = strainmap.matrixF();
      rotF = tF;
      U = StrainConverter::right_stretch_tensor(tF);
      if(m_rotate_flag) {
        rotF = U;
      }

      //make tstruc an un-rotated, un-strained version of struc
      // We modify the deformed structure so that its lattice is a deformed version of the nearest ideal lattice
      // Math: tF = R*U -> R=tF*U.inv() -> R.inv = U*tF.inv()
      tstruc.deform(U * tF.inverse());
      tstruc.lat_column_mat = imposed_lat.lat_column_mat();

      //tstruc.set_lattice(Lattice(imposed_lat.lat_column_mat()*strainmap.matrixN()), FRAC);
      //tstruc.set_lattice(imposed_lat, CART);

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
      basis_cost = bw * ConfigMapping::basis_cost(tdof, struc.n_mol());

      tot_cost = strain_cost + basis_cost;

      if(tot_cost < best_cost) {
        best_cost = tot_cost - m_tol;
        swap(best_assignment, assignment);
        //best_strain_cost = strain_cost;
        //best_basis_cost = basis_cost;
        cart_op = rotF * tF.inverse();
        //best_trans = strainmap.matrixN();
        std::swap(mapped_configdof, tdof);
        mapped_lat = imposed_lat;
      }
      // This finds first decomposition:
      //            struc.lattice() =  deformation*(supercell_list[import_scel_index].lattice())*equiv_mat
      //   that has cost function less than best_cost
      strain_cost = lw * strainmap.next_mapping_better_than(best_cost).strain_cost();
    }
    return true;
  }

  //*******************************************************************************************
  const std::vector<Lattice> &ConfigMapper::_lattices_of_vol(Index prim_vol) const {

    //If you specified that you wanted certain lattices, return those, otherwise do the
    //usual enumeration
    if(this->lattices_are_forced()) {
      //This may very well return an empty vector, saving painful time enumerating things
      return m_forced_superlat_map[prim_vol];

    }

    if(!valid_index(prim_vol)) {
      throw std::runtime_error("Cannot enumerate lattice of volume " + std::to_string(prim_vol) + ", which is out of bounds.\n");
    }

    //If we already have candidate lattices for the given volume, return those
    auto it = m_superlat_map.find(prim_vol);
    if(it != m_superlat_map.end())
      return it->second;

    //We don't have any lattices for the provided volume, enumerate them all!!!
    std::vector<Lattice> lat_vec;
    SupercellEnumerator<Lattice> enumerator(
      primclex().prim().lattice(),
      primclex().prim().point_group(),
      ScelEnumProps(prim_vol, prim_vol + 1));

    for(auto it = enumerator.begin(); it != enumerator.end(); ++it) {
      Lattice canon_lat = *it;
      if(!canon_lat.is_canonical(primclex().prim().point_group())) {
        canon_lat = canon_lat.canonical_form(primclex().prim().point_group());
      }
      lat_vec.push_back(canon_lat);
    }


    //Save the lattices we enumerated to the map, and return their value
    return m_superlat_map[prim_vol] = std::move(lat_vec);

  }

  //*******************************************************************************************
  const std::vector<Lattice> ConfigMapper::_lattices_of_vol_restricted(Index prim_vol, std::unordered_set<Eigen::Matrix3i, HermiteHash> &ref_hermites) const {
    //If you specified that you wanted certain lattices, return those, otherwise do the
    //usual enumeration
    if(this->lattices_are_forced()) {
      //This may very well return an empty vector, saving painful time enumerating things
      return m_forced_superlat_map[prim_vol];

    }

    if(!valid_index(prim_vol)) {
      throw std::runtime_error("Cannot enumerate lattice of volume " + std::to_string(prim_vol) + ", which is out of bounds.\n");
    }

    //If we already have candidate lattices for the given volume, return those
    auto it = m_superlat_map.find(prim_vol);
    if(it != m_superlat_map.end())
      return it->second;

    //We don't have any lattices for the provided volume, enumerate them all!!!
    std::vector<Lattice> lat_vec;
    SupercellEnumerator<Lattice> enumerator(
      primclex().prim().lattice(),
      primclex().prim().point_group(),
      ScelEnumProps(prim_vol, prim_vol + 1));

    //Save all the lattices we enumerate in their canonical form
    //std::cout << "size of enumerator" << std::distance(enumerator.begin(),enumerator.end());
    //std::cout << " size of equivalent set " << ref_hermites.size() << std::endl;
    for(auto it = enumerator.begin(); it != enumerator.end(); ++it) {

      if(std::none_of(ref_hermites.begin(), ref_hermites.end(), [&](const Eigen::Matrix3i x)->bool{
      return hermite_adjacency(iround(primclex().prim().lattice().inv_lat_column_mat() *
                                      it->lat_column_mat()), x);
      })) {
        continue;
      }
      Lattice canon_lat = *it;
      if(!canon_lat.is_canonical(primclex().prim().point_group())) {
        canon_lat = canon_lat.canonical_form(primclex().prim().point_group());
      }
      lat_vec.push_back(canon_lat);
    }


    //return the lattices within range of the reference hermite normal form of approximate transf_mat
    return lat_vec;

  }

  bool ConfigMapper::hermite_adjacency(const Eigen::Matrix3i test, const Eigen::Matrix3i ref) const {
    auto func = [](int a, int b) {
      return std::abs(a - b) < 3;
    };
    bool all_checks = func(test(0, 0), ref(0, 0)) && func(test(1, 1), ref(1, 1)) && func(test(2, 2), ref(2, 2));
    return all_checks;
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
     * Costs are calculated in context of the lattice.
     */
    //****************************************************************************************************************

    bool calc_cost_matrix(const Supercell &scel,
                          const SimpleStructure &rstruc,
                          const Coordinate &trans,
                          const Eigen::Matrix3d &metric,
                          Eigen::MatrixXd &cost_matrix) {

      if(rstruc.n_mol() > scel.num_sites())
        return false;
      double inf = 10E10;
      //if(cost_matrix.rows()!=scel.num_sites() || cost_matrix.cols()!=scel.num_sites())
      cost_matrix = Eigen::MatrixXd::Constant(scel.num_sites(), scel.num_sites(), inf);
      Index inf_counter;
      // loop through all the sites of the structure
      Index j = 0;

      // cart-to-frac conversion
      Eigen::Matrix3d c2f = rstruc.lat_column_mat.inverse();
      for(; j < rstruc.n_mol(); j++) {
        Coordinate current_relaxed_coord(c2f * rstruc.mol_info.coords.col(j), scel.lattice(), FRAC);
        current_relaxed_coord.cart() += trans.cart();
        // loop through all the sites in the supercell
        inf_counter = 0;
        for(Index i = 0; i < scel.num_sites(); i++) {

          // Check if relaxed atom j is allowed on site i
          // If so, populate cost_matrix normally
          if(scel.prim().basis()[scel.sublat(i)].contains(rstruc.mol_info.names[j])) {
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
          if(scel.prim().basis()[scel.sublat(i)].contains("Va")) {
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
                          const SimpleStructure &rstruc,
                          const Coordinate &trans,
                          const Eigen::Matrix3d &metric,
                          Eigen::MatrixXd &cost_matrix) {


      double inf = 10E10;
      const Supercell &scel(config.supercell());
      //if(cost_matrix.rows()!=scel.num_sites() || cost_matrix.cols()!=scel.num_sites())
      cost_matrix = Eigen::MatrixXd::Constant(scel.num_sites(), scel.num_sites(), inf);
      Index inf_counter;

      // cart-to-frac conversion
      Eigen::Matrix3d c2f = rstruc.lat_column_mat.inverse();

      // loop through all the sites of the structure
      Index j;
      for(j = 0; j < rstruc.n_mol(); j++) {
        Coordinate current_relaxed_coord(c2f * rstruc.mol_info.coords.col(j), scel.lattice(), FRAC);
        current_relaxed_coord.cart() += trans.cart();
        // loop through all the sites in the supercell
        inf_counter = 0;
        for(Index i = 0; i < scel.num_sites(); i++) {

          // Check if relaxed atom j is allowed on site i
          // If so, populate cost_matrix normally
          if(config.mol(i).name() == rstruc.mol_info.names[j]) {
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
          if(config.mol(i).name() == "Va") {
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
                            SimpleStructure rstruc,
                            MappedConfig &config_dof,
                            std::vector<Index> &best_assignments,
                            const bool translate_flag,
                            const double _tol) {
      Eigen::Matrix3d deformation = rstruc.lat_column_mat * scel.lattice().inv_lat_column_mat();
      // un-deform rstruc
      rstruc.deform(deformation.inverse());
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
                            SimpleStructure rstruc,
                            MappedConfig &config_dof,
                            std::vector<Index> &best_assignments,
                            const bool translate_flag,
                            const double _tol) {
      const Lattice &mapped_lat(config.supercell().lattice());
      Eigen::Matrix3d deformation = rstruc.lat_column_mat * mapped_lat.inv_lat_column_mat();
      // un-deform rstruc
      rstruc.deform(deformation.inverse());
      return preconditioned_struc_to_configdof(config, rstruc, deformation, config_dof, best_assignments, translate_flag, _tol);
    }

    //***************************************************************************************************
    // Options:
    //
    // translate_flag = true means that the rigid-translations are removed. (typically this option should be used)
    //
    // translate_flag = false means that rigid translations are not considered. (probably don't want to use this since structures should be considered equal if they are related by a rigid translation).
    //
    // preconditioned_struc_to_configdof is same as struc_to_configdof, except 'rstruc' is de-rotated and de-strained.
    // Any deformation is instead specified by 'deformation'
    bool preconditioned_struc_to_configdof(const Supercell &scel,
                                           const SimpleStructure &rstruc,
                                           const Eigen::Matrix3d &deformation,
                                           MappedConfig &config_dof,
                                           std::vector<Index> &best_assignments,
                                           const bool translate_flag,
                                           const double _tol) {
      // clear config_dof and set its deformation
      //config_dof.clear();

      config_dof.deformation = deformation;
      Eigen::Matrix3d metric(deformation.transpose()*deformation);
      //Initialize everything

      Eigen::MatrixXd cost_matrix;
      std::vector<Index> optimal_assignments;
      //BasicStructure<Site> best_ideal_struc(rstruc);
      Coordinate best_trans(scel.lattice());

      double min_mean = 10E10;

      // We want to get rid of translations.
      // define translation such that:
      //    IDEAL = RELAXED + translation
      // and use it when calculating cost matrix

      Index num_translations(1);

      if(rstruc.n_mol())
        num_translations += scel.prim().basis().size();

      for(Index n = 0; n < num_translations; n++) {
        double mean;

        //shift_struc has **ideal lattice**
        //BasicStructure<Site> shift_struc(rstruc);


        if(n > 0 && !scel.prim().basis()[n - 1].contains(rstruc.mol_info.names[0]))
          continue;

        Coordinate translation(scel.prim().lattice());

        // Always try the non-translated case (n==0), in case it gives best result
        // Also try translating first basis atom onto each chemically compatible site of PRIM (n>0)
        if(n > 0) {
          translation.cart() = scel.coord((n - 1) * scel.volume()).const_cart() - rstruc.mol_info.coords.col(0);
          translation.voronoi_within();
        }

        if(!ConfigMap_impl::calc_cost_matrix(scel, rstruc, translation, metric, cost_matrix)) {
          /// Indicates that structure is incompatible with supercell, so return false
          return false;
        }

        // The mapping routine is called here
        mean = hungarian_method(cost_matrix, optimal_assignments, _tol);

        // if optimal_assignments is smaller than rstruc.basis().size(), then rstruc is incompattible with supercell
        // (optimal_assignments.size()==0 if the hungarian routine detects an incompatibility)
        if(optimal_assignments.size() < rstruc.n_mol()) {
          return false;
        }


        // add small penalty (~_tol) for larger translation distances, so that shortest equivalent translation is used
        mean += _tol * translation.const_cart().norm() / 10.0;

        if(mean < min_mean) {

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
      config_dof.displacement.setZero(3, scel.num_sites());

      Eigen::Vector3d avg_disp(0, 0, 0);

      Coordinate disp_coord(scel.lattice());

      // Populate displacements given as the difference in the Coordinates
      // as described by best_assignments.
      for(Index i = 0; i < best_assignments.size(); i++) {

        // If we are dealing with a vacancy, its displacment must be zero.
        //if(best_assignments(i) >= rstruc.n_mol()) {
        //  --DO NOTHING--
        //}


        // Using min_dist routine to calculate the displacement vector that corresponds
        // to the distance used in the Cost Matrix and Hungarian Algorithm
        // The method returns the displacement vector pointing from the
        // IDEAL coordinate to the RELAXED coordinate
        if(best_assignments[i] < rstruc.n_mol()) {

          Coordinate relaxed_coord(rstruc.mol_info.coords.col(best_assignments[i]) + best_trans.const_cart(), scel.lattice(), CART);

          relaxed_coord.min_dist(scel.coord(i), disp_coord);
          config_dof.disp(i) = disp_coord.const_cart();

          avg_disp += config_dof.disp(i);
        }
      }

      avg_disp /= max(double(rstruc.n_mol()), 1.);


      // End of filling displacements


      // Make the assignment bitstring
      //
      // Loop through all supercell sites
      config_dof.occupation.setZero(scel.num_sites());
      std::string rel_basis_atom;
      for(Index i = 0; i < best_assignments.size(); i++) {
        // subtract off average displacement
        if(best_assignments[i] < rstruc.n_mol()) {
          config_dof.disp(i) -= avg_disp;
          // suppress small ugly numbers.
          for(Index j = 0; j < 3; j++) {
            if(almost_zero(config_dof.disp(i)[j], 1e-8))
              config_dof.disp(i)[j] = 0;
          }

          //Record basis atom
          rel_basis_atom = rstruc.mol_info.names[best_assignments[i]];
        }
        else {
          // Any value of the assignment vector larger than the number
          // of sites in the relaxed structure is by construction
          // specified as a vacancy.
          rel_basis_atom = "Va";
        }

        // set occupant and check for errors
        if(!scel.prim().basis()[scel.sublat(i)].contains(rel_basis_atom, config_dof.occupation[i])) {

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
                                           const SimpleStructure &rstruc,
                                           const Eigen::Matrix3d &deformation,
                                           MappedConfig &config_dof,
                                           std::vector<Index> &best_assignments,
                                           const bool translate_flag,
                                           const double _tol) {

      const Supercell &scel(config.supercell());

      // clear config_dof and set its deformation
      //config_dof.clear();

      config_dof.deformation = deformation;
      Eigen::Matrix3d metric(deformation.transpose()*deformation);
      //Initialize everything

      Eigen::MatrixXd cost_matrix;
      std::vector<Index> optimal_assignments;
      //BasicStructure<Site> best_ideal_struc(rstruc);
      Coordinate best_trans(scel.lattice());

      double min_mean = 10E10;

      // We want to get rid of translations.
      // define translation such that:
      //    IDEAL = RELAXED + translation
      // and use it when calculating cost matrix

      Index num_translations(1);

      num_translations += rstruc.n_mol();

      //num_translations = rstruc.n_mol();
      for(Index n = 0; n < num_translations; n++) {
        double mean;

        if(n > 0 && config.mol(0).name() != rstruc.mol_info.names[n - 1])
          continue;

        Coordinate translation(scel.lattice());

        // Always try the non-translated case (n==0), in case it gives best result
        // Also try translating first basis atom onto each chemically compatible site of PRIM (n>0)
        if(n > 0) {
          translation.cart() = scel.coord(0).const_cart() - rstruc.mol_info.coords.col(n - 1);
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

        // if optimal_assignments is smaller than rstruc.n_mol(), then rstruc is incompattible
        // with the supercell (optimal_assignments.size()==0 if the hungarian routine detects an incompatibility)
        if(optimal_assignments.size() < rstruc.n_mol()) {
          return false;
        }

        // add small penalty (~_tol) for larger translation distances, so that shortest equivalent translation is used
        mean += _tol * translation.const_cart().norm() / 10.0;
        if(mean < min_mean) {

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
      config_dof.displacement.setZero(3, scel.num_sites());

      Eigen::Vector3d avg_disp(0, 0, 0);

      Coordinate disp_coord(scel.lattice());

      // Populate displacements given as the difference in the Coordinates
      // as described by best_assignments.
      for(Index i = 0; i < best_assignments.size(); i++) {

        // If we are dealing with a vacancy, its displacment must be zero.
        //if(best_assignments(i) >= rstruc.n_mol()) {
        //  --DO NOTHING--
        //}


        // Using min_dist routine to calculate the displacement vector that corresponds
        // to the distance used in the Cost Matrix and Hungarian Algorithm
        // The method returns the displacement vector pointing from the
        // IDEAL coordinate to the RELAXED coordinate
        if(best_assignments[i] < rstruc.n_mol()) {

          Coordinate relaxed_coord(rstruc.mol_info.coords.col(best_assignments[i]) + best_trans.const_cart(), scel.lattice(), CART);

          relaxed_coord.min_dist(scel.coord(i), disp_coord);

          config_dof.disp(i) = disp_coord.const_cart();

          avg_disp += config_dof.disp(i);
        }
      }

      avg_disp /= max(double(rstruc.n_mol()), 1.);


      // End of filling displacements


      // Make the assignment bitstring
      //
      // Loop through all supercell sites
      config_dof.occupation.setZero(scel.num_sites());
      std::string rel_basis_atom;
      for(Index i = 0; i < best_assignments.size(); i++) {
        // subtract off average displacement (non-vacant sites only)
        // Any value of the assignment vector larger than the number
        // of sites in the relaxed structure is by construction
        // specified as a vacancy.
        if(best_assignments[i] < rstruc.n_mol()) {
          config_dof.disp(i) -= avg_disp;
          // suppress small ugly numbers.
          for(Index j = 0; j < 3; j++) {
            if(almost_zero(config_dof.disp(i)[j], 1e-8))
              config_dof.disp(i)[j] = 0;
          }
        }
        config_dof.occupation[i] = config.occ(i);
      }

      return true;
    }

  }
}
