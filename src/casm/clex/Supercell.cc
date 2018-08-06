#include "casm/clex/Supercell.hh"

#include <math.h>
#include <map>
#include <vector>
#include <stdlib.h>

//#include "casm/clusterography/HopCluster.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/ConfigIterator.hh"
#include "casm/clex/Clexulator.hh"

namespace CASM {

  //Given a Site and tolerance, return linear index into Configuration
  //   This may be slow, first converts Site -> UnitCellCoord,
  //   then finds UnitCellCoord in config_index_to_bijk
  Index Supercell::get_linear_index(const Site &site, double tol) const {
    //std::cout << "site: " << site << "  UCC: " << get_prim().get_unit_cell_coord(site, tol) << std::endl;
    Site tsite(site);
    tsite.within();
    return find(get_prim().get_unit_cell_coord(tsite, tol));
  };

  /*****************************************************************/

  //Given a Coordinate and tolerance, return linear index into Configuration
  //   This may be slow, first converts Coordinate -> UnitCellCoord,
  //   then finds UnitCellCoord in config_index_to_bijk
  Index Supercell::get_linear_index(const Coordinate &coord, double tol) const {
    //std::cout << "coord: " << coord << "  UCC: " << get_prim().get_unit_cell_coord(coord, tol) << std::endl;
    Coordinate tcoord(coord);
    tcoord.within();
    return find(get_prim().get_unit_cell_coord(tcoord, tol));
  };

  /*****************************************************************/

  Index Supercell::find(const UnitCellCoord &bijk) const {
    return bijk[0] * volume() + m_prim_grid.find(bijk);
  }

  /*****************************************************************/

  Coordinate Supercell::coord(const UnitCellCoord &bijk) const {
    Coordinate tcoord(m_prim_grid.coord(bijk, SCEL));
    tcoord.cart() += (*primclex).get_prim().basis[bijk[0]].cart();
    return tcoord;
  };

  /*****************************************************************/

  Coordinate Supercell::coord(Index l) const {
    Coordinate tcoord(m_prim_grid.coord(l % volume(), SCEL));
    tcoord.cart() += (*primclex).get_prim().basis[get_b(l)].cart();
    return tcoord;
  };

  /*****************************************************************/

  ReturnArray<int> Supercell::max_allowed_occupation() const {
    Array<int> max_allowed;

    // Figures out the maximum number of occupants in each basis site, to initialize counter with
    for(Index i = 0; i < get_prim().basis.size(); i++) {
      max_allowed.append(Array<int>(volume(), get_prim().basis[i].site_occupant().size() - 1));
    }
    //std::cout << "max_allowed_occupation is:  " << max_allowed << "\n\n";
    return max_allowed;
  }

  /*****************************************************************/

  const Structure &Supercell::get_prim() const {
    return primclex->get_prim();
  }

  /// \brief Returns the SuperNeighborList
  const SuperNeighborList &Supercell::nlist() const {

    // if any additions to the prim nlist, must update the super nlist
    if(get_primclex().nlist().size() != m_nlist_size_at_construction) {
      m_nlist.unique().reset();
    }

    // lazy construction of neighbor list
    if(!m_nlist) {
      m_nlist_size_at_construction = get_primclex().nlist().size();
      m_nlist = notstd::make_cloneable<SuperNeighborList>(
                  m_prim_grid,
                  get_primclex().nlist()
                );
    }
    return *m_nlist;
  };

  /*****************************************************************/

  // begin and end iterators for iterating over configurations
  Supercell::config_iterator Supercell::config_begin() {
    return config_iterator(primclex, m_id, 0);
  }

  Supercell::config_iterator Supercell::config_end() {
    return ++config_iterator(primclex, m_id, config_list.size() - 1);
  }

  // begin and end const_iterators for iterating over configurations
  Supercell::config_const_iterator Supercell::config_cbegin() const {
    return config_const_iterator(primclex, m_id, 0);
  }

  Supercell::config_const_iterator Supercell::config_cend() const {
    return ++config_const_iterator(primclex, m_id, config_list.size() - 1);
  }

  /// \brief Return supercell name
  ///
  /// - If lattice is the canonical equivalent, then return 'SCELV_A_B_C_D_E_F'
  /// - Else, return 'SCELV_A_B_C_D_E_F.$FG_INDEX', where $FG_INDEX is the index of the first
  ///   symmetry operation in the primitive structure's factor group such that the lattice
  ///   is equivalent to `apply(fg_op, canonical equivalent)`
  std::string Supercell::get_name() const {
    if(m_name.empty()) {
      _generate_name();
    }
    return m_name;
  };

  /*****************************************************************/

  const SymGroup &Supercell::factor_group() const {
    if(!m_factor_group.size())
      generate_factor_group();
    return m_factor_group;
  }

  /*****************************************************************/

  // permutation_symrep() populates permutation symrep if needed
  const Permutation &Supercell::factor_group_permute(Index i) const {
    return *(permutation_symrep().get_permutation(factor_group()[i]));
  }
  /*****************************************************************/

  // PrimGrid populates translation permutations if needed
  const Permutation &Supercell::translation_permute(Index i) const {
    return m_prim_grid.translation_permutation(i);
  }

  /*****************************************************************/

  // PrimGrid populates translation permutations if needed
  const Array<Permutation> &Supercell::translation_permute() const {
    return m_prim_grid.translation_permutations();
  }

  /*****************************************************************/

  /// \brief Begin iterator over translation permutations
  Supercell::permute_const_iterator Supercell::translate_begin() const {
    return permute_begin();
  }

  /*****************************************************************/

  /// \brief End iterator over translation permutations
  Supercell::permute_const_iterator Supercell::translate_end() const {
    return permute_begin().begin_next_fg_op();
  }

  /*****************************************************************/
  /* //Example usage case:
   *  Supercell my_supercell;
   *  Configuration my_config(my_supercell, configuration_info);
   *  ConfigDoF my_dof=my_config.configdof();
   *  my_dof.is_canonical(my_supercell.permute_begin(),my_supercell.permute_end());
   */
  Supercell::permute_const_iterator Supercell::permute_begin() const {
    return permute_it(0, 0); // starting indices
  }

  /*****************************************************************/

  Supercell::permute_const_iterator Supercell::permute_end() const {
    return permute_it(factor_group().size(), 0); // one past final indices
  }

  /*****************************************************************/

  Supercell::permute_const_iterator Supercell::permute_it(Index fg_index, Index trans_index) const {
    return permute_const_iterator(SymGroupRep::RemoteHandle(factor_group(), permutation_symrep_ID()),
                                  m_prim_grid,
                                  fg_index, trans_index); // one past final indices
  }

  /*****************************************************************/

  //Printing config_index_to_bijk
  void Supercell::print_bijk(std::ostream &stream) {
    for(Index i = 0; i < num_sites(); i++) {
      stream << uccoord(i);
    }
  }

  //*******************************************************************************
  /**
   *   enumerate_perturb_configurations, using filename of 'background' structure
   */
  //*******************************************************************************
  void Supercell::enumerate_perturb_configurations(const std::string &background, fs::path CSPECS, double tol, bool verbose, bool print) {
    Structure background_struc;
    fs::ifstream file(background);
    background_struc.read(file);
    enumerate_perturb_configurations(background_struc, CSPECS, tol, verbose, print);
  }

  //*******************************************************************************
  /**
   *   enumerate_perturb_configurations, using 'config' Configuration and 'CSPECS'
   *     to generate the 'background_config' and 'background_tree'.
   *     The factor group of the decorated config is used to generate the orbitree
   */
  //*******************************************************************************
  void Supercell::enumerate_perturb_configurations(Configuration background_config, fs::path CSPECS, double tol, bool verbose, bool print) {
    // Algorithm:
    // 1) generate orbitree in background
    // 2) generate background config
    // 3) for each orbit:
    //      perturb background config with decorated prototype cluster
    //        check if in config list
    // NOTE: This can be done much faster using permutation arithmetic
    if(verbose)   std::cout << "begin enumerate_perturb_configurations" << std::endl;

    // should generate the background_tree from the supercell-sized background structure
    //   this gets the right symmetry for the combination of perturbation and supercell shape

    if(verbose)   std::cout << "Generate background structure" << std::endl;
    Structure background_scel = superstructure(background_config);

    // generate the background config & orbitree
    //   std::cout << "generate background config and orbitree" << std::endl;
    SiteOrbitree background_tree(background_scel.lattice(), primclex->crystallography_tol());

    //fs::ifstream cspecsfile(CSPECS);
    //background_tree.read_CSPECS(cspecsfile);
    //cspecsfile.close();

    jsonParser json(CSPECS);
    background_tree.min_num_components = 2;
    background_tree.min_length = CASM::TOL;

    background_tree.max_length.clear();
    auto update_max_length = [&](int branch, double max_length) {
      while(branch > background_tree.max_length.size() - 1) {
        background_tree.max_length.push_back(0.0);
      }
      background_tree.max_length[branch] = max_length;
    };

    for(auto it = json["orbit_branch_specs"].cbegin(); it != json["orbit_branch_specs"].cend(); ++it) {
      update_max_length(std::stoi(it.name()), it->find("max_length")->get<double>());
    }
    background_tree.max_num_sites = background_tree.max_length.size() - 1;



    if(verbose)   std::cout << "Generate background orbitree" << std::endl;
    background_tree.generate_orbitree(background_scel);

    if(verbose) std::cout << "background_config: " << background_config.name() << std::endl;

    // for now, don't do anything with these here
    Array< Array< Array<Index> > > perturb_config_index;
    Array< Array< Array<permute_const_iterator> > > perturb_config_symop_index;

    if(verbose)   std::cout << "Enumerate perturb configurations" << std::endl;

    jsonParser jsonsrc = jsonParser::object();
    jsonsrc["supercell_name"] = get_name();
    jsonsrc["configid"] = background_config.get_id();

    enumerate_perturb_configurations(background_config, background_tree, perturb_config_index, perturb_config_symop_index, jsonsrc, tol);

    if(verbose) {
      for(Index nb = 0; nb < perturb_config_index.size(); nb++) {
        std::cout << "    Branch: " << nb << std::endl;

        for(Index no = 0; no < perturb_config_index[nb].size(); no++) {
          std::cout << "      Orbit: " << no << std::endl;
          background_tree.prototype(nb, no).print_decorated_sites(std::cout, 8, '\n');

          for(Index nd = 0; nd < perturb_config_index[nb][no].size(); nd++) {
            std::cout << "        config_index: " << perturb_config_index[nb][no][nd] << std::endl;
          }
        }
      }
    }

    if(print) {
      if(verbose)   std::cout << "Print info" << std::endl;

      // write in supercells/scel_name/config_name.perturb
      try {
        fs::create_directory("training_data");
      }
      catch(const fs::filesystem_error &ex) {
        std::cerr << "Error in Supercell::enumerate_perturb_configurations()." << std::endl;
        std::cerr << ex.what() << std::endl;
      }

      try {
        fs::create_directory(get_path());
      }
      catch(const fs::filesystem_error &ex) {
        std::cerr << "Error in Supercell::enumerate_perturb_configurations()." << std::endl;
        std::cerr << ex.what() << std::endl;
      }

      //const fs::path config_path = background_config.get_path() += ".perturb";    boost version clash;
      std::string pathstr = background_config.get_path().filename().string() + ".perturb";      //make string for only the filename (myfile.perturb)
      const fs::path config_path = background_config.get_path().remove_filename() /= pathstr;   //remve myfile from path and add myfile.perturb instead


      try {
        fs::create_directory(config_path);

        // write CSPECS, FCLUST, PERTURB.json
        //  - overwrite if necessary

        // write CSPECS
        {
          if(fs::exists(config_path / "CSPECS"))
            fs::remove(config_path / "CSPECS");
          fs::copy(fs::path(CSPECS), config_path / "CSPECS");
        }

        // write CLUST
        {
          if(fs::exists(config_path / "CLUST"))
            fs::remove(config_path / "CLUST");

          background_tree.write_proto_clust((config_path / "CLUST").string());

        }

        // write FCLUST
        {
          if(fs::exists(config_path / "FCLUST"))
            fs::remove(config_path / "FCLUST");

          background_tree.write_full_clust((config_path / "FCLUST").string());

        }

        // write PERTURB.json
        {
          if(fs::exists(config_path / "PERTURB.json"))
            fs::remove(config_path / "PERTURB.json");

          fs::ofstream file(config_path / "PERTURB.json");

          print_PERTURB_json(file, background_config, perturb_config_index, perturb_config_symop_index, false);
        }

      }
      catch(const fs::filesystem_error &ex) {
        std::cerr << "Error in Supercell::enumerate_perturb_configurations()." << std::endl;
        std::cerr << ex.what() << std::endl;
      }
    }

    if(verbose)   std::cout << "finish enumerate_perturb_configurations" << std::endl;

  }

  //*******************************************************************************
  /**
   *   enumerate_perturb_configurations, using 'background' Structure and 'CSPECS'
   *     to generate the 'background_config' and 'background_tree'.
   */
  //*******************************************************************************
  void Supercell::enumerate_perturb_configurations(const Structure &background, fs::path CSPECS, double tol, bool verbose, bool print) {

    Configuration background_config = configuration(background);
    enumerate_perturb_configurations(background_config, CSPECS, tol, verbose, print);

  };

  //*******************************************************************************
  /**
   *   Enumerate configurations that are perturbations of a 'background_config'.
   *     The 'perturbed' configurations differ from the 'background' structure by
   *     clusters in the 'background_tree'.  'tol' provides a tolerance for mapping
   *     the clusters to Configuration sites.
   *
   *   Enumerated configurations are added to 'Supercell::config_list' if they
   *     do not already exist there, using the 'permute_group' to check for equivalents.
   *
   *   Array< Array< Array<int> > > config_indices contains the mapping of [branch][orbit][decor] to config_list index
   *   Array< Array< Array<int> > > config_symop contains the index of the symop which mapped the config to canonical form
   *
   *   jsonsrc is a jsonParser (object type) describing the source of the enumerate configurations
   *
   */
  //*******************************************************************************
  void Supercell::enumerate_perturb_configurations(Configuration background_config,
                                                   const SiteOrbitree &background_tree,
                                                   Array< Array< Array<Index> > > &config_index,
                                                   Array< Array< Array<permute_const_iterator> > > &config_symop_index,
                                                   jsonParser &jsonsrc,
                                                   double tol) {

    //std::cout << "begin enumerate_perturb_configurations() ****" << std::endl;

    /* primitive pointer no longer exists
    if((*primclex).get_prim().lattice.primitive != background_tree.lattice.primitive) {
      std::cerr << "Error in Supercell::enumerate_perturb_configurations." << std::endl;
      std::cerr << "  'background_tree' lattice primitive is not primclex->prim lattice primitive" << std::endl;
      exit(1);
    }
    */

    // perturb_list.json: json["background_id"]["branch"]["orbit"][{decor = [linear_indices], id = config_id}]


    Configuration config = background_config;
    config.set_selected(false);

    // variables used for generating perturb configs
    Array< Array<int> > decor_map;
    Array<int> linear_indices;
    Array<int> orig_occ;
    Index index;
    permute_const_iterator permute_it;

    config_index.resize(background_tree.size());
    config_symop_index.resize(background_tree.size());


    // for each branch in 'background_tree'
    //std::cout << "loop over background_tree" << std::endl;
    for(Index nb = 0; nb < background_tree.size(); nb++) {

      //std::cout << "branch " << nb << std::endl;
      config_index[nb].resize(background_tree[nb].size());
      config_symop_index[nb].resize(background_tree[nb].size());

      // for each orbit
      for(Index no = 0; no < background_tree[nb].size(); no++) {
        //std::cout << "\n\n---------------------" << std::endl;
        //std::cout << "branch: " << nb << "  orbit " << no << std::endl;


        // get decor_map for prototype
        //std::cout << "get decor_map" << std::endl;
        decor_map = background_tree[nb][no].prototype.get_full_decor_map();

        // determine linear_index for cluster sites
        //std::cout << "get linear_indices and orig_occ" << std::endl;
        linear_indices.clear();
        orig_occ.clear();
        for(Index i = 0; i < background_tree[nb][no].prototype.size(); i++) {
          //std::cout << "  Site: " << i << "  :: " << background_tree[nb][no].prototype[i] << std::endl;
          linear_indices.push_back(get_linear_index(Coordinate(background_tree[nb][no].prototype[i]), tol));
          //std::cout << "    linear_index: " << linear_indices.back() << std::endl;
          //std::cout << "    frac_coord: " << frac_coord( config_index_to_bijk[linear_indices.back()]) << std::endl;
          orig_occ.push_back(config.occ(linear_indices[i]));

        }

        //Generate new clusters with different decorations using decor_map
        //std::cout << "decorate" << std::endl;
        for(Index i = 0; i < decor_map.size(); i++) {
          //std::cout << "decor_map " << i << ": " << decor_map[i] <<  std::endl;
          // set occupants
          for(Index j = 0; j < decor_map[i].size(); j++) {
            config.set_occ(linear_indices[j], decor_map[i][j]);
          }

          // At this point, 'config' is the perturbed config (using prototype & decor_map[i])

          jsonsrc["perturbation"].put_obj();
          jsonsrc["perturbation"]["branch"] = nb;
          jsonsrc["perturbation"]["orbit"] = no;
          jsonsrc["perturbation"]["decor"] = decor_map[i];

          config.set_source(jsonsrc);
          add_config(config, index, permute_it);
          //std::cout << "add_config: " << index << "  result: " << result << "  nb: " << nb << "  no: " << no << std::endl;


          config_index[nb][no].push_back(index);
          config_symop_index[nb][no].push_back(permute_it);

          //std::cout << "next" << std::endl << std::endl;
        }

        // reset 'config' to original occupants
        //std::cout << "reset background" << std::endl;
        for(Index i = 0; i < orig_occ.size(); i++) {
          config.set_occ(linear_indices[i], orig_occ[i]);
        }

        //std::cout << "next orbit" << std::endl;

      }
    }

    //std::cout << "finish enumerate_perturb_configurations() ****" << std::endl;

  };

  //*******************************************************************************
  /**
   *   Checks if the Configuration 'config' is contained in Supercell::config_list.
   *     Only checks Configuration::occupation for equivalence.
   *     Does not check for symmetrically equivalent Configurations, so put your
   *     'config' in canonical form first.
   */
  //*******************************************************************************
  bool Supercell::contains_config(const Configuration &config) const {
    Index index;
    return contains_config(config, index);
  };

  //*******************************************************************************
  /**
   *   Checks if the Configuration 'config' is contained in Supercell::config_list.
   *     Only checks Configuration::configdof for equivalence.
   *     Does not check for symmetrically equivalent Configurations, so put your
   *     'config' in canonical form first.
   *
   *   If equivalent found, 'index' contains it's index into config_list, else
   *     'index' = config_list.size().
   */
  //*******************************************************************************
  bool Supercell::contains_config(const Configuration &config, Index &index) const {
    auto res = m_config_map.find(&config);
    if(res == m_config_map.end()) {
      index = config_list.size();
      return false;
    }
    index = res->second;
    return true;
  };

  //*******************************************************************************
  Supercell::config_const_iterator Supercell::find(const Configuration &config) const {
    auto res = m_config_map.find(&config);
    if(res == m_config_map.end()) {
      return config_cend();
    }
    return config_const_iterator(&get_primclex(), get_id(), res->second);
  }

  //*******************************************************************************
  /**
   *   Converts 'config' to canonical form, then adds to config_list if not already
   *     present. Location in config_list is stored in 'index'.
   *     Permutation that resulted in canonical form is stored in 'permute_it'.
   *     Return 'true' if new config, 'false' otherwise.
   *
   *   Might want to rewrite without using new canon_config for memory/speed issues.
   */
  //*******************************************************************************

  bool Supercell::add_config(const Configuration &config) {
    Index index;
    Supercell::permute_const_iterator permute_it;
    return add_config(config, index, permute_it);
  }

  bool Supercell::add_config(const Configuration &config, Index &index, Supercell::permute_const_iterator &permute_it) {
    // 'permut_it' stores operation that takes 'config' to its canonical form
    permute_it = config.to_canonical();

    // std::cout << "    config: " << config.occupation() << std::endl;
    // std::cout << "     canon: " << canon_config.occupation() << std::endl;

    return add_canon_config(copy_apply(permute_it, config), index);
  }

  //*******************************************************************************
  /**
   *   Assumes 'canon_config' is in canonical form, adds to config_list if not already there.
   *     Location in config_list is stored in 'index'.
   */
  //*******************************************************************************
  bool Supercell::add_canon_config(const Configuration &canon_config, Index &index) {

    // Add 'canon_config' to 'config_list' if it doesn't already exist
    //   store it's index into 'config_list' in 'config_list_index'
    //std::cout << "check if canon_config is in config_list" << std::endl;
    if(!contains_config(canon_config, index)) {
      //std::cout << "new config" << std::endl;
      _add_canon_config(canon_config);
      return true;
      //std::cout << "    added" << std::endl;
    }
    else {
      config_list[index].push_back_source(canon_config.source());
    }
    return false;
  }

  //*******************************************************************************

  /// \brief Insert a configuration that may be non-canonical
  std::pair<Supercell::config_const_iterator, bool>
  Supercell::insert_config(const Configuration &config) {
    return insert_canon_config(config.canonical_form());
  }

  //*******************************************************************************

  /// \brief Insert a configuration that is known to be canonical
  std::pair<Supercell::config_const_iterator, bool>
  Supercell::insert_canon_config(const Configuration &canon_config) {
    Index index;
    std::pair<Supercell::config_const_iterator, bool> res;
    res.second = false;
    if(!contains_config(canon_config, index)) {
      _add_canon_config(canon_config);
      index = config_list.size() - 1;
      res.second = true;
    }
    res.first = config_const_iterator(&get_primclex(), get_id(), index);
    return res;
  }

  //*******************************************************************************
  /**
   *   Adds to config_list, assuming 'canon_config' is in canonical form and not
   *   already there
   */
  //*******************************************************************************
  void Supercell::_add_canon_config(const Configuration &canon_config) {

    if(this != &canon_config.get_supercell()) {
      throw std::runtime_error("Error adding Configuration to Supercell: Supercell mismatch");
    }
    //std::cout << "new config" << std::endl;
    config_list.push_back(canon_config);
    config_list.back().set_id(config_list.size() - 1);
    m_config_map.insert(
      std::make_pair(&config_list.back(),
                     boost::lexical_cast<Index>(config_list.back().get_id())));
    config_list.back().set_selected(false);
  }

  //*******************************************************************************

  void Supercell::read_config_list(const jsonParser &json) {

    // Provide an error check
    if(config_list.size() != 0) {
      std::cerr << "Error in Supercell::read_configuration." << std::endl;
      std::cerr << "  config_list.size() != 0, only use this once" << std::endl;
      exit(1);
    }

    if(!json.contains("supercells")) {
      return;
    }

    if(!json["supercells"].contains(get_name())) {
      return;
    }

    // Read all configurations for this supercell. They should be numbered sequentially, so read until not found.
    Index configid = 0;
    while(true) {
      std::stringstream ss;
      ss << configid;

      if(json["supercells"][get_name()].contains(ss.str())) {
        config_list.push_back(Configuration(json, *this, configid));
        m_config_map.insert(
          std::make_pair(&config_list.back(),
                         boost::lexical_cast<Index>(config_list.back().get_id())));
      }
      else {
        return;
      }
      configid++;
    }
  }


  //*******************************************************************************

  //Copy constructor is needed for proper initialization of m_prim_grid
  Supercell::Supercell(const Supercell &RHS) :
    primclex(RHS.primclex),
    real_super_lattice(RHS.real_super_lattice),
    recip_prim_lattice(RHS.recip_prim_lattice),
    m_prim_grid((*primclex).get_prim().lattice(), real_super_lattice, (*primclex).get_prim().basis.size()),
    recip_grid(recip_prim_lattice, (*primclex).get_prim().lattice().get_reciprocal()),
    m_name(RHS.m_name),
    m_nlist(RHS.m_nlist),
    m_canonical(nullptr),
    config_list(RHS.config_list),
    transf_mat(RHS.transf_mat),
    scaling(RHS.scaling),
    m_id(RHS.m_id) {
  }

  //*******************************************************************************

  Supercell::Supercell(PrimClex *_prim, const Eigen::Ref<const Eigen::Matrix3i> &transf_mat_init) :
    primclex(_prim),
    real_super_lattice((*primclex).get_prim().lattice().lat_column_mat() * transf_mat_init.cast<double>()),
    recip_prim_lattice(real_super_lattice.get_reciprocal()),
    m_prim_grid((*primclex).get_prim().lattice(), real_super_lattice, (*primclex).get_prim().basis.size()),
    recip_grid(recip_prim_lattice, (*primclex).get_prim().lattice().get_reciprocal()),
    m_canonical(nullptr),
    transf_mat(transf_mat_init) {
    scaling = 1.0;
    //    fill_reciprocal_supercell();
  }

  //*******************************************************************************

  Supercell::Supercell(PrimClex *_prim, const Lattice &superlattice) :
    primclex(_prim),
    //real_super_lattice((get_prim()).lattice().lat_column_mat()*transf_mat),
    real_super_lattice(superlattice),
    recip_prim_lattice(real_super_lattice.get_reciprocal()),
    m_prim_grid((*primclex).get_prim().lattice(), real_super_lattice, (*primclex).get_prim().basis.size()),
    recip_grid(recip_prim_lattice, (*primclex).get_prim().lattice().get_reciprocal()),
    m_canonical(nullptr),
    transf_mat(primclex->calc_transf_mat(superlattice)) {
    /*std::cerr << "IN SUPERCELL CONSTRUCTOR:\n"
              << "transf_mat is\n" << transf_mat << '\n'
              << "prim lattice is \n";
    (*primclex).get_prim().lattice.print(std::cerr);
    std::cerr << "\nSupercell Lattice is\n";
    real_super_lattice.print(std::cerr);
    std::cerr << "\nORIGINAL Supercell Lattice is\n";
    superlattice.print(std::cerr);

    std::cerr << "\nlat_column_mat() is\n" << (*primclex).get_prim().lattice.lat_column_mat()
              << "\n and product with transf_mat is \n" << (*primclex).get_prim().lattice.lat_column_mat()*transf_mat << "\n";
    */
    scaling = 1.0;

  }

  //*******************************************************************************
  /**
   * Run through every selected Configuration in *this and call write() on it. This will
   * update all the JSON files and also rewrite POS, DoF etc. Meant for when
   * you calculated some properties (e.g. formation energies or correlations) and
   * want it outputted, but didn't generate any new configurations.
   */

  jsonParser &Supercell::write_config_list(jsonParser &json) {
    for(Index c = 0; c < config_list.size(); c++) {
      config_list[c].write(json);
    }
    return json;
  }


  //*******************************************************************************
  /**
   *   Print the PERTURB file for perturbations enumerated around a Configuration
   *
   *   If 'print_config_name' == true, print config using Configuration::name()
   *   If 'print_config_name' == false, print config using config_index
   *
   */
  //*******************************************************************************
  void Supercell::print_PERTURB_json(std::ofstream &file,
                                     const Configuration &background_config,
                                     const Array< Array< Array<Index > > > &perturb_config_index,
                                     const Array< Array< Array<permute_const_iterator> > > &perturb_config_symop_index,
                                     bool print_config_name) const {

    jsonParser json = jsonParser::object();

    json["supercell_name"] = get_name();
    if(print_config_name) {
      json["config"] = background_config.name();
    }
    else {
      json["configid"] = background_config.get_id();
    }
    json["perturbations"] = jsonParser::array();

    for(Index nb = 0; nb < perturb_config_index.size(); nb++) {
      for(Index no = 0; no < perturb_config_index[nb].size(); no++) {
        for(Index nd = 0; nd < perturb_config_index[nb][no].size(); nd++) {

          jsonParser jsonobj = jsonParser::object();

          jsonobj["orbitbranch"] = nb;
          jsonobj["orbit"] = no;

          if(print_config_name) {
            jsonobj["config"] = get_config(perturb_config_index[nb][no][nd]).name();
          }
          else {
            jsonobj["configid"] = perturb_config_index[nb][no][nd];
          }
          jsonobj["symop"] = perturb_config_symop_index[nb][no][nd];

          json["perturbations"].push_back(jsonobj);

        }
      }
    }

    json.print(file);

  }

  //***********************************************************

  void Supercell::generate_factor_group()const {
    real_super_lattice.find_invariant_subgroup(get_prim().factor_group(), m_factor_group);
    m_factor_group.set_lattice(real_super_lattice);
    return;
  }

  //***********************************************************

  void Supercell::generate_permutations()const {
    if(!m_perm_symrep_ID.empty()) {
      std::cerr << "WARNING: In Supercell::generate_permutations(), but permutations data already exists.\n"
                << "         It will be overwritten.\n";
    }
    m_perm_symrep_ID = m_prim_grid.make_permutation_representation(factor_group(), get_prim().basis_permutation_symrep_ID());
    //m_trans_permute = m_prim_grid.make_translation_permutations(basis_size()); <--moved to PrimGrid

    /*
      std::cerr << "For SCEL " << " -- " << get_name() << " Translation Permutations are:\n";
      for(int i = 0; i < m_trans_permute.size(); i++)
      std::cerr << i << ":   " << m_trans_permute[i].perm_array() << "\n";

      std::cerr << "For SCEL " << " -- " << get_name() << " factor_group Permutations are:\n";
      for(int i = 0; i < m_factor_group.size(); i++){
    std::cerr << "Operation " << i << ":\n";
    m_factor_group[i].print(std::cerr,FRAC);
    std::cerr << '\n';
    std::cerr << i << ":   " << m_factor_group[i].get_permutation_rep(m_perm_symrep_ID)->perm_array() << '\n';

    }
    std:: cerr << "End permutations for SCEL " << get_name() << '\n';
    */

    return;
  }

  //***********************************************************

  void Supercell::_generate_name() const {
    //std::cout << "begin _generate_name()" << std::endl;
    m_name = CASM::generate_name(transf_mat);
  }

  //***********************************************************

  fs::path Supercell::get_path() const {
    return get_primclex().get_path() / "training_data" / get_name();
  }

  /*
   * Run through the configuration list and count how many of them
   * have been selected then return value.
   */

  Index Supercell::amount_selected() const {
    Index amount_selected = 0;
    for(Index c = 0; c < config_list.size(); c++) {
      if(config_list[c].selected()) {
        amount_selected++;
      }
    }
    return amount_selected;
  }

  //***********************************************************

  bool Supercell::is_canonical() const {
    return get_real_super_lattice().is_canonical(
             get_prim().point_group(),
             get_primclex().crystallography_tol());
  }

  //***********************************************************

  SymOp Supercell::to_canonical() const {
    return get_real_super_lattice().to_canonical(
             get_prim().point_group(),
             get_primclex().crystallography_tol());
  }

  //***********************************************************

  SymOp Supercell::from_canonical() const {
    return get_real_super_lattice().from_canonical(
             get_prim().point_group(),
             get_primclex().crystallography_tol());
  }

  //***********************************************************

  Supercell &Supercell::canonical_form() const {
    if(!m_canonical) {
      m_canonical = &get_primclex().get_supercell(
                      get_primclex().add_supercell(get_real_super_lattice()));
    }
    return *m_canonical;
  }

  //***********************************************************
  /**  Check if a Structure fits in this Supercell
   *  - Checks that 'structure'.lattice is supercell of 'real_super_lattice'
   *  - Does *NOT* check basis sites
   */
  //***********************************************************
  bool Supercell::is_supercell_of(const Structure &structure) const {
    Eigen::Matrix3d mat;
    return is_supercell_of(structure, mat);
  };

  //***********************************************************
  /**  Check if a Structure fits in this Supercell
   *  - Checks that 'structure'.lattice is supercell of 'real_super_lattice'
   *  - Does *NOT* check basis sites
   */
  //***********************************************************
  bool Supercell::is_supercell_of(const Structure &structure, Eigen::Matrix3d &mat) const {
    Structure tstruct = structure;
    SymGroup point_group;
    tstruct.lattice().generate_point_group(point_group);
    //if(real_super_lattice.is_supercell_of(tstruct.lattice, tstruct.factor_group().point_group(), mat)) {

    if(real_super_lattice.is_supercell_of(tstruct.lattice(), point_group, mat)) {

      return true;
    }
    return false;
  };

  //***********************************************************
  /**  Generate a Configuration from a Structure
   *  - Generally expected the user will first call
   *      Supercell::is_supercell_of(const Structure &structure, Matrix3<double> multimat)
   *  - tested OK for perfect prim coordinates, not yet tested with relaxed coordinates using 'tol'
   */
  //***********************************************************
  Configuration Supercell::configuration(const BasicStructure<Site> &structure_to_config, double tol) {
    //Because the user is a fool and the supercell may not be a supercell (This still doesn't check the basis!)
    Eigen::Matrix3d transmat;
    if(!structure_to_config.lattice().is_supercell_of(get_prim().lattice(), get_prim().factor_group(), transmat)) {
      std::cerr << "ERROR in Supercell::configuration" << std::endl;
      std::cerr << "The provided structure is not a supercell of the PRIM. Tranformation matrix was:" << std::endl;
      std::cerr << transmat << std::endl;
      exit(881);
    }

    std::cerr << "WARNING in Supercell::config(): This routine has not been tested on relaxed structures using 'tol'" << std::endl;
    //std::cout << "begin config()" << std::endl;
    //std::cout << "  mat:\n" << mat << std::endl;

    const Structure &prim = (*primclex).get_prim();

    // create a 'superstruc' that fills '*this'
    BasicStructure<Site> superstruc = structure_to_config.create_superstruc(real_super_lattice);

    //std::cout << "superstruc:\n";
    //superstruc.print(std::cout);
    //std::cout << " " << std::endl;

    // Set the occuation state of a Configuration from superstruc
    //   Allow Va on sites where Va are allowed
    //   Do not allow interstitials, print an error message and exit
    Configuration config(*this);

    // Initially set occupation to -1 (for unknown) on every site
    config.set_occupation(Array<int>(num_sites(), -1));

    Index linear_index, b;
    int val;

    // For each site in superstruc, set occ index
    for(Index i = 0; i < superstruc.basis.size(); i++) {
      //std::cout << "i: " << i << "  basis: " << superstruc.basis[i] << std::endl;
      linear_index = get_linear_index(Coordinate(superstruc.basis[i]), tol);
      b = get_b(linear_index);

      // check that we're not over-writing something already set
      if(config.occ(linear_index) != -1) {
        std::cerr << "Error in Supercell::config." << std::endl;
        std::cerr << "  Adding a second atom on site: linear index: " << linear_index << " bijk: " << uccoord(linear_index) << std::endl;
        exit(1);
      }

      // check that the Molecule in superstruc is allowed on the site in 'prim'
      if(!prim.basis[b].contains(superstruc.basis[i].occ_name(), val)) {
        std::cerr << "Error in Supercell::config." << std::endl;
        std::cerr << "  The molecule: " << superstruc.basis[i].occ_name() << " is not allowed on basis site " << b << " of the Supercell prim." << std::endl;
        exit(1);
      }
      config.set_occ(linear_index, val);
    }

    // Check that vacant sites are allowed
    for(Index i = 0; i < config.size(); i++) {
      if(config.occ(i) == -1) {
        b = get_b(i);

        if(prim.basis[b].contains("Va", val)) {
          config.set_occ(i, val);
        }
        else {
          std::cerr << "Error in Supercell::config." << std::endl;
          std::cerr << "  Missing atom.  Vacancies are not allowed on the site: " << uccoord(i) << std::endl;
          exit(1);
        }
      }
    }

    return config;

  };

  //***********************************************************
  /**  Returns a Structure equivalent to the Supercell
   *  - basis sites are ordered to agree with Supercell::config_index_to_bijk
   *  - occupation set to prim default, not curr_state
   */
  //***********************************************************
  Structure Supercell::superstructure() const {
    // create a 'superstruc' that fills '*this'
    Structure superstruc = (*primclex).get_prim().create_superstruc(real_super_lattice);

    Index linear_index;
    // sort basis sites so that they agree with config_index_to_bijk
    //   This sorting may not be necessary,
    //   but it depends on how we construct the config_index_to_bijk,
    //   so I'll leave it in for now just to be safe
    for(Index i = 0; i < superstruc.basis.size(); i++) {
      linear_index = get_linear_index(superstruc.basis[i]);
      superstruc.basis.swap_elem(i, linear_index);
    }

    //superstruc.reset();

    //set_site_internals() is better than Structure::reset(), because
    //it doesn't destroy all the info that
    //Structure::create_superstruc makes efficiently
    superstruc.set_site_internals();
    return superstruc;

  }

  //***********************************************************
  /**  Returns a Structure equivalent to the Supercell
   *  - basis sites are ordered to agree with Supercell::config_index_to_bijk
   *  - occupation set to config
   *  - prim set to (*primclex).prim
   */
  //***********************************************************
  Structure Supercell::superstructure(const Configuration &config) const {
    // create a 'superstruc' that fills '*this'
    Structure superstruc = superstructure();

    // set basis site occupants
    for(Index i = 0; i < superstruc.basis.size(); i++) {
      superstruc.basis[i].set_occ_value(config.occ(i));
    }

    // setting the occupation changes symmetry properties, so must reset
    superstruc.reset();

    return superstruc;

  }

  /**
   * This is a safer version that takes an Index instead of an actual Configuration.
   * It might be better to have the version that takes a Configuration private,
   * that way you can't pass it anything that's incompatible.
   */

  Structure Supercell::superstructure(Index config_index) const {
    if(config_index >= config_list.size()) {
      std::cerr << "ERROR in Supercell::superstructure" << std::endl;
      std::cerr << "Requested superstructure of configuration with index " << config_index << " but there are only " << config_list.size() << " configurations" << std::endl;
      exit(185);
    }
    return superstructure(config_list[config_index]);
  }

  //***********************************************************
  /**  Returns an Array<int> consistent with
   *     Configuration::occupation that is all vacancies.
   *     A site which can not contain a vacancy is set to -1.
   */
  //***********************************************************
  ReturnArray<int> Supercell::vacant() const {
    Array<int> occupation = Array<int>(num_sites(), -1);
    int b, index;
    for(Index i = 0; i < num_sites(); i++) {
      b = get_b(i);
      if(get_prim().basis[b].contains("Va", index)) {
        occupation[i] = index;
      }
    }
    return occupation;
  }

  //**********************************************************
  /** Initialize a nx3 vector of real space coords,
      corresponding to the real space grid points
  **/
  //**********************************************************
  Eigen::MatrixXd Supercell::real_coordinates() const {
    Eigen::MatrixXd real_coords(volume(), 3);
    for(int i = 0; i < volume(); i++) {
      Coordinate temp_real_point = m_prim_grid.coord(i, SCEL);
      temp_real_point.within(); //should this also be voronoi within?
      real_coords.row(i) = temp_real_point.const_cart().transpose();
    }
    return real_coords;
  }

  //**********************************************************
  /**   Initialize a nx3 vector of kpoints that we are interested in
   by using kpoints that are listed in recip_grid. We will only
   generate k-points in the IBZ, so always call voronoi_within on the
   coordinate
  */
  //**********************************************************
  Eigen::MatrixXd Supercell::recip_coordinates() const {
    Eigen::MatrixXd kpoint_coords(volume(), 3);
    Lattice temp_recip_lattice = (*primclex).get_prim().lattice().get_reciprocal();
    for(int i = 0; i < volume(); i++) {
      //      std::cout<<"UCC:"<<recip_grid.uccoord(i);
      Coordinate temp_kpoint = recip_grid.coord(i, PRIM);
      temp_kpoint.set_lattice(temp_recip_lattice, CART);
      temp_kpoint.within(); //This is temporary, should be replaced by a call to voronoi_within()

      kpoint_coords.row(i) = temp_kpoint.const_cart().transpose();
    }
    return kpoint_coords;
  }

  Array< bool > Supercell::is_commensurate_kpoint(const Eigen::MatrixXd &recip_coordinates, double tol) {
    Eigen::MatrixXd recip_frac_coords = recip_coordinates * recip_prim_lattice.inv_lat_column_mat();
    //    std::cout<<"Recip Frac Coords"<<std::endl<<recip_frac_coords<<std::endl;
    Array<bool> is_commensurate(recip_coordinates.rows(), true);
    for(int i = 0; i < recip_frac_coords.rows(); i++) {
      for(int j = 0; j < recip_frac_coords.cols(); j++) {
        if(std::abs(round(recip_frac_coords(i, j)) - recip_frac_coords(i, j)) > tol) {
          is_commensurate[i] = false;
          break;
        }
      }
    }
    return is_commensurate;
  }

  void Supercell::generate_fourier_matrix() {
    generate_fourier_matrix(real_coordinates(), recip_coordinates(), true);
  }

  void Supercell::generate_fourier_matrix(const Eigen::MatrixXd &real_coordinates, const Eigen::MatrixXd &recip_coordinates) {
    generate_fourier_matrix(real_coordinates, recip_coordinates, false);
  }

  void Supercell::generate_fourier_matrix(const Eigen::MatrixXd &real_coordinates, const Eigen::MatrixXd &recip_coordinates, const bool &override) {
    //Validate the input matrices
    if((real_coordinates.cols() != 3) || recip_coordinates.cols() != 3) {
      std::cerr << "ERROR in generate_fourier_matrix, your matrices are incorrectly initialized" << std::endl;
      std::cerr << "QUITTING" << std::endl;
      exit(666);
    }
    //Check if m_k_mesh is already full
    if(m_k_mesh.rows() != 0 || m_k_mesh.cols() != 0) {
      std::cerr << "WARNING in Supercell::generate_fourier_matrix. You already have a k-mesh in this Supercell"
                << " It will be overwritten" << std::endl;
    }
    m_k_mesh = recip_coordinates;
    // Setup the size of _fourier_matrix which is only the -i*r*k'
    Eigen::MatrixXcd _fourier_matrix(real_coordinates.rows(), m_k_mesh.rows());
    std::complex<double> pre_factor(0, -1);
    _fourier_matrix = pre_factor * real_coordinates * m_k_mesh.transpose();
    //Exponentiate every element of _fourier matrix and store in m_fourier_matrix
    m_fourier_matrix = _fourier_matrix.array().exp();
    //Find all those k-point that are not commensurate with this supercell and
    // set those columns of the fourier_matrix to be zeros
    Array<bool> is_commensurate(recip_coordinates.rows(), true);
    if(!override) {
      is_commensurate = is_commensurate_kpoint(recip_coordinates);
      for(int i = 0; i < m_fourier_matrix.cols(); i++) {
        if(!is_commensurate[i]) {
          m_fourier_matrix.col(i) = Eigen::MatrixXcd::Zero(m_fourier_matrix.rows(), 1);
        }
      }
    }
    generate_phase_factor((*primclex).shift_vectors(), is_commensurate, override);
    //    std::cout<<"Fourier Matrix:"<<std::endl<<m_fourier_matrix<<std::endl;
  }

  void Supercell::generate_phase_factor(const Eigen::MatrixXd &shift_vectors, const Array<bool> &is_commensurate, const bool &override) {
    Eigen::MatrixXcd _phase_factor(basis_size(), m_k_mesh.rows());
    std::complex<double> pre_factor(0, -1);
    //    std::cout<<"Shift vectors"<<std::endl<<shift_vectors<<std::endl;
    _phase_factor = pre_factor * shift_vectors * m_k_mesh.transpose();
    m_phase_factor = _phase_factor.array().exp();
    //Zero out all the rows of m_phase_factor that have k-points that are
    //not commensurate with the Supercell
    if(!override) {
      for(int i = 0; i < m_phase_factor.cols(); i++) {
        if(!is_commensurate[i]) {
          m_phase_factor.col(i) = Eigen::MatrixXcd::Zero(m_phase_factor.rows(), 1);
        }
      }
    }
    //std::cout<<"Phase factors:"<<std::endl<<m_phase_factor<<std::endl;
  }

  void Supercell::populate_structure_factor() {
    if(m_fourier_matrix.rows() == 0 || m_fourier_matrix.cols() == 0 || m_phase_factor.rows() == 0 || m_phase_factor.cols() == 0) {
      generate_fourier_matrix();
    }
    for(Index i = 0; i < config_list.size(); i++) {
      populate_structure_factor(i);
    }
    return;
  }

  void Supercell::populate_structure_factor(const Index &config_index) {
    if(m_fourier_matrix.rows() == 0 || m_fourier_matrix.cols() == 0 || m_phase_factor.rows() == 0 || m_phase_factor.cols() == 0) {
      generate_fourier_matrix();
    }
    config_list[config_index].calc_struct_fact();
    return;
  }

  bool Supercell::operator<(const Supercell &B) const {
    if(&get_primclex() != &B.get_primclex()) {
      throw std::runtime_error(
        "Error using Supercell::operator<(const Supercell& B): "
        "Only Supercell with the same PrimClex may be compared this way.");
    }
    if(volume() != B.volume()) {
      return volume() < B.volume();
    }
    return get_real_super_lattice() < B.get_real_super_lattice();
  }

  bool Supercell::_eq(const Supercell &B) const {
    if(&get_primclex() != &B.get_primclex()) {
      throw std::runtime_error(
        "Error using Supercell::operator==(const Supercell& B): "
        "Only Supercell with the same PrimClex may be compared this way.");
    }
    return get_transf_mat() == B.get_transf_mat();
  }

  Supercell &apply(const SymOp &op, Supercell &scel) {
    return scel = copy_apply(op, scel);
  }

  Supercell copy_apply(const SymOp &op, const Supercell &scel) {
    return Supercell(&scel.get_primclex(), copy_apply(op, scel.get_real_super_lattice()));
  }


  std::string generate_name(const Eigen::Matrix3i &transf_mat) {
    std::string name_str;

    Eigen::Matrix3i H = hermite_normal_form(transf_mat).first;
    name_str = "SCEL";
    std::stringstream tname;
    //Consider using a for loop with HermiteCounter_impl::_canonical_unroll here
    tname << H(0, 0)*H(1, 1)*H(2, 2) << "_" << H(0, 0) << "_" << H(1, 1) << "_" << H(2, 2) << "_" << H(1, 2) << "_" << H(0, 2) << "_" << H(0, 1);
    name_str.append(tname.str());

    return name_str;
  }

}
