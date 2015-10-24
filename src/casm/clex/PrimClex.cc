#include "casm/clex/PrimClex.hh"

#include <boost/algorithm/string.hpp>

#include "casm/clex/ConfigIterator.hh"
#include "casm/clex/ECIContainer.hh"
#include "casm/clusterography/jsonClust.hh"
#include "casm/system/RuntimeLibrary.hh"
#include "casm/casm_io/SafeOfstream.hh"


namespace CASM {



  //*******************************************************************************************
  //                                **** Constructors ****
  //*******************************************************************************************
  /// Initial construction of a PrimClex, from a primitive Structure
  PrimClex::PrimClex(const Structure &_prim) :
    prim(_prim),
    global_orbitree(_prim.lattice()) {
    //prim.generate_factor_group();
    return;
  };


  //*******************************************************************************************
  /// Construct PrimClex from existing CASM project directory
  ///  - read PrimClex and directory structure to generate all its Supercells and Configurations, etc.
  PrimClex::PrimClex(const fs::path &_root, std::ostream &sout):
    root(_root),
    m_dir(_root),
    m_settings(_root),
    prim(read_prim(m_dir.prim())),
    global_orbitree(prim.lattice()) {

    //prim.generate_factor_group();
    bool any_print = false;

    // here add stuff to read directory structure...

    // read .casmroot current settings
    try {
      jsonParser settings(root / ".casm" / "project_settings.json");
      from_json(curr_property, settings["curr_properties"]);
      from_json(curr_clex, settings["curr_clex"]);
      from_json(curr_calctype, settings["curr_calctype"]);
      from_json(curr_ref, settings["curr_ref"]);
      from_json(curr_bset, settings["curr_bset"]);
      from_json(curr_eci, settings["curr_eci"]);
      settings.get_else(compile_options, "compile_options", RuntimeLibrary::default_compile_options());
      settings.get_else(so_options, "so_options", RuntimeLibrary::default_so_options());
      from_json(m_name, settings["name"]);
    }
    catch(std::exception &e) {
      std::cerr << "Error in PrimClex::PrimClex(const fs::path &_root, std::ostream &sout) reading .casmroot" << std::endl;
      std::cerr << e.what() << std::endl;
      exit(1);
    }

    // read param composition
    auto comp_axes = m_dir.composition_axes(m_settings.calctype(), m_settings.ref());
    if(fs::is_regular_file(comp_axes)) {
      sout << "  Read " << comp_axes << std::endl;

      CompositionAxes opt(comp_axes);

      if(opt.has_current_axes) {
        m_has_composition_axes = true;
        m_comp_converter = opt.curr;
      }
    }

    // read supercells
    if(fs::is_regular_file(root / "training_data" / "SCEL")) {

      any_print = true;
      sout << "  Read " << root / "training_data" / "SCEL" << std::endl;
      fs::ifstream scel(root / "training_data" / "SCEL");
      read_supercells(scel);

    }

    // read config_list
    if(fs::is_regular_file(get_config_list_path())) {

      any_print = true;
      sout << "  Read " << get_config_list_path() << std::endl;
      read_config_list();
    }

    if(any_print)
      sout << "\n" << std::flush;


  }


  //*******************************************************************************************
  // **** Accessors ****
  //*******************************************************************************************

  /// Return project name
  std::string PrimClex::name() const {
    return m_name;
  }

  //*******************************************************************************************
  // ** Directory path accessors **
  //*******************************************************************************************
  /// Return casm project directory path
  fs::path PrimClex::get_path() const {
    return root;
  }

  //*******************************************************************************************
  /// Return supercell directory path
  fs::path PrimClex::get_path(const Index &scel_index) const {
    return root / "training_data" / supercell_list[scel_index].get_name();
  }

  //*******************************************************************************************
  /// Return configuration directory path
  fs::path PrimClex::get_path(const Index &scel_index, const Index &config_index) const {
    return get_path(scel_index) / supercell_list[scel_index].get_config(config_index).get_id();
  }

  //*******************************************************************************************
  /// Return config_list.json file path
  fs::path PrimClex::get_config_list_path() const {
    return root / ".casm" / "config_list.json";
  }


  // ** Current settings accessors **

  //*******************************************************************************************
  /// Return current property settings
  const std::vector<std::string> &PrimClex::get_curr_property() const {
    return curr_property;
  }

  //*******************************************************************************************
  /// Return current clex settings
  std::string PrimClex::get_curr_clex() const {
    return curr_clex;
  }

  //*******************************************************************************************
  /// Return current calctype setting
  std::string PrimClex::get_curr_calctype() const {
    return curr_calctype;
  }

  //*******************************************************************************************
  /// Return current reference setting
  std::string PrimClex::get_curr_ref() const {
    return curr_ref;
  }

  //*******************************************************************************************
  /// Return cluster settings
  std::string PrimClex::get_curr_bset() const {
    return curr_bset;
  }

  //*******************************************************************************************
  /// Return current global clexulator name
  std::string PrimClex::get_curr_clexulator() const {
    return name() + "_Clexulator";
  }

  //*******************************************************************************************
  /// Return current eci settings
  std::string PrimClex::get_curr_eci() const {
    return curr_eci;
  }

  //*******************************************************************************************
  /// Return compiler options
  std::string PrimClex::get_compile_options() const {
    return compile_options;
  }

  //*******************************************************************************************
  /// Return shared library options
  std::string PrimClex::get_so_options() const {
    return so_options;
  }

  // ** Composition accessors **

  //*******************************************************************************************
  /// const Access CompositionConverter object
  bool PrimClex::has_composition_axes() const {
    return m_has_composition_axes;
  }

  //*******************************************************************************************
  /// const Access CompositionConverter object
  const CompositionConverter &PrimClex::composition_axes() const {
    return m_comp_converter;
  }



  // ** Prim and Orbitree accessors **

  //*******************************************************************************************
  /// const Access to primitive Structure
  const Structure &PrimClex::get_prim() const {
    return prim;
  }

  //*******************************************************************************************
  /// const Access to global orbitree
  const SiteOrbitree &PrimClex::get_global_orbitree() const {
    return global_orbitree;
  };

  //*******************************************************************************************

  const Array<UnitCellCoord> &PrimClex::get_prim_nlist() const {
    return prim_nlist;
  }

  // ** Supercell and Configuration accessors **

  //*******************************************************************************************
  /// const Access entire supercell_list
  const boost::container::stable_vector<Supercell> &PrimClex::get_supercell_list() const {
    return supercell_list;
  };

  //*******************************************************************************************
  /// const Access supercell by index
  const Supercell &PrimClex::get_supercell(Index i) const {
    return supercell_list[i];
  };

  //*******************************************************************************************
  /// Access supercell by index
  Supercell &PrimClex::get_supercell(Index i) {
    return supercell_list[i];
  };

  //*******************************************************************************************
  /// const Access supercell by name
  const Supercell &PrimClex::get_supercell(std::string scellname) const {
    Index index;
    if(!contains_supercell(scellname, index)) {
      std::cout << "Error in PrimClex::get_supercell(std::string scellname) const." << std::endl;
      std::cout << "  supercell '" << scellname << "' not found." << std::endl;
      exit(1);
    }
    return supercell_list[index];
  };

  //*******************************************************************************************
  /// Access supercell by name
  Supercell &PrimClex::get_supercell(std::string scellname) {
    Index index;
    if(!contains_supercell(scellname, index)) {
      std::cout << "Error in PrimClex::get_supercell(std::string scellname)." << std::endl;
      std::cout << "  supercell '" << scellname << "' not found." << std::endl;
      exit(1);
    }
    return supercell_list[index];
  };

  //*******************************************************************************************
  /// access configuration by name (of the form "scellname/[NUMBER]", e.g., ("SCEL1_1_1_1_0_0_0/0")
  const Configuration &PrimClex::configuration(const std::string &configname) const {
    std::vector<std::string> splt_vec;
    boost::split(splt_vec, configname, boost::is_any_of("/"), boost::token_compress_on);
    Index config_ind;
    if(splt_vec.size() != 2) {
      std::cerr << "CRITICAL ERROR: In PrimClex::configuration(), cannot locate configuration " << configname << "\n"
                << "                Exiting...\n";
      assert(0);
      exit(1);
    }

    try {
      config_ind = boost::lexical_cast<Index>(splt_vec[1]);
    }
    catch(boost::bad_lexical_cast &) {
      std::cerr << "CRITICAL ERROR: In PrimClex::configuration(), malformed input:" << configname << "\n"
                << "                Exiting...\n";
      assert(0);
      exit(1);
    }


    return get_supercell(splt_vec[0]).get_config(config_ind);
  }

  //*******************************************************************************************

  Configuration &PrimClex::configuration(const std::string &configname) {
    std::vector<std::string> splt_vec;
    boost::split(splt_vec, configname, boost::is_any_of("/"), boost::token_compress_on);

    Index config_ind;
    if(splt_vec.size() != 2) {
      std::cerr << "CRITICAL ERROR: In PrimClex::configuration(), cannot locate configuration " << configname << "\n"
                << "                Exiting...\n";
      assert(0);
      exit(1);
    }

    try {
      config_ind = boost::lexical_cast<Index>(splt_vec[1]);
    }
    catch(boost::bad_lexical_cast &) {
      std::cerr << "CRITICAL ERROR: In PrimClex::configuration(), malformed input:" << configname << "\n"
                << "                Exiting...\n";
      assert(0);
      exit(1);
    }

    return get_supercell(splt_vec[0]).get_config(config_ind);
  }

  //*******************************************************************************************

  /// Configuration iterator: begin
  PrimClex::config_iterator PrimClex::config_begin() {
    if(supercell_list.size() == 0 || supercell_list[0].get_config_list().size() > 0)
      return config_iterator(this, 0, 0);
    return ++config_iterator(this, 0, 0);
  }

  //*******************************************************************************************
  /// Configuration iterator: end
  PrimClex::config_iterator PrimClex::config_end() {
    return config_iterator(this, supercell_list.size(), 0);
  }

  //*******************************************************************************************
  /// const Configuration iterator: begin
  PrimClex::config_const_iterator PrimClex::config_cbegin() const {
    if(supercell_list.size() == 0 || supercell_list[0].get_config_list().size() > 0)
      return config_const_iterator(this, 0, 0);
    return ++config_const_iterator(this, 0, 0);
  }

  //*******************************************************************************************
  /// const Configuration iterator: end
  PrimClex::config_const_iterator PrimClex::config_cend() const {
    return config_const_iterator(this, supercell_list.size(), 0);
  }


  //*******************************************************************************************
  /// Configuration iterator: begin
  PrimClex::config_iterator PrimClex::selected_config_begin() {
    //std::cout << "BEGINNING SELECTED CONFIG ITERATOR\n"
    //          << "supercell_list.size() is " << supercell_list.size() << "\n";

    if(supercell_list.size() == 0 || (supercell_list[0].get_config_list().size() > 0 && supercell_list[0].get_config(0).selected()))
      return config_iterator(this, 0, 0, true);
    return ++config_iterator(this, 0, 0, true);
  }

  //*******************************************************************************************
  /// Configuration iterator: end
  PrimClex::config_iterator PrimClex::selected_config_end() {
    return config_iterator(this, supercell_list.size(), 0, true);
  }

  //*******************************************************************************************
  /// const Configuration iterator: begin
  PrimClex::config_const_iterator PrimClex::selected_config_cbegin() const {
    if(supercell_list.size() == 0 || (supercell_list[0].get_config_list().size() > 0 && supercell_list[0].get_config(0).selected()))
      return config_const_iterator(this, 0, 0, true);
    return ++config_const_iterator(this, 0, 0, true);
  }

  //*******************************************************************************************
  /// const Configuration iterator: end
  PrimClex::config_const_iterator PrimClex::selected_config_cend() const {
    return config_const_iterator(this, supercell_list.size(), 0, true);
  }


  // ** Neighbor list accessors **
  //*******************************************************************************************
  /// Return number of sites in primitive cell neighborhood
  Index PrimClex::get_nlist_size() const {
    return prim_nlist.size();
  }

  //*******************************************************************************************
  /// Returns UnitCellCoord, 'delta', indicating where neighbor site 'nlist_index' is located in the neighborhood
  ///    neighbor_unit_cell_coord = UnitCellCoord(b,i,j,k) + this->get_nlist_uccoord(nlist_index)
  const UnitCellCoord &PrimClex::get_nlist_uccoord(Index nlist_index) const {
    return prim_nlist[nlist_index];
  }


  //*******************************************************************************************
  // **** IO ****
  //*******************************************************************************************
  /**
   * Re-write config_list.json, updating all the data
   */

  void PrimClex::write_config_list() {

    if(supercell_list.size() == 0) {
      fs::remove(get_config_list_path());
      return;
    }

    jsonParser json;

    if(fs::exists(get_config_list_path())) {
      json.read(get_config_list_path());
    }
    else {
      json.put_obj();
    }

    for(Index s = 0; s < supercell_list.size(); s++) {
      supercell_list[s].write_config_list(json);
    }

    SafeOfstream file;
    file.open(get_config_list_path());
    json.print(file.ofstream());
    file.close();

    return;
  }


  // **** Operators ****


  // **** Unorganized mess of functions ... ****

  //***********************************************************
  /**
   * Calculates the distances between all sites of a cluster from
   * the first site in the cluster of all clusters in the orbitree,
   * then stores this "delta" in the primitive neighbor list if it
   * hasn't been added already.
   * In the process, it updates the value
   * of nlist_ind() of every site in the given orbitree to
   * hold the index of where said site was placed in the neighbor
   * list, relative to the first site of the cluster.
   *
   * Example:
   * A triple cluster made from
   * Cluster site       |       Site in primitive       |       nlist_ind
   * 0                          2                               2
   * 1                          6                               7
   * 2                          1                               8
   *
   * Cluster site is just the index of the site inside the cluster (which is of type Array<Site>)
   * Site in primitive tells you what site in the basis of your primitive your cluster site corresponds to (inconsequential here)
   * nlist_ind is what you use to read prim_nlist and gets set in this routine.
   *
   * The entire prim_nlist is based off of the first site in the cluster. To get the
   * distance between cluster site 0 and 1, you look at prim_nlist[7], to get
   * the distance between cluster site 0 and 2, you look at prim_nlist[8].
   *
   */
  //***********************************************************

  void PrimClex::append_to_nlist(SiteOrbitree &new_tree) {
    // We include all the primitive-cell neighbors by default at the front of the prim_nlist
    // This might be overkill in the long term (e.g., for large primitive cells where only a few sites have DoFs)
    // In that case, we may change tactics and construct a separate list of nlist_inds that correspond to prim-cell neighbors
    // We need this information in some form for calculating differential correlations
    if(!prim_nlist.size()) {
      for(Index b = 0; b < prim.basis.size(); b++)
        prim_nlist.push_back(UnitCellCoord(b, 0, 0, 0));
    }

    Array<Index> clust_nlist_inds;

    //For each site we encounter we access the appropriate slot in the neighbor list and append all other sites
    //to it in the form of UnitCellCoords

    //branches
    for(Index i = 0; i < new_tree.size(); i++) {
      //orbits
      for(Index j = 0; j < new_tree[i].size(); j++) {
        //clusters
        for(Index k = 0; k < new_tree[i][j].size(); k++) {
          //tuccl corresponds to the first site in the cluster
          UnitCellCoord tuccl = prim.get_unit_cell_coord(new_tree[i][j][k][0]);

          clust_nlist_inds.resize(new_tree[i][j][k].size());
          //neighbor sites
          for(Index b = 0; b < new_tree[i][j][k].size(); b++) {
            //tuccb corresponds to a site that neighbors tuccl
            UnitCellCoord tuccb = prim.get_unit_cell_coord(new_tree[i][j][k][b]);
            UnitCellCoord delta = tuccb - tuccl;

            clust_nlist_inds[b] = prim_nlist.find(delta);

            //If delta is not already in prim_nlist, add it (the value of clust_nlist_inds[b] makes sense after the push_back)
            if(clust_nlist_inds[b] == prim_nlist.size()) {
              prim_nlist.push_back(delta);
            }
          }
          //Update nlist_inds of the cluster sites so they know where they live in the neighbor list
          new_tree[i][j][k].set_nlist_inds(clust_nlist_inds);
        }
      }
    }

    return;
  };

  //***********************************************************
  /**
   * Calculates the distances between all sites from each other
   * in every cluster of the orbitree, then stores this "delta"
   * in the primitive neighbor list if it hasn't been added already.
   * In the process, it updates the value
   * of nlist_ind of every site in the given orbitree to
   * hold the index of where said site was placed in the neighbor
   * list, relative to the first site of the cluster.
   *
   * This routine is essentially the same as PrimClex::append_to_nlist,
   * except it compares all sites in a cluster to all other sites
   * instead of just the first one (note additional loop over 'l').
   * The result is a neighbor list that contains transitional "delta",
   * without using flowertrees. It will be equal in size or larger than
   * neighbor list generated by PrimClex::append_to_nlist
   */
  //***********************************************************

  void PrimClex::append_to_nlist_perm(SiteOrbitree &new_tree) {

    // We include all the primitive-cell neighbors by default at the front of the prim_nlist
    // This might be overkill in the long term (e.g., for large primitive cells where only a few sites have DoFs)
    // In that case, we may change tactics and construct a separate list of nlist_inds that correspond to prim-cell neighbors
    // We need this information in some form for calculating differential correlations
    if(!prim_nlist.size()) {
      for(Index b = 0; b < prim.basis.size(); b++)
        prim_nlist.push_back(UnitCellCoord(b, 0, 0, 0));
    }



    //For each site we encounter we access the appropriate slot in the neighbor list and append all other sites
    //to it in the form of UnitCellCoords

    Array<Index> clust_nlist_inds;

    //branches
    for(Index i = 0; i < new_tree.size(); i++) {
      //orbits
      for(Index j = 0; j < new_tree[i].size(); j++) {
        //clusters
        for(Index k = 0; k < new_tree[i][j].size(); k++) {

          clust_nlist_inds.resize(new_tree[i][j][k].size());

          //sites
          for(Index l = 0; l < new_tree[i][j][k].size(); l++) {
            //tuccl corresponds to a particular site we're looking at
            UnitCellCoord tuccl = prim.get_unit_cell_coord(new_tree[i][j][k][l]);
            //neighbor sites
            for(Index b = 0; b < new_tree[i][j][k].size(); b++) {
              //tuccb corresponds to a site that neighbors tuccl
              UnitCellCoord tuccb = prim.get_unit_cell_coord(new_tree[i][j][k][b]);
              UnitCellCoord delta = tuccb - tuccl;

              clust_nlist_inds[b] = prim_nlist.find(delta);

              //If delta is not already in prim_nlist, add it (the value of clust_nlist_inds[b] makes sense after the push_back)
              if(clust_nlist_inds[b] == prim_nlist.size()) {
                prim_nlist.push_back(delta);
              }
            }

            // we set the first set of nlist_inds as the current indices
            if(l == 0) {
              new_tree[i][j][k].set_nlist_inds(clust_nlist_inds);
            }

            // save any other unique ones that we find
            Index t(0);
            for(t = 0; t < new_tree[i][j][k].trans_nlists().size(); t++) {
              if(clust_nlist_inds.all_in(new_tree[i][j][k].trans_nlist(t)))
                break;
            }
            if(t == new_tree[i][j][k].trans_nlists().size())
              new_tree[i][j][k].add_trans_nlist(clust_nlist_inds);
          }
        }
      }
    }

    return;
  };

  //*******************************************************************************************
  /**
   * Begin by constructing the flowertrees from the global
   * orbitree, then use every flowertree to get a full neighbor
   * list and updated nlist_ind of the sites in the flowertree. Then
   * update nlist_ind of each site the global orbitree, which won't add
   * any new sites to the neighbor list.
   *
   * You might not want to actually use this routine if you're
   * just interested in a global cluster expansion. You could
   * generate fewer neighbors (no flowertrees) and have all the
   * information you need. This routine would generate the
   * FULL neighbor list for every orbitree that lives in
   * the PrimClex.
   */
  //*******************************************************************************************

  void PrimClex::generate_full_nlist() {
    if(prim_nlist.size()) {
      std::cerr << "WARNING: In PrimClex::generate_full_nlist(), prim_nlist was previously initialized and is about to be overwritten!!\n";
      prim_nlist.clear();
    }

    for(Index b = 0; b < get_prim().basis.size(); b++) {
      if(get_prim().basis[b].site_occupant().size() > 1) {
        m_dof_manager.add_dof(get_prim().basis[b].site_occupant().type_name());
        break;
      }
    }
    for(Index b = 0; b < get_prim().basis.size(); b++) {
      for(Index i = 0; i < get_prim().basis[i].displacement().size(); i++)
        m_dof_manager.add_dof(get_prim().basis[b].displacement()[i].type_name());
    }

    //Update nlist_ind for global orbitree
    append_to_nlist_perm(global_orbitree);

    //Populate flowertrees for calculating correlations later and update indices of the sites
    /*prim.generate_flowertrees(global_orbitree, flowertrees);

    std::cout << "nlist size before flowertrees: " << prim_nlist.size() << " \n";
    for(Index q = 0; q < flowertrees.size(); q++) {
      append_to_nlist(flowertrees[q]);
      //flowertrees[q].print_full_clust(std::cout);
      //std::cout << std::endl;
    }
    std::cout << "nlist size after flowertrees: " << prim_nlist.size() << " \n";
    */

    m_dof_manager.resize_neighborhood(prim_nlist.size());

    // We can add more as needed
    m_dof_manager.register_dofs(global_orbitree);
    //std::ofstream clex_stream("clexulator.cc");
    //print_clexulator(clex_stream, "MyClexulator", global_orbitree);
    //clex_stream.close();

    return;
  };

  //*******************************************************************************************
  /**
   * Makes flowertrees out of the global orbitree that lives
   * in PrimClex
   */
  //*******************************************************************************************

  void PrimClex::populate_flowertrees() {
    //Populate flowertrees in one swoop
    prim.generate_flowertrees(global_orbitree, flowertrees);
    return;
  };
  //\John G 070713


  //John G 011013
  //*******************************************************************************************
  /**
   * Generates the global orbitree that lives
   * in PrimClex using the given CSPECS
   */
  //*******************************************************************************************
  void PrimClex::generate_global_orbitree(const fs::path &cspecs) {

    std::cout  << "  Read " << cspecs << std::endl << std::endl;

    fs::ifstream in_clust(cspecs);
    global_orbitree.read_CSPECS(in_clust);
    global_orbitree.min_num_components = 2;
    global_orbitree.min_length = 0.0001;
    in_clust.close();
    global_orbitree.generate_orbitree(prim);
    global_orbitree.collect_basis_info(prim);
    //std::cout << "----------------------------------------------\n\t\tGLOBAL ORBITREE\n----------------------------------------------\n ";
    //global_orbitree.print_full_clust(std::cout);
    return;
  };
  //\John G 011013

  //*******************************************************************************************
  void PrimClex::read_global_orbitree(const std::string &fclust_path) {
    read_global_orbitree(fs::path(fclust_path));
  }

  //*******************************************************************************************
  void PrimClex::read_global_orbitree(const fs::path &fclust_path) {

    global_orbitree.min_num_components = 2;     //What if we want other things?
    global_orbitree.min_length = 0.0001;
    from_json(jsonHelper(global_orbitree, prim), jsonParser(fclust_path.string()));
    global_orbitree.collect_basis_info(prim);

  }

  //*******************************************************************************************
  void PrimClex::read_custom_clusters(const std::string &custom_clusters, bool subclusters) {
    read_custom_clusters(fs::path(custom_clusters), subclusters);
  }
  //*******************************************************************************************
  void PrimClex::read_custom_clusters(const fs::path &custom_clusters, bool subclusters) {
    //Empty Orbitrees seem to have nothing set, I dont think this is a
    //safe way to create empty objects. It would be good if the
    //default constructor set them to 'invalid' values
    if(global_orbitree.size() == 0) {
      global_orbitree.min_num_components = 2;
      global_orbitree.min_length = 0.0001;
    }
    global_orbitree.read_custom_clusters_from_json(jsonParser(custom_clusters), prim, prim.factor_group(), subclusters);
  }

  //*******************************************************************************************
  /**  GENERATE_SUPERCELLS
   *   Generates all unique supercells between volStart and
   *   volEnd of PRIM volumes. Call this routine before you
   *   enumerate configurations.
   *
   *  ARN 100213
   */
  //*******************************************************************************************
  void PrimClex::generate_supercells(int volStart, int volEnd, bool verbose) {
    Array < Lattice > supercell_lattices;
    prim.lattice().generate_supercells(supercell_lattices, prim.factor_group(), volEnd, volStart);    //point_group?
    for(Index i = 0; i < supercell_lattices.size(); i++) {
      Index list_size = supercell_list.size();
      Index index = add_canonical_supercell(supercell_lattices[i]);
      if(supercell_list.size() != list_size) {
        std::cout << "  Generated: " << supercell_list[index].get_name() << "\n";
      }
      else {
        std::cout << "  Generated: " << supercell_list[index].get_name() << " (already existed)\n";
      }
    }

    //std::cout << supercell_lattices.size() << " supercells were generated\n";

  }
  //*******************************************************************************************
  void PrimClex::generate_supercell_nlists() {
    for(Index i = 0; i < supercell_list.size(); i++) {
      supercell_list[i].generate_neighbor_list();
    }
  }

  //*******************************************************************************************
  /**
   *  add a supercell if it doesn't already exist
   *    return the index of the supercell in supercell_list
   *
   *    TODO:
   *    Check to see if superlat is linear combination of previously
   *    enumerated Supercell. You'll have to transform the coordinates
   *    if you're trying to import something with a supercell that
   *    already exists. Should return operations involved.
   *    If supercell doesn't exist add it as it is without transforming
   */
  //*******************************************************************************************
  Index PrimClex::add_canonical_supercell(const Lattice &superlat) {

    if(!superlat.is_supercell_of(prim.lattice())) {
      std::cerr << "ERROR in PrimClex::add_canoical_supercell()." << std::endl
                << "  'superlat' lattice primitive is not prim lattice primitive" << std::endl
                << "  So that everything plays nicely, make sure the supercell lattice" << std::endl
                << "  is based on PrimClex::prim." << std::endl;
      exit(1);
    }

    // temporary version, just check transf_mat equivalence
    // Does this check for equivalent supercells with different transformation matrices? Seems like it should
    // Insert second loop that goes over a symmetry operation list and applies it to the transformation matrix
    Supercell scel(this, superlat);
    scel.set_id(supercell_list.size());
    for(Index i = 0; i < supercell_list.size(); i++) {
      if(supercell_list[i].get_transf_mat() == scel.get_transf_mat())
        return i;
    }

    // if not already existing, add it
    supercell_list.push_back(scel);
    return supercell_list.size() - 1;
  }
  //*******************************************************************************************
  /**
   *  add a supercell if it doesn't already exist
   *  return the index of the supercell in supercell_list
   *
   *  Unlike the initial version above this, this routine
   *  will use an Array of SymOps (probably the point group of the lattice)
   *  to determine whether the provided superlat is equivalent to a Supercell in the
   *  PrimClex list. It will then populate the transformation matrix that
   *  maps the passed lattice onto the one on the list (i.e. the one on on
   *  the list is treated as 'primitive'). If it doesn't map, then the
   *  matrix gets flattened to zeros.
   *
   *  This routine should replace the original one once it converges
   *  to something people can agree on.
   */
  //*******************************************************************************************
  Index PrimClex::add_supercell(const Lattice &superlat) {
    return add_canonical_supercell(niggli(superlat, prim.point_group(), tol()));

  }

  //*******************************************************************************************
  /**  ENUMERATE_ALL_CONFIGURATIONS
   *   Loops through all the supercells in supercell_list and calls
   *   the enumerate_configuration routines on all of them
   */
  //*******************************************************************************************
  void PrimClex::enumerate_all_configurations() {
    for(Index i = 0; i < supercell_list.size(); i++) {
      supercell_list[i].enumerate_all_occupation_configurations();
    }
  }

  //*******************************************************************************************
  void PrimClex::print_supercells() const {


    // also write supercells/supercell_path directories with LAT files
    try {
      fs::create_directory("training_data");
    }
    catch(const fs::filesystem_error &ex) {
      std::cerr << "Error in PrimClex::print_supercells()." << std::endl;
      std::cerr << ex.what() << std::endl;
    }

    BP::BP_Write scelfile((fs::path("training_data") / "SCEL").string());
    scelfile.newfile();

    print_supercells(scelfile.get_ostream());

    for(Index i = 0; i < supercell_list.size(); i++) {
      try {
        fs::create_directory(supercell_list[i].get_path());

        fs::path latpath = supercell_list[i].get_path() / "LAT";
        if(!fs::exists(latpath)) {
          BP::BP_Write latfile(latpath.string());
          latfile.newfile();
          supercell_list[i].get_real_super_lattice().print(latfile.get_ostream());
        }
      }
      catch(const fs::filesystem_error &ex) {
        std::cerr << "Error in PrimClex::print_supercells()." << std::endl;
        std::cerr << ex.what() << std::endl;
      }
    }
  }

  //*******************************************************************************************
  void PrimClex::print_supercells(std::ostream &stream) const {
    for(Index i = 0; i < supercell_list.size(); i++) {
      stream << "Supercell Name: '" << supercell_list[i].get_name() << "' Number: " << i << " Volume: " << supercell_list[i].get_transf_mat().determinant() << "\n";
      stream << "Supercell Transformation Matrix: \n";
      stream << supercell_list[i].get_transf_mat();
      stream << "\n";
    }
  }

  //*******************************************************************************************
  void PrimClex::read_supercells(std::istream &stream) {
    // expect a file with format:
    //
    // Supercell Number: 0 Volume: 1
    // Supercell Transformation Matrix:
    //  1 0 0
    //  0 1 0
    //  0 0 1
    //
    // Supercell Number: 1 Volume: 2
    // Supercell Transformation Matrix:
    //  1 0 -1
    //  0 1 0
    //  0 0 2

    Matrix3<double> mat;

    std::string s;
    while(!stream.eof()) {
      std::getline(stream, s);
      if(s[0] == 'S') {
        std::getline(stream, s);
        stream >> mat;

        //Lattice lat(prim.lattice.coord_trans(CASM::FRAC)*mat);
        Lattice lat(prim.lattice().coord_trans(FRAC)*mat);
        add_canonical_supercell(lat);
      }
    }
  }

  //*******************************************************************************************
  /**
   * Use the global_orbitree of *this to populate the global correlations in the
   * specified configuration config_index of the given supercell.
   * The method will match the size of the correlations in the configuration
   * to said orbitree and pass it down to methods of orbitree, filling up
   * values.
   */
  //*******************************************************************************************
  /*
    void PrimClex::populate_global_correlations(const Index &scell_index, const Index &config_index) {
      supercell_list[scell_index].populate_correlations(global_clexulator, config_index);
      return;
    }

    void PrimClex::populate_global_correlations(const Index &scell_index) {
      supercell_list[scell_index].populate_correlations(global_clexulator);
      return;
    }

    void PrimClex::generate_global_correlations() {
      for(Index i = 0; i < supercell_list.size(); i++) {
        supercell_list[i].populate_correlations(global_clexulator);
      }
      return;
    }

    void PrimClex::print_correlations(std::string corrFileName) {
      std::ofstream corrFile;
      corrFile.open(corrFileName.c_str());
      for(Index i = 0; i < supercell_list.size(); i++) {
        supercell_list[i].print_clex_correlations(corrFile);
      }
    }
  */
  //*******************************************************************************************
  /*
   * Just calls print_global_correlations_simple on every supercell
   * that lives in *this PrimClex, which in turns calls the required print
   * routine on every one of its configurations. Meant for creating old CASM style
   * corr.in files.
   */
  //*******************************************************************************************
  /*
    void PrimClex::print_global_correlations_simple(std::ostream &corrstream) const {
      for(Index s = 0; s < supercell_list.size(); s++) {
        supercell_list[s].print_global_correlations_simple(corrstream);
      }
      return;
    }

    */

  //*******************************************************************************************
  //Construct the basis set for each site
  void PrimClex::populate_basis_tables(const char &basis_type) {
    prim.fill_occupant_bases(basis_type);
  }

  //*******************************************************************************************
  void PrimClex::populate_cluster_basis_function_tables() {
    global_orbitree.generate_clust_bases();
    //global_orbitree.fill_discrete_bases_tensors();
  }

  //*******************************************************************************************
  // void PrimClex::read_relaxed_structure(Index superNum, Index configNum) {
  //   supercell_list[superNum].read_relaxed_structure(configNum, prim.lattice);
  // }

  //*******************************************************************************************
  // void PrimClex::collect_clex_relaxations() {
  //   for(Index i = 0; i < supercell_list.size(); i++) {
  //     supercell_list[i].read_clex_relaxations(prim.lattice);
  //   }
  // }

  //*******************************************************************************************
  /**
   *   Read the config_list file at 'file_name' and adds config
   *   to supercell, assuming it is already canonical.
   *   If 'print_dirs', call Supercell::print_clex_configuration()
   *   for each config to be made.
   */
  //*******************************************************************************************
  void PrimClex::read_config_list() {

    jsonParser json(get_config_list_path());

    if(!json.contains("supercells")) {
      return;
    }

    for(Index i = 0; i < supercell_list.size(); i++) {
      supercell_list[i].read_config_list(json);
    }
  }

  //*******************************************************************************************
  /*
   * Run through all the supercells and add up how many configurations are selected
   * in each one of them.
   */
  //*******************************************************************************************

  int PrimClex::amount_selected() const {
    int amount_selected = 0;
    for(Index s = 0; s < supercell_list.size(); s++) {
      amount_selected += supercell_list[s].amount_selected();
    }
    return amount_selected;
  }

  //*******************************************************************************************
  bool PrimClex::contains_supercell(std::string scellname, Index &index) const {
    for(Index i = 0; i < supercell_list.size(); i++) {
      if(supercell_list[i].get_name() == scellname) {
        index = i;
        return true;
      }
    }
    index = supercell_list.size();
    return false;

  };

  //*******************************************************************************************
  Matrix3<int> PrimClex::calc_transf_mat(const Lattice &superlat) const {
    Matrix3<int> tmp_transf_mat;
    Matrix3<double> ttrans = prim.lattice().coord_trans(FRAC).inverse() * superlat.coord_trans(FRAC);
    if(ttrans.is_integer()) {
      for(int i = 0; i < 3; i++) {
        for(int j = 0; j < 3; j++) {
          tmp_transf_mat(i, j) = round(ttrans(i, j));
        }
      }
    }
    else {
      std::cerr << "Error in PrimClex::calc_transf_mat(const Lattice &superlat)" << std::endl
                << "  Bad supercell, the transformation matrix is not integer. Exiting!" << std::endl;
      exit(1);
    }
    return tmp_transf_mat;
  }

  //*******************************************************************************************
  /// \brief Sets the composition axes
  ///
  /// Also:
  /// - updates all configuration references,
  /// - writes the updated configuration info
  /// - does not update composition_axes file
  void PrimClex::set_composition_axes(const CompositionConverter &_converter) {

    m_comp_converter = _converter;
    m_has_composition_axes = true;

    // We need some way to make sure param compositions get updated...
    generate_references();

  }

  //*******************************************************************************************
  /// Delete 'properties.ref_state.X.json' files,
  /// Then call 'clear_reference_properties'
  void PrimClex::clear_reference_states() {
    for(int i = 0; i < composition_axes().independent_compositions() + 1; i++) {
      fs::remove(dir().ref_state(settings().calctype(), settings().ref(), i));
    }
    generate_references();
  }

  //*******************************************************************************************
  /// Sets the root reference state to be the calculated properties of the chosen config
  /// Does not call 'generate_references' so written files will be out-of-date until you do so!!!
  void PrimClex::set_reference_state(int refid, const Configuration &config) {

    //std::cout << "begin set_reference_state()" << std::endl;

    std::string reason_invalid;
    if(!valid_reference_state(refid, config, reason_invalid)) {
      std::cerr << "Error in PrimClex::set_reference_state." << std::endl
                << "  Could not use  SCEL: " << config.get_supercell().get_name() << " ConfigID: " << config.get_id() << "  for reference state: " << refid << std::endl
                << "  " << reason_invalid << std::endl;
      exit(1);
    }

    jsonParser json;

    json["supercell_name"] = config.get_supercell().get_name();
    json["configid"] = config.get_id();

    json["param_composition"] = config.get_param_composition();
    json["ref_state"] = config.calc_properties();

    //std::cout << "Set refid: " << refid << std::endl;
    //std::cout << "Reference State:\n" << json << "\n" << std::endl;

    json.write(dir().ref_state(settings().calctype(), settings().ref(), refid));

    //std::cout << "finish set_reference_state()" << std::endl;
  }

  //*******************************************************************************************
  /// Check that it is valid to use 'config' as reference state 'refid', returns bool and if false, sets 'reason_invalid'
  ///   Currently checks:
  ///     1) that the necessary properties have been calculated,
  ///     2) that the same Configuration is not being used twice
  ///   Needs to check that reference states span composition space
  bool PrimClex::valid_reference_state(int refid, const Configuration &config, std::string &reason_invalid) const {

    // Check that the Configuration has all the curr_property calculated
    for(Index i = 0; i < get_curr_property().size(); i++) {
      if(!config.calc_properties().contains(get_curr_property()[i])) {
        reason_invalid = "You are attempting to use a Configuration for which the property '" + get_curr_property()[i] + "' has not been calculated.";
        return false;
      }
    }

    // Check that a Configuration is not being used for multiple reference states
    for(int i = 0; i < composition_axes().independent_compositions() + 1; i++) {
      if(i == refid)
        continue;

      if(!fs::exists(dir().ref_state(settings().calctype(), settings().ref(), i)))
        continue;

      jsonParser json(dir().ref_state(settings().calctype(), settings().ref(), i));

      if(json["configid"] == "custom" || json["supercell_name"] == "custom")
        continue;

      if(json["configid"] == config.get_id() && json["supercell_name"] == config.get_supercell().get_name()) {
        reason_invalid =  "You are attempting to use the same Configuration for >1 reference state and I can't allow that.";
        return false;
      }
    }

    // Here I should check that references span the necessary composition space, but I'm not yet

    //   First read all reference states that have already been set, then check that this new one is not a linear combination of those

    return true;
  }

  //*******************************************************************************************
  /// find calculated configurations closest to
  /// [0, 0, 0, ...], [1, 0, 0, ...], [0, 1, 0, ...], [0, 0, 1, ...], ...
  /// and set them as the root reference states, also calls generate_references
  /// Clears refrence states and properties whether or not it succeeds
  void PrimClex::set_reference_state_auto() {

    //std::cout << "begin set_reference_state_auto()" << std::endl;

    /// find calculated configurations closest to
    /// [0, 0, 0, ...], [1, 0, 0, ...], [0, 1, 0, ...], [0, 0, 1, ...], ...

    // Clear current references
    clear_reference_states();

    int Naxes = composition_axes().independent_compositions();

    //std::cout << "Naxes: " << Naxes << std::endl;

    Eigen::VectorXd target = Eigen::VectorXd::Zero(Naxes);

    //std::cout << "Ref: " << 0 << "  target: " << target.transpose() << std::endl;
    set_reference_state(0, closest_calculated_config(target));

    for(int i = 0; i < Naxes; i++) {

      target(i) = 1.0;
      //std::cout << "Ref: " << i + 1 << "  target: " << target.transpose() << std::endl;
      set_reference_state(i + 1, closest_calculated_config(target));
      target(i) = 0.0;
    }

    //std::cout << "generate references" << std::endl;
    generate_references();

    //std::cout << "finish set_reference_state_auto()" << std::endl;

  }


  //*******************************************************************************************
  /// Clear 'reference' and 'delta' properties from all Configurations
  /// Re-write all Configurations, updating:
  ///   param_composition.json
  ///   properties.calc.json
  ///   properties.ref.json
  ///   properties.delta.json
  void PrimClex::generate_references() {

    //std::cout << "begin PrimClex::generate_references()" << std::endl;

    for(config_iterator it = config_begin(); it != config_end(); ++it) {
      it->generate_reference();
    }

    write_config_list();

    //std::cout << "finish PrimClex::generate_references()" << std::endl;

  }


  //*******************************************************************************************
  /// private:

  //*******************************************************************************************
  /// Return the configuration closest in param_composition to the target_param_comp
  ///   Tie break returns configuration in smallest supercell (first found at that size)
  const Configuration &PrimClex::closest_calculated_config(const Eigen::VectorXd &target_param_comp) const {

    //std::cout << "begin closest_calculated_config()" << std::endl;

    /// return reference to Configuration with param_comp closest to target_param_comp
    ///   tie break goes to first Configuration with fewest atoms
    ///
    ///   must be Configurations for which the curr_properties have been calculated

    Eigen::VectorXd param_comp;
    Eigen::VectorXd closest_comp;

    double curr_dist;
    double close_dist = -1;
    Index close_super = -1;
    Index close_config;
    Index close_size;

    for(Index i = 0; i < supercell_list.size(); i++) {
      for(Index j = 0; j < supercell_list[i].get_config_list().size(); j++) {

        const Configuration &config = supercell_list[i].get_config(j);

        // check if the config has been calculated
        //std::cout << "\n\nScell: " << supercell_list[i].get_name() << "  Config: " << j << "\n" << config.calc_properties() << std::endl;

        if(config.calc_properties().size() == 0)
          continue;

        curr_dist = (target_param_comp - config.get_param_composition()).norm();

        if(!valid_index(close_super) ||
           (almost_equal(curr_dist, close_dist, TOL) && config.size() < close_size) ||
           (curr_dist < close_dist)) {
          close_super = i;
          close_config = j;
          close_dist = curr_dist;
          close_size = config.size();
        }

      } // for j
    } // for i

    if(!valid_index(close_super)) {
      std::cerr << "Error in PrimClex::closest_calculated_config" << std::endl
                << "  Could not find a calculated Configuration." << std::endl;
      exit(1);
    }

    //std::cout << "Closest Calculated:  SCEL: " << supercell_list[close_super].get_name() << "  Config: " <<
    //          close_config << "  ParamComp: " << supercell_list[close_super].get_config(close_config).get_param_composition().transpose() << std::endl;

    //std::cout << "finish closest_calculated_config()" << std::endl;

    return supercell_list[close_super].get_config(close_config);
  }


  //*******************************************************************************************
  Eigen::MatrixXd PrimClex::shift_vectors() const {
    Eigen::MatrixXd tau(prim.basis.size(), 3);
    for(int i = 0; i < prim.basis.size(); i++) {
      for(int j = 0; j < 3; j++)
        tau(i, j) = prim.basis[i].get(j, CART);
    }
    return tau;
  }

  //*******************************************************************************************
  void PrimClex::populate_structure_factor() {
    for(Index i = 0; i < supercell_list.size(); i++)
      populate_structure_factor(i);
    return;
  }

  //*******************************************************************************************
  void PrimClex::populate_structure_factor(const Index &scell_index) {
    supercell_list[scell_index].populate_structure_factor();
    return;
  }

  //*******************************************************************************************
  void PrimClex::populate_structure_factor(const Index &scell_index, const Index &config_index) {
    supercell_list[scell_index].populate_structure_factor(config_index);
    return;
  }

  //*******************************************************************************************
  Clexulator PrimClex::global_clexulator() const {
    if(!m_global_clexulator.initialized()) {
      if(!fs::exists(dir().clexulator_src(settings().name(), settings().bset()))) {
        throw std::runtime_error(
          std::string("Error loading clexulator ") + settings().bset() + ". No basis functions exist.");
      }
      m_global_clexulator = Clexulator(settings().global_clexulator(),
                                       dir().clexulator_dir(settings().bset()),
                                       settings().compile_options(),
                                       settings().so_options());
    }
    return m_global_clexulator;
  }

  //*******************************************************************************************
  ECIContainer PrimClex::global_eci(std::string clex_name)const {
    return ECIContainer(dir().eci_out(clex_name,
                                      settings().calctype(),
                                      settings().ref(),
                                      settings().bset(),
                                      settings().eci()));
  }

  //*******************************************************************************************
  /// \brief Make orbitree. For now specifically global.
  ///
  /// \param prim Primitive Structure. non-const due to Structure::set_site_internals.
  /// \param json bspecs.json file
  SiteOrbitree make_orbitree(Structure &prim, const jsonParser &json) {

    try {

      SiteOrbitree tree(prim.lattice());

      // --- first generate global orbitree -------------


      //global_orbitree.read_CSPECS(in_clust);
      tree.min_num_components = 2;
      tree.min_length = CASM::TOL;

      tree.max_length.clear();
      auto update_max_length = [&](int branch, double max_length) {
        while(branch > tree.max_length.size() - 1) {
          tree.max_length.push_back(0.0);
        }
        tree.max_length[branch] = max_length;
      };

      for(auto it = json["orbit_branch_specs"].cbegin(); it != json["orbit_branch_specs"].cend(); ++it) {
        update_max_length(std::stoi(it.name()), it->find("max_length")->get<double>());
      }
      tree.max_num_sites = tree.max_length.size() - 1;

      tree.generate_orbitree(prim);

      // --- then add custom orbits --------------------

      if(json.contains("orbit_specs")) {
        bool verbose = false;
        tree.read_custom_clusters_from_json(json["orbit_specs"], prim, prim.factor_group(), verbose);
      }

      return tree;
    }
    catch(...) {
      std::cerr << "Error in make_orbitree, with JSON input: \n";
      std::cerr << json << "\n";
      throw;
    }

  }

  //*******************************************************************************************
  /// \brief Expand a neighbor list to include neighborhood of another SiteOrbitree
  /**
   * Calculates the distances between all sites from each other
   * in every cluster of the orbitree, then stores this "delta"
   * in the primitive neighbor list if it hasn't been added already.
   * In the process, it updates the value
   * of nlist_ind of every site in the given orbitree to
   * hold the index of where said site was placed in the neighbor
   * list, relative to the first site of the cluster.
   */
  //*******************************************************************************************
  void expand_nlist(const Structure &prim, SiteOrbitree &tree, Array<UnitCellCoord> &nlist) {

    // We include all the primitive-cell neighbors by default at the front of the prim_nlist
    // This might be overkill in the long term (e.g., for large primitive cells where only a few sites have DoFs)
    // In that case, we may change tactics and construct a separate list of nlist_inds that correspond to prim-cell neighbors
    // We need this information in some form for calculating differential correlations
    if(!nlist.size()) {
      for(Index b = 0; b < prim.basis.size(); b++)
        nlist.push_back(UnitCellCoord(b, 0, 0, 0));
    }

    //For each site we encounter we access the appropriate slot in the neighbor list and append all other sites
    //to it in the form of UnitCellCoords

    Array<Index> clust_nlist_inds;

    //branches
    for(Index i = 0; i < tree.size(); i++) {
      //orbits
      for(Index j = 0; j < tree[i].size(); j++) {
        //clusters
        for(Index k = 0; k < tree[i][j].size(); k++) {

          clust_nlist_inds.resize(tree[i][j][k].size());

          //sites
          for(Index l = 0; l < tree[i][j][k].size(); l++) {
            //tuccl corresponds to a particular site we're looking at
            UnitCellCoord tuccl = prim.get_unit_cell_coord(tree[i][j][k][l]);
            //neighbor sites
            for(Index b = 0; b < tree[i][j][k].size(); b++) {
              //tuccb corresponds to a site that neighbors tuccl
              UnitCellCoord tuccb = prim.get_unit_cell_coord(tree[i][j][k][b]);
              UnitCellCoord delta = tuccb - tuccl;

              clust_nlist_inds[b] = nlist.find(delta);

              //If delta is not already in nlist, add it (the value of clust_nlist_inds[b] makes sense after the push_back)
              if(clust_nlist_inds[b] == nlist.size()) {
                nlist.push_back(delta);
              }
            }

            // we set the first set of nlist_inds as the current indices
            if(l == 0) {
              tree[i][j][k].set_nlist_inds(clust_nlist_inds);
            }

            // save any other unique ones that we find
            Index t(0);
            for(t = 0; t < tree[i][j][k].trans_nlists().size(); t++) {
              if(clust_nlist_inds.all_in(tree[i][j][k].trans_nlist(t)))
                break;
            }
            if(t == tree[i][j][k].trans_nlists().size())
              tree[i][j][k].add_trans_nlist(clust_nlist_inds);
          }
        }
      }
    }
  }

  //*******************************************************************************************
  /// \brief Print clexulator
  void print_clexulator(const Structure &prim,
                        SiteOrbitree &tree,
                        const Array<UnitCellCoord> &nlist,
                        std::string class_name,
                        std::ostream &stream) {

    DoFManager dof_manager;

    for(Index b = 0; b < prim.basis.size(); b++) {
      if(prim.basis[b].site_occupant().size() > 1) {
        dof_manager.add_dof(prim.basis[b].site_occupant().type_name());
        break;
      }
    }
    for(Index b = 0; b < prim.basis.size(); b++) {
      for(Index i = 0; i < prim.basis[i].displacement().size(); i++)
        dof_manager.add_dof(prim.basis[b].displacement()[i].type_name());
    }

    dof_manager.resize_neighborhood(nlist.size());

    // We can add more as needed
    dof_manager.register_dofs(tree);

    Index N_corr(tree.basis_set_size());
    std::stringstream private_def_stream, public_def_stream, interface_imp_stream, bfunc_imp_stream;

    std::string uclass_name;
    for(Index i = 0; i < class_name.size(); i++)
      uclass_name.push_back(std::toupper(class_name[i]));

    std::string indent(2, ' ');
    private_def_stream <<

                       indent << "  /// \\brief Clone the Clexulator\n" <<
                       indent << "  virtual " << class_name << "* _clone() const override {\n" <<
                       indent << "    return new " << class_name << "(*this);\n" <<
                       indent << "  }\n\n" <<

                       indent << "  // typedef for method pointers\n" <<
                       indent << "  typedef double (" << class_name << "::*BasisFuncPtr)() const;\n\n" <<

                       indent << "  // typedef for method pointers\n" <<
                       indent << "  typedef double (" << class_name << "::*DeltaBasisFuncPtr)(int, int) const;\n\n" <<

                       indent << "  // array of pointers to member functions for calculating basis functions\n" <<
                       indent << "  BasisFuncPtr m_orbit_func_list[" << N_corr << "];\n\n" <<

                       indent << "  // array of pointers to member functions for calculating flower functions\n" <<
                       indent << "  BasisFuncPtr m_flower_func_lists[" << prim.basis.size() << "][" << N_corr << "];\n\n" <<

                       /**** for separate 1D method pointer lists:
                       indent << "  BasisFuncPtr
                       for(Index i = 0; i < prim.basis.size(); i++) {
                         private_def_stream << " m_flower_func_at_" << i << "_list[" << N_corr << "]";
                         if(i + 1 < prim.basis.size())
                           private_def_stream << ',';
                       }
                       private_def_stream << ";\n\n" <<
                       **/

                       indent << "  // array of pointers to member functions for calculating DELTA flower functions\n" <<
                       indent << "  DeltaBasisFuncPtr m_delta_func_lists[" << prim.basis.size() << "][" << N_corr << "];\n\n";

    /**** for separate 1D method pointer lists:
    indent << "  DeltaBasisFuncPtr";

    for(Index i = 0; i < prim.basis.size(); i++) {
    private_def_stream << " m_delta_func_at_" << i << "_list[" << N_corr << "]";
    if(i + 1 < prim.basis.size())
    private_def_stream << ',';
    }
    private_def_stream << ";\n\n";
    **/

    dof_manager.print_clexulator_member_definitions(private_def_stream, prim, indent + "  ");

    // perhaps some more correlation calculating options here
    dof_manager.print_clexulator_private_method_definitions(private_def_stream, prim, indent + "  ");

    private_def_stream <<
                       indent << "  //default functions for basis function evaluation \n" <<
                       indent << "  double zero_func() const{ return 0.0;};\n" <<
                       indent << "  double zero_func(int,int) const{ return 0.0;};\n\n";

    public_def_stream <<
                      indent << "  " << class_name << "();\n\n" <<
                      indent << "  ~" << class_name << "();\n\n" <<

                      indent << "  /// \\brief Clone the " << class_name << "\n" <<
                      indent << "  std::unique_ptr<" << class_name << "> clone() const { \n" <<
                      indent << "    return std::unique_ptr<" << class_name << ">(_clone()); \n" <<
                      indent << "  }\n\n" <<

                      indent << "  /// \\brief Calculate contribution to global correlations from one unit cell\n" <<
                      indent << "  void calc_global_corr_contribution(double *corr_begin) const override;\n\n" <<

                      indent << "  /// \\brief Calculate contribution to select global correlations from one unit cell\n" <<
                      indent << "  void calc_restricted_global_corr_contribution(double *corr_begin, size_type const* ind_list_begin, size_type const* ind_list_end) const override;\n\n" <<

                      indent << "  /// \\brief Calculate point correlations about basis site 'b_index'\n" <<
                      indent << "  void calc_point_corr(int b_index, double *corr_begin) const override;\n\n" <<

                      indent << "  /// \\brief Calculate select point correlations about basis site 'b_index'\n" <<
                      indent << "  void calc_restricted_point_corr(int b_index, double *corr_begin, size_type const* ind_list_begin, size_type const* ind_list_end) const override;\n\n" <<

                      indent << "  /// \\brief Calculate the change in point correlations due to changing an occupant\n" <<
                      indent << "  void calc_delta_point_corr(int b_index, int occ_i, int occ_f, double *corr_begin) const override;\n\n" <<

                      indent << "  /// \\brief Calculate the change in select point correlations due to changing an occupant\n" <<
                      indent << "  void calc_restricted_delta_point_corr(int b_index, int occ_i, int occ_f, double *corr_begin, size_type const* ind_list_begin, size_type const* ind_list_end) const override;\n\n";

    dof_manager.print_clexulator_public_method_definitions(public_def_stream, prim, indent + "  ");


    //linear function index
    Index lf = 0, tlf;

    Array<FunctionVisitor *> labelers(dof_manager.get_function_label_visitors());
    //std::cout << "Initialized " << labelers.size() << " labelers \n";

    Array<std::string> orbit_method_names(N_corr);
    Array<Array<std::string> > flower_method_names(prim.basis.size(), Array<std::string>(N_corr));
    //Array< Array<Array<std::string> > > dflower_method_names(N_corr, Array<Array<std::string> >(prim.basis.size()));

    //this is very configuration-centric
    Array<Array<std::string> > dflower_method_names(prim.basis.size(), Array<std::string>(N_corr));

    // temporary storage for formula
    Array<std::string> formulae, tformulae;

    bool make_newline(false);

    //loop over orbits
    for(Index np = 0; np < tree.size(); np++) {
      for(Index no = 0; no < tree[np].size(); no++) {
        if(np == 0)
          bfunc_imp_stream <<
                           indent << "// Basis functions for empty cluster:\n";
        else {
          bfunc_imp_stream <<
                           indent << "/**** Basis functions for orbit " << np << ", " << no << "****\n";
          tree[np][no].prototype.print(bfunc_imp_stream, '\n');
          bfunc_imp_stream << "****/\n";
        }

        formulae = tree[np][no].orbit_function_cpp_strings(labelers);
        tlf = formulae.size();

        make_newline = false;
        for(Index nf = 0; nf < formulae.size(); nf++) {
          if(!formulae[nf].size())
            continue;
          make_newline = true;
          orbit_method_names[lf + nf] = "eval_bfunc_" + std::to_string(np) + "_" + std::to_string(no) + "_" + std::to_string(nf);
          private_def_stream <<
                             indent << "  double " << orbit_method_names[lf + nf] << "() const;\n";

          bfunc_imp_stream <<
                           indent << "double " << class_name << "::" << orbit_method_names[lf + nf] << "() const{\n" <<
                           indent << "  return " << formulae[nf] << ";\n" <<
                           indent << "}\n";
        }
        if(make_newline) {
          bfunc_imp_stream << '\n';
          private_def_stream << '\n';
        }
        make_newline = false;

        // loop over flowers (i.e., basis sites of prim)
        for(Index nb = 0; nb < prim.basis.size(); nb++) {
          formulae = tree[np][no].flower_function_cpp_strings(labelers, nb);
          for(Index nf = 0; nf < formulae.size(); nf++) {
            if(!formulae[nf].size())
              continue;
            make_newline = true;
            flower_method_names[nb][lf + nf] = "site_eval_at_" + std::to_string(nb) + "_bfunc_" + std::to_string(np) + "_" + std::to_string(no) + "_" + std::to_string(nf);
            private_def_stream <<
                               indent << "  double " << flower_method_names[nb][lf + nf] << "() const;\n";

            bfunc_imp_stream <<
                             indent << "double " << class_name << "::" << flower_method_names[nb][lf + nf] << "() const{\n" <<
                             indent << "  return " << formulae[nf] << ";\n" <<
                             indent << "}\n";

            //dflower_method_names[nb][lf + nf].resize(prim.basis[nb].occupant_basis().size());
          }
          if(make_newline) {
            bfunc_imp_stream << '\n';
            private_def_stream << '\n';
          }
          make_newline = false;

          // Very configuration-centric -> Find a way to move this block to OccupationDoFEnvironment:
          formulae.resize(formulae.size(), std::string());
          // loop over site basis functions
          for(Index nsbf = 0; nsbf < prim.basis[nb].occupant_basis().size(); nsbf++) {
            std::string delta_prefix = "(m_occ_func_" + std::to_string(nb) + "_" + std::to_string(nsbf) + "[occ_f] - m_occ_func_" + std::to_string(nb) + "_" + std::to_string(nsbf) + "[occ_i])";

            tformulae = tree[np][no].delta_occfunc_flower_function_cpp_strings(labelers, nb, nsbf);
            for(Index nf = 0; nf < tformulae.size(); nf++) {
              if(!tformulae[nf].size())
                continue;

              if(formulae[nf].size())
                formulae[nf] += " + ";

              formulae[nf] += delta_prefix;

              if(tformulae[nf] == "1" || tformulae[nf] == "(1)")
                continue;

              formulae[nf] += "*";
              formulae[nf] += tformulae[nf];
              //formulae[nf] += ")";
            }
          }
          for(Index nf = 0; nf < formulae.size(); nf++) {
            if(!formulae[nf].size())
              continue;
            make_newline = true;

            dflower_method_names[nb][lf + nf] = "delta_site_eval_at_" + std::to_string(nb) + "_bfunc_" + std::to_string(np) + "_" + std::to_string(no) + "_" + std::to_string(nf);
            private_def_stream <<
                               indent << "  double " << dflower_method_names[nb][lf + nf] << "(int occ_i, int occ_f) const;\n";

            bfunc_imp_stream <<
                             indent << "double " << class_name << "::" << dflower_method_names[nb][lf + nf] << "(int occ_i, int occ_f) const{\n" <<
                             indent << "  return " << formulae[nf] << ";\n" <<
                             indent << "}\n";
          }
          if(make_newline) {
            bfunc_imp_stream << '\n';
            private_def_stream << '\n';
          }
          make_newline = false;
        }
        // \End Configuration specific part

        lf += tlf;
      }
    }//Finished writing method definitions and implementations for basis functions

    //clean up:
    for(Index nl = 0; nl < labelers.size(); nl++)
      delete labelers[nl];
    labelers.clear();


    // Write constructor
    interface_imp_stream <<
                         indent << class_name << "::" << class_name << "() :\n" <<
                         indent << "  Clexulator_impl::Base(" << nlist.size() << ", " << N_corr << ") {\n";

    dof_manager.print_to_clexulator_constructor(interface_imp_stream, prim, indent + "  ");

    for(Index nf = 0; nf < orbit_method_names.size(); nf++) {
      if(orbit_method_names[nf].size() == 0)
        interface_imp_stream <<
                             indent << "  m_orbit_func_list[" << nf << "] = &" << class_name << "::zero_func;\n";
      else
        interface_imp_stream <<
                             indent << "  m_orbit_func_list[" << nf << "] = &" << class_name << "::" << orbit_method_names[nf] << ";\n";
    }
    interface_imp_stream << "\n\n";

    for(Index nb = 0; nb < flower_method_names.size(); nb++) {
      for(Index nf = 0; nf < flower_method_names[nb].size(); nf++) {
        if(flower_method_names[nb][nf].size() == 0)
          interface_imp_stream <<
                               indent << "  m_flower_func_lists[" << nb << "][" << nf << "] = &" << class_name << "::zero_func;\n";
        else
          interface_imp_stream <<
                               indent << "  m_flower_func_lists[" << nb << "][" << nf << "] = &" << class_name << "::" << flower_method_names[nb][nf] << ";\n";
      }
      interface_imp_stream << "\n\n";
    }

    for(Index nb = 0; nb < dflower_method_names.size(); nb++) {
      for(Index nf = 0; nf < dflower_method_names[nb].size(); nf++) {
        if(dflower_method_names[nb][nf].size() == 0)
          interface_imp_stream <<
                               indent << "  m_delta_func_lists[" << nb << "][" << nf << "] = &" << class_name << "::zero_func;\n";
        else
          interface_imp_stream <<
                               indent << "  m_delta_func_lists[" << nb << "][" << nf << "] = &" << class_name << "::" << dflower_method_names[nb][nf] << ";\n";
      }
      interface_imp_stream << "\n\n";
    }
    interface_imp_stream <<
                         indent << "}\n\n";

    // Write destructor

    interface_imp_stream <<
                         indent << class_name << "::~" << class_name << "(){\n" <<

                         indent << "  //nothing here for now\n" <<

                         indent << "}\n\n";

    // Write evaluation methods

    interface_imp_stream <<
                         indent << "/// \\brief Calculate contribution to global correlations from one unit cell\n" <<
                         indent << "void " << class_name << "::calc_global_corr_contribution(double *corr_begin) const {\n" <<
                         indent << "  for(size_type i=0; i<corr_size(); i++){\n" <<
                         indent << "    *(corr_begin+i) = (this->*m_orbit_func_list[i])();\n" <<
                         indent << "  }\n" <<
                         indent << "}\n\n" <<

                         indent << "/// \\brief Calculate contribution to select global correlations from one unit cell\n" <<
                         indent << "void " << class_name << "::calc_restricted_global_corr_contribution(double *corr_begin, size_type const* ind_list_begin, size_type const* ind_list_end) const {\n" <<
                         indent << "  for(; ind_list_begin<ind_list_end; ind_list_begin++){\n" <<
                         indent << "    *(corr_begin+*ind_list_begin) = (this->*m_orbit_func_list[*ind_list_begin])();\n" <<
                         indent << "  }\n" <<
                         indent << "}\n\n" <<

                         indent << "/// \\brief Calculate point correlations about basis site 'b_index'\n" <<
                         indent << "void " << class_name << "::calc_point_corr(int b_index, double *corr_begin) const {\n" <<
                         indent << "  for(size_type i=0; i<corr_size(); i++){\n" <<
                         indent << "    *(corr_begin+i) = (this->*m_flower_func_lists[b_index][i])();\n" <<
                         indent << "  }\n" <<
                         indent << "}\n\n" <<

                         indent << "/// \\brief Calculate select point correlations about basis site 'b_index'\n" <<
                         indent << "void " << class_name << "::calc_restricted_point_corr(int b_index, double *corr_begin, size_type const* ind_list_begin, size_type const* ind_list_end) const {\n" <<
                         indent << "  for(; ind_list_begin<ind_list_end; ind_list_begin++){\n" <<
                         indent << "    *(corr_begin+*ind_list_begin) = (this->*m_flower_func_lists[b_index][*ind_list_begin])();\n" <<
                         indent << "  }\n" <<
                         indent << "}\n\n" <<

                         indent << "/// \\brief Calculate the change in point correlations due to changing an occupant\n" <<
                         indent << "void " << class_name << "::calc_delta_point_corr(int b_index, int occ_i, int occ_f, double *corr_begin) const {\n" <<
                         indent << "  for(size_type i=0; i<corr_size(); i++){\n" <<
                         indent << "    *(corr_begin+i) = (this->*m_delta_func_lists[b_index][i])(occ_i, occ_f);\n" <<
                         indent << "  }\n" <<
                         indent << "}\n\n" <<

                         indent << "/// \\brief Calculate the change in select point correlations due to changing an occupant\n" <<
                         indent << "void " << class_name << "::calc_restricted_delta_point_corr(int b_index, int occ_i, int occ_f, double *corr_begin, size_type const* ind_list_begin, size_type const* ind_list_end) const {\n" <<
                         indent << "  for(; ind_list_begin<ind_list_end; ind_list_begin++){\n" <<
                         indent << "    *(corr_begin+*ind_list_begin) = (this->*m_delta_func_lists[b_index][*ind_list_begin])(occ_i, occ_f);\n" <<
                         indent << "  }\n" <<
                         indent << "}\n\n";


    // PUT EVERYTHING TOGETHER
    stream <<
           "#include <cstddef>\n" <<
           "#include \"casm/clex/Clexulator.hh\"\n" <<
           "\n\n\n" <<
           "/****** CLEXULATOR CLASS FOR PRIM ******" << std::endl;

    prim.print(stream);

    stream <<
           "**/\n\n\n" <<

           "/// \\brief Returns a Clexulator_impl::Base* owning a " << class_name << "\n" <<
           "extern \"C\" CASM::Clexulator_impl::Base* make_" + class_name << "();\n\n" <<

           "namespace CASM {\n\n" <<


           indent << "class " << class_name << " : public Clexulator_impl::Base {\n\n" <<

           indent << "public:\n\n" <<
           public_def_stream.str() << "\n" <<

           indent << "private:\n\n" <<
           private_def_stream.str() << "\n" <<

           indent << "};\n\n" << // close class definition

           indent <<

           "//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n" <<

           interface_imp_stream.str() <<
           bfunc_imp_stream.str() <<
           "}\n\n\n" <<      // close namespace

           "extern \"C\" {\n" <<
           indent << "/// \\brief Returns a Clexulator_impl::Base* owning a " << class_name << "\n" <<
           indent << "CASM::Clexulator_impl::Base* make_" + class_name << "() {\n" <<
           indent << "  return new CASM::" + class_name + "();\n" <<
           indent << "}\n\n" <<
           "}\n" <<

           "\n";
    // EOF

    return;
  }


}

