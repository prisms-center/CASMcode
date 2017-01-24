#include "casm/CASM_global_definitions.hh"
#include "casm/clex/ConfigDoF.hh"
#include "casm/symmetry/PermuteIterator.hh"
#include "casm/clex/Correlation.hh"
#include "casm/clex/Clexulator.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/Supercell.hh"
#include "casm/clex/NeighborList.hh"


namespace CASM {

  ConfigDoF::ConfigDoF(Index _N, double _tol) :
    m_N(_N),
    m_occupation(m_N),
    m_deformation(Eigen::Matrix3d::Identity()),
    m_has_deformation(false),
    m_tol(_tol) {
  }

  //*******************************************************************************

  ConfigDoF::ConfigDoF(const Array<int> &_occ, double _tol):
    m_N(_occ.size()),
    m_occupation(_occ),
    m_deformation(Eigen::Matrix3d::Identity()),
    m_has_deformation(false),
    m_tol(_tol) {
  }

  //*******************************************************************************

  ConfigDoF::ConfigDoF(const Array<int> &_occ, const Eigen::MatrixXd &_disp, const Eigen::Matrix3d &_deformation, double _tol):
    m_N(_occ.size()),
    m_occupation(_occ),
    m_displacement(_disp),
    m_deformation(_deformation),
    m_has_deformation(false),
    m_tol(_tol) {
    if(size() != occupation().size() || size() != displacement().cols()) {
      std::cerr << "CRITICAL ERROR: Attempting to initialize ConfigDoF with data structures of incompatible size: \n"
                << "                occupation().size() = " << occupation().size() << " and displacement().size() = " << displacement().size() << "\n"
                << "                Exiting...\n";
      exit(1);
    }

  }

  //*******************************************************************************

  void ConfigDoF::clear() {
    m_N = 0;
    clear_occupation();
    clear_displacement();
    clear_deformation();
  }

  void ConfigDoF::clear_occupation() {
    _occupation().clear();
  }

  void ConfigDoF::clear_displacement() {
    _displacement() = Eigen::MatrixXd();
  }

  void ConfigDoF::clear_deformation() {
    _deformation() = Eigen::Matrix3d::Identity();
    m_has_deformation = false;
  }

  //*******************************************************************************

  void ConfigDoF::swap(ConfigDoF &RHS) {
    _deformation().swap(RHS._deformation());
    _occupation().swap(RHS._occupation());
    _displacement().swap(RHS._displacement());
    std::swap(m_N, RHS.m_N);
    std::swap(m_tol, RHS.m_tol);
    std::swap(m_has_deformation, RHS.m_has_deformation);
  }

  //*******************************************************************************

  void ConfigDoF::set_occupation(const Array<int> &new_occupation) {
    if(!m_N)
      m_N = new_occupation.size();

    if(m_N != new_occupation.size()) {
      std::cerr << "CRITICAL ERROR: In ConfigDoF::set_occupation(), attempting to set occupation to size " << new_occupation.size() << ",\n"
                << "                which does not match initialized size of ConfigDoF -> " << size() << "\n"
                << "                Exiting...\n";
      assert(0);
      exit(1);
    }

    _occupation() = new_occupation;

  }

  //*******************************************************************************

  void ConfigDoF::set_displacement(const displacement_matrix_t &new_displacement) {
    if(!m_N)
      m_N = new_displacement.cols();

    if(m_N != new_displacement.cols()) {
      std::cerr << "CRITICAL ERROR: In ConfigDoF::set_displacement(), attempting to set displacement to size " << new_displacement.cols() << ",\n"
                << "                which does not match initialized size of ConfigDoF -> " << size() << "\n"
                << "                Exiting...\n";
      assert(0);
      exit(1);
    }

    _displacement() = new_displacement;
  }

  //*******************************************************************************

  void ConfigDoF::set_deformation(const Eigen::Matrix3d &new_deformation) {
    m_has_deformation = true;
    _deformation() = new_deformation;
  }

  //*******************************************************************************

  jsonParser &ConfigDoF::to_json(jsonParser &json) const {
    json = jsonParser::object();
    if(occupation().size())
      json["occupation"] = occupation();
    if(displacement().size())
      json["displacement"] = displacement();
    if(has_deformation())
      json["deformation"] = deformation();

    return json;
  }

  //*******************************************************************************
  void ConfigDoF::from_json(const jsonParser &json) {
    clear();
    json.get_if(_occupation(), "occupation");
    m_N = occupation().size();

    json.get_if(_displacement(), "displacement");
    if(displacement().cols() && size() && displacement().cols() != size()) {
      std::cerr << "CRITICAL ERROR: In ConfigDoF::from_json(), parsing displacement having size " << displacement().cols() << ",\n"
                << "                which does not match initialized size of ConfigDoF -> " << size() << "\n"
                << "                Exiting...\n";
      assert(0);
      exit(1);
    }

    if(json.contains("deformation")) {
      CASM::from_json(_deformation(), json["deformation"]);
      m_has_deformation = true;
    }
  }

  //*******************************************************************************

  jsonParser &to_json(const ConfigDoF &value, jsonParser &json) {
    return value.to_json(json);
  }

  //*******************************************************************************

  void from_json(ConfigDoF &value, const jsonParser &json) {
    value.from_json(json);
  }

  //*******************************************************************************

  // Calculate transformed ConfigDoF from PermuteIterator via
  //   transformed_configdof = permute_iterator * configdof
  ConfigDoF &apply(const PermuteIterator &it, ConfigDoF &dof) {
    Eigen::Matrix3d fg_cart_op = it.sym_op().matrix();
    if(dof.has_deformation())
      dof.set_deformation(fg_cart_op * dof.deformation() * fg_cart_op.transpose());

    Permutation tperm(it.combined_permute());
    if(dof.occupation().size())
      dof.set_occupation(tperm * dof.occupation());

    if(dof.displacement().cols()) {
      Eigen::MatrixXd new_disp = fg_cart_op * dof.displacement();
      dof.set_displacement(Eigen::MatrixXd(3, dof.size()));
      for(Index i = 0; i < dof.size(); i++)
        dof.disp(i) = new_disp.col(tperm[i]);
    }
    return dof;
  }

  //*******************************************************************************

  void swap(ConfigDoF &A, ConfigDoF &B) {
    A.swap(B);
  }

  /// \brief Returns correlations using 'clexulator'. Supercell needs a correctly populated neighbor list.
  Correlation correlations(const ConfigDoF &configdof, const Supercell &scel, Clexulator &clexulator) {

    //Size of the supercell will be used for normalizing correlations to a per primitive cell value
    int scel_vol = scel.volume();

    Correlation correlations = Correlation::Zero(clexulator.corr_size());

    //Inform Clexulator of the bitstring

    //TODO: This will probably get more complicated with displacements and stuff
    clexulator.set_config_occ(configdof.occupation().begin());
    //mc_clexor.set_config_disp(mc_confdof.m_displacements.begin());   //or whatever
    //mc_clexor.set_config_strain(mc_confdof.m_strain.begin());   //or whatever

    //Holds contribution to global correlations from a particular neighborhood
    std::vector<double> tcorr(clexulator.corr_size(), 0.0);
    //std::vector<double> corr(clexulator.corr_size(), 0.0);

    for(int v = 0; v < scel_vol; v++) {

      //Point the Clexulator to the right neighborhood
      clexulator.set_nlist(scel.nlist().sites(v).data());

      //Fill up contributions
      clexulator.calc_global_corr_contribution(&tcorr[0]);

      //Add contributions to total correlations
      for(int i = 0; i < tcorr.size(); i++) {
        correlations[i] += tcorr[i];
      }

    }

    // normalize by supercell volume
    for(int i = 0; i < clexulator.corr_size(); i++) {
      correlations[i] /= (double) scel_vol;
    }

    return correlations;
  }

  /// \brief Returns correlations using 'clexulator'. Supercell needs a correctly populated neighbor list.
  Eigen::VectorXd correlations_vec(const ConfigDoF &configdof, const Supercell &scel, Clexulator &clexulator) {

    //Size of the supercell will be used for normalizing correlations to a per primitive cell value
    int scel_vol = scel.volume();

    Eigen::VectorXd correlations = Eigen::VectorXd::Zero(clexulator.corr_size());

    //Inform Clexulator of the bitstring

    //TODO: This will probably get more complicated with displacements and stuff
    clexulator.set_config_occ(configdof.occupation().begin());
    //mc_clexor.set_config_disp(mc_confdof.m_displacements.begin());   //or whatever
    //mc_clexor.set_config_strain(mc_confdof.m_strain.begin());   //or whatever

    //Holds contribution to global correlations from a particular neighborhood
    Eigen::VectorXd tcorr = correlations;
    //std::vector<double> corr(clexulator.corr_size(), 0.0);

    for(int v = 0; v < scel_vol; v++) {

      //Point the Clexulator to the right neighborhood
      clexulator.set_nlist(scel.nlist().sites(v).data());

      //Fill up contributions
      clexulator.calc_global_corr_contribution(tcorr.data());

      correlations += tcorr;

    }

    correlations /= (double) scel_vol;

    return correlations;
  }

  /// \brief Returns num_each_molecule[ molecule_type], where 'molecule_type' is ordered as Structure::get_struc_molecule()
  ReturnArray<int> get_num_each_molecule(const ConfigDoF &configdof, const Supercell &scel) {

    // [basis_site][site_occupant_index]
    auto convert = get_index_converter(scel.get_prim(), scel.get_prim().get_struc_molecule());

    // create an array to count the number of each molecule
    Array<int> num_each_molecule(scel.get_prim().get_struc_molecule().size(), 0);

    // count the number of each molecule
    for(Index i = 0; i < configdof.size(); i++) {
      num_each_molecule[ convert[ scel.get_b(i) ][ configdof.occ(i)] ]++;
    }

    return num_each_molecule;
  }

  /// \brief Returns num_each_molecule(molecule_type), where 'molecule_type' is ordered as Structure::get_struc_molecule()
  Eigen::VectorXi get_num_each_molecule_vec(const ConfigDoF &configdof, const Supercell &scel) {

    // [basis_site][site_occupant_index]
    auto convert = get_index_converter(scel.get_prim(), scel.get_prim().get_struc_molecule());

    // create an array to count the number of each molecule
    Eigen::VectorXi num_each_molecule = Eigen::VectorXi::Zero(scel.get_prim().get_struc_molecule().size());

    // count the number of each molecule
    for(Index i = 0; i < configdof.size(); i++) {
      num_each_molecule(convert[ scel.get_b(i) ][ configdof.occ(i)])++;
    }

    return num_each_molecule;
  }

  /// \brief Returns comp_n, the number of each molecule per primitive cell, ordered as Structure::get_struc_molecule()
  Eigen::VectorXd comp_n(const ConfigDoF &configdof, const Supercell &scel) {
    return get_num_each_molecule_vec(configdof, scel).cast<double>() / scel.volume();
  }

}

