#include "casm/clex/ConfigDoF.hh"

#include "casm/CASM_global_definitions.hh"
#include "casm/casm_io/json_io/container.hh"
#include "casm/symmetry/PermuteIterator.hh"
#include "casm/crystallography/Structure.hh"
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

  ConfigDoF::ConfigDoF(const std::vector<int> &_occ, double _tol):
    m_N(_occ.size()),
    m_occupation(_occ),
    m_deformation(Eigen::Matrix3d::Identity()),
    m_has_deformation(false),
    m_tol(_tol) {
  }

  //*******************************************************************************

  ConfigDoF::ConfigDoF(const std::vector<int> &_occ, const Eigen::MatrixXd &_disp, const Eigen::Matrix3d &_deformation, double _tol):
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
    occupation().clear();
  }

  void ConfigDoF::clear_displacement() {
    displacement() = Eigen::MatrixXd();
  }

  void ConfigDoF::clear_deformation() {
    deformation() = Eigen::Matrix3d::Identity();
    m_has_deformation = false;
  }

  //*******************************************************************************

  void ConfigDoF::swap(ConfigDoF &RHS) {
    deformation().swap(RHS.deformation());
    occupation().swap(RHS.occupation());
    displacement().swap(RHS.displacement());
    std::swap(m_N, RHS.m_N);
    std::swap(m_tol, RHS.m_tol);
    std::swap(m_has_deformation, RHS.m_has_deformation);
  }

  //*******************************************************************************

  void ConfigDoF::set_occupation(const std::vector<int> &new_occupation) {
    if(!m_N)
      m_N = new_occupation.size();

    if(m_N != new_occupation.size()) {
      std::cerr << "CRITICAL ERROR: In ConfigDoF::set_occupation(), attempting to set occupation to size " << new_occupation.size() << ",\n"
                << "                which does not match initialized size of ConfigDoF -> " << size() << "\n"
                << "                Exiting...\n";
      assert(0);
      exit(1);
    }

    occupation() = new_occupation;

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

    displacement() = new_displacement;
  }

  //*******************************************************************************

  void ConfigDoF::set_deformation(const Eigen::Matrix3d &new_deformation) {
    m_has_deformation = true;
    deformation() = new_deformation;
  }

  //*******************************************************************************

  // Calculate transformed ConfigDoF from PermuteIterator via
  //   *this = permute_iterator * (*this)
  ConfigDoF &ConfigDoF::apply_sym(const PermuteIterator &it) {
    Eigen::Matrix3d fg_cart_op = it.sym_op().matrix();
    if(has_deformation()) {
      set_deformation(fg_cart_op * deformation() * fg_cart_op.transpose());
    }
    Permutation tperm(it.combined_permute());
    if(occupation().size()) {
      set_occupation(tperm * occupation());
    }
    if(displacement().cols()) {
      Eigen::MatrixXd new_disp = fg_cart_op * displacement();
      set_displacement(Eigen::MatrixXd(3, size()));
      for(Index i = 0; i < size(); i++)
        disp(i) = new_disp.col(tperm[i]);
    }

    return *this;
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
    json.get_if(occupation(), "occupation");
    m_N = occupation().size();

    json.get_if(displacement(), "displacement");
    if(displacement().cols() && size() && displacement().cols() != size()) {
      std::cerr << "CRITICAL ERROR: In ConfigDoF::from_json(), parsing displacement having size " << displacement().cols() << ",\n"
                << "                which does not match initialized size of ConfigDoF -> " << size() << "\n"
                << "                Exiting...\n";
      assert(0);
      exit(1);
    }

    if(json.contains("deformation")) {
      CASM::from_json(deformation(), json["deformation"]);
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

  void swap(ConfigDoF &A, ConfigDoF &B) {
    A.swap(B);
  }

  /// \brief Returns correlations using 'clexulator'. Supercell needs a correctly populated neighbor list.
  Eigen::VectorXd correlations(const ConfigDoF &configdof, const Supercell &scel, Clexulator &clexulator) {

    //Size of the supercell will be used for normalizing correlations to a per primitive cell value
    int scel_vol = scel.volume();

    Eigen::VectorXd correlations = Eigen::VectorXd::Zero(clexulator.corr_size());

    //Inform Clexulator of the bitstring

    //Holds contribution to global correlations from a particular neighborhood
    Eigen::VectorXd tcorr = correlations;
    //std::vector<double> corr(clexulator.corr_size(), 0.0);

    for(int v = 0; v < scel_vol; v++) {
      //Fill up contributions
      clexulator.calc_global_corr_contribution(configdof,
                                               scel.nlist().sites(v).data(),
                                               tcorr.data());

      correlations += tcorr;
    }

    correlations /= (double) scel_vol;

    return correlations;
  }

  /// \brief Returns num_each_molecule(molecule_type), where 'molecule_type' is ordered as Structure::struc_molecule()
  Eigen::VectorXi num_each_molecule(const ConfigDoF &configdof, const Supercell &scel) {

    // [basis_site][site_occupant_index]
    auto convert = make_index_converter(scel.prim(), scel.prim().struc_molecule());

    // create an array to count the number of each molecule
    Eigen::VectorXi num_each_molecule = Eigen::VectorXi::Zero(scel.prim().struc_molecule().size());

    // count the number of each molecule
    for(Index i = 0; i < configdof.size(); i++) {
      num_each_molecule(convert[ scel.sublat(i) ][ configdof.occ(i)])++;
    }

    return num_each_molecule;
  }

  /// \brief Returns comp_n, the number of each molecule per primitive cell, ordered as Structure::struc_molecule()
  Eigen::VectorXd comp_n(const ConfigDoF &configdof, const Supercell &scel) {
    return num_each_molecule(configdof, scel).cast<double>() / scel.volume();
  }

}

