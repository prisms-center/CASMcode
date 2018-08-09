#include "casm/clex/ChemicalReference.hh"
#include "casm/clex/ConfigIO.hh"
#include "casm/clex/Configuration.hh"
#include "casm/clex/ConfigIterator.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/CompositionConverter.hh"

#include "casm/misc/algorithm.hh"

namespace CASM {

  /// \brief Construct using the results of n(config), and e(config)
  ///
  /// The element of n(config) corresponding to vacancies is ignored, the
  /// rest are mapped using 'species_num[molecule_name[i]] = n(i)'
  ChemicalReferenceState::ChemicalReferenceState(
    const Configuration &config,
    std::function<Eigen::VectorXd(Configuration)> n,
    std::function<double(Configuration)> e) {

    auto names = config.get_prim().get_struc_molecule_name();
    auto vec = n(config);

    if(vec.size() != names.size()) {
      std::cerr << "Error constructing ChemicalReferenceState: Number of species != vec.size()\n";
      std::cerr << "  get_struc_molecule_name: " << names << "\n";
      std::cerr << "  n(config): " << vec << std::endl;
      throw std::runtime_error("Error constructing ChemicalReferenceState: Number of species != vec.size()");
    }

    for(int i = 0; i < names.size(); ++i) {
      if(!is_vacancy(names[i])) {
        species_num[names[i]] = vec(i);
      }
    }
    energy_per_species = e(config);
  }


  // --- ChemicalReference implementations -----------

  const std::string ChemicalReference::Name = "chem_ref";

  const std::string ChemicalReference::Desc = "Returns a reference energy as interpolated via a composition-energy hyperplane.";

  namespace {

    Eigen::MatrixXd _species_frac_matrix(const PrimClex &primclex,
                                         const std::vector<std::string> &ref_config) {
      Eigen::MatrixXd _N(primclex.composition_axes().components().size(), ref_config.size());
      for(int i = 0; i < ref_config.size(); ++i) {
        _N.col(i) = species_frac(primclex.configuration(ref_config[i]));
      }
      return _N;
    }

    Eigen::MatrixXd _species_frac_space(const Eigen::MatrixXd &_N) {
      // The input space
      Eigen::MatrixXd N(_N.rows(), _N.cols() - 1);
      for(int i = 0; i < N.cols(); i++) {
        N.col(i) = _N.col(i + 1) - _N.col(0);
      }
      return N;
    }

    int _rank(const Eigen::MatrixXd &N,
              double lin_alg_tol) {
      auto Qr = N.transpose().fullPivHouseholderQr();
      Qr.setThreshold(lin_alg_tol);
      return Qr.rank();
    }

    // Eliminate ref config until the rank of the defined species_frac space
    // is the same as the number of ref config.
    void _prune_ref_config(const PrimClex &primclex,
                           std::vector<std::string> &ref_config,
                           double lin_alg_tol) {
      // Each column of _N is the species_frac of the corresponding ref config
      Eigen::MatrixXd _N = _species_frac_matrix(primclex, ref_config);

      // Contains vectors spanning the space defined by the ref config
      // N has _N.cols() - 1 columns, with N.col(i) being _N.col(i) - _N.col(0)
      Eigen::MatrixXd N = _species_frac_space(_N);

      while(_rank(N, lin_alg_tol) != N.cols()) {

        // remove ref_config whose summed distance in species_frac space to all
        // others is the minimum
        double min_dist = std::numeric_limits<double>::max();
        double min_dist_ref = -1;
        for(int i = 0; i < ref_config.size(); ++i) {
          double tot_dist = 0.0;
          for(int j = 0; j < ref_config.size(); ++j) {
            if(i != j) {
              tot_dist += (_N.col(i) - _N.col(j)).norm();
            }
          }
          if(tot_dist < min_dist) {
            min_dist_ref = i;
            min_dist = tot_dist;
          }
        }

        ref_config.erase(ref_config.begin() + min_dist_ref);
        _N = _species_frac_matrix(primclex, ref_config);
        N = _species_frac_space(_N);
      }

    }

  }

  /// \brief Convert a set of ChemicalReferenceState to a hyperplane, including checks
  Eigen::VectorXd ChemicalReference::_calc_hyperplane(
    const Structure &prim,
    const std::vector<std::string> &struc_mol_name,
    Eigen::MatrixXd _N,
    Eigen::VectorXd E,
    double tol) {

    // --- convert input compositions to atom_frac

    Index Va_index = find_index_if(struc_mol_name, [ = ](const std::string & str) {
      return is_vacancy(str);
    });
    bool has_Va = (Va_index != struc_mol_name.size());

    if(has_Va) {
      _N.row(Va_index) = Eigen::VectorXd::Zero(_N.cols());
    }
    for(int i = 0; i < _N.cols(); i++) {
      _N.col(i) /= _N.col(i).sum();
    }

    // --- check that the input space is full rank (excluding Va) --------------

    // The input space
    Eigen::MatrixXd N = _species_frac_space(_N);

    int r = _rank(N, tol);

    if(r != N.cols()) {
      std::cerr << "Error in ChemicalReference::hyperplane " << std::endl;
      std::cerr << "Input space (column vectors of atom_frac):\n" << N << std::endl;
      std::cerr << "Rows correspond to: " << jsonParser(struc_mol_name) << std::endl;
      std::cerr << "Input space rank: " << r << std::endl;
      throw std::runtime_error("Error in ChemicalReference::hyperplane: Too many reference states specified");
    }

    //  --- check that the input space spans the prim space --------------------

    // get the prim composition space (column vectors are comp_n)
    Eigen::MatrixXd C = composition_space(prim, tol);

    // get Molecule allowed in prim, and how many there are
    std::vector<Molecule> struc_mol = prim.get_struc_molecule();
    for(int i = 0; i < struc_mol.size(); i++) {
      if(!is_molecule_name(struc_mol[i], struc_mol_name[i])) {
        std::cerr << "Error in ChemicalReference::hyperplane " << std::endl;
        std::cerr << "Initial struc_mol_name must be in same order as prim.get_struc_molecule()" << std::endl;
        throw std::runtime_error("Error in ChemicalReference::hyperplane: Molecule name mismatch");
      }
    }
    Index prim_N_mol = struc_mol.size();

    // ensure that there is a solution X for:
    //             C = N.topRows(prim_N_mol) * X
    //   (prim_space = input_space involving prim species * X)
    Eigen::MatrixXd X = N.topRows(prim_N_mol).fullPivHouseholderQr().solve(C);

    double relative_error = (N.topRows(prim_N_mol) * X - C).norm() / C.norm();

    if(relative_error > tol) {
      std::cerr << "Error in ChemicalReference::hyperplane " << std::endl;
      std::cerr << "Input space does not span the composition space of your prim." << std::endl;

      std::cerr << "Input space (column vectors in atom_frac space):\n" << N.topRows(prim_N_mol) << std::endl;
      std::cerr << "End members:\n" << end_members(prim) << std::endl;
      std::cerr << "Prim space (column vectors in atom_frac space):\n" << C << std::endl;
      std::cerr << "Rows correspond to: " << jsonParser(struc_mol_name) << std::endl;
      std::cerr << "X, prim_space = input_space*X: \n" << X << std::endl;
      std::cerr << "input_space*X: \n" << N.topRows(prim_N_mol)*X << std::endl;
      std::cerr << "relative_error: " << relative_error << std::endl;
      std::cerr << "tol: " << tol << std::endl;

      throw std::runtime_error("Error in ChemicalReference::hyperplane: Input space does not span prim space");
    }


    // --- solve N.transpose() * P = E, for P, the hyperplane reference --------

    Eigen::VectorXd P = _N.transpose().fullPivHouseholderQr().solve(E);

    relative_error = (_N.transpose() * P - E).norm() / E.norm();

    if(relative_error > tol) {
      std::cerr << "Error in ChemicalReference::hyperplane " << std::endl;
      std::cerr << "Could not solve for hyperplane reference." << std::endl;

      std::cerr << "Input space (column vectors in atom_frac space), N:\n" << N << std::endl;
      std::cerr << "Rows correspond to: " << jsonParser(struc_mol_name) << std::endl;
      std::cerr << "Solve: _N.transpose()*P = E" << std::endl;
      std::cerr << "_N.transpose():\n" << _N.transpose() << std::endl;
      std::cerr << "Reference state energies, E:\n" << E << std::endl;
      std::cerr << "P:\n" << P.transpose() << std::endl;
      std::cerr << "relative_err: " << relative_error << std::endl;
      std::cerr << "tol: " << tol << std::endl;

      throw std::runtime_error("Error in ChemicalReference::hyperplane: Could not solve for hyperplane reference");
    }

    if(has_Va && !almost_zero(P(Va_index))) {
      std::cerr << "Error in ChemicalReference::hyperplane " << std::endl;
      std::cerr << "Non-zero pure Va reference: " << P.transpose() << std::endl;
      std::cerr << "Elements correspond to: " << jsonParser(struc_mol_name) << std::endl;

      throw std::runtime_error("Error in ChemicalReference::hyperplane: Input space does not span prim space");
    }

    return P.head(prim_N_mol);
  }

  /// \brief Automatically set ChemicalReference using calculated Configurations
  ///        with 'extreme' compositions
  ChemicalReference auto_chemical_reference(const PrimClex &primclex, double lin_alg_tol) {

    auto closest_calculated_config = [&](const Eigen::VectorXd & target) {

      // return name of Configuration with param_comp closest to target_param_comp
      //   tie break goes to first Configuration with fewest atoms
      //
      //   must be Configurations for which the relaxed_energy has been calculated
      auto begin = primclex.config_begin();
      auto end = primclex.config_end();
      auto res = end;
      double close_dist = std::numeric_limits<double>::max();
      for(auto it = begin; it != end; ++it) {
        if(!it->calc_properties().contains("relaxed_energy")) {
          continue;
        }
        double curr_dist = (target - comp(*it)).norm();
        if((res == end) ||
           (almost_equal(curr_dist, close_dist, TOL) && it->size() < res->size()) ||
           (curr_dist < close_dist)) {
          res = it;
          close_dist = curr_dist;
        }
      }

      if(res == end) {
        std::cerr << "Error in auto_chemical_reference: Could not find any configurations\n";
        throw std::runtime_error("Error in auto_chemical_reference: Could not find any configurations");
      }
      return res->name();
    };



    // get number of reference states needed minus one... we'll also look for the 'origin'
    int Naxes = primclex.composition_axes().independent_compositions();
    Eigen::VectorXd target = Eigen::VectorXd::Zero(Naxes);

    // get 'origin' configuration
    std::vector<std::string> ref_config;
    ref_config.push_back(closest_calculated_config(target));

    for(int i = 0; i < Naxes; i++) {

      target(i) = 1.0;

      // get end member configurations
      std::string configname = closest_calculated_config(target);

      // ensure no repeats
      if(std::find(ref_config.begin(), ref_config.end(), configname) != ref_config.end()) {
        std::cerr << "Error in auto_chemical_reference: Could not find enough calculated configurations\n";
        throw std::runtime_error("Error in auto_chemical_reference: Could not find enough calculated configurations");
      }

      ref_config.push_back(configname);

      target(i) = 0.0;
    }

    // make sure there are the right number of references in species_frac space
    // (may be 1 extra if pure Va configuration is possible)
    _prune_ref_config(primclex, ref_config, lin_alg_tol);

    // construct ChemicalReferenceState
    std::vector<ChemicalReferenceState> ref_states;
    for(auto it = ref_config.begin(); it != ref_config.end(); ++it) {
      ref_states.emplace_back(primclex.configuration(*it),
                              ConfigIO::SpeciesFrac(),
                              ConfigIO::relaxed_energy_per_species());
    }

    return ChemicalReference(primclex.get_prim(), ref_states.begin(), ref_states.end(), lin_alg_tol);
  }

}
