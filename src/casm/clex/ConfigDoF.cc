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
    m_is_strained(false),
    m_tol(_tol) {
  }

  //*******************************************************************************

  ConfigDoF::ConfigDoF(const Array<int> &_occ, double _tol):
    m_N(_occ.size()),
    m_occupation(_occ),
    m_deformation(Eigen::Matrix3d::Identity()),
    m_is_strained(false),
    m_tol(_tol) {
  }

  //*******************************************************************************

  ConfigDoF::ConfigDoF(const Array<int> &_occ, const Eigen::MatrixXd &_disp, const Eigen::Matrix3d &_deformation, double _tol):
    m_N(_occ.size()),
    m_occupation(_occ),
    m_displacement(_disp),
    m_deformation(_deformation),
    m_is_strained(false),
    m_tol(_tol) {
    if(size() != occupation().size() || size() != displacement().cols()) {
      std::cerr << "CRITICAL ERROR: Attempting to initialize ConfigDoF with data structures of incompatible size: \n"
                << "                occupation().size() = " << occupation().size() << " and displacement().size() = " << displacement().size() << "\n"
                << "                Exiting...\n";
      exit(1);
    }

  }

  //*******************************************************************************

  bool ConfigDoF::operator==(const ConfigDoF &RHS)const {
    if(occupation() != RHS.occupation())
      return false;

    double _tol = max(tol(), RHS.tol());
    // floating-point comparison tolerance is _tol

    if(!has_displacement()) {
      if(RHS.has_displacement() && !almost_zero(RHS.displacement(), _tol)) {
        //std::cout << "FALSE 1: " << _tol << "\n" << RHS.displacement() << "\n";
        return false;
      }
    }
    //    has_displacement()==true;
    else if(!RHS.has_displacement()) {
      if(!almost_zero(displacement(), _tol)) {
        //std::cout << "FALSE 2: " << _tol << "\n" << displacement() << "\n";
        return false;
      }
    }
    //    has_displacement()==true && RHS.has_displacement==true;
    else if(!almost_equal(displacement(), RHS.displacement(), _tol)) {
      //std::cout << "FALSE 3: " << _tol << "\n" << RHS.displacement() << "\n\n" << displacement() << "\n";;
      return false;
    }

    if(is_strained() || RHS.is_strained()) {
      //Deformation defaults to identity--should work in all cases
      if(!Eigen::almost_equal(deformation(), RHS.deformation(), _tol)) {
        //std::cout << "FALSE 4: " << _tol << "\n" << RHS.deformation() << "\n\n" << deformation() << "\n";;
        return false;

      }
    }

    return true;
  }

  //*******************************************************************************

  void ConfigDoF::clear() {
    m_N = 0;
    _occupation().clear();
    _displacement() = Eigen::MatrixXd();
    _deformation() = Eigen::Matrix3d::Identity();
    m_is_strained = false;
  }

  //*******************************************************************************

  void ConfigDoF::swap(ConfigDoF &RHS) {
    _deformation().swap(RHS._deformation());
    _occupation().swap(RHS._occupation());
    _displacement().swap(RHS._displacement());
    std::swap(m_N, RHS.m_N);
    std::swap(m_tol, RHS.m_tol);
    std::swap(m_is_strained, RHS.m_is_strained);
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
    m_is_strained = true;
    _deformation() = new_deformation;
  }

  //*******************************************************************************

  bool ConfigDoF::is_primitive(PermuteIterator it_begin, double _tol) const {

    Index i, j;
    Index permute_ind;
    PermuteIterator it_next_fg(it_begin.begin_next_fg_op());

    displacement_matrix_t new_disp;

    ++it_begin; // skip very first operation (zero translation)

    // loop over all non-zero translations, trying to find one that leaves configdof unchanged
    // if such a translation exists, the configdof is not primitive
    for(; it_begin != it_next_fg; ++it_begin) {
      bool proceed_to_next = false;

      //check if translation maps occupation onto itself
      if(occupation().size()) {
        for(i = 0; i < size(); i++) {
          permute_ind = it_begin.permute_ind(i);
          if(occ(permute_ind) != occ(i)) {
            proceed_to_next = true;
            break;
          }
        }
      }

      if(proceed_to_next) {
        continue;
      }

      //check if translation maps displacement field onto itself
      if(displacement().cols()) {
        for(i = 0; i < size(); i++) {
          for(j = 0; j < 3; j++) {
            permute_ind = it_begin.permute_ind(i);
            if(almost_equal(disp(permute_ind)[j], disp(i)[j], _tol)) {
              proceed_to_next = true;
              break;
            }
          }
          if(proceed_to_next)
            break;
        }
      }

      // proceed_to_next==false implies we reached end of checks without finding a difference in the configuration
      if(!proceed_to_next)
        return false;
    }//\End Translation loop

    return true;

  };

  //*******************************************************************************

  std::vector<PermuteIterator> ConfigDoF::factor_group(PermuteIterator it_begin, PermuteIterator it_end, double _tol) const {
    std::vector<PermuteIterator> factorgroup;

    Index i, j;
    Index permute_ind;
    PermuteIterator it_next_fg;
    Eigen::Matrix3d fg_cart_op;
    Eigen::Matrix3d new_F;

    displacement_matrix_t new_disp;

    while(it_begin != it_end) {
      fg_cart_op = it_begin.sym_op().matrix();
      it_next_fg = it_begin.begin_next_fg_op();

      bool proceed_to_next = false;

      // check for canonical strain first, since it is fastest
      if(is_strained()) {
        new_F = fg_cart_op * deformation() * fg_cart_op.transpose();
        for(i = 0; i < 3; i++) {
          for(j = 0; j < 3; j++) {
            if(!almost_equal(new_F(i, j), deformation()(i, j), _tol)) {
              proceed_to_next = true;
              break;
            }
          }
          if(proceed_to_next)
            break;
        }
      }

      if(proceed_to_next) {
        it_begin = it_next_fg;
        continue;
      }


      // check for canonical displacement
      if(displacement().cols()) {
        new_disp = fg_cart_op * displacement();
      }

      while(it_begin != it_next_fg) {
        proceed_to_next = false;

        if(occupation().size()) {      // check for canonical occupation
          for(i = 0; i < size(); i++) {
            permute_ind = it_begin.permute_ind(i);
            if(occ(permute_ind) != occ(i)) {
              proceed_to_next = true;
              break;
            }
          }
        } //\End occupation check

        if(proceed_to_next) {
          ++it_begin;
          continue;
        }

        // check for canonical displacement
        if(displacement().cols()) {
          new_disp = fg_cart_op * displacement();

          for(i = 0; i < size(); i++) {
            for(j = 0; j < 3; j++) {
              permute_ind = it_begin.permute_ind(i);
              if(!almost_equal(new_disp(j, permute_ind), disp(i)[j], _tol)) {
                proceed_to_next = true;
                break;
              }
            }
            if(proceed_to_next)
              break;
          }
        }


        if(!proceed_to_next) {
          factorgroup.push_back(it_begin);
        }

        ++it_begin;


      }
      // at end of loop, we have it_begin = it_next_fg;
    }

    return factorgroup;

  };

  //*******************************************************************************
  // ROUTINES FOR GETTING/CHECKING CANONICAL FORM
  //*******************************************************************************
  /**
   *   The canonical form is the ConfigDoF, canon_config, such that any symmetrically-equivalent ConfigDoF, equiv_config, satisfies:
   *      * equiv_config.deformation()(i,j) < canon_config.deformation()(i,j)
   *        AND equiv_config.deformation()(k,l) == canon_config.deformation()(k,l)
   *            :: for [all 'l' when 'k'<'i'] and for ['l'<'j' when 'k'=='i']
   *      * IF equiv_config.deformation()(i,j) == canon_config.deformation()(i,j) for all 'i' and 'j'
   *        THEN equiv_config.occ(l) < canon_config(l)
   *             AND equiv_config.occ(m) == canon_config(m)
   *                 :: for [all 'm'<'l']
   *      * IF equiv_config.occ(m) == canon_config(m) for all 'm'
   *        THEN equiv_config.disp(s)[v] < equiv_config.disp(s)[v]
   *             AND equiv_config.disp(t)[w] < equiv_config.disp(t)[w]
   *                 :: for [all 'w' when 't'<'s'] and for ['w'<'v' when 't'=='s']
   *
   *   The canonical configuration is equivalent to forming a vector (config_vec) by concatenating:
   *      * deformation(), unrolled by row (i.e., deformation()(1,2) comes before deformation(2,1))
   *      * occupation()
   *      * displacement(), unrolled by column (i.e., displacement()(2,1) comes before displacement()(1,2) )
   *   And then finding the symmetrically equivalent config_vec whose elements have the highest attainable lexicographic order.
   */
  //*******************************************************************************
  /// Check if Configuration is canonical w.r.t. supercell factor group specified by it_begin && it_end
  /// tolerance '_tol' is used for fuzzy comparisons

  bool ConfigDoF::is_canonical(PermuteIterator it_begin, PermuteIterator it_end, double tol) const {
    return _is_canonical(it_begin, it_end, nullptr, tol);
  }

  //*******************************************************************************
  /// Check if Configuration is canonical w.r.t. supercell factor group specified by it_begin && it_end
  /// tolerance '_tol' is used for fuzzy comparisons
  /// populates factorgroup at the same time (only when is_canonical evaluates to 'true'), since algorithms are so similar

  bool ConfigDoF::is_canonical(PermuteIterator it_begin, PermuteIterator it_end, std::vector<PermuteIterator> &factorgroup, double tol) const {
    return _is_canonical(it_begin, it_end, &factorgroup, tol);
  }

  //*******************************************************************************
  /**
   *   Returns an equivalent Configuraiton in canonical form
   *      (largest-valued equivalent bitstring).
   */
  //*******************************************************************************

  ConfigDoF ConfigDoF::canonical_form(PermuteIterator it_begin, PermuteIterator it_end, double _tol) const {
    PermuteIterator it_canon;
    return _canonical_form(it_begin, it_end, it_canon, nullptr, _tol);
  }

  //*******************************************************************************
  /**
   *   Returns an equivalent Configuraiton in canonical form
   *      (largest-valued equivalent bitstring).
   *      Permutation that resulted in canonical form is stored in 'it_canon'.
   */
  //*******************************************************************************

  ConfigDoF ConfigDoF::canonical_form(PermuteIterator it_begin, PermuteIterator it_end,
                                      PermuteIterator &it_canon, double tol) const {
    return _canonical_form(it_begin, it_end, it_canon, nullptr, tol);
  }

  //*******************************************************************************
  /**
   *   Returns an equivalent Configuraiton in canonical form
   *      (largest-valued equivalent bitstring).
   *      Permutation that resulted in canonical form is stored in 'it_canon'.
   *      The 'factorgroup' of the ConfigDof (Array of permutations that leave the canonical form unchanged)
   *           is populated during the calculation
   */
  //*******************************************************************************

  ConfigDoF ConfigDoF::canonical_form(PermuteIterator it_begin, PermuteIterator it_end,
                                      PermuteIterator &it_canon, std::vector<PermuteIterator> &factorgroup, double tol) const {
    return _canonical_form(it_begin, it_end, it_canon, &factorgroup, tol);
  }

  //*******************************************************************************
  /// This version calculates the factor group of the configuration, but only if it is canonical (i.e., returns true), since loop terminates
  /// early otherwise.  This private method uses the pointer fg_ptr so that we only need one implementation for the various different public methods above
  bool ConfigDoF::_is_canonical(PermuteIterator it_begin, PermuteIterator it_end, std::vector<PermuteIterator> *fg_ptr, double _tol) const {

    if(fg_ptr)
      fg_ptr->clear();

    Index i, j;
    Index permute_ind;
    PermuteIterator it_next_fg;
    Eigen::Matrix3d fg_cart_op;
    Eigen::Matrix3d new_F;

    displacement_matrix_t new_disp;

    while(it_begin != it_end) {
      bool skip_to_next_op = false;
      fg_cart_op = it_begin.sym_op().matrix();
      it_next_fg = it_begin.begin_next_fg_op();

      // check for canonical strain first, since it is fastest
      if(is_strained()) {
        new_F = fg_cart_op * deformation() * fg_cart_op.transpose();
        for(i = 0; i < 3; i++) {
          for(j = 0; j < 3; j++) {
            if(new_F(i, j) > deformation()(i, j) + _tol) {
              if(fg_ptr)
                fg_ptr->clear();
              return false;
            }
            else if(new_F(i, j) < deformation()(i, j) - _tol) {
              skip_to_next_op = true;
              break;
            }
          }
          if(skip_to_next_op)
            break;
        }
      }

      if(skip_to_next_op) {
        it_begin = it_next_fg;
        continue;
      }

      if(displacement().cols())
        new_disp = fg_cart_op * displacement();

      // check for canonical occupation
      while(it_begin != it_next_fg) {
        skip_to_next_op = false;
        if(occupation().size()) {
          it_begin = it_begin;
          for(i = 0; i < size(); i++) {
            permute_ind = it_begin.permute_ind(i);
            if(occ(permute_ind) > occ(i)) {
              if(fg_ptr)
                fg_ptr->clear();
              return false;
            }
            else if(occ(permute_ind) < occ(i)) {
              skip_to_next_op = true;
              break;
            }
          }
        }

        if(skip_to_next_op) {
          ++it_begin;
          continue;
        }

        // check for canonical displacement
        if(displacement().cols()) {
          for(i = 0; i < size(); i++) {
            for(j = 0; j < 3; j++) {
              permute_ind = it_begin.permute_ind(i);
              if(disp(permute_ind)[j] > disp(i)[j] + _tol) {
                if(fg_ptr)
                  fg_ptr->clear();
                return false;
              }
              else if(disp(permute_ind)[j] < disp(i)[j] - _tol) {
                skip_to_next_op = true; // for sake of consistency & calculating factorgroup
                break;
              }
            }
            if(skip_to_next_op)
              break;
          }
        }
        if(!skip_to_next_op && fg_ptr)
          fg_ptr->push_back(it_begin);
        ++it_begin;
      }//\End Translation loop

    }//\End factor_group loop

    return true;

  };

  //*******************************************************************************
  /**
   *   Returns an equivalent Configuraiton in canonical form
   *      (largest-valued equivalent bitstring).
   *      Does not treat displacements, strain, etc.
   *      Permutation that resulted in canonical form is stored in 'it_canon'.
   *
   *   Might want to rewrite using a provided 'canon_config' for memory/speed issues.
   */
  //*******************************************************************************

  ConfigDoF ConfigDoF::_canonical_form(PermuteIterator it_begin, PermuteIterator it_end,
                                       PermuteIterator &it_canon, std::vector<PermuteIterator> *fg_ptr, double _tol) const {
    // canonical form is 'largest-valued' configuration bitstring
    if(fg_ptr)
      fg_ptr->clear();

    Index i, j;
    Index permute_ind;
    it_canon = it_begin;
    PermuteIterator it_next_fg;
    Eigen::Matrix3d fg_cart_op;
    Eigen::Matrix3d new_F;

    ConfigDoF best_config(*this);

    displacement_matrix_t new_disp;

    while(it_begin != it_end) {
      fg_cart_op = it_begin.sym_op().matrix();
      it_next_fg = it_begin.begin_next_fg_op();

      bool skip_to_next_op = false;

      // check for canonical strain first, since it is fastest
      if(is_strained()) {
        new_F = fg_cart_op * deformation() * fg_cart_op.transpose();
        for(i = 0; i < 3; i++) {
          for(j = 0; j < 3; j++) {
            if(new_F(i, j) > best_config.deformation()(i, j) + _tol) {
              it_canon = it_begin;
              best_config = it_canon * (*this);
              //if we were filling up the factor group previously, we need to clear it and start over:
              if(fg_ptr) {
                fg_ptr->clear();
                fg_ptr->push_back(it_canon);
              }

              // We will still do translation loop, but will do so using the next translation
              // this is safe even if there is only one translation operationl
              //    (++it_begin == it_next_fg, so we won't even enter next while loop)
              ++it_begin;
              break;
            }
            else if(new_F(i, j) < best_config.deformation()(i, j) - _tol) {
              skip_to_next_op = true;
              break;
            }
          }
          if(j < 3)
            break;
        }
      }//end if(is_strained())

      if(skip_to_next_op) {
        it_begin = it_next_fg;
        continue;
      }


      // check for canonical displacement
      if(has_displacement()) {
        new_disp = fg_cart_op * displacement();
      }

      while(it_begin != it_next_fg) {
        skip_to_next_op = false;

        if(has_occupation()) {      // check for canonical occupation
          for(i = 0; i < size(); i++) {
            permute_ind = it_begin.permute_ind(i);
            if(best_config.occ(i) < occ(permute_ind)) {
              skip_to_next_op = true;
              it_canon = it_begin;
              best_config = it_canon * (*this);
              //if we were filling up the factor group previously, we need to clear it and start over:
              if(fg_ptr) {
                fg_ptr->clear();
                fg_ptr->push_back(it_canon);
              }

              break;
            }
            else if(occ(permute_ind) < best_config.occ(i)) {
              skip_to_next_op = true;
              break;
            }
          }
        } //\End occupation check

        if(skip_to_next_op) {
          ++it_begin;
          continue;
        }

        // check for canonical displacement
        if(has_displacement()) {
          //new_disp = fg_cart_op * displacement();

          for(i = 0; i < size(); i++) {
            for(j = 0; j < 3; j++) {
              permute_ind = it_begin.permute_ind(i);
              if(best_config.disp(i)[j] < new_disp(j, permute_ind) - _tol) {
                skip_to_next_op = true; // not necessary, but in case we add other DoFs...
                it_canon = it_begin;
                best_config = it_canon * (*this);
                //if we were filling up the factor group previously, we need to clear it and start over:
                if(fg_ptr) {
                  fg_ptr->clear();
                  fg_ptr->push_back(it_canon);
                }

                break;
              }
              else if(new_disp(j, permute_ind) < best_config.disp(i)[j] - _tol) {
                skip_to_next_op = true;
                break;
              }
            }
            if(skip_to_next_op)
              break;
          }
        }

        // We get the factor_group here with minimal overhead:
        if(!skip_to_next_op && fg_ptr) { // if this is the case, then it_begin has the same effect as it_canon
          fg_ptr->push_back(it_begin);
        }

        ++it_begin;


      }
      // at end of loop, we have it_begin = it_next_fg;
    }

    if(fg_ptr) {
      //recycle it_begin, since we're done:
      it_begin = it_canon.inverse();
      for(Index i = 0; i < fg_ptr->size(); i++) {
        fg_ptr->at(i) = it_begin * (fg_ptr->at(i));
      }
    }
    return best_config;

  };

  //*******************************************************************************

  jsonParser &ConfigDoF::to_json(jsonParser &json) const {
    json = jsonParser::object();
    if(occupation().size())
      json["occupation"] = occupation();
    if(displacement().size())
      json["displacement"] = displacement();
    if(is_strained())
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
      m_is_strained = true;
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

  ConfigDoF operator*(const PermuteIterator &it, const ConfigDoF &dof) {
    ConfigDoF tconfig(dof.size());
    Eigen::Matrix3d fg_cart_op = it.sym_op().matrix();
    if(dof.is_strained())
      tconfig.set_deformation(fg_cart_op * dof.deformation() * fg_cart_op.transpose());

    Permutation tperm(*it);
    if(dof.occupation().size())
      tconfig.set_occupation(tperm * dof.occupation());

    if(dof.displacement().cols()) {
      Eigen::MatrixXd new_disp = fg_cart_op * dof.displacement();
      tconfig.set_displacement(Eigen::MatrixXd(3, dof.size()));
      for(Index i = 0; i < dof.size(); i++)
        tconfig.disp(i) = new_disp.col(tperm[i]);
    }

    return tconfig;
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

