#include "casm/clex/CompositionConverter.hh"

#include "casm/crystallography/Structure.hh"

// currently still relies on this for getting standard composition axes
#include "casm/clex/ParamComposition.hh"
#include "casm/casm_io/json_io/container.hh"

namespace CASM {

  /// \brief The dimensionality of the composition space
  ///
  /// Examples:
  /// - ZrOa: 1
  /// - AaB(1-a): 1
  /// - AaBbC(1-a-b): 2
  CompositionConverter::size_type CompositionConverter::independent_compositions() const {
    return m_to_x.rows();
  }

  /// \brief Composition variable names: "a", "b", ...
  std::string CompositionConverter::comp_var(size_type i) {
    return std::string(1, (char)(i + (int) 'a'));
  }


  /// \brief The order of components in mol composition vectors
  std::vector<std::string> CompositionConverter::components() const {
    return m_components;
  }

  /// \brief The mol composition of the parameteric composition axes origin
  ///
  /// - Matches order from components()
  ///
  Eigen::VectorXd CompositionConverter::origin() const {
    return m_origin;
  }

  /// \brief The mol composition of the parameteric composition axes end members
  ///
  /// - Matches order from components()
  ///
  Eigen::VectorXd CompositionConverter::end_member(size_type i) const {
    return m_end_members.col(i);
  }
  
  /// \brief Return the matrix Mij = dx_i/dn_j
  Eigen::MatrixXd CompositionConverter::dparam_dmol() const {
    return m_to_x;
  }
  
  /// \brief Return the matrix Mij = dn_i/dx_j
  Eigen::MatrixXd CompositionConverter::dmol_dparam() const {
    return m_to_n;
  }

  /// \brief Convert number of atoms per prim, 'n' to parametric composition 'x'
  ///
  /// \param n mol composition, matches order from components()
  ///
  Eigen::VectorXd CompositionConverter::param_composition(const Eigen::VectorXd &n) const {
    return m_to_x * (n - m_origin);
  }

  /// \brief Convert change in number of atoms per prim, 'dn' to change in parametric composition 'dx'
  ///
  /// \param dn mol composition, matches order from components()
  ///
  Eigen::VectorXd CompositionConverter::dparam_composition(const Eigen::VectorXd &dn) const {
    return m_to_x * dn;
  }

  /// \brief Convert parametric composition, 'x', to number of atoms per prim, 'n'
  ///
  /// - Matches order from components()
  ///
  Eigen::VectorXd CompositionConverter::mol_composition(const Eigen::VectorXd &x) const {
    return m_origin + m_to_n * x;
  }

  /// \brief Convert change in parametric composition, 'dx', to change in number of atoms per prim, 'dn'
  ///
  /// - Matches order from components()
  ///
  Eigen::VectorXd CompositionConverter::dmol_composition(const Eigen::VectorXd &dx) const {
    return m_to_n * dx;
  }

  /// \brief Convert dG/dn to dG/dx
  ///
  /// Notes:
  ///
  /// For conversion from dG/dx to dG/dn:
  /// - dG/dni = sum_j dG/dxj * dxj/dni
  /// - mu_n = m_to_mu_n * mu_x, where m_to_mu_n(i,j) = dxj/dni
  /// - Also, x = m_to_x*(n - origin)
  /// - so, dxj/dni = m_to_x(j,i)
  /// - that means m_to_mu_n == m_to_x.transpose()
  ///
  /// Similarly, for conversion from dG/dn to dG/dx:
  /// - mu_x = m_to_mu_x*mu_n, where m_to_mu_x(i,j) = dnj/dxi
  /// - and n = m_origin + m_to_n * x,
  /// - so dnj/dxi = m_to_n(j,i)
  /// - that means m_to_mu_x == m_to_n.transpose()
  ///
  Eigen::VectorXd CompositionConverter::param_chem_pot(const Eigen::VectorXd chem_pot) const {
    return m_to_n.transpose() * chem_pot;
  }

  /// \brief Convert dG/dx to dG/dn
  Eigen::VectorXd CompositionConverter::chem_pot(const Eigen::VectorXd param_chem_pot) const {
    return m_to_x.transpose() * param_chem_pot;
  }


  /// \brief Return formula for x->n
  std::string CompositionConverter::mol_formula() const {

    // n = origin + m_to_n * x

    int composition_var = (int)'a';
    std::stringstream tstr;

    // for each molecule:
    for(int i = 0; i < m_components.size(); i++) {

      bool first_char = true;

      // print mol name 'A('
      tstr << m_components[i] << "(";

      // constant term from origin
      // print 'x' if x != 0
      if(!almost_zero(m_origin(i))) {
        first_char = false;
        tstr << m_origin(i);
      }

      // terms from m_to_n columns
      for(int j = 0; j < m_to_n.cols(); j++) {

        // print nothing if x == 0
        if(almost_zero(m_to_n(i, j))) {
          continue;
        }

        // print '+' if x>0 mid-expression
        if(!first_char && m_to_n(i, j) > 0)
          tstr << '+';
        // print '-' if x<0
        else if(m_to_n(i, j) < 0)
          tstr << '-';
        // print absolute value of x if |x|!=1
        if(!almost_equal(std::abs(m_to_n(i, j)), 1.0))
          tstr << std::abs(m_to_n(i, j));
        //print variable ('a','b',etc...)
        tstr << (char)(composition_var + j);

        first_char = false;
      }

      // close ')'
      tstr << ")";
    }

    return tstr.str();

  }

  /// \brief Return formula for n->x
  std::string CompositionConverter::param_formula() const {

    // x_i = m_to_x*(n - origin) = m_to_x*n - m_to_x*origin

    int composition_var = (int)'a';
    std::stringstream tstr;

    Eigen::VectorXd v = -m_to_x * m_origin;

    // for each independent composition:
    for(int i = 0; i < independent_compositions(); i++) {

      bool first_char = true;

      // print mol name 'a('
      tstr << (char)(composition_var + i) << "(";

      // constant term from origin
      // print 'x' if x != 0
      if(!almost_zero(v(i))) {
        first_char = false;
        tstr << v(i);
      }

      // terms from m_to_x columns
      for(int j = 0; j < m_to_x.cols(); j++) {

        double coeff = m_to_x(i, j);

        // print nothing if n == 0
        if(almost_zero(coeff)) {
          continue;
        }

        // print 'A' or '+A' if n == 1
        if(almost_zero(coeff - 1)) {
          if(!first_char) {
            tstr << '+';
          }
          tstr << m_components[j];
        }

        // print '-A' if n == -1
        else if(almost_zero(coeff + 1)) {
          tstr << '-' << m_components[j];
        }

        // print 'nA' or '+nA' if n > 0
        else if(coeff > 0) {
          if(!first_char) {
            tstr << '+';
          }
          tstr << coeff << m_components[j];
        }

        // print '-nA' if n < 0
        else {
          tstr << coeff << m_components[j];
        }

        first_char = false;
      }

      // close ')'
      tstr << ")";
    }

    return tstr.str();

  }

  /// \brief Return formula for param_chem_pot->chem_pot
  std::string CompositionConverter::chem_pot_formula(int indent) const {
    // chem_pot = m_to_x.transpose() * param_chem_pot;
    std::stringstream ss;
    Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", std::string(' ', indent + 2) + "[", "]");
    ss << std::string(' ', indent) << "chem_pot = X * param_chem_pot, where X = \n" << m_to_x.transpose().format(CleanFmt);
    return ss.str();
  }

  /// \brief Return formula for chem_pot->param_chem_pot
  std::string CompositionConverter::param_chem_pot_formula(int indent) const {
    // param_chem_pot = m_to_n.transpose() * chem_pot;
    std::stringstream ss;
    Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", std::string(' ', indent + 2) + "[", "]");
    ss << std::string(' ', indent) << "param_chem_pot = X * chem_pot, where X = \n" << m_to_n.transpose().format(CleanFmt);
    return ss.str();
  }


  /// \brief Return formula for origin
  std::string CompositionConverter::origin_formula() const {
    return _n_formula(m_origin);
  }

  /// \brief Return formula for end member
  std::string CompositionConverter::end_member_formula(size_type i) const {
    return _n_formula(m_end_members.col(i));
  }


  /// \brief Helps make variadic template constructor possible
  void CompositionConverter::_add_end_member(Eigen::VectorXd _end_member) {

    _check_size(_end_member);

    int r = m_components.size();
    int c = m_end_members.cols();

    Eigen::MatrixXd tmp(r, c + 1);
    tmp.leftCols(c) = m_end_members;
    tmp.col(c) = _end_member;
    m_end_members = tmp;
  }

  /// \brief Check that origin and end member vectors have same size as the number of components
  void CompositionConverter::_check_size(const Eigen::VectorXd &vec) const {
    if(m_components.size() != vec.size()) {
      throw std::runtime_error(
        std::string("Error in CompositionConverter: origin or end member vector size does not match components size."));
    }
  }

  /// \brief Calculate conversion matrices m_to_n and m_to_x
  void CompositionConverter::_calc_conversion_matrices() {

    // calculate m_to_x and m_to_n:
    //
    // x = m_to_x*(n - origin)
    // n = origin + m_to_n*x

    // end_members.col(i) corresponds to x such that x[i] = 1, x[j!=i] = 0,
    //  -> end_members.col(i) = origin + m_to_n.col(i)

    // r > c
    int r = m_components.size();
    int c = m_end_members.cols();

    m_to_n = Eigen::MatrixXd(r, c);
    for(int i = 0; i < c; i++) {
      m_to_n.col(i) = m_end_members.col(i) - m_origin;
    }

    // x = m_to_x*(n - origin)
    //   -> x = m_to_x*(origin + m_to_n*x - origin)
    //   -> x = m_to_x*m_to_n*x
    Eigen::MatrixXd tmp = Eigen::MatrixXd::Identity(r, r);
    tmp.leftCols(c) = m_to_n;

    // to ensure tmp is full rank, fill right cols with vectors orthogonal to the m_to_n space
    Eigen::FullPivHouseholderQR<Eigen::MatrixXd> qr(m_to_n);
    Eigen::MatrixXd Q = qr.matrixQ();
    tmp.rightCols(r - c) = Q.rightCols(r - c);

    // get the m_to_x matrix
    m_to_x = tmp.inverse().topRows(c);
  }

  /// \brief Return formula for 'n'
  std::string CompositionConverter::_n_formula(const Eigen::VectorXd &vec) const {

    std::stringstream tstr;

    // for each molecule:
    for(int i = 0; i < vec.size(); i++) {

      // print 'A' if x == 1
      if(almost_zero(vec(i) - 1)) {
        tstr << m_components[i];
      }
      // print 'A(x)' if x != 0
      else if(!almost_zero(vec(i))) {
        tstr << m_components[i] << "(" << vec(i) << ")";
      }

    }

    return tstr.str();
  }

  /// \brief Pretty-print map of name/CompositionConverter pairs
  ///
  /// \param stream Output stream
  /// \param map Map of name/CompositionConverter pairs
  /// \param name Name for this set of composition axes
  ///
  void display_composition_axes(std::ostream &stream, const std::map<std::string, CompositionConverter> &map) {

    if(map.size() == 0) {
      return;
    }

    stream << std::setw(10) << "KEY" << " ";
    stream << std::setw(10) << "ORIGIN" << " ";
    for(int i = 0; i < map.begin()->second.independent_compositions(); i++) {
      stream << std::setw(10) << (char)((int)'a' + i) << " ";
    }
    stream << "    ";
    stream << "GENERAL FORMULA";
    stream << std::endl;

    stream << std::setw(10) << "  ---" << " ";
    stream << std::setw(10) << "  ---" << " ";
    for(int i = 0; i < map.begin()->second.independent_compositions(); i++) {
      stream << std::setw(10) << "  ---" << " ";
    }
    stream << "    ";
    stream << "---" << std::endl;

    for(auto it = map.cbegin(); it != map.cend(); ++it) {
      stream << std::setw(10) << it->first << " ";
      stream << std::setw(10) << it->second.origin_formula() << " ";
      for(int i = 0; i < it->second.independent_compositions(); ++i) {
        stream << std::setw(10) << it->second.end_member_formula(i) << " ";
      }
      stream << "    ";
      stream << std::setw(10) << it->second.mol_formula() << "\n";
    }
  }

  /// \brief Serialize CompositionConverter to JSON
  jsonParser &to_json(const CompositionConverter &f, jsonParser &json) {

    json = jsonParser::object();
    json["components"] = f.components();
    json["independent_compositions"] = f.independent_compositions();
    json["origin"] = f.origin();
    int var = (int) 'a';
    for(int i = 0; i < f.independent_compositions(); i++) {
      json[std::string(1, (char)(var + i)) ] = f.end_member(i);
    }
    json["mol_formula"] = f.mol_formula();
    json["param_formula"] = f.param_formula();

    return json;
  }

  /// \brief Deserialize CompositionConverter from JSON
  void from_json(CompositionConverter &f, const jsonParser &json) {

    try {

      std::vector<std::string> components;
      Eigen::VectorXd origin;

      int independent_compositions;

      from_json(components, json["components"]);
      from_json(origin, json["origin"]);

      from_json(independent_compositions, json["independent_compositions"]);
      Eigen::MatrixXd end_members(components.size(), independent_compositions);
      int var = (int) 'a';
      Eigen::VectorXd tvec;
      for(int i = 0; i < independent_compositions; i++) {
        from_json(tvec, json[std::string(1, (char)(var + i))]);
        end_members.col(i) = tvec;
      }

      f = CompositionConverter(components.begin(), components.end(), origin, end_members);
      return;
    }
    catch(...) {
      /// re-throw exceptions
      throw;
    }
  }

  /// \brief Generate a column matrix containing all the possible molecular end members
  ///
  /// \param prim A Structure to find the end members of (does
  ///             not check if it is actually primitive).
  ///
  /// - Each column corresponds to a point in composition space (specifying number
  ///   of each Molecule per prim)
  /// - Each row corresponds to a Molecule, ordered as from Structure::get_struc_molecule,
  ///   with units number Molecule / prim
  Eigen::MatrixXd end_members(const Structure &prim) {
    ParamComposition param_comp(prim);
    param_comp.generate_components();
    param_comp.generate_sublattice_map();
    param_comp.generate_prim_end_members();
    return param_comp.get_prim_end_members().transpose();
  }

  /// \brief Return the composition space of a Structure
  ///
  /// \param prim A Structure to find the standard composition space for (does
  ///             not check if it is actually primitive).
  /// \param tol tolerance for checking rank (default 1e-14)
  ///
  /// - Each column corresponds to an orthogonal vector in composition space
  /// - Each row corresponds to a Molecule, ordered as from Structure::get_struc_molecule,
  ///   with units number Molecule / prim
  Eigen::MatrixXd composition_space(const Structure &prim, double tol) {
    auto Qr = end_members(prim).fullPivHouseholderQr();
    Qr.setThreshold(tol);
    auto Q = Qr.matrixQ();
    return Q.leftCols(Qr.rank());
  }

  /// \brief Return the null composition space of a Structure
  ///
  /// \param prim A Structure to find the standard composition space for (does
  ///             not check if it is actually primitive).
  /// \param tol tolerance for checking rank (default 1e-14)
  ///
  /// - Each column corresponds to an orthogonal vector in composition space
  /// - Each row corresponds to a Molecule, ordered as from Structure::get_struc_molecule,
  ///   with units number Molecule / prim
  Eigen::MatrixXd null_composition_space(const Structure &prim, double tol) {
    auto Qr = end_members(prim).fullPivHouseholderQr();
    Qr.setThreshold(tol);
    auto Q = Qr.matrixQ();
    return Q.rightCols(Q.cols() - Qr.rank());
  }

}
