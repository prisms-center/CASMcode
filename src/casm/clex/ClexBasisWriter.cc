namespace CASM {
  namespace ClexBasisWriter_impl {

    std::string ClexBasis::clexulator_member_definitions(std::string const &class_name,
                                                         ClexBasis const &clex,
                                                         std::vector<std::unique_ptr<OrbitFunctionWriter> > const &orbit_func_writers;
                                                         std::string const &indent) {
      std::stringstream ss;
      ss <<

         indent << "/// \\brief Clone the Clexulator\n" <<
         indent << "virtual " << class_name << "* _clone() const override {\n" <<
         indent << "  return new " << class_name << "(*this);\n" <<
         indent << "}\n\n";

      for(auto const &writer_ptr : orbit_func_writers) {
        writer_ptr->print_typedefs(ss, class_name, indent);
      }

      for(auto const &writer_ptr : orbit_func_writers) {
        writer_ptr->print_eval_table_definitions(ss, class_name, clex, indent);
      }


      indent << "// typedef for method pointers\n" <<
             indent << "typedef double (" << class_name << "::*BasisFuncPtr)() const;\n\n" <<

             indent << "// typedef for method pointers\n" <<
             indent << "typedef double (" << class_name << "::*DeltaBasisFuncPtr)(int, int) const;\n\n" <<

             indent << "// array of pointers to member functions for calculating basis functions\n" <<
             indent << "BasisFuncPtr m_orbit_func_table[" << N_corr << "];\n\n" <<

             indent << "// array of pointers to member functions for calculating flower functions\n" <<
             indent << "BasisFuncPtr m_flower_func_table[" << Nsublat << "][" << N_corr << "];\n\n" <<

             indent << "// array of pointers to member functions for calculating DELTA flower functions\n" <<
             indent << "DeltaBasisFuncPtr m_delta_func_table[" << Nsublat << "][" << N_corr << "];\n\n";

      auto it(site_bases().begin()), end_it(site_bases().end());
      for(; it != end_it; ++it) {
        DoFType::traits(it->first).print_clexulator_member_definitions(ss, it->second(), indent);
      }
      return ss.str();
    }

    //*******************************************************************************************

    std::string ClexBasis::clexulator_private_method_definitions(std::string const &class_name,
                                                                 Index N_corr,
                                                                 Index N_sublat,
                                                                 std::string const &indent) {

      std::stringstream ss;
      ss <<
         indent << "  /// \\brief Clone the " << class_name << "\n" <<
         indent << "  std::unique_ptr<" << class_name << "> clone() const { \n" <<
         indent << "    return std::unique_ptr<" << class_name << ">(_clone()); \n" <<
         indent << "  }\n\n" <<

         indent << "  /// \\brief Calculate contribution to global correlations from one unit cell\n" <<
         indent << "  void _calc_global_corr_contribution(double *corr_begin) const override;\n\n" <<

         indent << "  /// \\brief Calculate contribution to select global correlations from one unit cell\n" <<
         indent << "  void _calc_restricted_global_corr_contribution(double *corr_begin, size_type const* ind_list_begin, size_type const* ind_list_end) const override;\n\n" <<

         indent << "  /// \\brief Calculate point correlations about basis site 'b_index'\n" <<
         indent << "  void _calc_point_corr(int b_index, double *corr_begin) const override;\n\n" <<

         indent << "  /// \\brief Calculate select point correlations about basis site 'b_index'\n" <<
         indent << "  void _calc_restricted_point_corr(int b_index, double *corr_begin, size_type const* ind_list_begin, size_type const* ind_list_end) const override;\n\n" <<

         indent << "  /// \\brief Calculate the change in point correlations due to changing an occupant\n" <<
         indent << "  void _calc_delta_point_corr(int b_index, int occ_i, int occ_f, double *corr_begin) const override;\n\n" <<

         indent << "  /// \\brief Calculate the change in select point correlations due to changing an occupant\n" <<
         indent << "  void _calc_restricted_delta_point_corr(int b_index, int occ_i, int occ_f, double *corr_begin, size_type const* ind_list_begin, size_type const* ind_list_end) const override;\n\n";


      auto it(site_bases().begin()), end_it(site_bases().end());
      for(; it != end_it; ++it) {
        DoFType::traits(it->first).print_clexulator_private_method_definitions(ss, it->second(), indent);
      }

      ss <<
         indent << "  //default functions for basis function evaluation \n" <<
         indent << "  double zero_func() const{ return 0.0;};\n" <<
         indent << "  double zero_func(int,int) const{ return 0.0;};\n\n";

      return ss.str();
    }

    //*******************************************************************************************

    std::string ClexBasis::clexulator_public_method_definitions(std::string const &class_name,
                                                                Index N_corr,
                                                                Index N_sublat,
                                                                std::string const &indent) {
      std::stringstream ss;
      ss <<
         indent << "  " << class_name << "();\n\n" <<
         indent << "  ~" << class_name << "();\n\n";

      auto it(site_bases().begin()), end_it(site_bases().end());
      for(; it != end_it; ++it) {
        DoFType::traits(it->first).print_clexulator_public_method_definitions(ss, it->second(), indent);
      }
      return ss.str();
    }

    //*******************************************************************************************

    std::tuple<std::string, std::string> clexulator_orbit_function_strings(ClexBasis::BSetOrbit const &_bset_orbit,
                                                                           OrbitType const &_clust_orbit,
                                                                           std::function<std::string(Index, Index)> method_namer,
                                                                           PrimNeighborList &_nlist,
                                                                           std::vector<std::unique_ptr<FunctionVisitor> > const &labelers) {

      std::stringstream bfunc_def_stream;
      std::stringstream bfunc_imp_stream;
      bool make_newline = false;

      std::vector<std::string>formulae = orbit_function_cpp_strings(_bset_orbit, _clust_orbit, _nlist, labelers);
      std::string method_name;
      make_newline = false;
      for(Index nf = 0; nf < formulae.size(); nf++) {
        if(!formulae[nf].size())
          continue;
        make_newline = true;

        method_name = method_namer(0, nf);
        bfunc_def_stream <<
                         indent << "  double " << method_name << "() const;\n";

        bfunc_imp_stream <<
                         indent << "double " << class_name << "::" << method_name << "() const{\n" <<
                         indent << "  return " << formulae[nf] << ";\n" <<
                         indent << "}\n";
      }
      if(make_newline) {
        bfunc_imp_stream << '\n';
        bfunct_def_stream << '\n';
      }
      return std::tuple<std::string, std::string>(bfunc_def_stream.str(), bfunc_imp_stream.str());
    }

    //*******************************************************************************************

    std::tuple<std::string, std::string> clexulator_flower_function_strings(ClexBasis::BSetOrbit const &_bset_orbit,
                                                                            OrbitType const &_clust_orbit,
                                                                            std::function<std::string(Index, Index)> method_namer,
                                                                            PrimNeighborList &_nlist,
                                                                            std::vector<std::unique_ptr<FunctionVisitor> > const &labelers) {
      std::stringstream bfunc_def_stream, bfunc_imp_stream;
      std::string method_name;
      // loop over flowers (i.e., basis sites of prim)
      auto it(_nlist.sublat_indices().begin()), end_it(_nlist.sublat_indices().end());
      for(; it != end_it; ++it) {
        Index nb = *it;
        auto nlist_index = find_index(_nlist.sublat_indices(), nb);

        std::map< UnitCell, std::vector< std::string > > formulae_pack = flower_function_cpp_strings(_bset_orbit,
                                                                         UnaryIdentity<BasisSet>(),
                                                                         _clust_orbit,
                                                                         _nlist,
                                                                         labelers,
                                                                         sublat_index);

        if(formulae_pack.size() > 1) {
          throw std::runtime_error("Printing clexulators is only supported for fully periodic systems!\n");
        }
        if(formulae_pack.size() == 0)
          continue;

        std::vector<std::string> const &formulae(formulae_pack.begin()->second);

        for(Index nf = 0; nf < formulae.size(); nf++) {
          if(!formulae[nf].size())
            continue;
          make_newline = true;

          method_name = method_namer(nb, nf);
          bfunc_def_stream <<
                           indent << "  double " << method_name << "() const;\n";

          bfunc_imp_stream <<
                           indent << "double " << class_name << "::" << method_name << "() const{\n" <<
                           indent << "  return " << formulae[nf] << ";\n" <<
                           indent << "}\n";

        }
      }
      if(make_newline) {
        bfunc_imp_stream << '\n';
        bfunc_def_stream << '\n';
      }

      return std::tuple<std::string, std::string>(bfunc_def_stream.str(), bfunc_imp_stream.str());
    }


    //*******************************************************************************************

    std::tuple<std::string, std::string, std::string> clexulator_dflower_function_strings(ClexBasis::BSetOrbit const &_bset_orbit,
        OrbitType const &_clust_orbit,
        std::function<std::string(Index, Index)> method_namer,
        PrimNeighborList &_nlist,
        std::vector<std::unique_ptr<FunctionVisitor> > const &labelers,
        FunctionVisitor const &prefactor_labeler) {
      std::stringstream bfunc_def_stream, bfunc_imp_stream;

      // loop over flowers (i.e., basis sites of prim)
      auto it(_nlist.sublat_indices().begin()), end_it(_nlist.sublat_indices().end());
      for(; it != end_it; ++it) {
        Index nb = *it;
        auto nlist_index = find_index(_nlist.sublat_indices(), nb);
        std::vector<std::string> formulae(_bset_orbit[0].size());

        // loop over site basis functions
        BasisSet site_basis(site_bases()[dof_key][nb]);
        site_basis.accept(site_func_labeler);
        for(Index nsbf = 0; nsbf < site_basis.size(); nsbf++) {
          Function const *func_ptr = site_basis[nsbf];

          auto get_quotient_bset = [func_ptr](BasisSet const & bset)->BasisSet{
            return bset.poly_quotient_set(func_ptr);
          };

          std::map< UnitCell, std::vector< std::string > > formulae_pack = flower_function_cpp_strings(_bset_orbit,
                                                                           get_quotient_bset,
                                                                           _clust_orbit,
                                                                           _nlist,
                                                                           labelers,
                                                                           sublat_index);

          if(formulae_pack.size() > 1) {
            throw std::runtime_error("Printing clexulators is only supported for fully periodic systems!\n");
          }
          if(formulae_pack.size() == 0)
            continue;

          std::vector<std::string> const &tformulae(formulae_pack.begin()->second);

          for(Index nf = 0; nf < tformulae.size(); nf++) {
            if(!tformulae[nf].size())
              continue;

            if(formulae[nf].size())
              formulae[nf] += " + ";

            formulae[nf] += site_basis[nsbf]->formula();

            if(tformulae[nf] == "1" || tformulae[nf] == "(1)")
              continue;

            formulae[nf] += "*";
            formulae[nf] += tformulae[nf];
            //formulae[nf] += ")";
          }
        }

        std::string method_name;
        for(Index nf = 0; nf < formulae.size(); nf++) {
          if(!formulae[nf].size()) {
            continue;
          }
          make_newline = true;

          method_name = method_namer(nb, nf);

          bfunc_def_stream <<
                           indent << "  double " << bfunc_name << "(int occ_i, int occ_f) const;\n";

          bfunc_imp_stream <<
                           indent << "double " << class_name << "::" << bfunc_name << "(int occ_i, int occ_f) const{\n" <<
                           indent << "  return " << formulae[nf] << ";\n" <<
                           indent << "}\n";
        }
        if(make_newline) {
          bfunc_imp_stream << '\n';
          bfunc_def_stream << '\n';
        }
        make_newline = false;
        // \End Configuration specific part
      }//\End loop over flowers
      return std::tuple<std::string, std::string, std::string>(bfunc_name_stream.str(), bfunc_def_stream.str(), bfunc_imp_stream.str());
    }

    //*******************************************************************************************
    std::string clexulator_interface_implementation(std::string const &class_name, PrimNeighborList &_nlist, std::string const &indent) {
      std::stringstream ss;
      // Write constructor
      ss << clexulator_constructor_implementation(class_name, _nlist, indent);

      // Write destructor
      ss <<
         indent << class_name << "::~" << class_name << "(){\n" <<

         indent << "  //nothing here for now\n" <<

         indent << "}\n\n";

      // Write evaluation methods

      ss <<
         indent << "/// \\brief Calculate contribution to global correlations from one unit cell\n" <<
         indent << "void " << class_name << "::_calc_global_corr_contribution(double *corr_begin) const {\n" <<
         indent << "  for(size_type i=0; i<corr_size(); i++){\n" <<
         indent << "    *(corr_begin+i) = (this->*m_orbit_func_table[i])();\n" <<
         indent << "  }\n" <<
         indent << "}\n\n" <<

         indent << "/// \\brief Calculate contribution to select global correlations from one unit cell\n" <<
         indent << "void " << class_name << "::_calc_restricted_global_corr_contribution(double *corr_begin, size_type const* ind_list_begin, size_type const* ind_list_end) const {\n" <<
         indent << "  for(; ind_list_begin<ind_list_end; ind_list_begin++){\n" <<
         indent << "    *(corr_begin+*ind_list_begin) = (this->*m_orbit_func_table[*ind_list_begin])();\n" <<
         indent << "  }\n" <<
         indent << "}\n\n" <<

         indent << "/// \\brief Calculate point correlations about basis site 'b_index'\n" <<
         indent << "void " << class_name << "::_calc_point_corr(int b_index, double *corr_begin) const {\n" <<
         indent << "  for(size_type i=0; i<corr_size(); i++){\n" <<
         indent << "    *(corr_begin+i) = (this->*m_flower_func_table[b_index][i])();\n" <<
         indent << "  }\n" <<
         indent << "}\n\n" <<

         indent << "/// \\brief Calculate select point correlations about basis site 'b_index'\n" <<
         indent << "void " << class_name << "::_calc_restricted_point_corr(int b_index, double *corr_begin, size_type const* ind_list_begin, size_type const* ind_list_end) const {\n" <<
         indent << "  for(; ind_list_begin<ind_list_end; ind_list_begin++){\n" <<
         indent << "    *(corr_begin+*ind_list_begin) = (this->*m_flower_func_table[b_index][*ind_list_begin])();\n" <<
         indent << "  }\n" <<
         indent << "}\n\n" <<

         indent << "/// \\brief Calculate the change in point correlations due to changing an occupant\n" <<
         indent << "void " << class_name << "::_calc_delta_point_corr(int b_index, int occ_i, int occ_f, double *corr_begin) const {\n" <<
         indent << "  for(size_type i=0; i<corr_size(); i++){\n" <<
         indent << "    *(corr_begin+i) = (this->*m_delta_func_table[b_index][i])(occ_i, occ_f);\n" <<
         indent << "  }\n" <<
         indent << "}\n\n" <<

         indent << "/// \\brief Calculate the change in select point correlations due to changing an occupant\n" <<
         indent << "void " << class_name << "::_calc_restricted_delta_point_corr(int b_index, int occ_i, int occ_f, double *corr_begin, size_type const* ind_list_begin, size_type const* ind_list_end) const {\n" <<
         indent << "  for(; ind_list_begin<ind_list_end; ind_list_begin++){\n" <<
         indent << "    *(corr_begin+*ind_list_begin) = (this->*m_delta_func_table[b_index][*ind_list_begin])(occ_i, occ_f);\n" <<
         indent << "  }\n" <<
         indent << "}\n\n";
    }

    //*******************************************************************************************
    std::string clexulator_constructor_implementation(std::string const &class_name,
                                                      PrimNeighborList &_nlist,
                                                      std::vector<std::string> const &orbit_method_names,
                                                      std::vector< std::vector<std::string> > const &flower_method_names,
                                                      std::vector< std::vector<std::string> > const &dflower_method_names,
                                                      std::string const &indent) {
      std::strinstream ss;
      // Write constructor
      ss <<
         indent << class_name << "::" << class_name << "() :\n" <<
         indent << "  Clexulator_impl::Base(" << _nlist.size() << ", " << N_corr << ") {\n";

      auto it(site_bases().begin()), end_it(site_bases().end());
      for(; it != end_it; ++it) {
        DoFType::traits(it->first).print_to_clexulator_constructor(ss, it->second(), indent + "  ");
      }

      for(Index nf = 0; nf < orbit_method_names.size(); nf++) {
        ss <<
           indent << "  m_orbit_func_table[" << nf << "] = &" << class_name << "::" << orbit_method_names[nf] << ";\n";
      }
      ss << "\n\n";

      for(Index nb = 0; nb < flower_method_names.size(); nb++) {
        for(Index nf = 0; nf < flower_method_names[nb].size(); nf++) {
          ss <<
             indent << "  m_flower_func_table[" << nb << "][" << nf << "] = &" << class_name << "::" << flower_method_names[nb][nf] << ";\n";
        }
        ss << "\n\n";
      }

      for(Index nb = 0; nb < dflower_method_names.size(); nb++) {
        for(Index nf = 0; nf < dflower_method_names[nb].size(); nf++) {
          ss <<
             indent << "  m_delta_func_table[" << nb << "][" << nf << "] = &" << class_name << "::" << dflower_method_names[nb][nf] << ";\n";
        }
        ss << "\n\n";
      }



      // Write weight matrix used for the neighbor list
      PrimNeighborList::Matrix3Type W = _nlist.weight_matrix();
      ss << indent << "  m_weight_matrix.row(0) << " << W(0, 0) << ", " << W(0, 1) << ", " << W(0, 2) << ";\n";
      ss << indent << "  m_weight_matrix.row(1) << " << W(1, 0) << ", " << W(1, 1) << ", " << W(1, 2) << ";\n";
      ss << indent << "  m_weight_matrix.row(2) << " << W(2, 0) << ", " << W(2, 1) << ", " << W(2, 2) << ";\n\n";

      // Write neighborhood of UnitCellCoord
      // expand the _nlist to contain 'global_orbitree' (all that is needed for now)
      std::set<UnitCellCoord> nbors;
      neighborhood(std::inserter(nbors, nbors.begin()), _tree TOL);


      ss << indent << "  m_neighborhood = std::set<UnitCellCoord> {\n";
      auto it = nbors.begin();
      while(it != nbors.end()) {
        ss << indent << "    {UnitCellCoord("
           << it->sublat() << ", "
           << it->unitcell(0) << ", "
           << it->unitcell(1) << ", "
           << it->unitcell(2) << ")}";
        ++it;
        if(it != nbors.end()) {
          ss << ",";
        }
        ss << "\n";
      }
      ss << indent << "  };\n";
      ss << "\n\n";

      ss << indent <<  "  m_orbit_neighborhood.resize(corr_size());\n";
      Index lno = 0;
      for(Index no = 0; no < _tree.size(); ++no) {
        std::set<UnitCellCoord> orbit_nbors;
        orbit_neighborhood(std::inserter(orbit_nbors, orbit_nbors.begin()), _tree, nb, no, TOL);

        Index proto_index = lno;

        ss << indent << "  m_orbit_neighborhood[" << lno << "] = std::set<UnitCellCoord> {\n";
        auto it = orbit_nbors.begin();
        while(it != orbit_nbors.end()) {
          ss << indent << "    {UnitCellCoord("
             << it->sublat() << ", "
             << it->unitcell(0) << ", "
             << it->unitcell(1) << ", "
             << it->unitcell(2) << ")}";
          ++it;
          if(it != orbit_nbors.end()) {
            ss << ",";
          }
          ss << "\n";
        }
        ss << indent << "  };\n";
        ++lno;
        for(Index nf = 1; nf < _tree[no].clust_basis.size(); ++nf) {
          ss << indent << "  m_orbit_neighborhood[" << lno << "] = m_orbit_neighborhood[" << proto_index << "];\n";
          ++lno;
        }
        ss << "\n";

      }


      ss <<
         indent << "}\n\n";

      return ss.str();
    }

    //*******************************************************************************************
    /// \brief Print clexulator
    void print_clexulator(const Structure &prim,
                          std::vector<OrbitType > const &_tree
                          const PrimNeighborList &_nlist,
                          std::string class_name,
                          std::ostream &stream,
                          double xtal_tol) {

      enum FuncStringType {func_definition = 0, func_implementation = 1};

      Index Nsublat = prim.basis.size();

      Index N_corr(basis_set_size());
      std::stringstream private_def_stream, public_def_stream, interface_imp_stream, bfunc_imp_stream, bfunc_def_stream;

      std::string uclass_name;
      for(Index i = 0; i < class_name.size(); i++)
        uclass_name.push_back(std::toupper(class_name[i]));

      std::string indent(2, ' ');

      std::string private_definitions = clexulator_member_definitions(class_name, N_corr, Nsublat, indent + "  ");

      private_definitions += clexulator_private_method_definitions(class_name, N_corr, Nsublat, indent + "  ");

      std::string public_definitions = clexulator_public_method_definitions(class_name, N_corr, N_sublat, indent + "  ");


      //linear function index
      Index lf = 0;

      std::vector<std::unique_ptr(FunctionVisitor)> labelers(get_function_label_visitors());
      //std::cout << "Initialized " << labelers.size() << " labelers \n";

      std::vector<std::string> orbit_method_names(N_corr, "zero_func");

      std::vector<std::vector<std::string> > flower_method_names(Nsublat, std::vector<std::string>(N_corr, "zero_func"));

      //this is very configuration-centric
      std::vector<std::vector<std::string> > dflower_method_names(Nsublat, std::vector<std::string>(N_corr, "zero_func"));

      // temporary storage for formula
      std::vector<std::string> formulae, tformulae;

      bool make_newline(false);

      //loop over orbits
      for(Index no = 0; no < _tree.size(); no++) {

        if(_tree[no].prototype().size() == 0)
          bfunc_imp_stream <<
                           indent << "// Basis functions for empty cluster:\n";
        else {
          bfunc_imp_stream <<
                           indent << "\/**** Basis functions for orbit " << no << "****\n";
          _tree[no].prototype().print(bfunc_imp_stream, '\n');
          bfunc_imp_stream << "****/\n";
        }

        auto orbit_method_namer = [lf, no, &orbit_method_names](Index nb, Index nf)->std::string{
          orbit_method_names[lf + nf] = "eval_bfunc_" + std::to_string(orbit_ind) + "_" + std::to_string(nf);
          return orbit_method_names[nb][lf + nf];
        }

        std::tuple<std::string, std::string> orbit_functions = clexulator_orbit_function_strings(m_bset_tree[no],
                                                               _tree[no],
                                                               lf,
                                                               _nlist,
                                                               labelers);

        bfunc_def_stream << orbit_function<func_definition>();
        bfunc_imp_stream << orbit_function<func_implementation>();

        auto flower_method_namer = [lf, no, &flower_method_names](Index nb, Index nf)->std::string{
          flower_method_names[nb][lf + nf] = "site_eval_bfunc_"  + std::to_string(orbit_ind) + "_" + std::to_string(nf) + "_at_" + std::to_string(nb);
          return flower_method_names[nb][lf + nf];
        }

        std::tuple<std::string, std::string> flower_functions = clexulator_flower_function_strings(m_bset_tree[no],
                                                                _tree[no],
                                                                lf,
                                                                _nlist,
                                                                labelers);

        bfunc_def_stream << flower_function<func_definition>();
        bfunc_imp_stream << flower_function<func_implementation>();

        auto dflower_method_namer = [lf, no, &dflower_method_names](Index nb, Index nf)->std::string{
          dflower_method_names[nb][lf + nf] = "site_eval_bfunc_"  + std::to_string(orbit_ind) + "_" + std::to_string(nf) + "_at_" + std::to_string(nb);
          return dflower_method_names[nb][lf + nf];
        }

        OccFuncLabeler site_prefactor_labeler("(m_occ_func_%b_%f[occ_f] - m_occ_func_%b_%f[occ_i])");

        std::tuple<std::string, std::string> dflower_functions = clexulator_dflower_function_strings(m_bset_tree[no],
                                                                 _tree[no],
                                                                 lf,
                                                                 _nlist,
                                                                 labelers,
                                                                 site_prefactor_labeler);


        bfunc_def_stream << dflower_function<func_definition>();
        bfunc_imp_stream << dflower_function<func_implementation>();

        lf += m_bset_tree[no][0].size();
      }//Finished writing method definitions and implementations for basis functions


      std::string interface_implementation = clexulator_interface_implementation(class_name, _nlist, indent);

      // PUT EVERYTHING TOGETHER
      stream <<
             "#include <cstddef>\n" <<
             "#include \"casm/clex/Clexulator.hh\"\n" <<
             "\n\n\n" <<
             "\/****** CLEXULATOR CLASS FOR PRIM ******" << std::endl;

      jsonParser json;
      write_prim(prim, json, FRAC);
      stream << json;

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
}
