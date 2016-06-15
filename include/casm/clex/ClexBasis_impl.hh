#include "casm/basis_set/DoFTraits.hh"
#include "casm/crystallography/UnitCellCoord.hh"
#include "casm/clusterography/IntegralCluster.hh"
#include "casm/clex/NeighborList.hh"

namespace CASM {
  template<typename OrbitIteratorType>
  void ClexBasis::generate(OrbitIteratorType _orbit_begin,
                           OrbitIteratorType _orbit_end,
                           jsonParser const &_bspecs,
                           Index _max_poly_order /*= -1*/) {
    std::vector<DoFKey> dof_keys;
    _bspecs.get_if(dof_keys, "dofs");
    std::vector<DoFKey> global_keys;
    std::vector<DoFKey> local_keys;

    //separate local_args from global_args
    for(DoFKey const &key : dof_keys) {
      if(m_global_bases.find(key) != m_global_bases.end()) {
        global_keys.push_back(key);
      }
      else if(m_site_bases.find(key) != m_site_bases.end()) {
        local_keys.push_back(key);
      }
      else {
        throw std::runtime_error(std::string("Attempting to build Clex basis set, but missing degree of freedom \"") + key + "\n");
      }
    }
    m_bset_tree.resize(std::distance(_orbit_begin, _orbit_end));

    auto bset_it = m_bset_tree.begin();
    for(; _orbit_begin != _orbit_end; ++_orbit_begin, ++bset_it) {
      bset_it->reserve(_orbit_begin->size());
      bset_it->push_back(_construct_prototype_basis(*_orbit_begin,
                                                    local_keys,
                                                    global_keys,
                                                    -1/* polynomial_order */));
      for(Index j = 1; j < _orbit_begin->size(); j++) {
        bset_it->push_back((*(_orbit_begin->equivalence_map(j).first)) * (*bset_it)[0]);
      }
    }
  }


  //*************************************************
  // @param local_args[i][j] is BasisSet for i'th DoFspace at j'th site of cluster
  //                                site  0            site 1            site 2
  // DoFKey 0: displacement       {x0,y0,z0}       {x1, y1, z1}       {x2,y2,z2}
  // DoFKey 1: configuration       {pA0,pB0}            {}            {pA2, pB2}
  //
  // Step 1:  Get the kroenecker product of cluster permutation with DoF symrep
  //
  //  permutation  |    kronecker prod  |    DoF Symrep (e.g., x--y displacement)
  //   [ 0  1 ]              v                     [cos -sin]
  //   [ 1  0 ]             XkX                    [sin  cos]
  //
  //       [x0]      [  0    0   cos -sin ]   [x0]
  //       [y0]      [  0    0   sin  cos ]   [y0]
  //   S * [x1]  =   [ cos -sin   0    0  ]   [x1]
  //       [x2]      [ sin  cos   0    0  ]   [x2]
  //
  // ----------------------------------------------------------------------
  //
  // Step 2: mix-in @param global_args to get all_argsets
  //
  // GLOBAL ARGS
  //
  // strain             { e1, e2, e3, e4, e5, e6}
  // composition        { comp_a, comp_b }
  //
  // arg_subsets =   [ {x0,y0,z0,x1,y1,z1,x2,y2,z2},
  //                   {pA0,pB0,pA2,pB2},
  //                   {e1,e2,e3,e4,e5,e6},
  //                   {comp_a,comp_b}]
  //
  template<typename OrbitType>
  BasisSet ClexBasis::_construct_prototype_basis(OrbitType const &_orbit,
                                                 std::vector<DoFKey> const &local_keys,
                                                 std::vector<DoFKey> const &global_keys,
                                                 Index max_poly_order) const {
    //std::cout<<"In IntegralCluster::generate_clust_basis, the size of this cluster is:"<<size()<<std::endl;
    //std::cout<<"valid_index evaluates to:"<<valid_index(max_poly_order)<<std::endl;

    // Default polynomial order is cluster size
    if(!valid_index(max_poly_order))
      max_poly_order = _orbit.prototype().size();

    //std::cout<<"Max_poly_order "<<max_poly_order<<std::endl;

    std::vector<BasisSet const *> arg_subsets;
    for(DoFKey const &key : global_keys) {
      auto find_it = m_global_bases.find(key);
      if(find_it != m_global_bases.end())
        arg_subsets.push_back(&(find_it->second));
      else
        throw std::runtime_error("Unable to construct basis sets. No known global DoF: " + key + "\n");
    }

    // copy local site bases to
    std::vector<BasisSet> all_local;
    all_local.reserve(local_keys.size());

    //Loop over dof's
    for(DoFKey const &key : local_keys) {
      // Make copies of local arguments to ensure that they are distinguishable by their DoF_IDs
      // i.e., make copies in 'tlocal' and reset the DoF_IDs to {0,1,2,etc...}
      auto find_it = m_site_bases.find(key);
      if(find_it == m_site_bases.end())
        throw std::runtime_error("Unable to construct basis sets. No known local DoF: " + key + "\n");

      std::vector<BasisSet> const &arg_vec(find_it->second);
      std::vector<BasisSet> tlocal;
      tlocal.reserve(_orbit.prototype().size());
      std::vector<BasisSet const *> site_args(_orbit.prototype().size(), nullptr);
      //Loop over sites
      for(Index i = 0; i < _orbit.prototype().size(); i++) {
        if(arg_vec[_orbit.prototype()[i].sublat()].size()) {
          tlocal.push_back(arg_vec[_orbit.prototype()[i].sublat()]);
          tlocal.back().set_dof_IDs(std::vector<Index>(1, i));
          site_args[i] = &tlocal.back();
        }
      }
      all_local.push_back(ClexBasis_impl::construct_clust_dof_basis(_orbit.prototype(), site_args));
      if(all_local.back().size())
        arg_subsets.push_back(&(all_local.back()));
    }

    return m_basis_builder->build(_orbit.prototype(), arg_subsets, max_poly_order, 1);
  }

  //********************************************************************

  // Divide by multiplicity. Same result as evaluating correlations via orbitree.
  template<typename OrbitType>
  std::vector<std::string> orbit_function_cpp_strings(ClexBasis::BSetOrbit _bset_orbit, // used as temporary
                                                      OrbitType const &_clust_orbit,
                                                      PrimNeighborList &_nlist,
                                                      std::vector<FunctionVisitor *> const &labelers) {
    std::string prefix, suffix;
    std::vector<std::string> formulae(_bset_orbit[0].size(), std::string());
    if(_clust_orbit.size() > 1) {
      prefix = "(";
      suffix = ")/" + std::to_string(_clust_orbit.size()) + ".";
    }

    for(Index ne = 0; ne < _bset_orbit.size(); ne++) {
      std::vector<PrimNeighborList::Scalar> nbor_IDs =
        _nlist.neighbor_indices(_clust_orbit[ne].elements().begin(),
                                _clust_orbit[ne].elements().end());
      _bset_orbit[ne].set_dof_IDs(std::vector<Index>(nbor_IDs.begin(), nbor_IDs.end()));
      for(Index nl = 0; nl < labelers.size(); nl++)
        _bset_orbit[ne].accept(*labelers[nl]);
    }

    for(Index nf = 0; nf < _bset_orbit[0].size(); nf++) {
      for(Index ne = 0; ne < _bset_orbit.size(); ne++) {
        if(!_bset_orbit[ne][nf] || (_bset_orbit[ne][nf]->is_zero()))
          continue;

        if(formulae[nf].empty())
          formulae[nf] += prefix;
        else if((_bset_orbit[ne][nf]->formula())[0] != '-' && (_bset_orbit[ne][nf]->formula())[0] != '+')
          formulae[nf] += " + ";

        formulae[nf] += _bset_orbit[ne][nf]->formula();
      }

      if(!formulae[nf].empty())
        formulae[nf] += suffix;
    }
    return formulae;
  }

  //********************************************************************
  /// nlist_index is the index of the basis site in the neighbor list
  template<typename OrbitType>
  std::map< UnitCell, std::vector< std::string > > flower_function_cpp_strings(ClexBasis::BSetOrbit _bset_orbit, // used as temporary
                                                                               std::function<BasisSet(BasisSet const &)> _bset_transform,
                                                                               OrbitType const &_clust_orbit,
                                                                               PrimNeighborList &_nlist,
                                                                               std::vector<FunctionVisitor *> const &labelers,
                                                                               Index sublat_index) {


    typedef std::vector<std::string> Formulae;
    typedef std::map<UnitCell, std::vector<std::string> > TransFormulae;
    TransFormulae result;

    std::string prefix, suffix;
    std::set<UnitCellCoord> trans_set;
    for(IntegralCluster const &equiv : _clust_orbit) {
      for(UnitCellCoord const &site : equiv.elements()) {
        if(site.sublat() == sublat_index)
          trans_set.insert(site);
      }
    }

    std::map<UnitCellCoord, std::set<UnitCellCoord> > unique_trans = ClexBasis_impl::unique_ucc(trans_set.begin(),
                                                                     trans_set.end(),
                                                                     _clust_orbit.sym_compare());

    //normalize by multiplicity (by convention)
    if(_clust_orbit.size() > 1) {
      prefix = "(";
      suffix = ")/" + std::to_string(_clust_orbit.size()) + ".";
    }


    for(std::pair<UnitCellCoord, std::set<UnitCellCoord> > const &trans_orbit : unique_trans) {

      Formulae &formulae(result.emplace(std::make_pair(trans_orbit.first, Formulae(_bset_orbit[0].size()))).first->second);

      // loop over equivalent clusters
      for(Index ne = 0; ne < _clust_orbit.size(); ne++) {

        // loop over cluster sites
        auto it(_clust_orbit[ne].cbegin()), end_it(_clust_orbit[ne].cbegin());
        for(; it != end_it; ++it) {

          // Continue if the cluster site doesn't belong to the target sublattice, or if we are working on a different translation
          if(sublat_index != it -> sublat() || trans_orbit.second.find(*it) == trans_orbit.second.end())
            continue;

          IntegralCluster trans_clust = _clust_orbit[ne] - it->unitcell();

          std::vector<PrimNeighborList::Scalar> nbor_IDs =
            _nlist.neighbor_indices(trans_clust.elements().begin(),
                                    trans_clust.elements().end());
          _bset_orbit[ne].set_dof_IDs(std::vector<Index>(nbor_IDs.begin(), nbor_IDs.end()));

          BasisSet transformed_bset(_bset_transform(_bset_orbit[ne]));
          for(Index nl = 0; nl < labelers.size(); nl++)
            transformed_bset.accept(*labelers[nl]);

          for(Index nf = 0; nf < transformed_bset.size(); nf++) {
            if(!transformed_bset[nf] || (transformed_bset[nf]->is_zero()))
              continue;

            if(formulae[nf].empty())
              formulae[nf] += prefix;
            else if((transformed_bset[nf]->formula())[0] != '-' && (transformed_bset[nf]->formula())[0] != '+')
              formulae[nf] += " + ";

            formulae[nf] += transformed_bset[nf]->formula();

          }
        }
      }
    }

    // append suffix to all formulae
    TransFormulae::iterator trans_it = result.begin(),
                            trans_end = result.end();
    for(; trans_it != trans_end; ++trans_it) {
      Formulae &formulae(trans_it->second);
      for(Index nf = 0; nf < formulae.size(); nf++) {
        if(!formulae[nf].empty())
          formulae[nf] += suffix;
      }
    }

    return result;

  }

  //***********************************************
  template<typename OrbitType>
  void print_proto_clust_funcs(ClexBasis const &_clex_basis,
                               std::ostream &out,
                               BasicStructure<Site> const &_prim,
                               std::vector<OrbitType > const &_tree) {
    //Prints out all prototype clusters (CLUST file)

    out << "COORD_MODE = " << COORD_MODE::NAME() << std::endl << std::endl;

    auto const &site_bases(_clex_basis.site_bases());

    out.flags(std::ios::showpoint | std::ios::fixed | std::ios::left);
    out.precision(5);
    auto it = site_bases.begin(), end_it = site_bases.end();
    out << "Basis site definitions for DoF " << it->first << ".\n";
    Index nb;
    for(Index nb = 0; nb < (it->second).size(); ++nb) {
      out << "  Basis site " << nb + 1 << ":\n"
          << "  ";
      _prim.basis[nb].print(out);
      out << "\n";
      //HERE
      out << DoFType::traits(it->first).site_basis_description((it->second)[nb], _prim.basis[nb]);
    }

    out << "\n\n";
    Index nf = 0;
    for(Index i = 0; i < _tree.size(); i++) {
      if(i == 0 || _tree[i].prototype().size() != _tree[i - 1].prototype().size())
        out << "** " << _tree[i].prototype().size() << "-site clusters ** \n" << std::flush;

      out << "      ** Orbit " << i + 1 << " of " << _tree.size() << " **"
          << "  Points: " << _tree[i].prototype().size()
          << "  Mult: " << _tree[i].size()
          << "  MinLength: " << _tree[i].prototype().min_length()
          << "  MaxLength: " << _tree[i].prototype().max_length()
          << '\n' << std::flush;

      out << "            " << "Prototype" << " of " << _tree[i].size() << " Equivalent Clusters in Orbit " << i << '\n' << std::flush;
      print_clust_basis(out,
                        _clex_basis.bset_orbit(i)[0],
                        _tree[i].prototype(),
                        nf, 8, '\n');

      nf += _clex_basis.bset_orbit(i)[0].size();
      out << "\n\n" << std::flush;

      if(_tree[i].prototype().size() != 0) out << '\n' << std::flush;
    }
  }

  //*******************************************************************************************
  /// \brief Print clexulator
  /*
  void print_clexulator(const Structure &prim,
                        SiteOrbitree &tree,
                        const PrimNeighborList &nlist,
                        std::string class_name,
                        std::ostream &stream,
                        double xtal_tol) {

    set_nlist_ind(prim, tree, nlist, xtal_tol);

    DoFManager dof_manager;
    Index Nsublat = prim.basis.size();
    for(Index b = 0; b < Nsublat; b++) {
      if(prim.basis[b].site_occupant().size() > 1) {
        dof_manager.add_dof(prim.basis[b].site_occupant().type_name());
        break;
      }
    }
    for(Index b = 0; b < Nsublat; b++) {
      for(Index i = 0; i < prim.basis[i].displacement().size(); i++)
        dof_manager.add_dof(prim.basis[b].displacement()[i].type_name());
    }

    dof_manager.resize_neighborhood(nlist.size()*nlist.sublat_indices().size());

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
                       indent << "  BasisFuncPtr m_flower_func_lists[" << Nsublat << "][" << N_corr << "];\n\n" <<

                       //for separate 1D method pointer lists:
                       //indent << "  BasisFuncPtr
                       //for(Index i = 0; i < Nsublat; i++) {
                         //private_def_stream << " m_flower_func_at_" << i << "_list[" << N_corr << "]";
                         //if(i + 1 < Nsublat)
                           //private_def_stream << ',';
                       //}
                       //private_def_stream << ";\n\n" <<


                       indent << "  // array of pointers to member functions for calculating DELTA flower functions\n" <<
                       indent << "  DeltaBasisFuncPtr m_delta_func_lists[" << Nsublat << "][" << N_corr << "];\n\n";

  // for separate 1D method pointer lists:
  //indent << "  DeltaBasisFuncPtr";

  //for(Index i = 0; i < Nsublat; i++) {
  //private_def_stream << " m_delta_func_at_" << i << "_list[" << N_corr << "]";
  //if(i + 1 < Nsublat)
  //private_def_stream << ',';
  //}
  //private_def_stream << ";\n\n";


    dof_manager.print_clexulator_member_definitions(private_def_stream, tree, indent + "  ");

    // perhaps some more correlation calculating options here
    dof_manager.print_clexulator_private_method_definitions(private_def_stream, tree, indent + "  ");

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

    dof_manager.print_clexulator_public_method_definitions(public_def_stream, tree, indent + "  ");


    //linear function index
    Index lf = 0, tlf;

    Array<FunctionVisitor *> labelers(dof_manager.get_function_label_visitors());
    //std::cout << "Initialized " << labelers.size() << " labelers \n";

    Array<std::string> orbit_method_names(N_corr);
    Array<Array<std::string> > flower_method_names(Nsublat, Array<std::string>(N_corr));
    //Array< Array<Array<std::string> > > dflower_method_names(N_corr, Array<Array<std::string> >(Nsublat));

    //this is very configuration-centric
    Array<Array<std::string> > dflower_method_names(Nsublat, Array<std::string>(N_corr));

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
                           indent << "\/**** Basis functions for orbit " << np << ", " << no << "****\n";
          tree[np][no].prototype.print(bfunc_imp_stream, '\n');
          bfunc_imp_stream << "\*\*\*\*\/\n";
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
        const SiteOrbitBranch &asym_unit(tree.asym_unit());
        for(Index na = 0; na < asym_unit.size(); na++) {
          for(Index ne = 0; ne < asym_unit[na].size(); ne++) {
            Index nb = asym_unit[na][ne][0].basis_ind();
            auto nlist_index = find_index(nlist.sublat_indices(), nb);
            if(nlist_index != nlist.sublat_indices().size()) {
              formulae = tree[np][no].flower_function_cpp_strings(labelers, nlist_index);
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

              }
            }
            if(make_newline) {
              bfunc_imp_stream << '\n';
              private_def_stream << '\n';
            }
            make_newline = false;

            // Very configuration-centric -> Find a way to move this block to OccupationDoFEnvironment:
            formulae.resize(formulae.size(), std::string());
            // loop over site basis functions
            const BasisSet &site_basis(asym_unit[na][ne].clust_basis);
            for(Index nsbf = 0; nsbf < site_basis.size(); nsbf++) {
              std::string delta_prefix = "(m_occ_func_" + std::to_string(nb) + "_" + std::to_string(nsbf) + "[occ_f] - m_occ_func_" + std::to_string(nb) + "_" + std::to_string(nsbf) + "[occ_i])";

              if(nlist_index != nlist.sublat_indices().size()) {
                tformulae = tree[np][no].delta_occfunc_flower_function_cpp_strings(site_basis, labelers, nlist_index, nb, nsbf);
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
            // \End Configuration specific part
          }
        }//\End loop over flowers

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

    dof_manager.print_to_clexulator_constructor(interface_imp_stream, tree, indent + "  ");

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

    // Write weight matrix used for the neighbor list
    PrimNeighborList::Matrix3Type W = nlist.weight_matrix();
    interface_imp_stream << indent << "  m_weight_matrix.row(0) << " << W(0, 0) << ", " << W(0, 1) << ", " << W(0, 2) << ";\n";
    interface_imp_stream << indent << "  m_weight_matrix.row(1) << " << W(1, 0) << ", " << W(1, 1) << ", " << W(1, 2) << ";\n";
    interface_imp_stream << indent << "  m_weight_matrix.row(2) << " << W(2, 0) << ", " << W(2, 1) << ", " << W(2, 2) << ";\n\n";

    // Write neighborhood of UnitCellCoord
    // expand the nlist to contain 'global_orbitree' (all that is needed for now)
    std::set<UnitCellCoord> nbors;
    neighborhood(std::inserter(nbors, nbors.begin()), tree, prim, TOL);


  //for(auto it = nbors.begin(); it != nbors.end(); ++it) {
  // interface_imp_stream << indent << "  m_neighborhood.insert(UnitCellCoord("
  //<< it->sublat() << ", "
  //                       << it->unitcell(0) << ", "
  //                       << it->unitcell(1) << ", "
  //                       << it->unitcell(2) << "));\n";
  //}


    interface_imp_stream << indent << "  m_neighborhood = std::set<UnitCellCoord> {\n";
    auto it = nbors.begin();
    while(it != nbors.end()) {
      interface_imp_stream << indent << "    {UnitCellCoord("
                           << it->sublat() << ", "
                           << it->unitcell(0) << ", "
                           << it->unitcell(1) << ", "
                           << it->unitcell(2) << ")}";
      ++it;
      if(it != nbors.end()) {
        interface_imp_stream << ",";
      }
      interface_imp_stream << "\n";
    }
    interface_imp_stream << indent << "  };\n";
    interface_imp_stream << "\n\n";

    interface_imp_stream << indent <<  "  m_orbit_neighborhood.resize(corr_size());\n";
    Index lno = 0;
    for(Index nb = 0; nb < tree.size(); ++nb) {
      for(Index no = 0; no < tree[nb].size(); ++no) {
        std::set<UnitCellCoord> orbit_nbors;
        orbit_neighborhood(std::inserter(orbit_nbors, orbit_nbors.begin()), tree, prim, nb, no, TOL);

        Index proto_index = lno;

        //for(auto it = orbit_nbors.begin(); it != orbit_nbors.end(); ++it) {
        //interface_imp_stream << indent << "  m_orbit_neighborhood[" << lno << "].insert(UnitCellCoord("
        //                     << it->sublat() << ", "
        //                     << it->unitcell(0) << ", "
        //                     << it->unitcell(1) << ", "
        //                     << it->unitcell(2) << "));\n";
        //}

        interface_imp_stream << indent << "  m_orbit_neighborhood[" << lno << "] = std::set<UnitCellCoord> {\n";
        auto it = orbit_nbors.begin();
        while(it != orbit_nbors.end()) {
          interface_imp_stream << indent << "    {UnitCellCoord("
                               << it->sublat() << ", "
                               << it->unitcell(0) << ", "
                               << it->unitcell(1) << ", "
                               << it->unitcell(2) << ")}";
          ++it;
          if(it != orbit_nbors.end()) {
            interface_imp_stream << ",";
          }
          interface_imp_stream << "\n";
        }
        interface_imp_stream << indent << "  };\n";
        ++lno;
        for(Index nf = 1; nf < tree.prototype(nb, no).clust_basis.size(); ++nf) {
          interface_imp_stream << indent << "  m_orbit_neighborhood[" << lno << "] = m_orbit_neighborhood[" << proto_index << "];\n";
          ++lno;
        }
        interface_imp_stream << "\n";

      }
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
           "\/****** CLEXULATOR CLASS FOR PRIM ******" << std::endl;

    jsonParser json;
    write_prim(prim, json, FRAC);
    stream << json;

    stream <<
           "\*\*\/\n\n\n" <<

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
  */

  //***********************************************
  namespace ClexBasis_impl {
    template<typename UCCIterType, typename IntegralClusterSymCompareType>
    std::map<UnitCellCoord, std::set<UnitCellCoord> > unique_ucc(UCCIterType begin,
                                                                 UCCIterType end,
                                                                 IntegralClusterSymCompareType const &sym_compare) {
      std::map<UnitCellCoord, std::set<UnitCellCoord> > tresult;

      typedef IntegralCluster cluster_type;
      typedef Orbit<cluster_type, IntegralClusterSymCompareType> orbit_type;

      if(begin == end)
        return tresult;
      // store orbits as we find them
      std::set<orbit_type> orbits;

      SymGroup identity_group(begin->prim().factor_group.begin(), begin->prim().factor_group.begin()++);
      orbit_type empty_orbit(cluster_type(begin->prim()), identity_group, sym_compare);

      // by looping over each site in the grid,
      for(; begin != end; ++begin) {

        // create a test cluster from prototype
        cluster_type test(empty_orbit.prototype());

        // add the new site
        test.elements().push_back(*begin);


        // try to find test cluster in already found orbits
        auto it = find_orbit(orbits.begin(), orbits.end(), test);
        if(it != orbits.end()) {
          tresult[it->prototype().element[0]].insert(*begin);
          continue;
        }

        tresult[test.element(0)].insert(*begin);

        // if not yet found, use test to generate a new Orbit
        orbits.insert(orbit_type(test, identity_group, sym_compare));
      }

      std::map<UnitCellCoord, std::set<UnitCellCoord> > result;
      for(auto &_set_pair : tresult)
        result.emplace(std::make_pair(*_set_pair.second.begin(), std::move(_set_pair.second)));
      return result;

    }
  }
}
