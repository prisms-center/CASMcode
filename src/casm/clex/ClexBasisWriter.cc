#include "casm/basis_set/DoFTraits.hh"
#include "casm/basis_set/FunctionVisitor.hh"
#include "casm/clex/ClexBasis.hh"
#include "casm/clex/ClexBasisWriter_impl.hh"
#include "casm/clex/OrbitFunctionTraits.hh"
#include "casm/crystallography/BasicStructure.hh"

namespace CASM {

ClexBasisWriter::ClexBasisWriter(Structure const &_prim,
                                 PARAM_PACK_TYPE const &param_pack_type) {
  if (param_pack_type == PARAM_PACK_TYPE::DEFAULT) {
    _initialize(_prim, ParamPackMixIn::basic_mix_in());
  } else if (param_pack_type == PARAM_PACK_TYPE::DIFF) {
    _initialize(_prim, ParamPackMixIn::diff_mix_in());
  } else {
    throw std::runtime_error(
        "Error in ClexBasisWriter: Invalid param_pack_type");
  }
}

ClexBasisWriter::ClexBasisWriter(Structure const &_prim,
                                 ParamPackMixIn const &_param_pack_mix_in) {
  _initialize(_prim, _param_pack_mix_in);
}

void ClexBasisWriter::_initialize(Structure const &_prim,
                                  ParamPackMixIn const &_param_pack_mix_in) {
  m_param_pack_mix_in = _param_pack_mix_in.clone();

  auto doftypes = xtal::all_local_dof_types(_prim);
  for (auto const &doftype : doftypes) {
    auto cv = DoFType::traits(doftype).clust_function_visitors();
    auto sv = DoFType::traits(doftype).site_function_visitors();
    for (auto &e : cv) m_clust_visitors.push_back(std::move(e));
    for (auto &e : sv) m_site_visitors.push_back(std::move(e));
  }

  doftypes = xtal::global_dof_types(_prim);
  for (auto const &doftype : doftypes) {
    auto cv = DoFType::traits(doftype).clust_function_visitors();
    auto sv = DoFType::traits(doftype).site_function_visitors();
    for (auto &e : cv) m_clust_visitors.push_back(std::move(e));
    for (auto &e : sv) m_site_visitors.push_back(std::move(e));
  }
}

namespace ClexBasisWriter_impl {

std::string clexulator_member_declarations(
    std::string const &class_name, ClexBasis const &clex,
    ParamPackMixIn const &_param_pack_mix_in,
    std::vector<std::unique_ptr<OrbitFunctionTraits> > const
        &orbit_func_writers,
    Index N_flower, std::string const &indent) {
  Index N_corr = clex.n_functions();

  std::stringstream ss;

  for (auto const &writer_ptr : orbit_func_writers) {
    writer_ptr->print_typedefs(ss, class_name, indent);
  }

  for (auto const &writer_ptr : orbit_func_writers) {
    writer_ptr->print_eval_table_declarations(ss, class_name, clex, indent);
  }

  ss <<

      indent
     << "// ParamPack object, which stores temporary data for calculations\n"
     << indent << "mutable ParamPack m_params;\n\n";

  Index ispec = 0;
  for (auto const &specialization :
       _param_pack_mix_in.scalar_specializations()) {
    ss <<

        indent << "// typedef for method pointers of scalar type "
       << specialization.second << "\n"
       << indent << "typedef " << specialization.second << " (" << class_name
       << "::*BasisFuncPtr_" << ispec << ")() const;\n\n"
       <<

        indent << "// typedef for method pointers\n"
       << indent << "typedef " << specialization.second << " (" << class_name
       << "::*DeltaBasisFuncPtr_" << ispec << ")(int, int) const;\n\n"
       <<

        indent
       << "// array of pointers to member functions for calculating basis "
          "functions of scalar type "
       << specialization.second << "\n"
       << indent << "BasisFuncPtr_" << ispec << " m_orbit_func_table_" << ispec
       << "[" << N_corr << "];\n\n"
       <<

        indent
       << "// array of pointers to member functions for calculating flower "
          "functions of scalar type "
       << specialization.second << "\n"
       << indent << "BasisFuncPtr_" << ispec << " m_flower_func_table_" << ispec
       << "[" << N_flower << "][" << N_corr << "];\n\n"
       <<

        indent
       << "// array of pointers to member functions for calculating DELTA "
          "flower functions of scalar type "
       << specialization.second << "\n"
       << indent << "DeltaBasisFuncPtr_" << ispec << " m_delta_func_table_"
       << ispec << "[" << N_flower << "][" << N_corr << "];\n\n";

    ++ispec;
  }

  for (auto const &dof : clex.site_bases())
    ss << DoFType::traits(dof.first).clexulator_member_declarations_string(
        clex.prim(), dof.second, indent);

  for (auto const &dof : clex.global_bases())
    ss << DoFType::traits(dof.first).clexulator_member_declarations_string(
        clex.prim(), dof.second, indent);

  ss << indent << "//ClexParamPack allocation for evaluated correlations \n"
     << indent << "ParamPack::Key m_corr_param_key;\n";

  for (auto const &dof : clex.site_bases()) {
    std::vector<DoFType::ParamAllocation> allo =
        DoFType::traits(dof.first).param_pack_allocation(clex.prim(),
                                                         dof.second);
    if (allo.empty()) continue;
    ss << indent << "//ClexParamPack allocation for DoF " << dof.first << "\n";
    for (const auto &el : allo)
      ss << indent << "ParamPack::Key m_" << el.param_name << "_param_key;\n";
    ss << "\n";
  }

  for (auto const &dof : clex.global_bases()) {
    std::vector<DoFType::ParamAllocation> allo =
        DoFType::traits(dof.first).param_pack_allocation(clex.prim(),
                                                         dof.second);
    if (allo.empty()) continue;
    ss << indent << "//ClexParamPack allocation for DoF " << dof.first << "\n";
    for (const auto &el : allo)
      ss << indent << "ParamPack::Key m_" << el.param_name << "_param_key;\n";
    ss << "\n";
  }

  return ss.str();
}

//*******************************************************************************************

std::string clexulator_private_method_declarations(
    std::string const &class_name, ClexBasis const &clex,
    std::string const &indent) {
  std::stringstream ss;
  ss <<

      indent << "/// \\brief Clone the " << class_name << "\n"
     << indent << "Clexulator_impl::Base *_clone() const override {\n"
     << indent << "  return new " << class_name << "(*this);\n"
     << indent << "}\n\n"
     <<

      indent
     << "/// \\brief Calculate contribution to global correlations from one "
        "unit cell\n"
     << indent << "/// Result is recorded in ClexParamPack\n"
     << indent << "void _calc_global_corr_contribution() const override;\n\n"
     <<

      indent
     << "/// \\brief Calculate contribution to global correlations from one "
        "unit cell "
     << indent
     << "/// Result is recorded in double array starting at corr_begin\n"
     << indent
     << "void _calc_global_corr_contribution(double *corr_begin) const "
        "override;\n\n"
     <<

      indent
     << "/// \\brief Calculate contribution to select global correlations from "
        "one unit cell into ClexParamPack\n"
     << indent << "/// Result is recorded in ClexParamPack\n"
     << indent
     << "void _calc_restricted_global_corr_contribution(size_type const "
        "*ind_list_begin, size_type const *ind_list_end) const override;\n\n"
     <<

      indent
     << "/// \\brief Calculate contribution to select global correlations from "
        "one unit cell\n"
     << indent
     << "/// Result is recorded in double array starting at corr_begin\n"
     << indent
     << "void _calc_restricted_global_corr_contribution(double *corr_begin, "
        "size_type const *ind_list_begin, size_type const *ind_list_end) const "
        "override;\n\n"
     <<

      indent
     << "/// \\brief Calculate point correlations about neighbor site "
        "'nlist_ind'\n"
     << indent
     << "/// For global clexulators, 'nlist_ind' only ranges over sites in the "
        "cell\n"
     << indent
     << "/// For local clexulators, 'nlist_ind' ranges over all sites in the "
        "neighborhood\n"
     << indent << "/// Result is recorded in ClexParamPack\n"
     << indent << "void _calc_point_corr(int nlist_ind) const override;\n\n"
     <<

      indent
     << "/// \\brief Calculate point correlations about neighbor site "
        "'nlist_ind'\n"
     << indent
     << "/// For global clexulators, 'nlist_ind' only ranges over sites in the "
        "cell\n"
     << indent
     << "/// For local clexulators, 'nlist_ind' ranges over all sites in the "
        "neighborhood\n"
     << indent
     << "/// Result is recorded in double array starting at corr_begin\n"
     << indent
     << "void _calc_point_corr(int nlist_ind, double *corr_begin) const "
        "override;\n\n"
     <<

      indent
     << "/// \\brief Calculate select point correlations about neighbor site "
        "'nlist_ind'\n"
     << indent
     << "/// For global clexulators, 'nlist_ind' only ranges over sites in the "
        "cell\n"
     << indent
     << "/// For local clexulators, 'nlist_ind' ranges over all sites in the "
        "neighborhood\n"
     << indent << "/// Result is recorded in ClexParamPack\n"
     << indent
     << "void _calc_restricted_point_corr(int nlist_ind, size_type const "
        "*ind_list_begin, size_type const *ind_list_end) const override;\n\n"
     <<

      indent
     << "/// \\brief Calculate select point correlations about neighbor site "
        "'nlist_ind'\n"
     << indent
     << "/// For global clexulators, 'nlist_ind' only ranges over sites in the "
        "cell\n"
     << indent
     << "/// For local clexulators, 'nlist_ind' ranges over all sites in the "
        "neighborhood\n"
     << indent
     << "/// Result is recorded in double array starting at corr_begin\n"
     << indent
     << "void _calc_restricted_point_corr(int nlist_ind, double *corr_begin, "
        "size_type const *ind_list_begin, size_type const *ind_list_end) const "
        "override;\n\n"
     <<

      indent
     << "/// \\brief Calculate the change in point correlations due to "
        "changing an occupant at neighbor site 'nlist_ind'\n"
     << indent
     << "/// For global clexulators, 'nlist_ind' only ranges over sites in the "
        "cell\n"
     << indent
     << "/// For local clexulators, 'nlist_ind' ranges over all sites in the "
        "neighborhood\n"
     << indent << "/// Result is recorded in ClexParamPack\n"
     << indent
     << "void _calc_delta_point_corr(int nlist_ind, int occ_i, int occ_f) "
        "const override;\n\n"
     <<

      indent
     << "/// \\brief Calculate the change in point correlations due to "
        "changing an occupant at neighbor site 'nlist_ind'\n"
     << indent
     << "/// For global clexulators, 'nlist_ind' only ranges over sites in the "
        "cell\n"
     << indent
     << "/// For local clexulators, 'nlist_ind' ranges over all sites in the "
        "neighborhood\n"
     << indent
     << "/// Result is recorded in double array starting at corr_begin\n"
     << indent
     << "void _calc_delta_point_corr(int nlist_ind, int occ_i, int occ_f, "
        "double *corr_begin) const override;\n\n"
     <<

      indent
     << "/// \\brief Calculate the change in select point correlations due to "
        "changing an occupant at neighbor site 'nlist_ind'\n"
     << indent
     << "/// For global clexulators, 'nlist_ind' only ranges over sites in the "
        "cell\n"
     << indent
     << "/// For local clexulators, 'nlist_ind' ranges over all sites in the "
        "neighborhood\n"
     << indent << "/// Result is recorded in ClexParamPack\n"
     << indent
     << "void _calc_restricted_delta_point_corr(int nlist_ind, int occ_i, int "
        "occ_f, size_type const *ind_list_begin, size_type const "
        "*ind_list_end) const override;\n\n"
     <<

      indent
     << "/// \\brief Calculate the change in select point correlations due to "
        "changing an occupant at neighbor site 'nlist_ind'\n"
     << indent
     << "/// For global clexulators, 'nlist_ind' only ranges over sites in the "
        "cell\n"
     << indent
     << "/// For local clexulators, 'nlist_ind' ranges over all sites in the "
        "neighborhood\n"
     << indent
     << "/// Result is recorded in double array starting at corr_begin\n"
     << indent
     << "void _calc_restricted_delta_point_corr(int nlist_ind, int occ_i, int "
        "occ_f, double *corr_begin, size_type const *ind_list_begin, size_type "
        "const *ind_list_end) const override;\n\n"
     <<

      indent << "template<typename Scalar>\n"
     << indent << "void _global_prepare() const;\n\n"
     <<

      indent << "template<typename Scalar>\n"
     << indent << "void _point_prepare(int nlist_ind) const;\n\n";

  {
    auto it(clex.site_bases().begin()), end_it(clex.site_bases().end());
    for (; it != end_it; ++it) {
      ss << DoFType::traits(it->first)
                .clexulator_private_method_declarations_string(
                    clex.prim(), it->second, indent);
    }
  }

  {
    auto it(clex.global_bases().begin()), end_it(clex.global_bases().end());
    for (; it != end_it; ++it) {
      ss << DoFType::traits(it->first)
                .clexulator_private_method_declarations_string(
                    clex.prim(), it->second, indent);
    }
  }

  ss << indent << "//default functions for basis function evaluation\n"
     << indent << "template <typename Scalar>\n"
     << indent << "Scalar zero_func() const {\n"
     << indent << "  return Scalar(0.0);\n"
     << indent << "}\n\n"
     <<

      indent << "template <typename Scalar>\n"
     << indent << "Scalar zero_func(int, int) const {\n"
     << indent << "  return Scalar(0.0);\n"
     << indent << "}\n\n";

  return ss.str();
}

//*******************************************************************************************

std::string clexulator_public_method_declarations(std::string const &class_name,
                                                  ClexBasis const &clex,
                                                  std::string const &indent) {
  std::stringstream ss;
  ss << indent << class_name << "();\n\n"
     << indent << "~" << class_name << "();\n\n"
     <<

      indent << "ClexParamPack const &param_pack() const override {\n"
     << indent << "  return m_params;\n"
     << indent << "}\n\n"
     <<

      indent << "ClexParamPack &param_pack() override {\n"
     << indent << "  return m_params;\n"
     << indent << "}\n\n";

  {
    auto it(clex.site_bases().begin()), end_it(clex.site_bases().end());
    for (; it != end_it; ++it) {
      ss << DoFType::traits(it->first)
                .clexulator_public_method_declarations_string(
                    clex.prim(), it->second, indent);
    }
  }
  {
    auto it(clex.global_bases().begin()), end_it(clex.global_bases().end());
    for (; it != end_it; ++it) {
      ss << DoFType::traits(it->first)
                .clexulator_public_method_declarations_string(
                    clex.prim(), it->second, indent);
    }
  }
  return ss.str();
}

//*******************************************************************************************

std::string clexulator_interface_declaration(
    std::string const &class_name, ClexBasis const &clex,
    ParamPackMixIn const &_param_pack_mix_in, std::string const &indent) {
  std::stringstream ss;
  // Write destructor
  ss << indent << class_name << "::~" << class_name << "() {\n"
     << indent << "  //nothing here for now\n"
     << indent << "}\n\n";

  // Write evaluation methods
  auto specializations = _param_pack_mix_in.scalar_specializations();
  ss << indent
     << "/// \\brief Calculate contribution to global correlations from one "
        "unit cell\n"
     << indent << "void " << class_name
     << "::_calc_global_corr_contribution(double *corr_begin) const {\n"
     << indent << "  _calc_global_corr_contribution();\n"
     << indent << "  for(size_type i = 0; i < corr_size(); i++) {\n"
     << indent
     << "    *(corr_begin + i) = ParamPack::Val<double>::get(m_params, "
        "m_corr_param_key, i);\n"
     << indent << "  }\n"
     << indent << "}\n\n"
     <<

      //-----

      indent
     << "/// \\brief Calculate contribution to global correlations from one "
        "unit cell\n"
     << indent << "void " << class_name
     << "::_calc_global_corr_contribution() const {\n"
     << indent << "  m_params.pre_eval();\n";

  Index ispec = 0;

  for (auto const &specialization : specializations) {
    if (specializations.size() > 1) {
      ss << indent << "  ";
      if (ispec > 0) ss << "else ";
      ss << "if(m_params.eval_mode() == " << specialization.first << ")";
    }
    ss << indent << "  {\n"
       << indent << "    _global_prepare<" << specialization.second << ">();\n"
       << indent << "    for(size_type i = 0; i < corr_size(); i++) {\n"
       << indent << "      ParamPack::Val<" << specialization.second
       << ">::set(m_params, m_corr_param_key, i, (this->*m_orbit_func_table_"
       << ispec << "[i])());\n"
       << indent << "    }\n"
       << indent << "  }\n";
    ++ispec;
  }
  ss << indent << "  m_params.post_eval();\n"
     << indent << "}\n\n"
     <<

      //-----

      indent
     << "/// \\brief Calculate contribution to select global correlations from "
        "one unit cell\n"
     << indent << "void " << class_name
     << "::_calc_restricted_global_corr_contribution(double *corr_begin, "
        "size_type const *ind_list_begin, size_type const *ind_list_end) const "
        "{\n"
     << indent
     << "  _calc_restricted_global_corr_contribution(ind_list_begin, "
        "ind_list_end);\n"
     << indent << "  for(; ind_list_begin < ind_list_end; ind_list_begin++) {\n"
     << indent
     << "    *(corr_begin + *ind_list_begin) = "
        "ParamPack::Val<double>::get(m_params, m_corr_param_key, "
        "*ind_list_begin);\n"
     << indent << "  }\n"
     << indent << "}\n\n"
     <<

      //-----

      indent
     << "/// \\brief Calculate contribution to select global correlations from "
        "one unit cell\n"
     << indent << "void " << class_name
     << "::_calc_restricted_global_corr_contribution(size_type const "
        "*ind_list_begin, size_type const *ind_list_end) const {\n"
     << indent << "  m_params.pre_eval();\n";

  ispec = 0;
  for (auto const &specialization : specializations) {
    if (specializations.size() > 1) {
      ss << indent << "  ";
      if (ispec > 0) ss << "else ";
      ss << "if(m_params.eval_mode() == " << specialization.first << ")";
    }
    ss << indent << "  {\n"
       << indent << "    _global_prepare<" << specialization.second << ">();\n"
       << indent
       << "    for(; ind_list_begin < ind_list_end; ind_list_begin++) {\n"
       << indent << "      ParamPack::Val<" << specialization.second
       << ">::set(m_params, m_corr_param_key, *ind_list_begin, "
          "(this->*m_orbit_func_table_"
       << ispec << "[*ind_list_begin])());\n"
       << indent << "    }\n"
       << indent << "  }\n";
    ++ispec;
  }
  ss << indent << "  m_params.post_eval();\n"
     << indent << "}\n\n"
     <<

      //-----

      indent
     << "/// \\brief Calculate point correlations about basis site "
        "'nlist_ind'\n"
     << indent << "void " << class_name
     << "::_calc_point_corr(int nlist_ind, double *corr_begin) const {\n"
     << indent << "  _calc_point_corr(nlist_ind);\n"
     << indent << "  for(size_type i = 0; i < corr_size(); i++) {\n"
     << indent
     << "    *(corr_begin + i) = ParamPack::Val<double>::get(m_params, "
        "m_corr_param_key, i);\n"
     << indent << "  }\n"
     << indent << "}\n\n"
     <<

      //-----

      indent
     << "/// \\brief Calculate point correlations about basis site "
        "'nlist_ind'\n"
     << indent << "void " << class_name
     << "::_calc_point_corr(int nlist_ind) const {\n"
     << indent << "  m_params.pre_eval();\n";

  ispec = 0;
  for (auto const &specialization : specializations) {
    if (specializations.size() > 1) {
      ss << indent << "  ";
      if (ispec > 0) ss << "else ";
      ss << "if(m_params.eval_mode() == " << specialization.first << ")";
    }
    ss << indent << "  {\n"
       << indent << "    _point_prepare<" << specialization.second
       << ">(nlist_ind);\n"
       << indent << "    for(size_type i = 0; i < corr_size(); i++) {\n"
       << indent << "      ParamPack::Val<" << specialization.second
       << ">::set(m_params, m_corr_param_key, i, (this->*m_flower_func_table_"
       << ispec << "[nlist_ind][i])());\n"
       << indent << "    }\n"
       << indent << "  }\n";
    ++ispec;
  }
  ss << indent << "  m_params.post_eval();\n"
     << indent << "}\n\n"
     <<

      //-----

      indent
     << "/// \\brief Calculate select point correlations about basis site "
        "'nlist_ind'\n"
     << indent << "void " << class_name
     << "::_calc_restricted_point_corr(int nlist_ind, double *corr_begin, "
        "size_type const *ind_list_begin, size_type const *ind_list_end) const "
        "{\n"
     << indent
     << "  _calc_restricted_point_corr(nlist_ind, ind_list_begin, "
        "ind_list_end);\n"
     << indent << "  for(; ind_list_begin < ind_list_end; ind_list_begin++) {\n"
     << indent
     << "    *(corr_begin + *ind_list_begin) = "
        "ParamPack::Val<double>::get(m_params, m_corr_param_key, "
        "*ind_list_begin);\n"
     << indent << "  }\n"
     << indent << "}\n\n"
     <<

      //-----

      indent
     << "/// \\brief Calculate select point correlations about basis site "
        "'nlist_ind'\n"
     << indent << "void " << class_name
     << "::_calc_restricted_point_corr(int nlist_ind, size_type const "
        "*ind_list_begin, size_type const *ind_list_end) const {\n"
     << indent << "  m_params.pre_eval();\n";

  ispec = 0;
  for (auto const &specialization : specializations) {
    if (specializations.size() > 1) {
      ss << indent << "  ";
      if (ispec > 0) ss << "else ";
      ss << "if(m_params.eval_mode() == " << specialization.first << ")";
    }
    ss << indent << "  {\n"
       << indent << "    _point_prepare<" << specialization.second
       << ">(nlist_ind);\n"
       << indent
       << "    for(; ind_list_begin < ind_list_end; ind_list_begin++) {\n"
       << indent << "      ParamPack::Val<" << specialization.second
       << ">::set(m_params, m_corr_param_key, *ind_list_begin, "
          "(this->*m_flower_func_table_"
       << ispec << "[nlist_ind][*ind_list_begin])());\n"
       << indent << "    }\n"
       << indent << "  }\n";
    ++ispec;
  }
  ss << indent << "  m_params.post_eval();\n"
     << indent << "}\n\n"
     <<

      //-----

      indent
     << "/// \\brief Calculate the change in point correlations due to "
        "changing an occupant\n"
     << indent << "void " << class_name
     << "::_calc_delta_point_corr(int nlist_ind, int occ_i, int occ_f, double "
        "*corr_begin) const {\n"
     << indent << "  _calc_delta_point_corr(nlist_ind, occ_i, occ_f);\n"
     << indent << "  for(size_type i = 0; i < corr_size(); i++) {\n"
     << indent
     << "    *(corr_begin + i) = ParamPack::Val<double>::get(m_params, "
        "m_corr_param_key, i);\n"
     << indent << "  }\n"
     << indent << "}\n\n"
     <<

      //-----

      indent
     << "/// \\brief Calculate the change in point correlations due to "
        "changing an occupant\n"
     << indent << "void " << class_name
     << "::_calc_delta_point_corr(int nlist_ind, int occ_i, int occ_f) const "
        "{\n"
     << indent << "  m_params.pre_eval();\n";

  ispec = 0;
  for (auto const &specialization : specializations) {
    if (specializations.size() > 1) {
      ss << indent << "  ";
      if (ispec > 0) ss << "else ";
      ss << "if(m_params.eval_mode() == " << specialization.first << ")";
    }
    ss << indent << "  {\n"
       << indent << "    _point_prepare<" << specialization.second
       << ">(nlist_ind);\n"
       << indent << "   for(size_type i = 0; i < corr_size(); i++) {\n"
       << indent << "      ParamPack::Val<" << specialization.second
       << ">::set(m_params, m_corr_param_key, i, (this->*m_delta_func_table_"
       << ispec << "[nlist_ind][i])(occ_i, occ_f));\n"
       << indent << "    }\n"
       << indent << "  }\n";
    ++ispec;
  }
  ss << indent << "  m_params.post_eval();\n"
     << indent << "}\n\n"
     <<

      //-----

      indent
     << "/// \\brief Calculate the change in select point correlations due to "
        "changing an occupant\n"
     << indent << "void " << class_name
     << "::_calc_restricted_delta_point_corr(int nlist_ind, int occ_i, int "
        "occ_f, double *corr_begin, size_type const *ind_list_begin, size_type "
        "const *ind_list_end) const {\n"
     << indent
     << "  _calc_restricted_delta_point_corr(nlist_ind, occ_i, occ_f, "
        "ind_list_begin, ind_list_end);\n"
     << indent << "  for(; ind_list_begin < ind_list_end; ind_list_begin++) {\n"
     << indent
     << "    *(corr_begin + *ind_list_begin) = "
        "ParamPack::Val<double>::get(m_params, m_corr_param_key, "
        "*ind_list_begin);\n"
     << indent << "  }\n"
     << indent << "}\n\n"
     <<

      //-----

      indent
     << "/// \\brief Calculate the change in select point correlations due to "
        "changing an occupant\n"
     << indent << "void " << class_name
     << "::_calc_restricted_delta_point_corr(int nlist_ind, int occ_i, int "
        "occ_f, size_type const *ind_list_begin, size_type const "
        "*ind_list_end) const {\n"
     << indent << "  m_params.pre_eval();\n";

  ispec = 0;
  for (auto const &specialization : specializations) {
    if (specializations.size() > 1) {
      ss << indent << "  ";
      if (ispec > 0) ss << "else ";
      ss << "if(m_params.eval_mode() == " << specialization.first << ")";
    }
    ss << indent << "  {\n"
       << indent << "    _point_prepare<" << specialization.second
       << ">(nlist_ind);\n"
       << indent
       << "    for(; ind_list_begin < ind_list_end; ind_list_begin++) {\n"
       << indent << "      ParamPack::Val<" << specialization.second
       << ">::set(m_params, m_corr_param_key, *ind_list_begin, "
          "(this->*m_delta_func_table_"
       << ispec << "[nlist_ind][*ind_list_begin])(occ_i, occ_f));\n"
       << indent << "    }\n"
       << indent << "  }\n";
    ++ispec;
  }
  ss << indent << "  m_params.post_eval();\n" << indent << "}\n\n";

  return ss.str();
}

}  // namespace ClexBasisWriter_impl
}  // namespace CASM
