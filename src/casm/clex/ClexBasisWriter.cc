#include "casm/clex/ClexBasisWriter.hh"
#include "casm/clex/ClexBasis.hh"
#include "casm/clex/OrbitFunctionTraits.hh"
#include "casm/basis_set/DoFTraits.hh"
#include "casm/basis_set/FunctionVisitor.hh"
namespace CASM {

  ClexBasisWriter::ClexBasisWriter(Structure const &_prim, ParamPackMixIn const *parampack_mix_in /*= nullptr*/) {
    //throw std::runtime_error("Error: print_clexulator is being re-implemented");

  }

  namespace ClexBasisWriter_impl {

    std::string clexulator_member_declarations(std::string const &class_name,
                                               ClexBasis const &clex,
                                               std::vector<std::unique_ptr<OrbitFunctionTraits> > const &orbit_func_writers,
                                               std::string const &indent) {
      Index N_corr = clex.n_functions();
      Index N_sublat = clex.n_sublat();
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
        writer_ptr->print_eval_table_declarations(ss, class_name, clex, indent);
      }


      ss <<
         indent << "// typedef for method pointers\n" <<
         indent << "typedef double (" << class_name << "::*BasisFuncPtr)() const;\n\n" <<

         indent << "// typedef for method pointers\n" <<
         indent << "typedef double (" << class_name << "::*DeltaBasisFuncPtr)(int, int) const;\n\n" <<

         indent << "// array of pointers to member functions for calculating basis functions\n" <<
         indent << "BasisFuncPtr m_orbit_func_table[" << N_corr << "];\n\n" <<

         indent << "// array of pointers to member functions for calculating flower functions\n" <<
         indent << "BasisFuncPtr m_flower_func_table[" << N_sublat << "][" << N_corr << "];\n\n" <<

         indent << "// array of pointers to member functions for calculating DELTA flower functions\n" <<
         indent << "DeltaBasisFuncPtr m_delta_func_table[" << N_sublat << "][" << N_corr << "];\n\n";

      auto it(clex.site_bases().begin()), end_it(clex.site_bases().end());
      for(; it != end_it; ++it) {
        ss << DoFType::traits(it->first).clexulator_member_declarations_string(clex.prim(), it->second, indent);
      }
      return ss.str();
    }

    //*******************************************************************************************

    std::string clexulator_private_method_declarations(std::string const &class_name,
                                                       ClexBasis const &clex,
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

         indent << "  /// \\brief Calculate point correlations about neighbor site 'n_index'\n" <<
         indent << "  /// For global clexulators, 'n_index' only ranges over sites in the cell\n" <<
         indent << "  /// For local clexulators, 'n_index' ranges over all sites in the neighborhood\n" <<
         indent << "  void _calc_point_corr(int n_index, double *corr_begin) const override;\n\n" <<

         indent << "  /// \\brief Calculate select point correlations about neighbor site 'n_index'\n" <<
         indent << "  /// For global clexulators, 'n_index' only ranges over sites in the cell\n" <<
         indent << "  /// For local clexulators, 'n_index' ranges over all sites in the neighborhood\n" <<
         indent << "  void _calc_restricted_point_corr(int n_index, double *corr_begin, size_type const* ind_list_begin, size_type const* ind_list_end) const override;\n\n" <<

         indent << "  /// \\brief Calculate the change in point correlations due to changing an occupant at neighbor site 'n_index'\n" <<
         indent << "  /// For global clexulators, 'n_index' only ranges over sites in the cell\n" <<
         indent << "  /// For local clexulators, 'n_index' ranges over all sites in the neighborhood\n" <<
         indent << "  void _calc_delta_point_corr(int n_index, int occ_i, int occ_f, double *corr_begin) const override;\n\n" <<

         indent << "  /// \\brief Calculate the change in select point correlations due to changing an occupant at neighbor site 'n_index'\n" <<
         indent << "  /// For global clexulators, 'n_index' only ranges over sites in the cell\n" <<
         indent << "  /// For local clexulators, 'n_index' ranges over all sites in the neighborhood\n" <<
         indent << "  void _calc_restricted_delta_point_corr(int b_index, int occ_i, int occ_f, double *corr_begin, size_type const* ind_list_begin, size_type const* ind_list_end) const override;\n\n" <<

         indent << "  void _global_prepare() const override;\n\n" <<

         indent << "  void _point_prepare(int neighbor_ind) const override;\n";



      auto it(clex.site_bases().begin()), end_it(clex.site_bases().end());
      for(; it != end_it; ++it) {
        ss << DoFType::traits(it->first).clexulator_private_method_declarations_string(clex.prim(), it->second, indent);
      }

      ss <<
         indent << "  //default functions for basis function evaluation \n" <<
         indent << "  double zero_func() const{ return 0.0;};\n" <<
         indent << "  double zero_func(int,int) const{ return 0.0;};\n\n";

      return ss.str();
    }

    //*******************************************************************************************

    std::string clexulator_public_method_declarations(std::string const &class_name,
                                                      ClexBasis const &clex,
                                                      std::string const &indent) {
      std::stringstream ss;
      ss <<
         indent << "  " << class_name << "();\n\n" <<
         indent << "  ~" << class_name << "();\n\n";

      auto it(clex.site_bases().begin()), end_it(clex.site_bases().end());
      for(; it != end_it; ++it) {
        ss << DoFType::traits(it->first).clexulator_public_method_declarations_string(clex.prim(), it->second, indent);
      }
      return ss.str();
    }


    //*******************************************************************************************

    std::string clexulator_interface_implementation(std::string const &class_name,
                                                    ClexBasis const &clex,
                                                    std::string const &indent) {
      std::stringstream ss;
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


      return ss.str();
    }

  }//\namespace ClexBasisWriter_impl
}//\namespace CASM
