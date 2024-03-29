#include <cstddef>

#include "casm/clex/Clexulator.hh"

/****** CLEXULATOR CLASS FOR PRIM ******
{
  "basis" : [
    {
      "coordinate" : [ 0.000000000000, 0.000000000000, 0.000000000000 ],
      "occupant_dof" : [ "A", "B", "C" ]
    }
  ],
  "coordinate_mode" : "Fractional",
  "lattice_vectors" : [
    [ 2.000000000000, 2.000000000000, 0.000000000000 ],
    [ 0.000000000000, 2.000000000000, 2.000000000000 ],
    [ 2.000000000000, 0.000000000000, 2.000000000000 ]
  ],
  "title" : "Test"
}**/

/// \brief Returns a clexulator::BaseClexulator* owning a test_Clexulator
extern "C" CASM::clexulator::BaseClexulator *make_test_Clexulator();

namespace CASM {

class test_Clexulator : public clexulator::BaseClexulator {
 public:
  test_Clexulator();

  ~test_Clexulator();

  /// \brief Clone the test_Clexulator
  std::unique_ptr<test_Clexulator> clone() const {
    return std::unique_ptr<test_Clexulator>(_clone());
  }

  /// \brief Calculate contribution to global correlations from one unit cell
  void calc_global_corr_contribution(double *corr_begin) const override;

  /// \brief Calculate contribution to select global correlations from one unit
  /// cell
  void calc_restricted_global_corr_contribution(
      double *corr_begin, size_type const *ind_list_begin,
      size_type const *ind_list_end) const override;

  /// \brief Calculate point correlations about basis site 'b_index'
  void calc_point_corr(int b_index, double *corr_begin) const override;

  /// \brief Calculate select point correlations about basis site 'b_index'
  void calc_restricted_point_corr(int b_index, double *corr_begin,
                                  size_type const *ind_list_begin,
                                  size_type const *ind_list_end) const override;

  /// \brief Calculate the change in point correlations due to changing an
  /// occupant
  void calc_delta_point_corr(int b_index, int occ_i, int occ_f,
                             double *corr_begin) const override;

  /// \brief Calculate the change in select point correlations due to changing
  /// an occupant
  void calc_restricted_delta_point_corr(
      int b_index, int occ_i, int occ_f, double *corr_begin,
      size_type const *ind_list_begin,
      size_type const *ind_list_end) const override;

 private:
  /// \brief Clone the Clexulator
  virtual test_Clexulator *_clone() const override {
    return new test_Clexulator(*this);
  }

  // typedef for method pointers
  typedef double (test_Clexulator::*BasisFuncPtr)() const;

  // typedef for method pointers
  typedef double (test_Clexulator::*DeltaBasisFuncPtr)(int, int) const;

  // array of pointers to member functions for calculating basis functions
  BasisFuncPtr m_orbit_func_list[75];

  // array of pointers to member functions for calculating flower functions
  BasisFuncPtr m_flower_func_lists[1][75];

  // array of pointers to member functions for calculating DELTA flower
  // functions
  DeltaBasisFuncPtr m_delta_func_lists[1][75];

  // Occupation Function tables for basis sites in asymmetric unit 0:
  //   - basis site 0:
  double m_occ_func_0_0[3];
  double m_occ_func_0_1[3];

  // Occupation Function accessors for basis site 0:
  const double &occ_func_0_0(const int &nlist_ind) const {
    return m_occ_func_0_0[*(m_occ_ptr + *(m_nlist_ptr + nlist_ind))];
  }
  const double &occ_func_0_1(const int &nlist_ind) const {
    return m_occ_func_0_1[*(m_occ_ptr + *(m_nlist_ptr + nlist_ind))];
  }

  // default functions for basis function evaluation
  double zero_func() const { return 0.0; };
  double zero_func(int, int) const { return 0.0; };

  double eval_bfunc_0_0_0() const;

  double eval_bfunc_1_0_0() const;
  double eval_bfunc_1_0_1() const;

  double site_eval_at_0_bfunc_1_0_0() const;
  double site_eval_at_0_bfunc_1_0_1() const;

  double delta_site_eval_at_0_bfunc_1_0_0(int occ_i, int occ_f) const;
  double delta_site_eval_at_0_bfunc_1_0_1(int occ_i, int occ_f) const;

  double eval_bfunc_2_0_0() const;
  double eval_bfunc_2_0_1() const;
  double eval_bfunc_2_0_2() const;

  double site_eval_at_0_bfunc_2_0_0() const;
  double site_eval_at_0_bfunc_2_0_1() const;
  double site_eval_at_0_bfunc_2_0_2() const;

  double delta_site_eval_at_0_bfunc_2_0_0(int occ_i, int occ_f) const;
  double delta_site_eval_at_0_bfunc_2_0_1(int occ_i, int occ_f) const;
  double delta_site_eval_at_0_bfunc_2_0_2(int occ_i, int occ_f) const;

  double eval_bfunc_2_1_0() const;
  double eval_bfunc_2_1_1() const;
  double eval_bfunc_2_1_2() const;

  double site_eval_at_0_bfunc_2_1_0() const;
  double site_eval_at_0_bfunc_2_1_1() const;
  double site_eval_at_0_bfunc_2_1_2() const;

  double delta_site_eval_at_0_bfunc_2_1_0(int occ_i, int occ_f) const;
  double delta_site_eval_at_0_bfunc_2_1_1(int occ_i, int occ_f) const;
  double delta_site_eval_at_0_bfunc_2_1_2(int occ_i, int occ_f) const;

  double eval_bfunc_2_2_0() const;
  double eval_bfunc_2_2_1() const;
  double eval_bfunc_2_2_2() const;

  double site_eval_at_0_bfunc_2_2_0() const;
  double site_eval_at_0_bfunc_2_2_1() const;
  double site_eval_at_0_bfunc_2_2_2() const;

  double delta_site_eval_at_0_bfunc_2_2_0(int occ_i, int occ_f) const;
  double delta_site_eval_at_0_bfunc_2_2_1(int occ_i, int occ_f) const;
  double delta_site_eval_at_0_bfunc_2_2_2(int occ_i, int occ_f) const;

  double eval_bfunc_2_3_0() const;
  double eval_bfunc_2_3_1() const;
  double eval_bfunc_2_3_2() const;

  double site_eval_at_0_bfunc_2_3_0() const;
  double site_eval_at_0_bfunc_2_3_1() const;
  double site_eval_at_0_bfunc_2_3_2() const;

  double delta_site_eval_at_0_bfunc_2_3_0(int occ_i, int occ_f) const;
  double delta_site_eval_at_0_bfunc_2_3_1(int occ_i, int occ_f) const;
  double delta_site_eval_at_0_bfunc_2_3_2(int occ_i, int occ_f) const;

  double eval_bfunc_2_4_0() const;
  double eval_bfunc_2_4_1() const;
  double eval_bfunc_2_4_2() const;

  double site_eval_at_0_bfunc_2_4_0() const;
  double site_eval_at_0_bfunc_2_4_1() const;
  double site_eval_at_0_bfunc_2_4_2() const;

  double delta_site_eval_at_0_bfunc_2_4_0(int occ_i, int occ_f) const;
  double delta_site_eval_at_0_bfunc_2_4_1(int occ_i, int occ_f) const;
  double delta_site_eval_at_0_bfunc_2_4_2(int occ_i, int occ_f) const;

  double eval_bfunc_2_5_0() const;
  double eval_bfunc_2_5_1() const;
  double eval_bfunc_2_5_2() const;

  double site_eval_at_0_bfunc_2_5_0() const;
  double site_eval_at_0_bfunc_2_5_1() const;
  double site_eval_at_0_bfunc_2_5_2() const;

  double delta_site_eval_at_0_bfunc_2_5_0(int occ_i, int occ_f) const;
  double delta_site_eval_at_0_bfunc_2_5_1(int occ_i, int occ_f) const;
  double delta_site_eval_at_0_bfunc_2_5_2(int occ_i, int occ_f) const;

  double eval_bfunc_3_0_0() const;
  double eval_bfunc_3_0_1() const;
  double eval_bfunc_3_0_2() const;
  double eval_bfunc_3_0_3() const;

  double site_eval_at_0_bfunc_3_0_0() const;
  double site_eval_at_0_bfunc_3_0_1() const;
  double site_eval_at_0_bfunc_3_0_2() const;
  double site_eval_at_0_bfunc_3_0_3() const;

  double delta_site_eval_at_0_bfunc_3_0_0(int occ_i, int occ_f) const;
  double delta_site_eval_at_0_bfunc_3_0_1(int occ_i, int occ_f) const;
  double delta_site_eval_at_0_bfunc_3_0_2(int occ_i, int occ_f) const;
  double delta_site_eval_at_0_bfunc_3_0_3(int occ_i, int occ_f) const;

  double eval_bfunc_3_1_0() const;
  double eval_bfunc_3_1_1() const;
  double eval_bfunc_3_1_2() const;
  double eval_bfunc_3_1_3() const;
  double eval_bfunc_3_1_4() const;
  double eval_bfunc_3_1_5() const;

  double site_eval_at_0_bfunc_3_1_0() const;
  double site_eval_at_0_bfunc_3_1_1() const;
  double site_eval_at_0_bfunc_3_1_2() const;
  double site_eval_at_0_bfunc_3_1_3() const;
  double site_eval_at_0_bfunc_3_1_4() const;
  double site_eval_at_0_bfunc_3_1_5() const;

  double delta_site_eval_at_0_bfunc_3_1_0(int occ_i, int occ_f) const;
  double delta_site_eval_at_0_bfunc_3_1_1(int occ_i, int occ_f) const;
  double delta_site_eval_at_0_bfunc_3_1_2(int occ_i, int occ_f) const;
  double delta_site_eval_at_0_bfunc_3_1_3(int occ_i, int occ_f) const;
  double delta_site_eval_at_0_bfunc_3_1_4(int occ_i, int occ_f) const;
  double delta_site_eval_at_0_bfunc_3_1_5(int occ_i, int occ_f) const;

  double eval_bfunc_3_2_0() const;
  double eval_bfunc_3_2_1() const;
  double eval_bfunc_3_2_2() const;
  double eval_bfunc_3_2_3() const;
  double eval_bfunc_3_2_4() const;
  double eval_bfunc_3_2_5() const;

  double site_eval_at_0_bfunc_3_2_0() const;
  double site_eval_at_0_bfunc_3_2_1() const;
  double site_eval_at_0_bfunc_3_2_2() const;
  double site_eval_at_0_bfunc_3_2_3() const;
  double site_eval_at_0_bfunc_3_2_4() const;
  double site_eval_at_0_bfunc_3_2_5() const;

  double delta_site_eval_at_0_bfunc_3_2_0(int occ_i, int occ_f) const;
  double delta_site_eval_at_0_bfunc_3_2_1(int occ_i, int occ_f) const;
  double delta_site_eval_at_0_bfunc_3_2_2(int occ_i, int occ_f) const;
  double delta_site_eval_at_0_bfunc_3_2_3(int occ_i, int occ_f) const;
  double delta_site_eval_at_0_bfunc_3_2_4(int occ_i, int occ_f) const;
  double delta_site_eval_at_0_bfunc_3_2_5(int occ_i, int occ_f) const;

  double eval_bfunc_3_3_0() const;
  double eval_bfunc_3_3_1() const;
  double eval_bfunc_3_3_2() const;
  double eval_bfunc_3_3_3() const;
  double eval_bfunc_3_3_4() const;
  double eval_bfunc_3_3_5() const;
  double eval_bfunc_3_3_6() const;
  double eval_bfunc_3_3_7() const;

  double site_eval_at_0_bfunc_3_3_0() const;
  double site_eval_at_0_bfunc_3_3_1() const;
  double site_eval_at_0_bfunc_3_3_2() const;
  double site_eval_at_0_bfunc_3_3_3() const;
  double site_eval_at_0_bfunc_3_3_4() const;
  double site_eval_at_0_bfunc_3_3_5() const;
  double site_eval_at_0_bfunc_3_3_6() const;
  double site_eval_at_0_bfunc_3_3_7() const;

  double delta_site_eval_at_0_bfunc_3_3_0(int occ_i, int occ_f) const;
  double delta_site_eval_at_0_bfunc_3_3_1(int occ_i, int occ_f) const;
  double delta_site_eval_at_0_bfunc_3_3_2(int occ_i, int occ_f) const;
  double delta_site_eval_at_0_bfunc_3_3_3(int occ_i, int occ_f) const;
  double delta_site_eval_at_0_bfunc_3_3_4(int occ_i, int occ_f) const;
  double delta_site_eval_at_0_bfunc_3_3_5(int occ_i, int occ_f) const;
  double delta_site_eval_at_0_bfunc_3_3_6(int occ_i, int occ_f) const;
  double delta_site_eval_at_0_bfunc_3_3_7(int occ_i, int occ_f) const;

  double eval_bfunc_3_4_0() const;
  double eval_bfunc_3_4_1() const;
  double eval_bfunc_3_4_2() const;
  double eval_bfunc_3_4_3() const;
  double eval_bfunc_3_4_4() const;
  double eval_bfunc_3_4_5() const;
  double eval_bfunc_3_4_6() const;
  double eval_bfunc_3_4_7() const;

  double site_eval_at_0_bfunc_3_4_0() const;
  double site_eval_at_0_bfunc_3_4_1() const;
  double site_eval_at_0_bfunc_3_4_2() const;
  double site_eval_at_0_bfunc_3_4_3() const;
  double site_eval_at_0_bfunc_3_4_4() const;
  double site_eval_at_0_bfunc_3_4_5() const;
  double site_eval_at_0_bfunc_3_4_6() const;
  double site_eval_at_0_bfunc_3_4_7() const;

  double delta_site_eval_at_0_bfunc_3_4_0(int occ_i, int occ_f) const;
  double delta_site_eval_at_0_bfunc_3_4_1(int occ_i, int occ_f) const;
  double delta_site_eval_at_0_bfunc_3_4_2(int occ_i, int occ_f) const;
  double delta_site_eval_at_0_bfunc_3_4_3(int occ_i, int occ_f) const;
  double delta_site_eval_at_0_bfunc_3_4_4(int occ_i, int occ_f) const;
  double delta_site_eval_at_0_bfunc_3_4_5(int occ_i, int occ_f) const;
  double delta_site_eval_at_0_bfunc_3_4_6(int occ_i, int occ_f) const;
  double delta_site_eval_at_0_bfunc_3_4_7(int occ_i, int occ_f) const;

  double eval_bfunc_4_0_0() const;
  double eval_bfunc_4_0_1() const;
  double eval_bfunc_4_0_2() const;
  double eval_bfunc_4_0_3() const;
  double eval_bfunc_4_0_4() const;
  double eval_bfunc_4_0_5() const;
  double eval_bfunc_4_0_6() const;
  double eval_bfunc_4_0_7() const;
  double eval_bfunc_4_0_8() const;
  double eval_bfunc_4_0_9() const;
  double eval_bfunc_4_0_10() const;
  double eval_bfunc_4_0_11() const;

  double site_eval_at_0_bfunc_4_0_0() const;
  double site_eval_at_0_bfunc_4_0_1() const;
  double site_eval_at_0_bfunc_4_0_2() const;
  double site_eval_at_0_bfunc_4_0_3() const;
  double site_eval_at_0_bfunc_4_0_4() const;
  double site_eval_at_0_bfunc_4_0_5() const;
  double site_eval_at_0_bfunc_4_0_6() const;
  double site_eval_at_0_bfunc_4_0_7() const;
  double site_eval_at_0_bfunc_4_0_8() const;
  double site_eval_at_0_bfunc_4_0_9() const;
  double site_eval_at_0_bfunc_4_0_10() const;
  double site_eval_at_0_bfunc_4_0_11() const;

  double delta_site_eval_at_0_bfunc_4_0_0(int occ_i, int occ_f) const;
  double delta_site_eval_at_0_bfunc_4_0_1(int occ_i, int occ_f) const;
  double delta_site_eval_at_0_bfunc_4_0_2(int occ_i, int occ_f) const;
  double delta_site_eval_at_0_bfunc_4_0_3(int occ_i, int occ_f) const;
  double delta_site_eval_at_0_bfunc_4_0_4(int occ_i, int occ_f) const;
  double delta_site_eval_at_0_bfunc_4_0_5(int occ_i, int occ_f) const;
  double delta_site_eval_at_0_bfunc_4_0_6(int occ_i, int occ_f) const;
  double delta_site_eval_at_0_bfunc_4_0_7(int occ_i, int occ_f) const;
  double delta_site_eval_at_0_bfunc_4_0_8(int occ_i, int occ_f) const;
  double delta_site_eval_at_0_bfunc_4_0_9(int occ_i, int occ_f) const;
  double delta_site_eval_at_0_bfunc_4_0_10(int occ_i, int occ_f) const;
  double delta_site_eval_at_0_bfunc_4_0_11(int occ_i, int occ_f) const;

  double eval_bfunc_4_1_0() const;
  double eval_bfunc_4_1_1() const;
  double eval_bfunc_4_1_2() const;
  double eval_bfunc_4_1_3() const;
  double eval_bfunc_4_1_4() const;
  double eval_bfunc_4_1_5() const;
  double eval_bfunc_4_1_6() const;
  double eval_bfunc_4_1_7() const;
  double eval_bfunc_4_1_8() const;
  double eval_bfunc_4_1_9() const;

  double site_eval_at_0_bfunc_4_1_0() const;
  double site_eval_at_0_bfunc_4_1_1() const;
  double site_eval_at_0_bfunc_4_1_2() const;
  double site_eval_at_0_bfunc_4_1_3() const;
  double site_eval_at_0_bfunc_4_1_4() const;
  double site_eval_at_0_bfunc_4_1_5() const;
  double site_eval_at_0_bfunc_4_1_6() const;
  double site_eval_at_0_bfunc_4_1_7() const;
  double site_eval_at_0_bfunc_4_1_8() const;
  double site_eval_at_0_bfunc_4_1_9() const;

  double delta_site_eval_at_0_bfunc_4_1_0(int occ_i, int occ_f) const;
  double delta_site_eval_at_0_bfunc_4_1_1(int occ_i, int occ_f) const;
  double delta_site_eval_at_0_bfunc_4_1_2(int occ_i, int occ_f) const;
  double delta_site_eval_at_0_bfunc_4_1_3(int occ_i, int occ_f) const;
  double delta_site_eval_at_0_bfunc_4_1_4(int occ_i, int occ_f) const;
  double delta_site_eval_at_0_bfunc_4_1_5(int occ_i, int occ_f) const;
  double delta_site_eval_at_0_bfunc_4_1_6(int occ_i, int occ_f) const;
  double delta_site_eval_at_0_bfunc_4_1_7(int occ_i, int occ_f) const;
  double delta_site_eval_at_0_bfunc_4_1_8(int occ_i, int occ_f) const;
  double delta_site_eval_at_0_bfunc_4_1_9(int occ_i, int occ_f) const;
};

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

test_Clexulator::test_Clexulator() : clexulator::BaseClexulator(177, 75) {
  m_occ_func_0_0[0] = 0.0000000000, m_occ_func_0_0[1] = 1.0000000000,
  m_occ_func_0_0[2] = 0.0000000000;

  m_occ_func_0_1[0] = 0.0000000000, m_occ_func_0_1[1] = 0.0000000000,
  m_occ_func_0_1[2] = 1.0000000000;

  m_orbit_func_list[0] = &test_Clexulator::eval_bfunc_0_0_0;
  m_orbit_func_list[1] = &test_Clexulator::eval_bfunc_1_0_0;
  m_orbit_func_list[2] = &test_Clexulator::eval_bfunc_1_0_1;
  m_orbit_func_list[3] = &test_Clexulator::eval_bfunc_2_0_0;
  m_orbit_func_list[4] = &test_Clexulator::eval_bfunc_2_0_1;
  m_orbit_func_list[5] = &test_Clexulator::eval_bfunc_2_0_2;
  m_orbit_func_list[6] = &test_Clexulator::eval_bfunc_2_1_0;
  m_orbit_func_list[7] = &test_Clexulator::eval_bfunc_2_1_1;
  m_orbit_func_list[8] = &test_Clexulator::eval_bfunc_2_1_2;
  m_orbit_func_list[9] = &test_Clexulator::eval_bfunc_2_2_0;
  m_orbit_func_list[10] = &test_Clexulator::eval_bfunc_2_2_1;
  m_orbit_func_list[11] = &test_Clexulator::eval_bfunc_2_2_2;
  m_orbit_func_list[12] = &test_Clexulator::eval_bfunc_2_3_0;
  m_orbit_func_list[13] = &test_Clexulator::eval_bfunc_2_3_1;
  m_orbit_func_list[14] = &test_Clexulator::eval_bfunc_2_3_2;
  m_orbit_func_list[15] = &test_Clexulator::eval_bfunc_2_4_0;
  m_orbit_func_list[16] = &test_Clexulator::eval_bfunc_2_4_1;
  m_orbit_func_list[17] = &test_Clexulator::eval_bfunc_2_4_2;
  m_orbit_func_list[18] = &test_Clexulator::eval_bfunc_2_5_0;
  m_orbit_func_list[19] = &test_Clexulator::eval_bfunc_2_5_1;
  m_orbit_func_list[20] = &test_Clexulator::eval_bfunc_2_5_2;
  m_orbit_func_list[21] = &test_Clexulator::eval_bfunc_3_0_0;
  m_orbit_func_list[22] = &test_Clexulator::eval_bfunc_3_0_1;
  m_orbit_func_list[23] = &test_Clexulator::eval_bfunc_3_0_2;
  m_orbit_func_list[24] = &test_Clexulator::eval_bfunc_3_0_3;
  m_orbit_func_list[25] = &test_Clexulator::eval_bfunc_3_1_0;
  m_orbit_func_list[26] = &test_Clexulator::eval_bfunc_3_1_1;
  m_orbit_func_list[27] = &test_Clexulator::eval_bfunc_3_1_2;
  m_orbit_func_list[28] = &test_Clexulator::eval_bfunc_3_1_3;
  m_orbit_func_list[29] = &test_Clexulator::eval_bfunc_3_1_4;
  m_orbit_func_list[30] = &test_Clexulator::eval_bfunc_3_1_5;
  m_orbit_func_list[31] = &test_Clexulator::eval_bfunc_3_2_0;
  m_orbit_func_list[32] = &test_Clexulator::eval_bfunc_3_2_1;
  m_orbit_func_list[33] = &test_Clexulator::eval_bfunc_3_2_2;
  m_orbit_func_list[34] = &test_Clexulator::eval_bfunc_3_2_3;
  m_orbit_func_list[35] = &test_Clexulator::eval_bfunc_3_2_4;
  m_orbit_func_list[36] = &test_Clexulator::eval_bfunc_3_2_5;
  m_orbit_func_list[37] = &test_Clexulator::eval_bfunc_3_3_0;
  m_orbit_func_list[38] = &test_Clexulator::eval_bfunc_3_3_1;
  m_orbit_func_list[39] = &test_Clexulator::eval_bfunc_3_3_2;
  m_orbit_func_list[40] = &test_Clexulator::eval_bfunc_3_3_3;
  m_orbit_func_list[41] = &test_Clexulator::eval_bfunc_3_3_4;
  m_orbit_func_list[42] = &test_Clexulator::eval_bfunc_3_3_5;
  m_orbit_func_list[43] = &test_Clexulator::eval_bfunc_3_3_6;
  m_orbit_func_list[44] = &test_Clexulator::eval_bfunc_3_3_7;
  m_orbit_func_list[45] = &test_Clexulator::eval_bfunc_3_4_0;
  m_orbit_func_list[46] = &test_Clexulator::eval_bfunc_3_4_1;
  m_orbit_func_list[47] = &test_Clexulator::eval_bfunc_3_4_2;
  m_orbit_func_list[48] = &test_Clexulator::eval_bfunc_3_4_3;
  m_orbit_func_list[49] = &test_Clexulator::eval_bfunc_3_4_4;
  m_orbit_func_list[50] = &test_Clexulator::eval_bfunc_3_4_5;
  m_orbit_func_list[51] = &test_Clexulator::eval_bfunc_3_4_6;
  m_orbit_func_list[52] = &test_Clexulator::eval_bfunc_3_4_7;
  m_orbit_func_list[53] = &test_Clexulator::eval_bfunc_4_0_0;
  m_orbit_func_list[54] = &test_Clexulator::eval_bfunc_4_0_1;
  m_orbit_func_list[55] = &test_Clexulator::eval_bfunc_4_0_2;
  m_orbit_func_list[56] = &test_Clexulator::eval_bfunc_4_0_3;
  m_orbit_func_list[57] = &test_Clexulator::eval_bfunc_4_0_4;
  m_orbit_func_list[58] = &test_Clexulator::eval_bfunc_4_0_5;
  m_orbit_func_list[59] = &test_Clexulator::eval_bfunc_4_0_6;
  m_orbit_func_list[60] = &test_Clexulator::eval_bfunc_4_0_7;
  m_orbit_func_list[61] = &test_Clexulator::eval_bfunc_4_0_8;
  m_orbit_func_list[62] = &test_Clexulator::eval_bfunc_4_0_9;
  m_orbit_func_list[63] = &test_Clexulator::eval_bfunc_4_0_10;
  m_orbit_func_list[64] = &test_Clexulator::eval_bfunc_4_0_11;
  m_orbit_func_list[65] = &test_Clexulator::eval_bfunc_4_1_0;
  m_orbit_func_list[66] = &test_Clexulator::eval_bfunc_4_1_1;
  m_orbit_func_list[67] = &test_Clexulator::eval_bfunc_4_1_2;
  m_orbit_func_list[68] = &test_Clexulator::eval_bfunc_4_1_3;
  m_orbit_func_list[69] = &test_Clexulator::eval_bfunc_4_1_4;
  m_orbit_func_list[70] = &test_Clexulator::eval_bfunc_4_1_5;
  m_orbit_func_list[71] = &test_Clexulator::eval_bfunc_4_1_6;
  m_orbit_func_list[72] = &test_Clexulator::eval_bfunc_4_1_7;
  m_orbit_func_list[73] = &test_Clexulator::eval_bfunc_4_1_8;
  m_orbit_func_list[74] = &test_Clexulator::eval_bfunc_4_1_9;

  m_flower_func_lists[0][0] = &test_Clexulator::zero_func;
  m_flower_func_lists[0][1] = &test_Clexulator::site_eval_at_0_bfunc_1_0_0;
  m_flower_func_lists[0][2] = &test_Clexulator::site_eval_at_0_bfunc_1_0_1;
  m_flower_func_lists[0][3] = &test_Clexulator::site_eval_at_0_bfunc_2_0_0;
  m_flower_func_lists[0][4] = &test_Clexulator::site_eval_at_0_bfunc_2_0_1;
  m_flower_func_lists[0][5] = &test_Clexulator::site_eval_at_0_bfunc_2_0_2;
  m_flower_func_lists[0][6] = &test_Clexulator::site_eval_at_0_bfunc_2_1_0;
  m_flower_func_lists[0][7] = &test_Clexulator::site_eval_at_0_bfunc_2_1_1;
  m_flower_func_lists[0][8] = &test_Clexulator::site_eval_at_0_bfunc_2_1_2;
  m_flower_func_lists[0][9] = &test_Clexulator::site_eval_at_0_bfunc_2_2_0;
  m_flower_func_lists[0][10] = &test_Clexulator::site_eval_at_0_bfunc_2_2_1;
  m_flower_func_lists[0][11] = &test_Clexulator::site_eval_at_0_bfunc_2_2_2;
  m_flower_func_lists[0][12] = &test_Clexulator::site_eval_at_0_bfunc_2_3_0;
  m_flower_func_lists[0][13] = &test_Clexulator::site_eval_at_0_bfunc_2_3_1;
  m_flower_func_lists[0][14] = &test_Clexulator::site_eval_at_0_bfunc_2_3_2;
  m_flower_func_lists[0][15] = &test_Clexulator::site_eval_at_0_bfunc_2_4_0;
  m_flower_func_lists[0][16] = &test_Clexulator::site_eval_at_0_bfunc_2_4_1;
  m_flower_func_lists[0][17] = &test_Clexulator::site_eval_at_0_bfunc_2_4_2;
  m_flower_func_lists[0][18] = &test_Clexulator::site_eval_at_0_bfunc_2_5_0;
  m_flower_func_lists[0][19] = &test_Clexulator::site_eval_at_0_bfunc_2_5_1;
  m_flower_func_lists[0][20] = &test_Clexulator::site_eval_at_0_bfunc_2_5_2;
  m_flower_func_lists[0][21] = &test_Clexulator::site_eval_at_0_bfunc_3_0_0;
  m_flower_func_lists[0][22] = &test_Clexulator::site_eval_at_0_bfunc_3_0_1;
  m_flower_func_lists[0][23] = &test_Clexulator::site_eval_at_0_bfunc_3_0_2;
  m_flower_func_lists[0][24] = &test_Clexulator::site_eval_at_0_bfunc_3_0_3;
  m_flower_func_lists[0][25] = &test_Clexulator::site_eval_at_0_bfunc_3_1_0;
  m_flower_func_lists[0][26] = &test_Clexulator::site_eval_at_0_bfunc_3_1_1;
  m_flower_func_lists[0][27] = &test_Clexulator::site_eval_at_0_bfunc_3_1_2;
  m_flower_func_lists[0][28] = &test_Clexulator::site_eval_at_0_bfunc_3_1_3;
  m_flower_func_lists[0][29] = &test_Clexulator::site_eval_at_0_bfunc_3_1_4;
  m_flower_func_lists[0][30] = &test_Clexulator::site_eval_at_0_bfunc_3_1_5;
  m_flower_func_lists[0][31] = &test_Clexulator::site_eval_at_0_bfunc_3_2_0;
  m_flower_func_lists[0][32] = &test_Clexulator::site_eval_at_0_bfunc_3_2_1;
  m_flower_func_lists[0][33] = &test_Clexulator::site_eval_at_0_bfunc_3_2_2;
  m_flower_func_lists[0][34] = &test_Clexulator::site_eval_at_0_bfunc_3_2_3;
  m_flower_func_lists[0][35] = &test_Clexulator::site_eval_at_0_bfunc_3_2_4;
  m_flower_func_lists[0][36] = &test_Clexulator::site_eval_at_0_bfunc_3_2_5;
  m_flower_func_lists[0][37] = &test_Clexulator::site_eval_at_0_bfunc_3_3_0;
  m_flower_func_lists[0][38] = &test_Clexulator::site_eval_at_0_bfunc_3_3_1;
  m_flower_func_lists[0][39] = &test_Clexulator::site_eval_at_0_bfunc_3_3_2;
  m_flower_func_lists[0][40] = &test_Clexulator::site_eval_at_0_bfunc_3_3_3;
  m_flower_func_lists[0][41] = &test_Clexulator::site_eval_at_0_bfunc_3_3_4;
  m_flower_func_lists[0][42] = &test_Clexulator::site_eval_at_0_bfunc_3_3_5;
  m_flower_func_lists[0][43] = &test_Clexulator::site_eval_at_0_bfunc_3_3_6;
  m_flower_func_lists[0][44] = &test_Clexulator::site_eval_at_0_bfunc_3_3_7;
  m_flower_func_lists[0][45] = &test_Clexulator::site_eval_at_0_bfunc_3_4_0;
  m_flower_func_lists[0][46] = &test_Clexulator::site_eval_at_0_bfunc_3_4_1;
  m_flower_func_lists[0][47] = &test_Clexulator::site_eval_at_0_bfunc_3_4_2;
  m_flower_func_lists[0][48] = &test_Clexulator::site_eval_at_0_bfunc_3_4_3;
  m_flower_func_lists[0][49] = &test_Clexulator::site_eval_at_0_bfunc_3_4_4;
  m_flower_func_lists[0][50] = &test_Clexulator::site_eval_at_0_bfunc_3_4_5;
  m_flower_func_lists[0][51] = &test_Clexulator::site_eval_at_0_bfunc_3_4_6;
  m_flower_func_lists[0][52] = &test_Clexulator::site_eval_at_0_bfunc_3_4_7;
  m_flower_func_lists[0][53] = &test_Clexulator::site_eval_at_0_bfunc_4_0_0;
  m_flower_func_lists[0][54] = &test_Clexulator::site_eval_at_0_bfunc_4_0_1;
  m_flower_func_lists[0][55] = &test_Clexulator::site_eval_at_0_bfunc_4_0_2;
  m_flower_func_lists[0][56] = &test_Clexulator::site_eval_at_0_bfunc_4_0_3;
  m_flower_func_lists[0][57] = &test_Clexulator::site_eval_at_0_bfunc_4_0_4;
  m_flower_func_lists[0][58] = &test_Clexulator::site_eval_at_0_bfunc_4_0_5;
  m_flower_func_lists[0][59] = &test_Clexulator::site_eval_at_0_bfunc_4_0_6;
  m_flower_func_lists[0][60] = &test_Clexulator::site_eval_at_0_bfunc_4_0_7;
  m_flower_func_lists[0][61] = &test_Clexulator::site_eval_at_0_bfunc_4_0_8;
  m_flower_func_lists[0][62] = &test_Clexulator::site_eval_at_0_bfunc_4_0_9;
  m_flower_func_lists[0][63] = &test_Clexulator::site_eval_at_0_bfunc_4_0_10;
  m_flower_func_lists[0][64] = &test_Clexulator::site_eval_at_0_bfunc_4_0_11;
  m_flower_func_lists[0][65] = &test_Clexulator::site_eval_at_0_bfunc_4_1_0;
  m_flower_func_lists[0][66] = &test_Clexulator::site_eval_at_0_bfunc_4_1_1;
  m_flower_func_lists[0][67] = &test_Clexulator::site_eval_at_0_bfunc_4_1_2;
  m_flower_func_lists[0][68] = &test_Clexulator::site_eval_at_0_bfunc_4_1_3;
  m_flower_func_lists[0][69] = &test_Clexulator::site_eval_at_0_bfunc_4_1_4;
  m_flower_func_lists[0][70] = &test_Clexulator::site_eval_at_0_bfunc_4_1_5;
  m_flower_func_lists[0][71] = &test_Clexulator::site_eval_at_0_bfunc_4_1_6;
  m_flower_func_lists[0][72] = &test_Clexulator::site_eval_at_0_bfunc_4_1_7;
  m_flower_func_lists[0][73] = &test_Clexulator::site_eval_at_0_bfunc_4_1_8;
  m_flower_func_lists[0][74] = &test_Clexulator::site_eval_at_0_bfunc_4_1_9;

  m_delta_func_lists[0][0] = &test_Clexulator::zero_func;
  m_delta_func_lists[0][1] = &test_Clexulator::delta_site_eval_at_0_bfunc_1_0_0;
  m_delta_func_lists[0][2] = &test_Clexulator::delta_site_eval_at_0_bfunc_1_0_1;
  m_delta_func_lists[0][3] = &test_Clexulator::delta_site_eval_at_0_bfunc_2_0_0;
  m_delta_func_lists[0][4] = &test_Clexulator::delta_site_eval_at_0_bfunc_2_0_1;
  m_delta_func_lists[0][5] = &test_Clexulator::delta_site_eval_at_0_bfunc_2_0_2;
  m_delta_func_lists[0][6] = &test_Clexulator::delta_site_eval_at_0_bfunc_2_1_0;
  m_delta_func_lists[0][7] = &test_Clexulator::delta_site_eval_at_0_bfunc_2_1_1;
  m_delta_func_lists[0][8] = &test_Clexulator::delta_site_eval_at_0_bfunc_2_1_2;
  m_delta_func_lists[0][9] = &test_Clexulator::delta_site_eval_at_0_bfunc_2_2_0;
  m_delta_func_lists[0][10] =
      &test_Clexulator::delta_site_eval_at_0_bfunc_2_2_1;
  m_delta_func_lists[0][11] =
      &test_Clexulator::delta_site_eval_at_0_bfunc_2_2_2;
  m_delta_func_lists[0][12] =
      &test_Clexulator::delta_site_eval_at_0_bfunc_2_3_0;
  m_delta_func_lists[0][13] =
      &test_Clexulator::delta_site_eval_at_0_bfunc_2_3_1;
  m_delta_func_lists[0][14] =
      &test_Clexulator::delta_site_eval_at_0_bfunc_2_3_2;
  m_delta_func_lists[0][15] =
      &test_Clexulator::delta_site_eval_at_0_bfunc_2_4_0;
  m_delta_func_lists[0][16] =
      &test_Clexulator::delta_site_eval_at_0_bfunc_2_4_1;
  m_delta_func_lists[0][17] =
      &test_Clexulator::delta_site_eval_at_0_bfunc_2_4_2;
  m_delta_func_lists[0][18] =
      &test_Clexulator::delta_site_eval_at_0_bfunc_2_5_0;
  m_delta_func_lists[0][19] =
      &test_Clexulator::delta_site_eval_at_0_bfunc_2_5_1;
  m_delta_func_lists[0][20] =
      &test_Clexulator::delta_site_eval_at_0_bfunc_2_5_2;
  m_delta_func_lists[0][21] =
      &test_Clexulator::delta_site_eval_at_0_bfunc_3_0_0;
  m_delta_func_lists[0][22] =
      &test_Clexulator::delta_site_eval_at_0_bfunc_3_0_1;
  m_delta_func_lists[0][23] =
      &test_Clexulator::delta_site_eval_at_0_bfunc_3_0_2;
  m_delta_func_lists[0][24] =
      &test_Clexulator::delta_site_eval_at_0_bfunc_3_0_3;
  m_delta_func_lists[0][25] =
      &test_Clexulator::delta_site_eval_at_0_bfunc_3_1_0;
  m_delta_func_lists[0][26] =
      &test_Clexulator::delta_site_eval_at_0_bfunc_3_1_1;
  m_delta_func_lists[0][27] =
      &test_Clexulator::delta_site_eval_at_0_bfunc_3_1_2;
  m_delta_func_lists[0][28] =
      &test_Clexulator::delta_site_eval_at_0_bfunc_3_1_3;
  m_delta_func_lists[0][29] =
      &test_Clexulator::delta_site_eval_at_0_bfunc_3_1_4;
  m_delta_func_lists[0][30] =
      &test_Clexulator::delta_site_eval_at_0_bfunc_3_1_5;
  m_delta_func_lists[0][31] =
      &test_Clexulator::delta_site_eval_at_0_bfunc_3_2_0;
  m_delta_func_lists[0][32] =
      &test_Clexulator::delta_site_eval_at_0_bfunc_3_2_1;
  m_delta_func_lists[0][33] =
      &test_Clexulator::delta_site_eval_at_0_bfunc_3_2_2;
  m_delta_func_lists[0][34] =
      &test_Clexulator::delta_site_eval_at_0_bfunc_3_2_3;
  m_delta_func_lists[0][35] =
      &test_Clexulator::delta_site_eval_at_0_bfunc_3_2_4;
  m_delta_func_lists[0][36] =
      &test_Clexulator::delta_site_eval_at_0_bfunc_3_2_5;
  m_delta_func_lists[0][37] =
      &test_Clexulator::delta_site_eval_at_0_bfunc_3_3_0;
  m_delta_func_lists[0][38] =
      &test_Clexulator::delta_site_eval_at_0_bfunc_3_3_1;
  m_delta_func_lists[0][39] =
      &test_Clexulator::delta_site_eval_at_0_bfunc_3_3_2;
  m_delta_func_lists[0][40] =
      &test_Clexulator::delta_site_eval_at_0_bfunc_3_3_3;
  m_delta_func_lists[0][41] =
      &test_Clexulator::delta_site_eval_at_0_bfunc_3_3_4;
  m_delta_func_lists[0][42] =
      &test_Clexulator::delta_site_eval_at_0_bfunc_3_3_5;
  m_delta_func_lists[0][43] =
      &test_Clexulator::delta_site_eval_at_0_bfunc_3_3_6;
  m_delta_func_lists[0][44] =
      &test_Clexulator::delta_site_eval_at_0_bfunc_3_3_7;
  m_delta_func_lists[0][45] =
      &test_Clexulator::delta_site_eval_at_0_bfunc_3_4_0;
  m_delta_func_lists[0][46] =
      &test_Clexulator::delta_site_eval_at_0_bfunc_3_4_1;
  m_delta_func_lists[0][47] =
      &test_Clexulator::delta_site_eval_at_0_bfunc_3_4_2;
  m_delta_func_lists[0][48] =
      &test_Clexulator::delta_site_eval_at_0_bfunc_3_4_3;
  m_delta_func_lists[0][49] =
      &test_Clexulator::delta_site_eval_at_0_bfunc_3_4_4;
  m_delta_func_lists[0][50] =
      &test_Clexulator::delta_site_eval_at_0_bfunc_3_4_5;
  m_delta_func_lists[0][51] =
      &test_Clexulator::delta_site_eval_at_0_bfunc_3_4_6;
  m_delta_func_lists[0][52] =
      &test_Clexulator::delta_site_eval_at_0_bfunc_3_4_7;
  m_delta_func_lists[0][53] =
      &test_Clexulator::delta_site_eval_at_0_bfunc_4_0_0;
  m_delta_func_lists[0][54] =
      &test_Clexulator::delta_site_eval_at_0_bfunc_4_0_1;
  m_delta_func_lists[0][55] =
      &test_Clexulator::delta_site_eval_at_0_bfunc_4_0_2;
  m_delta_func_lists[0][56] =
      &test_Clexulator::delta_site_eval_at_0_bfunc_4_0_3;
  m_delta_func_lists[0][57] =
      &test_Clexulator::delta_site_eval_at_0_bfunc_4_0_4;
  m_delta_func_lists[0][58] =
      &test_Clexulator::delta_site_eval_at_0_bfunc_4_0_5;
  m_delta_func_lists[0][59] =
      &test_Clexulator::delta_site_eval_at_0_bfunc_4_0_6;
  m_delta_func_lists[0][60] =
      &test_Clexulator::delta_site_eval_at_0_bfunc_4_0_7;
  m_delta_func_lists[0][61] =
      &test_Clexulator::delta_site_eval_at_0_bfunc_4_0_8;
  m_delta_func_lists[0][62] =
      &test_Clexulator::delta_site_eval_at_0_bfunc_4_0_9;
  m_delta_func_lists[0][63] =
      &test_Clexulator::delta_site_eval_at_0_bfunc_4_0_10;
  m_delta_func_lists[0][64] =
      &test_Clexulator::delta_site_eval_at_0_bfunc_4_0_11;
  m_delta_func_lists[0][65] =
      &test_Clexulator::delta_site_eval_at_0_bfunc_4_1_0;
  m_delta_func_lists[0][66] =
      &test_Clexulator::delta_site_eval_at_0_bfunc_4_1_1;
  m_delta_func_lists[0][67] =
      &test_Clexulator::delta_site_eval_at_0_bfunc_4_1_2;
  m_delta_func_lists[0][68] =
      &test_Clexulator::delta_site_eval_at_0_bfunc_4_1_3;
  m_delta_func_lists[0][69] =
      &test_Clexulator::delta_site_eval_at_0_bfunc_4_1_4;
  m_delta_func_lists[0][70] =
      &test_Clexulator::delta_site_eval_at_0_bfunc_4_1_5;
  m_delta_func_lists[0][71] =
      &test_Clexulator::delta_site_eval_at_0_bfunc_4_1_6;
  m_delta_func_lists[0][72] =
      &test_Clexulator::delta_site_eval_at_0_bfunc_4_1_7;
  m_delta_func_lists[0][73] =
      &test_Clexulator::delta_site_eval_at_0_bfunc_4_1_8;
  m_delta_func_lists[0][74] =
      &test_Clexulator::delta_site_eval_at_0_bfunc_4_1_9;

  m_weight_matrix.row(0) << 2, 1, 1;
  m_weight_matrix.row(1) << 1, 2, 1;
  m_weight_matrix.row(2) << 1, 1, 2;

  m_neighborhood.insert(UnitCellCoord(0, -3, 0, 0));
  m_neighborhood.insert(UnitCellCoord(0, -3, 0, 3));
  m_neighborhood.insert(UnitCellCoord(0, -3, 1, 1));
  m_neighborhood.insert(UnitCellCoord(0, -3, 3, 0));
  m_neighborhood.insert(UnitCellCoord(0, -2, 0, 0));
  m_neighborhood.insert(UnitCellCoord(0, -2, 0, 1));
  m_neighborhood.insert(UnitCellCoord(0, -2, 0, 2));
  m_neighborhood.insert(UnitCellCoord(0, -2, 1, 0));
  m_neighborhood.insert(UnitCellCoord(0, -2, 1, 1));
  m_neighborhood.insert(UnitCellCoord(0, -2, 2, 0));
  m_neighborhood.insert(UnitCellCoord(0, -1, -1, -1));
  m_neighborhood.insert(UnitCellCoord(0, -1, -1, 0));
  m_neighborhood.insert(UnitCellCoord(0, -1, -1, 1));
  m_neighborhood.insert(UnitCellCoord(0, -1, -1, 2));
  m_neighborhood.insert(UnitCellCoord(0, -1, -1, 3));
  m_neighborhood.insert(UnitCellCoord(0, -1, 0, -1));
  m_neighborhood.insert(UnitCellCoord(0, -1, 0, 0));
  m_neighborhood.insert(UnitCellCoord(0, -1, 0, 1));
  m_neighborhood.insert(UnitCellCoord(0, -1, 0, 2));
  m_neighborhood.insert(UnitCellCoord(0, -1, 1, -1));
  m_neighborhood.insert(UnitCellCoord(0, -1, 1, 0));
  m_neighborhood.insert(UnitCellCoord(0, -1, 1, 1));
  m_neighborhood.insert(UnitCellCoord(0, -1, 2, -1));
  m_neighborhood.insert(UnitCellCoord(0, -1, 2, 0));
  m_neighborhood.insert(UnitCellCoord(0, -1, 3, -1));
  m_neighborhood.insert(UnitCellCoord(0, 0, -3, 0));
  m_neighborhood.insert(UnitCellCoord(0, 0, -3, 3));
  m_neighborhood.insert(UnitCellCoord(0, 0, -2, 0));
  m_neighborhood.insert(UnitCellCoord(0, 0, -2, 1));
  m_neighborhood.insert(UnitCellCoord(0, 0, -2, 2));
  m_neighborhood.insert(UnitCellCoord(0, 0, -1, -1));
  m_neighborhood.insert(UnitCellCoord(0, 0, -1, 0));
  m_neighborhood.insert(UnitCellCoord(0, 0, -1, 1));
  m_neighborhood.insert(UnitCellCoord(0, 0, -1, 2));
  m_neighborhood.insert(UnitCellCoord(0, 0, 0, -3));
  m_neighborhood.insert(UnitCellCoord(0, 0, 0, -2));
  m_neighborhood.insert(UnitCellCoord(0, 0, 0, -1));
  m_neighborhood.insert(UnitCellCoord(0, 0, 0, 0));
  m_neighborhood.insert(UnitCellCoord(0, 0, 0, 1));
  m_neighborhood.insert(UnitCellCoord(0, 0, 0, 2));
  m_neighborhood.insert(UnitCellCoord(0, 0, 0, 3));
  m_neighborhood.insert(UnitCellCoord(0, 0, 1, -2));
  m_neighborhood.insert(UnitCellCoord(0, 0, 1, -1));
  m_neighborhood.insert(UnitCellCoord(0, 0, 1, 0));
  m_neighborhood.insert(UnitCellCoord(0, 0, 1, 1));
  m_neighborhood.insert(UnitCellCoord(0, 0, 2, -2));
  m_neighborhood.insert(UnitCellCoord(0, 0, 2, -1));
  m_neighborhood.insert(UnitCellCoord(0, 0, 2, 0));
  m_neighborhood.insert(UnitCellCoord(0, 0, 3, -3));
  m_neighborhood.insert(UnitCellCoord(0, 0, 3, 0));
  m_neighborhood.insert(UnitCellCoord(0, 1, -3, 1));
  m_neighborhood.insert(UnitCellCoord(0, 1, -2, 0));
  m_neighborhood.insert(UnitCellCoord(0, 1, -2, 1));
  m_neighborhood.insert(UnitCellCoord(0, 1, -1, -1));
  m_neighborhood.insert(UnitCellCoord(0, 1, -1, 0));
  m_neighborhood.insert(UnitCellCoord(0, 1, -1, 1));
  m_neighborhood.insert(UnitCellCoord(0, 1, 0, -2));
  m_neighborhood.insert(UnitCellCoord(0, 1, 0, -1));
  m_neighborhood.insert(UnitCellCoord(0, 1, 0, 0));
  m_neighborhood.insert(UnitCellCoord(0, 1, 0, 1));
  m_neighborhood.insert(UnitCellCoord(0, 1, 1, -3));
  m_neighborhood.insert(UnitCellCoord(0, 1, 1, -2));
  m_neighborhood.insert(UnitCellCoord(0, 1, 1, -1));
  m_neighborhood.insert(UnitCellCoord(0, 1, 1, 0));
  m_neighborhood.insert(UnitCellCoord(0, 1, 1, 1));
  m_neighborhood.insert(UnitCellCoord(0, 2, -2, 0));
  m_neighborhood.insert(UnitCellCoord(0, 2, -1, -1));
  m_neighborhood.insert(UnitCellCoord(0, 2, -1, 0));
  m_neighborhood.insert(UnitCellCoord(0, 2, 0, -2));
  m_neighborhood.insert(UnitCellCoord(0, 2, 0, -1));
  m_neighborhood.insert(UnitCellCoord(0, 2, 0, 0));
  m_neighborhood.insert(UnitCellCoord(0, 3, -3, 0));
  m_neighborhood.insert(UnitCellCoord(0, 3, -1, -1));
  m_neighborhood.insert(UnitCellCoord(0, 3, 0, -3));
  m_neighborhood.insert(UnitCellCoord(0, 3, 0, 0));

  m_orbit_neighborhood.resize(corr_size());

  m_orbit_neighborhood[1].insert(UnitCellCoord(0, 0, 0, 0));
  m_orbit_neighborhood[2] = m_orbit_neighborhood[1];

  m_orbit_neighborhood[3].insert(UnitCellCoord(0, -1, 0, 0));
  m_orbit_neighborhood[3].insert(UnitCellCoord(0, -1, 0, 1));
  m_orbit_neighborhood[3].insert(UnitCellCoord(0, -1, 1, 0));
  m_orbit_neighborhood[3].insert(UnitCellCoord(0, 0, -1, 0));
  m_orbit_neighborhood[3].insert(UnitCellCoord(0, 0, -1, 1));
  m_orbit_neighborhood[3].insert(UnitCellCoord(0, 0, 0, -1));
  m_orbit_neighborhood[3].insert(UnitCellCoord(0, 0, 0, 0));
  m_orbit_neighborhood[3].insert(UnitCellCoord(0, 0, 0, 1));
  m_orbit_neighborhood[3].insert(UnitCellCoord(0, 0, 1, -1));
  m_orbit_neighborhood[3].insert(UnitCellCoord(0, 0, 1, 0));
  m_orbit_neighborhood[3].insert(UnitCellCoord(0, 1, -1, 0));
  m_orbit_neighborhood[3].insert(UnitCellCoord(0, 1, 0, -1));
  m_orbit_neighborhood[3].insert(UnitCellCoord(0, 1, 0, 0));
  m_orbit_neighborhood[4] = m_orbit_neighborhood[3];
  m_orbit_neighborhood[5] = m_orbit_neighborhood[3];

  m_orbit_neighborhood[6].insert(UnitCellCoord(0, -1, -1, 1));
  m_orbit_neighborhood[6].insert(UnitCellCoord(0, -1, 1, -1));
  m_orbit_neighborhood[6].insert(UnitCellCoord(0, -1, 1, 1));
  m_orbit_neighborhood[6].insert(UnitCellCoord(0, 0, 0, 0));
  m_orbit_neighborhood[6].insert(UnitCellCoord(0, 1, -1, -1));
  m_orbit_neighborhood[6].insert(UnitCellCoord(0, 1, -1, 1));
  m_orbit_neighborhood[6].insert(UnitCellCoord(0, 1, 1, -1));
  m_orbit_neighborhood[7] = m_orbit_neighborhood[6];
  m_orbit_neighborhood[8] = m_orbit_neighborhood[6];

  m_orbit_neighborhood[9].insert(UnitCellCoord(0, -2, 0, 1));
  m_orbit_neighborhood[9].insert(UnitCellCoord(0, -2, 1, 0));
  m_orbit_neighborhood[9].insert(UnitCellCoord(0, -2, 1, 1));
  m_orbit_neighborhood[9].insert(UnitCellCoord(0, -1, -1, 0));
  m_orbit_neighborhood[9].insert(UnitCellCoord(0, -1, -1, 2));
  m_orbit_neighborhood[9].insert(UnitCellCoord(0, -1, 0, -1));
  m_orbit_neighborhood[9].insert(UnitCellCoord(0, -1, 0, 2));
  m_orbit_neighborhood[9].insert(UnitCellCoord(0, -1, 2, -1));
  m_orbit_neighborhood[9].insert(UnitCellCoord(0, -1, 2, 0));
  m_orbit_neighborhood[9].insert(UnitCellCoord(0, 0, -2, 1));
  m_orbit_neighborhood[9].insert(UnitCellCoord(0, 0, -1, -1));
  m_orbit_neighborhood[9].insert(UnitCellCoord(0, 0, -1, 2));
  m_orbit_neighborhood[9].insert(UnitCellCoord(0, 0, 0, 0));
  m_orbit_neighborhood[9].insert(UnitCellCoord(0, 0, 1, -2));
  m_orbit_neighborhood[9].insert(UnitCellCoord(0, 0, 1, 1));
  m_orbit_neighborhood[9].insert(UnitCellCoord(0, 0, 2, -1));
  m_orbit_neighborhood[9].insert(UnitCellCoord(0, 1, -2, 0));
  m_orbit_neighborhood[9].insert(UnitCellCoord(0, 1, -2, 1));
  m_orbit_neighborhood[9].insert(UnitCellCoord(0, 1, 0, -2));
  m_orbit_neighborhood[9].insert(UnitCellCoord(0, 1, 0, 1));
  m_orbit_neighborhood[9].insert(UnitCellCoord(0, 1, 1, -2));
  m_orbit_neighborhood[9].insert(UnitCellCoord(0, 1, 1, 0));
  m_orbit_neighborhood[9].insert(UnitCellCoord(0, 2, -1, -1));
  m_orbit_neighborhood[9].insert(UnitCellCoord(0, 2, -1, 0));
  m_orbit_neighborhood[9].insert(UnitCellCoord(0, 2, 0, -1));
  m_orbit_neighborhood[10] = m_orbit_neighborhood[9];
  m_orbit_neighborhood[11] = m_orbit_neighborhood[9];

  m_orbit_neighborhood[12].insert(UnitCellCoord(0, -2, 0, 0));
  m_orbit_neighborhood[12].insert(UnitCellCoord(0, -2, 0, 2));
  m_orbit_neighborhood[12].insert(UnitCellCoord(0, -2, 2, 0));
  m_orbit_neighborhood[12].insert(UnitCellCoord(0, 0, -2, 0));
  m_orbit_neighborhood[12].insert(UnitCellCoord(0, 0, -2, 2));
  m_orbit_neighborhood[12].insert(UnitCellCoord(0, 0, 0, -2));
  m_orbit_neighborhood[12].insert(UnitCellCoord(0, 0, 0, 0));
  m_orbit_neighborhood[12].insert(UnitCellCoord(0, 0, 0, 2));
  m_orbit_neighborhood[12].insert(UnitCellCoord(0, 0, 2, -2));
  m_orbit_neighborhood[12].insert(UnitCellCoord(0, 0, 2, 0));
  m_orbit_neighborhood[12].insert(UnitCellCoord(0, 2, -2, 0));
  m_orbit_neighborhood[12].insert(UnitCellCoord(0, 2, 0, -2));
  m_orbit_neighborhood[12].insert(UnitCellCoord(0, 2, 0, 0));
  m_orbit_neighborhood[13] = m_orbit_neighborhood[12];
  m_orbit_neighborhood[14] = m_orbit_neighborhood[12];

  m_orbit_neighborhood[15].insert(UnitCellCoord(0, -3, 1, 1));
  m_orbit_neighborhood[15].insert(UnitCellCoord(0, -1, -1, -1));
  m_orbit_neighborhood[15].insert(UnitCellCoord(0, -1, -1, 3));
  m_orbit_neighborhood[15].insert(UnitCellCoord(0, -1, 3, -1));
  m_orbit_neighborhood[15].insert(UnitCellCoord(0, 0, 0, 0));
  m_orbit_neighborhood[15].insert(UnitCellCoord(0, 1, -3, 1));
  m_orbit_neighborhood[15].insert(UnitCellCoord(0, 1, 1, -3));
  m_orbit_neighborhood[15].insert(UnitCellCoord(0, 1, 1, 1));
  m_orbit_neighborhood[15].insert(UnitCellCoord(0, 3, -1, -1));
  m_orbit_neighborhood[16] = m_orbit_neighborhood[15];
  m_orbit_neighborhood[17] = m_orbit_neighborhood[15];

  m_orbit_neighborhood[18].insert(UnitCellCoord(0, -3, 0, 0));
  m_orbit_neighborhood[18].insert(UnitCellCoord(0, -3, 0, 3));
  m_orbit_neighborhood[18].insert(UnitCellCoord(0, -3, 3, 0));
  m_orbit_neighborhood[18].insert(UnitCellCoord(0, 0, -3, 0));
  m_orbit_neighborhood[18].insert(UnitCellCoord(0, 0, -3, 3));
  m_orbit_neighborhood[18].insert(UnitCellCoord(0, 0, 0, -3));
  m_orbit_neighborhood[18].insert(UnitCellCoord(0, 0, 0, 0));
  m_orbit_neighborhood[18].insert(UnitCellCoord(0, 0, 0, 3));
  m_orbit_neighborhood[18].insert(UnitCellCoord(0, 0, 3, -3));
  m_orbit_neighborhood[18].insert(UnitCellCoord(0, 0, 3, 0));
  m_orbit_neighborhood[18].insert(UnitCellCoord(0, 3, -3, 0));
  m_orbit_neighborhood[18].insert(UnitCellCoord(0, 3, 0, -3));
  m_orbit_neighborhood[18].insert(UnitCellCoord(0, 3, 0, 0));
  m_orbit_neighborhood[19] = m_orbit_neighborhood[18];
  m_orbit_neighborhood[20] = m_orbit_neighborhood[18];

  m_orbit_neighborhood[21].insert(UnitCellCoord(0, -1, 0, 0));
  m_orbit_neighborhood[21].insert(UnitCellCoord(0, -1, 0, 1));
  m_orbit_neighborhood[21].insert(UnitCellCoord(0, -1, 1, 0));
  m_orbit_neighborhood[21].insert(UnitCellCoord(0, 0, -1, 0));
  m_orbit_neighborhood[21].insert(UnitCellCoord(0, 0, -1, 1));
  m_orbit_neighborhood[21].insert(UnitCellCoord(0, 0, 0, -1));
  m_orbit_neighborhood[21].insert(UnitCellCoord(0, 0, 0, 0));
  m_orbit_neighborhood[21].insert(UnitCellCoord(0, 0, 0, 1));
  m_orbit_neighborhood[21].insert(UnitCellCoord(0, 0, 1, -1));
  m_orbit_neighborhood[21].insert(UnitCellCoord(0, 0, 1, 0));
  m_orbit_neighborhood[21].insert(UnitCellCoord(0, 1, -1, 0));
  m_orbit_neighborhood[21].insert(UnitCellCoord(0, 1, 0, -1));
  m_orbit_neighborhood[21].insert(UnitCellCoord(0, 1, 0, 0));
  m_orbit_neighborhood[22] = m_orbit_neighborhood[21];
  m_orbit_neighborhood[23] = m_orbit_neighborhood[21];
  m_orbit_neighborhood[24] = m_orbit_neighborhood[21];

  m_orbit_neighborhood[25].insert(UnitCellCoord(0, -2, 0, 1));
  m_orbit_neighborhood[25].insert(UnitCellCoord(0, -2, 1, 0));
  m_orbit_neighborhood[25].insert(UnitCellCoord(0, -2, 1, 1));
  m_orbit_neighborhood[25].insert(UnitCellCoord(0, -1, -1, 0));
  m_orbit_neighborhood[25].insert(UnitCellCoord(0, -1, -1, 2));
  m_orbit_neighborhood[25].insert(UnitCellCoord(0, -1, 0, -1));
  m_orbit_neighborhood[25].insert(UnitCellCoord(0, -1, 0, 0));
  m_orbit_neighborhood[25].insert(UnitCellCoord(0, -1, 0, 1));
  m_orbit_neighborhood[25].insert(UnitCellCoord(0, -1, 0, 2));
  m_orbit_neighborhood[25].insert(UnitCellCoord(0, -1, 1, 0));
  m_orbit_neighborhood[25].insert(UnitCellCoord(0, -1, 2, -1));
  m_orbit_neighborhood[25].insert(UnitCellCoord(0, -1, 2, 0));
  m_orbit_neighborhood[25].insert(UnitCellCoord(0, 0, -2, 1));
  m_orbit_neighborhood[25].insert(UnitCellCoord(0, 0, -1, -1));
  m_orbit_neighborhood[25].insert(UnitCellCoord(0, 0, -1, 0));
  m_orbit_neighborhood[25].insert(UnitCellCoord(0, 0, -1, 1));
  m_orbit_neighborhood[25].insert(UnitCellCoord(0, 0, -1, 2));
  m_orbit_neighborhood[25].insert(UnitCellCoord(0, 0, 0, -1));
  m_orbit_neighborhood[25].insert(UnitCellCoord(0, 0, 0, 0));
  m_orbit_neighborhood[25].insert(UnitCellCoord(0, 0, 0, 1));
  m_orbit_neighborhood[25].insert(UnitCellCoord(0, 0, 1, -2));
  m_orbit_neighborhood[25].insert(UnitCellCoord(0, 0, 1, -1));
  m_orbit_neighborhood[25].insert(UnitCellCoord(0, 0, 1, 0));
  m_orbit_neighborhood[25].insert(UnitCellCoord(0, 0, 1, 1));
  m_orbit_neighborhood[25].insert(UnitCellCoord(0, 0, 2, -1));
  m_orbit_neighborhood[25].insert(UnitCellCoord(0, 1, -2, 0));
  m_orbit_neighborhood[25].insert(UnitCellCoord(0, 1, -2, 1));
  m_orbit_neighborhood[25].insert(UnitCellCoord(0, 1, -1, 0));
  m_orbit_neighborhood[25].insert(UnitCellCoord(0, 1, 0, -2));
  m_orbit_neighborhood[25].insert(UnitCellCoord(0, 1, 0, -1));
  m_orbit_neighborhood[25].insert(UnitCellCoord(0, 1, 0, 0));
  m_orbit_neighborhood[25].insert(UnitCellCoord(0, 1, 0, 1));
  m_orbit_neighborhood[25].insert(UnitCellCoord(0, 1, 1, -2));
  m_orbit_neighborhood[25].insert(UnitCellCoord(0, 1, 1, 0));
  m_orbit_neighborhood[25].insert(UnitCellCoord(0, 2, -1, -1));
  m_orbit_neighborhood[25].insert(UnitCellCoord(0, 2, -1, 0));
  m_orbit_neighborhood[25].insert(UnitCellCoord(0, 2, 0, -1));
  m_orbit_neighborhood[26] = m_orbit_neighborhood[25];
  m_orbit_neighborhood[27] = m_orbit_neighborhood[25];
  m_orbit_neighborhood[28] = m_orbit_neighborhood[25];
  m_orbit_neighborhood[29] = m_orbit_neighborhood[25];
  m_orbit_neighborhood[30] = m_orbit_neighborhood[25];

  m_orbit_neighborhood[31].insert(UnitCellCoord(0, -2, 0, 0));
  m_orbit_neighborhood[31].insert(UnitCellCoord(0, -2, 0, 2));
  m_orbit_neighborhood[31].insert(UnitCellCoord(0, -2, 2, 0));
  m_orbit_neighborhood[31].insert(UnitCellCoord(0, -1, 0, 0));
  m_orbit_neighborhood[31].insert(UnitCellCoord(0, -1, 0, 1));
  m_orbit_neighborhood[31].insert(UnitCellCoord(0, -1, 1, 0));
  m_orbit_neighborhood[31].insert(UnitCellCoord(0, 0, -2, 0));
  m_orbit_neighborhood[31].insert(UnitCellCoord(0, 0, -2, 2));
  m_orbit_neighborhood[31].insert(UnitCellCoord(0, 0, -1, 0));
  m_orbit_neighborhood[31].insert(UnitCellCoord(0, 0, -1, 1));
  m_orbit_neighborhood[31].insert(UnitCellCoord(0, 0, 0, -2));
  m_orbit_neighborhood[31].insert(UnitCellCoord(0, 0, 0, -1));
  m_orbit_neighborhood[31].insert(UnitCellCoord(0, 0, 0, 0));
  m_orbit_neighborhood[31].insert(UnitCellCoord(0, 0, 0, 1));
  m_orbit_neighborhood[31].insert(UnitCellCoord(0, 0, 0, 2));
  m_orbit_neighborhood[31].insert(UnitCellCoord(0, 0, 1, -1));
  m_orbit_neighborhood[31].insert(UnitCellCoord(0, 0, 1, 0));
  m_orbit_neighborhood[31].insert(UnitCellCoord(0, 0, 2, -2));
  m_orbit_neighborhood[31].insert(UnitCellCoord(0, 0, 2, 0));
  m_orbit_neighborhood[31].insert(UnitCellCoord(0, 1, -1, 0));
  m_orbit_neighborhood[31].insert(UnitCellCoord(0, 1, 0, -1));
  m_orbit_neighborhood[31].insert(UnitCellCoord(0, 1, 0, 0));
  m_orbit_neighborhood[31].insert(UnitCellCoord(0, 2, -2, 0));
  m_orbit_neighborhood[31].insert(UnitCellCoord(0, 2, 0, -2));
  m_orbit_neighborhood[31].insert(UnitCellCoord(0, 2, 0, 0));
  m_orbit_neighborhood[32] = m_orbit_neighborhood[31];
  m_orbit_neighborhood[33] = m_orbit_neighborhood[31];
  m_orbit_neighborhood[34] = m_orbit_neighborhood[31];
  m_orbit_neighborhood[35] = m_orbit_neighborhood[31];
  m_orbit_neighborhood[36] = m_orbit_neighborhood[31];

  m_orbit_neighborhood[37].insert(UnitCellCoord(0, -3, 1, 1));
  m_orbit_neighborhood[37].insert(UnitCellCoord(0, -2, 0, 1));
  m_orbit_neighborhood[37].insert(UnitCellCoord(0, -2, 1, 0));
  m_orbit_neighborhood[37].insert(UnitCellCoord(0, -2, 1, 1));
  m_orbit_neighborhood[37].insert(UnitCellCoord(0, -1, -1, -1));
  m_orbit_neighborhood[37].insert(UnitCellCoord(0, -1, -1, 0));
  m_orbit_neighborhood[37].insert(UnitCellCoord(0, -1, -1, 2));
  m_orbit_neighborhood[37].insert(UnitCellCoord(0, -1, -1, 3));
  m_orbit_neighborhood[37].insert(UnitCellCoord(0, -1, 0, -1));
  m_orbit_neighborhood[37].insert(UnitCellCoord(0, -1, 0, 0));
  m_orbit_neighborhood[37].insert(UnitCellCoord(0, -1, 0, 1));
  m_orbit_neighborhood[37].insert(UnitCellCoord(0, -1, 0, 2));
  m_orbit_neighborhood[37].insert(UnitCellCoord(0, -1, 1, 0));
  m_orbit_neighborhood[37].insert(UnitCellCoord(0, -1, 2, -1));
  m_orbit_neighborhood[37].insert(UnitCellCoord(0, -1, 2, 0));
  m_orbit_neighborhood[37].insert(UnitCellCoord(0, -1, 3, -1));
  m_orbit_neighborhood[37].insert(UnitCellCoord(0, 0, -2, 1));
  m_orbit_neighborhood[37].insert(UnitCellCoord(0, 0, -1, -1));
  m_orbit_neighborhood[37].insert(UnitCellCoord(0, 0, -1, 0));
  m_orbit_neighborhood[37].insert(UnitCellCoord(0, 0, -1, 1));
  m_orbit_neighborhood[37].insert(UnitCellCoord(0, 0, -1, 2));
  m_orbit_neighborhood[37].insert(UnitCellCoord(0, 0, 0, -1));
  m_orbit_neighborhood[37].insert(UnitCellCoord(0, 0, 0, 0));
  m_orbit_neighborhood[37].insert(UnitCellCoord(0, 0, 0, 1));
  m_orbit_neighborhood[37].insert(UnitCellCoord(0, 0, 1, -2));
  m_orbit_neighborhood[37].insert(UnitCellCoord(0, 0, 1, -1));
  m_orbit_neighborhood[37].insert(UnitCellCoord(0, 0, 1, 0));
  m_orbit_neighborhood[37].insert(UnitCellCoord(0, 0, 1, 1));
  m_orbit_neighborhood[37].insert(UnitCellCoord(0, 0, 2, -1));
  m_orbit_neighborhood[37].insert(UnitCellCoord(0, 1, -3, 1));
  m_orbit_neighborhood[37].insert(UnitCellCoord(0, 1, -2, 0));
  m_orbit_neighborhood[37].insert(UnitCellCoord(0, 1, -2, 1));
  m_orbit_neighborhood[37].insert(UnitCellCoord(0, 1, -1, 0));
  m_orbit_neighborhood[37].insert(UnitCellCoord(0, 1, 0, -2));
  m_orbit_neighborhood[37].insert(UnitCellCoord(0, 1, 0, -1));
  m_orbit_neighborhood[37].insert(UnitCellCoord(0, 1, 0, 0));
  m_orbit_neighborhood[37].insert(UnitCellCoord(0, 1, 0, 1));
  m_orbit_neighborhood[37].insert(UnitCellCoord(0, 1, 1, -3));
  m_orbit_neighborhood[37].insert(UnitCellCoord(0, 1, 1, -2));
  m_orbit_neighborhood[37].insert(UnitCellCoord(0, 1, 1, 0));
  m_orbit_neighborhood[37].insert(UnitCellCoord(0, 1, 1, 1));
  m_orbit_neighborhood[37].insert(UnitCellCoord(0, 2, -1, -1));
  m_orbit_neighborhood[37].insert(UnitCellCoord(0, 2, -1, 0));
  m_orbit_neighborhood[37].insert(UnitCellCoord(0, 2, 0, -1));
  m_orbit_neighborhood[37].insert(UnitCellCoord(0, 3, -1, -1));
  m_orbit_neighborhood[38] = m_orbit_neighborhood[37];
  m_orbit_neighborhood[39] = m_orbit_neighborhood[37];
  m_orbit_neighborhood[40] = m_orbit_neighborhood[37];
  m_orbit_neighborhood[41] = m_orbit_neighborhood[37];
  m_orbit_neighborhood[42] = m_orbit_neighborhood[37];
  m_orbit_neighborhood[43] = m_orbit_neighborhood[37];
  m_orbit_neighborhood[44] = m_orbit_neighborhood[37];

  m_orbit_neighborhood[45].insert(UnitCellCoord(0, -3, 0, 0));
  m_orbit_neighborhood[45].insert(UnitCellCoord(0, -3, 0, 3));
  m_orbit_neighborhood[45].insert(UnitCellCoord(0, -3, 3, 0));
  m_orbit_neighborhood[45].insert(UnitCellCoord(0, -2, 0, 0));
  m_orbit_neighborhood[45].insert(UnitCellCoord(0, -2, 0, 2));
  m_orbit_neighborhood[45].insert(UnitCellCoord(0, -2, 2, 0));
  m_orbit_neighborhood[45].insert(UnitCellCoord(0, -1, 0, 0));
  m_orbit_neighborhood[45].insert(UnitCellCoord(0, -1, 0, 1));
  m_orbit_neighborhood[45].insert(UnitCellCoord(0, -1, 1, 0));
  m_orbit_neighborhood[45].insert(UnitCellCoord(0, 0, -3, 0));
  m_orbit_neighborhood[45].insert(UnitCellCoord(0, 0, -3, 3));
  m_orbit_neighborhood[45].insert(UnitCellCoord(0, 0, -2, 0));
  m_orbit_neighborhood[45].insert(UnitCellCoord(0, 0, -2, 2));
  m_orbit_neighborhood[45].insert(UnitCellCoord(0, 0, -1, 0));
  m_orbit_neighborhood[45].insert(UnitCellCoord(0, 0, -1, 1));
  m_orbit_neighborhood[45].insert(UnitCellCoord(0, 0, 0, -3));
  m_orbit_neighborhood[45].insert(UnitCellCoord(0, 0, 0, -2));
  m_orbit_neighborhood[45].insert(UnitCellCoord(0, 0, 0, -1));
  m_orbit_neighborhood[45].insert(UnitCellCoord(0, 0, 0, 0));
  m_orbit_neighborhood[45].insert(UnitCellCoord(0, 0, 0, 1));
  m_orbit_neighborhood[45].insert(UnitCellCoord(0, 0, 0, 2));
  m_orbit_neighborhood[45].insert(UnitCellCoord(0, 0, 0, 3));
  m_orbit_neighborhood[45].insert(UnitCellCoord(0, 0, 1, -1));
  m_orbit_neighborhood[45].insert(UnitCellCoord(0, 0, 1, 0));
  m_orbit_neighborhood[45].insert(UnitCellCoord(0, 0, 2, -2));
  m_orbit_neighborhood[45].insert(UnitCellCoord(0, 0, 2, 0));
  m_orbit_neighborhood[45].insert(UnitCellCoord(0, 0, 3, -3));
  m_orbit_neighborhood[45].insert(UnitCellCoord(0, 0, 3, 0));
  m_orbit_neighborhood[45].insert(UnitCellCoord(0, 1, -1, 0));
  m_orbit_neighborhood[45].insert(UnitCellCoord(0, 1, 0, -1));
  m_orbit_neighborhood[45].insert(UnitCellCoord(0, 1, 0, 0));
  m_orbit_neighborhood[45].insert(UnitCellCoord(0, 2, -2, 0));
  m_orbit_neighborhood[45].insert(UnitCellCoord(0, 2, 0, -2));
  m_orbit_neighborhood[45].insert(UnitCellCoord(0, 2, 0, 0));
  m_orbit_neighborhood[45].insert(UnitCellCoord(0, 3, -3, 0));
  m_orbit_neighborhood[45].insert(UnitCellCoord(0, 3, 0, -3));
  m_orbit_neighborhood[45].insert(UnitCellCoord(0, 3, 0, 0));
  m_orbit_neighborhood[46] = m_orbit_neighborhood[45];
  m_orbit_neighborhood[47] = m_orbit_neighborhood[45];
  m_orbit_neighborhood[48] = m_orbit_neighborhood[45];
  m_orbit_neighborhood[49] = m_orbit_neighborhood[45];
  m_orbit_neighborhood[50] = m_orbit_neighborhood[45];
  m_orbit_neighborhood[51] = m_orbit_neighborhood[45];
  m_orbit_neighborhood[52] = m_orbit_neighborhood[45];

  m_orbit_neighborhood[53].insert(UnitCellCoord(0, -3, 1, 1));
  m_orbit_neighborhood[53].insert(UnitCellCoord(0, -2, 0, 1));
  m_orbit_neighborhood[53].insert(UnitCellCoord(0, -2, 1, 0));
  m_orbit_neighborhood[53].insert(UnitCellCoord(0, -2, 1, 1));
  m_orbit_neighborhood[53].insert(UnitCellCoord(0, -1, -1, -1));
  m_orbit_neighborhood[53].insert(UnitCellCoord(0, -1, -1, 0));
  m_orbit_neighborhood[53].insert(UnitCellCoord(0, -1, -1, 2));
  m_orbit_neighborhood[53].insert(UnitCellCoord(0, -1, -1, 3));
  m_orbit_neighborhood[53].insert(UnitCellCoord(0, -1, 0, -1));
  m_orbit_neighborhood[53].insert(UnitCellCoord(0, -1, 0, 0));
  m_orbit_neighborhood[53].insert(UnitCellCoord(0, -1, 0, 1));
  m_orbit_neighborhood[53].insert(UnitCellCoord(0, -1, 0, 2));
  m_orbit_neighborhood[53].insert(UnitCellCoord(0, -1, 1, 0));
  m_orbit_neighborhood[53].insert(UnitCellCoord(0, -1, 2, -1));
  m_orbit_neighborhood[53].insert(UnitCellCoord(0, -1, 2, 0));
  m_orbit_neighborhood[53].insert(UnitCellCoord(0, -1, 3, -1));
  m_orbit_neighborhood[53].insert(UnitCellCoord(0, 0, -2, 1));
  m_orbit_neighborhood[53].insert(UnitCellCoord(0, 0, -1, -1));
  m_orbit_neighborhood[53].insert(UnitCellCoord(0, 0, -1, 0));
  m_orbit_neighborhood[53].insert(UnitCellCoord(0, 0, -1, 1));
  m_orbit_neighborhood[53].insert(UnitCellCoord(0, 0, -1, 2));
  m_orbit_neighborhood[53].insert(UnitCellCoord(0, 0, 0, -1));
  m_orbit_neighborhood[53].insert(UnitCellCoord(0, 0, 0, 0));
  m_orbit_neighborhood[53].insert(UnitCellCoord(0, 0, 0, 1));
  m_orbit_neighborhood[53].insert(UnitCellCoord(0, 0, 1, -2));
  m_orbit_neighborhood[53].insert(UnitCellCoord(0, 0, 1, -1));
  m_orbit_neighborhood[53].insert(UnitCellCoord(0, 0, 1, 0));
  m_orbit_neighborhood[53].insert(UnitCellCoord(0, 0, 1, 1));
  m_orbit_neighborhood[53].insert(UnitCellCoord(0, 0, 2, -1));
  m_orbit_neighborhood[53].insert(UnitCellCoord(0, 1, -3, 1));
  m_orbit_neighborhood[53].insert(UnitCellCoord(0, 1, -2, 0));
  m_orbit_neighborhood[53].insert(UnitCellCoord(0, 1, -2, 1));
  m_orbit_neighborhood[53].insert(UnitCellCoord(0, 1, -1, 0));
  m_orbit_neighborhood[53].insert(UnitCellCoord(0, 1, 0, -2));
  m_orbit_neighborhood[53].insert(UnitCellCoord(0, 1, 0, -1));
  m_orbit_neighborhood[53].insert(UnitCellCoord(0, 1, 0, 0));
  m_orbit_neighborhood[53].insert(UnitCellCoord(0, 1, 0, 1));
  m_orbit_neighborhood[53].insert(UnitCellCoord(0, 1, 1, -3));
  m_orbit_neighborhood[53].insert(UnitCellCoord(0, 1, 1, -2));
  m_orbit_neighborhood[53].insert(UnitCellCoord(0, 1, 1, 0));
  m_orbit_neighborhood[53].insert(UnitCellCoord(0, 1, 1, 1));
  m_orbit_neighborhood[53].insert(UnitCellCoord(0, 2, -1, -1));
  m_orbit_neighborhood[53].insert(UnitCellCoord(0, 2, -1, 0));
  m_orbit_neighborhood[53].insert(UnitCellCoord(0, 2, 0, -1));
  m_orbit_neighborhood[53].insert(UnitCellCoord(0, 3, -1, -1));
  m_orbit_neighborhood[54] = m_orbit_neighborhood[53];
  m_orbit_neighborhood[55] = m_orbit_neighborhood[53];
  m_orbit_neighborhood[56] = m_orbit_neighborhood[53];
  m_orbit_neighborhood[57] = m_orbit_neighborhood[53];
  m_orbit_neighborhood[58] = m_orbit_neighborhood[53];
  m_orbit_neighborhood[59] = m_orbit_neighborhood[53];
  m_orbit_neighborhood[60] = m_orbit_neighborhood[53];
  m_orbit_neighborhood[61] = m_orbit_neighborhood[53];
  m_orbit_neighborhood[62] = m_orbit_neighborhood[53];
  m_orbit_neighborhood[63] = m_orbit_neighborhood[53];
  m_orbit_neighborhood[64] = m_orbit_neighborhood[53];

  m_orbit_neighborhood[65].insert(UnitCellCoord(0, -3, 0, 0));
  m_orbit_neighborhood[65].insert(UnitCellCoord(0, -3, 0, 3));
  m_orbit_neighborhood[65].insert(UnitCellCoord(0, -3, 3, 0));
  m_orbit_neighborhood[65].insert(UnitCellCoord(0, -2, 0, 0));
  m_orbit_neighborhood[65].insert(UnitCellCoord(0, -2, 0, 2));
  m_orbit_neighborhood[65].insert(UnitCellCoord(0, -2, 2, 0));
  m_orbit_neighborhood[65].insert(UnitCellCoord(0, -1, 0, 0));
  m_orbit_neighborhood[65].insert(UnitCellCoord(0, -1, 0, 1));
  m_orbit_neighborhood[65].insert(UnitCellCoord(0, -1, 1, 0));
  m_orbit_neighborhood[65].insert(UnitCellCoord(0, 0, -3, 0));
  m_orbit_neighborhood[65].insert(UnitCellCoord(0, 0, -3, 3));
  m_orbit_neighborhood[65].insert(UnitCellCoord(0, 0, -2, 0));
  m_orbit_neighborhood[65].insert(UnitCellCoord(0, 0, -2, 2));
  m_orbit_neighborhood[65].insert(UnitCellCoord(0, 0, -1, 0));
  m_orbit_neighborhood[65].insert(UnitCellCoord(0, 0, -1, 1));
  m_orbit_neighborhood[65].insert(UnitCellCoord(0, 0, 0, -3));
  m_orbit_neighborhood[65].insert(UnitCellCoord(0, 0, 0, -2));
  m_orbit_neighborhood[65].insert(UnitCellCoord(0, 0, 0, -1));
  m_orbit_neighborhood[65].insert(UnitCellCoord(0, 0, 0, 0));
  m_orbit_neighborhood[65].insert(UnitCellCoord(0, 0, 0, 1));
  m_orbit_neighborhood[65].insert(UnitCellCoord(0, 0, 0, 2));
  m_orbit_neighborhood[65].insert(UnitCellCoord(0, 0, 0, 3));
  m_orbit_neighborhood[65].insert(UnitCellCoord(0, 0, 1, -1));
  m_orbit_neighborhood[65].insert(UnitCellCoord(0, 0, 1, 0));
  m_orbit_neighborhood[65].insert(UnitCellCoord(0, 0, 2, -2));
  m_orbit_neighborhood[65].insert(UnitCellCoord(0, 0, 2, 0));
  m_orbit_neighborhood[65].insert(UnitCellCoord(0, 0, 3, -3));
  m_orbit_neighborhood[65].insert(UnitCellCoord(0, 0, 3, 0));
  m_orbit_neighborhood[65].insert(UnitCellCoord(0, 1, -1, 0));
  m_orbit_neighborhood[65].insert(UnitCellCoord(0, 1, 0, -1));
  m_orbit_neighborhood[65].insert(UnitCellCoord(0, 1, 0, 0));
  m_orbit_neighborhood[65].insert(UnitCellCoord(0, 2, -2, 0));
  m_orbit_neighborhood[65].insert(UnitCellCoord(0, 2, 0, -2));
  m_orbit_neighborhood[65].insert(UnitCellCoord(0, 2, 0, 0));
  m_orbit_neighborhood[65].insert(UnitCellCoord(0, 3, -3, 0));
  m_orbit_neighborhood[65].insert(UnitCellCoord(0, 3, 0, -3));
  m_orbit_neighborhood[65].insert(UnitCellCoord(0, 3, 0, 0));
  m_orbit_neighborhood[66] = m_orbit_neighborhood[65];
  m_orbit_neighborhood[67] = m_orbit_neighborhood[65];
  m_orbit_neighborhood[68] = m_orbit_neighborhood[65];
  m_orbit_neighborhood[69] = m_orbit_neighborhood[65];
  m_orbit_neighborhood[70] = m_orbit_neighborhood[65];
  m_orbit_neighborhood[71] = m_orbit_neighborhood[65];
  m_orbit_neighborhood[72] = m_orbit_neighborhood[65];
  m_orbit_neighborhood[73] = m_orbit_neighborhood[65];
  m_orbit_neighborhood[74] = m_orbit_neighborhood[65];
}

test_Clexulator::~test_Clexulator() {
  // nothing here for now
}

/// \brief Calculate contribution to global correlations from one unit cell
void test_Clexulator::calc_global_corr_contribution(double *corr_begin) const {
  for (size_type i = 0; i < corr_size(); i++) {
    *(corr_begin + i) = (this->*m_orbit_func_list[i])();
  }
}

/// \brief Calculate contribution to select global correlations from one unit
/// cell
void test_Clexulator::calc_restricted_global_corr_contribution(
    double *corr_begin, size_type const *ind_list_begin,
    size_type const *ind_list_end) const {
  for (; ind_list_begin < ind_list_end; ind_list_begin++) {
    *(corr_begin + *ind_list_begin) =
        (this->*m_orbit_func_list[*ind_list_begin])();
  }
}

/// \brief Calculate point correlations about basis site 'b_index'
void test_Clexulator::calc_point_corr(int b_index, double *corr_begin) const {
  for (size_type i = 0; i < corr_size(); i++) {
    *(corr_begin + i) = (this->*m_flower_func_lists[b_index][i])();
  }
}

/// \brief Calculate select point correlations about basis site 'b_index'
void test_Clexulator::calc_restricted_point_corr(
    int b_index, double *corr_begin, size_type const *ind_list_begin,
    size_type const *ind_list_end) const {
  for (; ind_list_begin < ind_list_end; ind_list_begin++) {
    *(corr_begin + *ind_list_begin) =
        (this->*m_flower_func_lists[b_index][*ind_list_begin])();
  }
}

/// \brief Calculate the change in point correlations due to changing an
/// occupant
void test_Clexulator::calc_delta_point_corr(int b_index, int occ_i, int occ_f,
                                            double *corr_begin) const {
  for (size_type i = 0; i < corr_size(); i++) {
    *(corr_begin + i) = (this->*m_delta_func_lists[b_index][i])(occ_i, occ_f);
  }
}

/// \brief Calculate the change in select point correlations due to changing an
/// occupant
void test_Clexulator::calc_restricted_delta_point_corr(
    int b_index, int occ_i, int occ_f, double *corr_begin,
    size_type const *ind_list_begin, size_type const *ind_list_end) const {
  for (; ind_list_begin < ind_list_end; ind_list_begin++) {
    *(corr_begin + *ind_list_begin) =
        (this->*m_delta_func_lists[b_index][*ind_list_begin])(occ_i, occ_f);
  }
}

// Basis functions for empty cluster:
double test_Clexulator::eval_bfunc_0_0_0() const { return (1); }

/**** Basis functions for orbit 1, 0****
#Points: 1
MaxLength: 0  MinLength: 0
 0.0000000   0.0000000   0.0000000 A B C
****/
double test_Clexulator::eval_bfunc_1_0_0() const { return (occ_func_0_0(0)); }
double test_Clexulator::eval_bfunc_1_0_1() const { return (occ_func_0_1(0)); }

double test_Clexulator::site_eval_at_0_bfunc_1_0_0() const {
  return (occ_func_0_0(0));
}
double test_Clexulator::site_eval_at_0_bfunc_1_0_1() const {
  return (occ_func_0_1(0));
}

double test_Clexulator::delta_site_eval_at_0_bfunc_1_0_0(int occ_i,
                                                         int occ_f) const {
  return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]);
}
double test_Clexulator::delta_site_eval_at_0_bfunc_1_0_1(int occ_i,
                                                         int occ_f) const {
  return (m_occ_func_0_1[occ_f] - m_occ_func_0_1[occ_i]);
}

/**** Basis functions for orbit 2, 0****
#Points: 2
MaxLength: 2.8284271  MinLength: 2.8284271
 0.0000000   0.0000000   0.0000000 A B C
 0.0000000   0.0000000  -1.0000000 A B C
****/
double test_Clexulator::eval_bfunc_2_0_0() const {
  return ((occ_func_0_0(0) * occ_func_0_0(6)) +
          (occ_func_0_0(0) * occ_func_0_0(4)) +
          (occ_func_0_0(0) * occ_func_0_0(8)) +
          (occ_func_0_0(0) * occ_func_0_0(3)) +
          (occ_func_0_0(0) * occ_func_0_0(1)) +
          (occ_func_0_0(0) * occ_func_0_0(11))) /
         6.0;
}
double test_Clexulator::eval_bfunc_2_0_1() const {
  return (((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(6) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(6))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(4) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(4))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(8) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(8))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(3) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(3))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(1) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(1))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(11) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(11)))) /
         6.0;
}
double test_Clexulator::eval_bfunc_2_0_2() const {
  return ((occ_func_0_1(0) * occ_func_0_1(6)) +
          (occ_func_0_1(0) * occ_func_0_1(4)) +
          (occ_func_0_1(0) * occ_func_0_1(8)) +
          (occ_func_0_1(0) * occ_func_0_1(3)) +
          (occ_func_0_1(0) * occ_func_0_1(1)) +
          (occ_func_0_1(0) * occ_func_0_1(11))) /
         6.0;
}

double test_Clexulator::site_eval_at_0_bfunc_2_0_0() const {
  return ((occ_func_0_0(0) * occ_func_0_0(6)) +
          (occ_func_0_0(7) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(4)) +
          (occ_func_0_0(9) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(8)) +
          (occ_func_0_0(5) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(3)) +
          (occ_func_0_0(10) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(1)) +
          (occ_func_0_0(12) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(11)) +
          (occ_func_0_0(2) * occ_func_0_0(0))) /
         6.0;
}
double test_Clexulator::site_eval_at_0_bfunc_2_0_1() const {
  return (((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(6) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(6))) +
          ((0.7071067812 * occ_func_0_1(7) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(7) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(4) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(4))) +
          ((0.7071067812 * occ_func_0_1(9) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(9) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(8) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(8))) +
          ((0.7071067812 * occ_func_0_1(5) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(5) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(3) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(3))) +
          ((0.7071067812 * occ_func_0_1(10) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(10) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(1) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(1))) +
          ((0.7071067812 * occ_func_0_1(12) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(12) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(11) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(11))) +
          ((0.7071067812 * occ_func_0_1(2) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(2) * occ_func_0_1(0)))) /
         6.0;
}
double test_Clexulator::site_eval_at_0_bfunc_2_0_2() const {
  return ((occ_func_0_1(0) * occ_func_0_1(6)) +
          (occ_func_0_1(7) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(4)) +
          (occ_func_0_1(9) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(8)) +
          (occ_func_0_1(5) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(3)) +
          (occ_func_0_1(10) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(1)) +
          (occ_func_0_1(12) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(11)) +
          (occ_func_0_1(2) * occ_func_0_1(0))) /
         6.0;
}

double test_Clexulator::delta_site_eval_at_0_bfunc_2_0_0(int occ_i,
                                                         int occ_f) const {
  return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) *
         ((occ_func_0_0(6)) + (occ_func_0_0(7)) + (occ_func_0_0(4)) +
          (occ_func_0_0(9)) + (occ_func_0_0(8)) + (occ_func_0_0(5)) +
          (occ_func_0_0(3)) + (occ_func_0_0(10)) + (occ_func_0_0(1)) +
          (occ_func_0_0(12)) + (occ_func_0_0(11)) + (occ_func_0_0(2))) /
         6.0;
}
double test_Clexulator::delta_site_eval_at_0_bfunc_2_0_1(int occ_i,
                                                         int occ_f) const {
  return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) *
             ((0.7071067812 * occ_func_0_1(6)) +
              (0.7071067812 * occ_func_0_1(7)) +
              (0.7071067812 * occ_func_0_1(4)) +
              (0.7071067812 * occ_func_0_1(9)) +
              (0.7071067812 * occ_func_0_1(8)) +
              (0.7071067812 * occ_func_0_1(5)) +
              (0.7071067812 * occ_func_0_1(3)) +
              (0.7071067812 * occ_func_0_1(10)) +
              (0.7071067812 * occ_func_0_1(1)) +
              (0.7071067812 * occ_func_0_1(12)) +
              (0.7071067812 * occ_func_0_1(11)) +
              (0.7071067812 * occ_func_0_1(2))) /
             6.0 +
         (m_occ_func_0_1[occ_f] - m_occ_func_0_1[occ_i]) *
             ((0.7071067812 * occ_func_0_0(6)) +
              (0.7071067812 * occ_func_0_0(7)) +
              (0.7071067812 * occ_func_0_0(4)) +
              (0.7071067812 * occ_func_0_0(9)) +
              (0.7071067812 * occ_func_0_0(8)) +
              (0.7071067812 * occ_func_0_0(5)) +
              (0.7071067812 * occ_func_0_0(3)) +
              (0.7071067812 * occ_func_0_0(10)) +
              (0.7071067812 * occ_func_0_0(1)) +
              (0.7071067812 * occ_func_0_0(12)) +
              (0.7071067812 * occ_func_0_0(11)) +
              (0.7071067812 * occ_func_0_0(2))) /
             6.0;
}
double test_Clexulator::delta_site_eval_at_0_bfunc_2_0_2(int occ_i,
                                                         int occ_f) const {
  return (m_occ_func_0_1[occ_f] - m_occ_func_0_1[occ_i]) *
         ((occ_func_0_1(6)) + (occ_func_0_1(7)) + (occ_func_0_1(4)) +
          (occ_func_0_1(9)) + (occ_func_0_1(8)) + (occ_func_0_1(5)) +
          (occ_func_0_1(3)) + (occ_func_0_1(10)) + (occ_func_0_1(1)) +
          (occ_func_0_1(12)) + (occ_func_0_1(11)) + (occ_func_0_1(2))) /
         6.0;
}

/**** Basis functions for orbit 2, 1****
#Points: 2
MaxLength: 4.0000000  MinLength: 4.0000000
 0.0000000   0.0000000   0.0000000 A B C
 1.0000000  -1.0000000  -1.0000000 A B C
****/
double test_Clexulator::eval_bfunc_2_1_0() const {
  return ((occ_func_0_0(0) * occ_func_0_0(16)) +
          (occ_func_0_0(0) * occ_func_0_0(18)) +
          (occ_func_0_0(0) * occ_func_0_0(14))) /
         3.0;
}
double test_Clexulator::eval_bfunc_2_1_1() const {
  return (((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(16) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(16))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(18) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(18))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(14) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(14)))) /
         3.0;
}
double test_Clexulator::eval_bfunc_2_1_2() const {
  return ((occ_func_0_1(0) * occ_func_0_1(16)) +
          (occ_func_0_1(0) * occ_func_0_1(18)) +
          (occ_func_0_1(0) * occ_func_0_1(14))) /
         3.0;
}

double test_Clexulator::site_eval_at_0_bfunc_2_1_0() const {
  return ((occ_func_0_0(0) * occ_func_0_0(16)) +
          (occ_func_0_0(15) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(18)) +
          (occ_func_0_0(13) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(14)) +
          (occ_func_0_0(17) * occ_func_0_0(0))) /
         3.0;
}
double test_Clexulator::site_eval_at_0_bfunc_2_1_1() const {
  return (((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(16) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(16))) +
          ((0.7071067812 * occ_func_0_1(15) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(15) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(18) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(18))) +
          ((0.7071067812 * occ_func_0_1(13) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(13) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(14) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(14))) +
          ((0.7071067812 * occ_func_0_1(17) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(17) * occ_func_0_1(0)))) /
         3.0;
}
double test_Clexulator::site_eval_at_0_bfunc_2_1_2() const {
  return ((occ_func_0_1(0) * occ_func_0_1(16)) +
          (occ_func_0_1(15) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(18)) +
          (occ_func_0_1(13) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(14)) +
          (occ_func_0_1(17) * occ_func_0_1(0))) /
         3.0;
}

double test_Clexulator::delta_site_eval_at_0_bfunc_2_1_0(int occ_i,
                                                         int occ_f) const {
  return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) *
         ((occ_func_0_0(16)) + (occ_func_0_0(15)) + (occ_func_0_0(18)) +
          (occ_func_0_0(13)) + (occ_func_0_0(14)) + (occ_func_0_0(17))) /
         3.0;
}
double test_Clexulator::delta_site_eval_at_0_bfunc_2_1_1(int occ_i,
                                                         int occ_f) const {
  return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) *
             ((0.7071067812 * occ_func_0_1(16)) +
              (0.7071067812 * occ_func_0_1(15)) +
              (0.7071067812 * occ_func_0_1(18)) +
              (0.7071067812 * occ_func_0_1(13)) +
              (0.7071067812 * occ_func_0_1(14)) +
              (0.7071067812 * occ_func_0_1(17))) /
             3.0 +
         (m_occ_func_0_1[occ_f] - m_occ_func_0_1[occ_i]) *
             ((0.7071067812 * occ_func_0_0(16)) +
              (0.7071067812 * occ_func_0_0(15)) +
              (0.7071067812 * occ_func_0_0(18)) +
              (0.7071067812 * occ_func_0_0(13)) +
              (0.7071067812 * occ_func_0_0(14)) +
              (0.7071067812 * occ_func_0_0(17))) /
             3.0;
}
double test_Clexulator::delta_site_eval_at_0_bfunc_2_1_2(int occ_i,
                                                         int occ_f) const {
  return (m_occ_func_0_1[occ_f] - m_occ_func_0_1[occ_i]) *
         ((occ_func_0_1(16)) + (occ_func_0_1(15)) + (occ_func_0_1(18)) +
          (occ_func_0_1(13)) + (occ_func_0_1(14)) + (occ_func_0_1(17))) /
         3.0;
}

/**** Basis functions for orbit 2, 2****
#Points: 2
MaxLength: 4.8989795  MinLength: 4.8989795
 0.0000000   0.0000000   0.0000000 A B C
 1.0000000   1.0000000   0.0000000 A B C
****/
double test_Clexulator::eval_bfunc_2_2_0() const {
  return ((occ_func_0_0(0) * occ_func_0_0(39)) +
          (occ_func_0_0(0) * occ_func_0_0(26)) +
          (occ_func_0_0(0) * occ_func_0_0(25)) +
          (occ_func_0_0(0) * occ_func_0_0(42)) +
          (occ_func_0_0(0) * occ_func_0_0(33)) +
          (occ_func_0_0(0) * occ_func_0_0(40)) +
          (occ_func_0_0(0) * occ_func_0_0(30)) +
          (occ_func_0_0(0) * occ_func_0_0(37)) +
          (occ_func_0_0(0) * occ_func_0_0(32)) +
          (occ_func_0_0(0) * occ_func_0_0(20)) +
          (occ_func_0_0(0) * occ_func_0_0(34)) +
          (occ_func_0_0(0) * occ_func_0_0(38))) /
         12.0;
}
double test_Clexulator::eval_bfunc_2_2_1() const {
  return (((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(39) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(39))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(26) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(26))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(25) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(25))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(42) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(42))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(33) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(33))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(40) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(40))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(30) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(30))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(37) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(37))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(32) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(32))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(20) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(20))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(34) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(34))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(38) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(38)))) /
         12.0;
}
double test_Clexulator::eval_bfunc_2_2_2() const {
  return ((occ_func_0_1(0) * occ_func_0_1(39)) +
          (occ_func_0_1(0) * occ_func_0_1(26)) +
          (occ_func_0_1(0) * occ_func_0_1(25)) +
          (occ_func_0_1(0) * occ_func_0_1(42)) +
          (occ_func_0_1(0) * occ_func_0_1(33)) +
          (occ_func_0_1(0) * occ_func_0_1(40)) +
          (occ_func_0_1(0) * occ_func_0_1(30)) +
          (occ_func_0_1(0) * occ_func_0_1(37)) +
          (occ_func_0_1(0) * occ_func_0_1(32)) +
          (occ_func_0_1(0) * occ_func_0_1(20)) +
          (occ_func_0_1(0) * occ_func_0_1(34)) +
          (occ_func_0_1(0) * occ_func_0_1(38))) /
         12.0;
}

double test_Clexulator::site_eval_at_0_bfunc_2_2_0() const {
  return ((occ_func_0_0(0) * occ_func_0_0(39)) +
          (occ_func_0_0(22) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(26)) +
          (occ_func_0_0(35) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(25)) +
          (occ_func_0_0(36) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(42)) +
          (occ_func_0_0(19) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(33)) +
          (occ_func_0_0(28) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(40)) +
          (occ_func_0_0(21) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(30)) +
          (occ_func_0_0(31) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(37)) +
          (occ_func_0_0(24) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(32)) +
          (occ_func_0_0(29) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(20)) +
          (occ_func_0_0(41) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(34)) +
          (occ_func_0_0(27) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(38)) +
          (occ_func_0_0(23) * occ_func_0_0(0))) /
         12.0;
}
double test_Clexulator::site_eval_at_0_bfunc_2_2_1() const {
  return (((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(39) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(39))) +
          ((0.7071067812 * occ_func_0_1(22) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(22) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(26) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(26))) +
          ((0.7071067812 * occ_func_0_1(35) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(35) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(25) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(25))) +
          ((0.7071067812 * occ_func_0_1(36) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(36) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(42) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(42))) +
          ((0.7071067812 * occ_func_0_1(19) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(19) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(33) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(33))) +
          ((0.7071067812 * occ_func_0_1(28) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(28) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(40) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(40))) +
          ((0.7071067812 * occ_func_0_1(21) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(21) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(30) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(30))) +
          ((0.7071067812 * occ_func_0_1(31) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(31) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(37) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(37))) +
          ((0.7071067812 * occ_func_0_1(24) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(24) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(32) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(32))) +
          ((0.7071067812 * occ_func_0_1(29) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(29) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(20) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(20))) +
          ((0.7071067812 * occ_func_0_1(41) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(41) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(34) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(34))) +
          ((0.7071067812 * occ_func_0_1(27) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(27) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(38) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(38))) +
          ((0.7071067812 * occ_func_0_1(23) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(23) * occ_func_0_1(0)))) /
         12.0;
}
double test_Clexulator::site_eval_at_0_bfunc_2_2_2() const {
  return ((occ_func_0_1(0) * occ_func_0_1(39)) +
          (occ_func_0_1(22) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(26)) +
          (occ_func_0_1(35) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(25)) +
          (occ_func_0_1(36) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(42)) +
          (occ_func_0_1(19) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(33)) +
          (occ_func_0_1(28) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(40)) +
          (occ_func_0_1(21) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(30)) +
          (occ_func_0_1(31) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(37)) +
          (occ_func_0_1(24) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(32)) +
          (occ_func_0_1(29) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(20)) +
          (occ_func_0_1(41) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(34)) +
          (occ_func_0_1(27) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(38)) +
          (occ_func_0_1(23) * occ_func_0_1(0))) /
         12.0;
}

double test_Clexulator::delta_site_eval_at_0_bfunc_2_2_0(int occ_i,
                                                         int occ_f) const {
  return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) *
         ((occ_func_0_0(39)) + (occ_func_0_0(22)) + (occ_func_0_0(26)) +
          (occ_func_0_0(35)) + (occ_func_0_0(25)) + (occ_func_0_0(36)) +
          (occ_func_0_0(42)) + (occ_func_0_0(19)) + (occ_func_0_0(33)) +
          (occ_func_0_0(28)) + (occ_func_0_0(40)) + (occ_func_0_0(21)) +
          (occ_func_0_0(30)) + (occ_func_0_0(31)) + (occ_func_0_0(37)) +
          (occ_func_0_0(24)) + (occ_func_0_0(32)) + (occ_func_0_0(29)) +
          (occ_func_0_0(20)) + (occ_func_0_0(41)) + (occ_func_0_0(34)) +
          (occ_func_0_0(27)) + (occ_func_0_0(38)) + (occ_func_0_0(23))) /
         12.0;
}
double test_Clexulator::delta_site_eval_at_0_bfunc_2_2_1(int occ_i,
                                                         int occ_f) const {
  return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) *
             ((0.7071067812 * occ_func_0_1(39)) +
              (0.7071067812 * occ_func_0_1(22)) +
              (0.7071067812 * occ_func_0_1(26)) +
              (0.7071067812 * occ_func_0_1(35)) +
              (0.7071067812 * occ_func_0_1(25)) +
              (0.7071067812 * occ_func_0_1(36)) +
              (0.7071067812 * occ_func_0_1(42)) +
              (0.7071067812 * occ_func_0_1(19)) +
              (0.7071067812 * occ_func_0_1(33)) +
              (0.7071067812 * occ_func_0_1(28)) +
              (0.7071067812 * occ_func_0_1(40)) +
              (0.7071067812 * occ_func_0_1(21)) +
              (0.7071067812 * occ_func_0_1(30)) +
              (0.7071067812 * occ_func_0_1(31)) +
              (0.7071067812 * occ_func_0_1(37)) +
              (0.7071067812 * occ_func_0_1(24)) +
              (0.7071067812 * occ_func_0_1(32)) +
              (0.7071067812 * occ_func_0_1(29)) +
              (0.7071067812 * occ_func_0_1(20)) +
              (0.7071067812 * occ_func_0_1(41)) +
              (0.7071067812 * occ_func_0_1(34)) +
              (0.7071067812 * occ_func_0_1(27)) +
              (0.7071067812 * occ_func_0_1(38)) +
              (0.7071067812 * occ_func_0_1(23))) /
             12.0 +
         (m_occ_func_0_1[occ_f] - m_occ_func_0_1[occ_i]) *
             ((0.7071067812 * occ_func_0_0(39)) +
              (0.7071067812 * occ_func_0_0(22)) +
              (0.7071067812 * occ_func_0_0(26)) +
              (0.7071067812 * occ_func_0_0(35)) +
              (0.7071067812 * occ_func_0_0(25)) +
              (0.7071067812 * occ_func_0_0(36)) +
              (0.7071067812 * occ_func_0_0(42)) +
              (0.7071067812 * occ_func_0_0(19)) +
              (0.7071067812 * occ_func_0_0(33)) +
              (0.7071067812 * occ_func_0_0(28)) +
              (0.7071067812 * occ_func_0_0(40)) +
              (0.7071067812 * occ_func_0_0(21)) +
              (0.7071067812 * occ_func_0_0(30)) +
              (0.7071067812 * occ_func_0_0(31)) +
              (0.7071067812 * occ_func_0_0(37)) +
              (0.7071067812 * occ_func_0_0(24)) +
              (0.7071067812 * occ_func_0_0(32)) +
              (0.7071067812 * occ_func_0_0(29)) +
              (0.7071067812 * occ_func_0_0(20)) +
              (0.7071067812 * occ_func_0_0(41)) +
              (0.7071067812 * occ_func_0_0(34)) +
              (0.7071067812 * occ_func_0_0(27)) +
              (0.7071067812 * occ_func_0_0(38)) +
              (0.7071067812 * occ_func_0_0(23))) /
             12.0;
}
double test_Clexulator::delta_site_eval_at_0_bfunc_2_2_2(int occ_i,
                                                         int occ_f) const {
  return (m_occ_func_0_1[occ_f] - m_occ_func_0_1[occ_i]) *
         ((occ_func_0_1(39)) + (occ_func_0_1(22)) + (occ_func_0_1(26)) +
          (occ_func_0_1(35)) + (occ_func_0_1(25)) + (occ_func_0_1(36)) +
          (occ_func_0_1(42)) + (occ_func_0_1(19)) + (occ_func_0_1(33)) +
          (occ_func_0_1(28)) + (occ_func_0_1(40)) + (occ_func_0_1(21)) +
          (occ_func_0_1(30)) + (occ_func_0_1(31)) + (occ_func_0_1(37)) +
          (occ_func_0_1(24)) + (occ_func_0_1(32)) + (occ_func_0_1(29)) +
          (occ_func_0_1(20)) + (occ_func_0_1(41)) + (occ_func_0_1(34)) +
          (occ_func_0_1(27)) + (occ_func_0_1(38)) + (occ_func_0_1(23))) /
         12.0;
}

/**** Basis functions for orbit 2, 3****
#Points: 2
MaxLength: 5.6568542  MinLength: 5.6568542
 0.0000000   0.0000000   0.0000000 A B C
 2.0000000   0.0000000   0.0000000 A B C
****/
double test_Clexulator::eval_bfunc_2_3_0() const {
  return ((occ_func_0_0(0) * occ_func_0_0(54)) +
          (occ_func_0_0(0) * occ_func_0_0(50)) +
          (occ_func_0_0(0) * occ_func_0_0(49)) +
          (occ_func_0_0(0) * occ_func_0_0(53)) +
          (occ_func_0_0(0) * occ_func_0_0(51)) +
          (occ_func_0_0(0) * occ_func_0_0(52))) /
         6.0;
}
double test_Clexulator::eval_bfunc_2_3_1() const {
  return (((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(54) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(54))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(50) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(50))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(49) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(49))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(53) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(53))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(51) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(51))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(52) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(52)))) /
         6.0;
}
double test_Clexulator::eval_bfunc_2_3_2() const {
  return ((occ_func_0_1(0) * occ_func_0_1(54)) +
          (occ_func_0_1(0) * occ_func_0_1(50)) +
          (occ_func_0_1(0) * occ_func_0_1(49)) +
          (occ_func_0_1(0) * occ_func_0_1(53)) +
          (occ_func_0_1(0) * occ_func_0_1(51)) +
          (occ_func_0_1(0) * occ_func_0_1(52))) /
         6.0;
}

double test_Clexulator::site_eval_at_0_bfunc_2_3_0() const {
  return ((occ_func_0_0(0) * occ_func_0_0(54)) +
          (occ_func_0_0(43) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(50)) +
          (occ_func_0_0(47) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(49)) +
          (occ_func_0_0(48) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(53)) +
          (occ_func_0_0(44) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(51)) +
          (occ_func_0_0(46) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(52)) +
          (occ_func_0_0(45) * occ_func_0_0(0))) /
         6.0;
}
double test_Clexulator::site_eval_at_0_bfunc_2_3_1() const {
  return (((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(54) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(54))) +
          ((0.7071067812 * occ_func_0_1(43) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(43) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(50) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(50))) +
          ((0.7071067812 * occ_func_0_1(47) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(47) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(49) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(49))) +
          ((0.7071067812 * occ_func_0_1(48) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(48) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(53) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(53))) +
          ((0.7071067812 * occ_func_0_1(44) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(44) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(51) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(51))) +
          ((0.7071067812 * occ_func_0_1(46) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(46) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(52) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(52))) +
          ((0.7071067812 * occ_func_0_1(45) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(45) * occ_func_0_1(0)))) /
         6.0;
}
double test_Clexulator::site_eval_at_0_bfunc_2_3_2() const {
  return ((occ_func_0_1(0) * occ_func_0_1(54)) +
          (occ_func_0_1(43) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(50)) +
          (occ_func_0_1(47) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(49)) +
          (occ_func_0_1(48) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(53)) +
          (occ_func_0_1(44) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(51)) +
          (occ_func_0_1(46) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(52)) +
          (occ_func_0_1(45) * occ_func_0_1(0))) /
         6.0;
}

double test_Clexulator::delta_site_eval_at_0_bfunc_2_3_0(int occ_i,
                                                         int occ_f) const {
  return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) *
         ((occ_func_0_0(54)) + (occ_func_0_0(43)) + (occ_func_0_0(50)) +
          (occ_func_0_0(47)) + (occ_func_0_0(49)) + (occ_func_0_0(48)) +
          (occ_func_0_0(53)) + (occ_func_0_0(44)) + (occ_func_0_0(51)) +
          (occ_func_0_0(46)) + (occ_func_0_0(52)) + (occ_func_0_0(45))) /
         6.0;
}
double test_Clexulator::delta_site_eval_at_0_bfunc_2_3_1(int occ_i,
                                                         int occ_f) const {
  return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) *
             ((0.7071067812 * occ_func_0_1(54)) +
              (0.7071067812 * occ_func_0_1(43)) +
              (0.7071067812 * occ_func_0_1(50)) +
              (0.7071067812 * occ_func_0_1(47)) +
              (0.7071067812 * occ_func_0_1(49)) +
              (0.7071067812 * occ_func_0_1(48)) +
              (0.7071067812 * occ_func_0_1(53)) +
              (0.7071067812 * occ_func_0_1(44)) +
              (0.7071067812 * occ_func_0_1(51)) +
              (0.7071067812 * occ_func_0_1(46)) +
              (0.7071067812 * occ_func_0_1(52)) +
              (0.7071067812 * occ_func_0_1(45))) /
             6.0 +
         (m_occ_func_0_1[occ_f] - m_occ_func_0_1[occ_i]) *
             ((0.7071067812 * occ_func_0_0(54)) +
              (0.7071067812 * occ_func_0_0(43)) +
              (0.7071067812 * occ_func_0_0(50)) +
              (0.7071067812 * occ_func_0_0(47)) +
              (0.7071067812 * occ_func_0_0(49)) +
              (0.7071067812 * occ_func_0_0(48)) +
              (0.7071067812 * occ_func_0_0(53)) +
              (0.7071067812 * occ_func_0_0(44)) +
              (0.7071067812 * occ_func_0_0(51)) +
              (0.7071067812 * occ_func_0_0(46)) +
              (0.7071067812 * occ_func_0_0(52)) +
              (0.7071067812 * occ_func_0_0(45))) /
             6.0;
}
double test_Clexulator::delta_site_eval_at_0_bfunc_2_3_2(int occ_i,
                                                         int occ_f) const {
  return (m_occ_func_0_1[occ_f] - m_occ_func_0_1[occ_i]) *
         ((occ_func_0_1(54)) + (occ_func_0_1(43)) + (occ_func_0_1(50)) +
          (occ_func_0_1(47)) + (occ_func_0_1(49)) + (occ_func_0_1(48)) +
          (occ_func_0_1(53)) + (occ_func_0_1(44)) + (occ_func_0_1(51)) +
          (occ_func_0_1(46)) + (occ_func_0_1(52)) + (occ_func_0_1(45))) /
         6.0;
}

/**** Basis functions for orbit 2, 4****
#Points: 2
MaxLength: 6.9282032  MinLength: 6.9282032
 0.0000000   0.0000000   0.0000000 A B C
 1.0000000   1.0000000   1.0000000 A B C
****/
double test_Clexulator::eval_bfunc_2_4_0() const {
  return ((occ_func_0_0(0) * occ_func_0_0(85)) +
          (occ_func_0_0(0) * occ_func_0_0(82)) +
          (occ_func_0_0(0) * occ_func_0_0(81)) +
          (occ_func_0_0(0) * occ_func_0_0(86))) /
         4.0;
}
double test_Clexulator::eval_bfunc_2_4_1() const {
  return (((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(85) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(85))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(82) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(82))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(81) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(81))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(86) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(86)))) /
         4.0;
}
double test_Clexulator::eval_bfunc_2_4_2() const {
  return ((occ_func_0_1(0) * occ_func_0_1(85)) +
          (occ_func_0_1(0) * occ_func_0_1(82)) +
          (occ_func_0_1(0) * occ_func_0_1(81)) +
          (occ_func_0_1(0) * occ_func_0_1(86))) /
         4.0;
}

double test_Clexulator::site_eval_at_0_bfunc_2_4_0() const {
  return ((occ_func_0_0(0) * occ_func_0_0(85)) +
          (occ_func_0_0(80) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(82)) +
          (occ_func_0_0(83) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(81)) +
          (occ_func_0_0(84) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(86)) +
          (occ_func_0_0(79) * occ_func_0_0(0))) /
         4.0;
}
double test_Clexulator::site_eval_at_0_bfunc_2_4_1() const {
  return (((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(85) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(85))) +
          ((0.7071067812 * occ_func_0_1(80) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(80) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(82) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(82))) +
          ((0.7071067812 * occ_func_0_1(83) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(83) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(81) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(81))) +
          ((0.7071067812 * occ_func_0_1(84) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(84) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(86) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(86))) +
          ((0.7071067812 * occ_func_0_1(79) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(79) * occ_func_0_1(0)))) /
         4.0;
}
double test_Clexulator::site_eval_at_0_bfunc_2_4_2() const {
  return ((occ_func_0_1(0) * occ_func_0_1(85)) +
          (occ_func_0_1(80) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(82)) +
          (occ_func_0_1(83) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(81)) +
          (occ_func_0_1(84) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(86)) +
          (occ_func_0_1(79) * occ_func_0_1(0))) /
         4.0;
}

double test_Clexulator::delta_site_eval_at_0_bfunc_2_4_0(int occ_i,
                                                         int occ_f) const {
  return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) *
         ((occ_func_0_0(85)) + (occ_func_0_0(80)) + (occ_func_0_0(82)) +
          (occ_func_0_0(83)) + (occ_func_0_0(81)) + (occ_func_0_0(84)) +
          (occ_func_0_0(86)) + (occ_func_0_0(79))) /
         4.0;
}
double test_Clexulator::delta_site_eval_at_0_bfunc_2_4_1(int occ_i,
                                                         int occ_f) const {
  return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) *
             ((0.7071067812 * occ_func_0_1(85)) +
              (0.7071067812 * occ_func_0_1(80)) +
              (0.7071067812 * occ_func_0_1(82)) +
              (0.7071067812 * occ_func_0_1(83)) +
              (0.7071067812 * occ_func_0_1(81)) +
              (0.7071067812 * occ_func_0_1(84)) +
              (0.7071067812 * occ_func_0_1(86)) +
              (0.7071067812 * occ_func_0_1(79))) /
             4.0 +
         (m_occ_func_0_1[occ_f] - m_occ_func_0_1[occ_i]) *
             ((0.7071067812 * occ_func_0_0(85)) +
              (0.7071067812 * occ_func_0_0(80)) +
              (0.7071067812 * occ_func_0_0(82)) +
              (0.7071067812 * occ_func_0_0(83)) +
              (0.7071067812 * occ_func_0_0(81)) +
              (0.7071067812 * occ_func_0_0(84)) +
              (0.7071067812 * occ_func_0_0(86)) +
              (0.7071067812 * occ_func_0_0(79))) /
             4.0;
}
double test_Clexulator::delta_site_eval_at_0_bfunc_2_4_2(int occ_i,
                                                         int occ_f) const {
  return (m_occ_func_0_1[occ_f] - m_occ_func_0_1[occ_i]) *
         ((occ_func_0_1(85)) + (occ_func_0_1(80)) + (occ_func_0_1(82)) +
          (occ_func_0_1(83)) + (occ_func_0_1(81)) + (occ_func_0_1(84)) +
          (occ_func_0_1(86)) + (occ_func_0_1(79))) /
         4.0;
}

/**** Basis functions for orbit 2, 5****
#Points: 2
MaxLength: 8.4852814  MinLength: 8.4852814
 0.0000000   0.0000000   0.0000000 A B C
 3.0000000   0.0000000   0.0000000 A B C
****/
double test_Clexulator::eval_bfunc_2_5_0() const {
  return ((occ_func_0_0(0) * occ_func_0_0(175)) +
          (occ_func_0_0(0) * occ_func_0_0(160)) +
          (occ_func_0_0(0) * occ_func_0_0(159)) +
          (occ_func_0_0(0) * occ_func_0_0(174)) +
          (occ_func_0_0(0) * occ_func_0_0(161)) +
          (occ_func_0_0(0) * occ_func_0_0(171))) /
         6.0;
}
double test_Clexulator::eval_bfunc_2_5_1() const {
  return (((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(175) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(175))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(160) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(160))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(159) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(159))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(174) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(174))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(161) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(161))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(171) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(171)))) /
         6.0;
}
double test_Clexulator::eval_bfunc_2_5_2() const {
  return ((occ_func_0_1(0) * occ_func_0_1(175)) +
          (occ_func_0_1(0) * occ_func_0_1(160)) +
          (occ_func_0_1(0) * occ_func_0_1(159)) +
          (occ_func_0_1(0) * occ_func_0_1(174)) +
          (occ_func_0_1(0) * occ_func_0_1(161)) +
          (occ_func_0_1(0) * occ_func_0_1(171))) /
         6.0;
}

double test_Clexulator::site_eval_at_0_bfunc_2_5_0() const {
  return ((occ_func_0_0(0) * occ_func_0_0(175)) +
          (occ_func_0_0(142) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(160)) +
          (occ_func_0_0(157) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(159)) +
          (occ_func_0_0(158) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(174)) +
          (occ_func_0_0(143) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(161)) +
          (occ_func_0_0(156) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(171)) +
          (occ_func_0_0(146) * occ_func_0_0(0))) /
         6.0;
}
double test_Clexulator::site_eval_at_0_bfunc_2_5_1() const {
  return (((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(175) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(175))) +
          ((0.7071067812 * occ_func_0_1(142) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(142) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(160) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(160))) +
          ((0.7071067812 * occ_func_0_1(157) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(157) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(159) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(159))) +
          ((0.7071067812 * occ_func_0_1(158) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(158) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(174) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(174))) +
          ((0.7071067812 * occ_func_0_1(143) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(143) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(161) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(161))) +
          ((0.7071067812 * occ_func_0_1(156) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(156) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(171) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(171))) +
          ((0.7071067812 * occ_func_0_1(146) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(146) * occ_func_0_1(0)))) /
         6.0;
}
double test_Clexulator::site_eval_at_0_bfunc_2_5_2() const {
  return ((occ_func_0_1(0) * occ_func_0_1(175)) +
          (occ_func_0_1(142) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(160)) +
          (occ_func_0_1(157) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(159)) +
          (occ_func_0_1(158) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(174)) +
          (occ_func_0_1(143) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(161)) +
          (occ_func_0_1(156) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(171)) +
          (occ_func_0_1(146) * occ_func_0_1(0))) /
         6.0;
}

double test_Clexulator::delta_site_eval_at_0_bfunc_2_5_0(int occ_i,
                                                         int occ_f) const {
  return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) *
         ((occ_func_0_0(175)) + (occ_func_0_0(142)) + (occ_func_0_0(160)) +
          (occ_func_0_0(157)) + (occ_func_0_0(159)) + (occ_func_0_0(158)) +
          (occ_func_0_0(174)) + (occ_func_0_0(143)) + (occ_func_0_0(161)) +
          (occ_func_0_0(156)) + (occ_func_0_0(171)) + (occ_func_0_0(146))) /
         6.0;
}
double test_Clexulator::delta_site_eval_at_0_bfunc_2_5_1(int occ_i,
                                                         int occ_f) const {
  return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) *
             ((0.7071067812 * occ_func_0_1(175)) +
              (0.7071067812 * occ_func_0_1(142)) +
              (0.7071067812 * occ_func_0_1(160)) +
              (0.7071067812 * occ_func_0_1(157)) +
              (0.7071067812 * occ_func_0_1(159)) +
              (0.7071067812 * occ_func_0_1(158)) +
              (0.7071067812 * occ_func_0_1(174)) +
              (0.7071067812 * occ_func_0_1(143)) +
              (0.7071067812 * occ_func_0_1(161)) +
              (0.7071067812 * occ_func_0_1(156)) +
              (0.7071067812 * occ_func_0_1(171)) +
              (0.7071067812 * occ_func_0_1(146))) /
             6.0 +
         (m_occ_func_0_1[occ_f] - m_occ_func_0_1[occ_i]) *
             ((0.7071067812 * occ_func_0_0(175)) +
              (0.7071067812 * occ_func_0_0(142)) +
              (0.7071067812 * occ_func_0_0(160)) +
              (0.7071067812 * occ_func_0_0(157)) +
              (0.7071067812 * occ_func_0_0(159)) +
              (0.7071067812 * occ_func_0_0(158)) +
              (0.7071067812 * occ_func_0_0(174)) +
              (0.7071067812 * occ_func_0_0(143)) +
              (0.7071067812 * occ_func_0_0(161)) +
              (0.7071067812 * occ_func_0_0(156)) +
              (0.7071067812 * occ_func_0_0(171)) +
              (0.7071067812 * occ_func_0_0(146))) /
             6.0;
}
double test_Clexulator::delta_site_eval_at_0_bfunc_2_5_2(int occ_i,
                                                         int occ_f) const {
  return (m_occ_func_0_1[occ_f] - m_occ_func_0_1[occ_i]) *
         ((occ_func_0_1(175)) + (occ_func_0_1(142)) + (occ_func_0_1(160)) +
          (occ_func_0_1(157)) + (occ_func_0_1(159)) + (occ_func_0_1(158)) +
          (occ_func_0_1(174)) + (occ_func_0_1(143)) + (occ_func_0_1(161)) +
          (occ_func_0_1(156)) + (occ_func_0_1(171)) + (occ_func_0_1(146))) /
         6.0;
}

/**** Basis functions for orbit 3, 0****
#Points: 3
MaxLength: 2.8284271  MinLength: 2.8284271
 0.0000000   0.0000000   0.0000000 A B C
 0.0000000   0.0000000  -1.0000000 A B C
 1.0000000   0.0000000  -1.0000000 A B C
****/
double test_Clexulator::eval_bfunc_3_0_0() const {
  return ((occ_func_0_0(0) * occ_func_0_0(6) * occ_func_0_0(11)) +
          (occ_func_0_0(0) * occ_func_0_0(4) * occ_func_0_0(6)) +
          (occ_func_0_0(0) * occ_func_0_0(3) * occ_func_0_0(8)) +
          (occ_func_0_0(0) * occ_func_0_0(10) * occ_func_0_0(12)) +
          (occ_func_0_0(0) * occ_func_0_0(5) * occ_func_0_0(10)) +
          (occ_func_0_0(0) * occ_func_0_0(4) * occ_func_0_0(5)) +
          (occ_func_0_0(0) * occ_func_0_0(9) * occ_func_0_0(12)) +
          (occ_func_0_0(0) * occ_func_0_0(2) * occ_func_0_0(7))) /
         8.0;
}
double test_Clexulator::eval_bfunc_3_0_1() const {
  return (((0.5773502692 * occ_func_0_1(0) * occ_func_0_0(6) *
                occ_func_0_0(11) +
            0.5773502692 * occ_func_0_0(0) * occ_func_0_1(6) *
                occ_func_0_0(11) +
            0.5773502692 * occ_func_0_0(0) * occ_func_0_0(6) *
                occ_func_0_1(11))) +
          ((0.5773502692 * occ_func_0_1(0) * occ_func_0_0(4) * occ_func_0_0(6) +
            0.5773502692 * occ_func_0_0(0) * occ_func_0_1(4) * occ_func_0_0(6) +
            0.5773502692 * occ_func_0_0(0) * occ_func_0_0(4) *
                occ_func_0_1(6))) +
          ((0.5773502692 * occ_func_0_1(0) * occ_func_0_0(3) * occ_func_0_0(8) +
            0.5773502692 * occ_func_0_0(0) * occ_func_0_1(3) * occ_func_0_0(8) +
            0.5773502692 * occ_func_0_0(0) * occ_func_0_0(3) *
                occ_func_0_1(8))) +
          ((0.5773502692 * occ_func_0_1(0) * occ_func_0_0(10) *
                occ_func_0_0(12) +
            0.5773502692 * occ_func_0_0(0) * occ_func_0_1(10) *
                occ_func_0_0(12) +
            0.5773502692 * occ_func_0_0(0) * occ_func_0_0(10) *
                occ_func_0_1(12))) +
          ((0.5773502692 * occ_func_0_1(0) * occ_func_0_0(5) *
                occ_func_0_0(10) +
            0.5773502692 * occ_func_0_0(0) * occ_func_0_1(5) *
                occ_func_0_0(10) +
            0.5773502692 * occ_func_0_0(0) * occ_func_0_0(5) *
                occ_func_0_1(10))) +
          ((0.5773502692 * occ_func_0_1(0) * occ_func_0_0(4) * occ_func_0_0(5) +
            0.5773502692 * occ_func_0_0(0) * occ_func_0_1(4) * occ_func_0_0(5) +
            0.5773502692 * occ_func_0_0(0) * occ_func_0_0(4) *
                occ_func_0_1(5))) +
          ((0.5773502692 * occ_func_0_1(0) * occ_func_0_0(9) *
                occ_func_0_0(12) +
            0.5773502692 * occ_func_0_0(0) * occ_func_0_1(9) *
                occ_func_0_0(12) +
            0.5773502692 * occ_func_0_0(0) * occ_func_0_0(9) *
                occ_func_0_1(12))) +
          ((0.5773502692 * occ_func_0_1(0) * occ_func_0_0(2) * occ_func_0_0(7) +
            0.5773502692 * occ_func_0_0(0) * occ_func_0_1(2) * occ_func_0_0(7) +
            0.5773502692 * occ_func_0_0(0) * occ_func_0_0(2) *
                occ_func_0_1(7)))) /
         8.0;
}
double test_Clexulator::eval_bfunc_3_0_2() const {
  return (((0.5773502692 * occ_func_0_1(0) * occ_func_0_1(6) *
                occ_func_0_0(11) +
            0.5773502692 * occ_func_0_1(0) * occ_func_0_0(6) *
                occ_func_0_1(11) +
            0.5773502692 * occ_func_0_0(0) * occ_func_0_1(6) *
                occ_func_0_1(11))) +
          ((0.5773502692 * occ_func_0_1(0) * occ_func_0_1(4) * occ_func_0_0(6) +
            0.5773502692 * occ_func_0_1(0) * occ_func_0_0(4) * occ_func_0_1(6) +
            0.5773502692 * occ_func_0_0(0) * occ_func_0_1(4) *
                occ_func_0_1(6))) +
          ((0.5773502692 * occ_func_0_1(0) * occ_func_0_1(3) * occ_func_0_0(8) +
            0.5773502692 * occ_func_0_1(0) * occ_func_0_0(3) * occ_func_0_1(8) +
            0.5773502692 * occ_func_0_0(0) * occ_func_0_1(3) *
                occ_func_0_1(8))) +
          ((0.5773502692 * occ_func_0_1(0) * occ_func_0_1(10) *
                occ_func_0_0(12) +
            0.5773502692 * occ_func_0_1(0) * occ_func_0_0(10) *
                occ_func_0_1(12) +
            0.5773502692 * occ_func_0_0(0) * occ_func_0_1(10) *
                occ_func_0_1(12))) +
          ((0.5773502692 * occ_func_0_1(0) * occ_func_0_1(5) *
                occ_func_0_0(10) +
            0.5773502692 * occ_func_0_1(0) * occ_func_0_0(5) *
                occ_func_0_1(10) +
            0.5773502692 * occ_func_0_0(0) * occ_func_0_1(5) *
                occ_func_0_1(10))) +
          ((0.5773502692 * occ_func_0_1(0) * occ_func_0_1(4) * occ_func_0_0(5) +
            0.5773502692 * occ_func_0_1(0) * occ_func_0_0(4) * occ_func_0_1(5) +
            0.5773502692 * occ_func_0_0(0) * occ_func_0_1(4) *
                occ_func_0_1(5))) +
          ((0.5773502692 * occ_func_0_1(0) * occ_func_0_1(9) *
                occ_func_0_0(12) +
            0.5773502692 * occ_func_0_1(0) * occ_func_0_0(9) *
                occ_func_0_1(12) +
            0.5773502692 * occ_func_0_0(0) * occ_func_0_1(9) *
                occ_func_0_1(12))) +
          ((0.5773502692 * occ_func_0_1(0) * occ_func_0_1(2) * occ_func_0_0(7) +
            0.5773502692 * occ_func_0_1(0) * occ_func_0_0(2) * occ_func_0_1(7) +
            0.5773502692 * occ_func_0_0(0) * occ_func_0_1(2) *
                occ_func_0_1(7)))) /
         8.0;
}
double test_Clexulator::eval_bfunc_3_0_3() const {
  return ((occ_func_0_1(0) * occ_func_0_1(6) * occ_func_0_1(11)) +
          (occ_func_0_1(0) * occ_func_0_1(4) * occ_func_0_1(6)) +
          (occ_func_0_1(0) * occ_func_0_1(3) * occ_func_0_1(8)) +
          (occ_func_0_1(0) * occ_func_0_1(10) * occ_func_0_1(12)) +
          (occ_func_0_1(0) * occ_func_0_1(5) * occ_func_0_1(10)) +
          (occ_func_0_1(0) * occ_func_0_1(4) * occ_func_0_1(5)) +
          (occ_func_0_1(0) * occ_func_0_1(9) * occ_func_0_1(12)) +
          (occ_func_0_1(0) * occ_func_0_1(2) * occ_func_0_1(7))) /
         8.0;
}

double test_Clexulator::site_eval_at_0_bfunc_3_0_0() const {
  return ((occ_func_0_0(0) * occ_func_0_0(6) * occ_func_0_0(11)) +
          (occ_func_0_0(7) * occ_func_0_0(0) * occ_func_0_0(12)) +
          (occ_func_0_0(2) * occ_func_0_0(1) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(4) * occ_func_0_0(6)) +
          (occ_func_0_0(9) * occ_func_0_0(0) * occ_func_0_0(8)) +
          (occ_func_0_0(7) * occ_func_0_0(5) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(3) * occ_func_0_0(8)) +
          (occ_func_0_0(10) * occ_func_0_0(0) * occ_func_0_0(11)) +
          (occ_func_0_0(5) * occ_func_0_0(2) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(10) * occ_func_0_0(12)) +
          (occ_func_0_0(3) * occ_func_0_0(0) * occ_func_0_0(9)) +
          (occ_func_0_0(1) * occ_func_0_0(4) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(5) * occ_func_0_0(10)) +
          (occ_func_0_0(8) * occ_func_0_0(0) * occ_func_0_0(11)) +
          (occ_func_0_0(3) * occ_func_0_0(2) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(4) * occ_func_0_0(5)) +
          (occ_func_0_0(9) * occ_func_0_0(0) * occ_func_0_0(7)) +
          (occ_func_0_0(8) * occ_func_0_0(6) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(9) * occ_func_0_0(12)) +
          (occ_func_0_0(4) * occ_func_0_0(0) * occ_func_0_0(10)) +
          (occ_func_0_0(1) * occ_func_0_0(3) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(2) * occ_func_0_0(7)) +
          (occ_func_0_0(11) * occ_func_0_0(0) * occ_func_0_0(12)) +
          (occ_func_0_0(6) * occ_func_0_0(1) * occ_func_0_0(0))) /
         8.0;
}
double test_Clexulator::site_eval_at_0_bfunc_3_0_1() const {
  return (((0.5773502692 * occ_func_0_1(0) * occ_func_0_0(6) *
                occ_func_0_0(11) +
            0.5773502692 * occ_func_0_0(0) * occ_func_0_1(6) *
                occ_func_0_0(11) +
            0.5773502692 * occ_func_0_0(0) * occ_func_0_0(6) *
                occ_func_0_1(11))) +
          ((0.5773502692 * occ_func_0_1(7) * occ_func_0_0(0) *
                occ_func_0_0(12) +
            0.5773502692 * occ_func_0_0(7) * occ_func_0_1(0) *
                occ_func_0_0(12) +
            0.5773502692 * occ_func_0_0(7) * occ_func_0_0(0) *
                occ_func_0_1(12))) +
          ((0.5773502692 * occ_func_0_1(2) * occ_func_0_0(1) * occ_func_0_0(0) +
            0.5773502692 * occ_func_0_0(2) * occ_func_0_1(1) * occ_func_0_0(0) +
            0.5773502692 * occ_func_0_0(2) * occ_func_0_0(1) *
                occ_func_0_1(0))) +
          ((0.5773502692 * occ_func_0_1(0) * occ_func_0_0(4) * occ_func_0_0(6) +
            0.5773502692 * occ_func_0_0(0) * occ_func_0_1(4) * occ_func_0_0(6) +
            0.5773502692 * occ_func_0_0(0) * occ_func_0_0(4) *
                occ_func_0_1(6))) +
          ((0.5773502692 * occ_func_0_1(9) * occ_func_0_0(0) * occ_func_0_0(8) +
            0.5773502692 * occ_func_0_0(9) * occ_func_0_1(0) * occ_func_0_0(8) +
            0.5773502692 * occ_func_0_0(9) * occ_func_0_0(0) *
                occ_func_0_1(8))) +
          ((0.5773502692 * occ_func_0_1(7) * occ_func_0_0(5) * occ_func_0_0(0) +
            0.5773502692 * occ_func_0_0(7) * occ_func_0_1(5) * occ_func_0_0(0) +
            0.5773502692 * occ_func_0_0(7) * occ_func_0_0(5) *
                occ_func_0_1(0))) +
          ((0.5773502692 * occ_func_0_1(0) * occ_func_0_0(3) * occ_func_0_0(8) +
            0.5773502692 * occ_func_0_0(0) * occ_func_0_1(3) * occ_func_0_0(8) +
            0.5773502692 * occ_func_0_0(0) * occ_func_0_0(3) *
                occ_func_0_1(8))) +
          ((0.5773502692 * occ_func_0_1(10) * occ_func_0_0(0) *
                occ_func_0_0(11) +
            0.5773502692 * occ_func_0_0(10) * occ_func_0_1(0) *
                occ_func_0_0(11) +
            0.5773502692 * occ_func_0_0(10) * occ_func_0_0(0) *
                occ_func_0_1(11))) +
          ((0.5773502692 * occ_func_0_1(5) * occ_func_0_0(2) * occ_func_0_0(0) +
            0.5773502692 * occ_func_0_0(5) * occ_func_0_1(2) * occ_func_0_0(0) +
            0.5773502692 * occ_func_0_0(5) * occ_func_0_0(2) *
                occ_func_0_1(0))) +
          ((0.5773502692 * occ_func_0_1(0) * occ_func_0_0(10) *
                occ_func_0_0(12) +
            0.5773502692 * occ_func_0_0(0) * occ_func_0_1(10) *
                occ_func_0_0(12) +
            0.5773502692 * occ_func_0_0(0) * occ_func_0_0(10) *
                occ_func_0_1(12))) +
          ((0.5773502692 * occ_func_0_1(3) * occ_func_0_0(0) * occ_func_0_0(9) +
            0.5773502692 * occ_func_0_0(3) * occ_func_0_1(0) * occ_func_0_0(9) +
            0.5773502692 * occ_func_0_0(3) * occ_func_0_0(0) *
                occ_func_0_1(9))) +
          ((0.5773502692 * occ_func_0_1(1) * occ_func_0_0(4) * occ_func_0_0(0) +
            0.5773502692 * occ_func_0_0(1) * occ_func_0_1(4) * occ_func_0_0(0) +
            0.5773502692 * occ_func_0_0(1) * occ_func_0_0(4) *
                occ_func_0_1(0))) +
          ((0.5773502692 * occ_func_0_1(0) * occ_func_0_0(5) *
                occ_func_0_0(10) +
            0.5773502692 * occ_func_0_0(0) * occ_func_0_1(5) *
                occ_func_0_0(10) +
            0.5773502692 * occ_func_0_0(0) * occ_func_0_0(5) *
                occ_func_0_1(10))) +
          ((0.5773502692 * occ_func_0_1(8) * occ_func_0_0(0) *
                occ_func_0_0(11) +
            0.5773502692 * occ_func_0_0(8) * occ_func_0_1(0) *
                occ_func_0_0(11) +
            0.5773502692 * occ_func_0_0(8) * occ_func_0_0(0) *
                occ_func_0_1(11))) +
          ((0.5773502692 * occ_func_0_1(3) * occ_func_0_0(2) * occ_func_0_0(0) +
            0.5773502692 * occ_func_0_0(3) * occ_func_0_1(2) * occ_func_0_0(0) +
            0.5773502692 * occ_func_0_0(3) * occ_func_0_0(2) *
                occ_func_0_1(0))) +
          ((0.5773502692 * occ_func_0_1(0) * occ_func_0_0(4) * occ_func_0_0(5) +
            0.5773502692 * occ_func_0_0(0) * occ_func_0_1(4) * occ_func_0_0(5) +
            0.5773502692 * occ_func_0_0(0) * occ_func_0_0(4) *
                occ_func_0_1(5))) +
          ((0.5773502692 * occ_func_0_1(9) * occ_func_0_0(0) * occ_func_0_0(7) +
            0.5773502692 * occ_func_0_0(9) * occ_func_0_1(0) * occ_func_0_0(7) +
            0.5773502692 * occ_func_0_0(9) * occ_func_0_0(0) *
                occ_func_0_1(7))) +
          ((0.5773502692 * occ_func_0_1(8) * occ_func_0_0(6) * occ_func_0_0(0) +
            0.5773502692 * occ_func_0_0(8) * occ_func_0_1(6) * occ_func_0_0(0) +
            0.5773502692 * occ_func_0_0(8) * occ_func_0_0(6) *
                occ_func_0_1(0))) +
          ((0.5773502692 * occ_func_0_1(0) * occ_func_0_0(9) *
                occ_func_0_0(12) +
            0.5773502692 * occ_func_0_0(0) * occ_func_0_1(9) *
                occ_func_0_0(12) +
            0.5773502692 * occ_func_0_0(0) * occ_func_0_0(9) *
                occ_func_0_1(12))) +
          ((0.5773502692 * occ_func_0_1(4) * occ_func_0_0(0) *
                occ_func_0_0(10) +
            0.5773502692 * occ_func_0_0(4) * occ_func_0_1(0) *
                occ_func_0_0(10) +
            0.5773502692 * occ_func_0_0(4) * occ_func_0_0(0) *
                occ_func_0_1(10))) +
          ((0.5773502692 * occ_func_0_1(1) * occ_func_0_0(3) * occ_func_0_0(0) +
            0.5773502692 * occ_func_0_0(1) * occ_func_0_1(3) * occ_func_0_0(0) +
            0.5773502692 * occ_func_0_0(1) * occ_func_0_0(3) *
                occ_func_0_1(0))) +
          ((0.5773502692 * occ_func_0_1(0) * occ_func_0_0(2) * occ_func_0_0(7) +
            0.5773502692 * occ_func_0_0(0) * occ_func_0_1(2) * occ_func_0_0(7) +
            0.5773502692 * occ_func_0_0(0) * occ_func_0_0(2) *
                occ_func_0_1(7))) +
          ((0.5773502692 * occ_func_0_1(11) * occ_func_0_0(0) *
                occ_func_0_0(12) +
            0.5773502692 * occ_func_0_0(11) * occ_func_0_1(0) *
                occ_func_0_0(12) +
            0.5773502692 * occ_func_0_0(11) * occ_func_0_0(0) *
                occ_func_0_1(12))) +
          ((0.5773502692 * occ_func_0_1(6) * occ_func_0_0(1) * occ_func_0_0(0) +
            0.5773502692 * occ_func_0_0(6) * occ_func_0_1(1) * occ_func_0_0(0) +
            0.5773502692 * occ_func_0_0(6) * occ_func_0_0(1) *
                occ_func_0_1(0)))) /
         8.0;
}
double test_Clexulator::site_eval_at_0_bfunc_3_0_2() const {
  return (((0.5773502692 * occ_func_0_1(0) * occ_func_0_1(6) *
                occ_func_0_0(11) +
            0.5773502692 * occ_func_0_1(0) * occ_func_0_0(6) *
                occ_func_0_1(11) +
            0.5773502692 * occ_func_0_0(0) * occ_func_0_1(6) *
                occ_func_0_1(11))) +
          ((0.5773502692 * occ_func_0_1(7) * occ_func_0_1(0) *
                occ_func_0_0(12) +
            0.5773502692 * occ_func_0_1(7) * occ_func_0_0(0) *
                occ_func_0_1(12) +
            0.5773502692 * occ_func_0_0(7) * occ_func_0_1(0) *
                occ_func_0_1(12))) +
          ((0.5773502692 * occ_func_0_1(2) * occ_func_0_1(1) * occ_func_0_0(0) +
            0.5773502692 * occ_func_0_1(2) * occ_func_0_0(1) * occ_func_0_1(0) +
            0.5773502692 * occ_func_0_0(2) * occ_func_0_1(1) *
                occ_func_0_1(0))) +
          ((0.5773502692 * occ_func_0_1(0) * occ_func_0_1(4) * occ_func_0_0(6) +
            0.5773502692 * occ_func_0_1(0) * occ_func_0_0(4) * occ_func_0_1(6) +
            0.5773502692 * occ_func_0_0(0) * occ_func_0_1(4) *
                occ_func_0_1(6))) +
          ((0.5773502692 * occ_func_0_1(9) * occ_func_0_1(0) * occ_func_0_0(8) +
            0.5773502692 * occ_func_0_1(9) * occ_func_0_0(0) * occ_func_0_1(8) +
            0.5773502692 * occ_func_0_0(9) * occ_func_0_1(0) *
                occ_func_0_1(8))) +
          ((0.5773502692 * occ_func_0_1(7) * occ_func_0_1(5) * occ_func_0_0(0) +
            0.5773502692 * occ_func_0_1(7) * occ_func_0_0(5) * occ_func_0_1(0) +
            0.5773502692 * occ_func_0_0(7) * occ_func_0_1(5) *
                occ_func_0_1(0))) +
          ((0.5773502692 * occ_func_0_1(0) * occ_func_0_1(3) * occ_func_0_0(8) +
            0.5773502692 * occ_func_0_1(0) * occ_func_0_0(3) * occ_func_0_1(8) +
            0.5773502692 * occ_func_0_0(0) * occ_func_0_1(3) *
                occ_func_0_1(8))) +
          ((0.5773502692 * occ_func_0_1(10) * occ_func_0_1(0) *
                occ_func_0_0(11) +
            0.5773502692 * occ_func_0_1(10) * occ_func_0_0(0) *
                occ_func_0_1(11) +
            0.5773502692 * occ_func_0_0(10) * occ_func_0_1(0) *
                occ_func_0_1(11))) +
          ((0.5773502692 * occ_func_0_1(5) * occ_func_0_1(2) * occ_func_0_0(0) +
            0.5773502692 * occ_func_0_1(5) * occ_func_0_0(2) * occ_func_0_1(0) +
            0.5773502692 * occ_func_0_0(5) * occ_func_0_1(2) *
                occ_func_0_1(0))) +
          ((0.5773502692 * occ_func_0_1(0) * occ_func_0_1(10) *
                occ_func_0_0(12) +
            0.5773502692 * occ_func_0_1(0) * occ_func_0_0(10) *
                occ_func_0_1(12) +
            0.5773502692 * occ_func_0_0(0) * occ_func_0_1(10) *
                occ_func_0_1(12))) +
          ((0.5773502692 * occ_func_0_1(3) * occ_func_0_1(0) * occ_func_0_0(9) +
            0.5773502692 * occ_func_0_1(3) * occ_func_0_0(0) * occ_func_0_1(9) +
            0.5773502692 * occ_func_0_0(3) * occ_func_0_1(0) *
                occ_func_0_1(9))) +
          ((0.5773502692 * occ_func_0_1(1) * occ_func_0_1(4) * occ_func_0_0(0) +
            0.5773502692 * occ_func_0_1(1) * occ_func_0_0(4) * occ_func_0_1(0) +
            0.5773502692 * occ_func_0_0(1) * occ_func_0_1(4) *
                occ_func_0_1(0))) +
          ((0.5773502692 * occ_func_0_1(0) * occ_func_0_1(5) *
                occ_func_0_0(10) +
            0.5773502692 * occ_func_0_1(0) * occ_func_0_0(5) *
                occ_func_0_1(10) +
            0.5773502692 * occ_func_0_0(0) * occ_func_0_1(5) *
                occ_func_0_1(10))) +
          ((0.5773502692 * occ_func_0_1(8) * occ_func_0_1(0) *
                occ_func_0_0(11) +
            0.5773502692 * occ_func_0_1(8) * occ_func_0_0(0) *
                occ_func_0_1(11) +
            0.5773502692 * occ_func_0_0(8) * occ_func_0_1(0) *
                occ_func_0_1(11))) +
          ((0.5773502692 * occ_func_0_1(3) * occ_func_0_1(2) * occ_func_0_0(0) +
            0.5773502692 * occ_func_0_1(3) * occ_func_0_0(2) * occ_func_0_1(0) +
            0.5773502692 * occ_func_0_0(3) * occ_func_0_1(2) *
                occ_func_0_1(0))) +
          ((0.5773502692 * occ_func_0_1(0) * occ_func_0_1(4) * occ_func_0_0(5) +
            0.5773502692 * occ_func_0_1(0) * occ_func_0_0(4) * occ_func_0_1(5) +
            0.5773502692 * occ_func_0_0(0) * occ_func_0_1(4) *
                occ_func_0_1(5))) +
          ((0.5773502692 * occ_func_0_1(9) * occ_func_0_1(0) * occ_func_0_0(7) +
            0.5773502692 * occ_func_0_1(9) * occ_func_0_0(0) * occ_func_0_1(7) +
            0.5773502692 * occ_func_0_0(9) * occ_func_0_1(0) *
                occ_func_0_1(7))) +
          ((0.5773502692 * occ_func_0_1(8) * occ_func_0_1(6) * occ_func_0_0(0) +
            0.5773502692 * occ_func_0_1(8) * occ_func_0_0(6) * occ_func_0_1(0) +
            0.5773502692 * occ_func_0_0(8) * occ_func_0_1(6) *
                occ_func_0_1(0))) +
          ((0.5773502692 * occ_func_0_1(0) * occ_func_0_1(9) *
                occ_func_0_0(12) +
            0.5773502692 * occ_func_0_1(0) * occ_func_0_0(9) *
                occ_func_0_1(12) +
            0.5773502692 * occ_func_0_0(0) * occ_func_0_1(9) *
                occ_func_0_1(12))) +
          ((0.5773502692 * occ_func_0_1(4) * occ_func_0_1(0) *
                occ_func_0_0(10) +
            0.5773502692 * occ_func_0_1(4) * occ_func_0_0(0) *
                occ_func_0_1(10) +
            0.5773502692 * occ_func_0_0(4) * occ_func_0_1(0) *
                occ_func_0_1(10))) +
          ((0.5773502692 * occ_func_0_1(1) * occ_func_0_1(3) * occ_func_0_0(0) +
            0.5773502692 * occ_func_0_1(1) * occ_func_0_0(3) * occ_func_0_1(0) +
            0.5773502692 * occ_func_0_0(1) * occ_func_0_1(3) *
                occ_func_0_1(0))) +
          ((0.5773502692 * occ_func_0_1(0) * occ_func_0_1(2) * occ_func_0_0(7) +
            0.5773502692 * occ_func_0_1(0) * occ_func_0_0(2) * occ_func_0_1(7) +
            0.5773502692 * occ_func_0_0(0) * occ_func_0_1(2) *
                occ_func_0_1(7))) +
          ((0.5773502692 * occ_func_0_1(11) * occ_func_0_1(0) *
                occ_func_0_0(12) +
            0.5773502692 * occ_func_0_1(11) * occ_func_0_0(0) *
                occ_func_0_1(12) +
            0.5773502692 * occ_func_0_0(11) * occ_func_0_1(0) *
                occ_func_0_1(12))) +
          ((0.5773502692 * occ_func_0_1(6) * occ_func_0_1(1) * occ_func_0_0(0) +
            0.5773502692 * occ_func_0_1(6) * occ_func_0_0(1) * occ_func_0_1(0) +
            0.5773502692 * occ_func_0_0(6) * occ_func_0_1(1) *
                occ_func_0_1(0)))) /
         8.0;
}
double test_Clexulator::site_eval_at_0_bfunc_3_0_3() const {
  return ((occ_func_0_1(0) * occ_func_0_1(6) * occ_func_0_1(11)) +
          (occ_func_0_1(7) * occ_func_0_1(0) * occ_func_0_1(12)) +
          (occ_func_0_1(2) * occ_func_0_1(1) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(4) * occ_func_0_1(6)) +
          (occ_func_0_1(9) * occ_func_0_1(0) * occ_func_0_1(8)) +
          (occ_func_0_1(7) * occ_func_0_1(5) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(3) * occ_func_0_1(8)) +
          (occ_func_0_1(10) * occ_func_0_1(0) * occ_func_0_1(11)) +
          (occ_func_0_1(5) * occ_func_0_1(2) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(10) * occ_func_0_1(12)) +
          (occ_func_0_1(3) * occ_func_0_1(0) * occ_func_0_1(9)) +
          (occ_func_0_1(1) * occ_func_0_1(4) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(5) * occ_func_0_1(10)) +
          (occ_func_0_1(8) * occ_func_0_1(0) * occ_func_0_1(11)) +
          (occ_func_0_1(3) * occ_func_0_1(2) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(4) * occ_func_0_1(5)) +
          (occ_func_0_1(9) * occ_func_0_1(0) * occ_func_0_1(7)) +
          (occ_func_0_1(8) * occ_func_0_1(6) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(9) * occ_func_0_1(12)) +
          (occ_func_0_1(4) * occ_func_0_1(0) * occ_func_0_1(10)) +
          (occ_func_0_1(1) * occ_func_0_1(3) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(2) * occ_func_0_1(7)) +
          (occ_func_0_1(11) * occ_func_0_1(0) * occ_func_0_1(12)) +
          (occ_func_0_1(6) * occ_func_0_1(1) * occ_func_0_1(0))) /
         8.0;
}

double test_Clexulator::delta_site_eval_at_0_bfunc_3_0_0(int occ_i,
                                                         int occ_f) const {
  return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) *
         ((occ_func_0_0(6) * occ_func_0_0(11)) +
          (occ_func_0_0(7) * occ_func_0_0(12)) +
          (occ_func_0_0(2) * occ_func_0_0(1)) +
          (occ_func_0_0(4) * occ_func_0_0(6)) +
          (occ_func_0_0(9) * occ_func_0_0(8)) +
          (occ_func_0_0(7) * occ_func_0_0(5)) +
          (occ_func_0_0(3) * occ_func_0_0(8)) +
          (occ_func_0_0(10) * occ_func_0_0(11)) +
          (occ_func_0_0(5) * occ_func_0_0(2)) +
          (occ_func_0_0(10) * occ_func_0_0(12)) +
          (occ_func_0_0(3) * occ_func_0_0(9)) +
          (occ_func_0_0(1) * occ_func_0_0(4)) +
          (occ_func_0_0(5) * occ_func_0_0(10)) +
          (occ_func_0_0(8) * occ_func_0_0(11)) +
          (occ_func_0_0(3) * occ_func_0_0(2)) +
          (occ_func_0_0(4) * occ_func_0_0(5)) +
          (occ_func_0_0(9) * occ_func_0_0(7)) +
          (occ_func_0_0(8) * occ_func_0_0(6)) +
          (occ_func_0_0(9) * occ_func_0_0(12)) +
          (occ_func_0_0(4) * occ_func_0_0(10)) +
          (occ_func_0_0(1) * occ_func_0_0(3)) +
          (occ_func_0_0(2) * occ_func_0_0(7)) +
          (occ_func_0_0(11) * occ_func_0_0(12)) +
          (occ_func_0_0(6) * occ_func_0_0(1))) /
         8.0;
}
double test_Clexulator::delta_site_eval_at_0_bfunc_3_0_1(int occ_i,
                                                         int occ_f) const {
  return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) *
             (((0.5773502692 * occ_func_0_1(6) * occ_func_0_0(11) +
                0.5773502692 * occ_func_0_0(6) * occ_func_0_1(11))) +
              ((0.5773502692 * occ_func_0_1(7) * occ_func_0_0(12) +
                0.5773502692 * occ_func_0_0(7) * occ_func_0_1(12))) +
              ((0.5773502692 * occ_func_0_1(2) * occ_func_0_0(1) +
                0.5773502692 * occ_func_0_0(2) * occ_func_0_1(1))) +
              ((0.5773502692 * occ_func_0_1(4) * occ_func_0_0(6) +
                0.5773502692 * occ_func_0_0(4) * occ_func_0_1(6))) +
              ((0.5773502692 * occ_func_0_1(9) * occ_func_0_0(8) +
                0.5773502692 * occ_func_0_0(9) * occ_func_0_1(8))) +
              ((0.5773502692 * occ_func_0_1(7) * occ_func_0_0(5) +
                0.5773502692 * occ_func_0_0(7) * occ_func_0_1(5))) +
              ((0.5773502692 * occ_func_0_1(3) * occ_func_0_0(8) +
                0.5773502692 * occ_func_0_0(3) * occ_func_0_1(8))) +
              ((0.5773502692 * occ_func_0_1(10) * occ_func_0_0(11) +
                0.5773502692 * occ_func_0_0(10) * occ_func_0_1(11))) +
              ((0.5773502692 * occ_func_0_1(5) * occ_func_0_0(2) +
                0.5773502692 * occ_func_0_0(5) * occ_func_0_1(2))) +
              ((0.5773502692 * occ_func_0_1(10) * occ_func_0_0(12) +
                0.5773502692 * occ_func_0_0(10) * occ_func_0_1(12))) +
              ((0.5773502692 * occ_func_0_1(3) * occ_func_0_0(9) +
                0.5773502692 * occ_func_0_0(3) * occ_func_0_1(9))) +
              ((0.5773502692 * occ_func_0_1(1) * occ_func_0_0(4) +
                0.5773502692 * occ_func_0_0(1) * occ_func_0_1(4))) +
              ((0.5773502692 * occ_func_0_1(5) * occ_func_0_0(10) +
                0.5773502692 * occ_func_0_0(5) * occ_func_0_1(10))) +
              ((0.5773502692 * occ_func_0_1(8) * occ_func_0_0(11) +
                0.5773502692 * occ_func_0_0(8) * occ_func_0_1(11))) +
              ((0.5773502692 * occ_func_0_1(3) * occ_func_0_0(2) +
                0.5773502692 * occ_func_0_0(3) * occ_func_0_1(2))) +
              ((0.5773502692 * occ_func_0_1(4) * occ_func_0_0(5) +
                0.5773502692 * occ_func_0_0(4) * occ_func_0_1(5))) +
              ((0.5773502692 * occ_func_0_1(9) * occ_func_0_0(7) +
                0.5773502692 * occ_func_0_0(9) * occ_func_0_1(7))) +
              ((0.5773502692 * occ_func_0_1(8) * occ_func_0_0(6) +
                0.5773502692 * occ_func_0_0(8) * occ_func_0_1(6))) +
              ((0.5773502692 * occ_func_0_1(9) * occ_func_0_0(12) +
                0.5773502692 * occ_func_0_0(9) * occ_func_0_1(12))) +
              ((0.5773502692 * occ_func_0_1(4) * occ_func_0_0(10) +
                0.5773502692 * occ_func_0_0(4) * occ_func_0_1(10))) +
              ((0.5773502692 * occ_func_0_1(1) * occ_func_0_0(3) +
                0.5773502692 * occ_func_0_0(1) * occ_func_0_1(3))) +
              ((0.5773502692 * occ_func_0_1(2) * occ_func_0_0(7) +
                0.5773502692 * occ_func_0_0(2) * occ_func_0_1(7))) +
              ((0.5773502692 * occ_func_0_1(11) * occ_func_0_0(12) +
                0.5773502692 * occ_func_0_0(11) * occ_func_0_1(12))) +
              ((0.5773502692 * occ_func_0_1(6) * occ_func_0_0(1) +
                0.5773502692 * occ_func_0_0(6) * occ_func_0_1(1)))) /
             8.0 +
         (m_occ_func_0_1[occ_f] - m_occ_func_0_1[occ_i]) *
             ((0.5773502692 * occ_func_0_0(6) * occ_func_0_0(11)) +
              (0.5773502692 * occ_func_0_0(7) * occ_func_0_0(12)) +
              (0.5773502692 * occ_func_0_0(2) * occ_func_0_0(1)) +
              (0.5773502692 * occ_func_0_0(4) * occ_func_0_0(6)) +
              (0.5773502692 * occ_func_0_0(9) * occ_func_0_0(8)) +
              (0.5773502692 * occ_func_0_0(7) * occ_func_0_0(5)) +
              (0.5773502692 * occ_func_0_0(3) * occ_func_0_0(8)) +
              (0.5773502692 * occ_func_0_0(10) * occ_func_0_0(11)) +
              (0.5773502692 * occ_func_0_0(5) * occ_func_0_0(2)) +
              (0.5773502692 * occ_func_0_0(10) * occ_func_0_0(12)) +
              (0.5773502692 * occ_func_0_0(3) * occ_func_0_0(9)) +
              (0.5773502692 * occ_func_0_0(1) * occ_func_0_0(4)) +
              (0.5773502692 * occ_func_0_0(5) * occ_func_0_0(10)) +
              (0.5773502692 * occ_func_0_0(8) * occ_func_0_0(11)) +
              (0.5773502692 * occ_func_0_0(3) * occ_func_0_0(2)) +
              (0.5773502692 * occ_func_0_0(4) * occ_func_0_0(5)) +
              (0.5773502692 * occ_func_0_0(9) * occ_func_0_0(7)) +
              (0.5773502692 * occ_func_0_0(8) * occ_func_0_0(6)) +
              (0.5773502692 * occ_func_0_0(9) * occ_func_0_0(12)) +
              (0.5773502692 * occ_func_0_0(4) * occ_func_0_0(10)) +
              (0.5773502692 * occ_func_0_0(1) * occ_func_0_0(3)) +
              (0.5773502692 * occ_func_0_0(2) * occ_func_0_0(7)) +
              (0.5773502692 * occ_func_0_0(11) * occ_func_0_0(12)) +
              (0.5773502692 * occ_func_0_0(6) * occ_func_0_0(1))) /
             8.0;
}
double test_Clexulator::delta_site_eval_at_0_bfunc_3_0_2(int occ_i,
                                                         int occ_f) const {
  return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) *
             ((0.5773502692 * occ_func_0_1(6) * occ_func_0_1(11)) +
              (0.5773502692 * occ_func_0_1(7) * occ_func_0_1(12)) +
              (0.5773502692 * occ_func_0_1(2) * occ_func_0_1(1)) +
              (0.5773502692 * occ_func_0_1(4) * occ_func_0_1(6)) +
              (0.5773502692 * occ_func_0_1(9) * occ_func_0_1(8)) +
              (0.5773502692 * occ_func_0_1(7) * occ_func_0_1(5)) +
              (0.5773502692 * occ_func_0_1(3) * occ_func_0_1(8)) +
              (0.5773502692 * occ_func_0_1(10) * occ_func_0_1(11)) +
              (0.5773502692 * occ_func_0_1(5) * occ_func_0_1(2)) +
              (0.5773502692 * occ_func_0_1(10) * occ_func_0_1(12)) +
              (0.5773502692 * occ_func_0_1(3) * occ_func_0_1(9)) +
              (0.5773502692 * occ_func_0_1(1) * occ_func_0_1(4)) +
              (0.5773502692 * occ_func_0_1(5) * occ_func_0_1(10)) +
              (0.5773502692 * occ_func_0_1(8) * occ_func_0_1(11)) +
              (0.5773502692 * occ_func_0_1(3) * occ_func_0_1(2)) +
              (0.5773502692 * occ_func_0_1(4) * occ_func_0_1(5)) +
              (0.5773502692 * occ_func_0_1(9) * occ_func_0_1(7)) +
              (0.5773502692 * occ_func_0_1(8) * occ_func_0_1(6)) +
              (0.5773502692 * occ_func_0_1(9) * occ_func_0_1(12)) +
              (0.5773502692 * occ_func_0_1(4) * occ_func_0_1(10)) +
              (0.5773502692 * occ_func_0_1(1) * occ_func_0_1(3)) +
              (0.5773502692 * occ_func_0_1(2) * occ_func_0_1(7)) +
              (0.5773502692 * occ_func_0_1(11) * occ_func_0_1(12)) +
              (0.5773502692 * occ_func_0_1(6) * occ_func_0_1(1))) /
             8.0 +
         (m_occ_func_0_1[occ_f] - m_occ_func_0_1[occ_i]) *
             (((0.5773502692 * occ_func_0_1(6) * occ_func_0_0(11) +
                0.5773502692 * occ_func_0_0(6) * occ_func_0_1(11))) +
              ((0.5773502692 * occ_func_0_1(7) * occ_func_0_0(12) +
                0.5773502692 * occ_func_0_0(7) * occ_func_0_1(12))) +
              ((0.5773502692 * occ_func_0_1(2) * occ_func_0_0(1) +
                0.5773502692 * occ_func_0_0(2) * occ_func_0_1(1))) +
              ((0.5773502692 * occ_func_0_1(4) * occ_func_0_0(6) +
                0.5773502692 * occ_func_0_0(4) * occ_func_0_1(6))) +
              ((0.5773502692 * occ_func_0_1(9) * occ_func_0_0(8) +
                0.5773502692 * occ_func_0_0(9) * occ_func_0_1(8))) +
              ((0.5773502692 * occ_func_0_1(7) * occ_func_0_0(5) +
                0.5773502692 * occ_func_0_0(7) * occ_func_0_1(5))) +
              ((0.5773502692 * occ_func_0_1(3) * occ_func_0_0(8) +
                0.5773502692 * occ_func_0_0(3) * occ_func_0_1(8))) +
              ((0.5773502692 * occ_func_0_1(10) * occ_func_0_0(11) +
                0.5773502692 * occ_func_0_0(10) * occ_func_0_1(11))) +
              ((0.5773502692 * occ_func_0_1(5) * occ_func_0_0(2) +
                0.5773502692 * occ_func_0_0(5) * occ_func_0_1(2))) +
              ((0.5773502692 * occ_func_0_1(10) * occ_func_0_0(12) +
                0.5773502692 * occ_func_0_0(10) * occ_func_0_1(12))) +
              ((0.5773502692 * occ_func_0_1(3) * occ_func_0_0(9) +
                0.5773502692 * occ_func_0_0(3) * occ_func_0_1(9))) +
              ((0.5773502692 * occ_func_0_1(1) * occ_func_0_0(4) +
                0.5773502692 * occ_func_0_0(1) * occ_func_0_1(4))) +
              ((0.5773502692 * occ_func_0_1(5) * occ_func_0_0(10) +
                0.5773502692 * occ_func_0_0(5) * occ_func_0_1(10))) +
              ((0.5773502692 * occ_func_0_1(8) * occ_func_0_0(11) +
                0.5773502692 * occ_func_0_0(8) * occ_func_0_1(11))) +
              ((0.5773502692 * occ_func_0_1(3) * occ_func_0_0(2) +
                0.5773502692 * occ_func_0_0(3) * occ_func_0_1(2))) +
              ((0.5773502692 * occ_func_0_1(4) * occ_func_0_0(5) +
                0.5773502692 * occ_func_0_0(4) * occ_func_0_1(5))) +
              ((0.5773502692 * occ_func_0_1(9) * occ_func_0_0(7) +
                0.5773502692 * occ_func_0_0(9) * occ_func_0_1(7))) +
              ((0.5773502692 * occ_func_0_1(8) * occ_func_0_0(6) +
                0.5773502692 * occ_func_0_0(8) * occ_func_0_1(6))) +
              ((0.5773502692 * occ_func_0_1(9) * occ_func_0_0(12) +
                0.5773502692 * occ_func_0_0(9) * occ_func_0_1(12))) +
              ((0.5773502692 * occ_func_0_1(4) * occ_func_0_0(10) +
                0.5773502692 * occ_func_0_0(4) * occ_func_0_1(10))) +
              ((0.5773502692 * occ_func_0_1(1) * occ_func_0_0(3) +
                0.5773502692 * occ_func_0_0(1) * occ_func_0_1(3))) +
              ((0.5773502692 * occ_func_0_1(2) * occ_func_0_0(7) +
                0.5773502692 * occ_func_0_0(2) * occ_func_0_1(7))) +
              ((0.5773502692 * occ_func_0_1(11) * occ_func_0_0(12) +
                0.5773502692 * occ_func_0_0(11) * occ_func_0_1(12))) +
              ((0.5773502692 * occ_func_0_1(6) * occ_func_0_0(1) +
                0.5773502692 * occ_func_0_0(6) * occ_func_0_1(1)))) /
             8.0;
}
double test_Clexulator::delta_site_eval_at_0_bfunc_3_0_3(int occ_i,
                                                         int occ_f) const {
  return (m_occ_func_0_1[occ_f] - m_occ_func_0_1[occ_i]) *
         ((occ_func_0_1(6) * occ_func_0_1(11)) +
          (occ_func_0_1(7) * occ_func_0_1(12)) +
          (occ_func_0_1(2) * occ_func_0_1(1)) +
          (occ_func_0_1(4) * occ_func_0_1(6)) +
          (occ_func_0_1(9) * occ_func_0_1(8)) +
          (occ_func_0_1(7) * occ_func_0_1(5)) +
          (occ_func_0_1(3) * occ_func_0_1(8)) +
          (occ_func_0_1(10) * occ_func_0_1(11)) +
          (occ_func_0_1(5) * occ_func_0_1(2)) +
          (occ_func_0_1(10) * occ_func_0_1(12)) +
          (occ_func_0_1(3) * occ_func_0_1(9)) +
          (occ_func_0_1(1) * occ_func_0_1(4)) +
          (occ_func_0_1(5) * occ_func_0_1(10)) +
          (occ_func_0_1(8) * occ_func_0_1(11)) +
          (occ_func_0_1(3) * occ_func_0_1(2)) +
          (occ_func_0_1(4) * occ_func_0_1(5)) +
          (occ_func_0_1(9) * occ_func_0_1(7)) +
          (occ_func_0_1(8) * occ_func_0_1(6)) +
          (occ_func_0_1(9) * occ_func_0_1(12)) +
          (occ_func_0_1(4) * occ_func_0_1(10)) +
          (occ_func_0_1(1) * occ_func_0_1(3)) +
          (occ_func_0_1(2) * occ_func_0_1(7)) +
          (occ_func_0_1(11) * occ_func_0_1(12)) +
          (occ_func_0_1(6) * occ_func_0_1(1))) /
         8.0;
}

/**** Basis functions for orbit 3, 1****
#Points: 3
MaxLength: 4.8989795  MinLength: 2.8284271
 0.0000000   0.0000000   0.0000000 A B C
 0.0000000  -1.0000000   1.0000000 A B C
 1.0000000   0.0000000   1.0000000 A B C
****/
double test_Clexulator::eval_bfunc_3_1_0() const {
  return ((occ_func_0_0(0) * occ_func_0_0(5) * occ_func_0_0(37)) +
          (occ_func_0_0(0) * occ_func_0_0(12) * occ_func_0_0(33)) +
          (occ_func_0_0(0) * occ_func_0_0(10) * occ_func_0_0(30)) +
          (occ_func_0_0(0) * occ_func_0_0(4) * occ_func_0_0(40)) +
          (occ_func_0_0(0) * occ_func_0_0(2) * occ_func_0_0(27)) +
          (occ_func_0_0(0) * occ_func_0_0(7) * occ_func_0_0(41)) +
          (occ_func_0_0(0) * occ_func_0_0(1) * occ_func_0_0(23)) +
          (occ_func_0_0(0) * occ_func_0_0(9) * occ_func_0_0(38)) +
          (occ_func_0_0(0) * occ_func_0_0(10) * occ_func_0_0(36)) +
          (occ_func_0_0(0) * occ_func_0_0(3) * occ_func_0_0(32)) +
          (occ_func_0_0(0) * occ_func_0_0(6) * occ_func_0_0(34)) +
          (occ_func_0_0(0) * occ_func_0_0(11) * occ_func_0_0(39)) +
          (occ_func_0_0(0) * occ_func_0_0(7) * occ_func_0_0(21)) +
          (occ_func_0_0(0) * occ_func_0_0(4) * occ_func_0_0(19)) +
          (occ_func_0_0(0) * occ_func_0_0(2) * occ_func_0_0(28)) +
          (occ_func_0_0(0) * occ_func_0_0(8) * occ_func_0_0(42)) +
          (occ_func_0_0(0) * occ_func_0_0(5) * occ_func_0_0(22)) +
          (occ_func_0_0(0) * occ_func_0_0(6) * occ_func_0_0(26)) +
          (occ_func_0_0(0) * occ_func_0_0(8) * occ_func_0_0(20)) +
          (occ_func_0_0(0) * occ_func_0_0(1) * occ_func_0_0(31)) +
          (occ_func_0_0(0) * occ_func_0_0(9) * occ_func_0_0(25)) +
          (occ_func_0_0(0) * occ_func_0_0(11) * occ_func_0_0(29)) +
          (occ_func_0_0(0) * occ_func_0_0(12) * occ_func_0_0(35)) +
          (occ_func_0_0(0) * occ_func_0_0(3) * occ_func_0_0(24))) /
         24.0;
}
double test_Clexulator::eval_bfunc_3_1_1() const {
  return (((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(5) *
                occ_func_0_0(37) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(5) *
                occ_func_0_0(37))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(12) *
                occ_func_0_0(33) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(12) *
                occ_func_0_0(33))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(10) *
                occ_func_0_0(30) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(10) *
                occ_func_0_0(30))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(4) *
                occ_func_0_0(40) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(4) *
                occ_func_0_0(40))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(2) *
                occ_func_0_0(27) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(2) *
                occ_func_0_0(27))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(7) *
                occ_func_0_0(41) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(7) *
                occ_func_0_0(41))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(1) *
                occ_func_0_0(23) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(1) *
                occ_func_0_0(23))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(9) *
                occ_func_0_0(38) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(9) *
                occ_func_0_0(38))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(10) *
                occ_func_0_0(36) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(10) *
                occ_func_0_0(36))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(3) *
                occ_func_0_0(32) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(3) *
                occ_func_0_0(32))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(6) *
                occ_func_0_0(34) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(6) *
                occ_func_0_0(34))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(11) *
                occ_func_0_0(39) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(11) *
                occ_func_0_0(39))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(7) *
                occ_func_0_0(21) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(7) *
                occ_func_0_0(21))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(4) *
                occ_func_0_0(19) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(4) *
                occ_func_0_0(19))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(2) *
                occ_func_0_0(28) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(2) *
                occ_func_0_0(28))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(8) *
                occ_func_0_0(42) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(8) *
                occ_func_0_0(42))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(5) *
                occ_func_0_0(22) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(5) *
                occ_func_0_0(22))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(6) *
                occ_func_0_0(26) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(6) *
                occ_func_0_0(26))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(8) *
                occ_func_0_0(20) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(8) *
                occ_func_0_0(20))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(1) *
                occ_func_0_0(31) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(1) *
                occ_func_0_0(31))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(9) *
                occ_func_0_0(25) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(9) *
                occ_func_0_0(25))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(11) *
                occ_func_0_0(29) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(11) *
                occ_func_0_0(29))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(12) *
                occ_func_0_0(35) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(12) *
                occ_func_0_0(35))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(3) *
                occ_func_0_0(24) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(3) *
                occ_func_0_0(24)))) /
         24.0;
}
double test_Clexulator::eval_bfunc_3_1_2() const {
  return ((occ_func_0_1(0) * occ_func_0_1(5) * occ_func_0_0(37)) +
          (occ_func_0_1(0) * occ_func_0_1(12) * occ_func_0_0(33)) +
          (occ_func_0_1(0) * occ_func_0_1(10) * occ_func_0_0(30)) +
          (occ_func_0_1(0) * occ_func_0_1(4) * occ_func_0_0(40)) +
          (occ_func_0_1(0) * occ_func_0_1(2) * occ_func_0_0(27)) +
          (occ_func_0_1(0) * occ_func_0_1(7) * occ_func_0_0(41)) +
          (occ_func_0_1(0) * occ_func_0_1(1) * occ_func_0_0(23)) +
          (occ_func_0_1(0) * occ_func_0_1(9) * occ_func_0_0(38)) +
          (occ_func_0_1(0) * occ_func_0_1(10) * occ_func_0_0(36)) +
          (occ_func_0_1(0) * occ_func_0_1(3) * occ_func_0_0(32)) +
          (occ_func_0_1(0) * occ_func_0_1(6) * occ_func_0_0(34)) +
          (occ_func_0_1(0) * occ_func_0_1(11) * occ_func_0_0(39)) +
          (occ_func_0_1(0) * occ_func_0_1(7) * occ_func_0_0(21)) +
          (occ_func_0_1(0) * occ_func_0_1(4) * occ_func_0_0(19)) +
          (occ_func_0_1(0) * occ_func_0_1(2) * occ_func_0_0(28)) +
          (occ_func_0_1(0) * occ_func_0_1(8) * occ_func_0_0(42)) +
          (occ_func_0_1(0) * occ_func_0_1(5) * occ_func_0_0(22)) +
          (occ_func_0_1(0) * occ_func_0_1(6) * occ_func_0_0(26)) +
          (occ_func_0_1(0) * occ_func_0_1(8) * occ_func_0_0(20)) +
          (occ_func_0_1(0) * occ_func_0_1(1) * occ_func_0_0(31)) +
          (occ_func_0_1(0) * occ_func_0_1(9) * occ_func_0_0(25)) +
          (occ_func_0_1(0) * occ_func_0_1(11) * occ_func_0_0(29)) +
          (occ_func_0_1(0) * occ_func_0_1(12) * occ_func_0_0(35)) +
          (occ_func_0_1(0) * occ_func_0_1(3) * occ_func_0_0(24))) /
         24.0;
}
double test_Clexulator::eval_bfunc_3_1_3() const {
  return ((occ_func_0_0(0) * occ_func_0_0(5) * occ_func_0_1(37)) +
          (occ_func_0_0(0) * occ_func_0_0(12) * occ_func_0_1(33)) +
          (occ_func_0_0(0) * occ_func_0_0(10) * occ_func_0_1(30)) +
          (occ_func_0_0(0) * occ_func_0_0(4) * occ_func_0_1(40)) +
          (occ_func_0_0(0) * occ_func_0_0(2) * occ_func_0_1(27)) +
          (occ_func_0_0(0) * occ_func_0_0(7) * occ_func_0_1(41)) +
          (occ_func_0_0(0) * occ_func_0_0(1) * occ_func_0_1(23)) +
          (occ_func_0_0(0) * occ_func_0_0(9) * occ_func_0_1(38)) +
          (occ_func_0_0(0) * occ_func_0_0(10) * occ_func_0_1(36)) +
          (occ_func_0_0(0) * occ_func_0_0(3) * occ_func_0_1(32)) +
          (occ_func_0_0(0) * occ_func_0_0(6) * occ_func_0_1(34)) +
          (occ_func_0_0(0) * occ_func_0_0(11) * occ_func_0_1(39)) +
          (occ_func_0_0(0) * occ_func_0_0(7) * occ_func_0_1(21)) +
          (occ_func_0_0(0) * occ_func_0_0(4) * occ_func_0_1(19)) +
          (occ_func_0_0(0) * occ_func_0_0(2) * occ_func_0_1(28)) +
          (occ_func_0_0(0) * occ_func_0_0(8) * occ_func_0_1(42)) +
          (occ_func_0_0(0) * occ_func_0_0(5) * occ_func_0_1(22)) +
          (occ_func_0_0(0) * occ_func_0_0(6) * occ_func_0_1(26)) +
          (occ_func_0_0(0) * occ_func_0_0(8) * occ_func_0_1(20)) +
          (occ_func_0_0(0) * occ_func_0_0(1) * occ_func_0_1(31)) +
          (occ_func_0_0(0) * occ_func_0_0(9) * occ_func_0_1(25)) +
          (occ_func_0_0(0) * occ_func_0_0(11) * occ_func_0_1(29)) +
          (occ_func_0_0(0) * occ_func_0_0(12) * occ_func_0_1(35)) +
          (occ_func_0_0(0) * occ_func_0_0(3) * occ_func_0_1(24))) /
         24.0;
}
double test_Clexulator::eval_bfunc_3_1_4() const {
  return (((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(5) *
                occ_func_0_1(37) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(5) *
                occ_func_0_1(37))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(12) *
                occ_func_0_1(33) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(12) *
                occ_func_0_1(33))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(10) *
                occ_func_0_1(30) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(10) *
                occ_func_0_1(30))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(4) *
                occ_func_0_1(40) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(4) *
                occ_func_0_1(40))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(2) *
                occ_func_0_1(27) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(2) *
                occ_func_0_1(27))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(7) *
                occ_func_0_1(41) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(7) *
                occ_func_0_1(41))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(1) *
                occ_func_0_1(23) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(1) *
                occ_func_0_1(23))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(9) *
                occ_func_0_1(38) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(9) *
                occ_func_0_1(38))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(10) *
                occ_func_0_1(36) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(10) *
                occ_func_0_1(36))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(3) *
                occ_func_0_1(32) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(3) *
                occ_func_0_1(32))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(6) *
                occ_func_0_1(34) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(6) *
                occ_func_0_1(34))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(11) *
                occ_func_0_1(39) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(11) *
                occ_func_0_1(39))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(7) *
                occ_func_0_1(21) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(7) *
                occ_func_0_1(21))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(4) *
                occ_func_0_1(19) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(4) *
                occ_func_0_1(19))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(2) *
                occ_func_0_1(28) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(2) *
                occ_func_0_1(28))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(8) *
                occ_func_0_1(42) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(8) *
                occ_func_0_1(42))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(5) *
                occ_func_0_1(22) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(5) *
                occ_func_0_1(22))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(6) *
                occ_func_0_1(26) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(6) *
                occ_func_0_1(26))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(8) *
                occ_func_0_1(20) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(8) *
                occ_func_0_1(20))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(1) *
                occ_func_0_1(31) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(1) *
                occ_func_0_1(31))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(9) *
                occ_func_0_1(25) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(9) *
                occ_func_0_1(25))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(11) *
                occ_func_0_1(29) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(11) *
                occ_func_0_1(29))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(12) *
                occ_func_0_1(35) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(12) *
                occ_func_0_1(35))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(3) *
                occ_func_0_1(24) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(3) *
                occ_func_0_1(24)))) /
         24.0;
}
double test_Clexulator::eval_bfunc_3_1_5() const {
  return ((occ_func_0_1(0) * occ_func_0_1(5) * occ_func_0_1(37)) +
          (occ_func_0_1(0) * occ_func_0_1(12) * occ_func_0_1(33)) +
          (occ_func_0_1(0) * occ_func_0_1(10) * occ_func_0_1(30)) +
          (occ_func_0_1(0) * occ_func_0_1(4) * occ_func_0_1(40)) +
          (occ_func_0_1(0) * occ_func_0_1(2) * occ_func_0_1(27)) +
          (occ_func_0_1(0) * occ_func_0_1(7) * occ_func_0_1(41)) +
          (occ_func_0_1(0) * occ_func_0_1(1) * occ_func_0_1(23)) +
          (occ_func_0_1(0) * occ_func_0_1(9) * occ_func_0_1(38)) +
          (occ_func_0_1(0) * occ_func_0_1(10) * occ_func_0_1(36)) +
          (occ_func_0_1(0) * occ_func_0_1(3) * occ_func_0_1(32)) +
          (occ_func_0_1(0) * occ_func_0_1(6) * occ_func_0_1(34)) +
          (occ_func_0_1(0) * occ_func_0_1(11) * occ_func_0_1(39)) +
          (occ_func_0_1(0) * occ_func_0_1(7) * occ_func_0_1(21)) +
          (occ_func_0_1(0) * occ_func_0_1(4) * occ_func_0_1(19)) +
          (occ_func_0_1(0) * occ_func_0_1(2) * occ_func_0_1(28)) +
          (occ_func_0_1(0) * occ_func_0_1(8) * occ_func_0_1(42)) +
          (occ_func_0_1(0) * occ_func_0_1(5) * occ_func_0_1(22)) +
          (occ_func_0_1(0) * occ_func_0_1(6) * occ_func_0_1(26)) +
          (occ_func_0_1(0) * occ_func_0_1(8) * occ_func_0_1(20)) +
          (occ_func_0_1(0) * occ_func_0_1(1) * occ_func_0_1(31)) +
          (occ_func_0_1(0) * occ_func_0_1(9) * occ_func_0_1(25)) +
          (occ_func_0_1(0) * occ_func_0_1(11) * occ_func_0_1(29)) +
          (occ_func_0_1(0) * occ_func_0_1(12) * occ_func_0_1(35)) +
          (occ_func_0_1(0) * occ_func_0_1(3) * occ_func_0_1(24))) /
         24.0;
}

double test_Clexulator::site_eval_at_0_bfunc_3_1_0() const {
  return ((occ_func_0_0(0) * occ_func_0_0(5) * occ_func_0_0(37)) +
          (occ_func_0_0(8) * occ_func_0_0(0) * occ_func_0_0(39)) +
          (occ_func_0_0(24) * occ_func_0_0(22) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(12) * occ_func_0_0(33)) +
          (occ_func_0_0(1) * occ_func_0_0(0) * occ_func_0_0(26)) +
          (occ_func_0_0(28) * occ_func_0_0(35) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(10) * occ_func_0_0(30)) +
          (occ_func_0_0(3) * occ_func_0_0(0) * occ_func_0_0(25)) +
          (occ_func_0_0(31) * occ_func_0_0(36) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(4) * occ_func_0_0(40)) +
          (occ_func_0_0(9) * occ_func_0_0(0) * occ_func_0_0(42)) +
          (occ_func_0_0(21) * occ_func_0_0(19) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(2) * occ_func_0_0(27)) +
          (occ_func_0_0(11) * occ_func_0_0(0) * occ_func_0_0(33)) +
          (occ_func_0_0(34) * occ_func_0_0(28) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(7) * occ_func_0_0(41)) +
          (occ_func_0_0(6) * occ_func_0_0(0) * occ_func_0_0(40)) +
          (occ_func_0_0(20) * occ_func_0_0(21) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(1) * occ_func_0_0(23)) +
          (occ_func_0_0(12) * occ_func_0_0(0) * occ_func_0_0(30)) +
          (occ_func_0_0(38) * occ_func_0_0(31) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(9) * occ_func_0_0(38)) +
          (occ_func_0_0(4) * occ_func_0_0(0) * occ_func_0_0(36)) +
          (occ_func_0_0(23) * occ_func_0_0(25) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(10) * occ_func_0_0(36)) +
          (occ_func_0_0(3) * occ_func_0_0(0) * occ_func_0_0(31)) +
          (occ_func_0_0(25) * occ_func_0_0(30) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(3) * occ_func_0_0(32)) +
          (occ_func_0_0(10) * occ_func_0_0(0) * occ_func_0_0(37)) +
          (occ_func_0_0(29) * occ_func_0_0(24) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(6) * occ_func_0_0(34)) +
          (occ_func_0_0(7) * occ_func_0_0(0) * occ_func_0_0(35)) +
          (occ_func_0_0(27) * occ_func_0_0(26) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(11) * occ_func_0_0(39)) +
          (occ_func_0_0(2) * occ_func_0_0(0) * occ_func_0_0(32)) +
          (occ_func_0_0(22) * occ_func_0_0(29) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(7) * occ_func_0_0(21)) +
          (occ_func_0_0(6) * occ_func_0_0(0) * occ_func_0_0(20)) +
          (occ_func_0_0(40) * occ_func_0_0(41) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(4) * occ_func_0_0(19)) +
          (occ_func_0_0(9) * occ_func_0_0(0) * occ_func_0_0(21)) +
          (occ_func_0_0(42) * occ_func_0_0(40) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(2) * occ_func_0_0(28)) +
          (occ_func_0_0(11) * occ_func_0_0(0) * occ_func_0_0(34)) +
          (occ_func_0_0(33) * occ_func_0_0(27) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(8) * occ_func_0_0(42)) +
          (occ_func_0_0(5) * occ_func_0_0(0) * occ_func_0_0(41)) +
          (occ_func_0_0(19) * occ_func_0_0(20) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(5) * occ_func_0_0(22)) +
          (occ_func_0_0(8) * occ_func_0_0(0) * occ_func_0_0(24)) +
          (occ_func_0_0(39) * occ_func_0_0(37) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(6) * occ_func_0_0(26)) +
          (occ_func_0_0(7) * occ_func_0_0(0) * occ_func_0_0(27)) +
          (occ_func_0_0(35) * occ_func_0_0(34) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(8) * occ_func_0_0(20)) +
          (occ_func_0_0(5) * occ_func_0_0(0) * occ_func_0_0(19)) +
          (occ_func_0_0(41) * occ_func_0_0(42) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(1) * occ_func_0_0(31)) +
          (occ_func_0_0(12) * occ_func_0_0(0) * occ_func_0_0(38)) +
          (occ_func_0_0(30) * occ_func_0_0(23) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(9) * occ_func_0_0(25)) +
          (occ_func_0_0(4) * occ_func_0_0(0) * occ_func_0_0(23)) +
          (occ_func_0_0(36) * occ_func_0_0(38) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(11) * occ_func_0_0(29)) +
          (occ_func_0_0(2) * occ_func_0_0(0) * occ_func_0_0(22)) +
          (occ_func_0_0(32) * occ_func_0_0(39) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(12) * occ_func_0_0(35)) +
          (occ_func_0_0(1) * occ_func_0_0(0) * occ_func_0_0(28)) +
          (occ_func_0_0(26) * occ_func_0_0(33) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(3) * occ_func_0_0(24)) +
          (occ_func_0_0(10) * occ_func_0_0(0) * occ_func_0_0(29)) +
          (occ_func_0_0(37) * occ_func_0_0(32) * occ_func_0_0(0))) /
         24.0;
}
double test_Clexulator::site_eval_at_0_bfunc_3_1_1() const {
  return (((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(5) *
                occ_func_0_0(37) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(5) *
                occ_func_0_0(37))) +
          ((0.7071067812 * occ_func_0_1(8) * occ_func_0_0(0) *
                occ_func_0_0(39) +
            0.7071067812 * occ_func_0_0(8) * occ_func_0_1(0) *
                occ_func_0_0(39))) +
          ((0.7071067812 * occ_func_0_1(24) * occ_func_0_0(22) *
                occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(24) * occ_func_0_1(22) *
                occ_func_0_0(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(12) *
                occ_func_0_0(33) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(12) *
                occ_func_0_0(33))) +
          ((0.7071067812 * occ_func_0_1(1) * occ_func_0_0(0) *
                occ_func_0_0(26) +
            0.7071067812 * occ_func_0_0(1) * occ_func_0_1(0) *
                occ_func_0_0(26))) +
          ((0.7071067812 * occ_func_0_1(28) * occ_func_0_0(35) *
                occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(28) * occ_func_0_1(35) *
                occ_func_0_0(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(10) *
                occ_func_0_0(30) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(10) *
                occ_func_0_0(30))) +
          ((0.7071067812 * occ_func_0_1(3) * occ_func_0_0(0) *
                occ_func_0_0(25) +
            0.7071067812 * occ_func_0_0(3) * occ_func_0_1(0) *
                occ_func_0_0(25))) +
          ((0.7071067812 * occ_func_0_1(31) * occ_func_0_0(36) *
                occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(31) * occ_func_0_1(36) *
                occ_func_0_0(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(4) *
                occ_func_0_0(40) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(4) *
                occ_func_0_0(40))) +
          ((0.7071067812 * occ_func_0_1(9) * occ_func_0_0(0) *
                occ_func_0_0(42) +
            0.7071067812 * occ_func_0_0(9) * occ_func_0_1(0) *
                occ_func_0_0(42))) +
          ((0.7071067812 * occ_func_0_1(21) * occ_func_0_0(19) *
                occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(21) * occ_func_0_1(19) *
                occ_func_0_0(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(2) *
                occ_func_0_0(27) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(2) *
                occ_func_0_0(27))) +
          ((0.7071067812 * occ_func_0_1(11) * occ_func_0_0(0) *
                occ_func_0_0(33) +
            0.7071067812 * occ_func_0_0(11) * occ_func_0_1(0) *
                occ_func_0_0(33))) +
          ((0.7071067812 * occ_func_0_1(34) * occ_func_0_0(28) *
                occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(34) * occ_func_0_1(28) *
                occ_func_0_0(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(7) *
                occ_func_0_0(41) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(7) *
                occ_func_0_0(41))) +
          ((0.7071067812 * occ_func_0_1(6) * occ_func_0_0(0) *
                occ_func_0_0(40) +
            0.7071067812 * occ_func_0_0(6) * occ_func_0_1(0) *
                occ_func_0_0(40))) +
          ((0.7071067812 * occ_func_0_1(20) * occ_func_0_0(21) *
                occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(20) * occ_func_0_1(21) *
                occ_func_0_0(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(1) *
                occ_func_0_0(23) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(1) *
                occ_func_0_0(23))) +
          ((0.7071067812 * occ_func_0_1(12) * occ_func_0_0(0) *
                occ_func_0_0(30) +
            0.7071067812 * occ_func_0_0(12) * occ_func_0_1(0) *
                occ_func_0_0(30))) +
          ((0.7071067812 * occ_func_0_1(38) * occ_func_0_0(31) *
                occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(38) * occ_func_0_1(31) *
                occ_func_0_0(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(9) *
                occ_func_0_0(38) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(9) *
                occ_func_0_0(38))) +
          ((0.7071067812 * occ_func_0_1(4) * occ_func_0_0(0) *
                occ_func_0_0(36) +
            0.7071067812 * occ_func_0_0(4) * occ_func_0_1(0) *
                occ_func_0_0(36))) +
          ((0.7071067812 * occ_func_0_1(23) * occ_func_0_0(25) *
                occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(23) * occ_func_0_1(25) *
                occ_func_0_0(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(10) *
                occ_func_0_0(36) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(10) *
                occ_func_0_0(36))) +
          ((0.7071067812 * occ_func_0_1(3) * occ_func_0_0(0) *
                occ_func_0_0(31) +
            0.7071067812 * occ_func_0_0(3) * occ_func_0_1(0) *
                occ_func_0_0(31))) +
          ((0.7071067812 * occ_func_0_1(25) * occ_func_0_0(30) *
                occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(25) * occ_func_0_1(30) *
                occ_func_0_0(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(3) *
                occ_func_0_0(32) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(3) *
                occ_func_0_0(32))) +
          ((0.7071067812 * occ_func_0_1(10) * occ_func_0_0(0) *
                occ_func_0_0(37) +
            0.7071067812 * occ_func_0_0(10) * occ_func_0_1(0) *
                occ_func_0_0(37))) +
          ((0.7071067812 * occ_func_0_1(29) * occ_func_0_0(24) *
                occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(29) * occ_func_0_1(24) *
                occ_func_0_0(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(6) *
                occ_func_0_0(34) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(6) *
                occ_func_0_0(34))) +
          ((0.7071067812 * occ_func_0_1(7) * occ_func_0_0(0) *
                occ_func_0_0(35) +
            0.7071067812 * occ_func_0_0(7) * occ_func_0_1(0) *
                occ_func_0_0(35))) +
          ((0.7071067812 * occ_func_0_1(27) * occ_func_0_0(26) *
                occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(27) * occ_func_0_1(26) *
                occ_func_0_0(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(11) *
                occ_func_0_0(39) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(11) *
                occ_func_0_0(39))) +
          ((0.7071067812 * occ_func_0_1(2) * occ_func_0_0(0) *
                occ_func_0_0(32) +
            0.7071067812 * occ_func_0_0(2) * occ_func_0_1(0) *
                occ_func_0_0(32))) +
          ((0.7071067812 * occ_func_0_1(22) * occ_func_0_0(29) *
                occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(22) * occ_func_0_1(29) *
                occ_func_0_0(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(7) *
                occ_func_0_0(21) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(7) *
                occ_func_0_0(21))) +
          ((0.7071067812 * occ_func_0_1(6) * occ_func_0_0(0) *
                occ_func_0_0(20) +
            0.7071067812 * occ_func_0_0(6) * occ_func_0_1(0) *
                occ_func_0_0(20))) +
          ((0.7071067812 * occ_func_0_1(40) * occ_func_0_0(41) *
                occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(40) * occ_func_0_1(41) *
                occ_func_0_0(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(4) *
                occ_func_0_0(19) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(4) *
                occ_func_0_0(19))) +
          ((0.7071067812 * occ_func_0_1(9) * occ_func_0_0(0) *
                occ_func_0_0(21) +
            0.7071067812 * occ_func_0_0(9) * occ_func_0_1(0) *
                occ_func_0_0(21))) +
          ((0.7071067812 * occ_func_0_1(42) * occ_func_0_0(40) *
                occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(42) * occ_func_0_1(40) *
                occ_func_0_0(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(2) *
                occ_func_0_0(28) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(2) *
                occ_func_0_0(28))) +
          ((0.7071067812 * occ_func_0_1(11) * occ_func_0_0(0) *
                occ_func_0_0(34) +
            0.7071067812 * occ_func_0_0(11) * occ_func_0_1(0) *
                occ_func_0_0(34))) +
          ((0.7071067812 * occ_func_0_1(33) * occ_func_0_0(27) *
                occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(33) * occ_func_0_1(27) *
                occ_func_0_0(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(8) *
                occ_func_0_0(42) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(8) *
                occ_func_0_0(42))) +
          ((0.7071067812 * occ_func_0_1(5) * occ_func_0_0(0) *
                occ_func_0_0(41) +
            0.7071067812 * occ_func_0_0(5) * occ_func_0_1(0) *
                occ_func_0_0(41))) +
          ((0.7071067812 * occ_func_0_1(19) * occ_func_0_0(20) *
                occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(19) * occ_func_0_1(20) *
                occ_func_0_0(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(5) *
                occ_func_0_0(22) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(5) *
                occ_func_0_0(22))) +
          ((0.7071067812 * occ_func_0_1(8) * occ_func_0_0(0) *
                occ_func_0_0(24) +
            0.7071067812 * occ_func_0_0(8) * occ_func_0_1(0) *
                occ_func_0_0(24))) +
          ((0.7071067812 * occ_func_0_1(39) * occ_func_0_0(37) *
                occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(39) * occ_func_0_1(37) *
                occ_func_0_0(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(6) *
                occ_func_0_0(26) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(6) *
                occ_func_0_0(26))) +
          ((0.7071067812 * occ_func_0_1(7) * occ_func_0_0(0) *
                occ_func_0_0(27) +
            0.7071067812 * occ_func_0_0(7) * occ_func_0_1(0) *
                occ_func_0_0(27))) +
          ((0.7071067812 * occ_func_0_1(35) * occ_func_0_0(34) *
                occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(35) * occ_func_0_1(34) *
                occ_func_0_0(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(8) *
                occ_func_0_0(20) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(8) *
                occ_func_0_0(20))) +
          ((0.7071067812 * occ_func_0_1(5) * occ_func_0_0(0) *
                occ_func_0_0(19) +
            0.7071067812 * occ_func_0_0(5) * occ_func_0_1(0) *
                occ_func_0_0(19))) +
          ((0.7071067812 * occ_func_0_1(41) * occ_func_0_0(42) *
                occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(41) * occ_func_0_1(42) *
                occ_func_0_0(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(1) *
                occ_func_0_0(31) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(1) *
                occ_func_0_0(31))) +
          ((0.7071067812 * occ_func_0_1(12) * occ_func_0_0(0) *
                occ_func_0_0(38) +
            0.7071067812 * occ_func_0_0(12) * occ_func_0_1(0) *
                occ_func_0_0(38))) +
          ((0.7071067812 * occ_func_0_1(30) * occ_func_0_0(23) *
                occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(30) * occ_func_0_1(23) *
                occ_func_0_0(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(9) *
                occ_func_0_0(25) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(9) *
                occ_func_0_0(25))) +
          ((0.7071067812 * occ_func_0_1(4) * occ_func_0_0(0) *
                occ_func_0_0(23) +
            0.7071067812 * occ_func_0_0(4) * occ_func_0_1(0) *
                occ_func_0_0(23))) +
          ((0.7071067812 * occ_func_0_1(36) * occ_func_0_0(38) *
                occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(36) * occ_func_0_1(38) *
                occ_func_0_0(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(11) *
                occ_func_0_0(29) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(11) *
                occ_func_0_0(29))) +
          ((0.7071067812 * occ_func_0_1(2) * occ_func_0_0(0) *
                occ_func_0_0(22) +
            0.7071067812 * occ_func_0_0(2) * occ_func_0_1(0) *
                occ_func_0_0(22))) +
          ((0.7071067812 * occ_func_0_1(32) * occ_func_0_0(39) *
                occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(32) * occ_func_0_1(39) *
                occ_func_0_0(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(12) *
                occ_func_0_0(35) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(12) *
                occ_func_0_0(35))) +
          ((0.7071067812 * occ_func_0_1(1) * occ_func_0_0(0) *
                occ_func_0_0(28) +
            0.7071067812 * occ_func_0_0(1) * occ_func_0_1(0) *
                occ_func_0_0(28))) +
          ((0.7071067812 * occ_func_0_1(26) * occ_func_0_0(33) *
                occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(26) * occ_func_0_1(33) *
                occ_func_0_0(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(3) *
                occ_func_0_0(24) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(3) *
                occ_func_0_0(24))) +
          ((0.7071067812 * occ_func_0_1(10) * occ_func_0_0(0) *
                occ_func_0_0(29) +
            0.7071067812 * occ_func_0_0(10) * occ_func_0_1(0) *
                occ_func_0_0(29))) +
          ((0.7071067812 * occ_func_0_1(37) * occ_func_0_0(32) *
                occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(37) * occ_func_0_1(32) *
                occ_func_0_0(0)))) /
         24.0;
}
double test_Clexulator::site_eval_at_0_bfunc_3_1_2() const {
  return ((occ_func_0_1(0) * occ_func_0_1(5) * occ_func_0_0(37)) +
          (occ_func_0_1(8) * occ_func_0_1(0) * occ_func_0_0(39)) +
          (occ_func_0_1(24) * occ_func_0_1(22) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_1(12) * occ_func_0_0(33)) +
          (occ_func_0_1(1) * occ_func_0_1(0) * occ_func_0_0(26)) +
          (occ_func_0_1(28) * occ_func_0_1(35) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_1(10) * occ_func_0_0(30)) +
          (occ_func_0_1(3) * occ_func_0_1(0) * occ_func_0_0(25)) +
          (occ_func_0_1(31) * occ_func_0_1(36) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_1(4) * occ_func_0_0(40)) +
          (occ_func_0_1(9) * occ_func_0_1(0) * occ_func_0_0(42)) +
          (occ_func_0_1(21) * occ_func_0_1(19) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_1(2) * occ_func_0_0(27)) +
          (occ_func_0_1(11) * occ_func_0_1(0) * occ_func_0_0(33)) +
          (occ_func_0_1(34) * occ_func_0_1(28) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_1(7) * occ_func_0_0(41)) +
          (occ_func_0_1(6) * occ_func_0_1(0) * occ_func_0_0(40)) +
          (occ_func_0_1(20) * occ_func_0_1(21) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_1(1) * occ_func_0_0(23)) +
          (occ_func_0_1(12) * occ_func_0_1(0) * occ_func_0_0(30)) +
          (occ_func_0_1(38) * occ_func_0_1(31) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_1(9) * occ_func_0_0(38)) +
          (occ_func_0_1(4) * occ_func_0_1(0) * occ_func_0_0(36)) +
          (occ_func_0_1(23) * occ_func_0_1(25) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_1(10) * occ_func_0_0(36)) +
          (occ_func_0_1(3) * occ_func_0_1(0) * occ_func_0_0(31)) +
          (occ_func_0_1(25) * occ_func_0_1(30) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_1(3) * occ_func_0_0(32)) +
          (occ_func_0_1(10) * occ_func_0_1(0) * occ_func_0_0(37)) +
          (occ_func_0_1(29) * occ_func_0_1(24) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_1(6) * occ_func_0_0(34)) +
          (occ_func_0_1(7) * occ_func_0_1(0) * occ_func_0_0(35)) +
          (occ_func_0_1(27) * occ_func_0_1(26) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_1(11) * occ_func_0_0(39)) +
          (occ_func_0_1(2) * occ_func_0_1(0) * occ_func_0_0(32)) +
          (occ_func_0_1(22) * occ_func_0_1(29) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_1(7) * occ_func_0_0(21)) +
          (occ_func_0_1(6) * occ_func_0_1(0) * occ_func_0_0(20)) +
          (occ_func_0_1(40) * occ_func_0_1(41) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_1(4) * occ_func_0_0(19)) +
          (occ_func_0_1(9) * occ_func_0_1(0) * occ_func_0_0(21)) +
          (occ_func_0_1(42) * occ_func_0_1(40) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_1(2) * occ_func_0_0(28)) +
          (occ_func_0_1(11) * occ_func_0_1(0) * occ_func_0_0(34)) +
          (occ_func_0_1(33) * occ_func_0_1(27) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_1(8) * occ_func_0_0(42)) +
          (occ_func_0_1(5) * occ_func_0_1(0) * occ_func_0_0(41)) +
          (occ_func_0_1(19) * occ_func_0_1(20) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_1(5) * occ_func_0_0(22)) +
          (occ_func_0_1(8) * occ_func_0_1(0) * occ_func_0_0(24)) +
          (occ_func_0_1(39) * occ_func_0_1(37) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_1(6) * occ_func_0_0(26)) +
          (occ_func_0_1(7) * occ_func_0_1(0) * occ_func_0_0(27)) +
          (occ_func_0_1(35) * occ_func_0_1(34) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_1(8) * occ_func_0_0(20)) +
          (occ_func_0_1(5) * occ_func_0_1(0) * occ_func_0_0(19)) +
          (occ_func_0_1(41) * occ_func_0_1(42) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_1(1) * occ_func_0_0(31)) +
          (occ_func_0_1(12) * occ_func_0_1(0) * occ_func_0_0(38)) +
          (occ_func_0_1(30) * occ_func_0_1(23) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_1(9) * occ_func_0_0(25)) +
          (occ_func_0_1(4) * occ_func_0_1(0) * occ_func_0_0(23)) +
          (occ_func_0_1(36) * occ_func_0_1(38) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_1(11) * occ_func_0_0(29)) +
          (occ_func_0_1(2) * occ_func_0_1(0) * occ_func_0_0(22)) +
          (occ_func_0_1(32) * occ_func_0_1(39) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_1(12) * occ_func_0_0(35)) +
          (occ_func_0_1(1) * occ_func_0_1(0) * occ_func_0_0(28)) +
          (occ_func_0_1(26) * occ_func_0_1(33) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_1(3) * occ_func_0_0(24)) +
          (occ_func_0_1(10) * occ_func_0_1(0) * occ_func_0_0(29)) +
          (occ_func_0_1(37) * occ_func_0_1(32) * occ_func_0_0(0))) /
         24.0;
}
double test_Clexulator::site_eval_at_0_bfunc_3_1_3() const {
  return ((occ_func_0_0(0) * occ_func_0_0(5) * occ_func_0_1(37)) +
          (occ_func_0_0(8) * occ_func_0_0(0) * occ_func_0_1(39)) +
          (occ_func_0_0(24) * occ_func_0_0(22) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_0(12) * occ_func_0_1(33)) +
          (occ_func_0_0(1) * occ_func_0_0(0) * occ_func_0_1(26)) +
          (occ_func_0_0(28) * occ_func_0_0(35) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_0(10) * occ_func_0_1(30)) +
          (occ_func_0_0(3) * occ_func_0_0(0) * occ_func_0_1(25)) +
          (occ_func_0_0(31) * occ_func_0_0(36) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_0(4) * occ_func_0_1(40)) +
          (occ_func_0_0(9) * occ_func_0_0(0) * occ_func_0_1(42)) +
          (occ_func_0_0(21) * occ_func_0_0(19) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_0(2) * occ_func_0_1(27)) +
          (occ_func_0_0(11) * occ_func_0_0(0) * occ_func_0_1(33)) +
          (occ_func_0_0(34) * occ_func_0_0(28) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_0(7) * occ_func_0_1(41)) +
          (occ_func_0_0(6) * occ_func_0_0(0) * occ_func_0_1(40)) +
          (occ_func_0_0(20) * occ_func_0_0(21) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_0(1) * occ_func_0_1(23)) +
          (occ_func_0_0(12) * occ_func_0_0(0) * occ_func_0_1(30)) +
          (occ_func_0_0(38) * occ_func_0_0(31) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_0(9) * occ_func_0_1(38)) +
          (occ_func_0_0(4) * occ_func_0_0(0) * occ_func_0_1(36)) +
          (occ_func_0_0(23) * occ_func_0_0(25) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_0(10) * occ_func_0_1(36)) +
          (occ_func_0_0(3) * occ_func_0_0(0) * occ_func_0_1(31)) +
          (occ_func_0_0(25) * occ_func_0_0(30) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_0(3) * occ_func_0_1(32)) +
          (occ_func_0_0(10) * occ_func_0_0(0) * occ_func_0_1(37)) +
          (occ_func_0_0(29) * occ_func_0_0(24) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_0(6) * occ_func_0_1(34)) +
          (occ_func_0_0(7) * occ_func_0_0(0) * occ_func_0_1(35)) +
          (occ_func_0_0(27) * occ_func_0_0(26) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_0(11) * occ_func_0_1(39)) +
          (occ_func_0_0(2) * occ_func_0_0(0) * occ_func_0_1(32)) +
          (occ_func_0_0(22) * occ_func_0_0(29) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_0(7) * occ_func_0_1(21)) +
          (occ_func_0_0(6) * occ_func_0_0(0) * occ_func_0_1(20)) +
          (occ_func_0_0(40) * occ_func_0_0(41) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_0(4) * occ_func_0_1(19)) +
          (occ_func_0_0(9) * occ_func_0_0(0) * occ_func_0_1(21)) +
          (occ_func_0_0(42) * occ_func_0_0(40) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_0(2) * occ_func_0_1(28)) +
          (occ_func_0_0(11) * occ_func_0_0(0) * occ_func_0_1(34)) +
          (occ_func_0_0(33) * occ_func_0_0(27) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_0(8) * occ_func_0_1(42)) +
          (occ_func_0_0(5) * occ_func_0_0(0) * occ_func_0_1(41)) +
          (occ_func_0_0(19) * occ_func_0_0(20) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_0(5) * occ_func_0_1(22)) +
          (occ_func_0_0(8) * occ_func_0_0(0) * occ_func_0_1(24)) +
          (occ_func_0_0(39) * occ_func_0_0(37) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_0(6) * occ_func_0_1(26)) +
          (occ_func_0_0(7) * occ_func_0_0(0) * occ_func_0_1(27)) +
          (occ_func_0_0(35) * occ_func_0_0(34) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_0(8) * occ_func_0_1(20)) +
          (occ_func_0_0(5) * occ_func_0_0(0) * occ_func_0_1(19)) +
          (occ_func_0_0(41) * occ_func_0_0(42) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_0(1) * occ_func_0_1(31)) +
          (occ_func_0_0(12) * occ_func_0_0(0) * occ_func_0_1(38)) +
          (occ_func_0_0(30) * occ_func_0_0(23) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_0(9) * occ_func_0_1(25)) +
          (occ_func_0_0(4) * occ_func_0_0(0) * occ_func_0_1(23)) +
          (occ_func_0_0(36) * occ_func_0_0(38) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_0(11) * occ_func_0_1(29)) +
          (occ_func_0_0(2) * occ_func_0_0(0) * occ_func_0_1(22)) +
          (occ_func_0_0(32) * occ_func_0_0(39) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_0(12) * occ_func_0_1(35)) +
          (occ_func_0_0(1) * occ_func_0_0(0) * occ_func_0_1(28)) +
          (occ_func_0_0(26) * occ_func_0_0(33) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_0(3) * occ_func_0_1(24)) +
          (occ_func_0_0(10) * occ_func_0_0(0) * occ_func_0_1(29)) +
          (occ_func_0_0(37) * occ_func_0_0(32) * occ_func_0_1(0))) /
         24.0;
}
double test_Clexulator::site_eval_at_0_bfunc_3_1_4() const {
  return (((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(5) *
                occ_func_0_1(37) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(5) *
                occ_func_0_1(37))) +
          ((0.7071067812 * occ_func_0_1(8) * occ_func_0_0(0) *
                occ_func_0_1(39) +
            0.7071067812 * occ_func_0_0(8) * occ_func_0_1(0) *
                occ_func_0_1(39))) +
          ((0.7071067812 * occ_func_0_1(24) * occ_func_0_0(22) *
                occ_func_0_1(0) +
            0.7071067812 * occ_func_0_0(24) * occ_func_0_1(22) *
                occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(12) *
                occ_func_0_1(33) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(12) *
                occ_func_0_1(33))) +
          ((0.7071067812 * occ_func_0_1(1) * occ_func_0_0(0) *
                occ_func_0_1(26) +
            0.7071067812 * occ_func_0_0(1) * occ_func_0_1(0) *
                occ_func_0_1(26))) +
          ((0.7071067812 * occ_func_0_1(28) * occ_func_0_0(35) *
                occ_func_0_1(0) +
            0.7071067812 * occ_func_0_0(28) * occ_func_0_1(35) *
                occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(10) *
                occ_func_0_1(30) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(10) *
                occ_func_0_1(30))) +
          ((0.7071067812 * occ_func_0_1(3) * occ_func_0_0(0) *
                occ_func_0_1(25) +
            0.7071067812 * occ_func_0_0(3) * occ_func_0_1(0) *
                occ_func_0_1(25))) +
          ((0.7071067812 * occ_func_0_1(31) * occ_func_0_0(36) *
                occ_func_0_1(0) +
            0.7071067812 * occ_func_0_0(31) * occ_func_0_1(36) *
                occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(4) *
                occ_func_0_1(40) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(4) *
                occ_func_0_1(40))) +
          ((0.7071067812 * occ_func_0_1(9) * occ_func_0_0(0) *
                occ_func_0_1(42) +
            0.7071067812 * occ_func_0_0(9) * occ_func_0_1(0) *
                occ_func_0_1(42))) +
          ((0.7071067812 * occ_func_0_1(21) * occ_func_0_0(19) *
                occ_func_0_1(0) +
            0.7071067812 * occ_func_0_0(21) * occ_func_0_1(19) *
                occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(2) *
                occ_func_0_1(27) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(2) *
                occ_func_0_1(27))) +
          ((0.7071067812 * occ_func_0_1(11) * occ_func_0_0(0) *
                occ_func_0_1(33) +
            0.7071067812 * occ_func_0_0(11) * occ_func_0_1(0) *
                occ_func_0_1(33))) +
          ((0.7071067812 * occ_func_0_1(34) * occ_func_0_0(28) *
                occ_func_0_1(0) +
            0.7071067812 * occ_func_0_0(34) * occ_func_0_1(28) *
                occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(7) *
                occ_func_0_1(41) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(7) *
                occ_func_0_1(41))) +
          ((0.7071067812 * occ_func_0_1(6) * occ_func_0_0(0) *
                occ_func_0_1(40) +
            0.7071067812 * occ_func_0_0(6) * occ_func_0_1(0) *
                occ_func_0_1(40))) +
          ((0.7071067812 * occ_func_0_1(20) * occ_func_0_0(21) *
                occ_func_0_1(0) +
            0.7071067812 * occ_func_0_0(20) * occ_func_0_1(21) *
                occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(1) *
                occ_func_0_1(23) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(1) *
                occ_func_0_1(23))) +
          ((0.7071067812 * occ_func_0_1(12) * occ_func_0_0(0) *
                occ_func_0_1(30) +
            0.7071067812 * occ_func_0_0(12) * occ_func_0_1(0) *
                occ_func_0_1(30))) +
          ((0.7071067812 * occ_func_0_1(38) * occ_func_0_0(31) *
                occ_func_0_1(0) +
            0.7071067812 * occ_func_0_0(38) * occ_func_0_1(31) *
                occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(9) *
                occ_func_0_1(38) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(9) *
                occ_func_0_1(38))) +
          ((0.7071067812 * occ_func_0_1(4) * occ_func_0_0(0) *
                occ_func_0_1(36) +
            0.7071067812 * occ_func_0_0(4) * occ_func_0_1(0) *
                occ_func_0_1(36))) +
          ((0.7071067812 * occ_func_0_1(23) * occ_func_0_0(25) *
                occ_func_0_1(0) +
            0.7071067812 * occ_func_0_0(23) * occ_func_0_1(25) *
                occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(10) *
                occ_func_0_1(36) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(10) *
                occ_func_0_1(36))) +
          ((0.7071067812 * occ_func_0_1(3) * occ_func_0_0(0) *
                occ_func_0_1(31) +
            0.7071067812 * occ_func_0_0(3) * occ_func_0_1(0) *
                occ_func_0_1(31))) +
          ((0.7071067812 * occ_func_0_1(25) * occ_func_0_0(30) *
                occ_func_0_1(0) +
            0.7071067812 * occ_func_0_0(25) * occ_func_0_1(30) *
                occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(3) *
                occ_func_0_1(32) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(3) *
                occ_func_0_1(32))) +
          ((0.7071067812 * occ_func_0_1(10) * occ_func_0_0(0) *
                occ_func_0_1(37) +
            0.7071067812 * occ_func_0_0(10) * occ_func_0_1(0) *
                occ_func_0_1(37))) +
          ((0.7071067812 * occ_func_0_1(29) * occ_func_0_0(24) *
                occ_func_0_1(0) +
            0.7071067812 * occ_func_0_0(29) * occ_func_0_1(24) *
                occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(6) *
                occ_func_0_1(34) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(6) *
                occ_func_0_1(34))) +
          ((0.7071067812 * occ_func_0_1(7) * occ_func_0_0(0) *
                occ_func_0_1(35) +
            0.7071067812 * occ_func_0_0(7) * occ_func_0_1(0) *
                occ_func_0_1(35))) +
          ((0.7071067812 * occ_func_0_1(27) * occ_func_0_0(26) *
                occ_func_0_1(0) +
            0.7071067812 * occ_func_0_0(27) * occ_func_0_1(26) *
                occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(11) *
                occ_func_0_1(39) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(11) *
                occ_func_0_1(39))) +
          ((0.7071067812 * occ_func_0_1(2) * occ_func_0_0(0) *
                occ_func_0_1(32) +
            0.7071067812 * occ_func_0_0(2) * occ_func_0_1(0) *
                occ_func_0_1(32))) +
          ((0.7071067812 * occ_func_0_1(22) * occ_func_0_0(29) *
                occ_func_0_1(0) +
            0.7071067812 * occ_func_0_0(22) * occ_func_0_1(29) *
                occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(7) *
                occ_func_0_1(21) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(7) *
                occ_func_0_1(21))) +
          ((0.7071067812 * occ_func_0_1(6) * occ_func_0_0(0) *
                occ_func_0_1(20) +
            0.7071067812 * occ_func_0_0(6) * occ_func_0_1(0) *
                occ_func_0_1(20))) +
          ((0.7071067812 * occ_func_0_1(40) * occ_func_0_0(41) *
                occ_func_0_1(0) +
            0.7071067812 * occ_func_0_0(40) * occ_func_0_1(41) *
                occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(4) *
                occ_func_0_1(19) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(4) *
                occ_func_0_1(19))) +
          ((0.7071067812 * occ_func_0_1(9) * occ_func_0_0(0) *
                occ_func_0_1(21) +
            0.7071067812 * occ_func_0_0(9) * occ_func_0_1(0) *
                occ_func_0_1(21))) +
          ((0.7071067812 * occ_func_0_1(42) * occ_func_0_0(40) *
                occ_func_0_1(0) +
            0.7071067812 * occ_func_0_0(42) * occ_func_0_1(40) *
                occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(2) *
                occ_func_0_1(28) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(2) *
                occ_func_0_1(28))) +
          ((0.7071067812 * occ_func_0_1(11) * occ_func_0_0(0) *
                occ_func_0_1(34) +
            0.7071067812 * occ_func_0_0(11) * occ_func_0_1(0) *
                occ_func_0_1(34))) +
          ((0.7071067812 * occ_func_0_1(33) * occ_func_0_0(27) *
                occ_func_0_1(0) +
            0.7071067812 * occ_func_0_0(33) * occ_func_0_1(27) *
                occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(8) *
                occ_func_0_1(42) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(8) *
                occ_func_0_1(42))) +
          ((0.7071067812 * occ_func_0_1(5) * occ_func_0_0(0) *
                occ_func_0_1(41) +
            0.7071067812 * occ_func_0_0(5) * occ_func_0_1(0) *
                occ_func_0_1(41))) +
          ((0.7071067812 * occ_func_0_1(19) * occ_func_0_0(20) *
                occ_func_0_1(0) +
            0.7071067812 * occ_func_0_0(19) * occ_func_0_1(20) *
                occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(5) *
                occ_func_0_1(22) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(5) *
                occ_func_0_1(22))) +
          ((0.7071067812 * occ_func_0_1(8) * occ_func_0_0(0) *
                occ_func_0_1(24) +
            0.7071067812 * occ_func_0_0(8) * occ_func_0_1(0) *
                occ_func_0_1(24))) +
          ((0.7071067812 * occ_func_0_1(39) * occ_func_0_0(37) *
                occ_func_0_1(0) +
            0.7071067812 * occ_func_0_0(39) * occ_func_0_1(37) *
                occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(6) *
                occ_func_0_1(26) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(6) *
                occ_func_0_1(26))) +
          ((0.7071067812 * occ_func_0_1(7) * occ_func_0_0(0) *
                occ_func_0_1(27) +
            0.7071067812 * occ_func_0_0(7) * occ_func_0_1(0) *
                occ_func_0_1(27))) +
          ((0.7071067812 * occ_func_0_1(35) * occ_func_0_0(34) *
                occ_func_0_1(0) +
            0.7071067812 * occ_func_0_0(35) * occ_func_0_1(34) *
                occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(8) *
                occ_func_0_1(20) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(8) *
                occ_func_0_1(20))) +
          ((0.7071067812 * occ_func_0_1(5) * occ_func_0_0(0) *
                occ_func_0_1(19) +
            0.7071067812 * occ_func_0_0(5) * occ_func_0_1(0) *
                occ_func_0_1(19))) +
          ((0.7071067812 * occ_func_0_1(41) * occ_func_0_0(42) *
                occ_func_0_1(0) +
            0.7071067812 * occ_func_0_0(41) * occ_func_0_1(42) *
                occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(1) *
                occ_func_0_1(31) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(1) *
                occ_func_0_1(31))) +
          ((0.7071067812 * occ_func_0_1(12) * occ_func_0_0(0) *
                occ_func_0_1(38) +
            0.7071067812 * occ_func_0_0(12) * occ_func_0_1(0) *
                occ_func_0_1(38))) +
          ((0.7071067812 * occ_func_0_1(30) * occ_func_0_0(23) *
                occ_func_0_1(0) +
            0.7071067812 * occ_func_0_0(30) * occ_func_0_1(23) *
                occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(9) *
                occ_func_0_1(25) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(9) *
                occ_func_0_1(25))) +
          ((0.7071067812 * occ_func_0_1(4) * occ_func_0_0(0) *
                occ_func_0_1(23) +
            0.7071067812 * occ_func_0_0(4) * occ_func_0_1(0) *
                occ_func_0_1(23))) +
          ((0.7071067812 * occ_func_0_1(36) * occ_func_0_0(38) *
                occ_func_0_1(0) +
            0.7071067812 * occ_func_0_0(36) * occ_func_0_1(38) *
                occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(11) *
                occ_func_0_1(29) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(11) *
                occ_func_0_1(29))) +
          ((0.7071067812 * occ_func_0_1(2) * occ_func_0_0(0) *
                occ_func_0_1(22) +
            0.7071067812 * occ_func_0_0(2) * occ_func_0_1(0) *
                occ_func_0_1(22))) +
          ((0.7071067812 * occ_func_0_1(32) * occ_func_0_0(39) *
                occ_func_0_1(0) +
            0.7071067812 * occ_func_0_0(32) * occ_func_0_1(39) *
                occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(12) *
                occ_func_0_1(35) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(12) *
                occ_func_0_1(35))) +
          ((0.7071067812 * occ_func_0_1(1) * occ_func_0_0(0) *
                occ_func_0_1(28) +
            0.7071067812 * occ_func_0_0(1) * occ_func_0_1(0) *
                occ_func_0_1(28))) +
          ((0.7071067812 * occ_func_0_1(26) * occ_func_0_0(33) *
                occ_func_0_1(0) +
            0.7071067812 * occ_func_0_0(26) * occ_func_0_1(33) *
                occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(3) *
                occ_func_0_1(24) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(3) *
                occ_func_0_1(24))) +
          ((0.7071067812 * occ_func_0_1(10) * occ_func_0_0(0) *
                occ_func_0_1(29) +
            0.7071067812 * occ_func_0_0(10) * occ_func_0_1(0) *
                occ_func_0_1(29))) +
          ((0.7071067812 * occ_func_0_1(37) * occ_func_0_0(32) *
                occ_func_0_1(0) +
            0.7071067812 * occ_func_0_0(37) * occ_func_0_1(32) *
                occ_func_0_1(0)))) /
         24.0;
}
double test_Clexulator::site_eval_at_0_bfunc_3_1_5() const {
  return ((occ_func_0_1(0) * occ_func_0_1(5) * occ_func_0_1(37)) +
          (occ_func_0_1(8) * occ_func_0_1(0) * occ_func_0_1(39)) +
          (occ_func_0_1(24) * occ_func_0_1(22) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(12) * occ_func_0_1(33)) +
          (occ_func_0_1(1) * occ_func_0_1(0) * occ_func_0_1(26)) +
          (occ_func_0_1(28) * occ_func_0_1(35) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(10) * occ_func_0_1(30)) +
          (occ_func_0_1(3) * occ_func_0_1(0) * occ_func_0_1(25)) +
          (occ_func_0_1(31) * occ_func_0_1(36) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(4) * occ_func_0_1(40)) +
          (occ_func_0_1(9) * occ_func_0_1(0) * occ_func_0_1(42)) +
          (occ_func_0_1(21) * occ_func_0_1(19) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(2) * occ_func_0_1(27)) +
          (occ_func_0_1(11) * occ_func_0_1(0) * occ_func_0_1(33)) +
          (occ_func_0_1(34) * occ_func_0_1(28) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(7) * occ_func_0_1(41)) +
          (occ_func_0_1(6) * occ_func_0_1(0) * occ_func_0_1(40)) +
          (occ_func_0_1(20) * occ_func_0_1(21) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(1) * occ_func_0_1(23)) +
          (occ_func_0_1(12) * occ_func_0_1(0) * occ_func_0_1(30)) +
          (occ_func_0_1(38) * occ_func_0_1(31) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(9) * occ_func_0_1(38)) +
          (occ_func_0_1(4) * occ_func_0_1(0) * occ_func_0_1(36)) +
          (occ_func_0_1(23) * occ_func_0_1(25) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(10) * occ_func_0_1(36)) +
          (occ_func_0_1(3) * occ_func_0_1(0) * occ_func_0_1(31)) +
          (occ_func_0_1(25) * occ_func_0_1(30) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(3) * occ_func_0_1(32)) +
          (occ_func_0_1(10) * occ_func_0_1(0) * occ_func_0_1(37)) +
          (occ_func_0_1(29) * occ_func_0_1(24) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(6) * occ_func_0_1(34)) +
          (occ_func_0_1(7) * occ_func_0_1(0) * occ_func_0_1(35)) +
          (occ_func_0_1(27) * occ_func_0_1(26) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(11) * occ_func_0_1(39)) +
          (occ_func_0_1(2) * occ_func_0_1(0) * occ_func_0_1(32)) +
          (occ_func_0_1(22) * occ_func_0_1(29) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(7) * occ_func_0_1(21)) +
          (occ_func_0_1(6) * occ_func_0_1(0) * occ_func_0_1(20)) +
          (occ_func_0_1(40) * occ_func_0_1(41) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(4) * occ_func_0_1(19)) +
          (occ_func_0_1(9) * occ_func_0_1(0) * occ_func_0_1(21)) +
          (occ_func_0_1(42) * occ_func_0_1(40) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(2) * occ_func_0_1(28)) +
          (occ_func_0_1(11) * occ_func_0_1(0) * occ_func_0_1(34)) +
          (occ_func_0_1(33) * occ_func_0_1(27) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(8) * occ_func_0_1(42)) +
          (occ_func_0_1(5) * occ_func_0_1(0) * occ_func_0_1(41)) +
          (occ_func_0_1(19) * occ_func_0_1(20) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(5) * occ_func_0_1(22)) +
          (occ_func_0_1(8) * occ_func_0_1(0) * occ_func_0_1(24)) +
          (occ_func_0_1(39) * occ_func_0_1(37) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(6) * occ_func_0_1(26)) +
          (occ_func_0_1(7) * occ_func_0_1(0) * occ_func_0_1(27)) +
          (occ_func_0_1(35) * occ_func_0_1(34) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(8) * occ_func_0_1(20)) +
          (occ_func_0_1(5) * occ_func_0_1(0) * occ_func_0_1(19)) +
          (occ_func_0_1(41) * occ_func_0_1(42) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(1) * occ_func_0_1(31)) +
          (occ_func_0_1(12) * occ_func_0_1(0) * occ_func_0_1(38)) +
          (occ_func_0_1(30) * occ_func_0_1(23) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(9) * occ_func_0_1(25)) +
          (occ_func_0_1(4) * occ_func_0_1(0) * occ_func_0_1(23)) +
          (occ_func_0_1(36) * occ_func_0_1(38) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(11) * occ_func_0_1(29)) +
          (occ_func_0_1(2) * occ_func_0_1(0) * occ_func_0_1(22)) +
          (occ_func_0_1(32) * occ_func_0_1(39) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(12) * occ_func_0_1(35)) +
          (occ_func_0_1(1) * occ_func_0_1(0) * occ_func_0_1(28)) +
          (occ_func_0_1(26) * occ_func_0_1(33) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(3) * occ_func_0_1(24)) +
          (occ_func_0_1(10) * occ_func_0_1(0) * occ_func_0_1(29)) +
          (occ_func_0_1(37) * occ_func_0_1(32) * occ_func_0_1(0))) /
         24.0;
}

double test_Clexulator::delta_site_eval_at_0_bfunc_3_1_0(int occ_i,
                                                         int occ_f) const {
  return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) *
         ((occ_func_0_0(5) * occ_func_0_0(37)) +
          (occ_func_0_0(8) * occ_func_0_0(39)) +
          (occ_func_0_0(24) * occ_func_0_0(22)) +
          (occ_func_0_0(12) * occ_func_0_0(33)) +
          (occ_func_0_0(1) * occ_func_0_0(26)) +
          (occ_func_0_0(28) * occ_func_0_0(35)) +
          (occ_func_0_0(10) * occ_func_0_0(30)) +
          (occ_func_0_0(3) * occ_func_0_0(25)) +
          (occ_func_0_0(31) * occ_func_0_0(36)) +
          (occ_func_0_0(4) * occ_func_0_0(40)) +
          (occ_func_0_0(9) * occ_func_0_0(42)) +
          (occ_func_0_0(21) * occ_func_0_0(19)) +
          (occ_func_0_0(2) * occ_func_0_0(27)) +
          (occ_func_0_0(11) * occ_func_0_0(33)) +
          (occ_func_0_0(34) * occ_func_0_0(28)) +
          (occ_func_0_0(7) * occ_func_0_0(41)) +
          (occ_func_0_0(6) * occ_func_0_0(40)) +
          (occ_func_0_0(20) * occ_func_0_0(21)) +
          (occ_func_0_0(1) * occ_func_0_0(23)) +
          (occ_func_0_0(12) * occ_func_0_0(30)) +
          (occ_func_0_0(38) * occ_func_0_0(31)) +
          (occ_func_0_0(9) * occ_func_0_0(38)) +
          (occ_func_0_0(4) * occ_func_0_0(36)) +
          (occ_func_0_0(23) * occ_func_0_0(25)) +
          (occ_func_0_0(10) * occ_func_0_0(36)) +
          (occ_func_0_0(3) * occ_func_0_0(31)) +
          (occ_func_0_0(25) * occ_func_0_0(30)) +
          (occ_func_0_0(3) * occ_func_0_0(32)) +
          (occ_func_0_0(10) * occ_func_0_0(37)) +
          (occ_func_0_0(29) * occ_func_0_0(24)) +
          (occ_func_0_0(6) * occ_func_0_0(34)) +
          (occ_func_0_0(7) * occ_func_0_0(35)) +
          (occ_func_0_0(27) * occ_func_0_0(26)) +
          (occ_func_0_0(11) * occ_func_0_0(39)) +
          (occ_func_0_0(2) * occ_func_0_0(32)) +
          (occ_func_0_0(22) * occ_func_0_0(29)) +
          (occ_func_0_0(7) * occ_func_0_0(21)) +
          (occ_func_0_0(6) * occ_func_0_0(20)) +
          (occ_func_0_0(40) * occ_func_0_0(41)) +
          (occ_func_0_0(4) * occ_func_0_0(19)) +
          (occ_func_0_0(9) * occ_func_0_0(21)) +
          (occ_func_0_0(42) * occ_func_0_0(40)) +
          (occ_func_0_0(2) * occ_func_0_0(28)) +
          (occ_func_0_0(11) * occ_func_0_0(34)) +
          (occ_func_0_0(33) * occ_func_0_0(27)) +
          (occ_func_0_0(8) * occ_func_0_0(42)) +
          (occ_func_0_0(5) * occ_func_0_0(41)) +
          (occ_func_0_0(19) * occ_func_0_0(20)) +
          (occ_func_0_0(5) * occ_func_0_0(22)) +
          (occ_func_0_0(8) * occ_func_0_0(24)) +
          (occ_func_0_0(39) * occ_func_0_0(37)) +
          (occ_func_0_0(6) * occ_func_0_0(26)) +
          (occ_func_0_0(7) * occ_func_0_0(27)) +
          (occ_func_0_0(35) * occ_func_0_0(34)) +
          (occ_func_0_0(8) * occ_func_0_0(20)) +
          (occ_func_0_0(5) * occ_func_0_0(19)) +
          (occ_func_0_0(41) * occ_func_0_0(42)) +
          (occ_func_0_0(1) * occ_func_0_0(31)) +
          (occ_func_0_0(12) * occ_func_0_0(38)) +
          (occ_func_0_0(30) * occ_func_0_0(23)) +
          (occ_func_0_0(9) * occ_func_0_0(25)) +
          (occ_func_0_0(4) * occ_func_0_0(23)) +
          (occ_func_0_0(36) * occ_func_0_0(38)) +
          (occ_func_0_0(11) * occ_func_0_0(29)) +
          (occ_func_0_0(2) * occ_func_0_0(22)) +
          (occ_func_0_0(32) * occ_func_0_0(39)) +
          (occ_func_0_0(12) * occ_func_0_0(35)) +
          (occ_func_0_0(1) * occ_func_0_0(28)) +
          (occ_func_0_0(26) * occ_func_0_0(33)) +
          (occ_func_0_0(3) * occ_func_0_0(24)) +
          (occ_func_0_0(10) * occ_func_0_0(29)) +
          (occ_func_0_0(37) * occ_func_0_0(32))) /
         24.0;
}
double test_Clexulator::delta_site_eval_at_0_bfunc_3_1_1(int occ_i,
                                                         int occ_f) const {
  return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) *
             ((0.7071067812 * occ_func_0_1(5) * occ_func_0_0(37)) +
              (0.7071067812 * occ_func_0_1(8) * occ_func_0_0(39)) +
              ((0.7071067812 * occ_func_0_1(24) * occ_func_0_0(22) +
                0.7071067812 * occ_func_0_0(24) * occ_func_0_1(22))) +
              (0.7071067812 * occ_func_0_1(12) * occ_func_0_0(33)) +
              (0.7071067812 * occ_func_0_1(1) * occ_func_0_0(26)) +
              ((0.7071067812 * occ_func_0_1(28) * occ_func_0_0(35) +
                0.7071067812 * occ_func_0_0(28) * occ_func_0_1(35))) +
              (0.7071067812 * occ_func_0_1(10) * occ_func_0_0(30)) +
              (0.7071067812 * occ_func_0_1(3) * occ_func_0_0(25)) +
              ((0.7071067812 * occ_func_0_1(31) * occ_func_0_0(36) +
                0.7071067812 * occ_func_0_0(31) * occ_func_0_1(36))) +
              (0.7071067812 * occ_func_0_1(4) * occ_func_0_0(40)) +
              (0.7071067812 * occ_func_0_1(9) * occ_func_0_0(42)) +
              ((0.7071067812 * occ_func_0_1(21) * occ_func_0_0(19) +
                0.7071067812 * occ_func_0_0(21) * occ_func_0_1(19))) +
              (0.7071067812 * occ_func_0_1(2) * occ_func_0_0(27)) +
              (0.7071067812 * occ_func_0_1(11) * occ_func_0_0(33)) +
              ((0.7071067812 * occ_func_0_1(34) * occ_func_0_0(28) +
                0.7071067812 * occ_func_0_0(34) * occ_func_0_1(28))) +
              (0.7071067812 * occ_func_0_1(7) * occ_func_0_0(41)) +
              (0.7071067812 * occ_func_0_1(6) * occ_func_0_0(40)) +
              ((0.7071067812 * occ_func_0_1(20) * occ_func_0_0(21) +
                0.7071067812 * occ_func_0_0(20) * occ_func_0_1(21))) +
              (0.7071067812 * occ_func_0_1(1) * occ_func_0_0(23)) +
              (0.7071067812 * occ_func_0_1(12) * occ_func_0_0(30)) +
              ((0.7071067812 * occ_func_0_1(38) * occ_func_0_0(31) +
                0.7071067812 * occ_func_0_0(38) * occ_func_0_1(31))) +
              (0.7071067812 * occ_func_0_1(9) * occ_func_0_0(38)) +
              (0.7071067812 * occ_func_0_1(4) * occ_func_0_0(36)) +
              ((0.7071067812 * occ_func_0_1(23) * occ_func_0_0(25) +
                0.7071067812 * occ_func_0_0(23) * occ_func_0_1(25))) +
              (0.7071067812 * occ_func_0_1(10) * occ_func_0_0(36)) +
              (0.7071067812 * occ_func_0_1(3) * occ_func_0_0(31)) +
              ((0.7071067812 * occ_func_0_1(25) * occ_func_0_0(30) +
                0.7071067812 * occ_func_0_0(25) * occ_func_0_1(30))) +
              (0.7071067812 * occ_func_0_1(3) * occ_func_0_0(32)) +
              (0.7071067812 * occ_func_0_1(10) * occ_func_0_0(37)) +
              ((0.7071067812 * occ_func_0_1(29) * occ_func_0_0(24) +
                0.7071067812 * occ_func_0_0(29) * occ_func_0_1(24))) +
              (0.7071067812 * occ_func_0_1(6) * occ_func_0_0(34)) +
              (0.7071067812 * occ_func_0_1(7) * occ_func_0_0(35)) +
              ((0.7071067812 * occ_func_0_1(27) * occ_func_0_0(26) +
                0.7071067812 * occ_func_0_0(27) * occ_func_0_1(26))) +
              (0.7071067812 * occ_func_0_1(11) * occ_func_0_0(39)) +
              (0.7071067812 * occ_func_0_1(2) * occ_func_0_0(32)) +
              ((0.7071067812 * occ_func_0_1(22) * occ_func_0_0(29) +
                0.7071067812 * occ_func_0_0(22) * occ_func_0_1(29))) +
              (0.7071067812 * occ_func_0_1(7) * occ_func_0_0(21)) +
              (0.7071067812 * occ_func_0_1(6) * occ_func_0_0(20)) +
              ((0.7071067812 * occ_func_0_1(40) * occ_func_0_0(41) +
                0.7071067812 * occ_func_0_0(40) * occ_func_0_1(41))) +
              (0.7071067812 * occ_func_0_1(4) * occ_func_0_0(19)) +
              (0.7071067812 * occ_func_0_1(9) * occ_func_0_0(21)) +
              ((0.7071067812 * occ_func_0_1(42) * occ_func_0_0(40) +
                0.7071067812 * occ_func_0_0(42) * occ_func_0_1(40))) +
              (0.7071067812 * occ_func_0_1(2) * occ_func_0_0(28)) +
              (0.7071067812 * occ_func_0_1(11) * occ_func_0_0(34)) +
              ((0.7071067812 * occ_func_0_1(33) * occ_func_0_0(27) +
                0.7071067812 * occ_func_0_0(33) * occ_func_0_1(27))) +
              (0.7071067812 * occ_func_0_1(8) * occ_func_0_0(42)) +
              (0.7071067812 * occ_func_0_1(5) * occ_func_0_0(41)) +
              ((0.7071067812 * occ_func_0_1(19) * occ_func_0_0(20) +
                0.7071067812 * occ_func_0_0(19) * occ_func_0_1(20))) +
              (0.7071067812 * occ_func_0_1(5) * occ_func_0_0(22)) +
              (0.7071067812 * occ_func_0_1(8) * occ_func_0_0(24)) +
              ((0.7071067812 * occ_func_0_1(39) * occ_func_0_0(37) +
                0.7071067812 * occ_func_0_0(39) * occ_func_0_1(37))) +
              (0.7071067812 * occ_func_0_1(6) * occ_func_0_0(26)) +
              (0.7071067812 * occ_func_0_1(7) * occ_func_0_0(27)) +
              ((0.7071067812 * occ_func_0_1(35) * occ_func_0_0(34) +
                0.7071067812 * occ_func_0_0(35) * occ_func_0_1(34))) +
              (0.7071067812 * occ_func_0_1(8) * occ_func_0_0(20)) +
              (0.7071067812 * occ_func_0_1(5) * occ_func_0_0(19)) +
              ((0.7071067812 * occ_func_0_1(41) * occ_func_0_0(42) +
                0.7071067812 * occ_func_0_0(41) * occ_func_0_1(42))) +
              (0.7071067812 * occ_func_0_1(1) * occ_func_0_0(31)) +
              (0.7071067812 * occ_func_0_1(12) * occ_func_0_0(38)) +
              ((0.7071067812 * occ_func_0_1(30) * occ_func_0_0(23) +
                0.7071067812 * occ_func_0_0(30) * occ_func_0_1(23))) +
              (0.7071067812 * occ_func_0_1(9) * occ_func_0_0(25)) +
              (0.7071067812 * occ_func_0_1(4) * occ_func_0_0(23)) +
              ((0.7071067812 * occ_func_0_1(36) * occ_func_0_0(38) +
                0.7071067812 * occ_func_0_0(36) * occ_func_0_1(38))) +
              (0.7071067812 * occ_func_0_1(11) * occ_func_0_0(29)) +
              (0.7071067812 * occ_func_0_1(2) * occ_func_0_0(22)) +
              ((0.7071067812 * occ_func_0_1(32) * occ_func_0_0(39) +
                0.7071067812 * occ_func_0_0(32) * occ_func_0_1(39))) +
              (0.7071067812 * occ_func_0_1(12) * occ_func_0_0(35)) +
              (0.7071067812 * occ_func_0_1(1) * occ_func_0_0(28)) +
              ((0.7071067812 * occ_func_0_1(26) * occ_func_0_0(33) +
                0.7071067812 * occ_func_0_0(26) * occ_func_0_1(33))) +
              (0.7071067812 * occ_func_0_1(3) * occ_func_0_0(24)) +
              (0.7071067812 * occ_func_0_1(10) * occ_func_0_0(29)) +
              ((0.7071067812 * occ_func_0_1(37) * occ_func_0_0(32) +
                0.7071067812 * occ_func_0_0(37) * occ_func_0_1(32)))) /
             24.0 +
         (m_occ_func_0_1[occ_f] - m_occ_func_0_1[occ_i]) *
             ((0.7071067812 * occ_func_0_0(5) * occ_func_0_0(37)) +
              (0.7071067812 * occ_func_0_0(8) * occ_func_0_0(39)) +
              (0.7071067812 * occ_func_0_0(12) * occ_func_0_0(33)) +
              (0.7071067812 * occ_func_0_0(1) * occ_func_0_0(26)) +
              (0.7071067812 * occ_func_0_0(10) * occ_func_0_0(30)) +
              (0.7071067812 * occ_func_0_0(3) * occ_func_0_0(25)) +
              (0.7071067812 * occ_func_0_0(4) * occ_func_0_0(40)) +
              (0.7071067812 * occ_func_0_0(9) * occ_func_0_0(42)) +
              (0.7071067812 * occ_func_0_0(2) * occ_func_0_0(27)) +
              (0.7071067812 * occ_func_0_0(11) * occ_func_0_0(33)) +
              (0.7071067812 * occ_func_0_0(7) * occ_func_0_0(41)) +
              (0.7071067812 * occ_func_0_0(6) * occ_func_0_0(40)) +
              (0.7071067812 * occ_func_0_0(1) * occ_func_0_0(23)) +
              (0.7071067812 * occ_func_0_0(12) * occ_func_0_0(30)) +
              (0.7071067812 * occ_func_0_0(9) * occ_func_0_0(38)) +
              (0.7071067812 * occ_func_0_0(4) * occ_func_0_0(36)) +
              (0.7071067812 * occ_func_0_0(10) * occ_func_0_0(36)) +
              (0.7071067812 * occ_func_0_0(3) * occ_func_0_0(31)) +
              (0.7071067812 * occ_func_0_0(3) * occ_func_0_0(32)) +
              (0.7071067812 * occ_func_0_0(10) * occ_func_0_0(37)) +
              (0.7071067812 * occ_func_0_0(6) * occ_func_0_0(34)) +
              (0.7071067812 * occ_func_0_0(7) * occ_func_0_0(35)) +
              (0.7071067812 * occ_func_0_0(11) * occ_func_0_0(39)) +
              (0.7071067812 * occ_func_0_0(2) * occ_func_0_0(32)) +
              (0.7071067812 * occ_func_0_0(7) * occ_func_0_0(21)) +
              (0.7071067812 * occ_func_0_0(6) * occ_func_0_0(20)) +
              (0.7071067812 * occ_func_0_0(4) * occ_func_0_0(19)) +
              (0.7071067812 * occ_func_0_0(9) * occ_func_0_0(21)) +
              (0.7071067812 * occ_func_0_0(2) * occ_func_0_0(28)) +
              (0.7071067812 * occ_func_0_0(11) * occ_func_0_0(34)) +
              (0.7071067812 * occ_func_0_0(8) * occ_func_0_0(42)) +
              (0.7071067812 * occ_func_0_0(5) * occ_func_0_0(41)) +
              (0.7071067812 * occ_func_0_0(5) * occ_func_0_0(22)) +
              (0.7071067812 * occ_func_0_0(8) * occ_func_0_0(24)) +
              (0.7071067812 * occ_func_0_0(6) * occ_func_0_0(26)) +
              (0.7071067812 * occ_func_0_0(7) * occ_func_0_0(27)) +
              (0.7071067812 * occ_func_0_0(8) * occ_func_0_0(20)) +
              (0.7071067812 * occ_func_0_0(5) * occ_func_0_0(19)) +
              (0.7071067812 * occ_func_0_0(1) * occ_func_0_0(31)) +
              (0.7071067812 * occ_func_0_0(12) * occ_func_0_0(38)) +
              (0.7071067812 * occ_func_0_0(9) * occ_func_0_0(25)) +
              (0.7071067812 * occ_func_0_0(4) * occ_func_0_0(23)) +
              (0.7071067812 * occ_func_0_0(11) * occ_func_0_0(29)) +
              (0.7071067812 * occ_func_0_0(2) * occ_func_0_0(22)) +
              (0.7071067812 * occ_func_0_0(12) * occ_func_0_0(35)) +
              (0.7071067812 * occ_func_0_0(1) * occ_func_0_0(28)) +
              (0.7071067812 * occ_func_0_0(3) * occ_func_0_0(24)) +
              (0.7071067812 * occ_func_0_0(10) * occ_func_0_0(29))) /
             24.0;
}
double test_Clexulator::delta_site_eval_at_0_bfunc_3_1_2(int occ_i,
                                                         int occ_f) const {
  return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) *
             ((occ_func_0_1(24) * occ_func_0_1(22)) +
              (occ_func_0_1(28) * occ_func_0_1(35)) +
              (occ_func_0_1(31) * occ_func_0_1(36)) +
              (occ_func_0_1(21) * occ_func_0_1(19)) +
              (occ_func_0_1(34) * occ_func_0_1(28)) +
              (occ_func_0_1(20) * occ_func_0_1(21)) +
              (occ_func_0_1(38) * occ_func_0_1(31)) +
              (occ_func_0_1(23) * occ_func_0_1(25)) +
              (occ_func_0_1(25) * occ_func_0_1(30)) +
              (occ_func_0_1(29) * occ_func_0_1(24)) +
              (occ_func_0_1(27) * occ_func_0_1(26)) +
              (occ_func_0_1(22) * occ_func_0_1(29)) +
              (occ_func_0_1(40) * occ_func_0_1(41)) +
              (occ_func_0_1(42) * occ_func_0_1(40)) +
              (occ_func_0_1(33) * occ_func_0_1(27)) +
              (occ_func_0_1(19) * occ_func_0_1(20)) +
              (occ_func_0_1(39) * occ_func_0_1(37)) +
              (occ_func_0_1(35) * occ_func_0_1(34)) +
              (occ_func_0_1(41) * occ_func_0_1(42)) +
              (occ_func_0_1(30) * occ_func_0_1(23)) +
              (occ_func_0_1(36) * occ_func_0_1(38)) +
              (occ_func_0_1(32) * occ_func_0_1(39)) +
              (occ_func_0_1(26) * occ_func_0_1(33)) +
              (occ_func_0_1(37) * occ_func_0_1(32))) /
             24.0 +
         (m_occ_func_0_1[occ_f] - m_occ_func_0_1[occ_i]) *
             ((occ_func_0_1(5) * occ_func_0_0(37)) +
              (occ_func_0_1(8) * occ_func_0_0(39)) +
              (occ_func_0_1(12) * occ_func_0_0(33)) +
              (occ_func_0_1(1) * occ_func_0_0(26)) +
              (occ_func_0_1(10) * occ_func_0_0(30)) +
              (occ_func_0_1(3) * occ_func_0_0(25)) +
              (occ_func_0_1(4) * occ_func_0_0(40)) +
              (occ_func_0_1(9) * occ_func_0_0(42)) +
              (occ_func_0_1(2) * occ_func_0_0(27)) +
              (occ_func_0_1(11) * occ_func_0_0(33)) +
              (occ_func_0_1(7) * occ_func_0_0(41)) +
              (occ_func_0_1(6) * occ_func_0_0(40)) +
              (occ_func_0_1(1) * occ_func_0_0(23)) +
              (occ_func_0_1(12) * occ_func_0_0(30)) +
              (occ_func_0_1(9) * occ_func_0_0(38)) +
              (occ_func_0_1(4) * occ_func_0_0(36)) +
              (occ_func_0_1(10) * occ_func_0_0(36)) +
              (occ_func_0_1(3) * occ_func_0_0(31)) +
              (occ_func_0_1(3) * occ_func_0_0(32)) +
              (occ_func_0_1(10) * occ_func_0_0(37)) +
              (occ_func_0_1(6) * occ_func_0_0(34)) +
              (occ_func_0_1(7) * occ_func_0_0(35)) +
              (occ_func_0_1(11) * occ_func_0_0(39)) +
              (occ_func_0_1(2) * occ_func_0_0(32)) +
              (occ_func_0_1(7) * occ_func_0_0(21)) +
              (occ_func_0_1(6) * occ_func_0_0(20)) +
              (occ_func_0_1(4) * occ_func_0_0(19)) +
              (occ_func_0_1(9) * occ_func_0_0(21)) +
              (occ_func_0_1(2) * occ_func_0_0(28)) +
              (occ_func_0_1(11) * occ_func_0_0(34)) +
              (occ_func_0_1(8) * occ_func_0_0(42)) +
              (occ_func_0_1(5) * occ_func_0_0(41)) +
              (occ_func_0_1(5) * occ_func_0_0(22)) +
              (occ_func_0_1(8) * occ_func_0_0(24)) +
              (occ_func_0_1(6) * occ_func_0_0(26)) +
              (occ_func_0_1(7) * occ_func_0_0(27)) +
              (occ_func_0_1(8) * occ_func_0_0(20)) +
              (occ_func_0_1(5) * occ_func_0_0(19)) +
              (occ_func_0_1(1) * occ_func_0_0(31)) +
              (occ_func_0_1(12) * occ_func_0_0(38)) +
              (occ_func_0_1(9) * occ_func_0_0(25)) +
              (occ_func_0_1(4) * occ_func_0_0(23)) +
              (occ_func_0_1(11) * occ_func_0_0(29)) +
              (occ_func_0_1(2) * occ_func_0_0(22)) +
              (occ_func_0_1(12) * occ_func_0_0(35)) +
              (occ_func_0_1(1) * occ_func_0_0(28)) +
              (occ_func_0_1(3) * occ_func_0_0(24)) +
              (occ_func_0_1(10) * occ_func_0_0(29))) /
             24.0;
}
double test_Clexulator::delta_site_eval_at_0_bfunc_3_1_3(int occ_i,
                                                         int occ_f) const {
  return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) *
             ((occ_func_0_0(5) * occ_func_0_1(37)) +
              (occ_func_0_0(8) * occ_func_0_1(39)) +
              (occ_func_0_0(12) * occ_func_0_1(33)) +
              (occ_func_0_0(1) * occ_func_0_1(26)) +
              (occ_func_0_0(10) * occ_func_0_1(30)) +
              (occ_func_0_0(3) * occ_func_0_1(25)) +
              (occ_func_0_0(4) * occ_func_0_1(40)) +
              (occ_func_0_0(9) * occ_func_0_1(42)) +
              (occ_func_0_0(2) * occ_func_0_1(27)) +
              (occ_func_0_0(11) * occ_func_0_1(33)) +
              (occ_func_0_0(7) * occ_func_0_1(41)) +
              (occ_func_0_0(6) * occ_func_0_1(40)) +
              (occ_func_0_0(1) * occ_func_0_1(23)) +
              (occ_func_0_0(12) * occ_func_0_1(30)) +
              (occ_func_0_0(9) * occ_func_0_1(38)) +
              (occ_func_0_0(4) * occ_func_0_1(36)) +
              (occ_func_0_0(10) * occ_func_0_1(36)) +
              (occ_func_0_0(3) * occ_func_0_1(31)) +
              (occ_func_0_0(3) * occ_func_0_1(32)) +
              (occ_func_0_0(10) * occ_func_0_1(37)) +
              (occ_func_0_0(6) * occ_func_0_1(34)) +
              (occ_func_0_0(7) * occ_func_0_1(35)) +
              (occ_func_0_0(11) * occ_func_0_1(39)) +
              (occ_func_0_0(2) * occ_func_0_1(32)) +
              (occ_func_0_0(7) * occ_func_0_1(21)) +
              (occ_func_0_0(6) * occ_func_0_1(20)) +
              (occ_func_0_0(4) * occ_func_0_1(19)) +
              (occ_func_0_0(9) * occ_func_0_1(21)) +
              (occ_func_0_0(2) * occ_func_0_1(28)) +
              (occ_func_0_0(11) * occ_func_0_1(34)) +
              (occ_func_0_0(8) * occ_func_0_1(42)) +
              (occ_func_0_0(5) * occ_func_0_1(41)) +
              (occ_func_0_0(5) * occ_func_0_1(22)) +
              (occ_func_0_0(8) * occ_func_0_1(24)) +
              (occ_func_0_0(6) * occ_func_0_1(26)) +
              (occ_func_0_0(7) * occ_func_0_1(27)) +
              (occ_func_0_0(8) * occ_func_0_1(20)) +
              (occ_func_0_0(5) * occ_func_0_1(19)) +
              (occ_func_0_0(1) * occ_func_0_1(31)) +
              (occ_func_0_0(12) * occ_func_0_1(38)) +
              (occ_func_0_0(9) * occ_func_0_1(25)) +
              (occ_func_0_0(4) * occ_func_0_1(23)) +
              (occ_func_0_0(11) * occ_func_0_1(29)) +
              (occ_func_0_0(2) * occ_func_0_1(22)) +
              (occ_func_0_0(12) * occ_func_0_1(35)) +
              (occ_func_0_0(1) * occ_func_0_1(28)) +
              (occ_func_0_0(3) * occ_func_0_1(24)) +
              (occ_func_0_0(10) * occ_func_0_1(29))) /
             24.0 +
         (m_occ_func_0_1[occ_f] - m_occ_func_0_1[occ_i]) *
             ((occ_func_0_0(24) * occ_func_0_0(22)) +
              (occ_func_0_0(28) * occ_func_0_0(35)) +
              (occ_func_0_0(31) * occ_func_0_0(36)) +
              (occ_func_0_0(21) * occ_func_0_0(19)) +
              (occ_func_0_0(34) * occ_func_0_0(28)) +
              (occ_func_0_0(20) * occ_func_0_0(21)) +
              (occ_func_0_0(38) * occ_func_0_0(31)) +
              (occ_func_0_0(23) * occ_func_0_0(25)) +
              (occ_func_0_0(25) * occ_func_0_0(30)) +
              (occ_func_0_0(29) * occ_func_0_0(24)) +
              (occ_func_0_0(27) * occ_func_0_0(26)) +
              (occ_func_0_0(22) * occ_func_0_0(29)) +
              (occ_func_0_0(40) * occ_func_0_0(41)) +
              (occ_func_0_0(42) * occ_func_0_0(40)) +
              (occ_func_0_0(33) * occ_func_0_0(27)) +
              (occ_func_0_0(19) * occ_func_0_0(20)) +
              (occ_func_0_0(39) * occ_func_0_0(37)) +
              (occ_func_0_0(35) * occ_func_0_0(34)) +
              (occ_func_0_0(41) * occ_func_0_0(42)) +
              (occ_func_0_0(30) * occ_func_0_0(23)) +
              (occ_func_0_0(36) * occ_func_0_0(38)) +
              (occ_func_0_0(32) * occ_func_0_0(39)) +
              (occ_func_0_0(26) * occ_func_0_0(33)) +
              (occ_func_0_0(37) * occ_func_0_0(32))) /
             24.0;
}
double test_Clexulator::delta_site_eval_at_0_bfunc_3_1_4(int occ_i,
                                                         int occ_f) const {
  return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) *
             ((0.7071067812 * occ_func_0_1(5) * occ_func_0_1(37)) +
              (0.7071067812 * occ_func_0_1(8) * occ_func_0_1(39)) +
              (0.7071067812 * occ_func_0_1(12) * occ_func_0_1(33)) +
              (0.7071067812 * occ_func_0_1(1) * occ_func_0_1(26)) +
              (0.7071067812 * occ_func_0_1(10) * occ_func_0_1(30)) +
              (0.7071067812 * occ_func_0_1(3) * occ_func_0_1(25)) +
              (0.7071067812 * occ_func_0_1(4) * occ_func_0_1(40)) +
              (0.7071067812 * occ_func_0_1(9) * occ_func_0_1(42)) +
              (0.7071067812 * occ_func_0_1(2) * occ_func_0_1(27)) +
              (0.7071067812 * occ_func_0_1(11) * occ_func_0_1(33)) +
              (0.7071067812 * occ_func_0_1(7) * occ_func_0_1(41)) +
              (0.7071067812 * occ_func_0_1(6) * occ_func_0_1(40)) +
              (0.7071067812 * occ_func_0_1(1) * occ_func_0_1(23)) +
              (0.7071067812 * occ_func_0_1(12) * occ_func_0_1(30)) +
              (0.7071067812 * occ_func_0_1(9) * occ_func_0_1(38)) +
              (0.7071067812 * occ_func_0_1(4) * occ_func_0_1(36)) +
              (0.7071067812 * occ_func_0_1(10) * occ_func_0_1(36)) +
              (0.7071067812 * occ_func_0_1(3) * occ_func_0_1(31)) +
              (0.7071067812 * occ_func_0_1(3) * occ_func_0_1(32)) +
              (0.7071067812 * occ_func_0_1(10) * occ_func_0_1(37)) +
              (0.7071067812 * occ_func_0_1(6) * occ_func_0_1(34)) +
              (0.7071067812 * occ_func_0_1(7) * occ_func_0_1(35)) +
              (0.7071067812 * occ_func_0_1(11) * occ_func_0_1(39)) +
              (0.7071067812 * occ_func_0_1(2) * occ_func_0_1(32)) +
              (0.7071067812 * occ_func_0_1(7) * occ_func_0_1(21)) +
              (0.7071067812 * occ_func_0_1(6) * occ_func_0_1(20)) +
              (0.7071067812 * occ_func_0_1(4) * occ_func_0_1(19)) +
              (0.7071067812 * occ_func_0_1(9) * occ_func_0_1(21)) +
              (0.7071067812 * occ_func_0_1(2) * occ_func_0_1(28)) +
              (0.7071067812 * occ_func_0_1(11) * occ_func_0_1(34)) +
              (0.7071067812 * occ_func_0_1(8) * occ_func_0_1(42)) +
              (0.7071067812 * occ_func_0_1(5) * occ_func_0_1(41)) +
              (0.7071067812 * occ_func_0_1(5) * occ_func_0_1(22)) +
              (0.7071067812 * occ_func_0_1(8) * occ_func_0_1(24)) +
              (0.7071067812 * occ_func_0_1(6) * occ_func_0_1(26)) +
              (0.7071067812 * occ_func_0_1(7) * occ_func_0_1(27)) +
              (0.7071067812 * occ_func_0_1(8) * occ_func_0_1(20)) +
              (0.7071067812 * occ_func_0_1(5) * occ_func_0_1(19)) +
              (0.7071067812 * occ_func_0_1(1) * occ_func_0_1(31)) +
              (0.7071067812 * occ_func_0_1(12) * occ_func_0_1(38)) +
              (0.7071067812 * occ_func_0_1(9) * occ_func_0_1(25)) +
              (0.7071067812 * occ_func_0_1(4) * occ_func_0_1(23)) +
              (0.7071067812 * occ_func_0_1(11) * occ_func_0_1(29)) +
              (0.7071067812 * occ_func_0_1(2) * occ_func_0_1(22)) +
              (0.7071067812 * occ_func_0_1(12) * occ_func_0_1(35)) +
              (0.7071067812 * occ_func_0_1(1) * occ_func_0_1(28)) +
              (0.7071067812 * occ_func_0_1(3) * occ_func_0_1(24)) +
              (0.7071067812 * occ_func_0_1(10) * occ_func_0_1(29))) /
             24.0 +
         (m_occ_func_0_1[occ_f] - m_occ_func_0_1[occ_i]) *
             ((0.7071067812 * occ_func_0_0(5) * occ_func_0_1(37)) +
              (0.7071067812 * occ_func_0_0(8) * occ_func_0_1(39)) +
              ((0.7071067812 * occ_func_0_1(24) * occ_func_0_0(22) +
                0.7071067812 * occ_func_0_0(24) * occ_func_0_1(22))) +
              (0.7071067812 * occ_func_0_0(12) * occ_func_0_1(33)) +
              (0.7071067812 * occ_func_0_0(1) * occ_func_0_1(26)) +
              ((0.7071067812 * occ_func_0_1(28) * occ_func_0_0(35) +
                0.7071067812 * occ_func_0_0(28) * occ_func_0_1(35))) +
              (0.7071067812 * occ_func_0_0(10) * occ_func_0_1(30)) +
              (0.7071067812 * occ_func_0_0(3) * occ_func_0_1(25)) +
              ((0.7071067812 * occ_func_0_1(31) * occ_func_0_0(36) +
                0.7071067812 * occ_func_0_0(31) * occ_func_0_1(36))) +
              (0.7071067812 * occ_func_0_0(4) * occ_func_0_1(40)) +
              (0.7071067812 * occ_func_0_0(9) * occ_func_0_1(42)) +
              ((0.7071067812 * occ_func_0_1(21) * occ_func_0_0(19) +
                0.7071067812 * occ_func_0_0(21) * occ_func_0_1(19))) +
              (0.7071067812 * occ_func_0_0(2) * occ_func_0_1(27)) +
              (0.7071067812 * occ_func_0_0(11) * occ_func_0_1(33)) +
              ((0.7071067812 * occ_func_0_1(34) * occ_func_0_0(28) +
                0.7071067812 * occ_func_0_0(34) * occ_func_0_1(28))) +
              (0.7071067812 * occ_func_0_0(7) * occ_func_0_1(41)) +
              (0.7071067812 * occ_func_0_0(6) * occ_func_0_1(40)) +
              ((0.7071067812 * occ_func_0_1(20) * occ_func_0_0(21) +
                0.7071067812 * occ_func_0_0(20) * occ_func_0_1(21))) +
              (0.7071067812 * occ_func_0_0(1) * occ_func_0_1(23)) +
              (0.7071067812 * occ_func_0_0(12) * occ_func_0_1(30)) +
              ((0.7071067812 * occ_func_0_1(38) * occ_func_0_0(31) +
                0.7071067812 * occ_func_0_0(38) * occ_func_0_1(31))) +
              (0.7071067812 * occ_func_0_0(9) * occ_func_0_1(38)) +
              (0.7071067812 * occ_func_0_0(4) * occ_func_0_1(36)) +
              ((0.7071067812 * occ_func_0_1(23) * occ_func_0_0(25) +
                0.7071067812 * occ_func_0_0(23) * occ_func_0_1(25))) +
              (0.7071067812 * occ_func_0_0(10) * occ_func_0_1(36)) +
              (0.7071067812 * occ_func_0_0(3) * occ_func_0_1(31)) +
              ((0.7071067812 * occ_func_0_1(25) * occ_func_0_0(30) +
                0.7071067812 * occ_func_0_0(25) * occ_func_0_1(30))) +
              (0.7071067812 * occ_func_0_0(3) * occ_func_0_1(32)) +
              (0.7071067812 * occ_func_0_0(10) * occ_func_0_1(37)) +
              ((0.7071067812 * occ_func_0_1(29) * occ_func_0_0(24) +
                0.7071067812 * occ_func_0_0(29) * occ_func_0_1(24))) +
              (0.7071067812 * occ_func_0_0(6) * occ_func_0_1(34)) +
              (0.7071067812 * occ_func_0_0(7) * occ_func_0_1(35)) +
              ((0.7071067812 * occ_func_0_1(27) * occ_func_0_0(26) +
                0.7071067812 * occ_func_0_0(27) * occ_func_0_1(26))) +
              (0.7071067812 * occ_func_0_0(11) * occ_func_0_1(39)) +
              (0.7071067812 * occ_func_0_0(2) * occ_func_0_1(32)) +
              ((0.7071067812 * occ_func_0_1(22) * occ_func_0_0(29) +
                0.7071067812 * occ_func_0_0(22) * occ_func_0_1(29))) +
              (0.7071067812 * occ_func_0_0(7) * occ_func_0_1(21)) +
              (0.7071067812 * occ_func_0_0(6) * occ_func_0_1(20)) +
              ((0.7071067812 * occ_func_0_1(40) * occ_func_0_0(41) +
                0.7071067812 * occ_func_0_0(40) * occ_func_0_1(41))) +
              (0.7071067812 * occ_func_0_0(4) * occ_func_0_1(19)) +
              (0.7071067812 * occ_func_0_0(9) * occ_func_0_1(21)) +
              ((0.7071067812 * occ_func_0_1(42) * occ_func_0_0(40) +
                0.7071067812 * occ_func_0_0(42) * occ_func_0_1(40))) +
              (0.7071067812 * occ_func_0_0(2) * occ_func_0_1(28)) +
              (0.7071067812 * occ_func_0_0(11) * occ_func_0_1(34)) +
              ((0.7071067812 * occ_func_0_1(33) * occ_func_0_0(27) +
                0.7071067812 * occ_func_0_0(33) * occ_func_0_1(27))) +
              (0.7071067812 * occ_func_0_0(8) * occ_func_0_1(42)) +
              (0.7071067812 * occ_func_0_0(5) * occ_func_0_1(41)) +
              ((0.7071067812 * occ_func_0_1(19) * occ_func_0_0(20) +
                0.7071067812 * occ_func_0_0(19) * occ_func_0_1(20))) +
              (0.7071067812 * occ_func_0_0(5) * occ_func_0_1(22)) +
              (0.7071067812 * occ_func_0_0(8) * occ_func_0_1(24)) +
              ((0.7071067812 * occ_func_0_1(39) * occ_func_0_0(37) +
                0.7071067812 * occ_func_0_0(39) * occ_func_0_1(37))) +
              (0.7071067812 * occ_func_0_0(6) * occ_func_0_1(26)) +
              (0.7071067812 * occ_func_0_0(7) * occ_func_0_1(27)) +
              ((0.7071067812 * occ_func_0_1(35) * occ_func_0_0(34) +
                0.7071067812 * occ_func_0_0(35) * occ_func_0_1(34))) +
              (0.7071067812 * occ_func_0_0(8) * occ_func_0_1(20)) +
              (0.7071067812 * occ_func_0_0(5) * occ_func_0_1(19)) +
              ((0.7071067812 * occ_func_0_1(41) * occ_func_0_0(42) +
                0.7071067812 * occ_func_0_0(41) * occ_func_0_1(42))) +
              (0.7071067812 * occ_func_0_0(1) * occ_func_0_1(31)) +
              (0.7071067812 * occ_func_0_0(12) * occ_func_0_1(38)) +
              ((0.7071067812 * occ_func_0_1(30) * occ_func_0_0(23) +
                0.7071067812 * occ_func_0_0(30) * occ_func_0_1(23))) +
              (0.7071067812 * occ_func_0_0(9) * occ_func_0_1(25)) +
              (0.7071067812 * occ_func_0_0(4) * occ_func_0_1(23)) +
              ((0.7071067812 * occ_func_0_1(36) * occ_func_0_0(38) +
                0.7071067812 * occ_func_0_0(36) * occ_func_0_1(38))) +
              (0.7071067812 * occ_func_0_0(11) * occ_func_0_1(29)) +
              (0.7071067812 * occ_func_0_0(2) * occ_func_0_1(22)) +
              ((0.7071067812 * occ_func_0_1(32) * occ_func_0_0(39) +
                0.7071067812 * occ_func_0_0(32) * occ_func_0_1(39))) +
              (0.7071067812 * occ_func_0_0(12) * occ_func_0_1(35)) +
              (0.7071067812 * occ_func_0_0(1) * occ_func_0_1(28)) +
              ((0.7071067812 * occ_func_0_1(26) * occ_func_0_0(33) +
                0.7071067812 * occ_func_0_0(26) * occ_func_0_1(33))) +
              (0.7071067812 * occ_func_0_0(3) * occ_func_0_1(24)) +
              (0.7071067812 * occ_func_0_0(10) * occ_func_0_1(29)) +
              ((0.7071067812 * occ_func_0_1(37) * occ_func_0_0(32) +
                0.7071067812 * occ_func_0_0(37) * occ_func_0_1(32)))) /
             24.0;
}
double test_Clexulator::delta_site_eval_at_0_bfunc_3_1_5(int occ_i,
                                                         int occ_f) const {
  return (m_occ_func_0_1[occ_f] - m_occ_func_0_1[occ_i]) *
         ((occ_func_0_1(5) * occ_func_0_1(37)) +
          (occ_func_0_1(8) * occ_func_0_1(39)) +
          (occ_func_0_1(24) * occ_func_0_1(22)) +
          (occ_func_0_1(12) * occ_func_0_1(33)) +
          (occ_func_0_1(1) * occ_func_0_1(26)) +
          (occ_func_0_1(28) * occ_func_0_1(35)) +
          (occ_func_0_1(10) * occ_func_0_1(30)) +
          (occ_func_0_1(3) * occ_func_0_1(25)) +
          (occ_func_0_1(31) * occ_func_0_1(36)) +
          (occ_func_0_1(4) * occ_func_0_1(40)) +
          (occ_func_0_1(9) * occ_func_0_1(42)) +
          (occ_func_0_1(21) * occ_func_0_1(19)) +
          (occ_func_0_1(2) * occ_func_0_1(27)) +
          (occ_func_0_1(11) * occ_func_0_1(33)) +
          (occ_func_0_1(34) * occ_func_0_1(28)) +
          (occ_func_0_1(7) * occ_func_0_1(41)) +
          (occ_func_0_1(6) * occ_func_0_1(40)) +
          (occ_func_0_1(20) * occ_func_0_1(21)) +
          (occ_func_0_1(1) * occ_func_0_1(23)) +
          (occ_func_0_1(12) * occ_func_0_1(30)) +
          (occ_func_0_1(38) * occ_func_0_1(31)) +
          (occ_func_0_1(9) * occ_func_0_1(38)) +
          (occ_func_0_1(4) * occ_func_0_1(36)) +
          (occ_func_0_1(23) * occ_func_0_1(25)) +
          (occ_func_0_1(10) * occ_func_0_1(36)) +
          (occ_func_0_1(3) * occ_func_0_1(31)) +
          (occ_func_0_1(25) * occ_func_0_1(30)) +
          (occ_func_0_1(3) * occ_func_0_1(32)) +
          (occ_func_0_1(10) * occ_func_0_1(37)) +
          (occ_func_0_1(29) * occ_func_0_1(24)) +
          (occ_func_0_1(6) * occ_func_0_1(34)) +
          (occ_func_0_1(7) * occ_func_0_1(35)) +
          (occ_func_0_1(27) * occ_func_0_1(26)) +
          (occ_func_0_1(11) * occ_func_0_1(39)) +
          (occ_func_0_1(2) * occ_func_0_1(32)) +
          (occ_func_0_1(22) * occ_func_0_1(29)) +
          (occ_func_0_1(7) * occ_func_0_1(21)) +
          (occ_func_0_1(6) * occ_func_0_1(20)) +
          (occ_func_0_1(40) * occ_func_0_1(41)) +
          (occ_func_0_1(4) * occ_func_0_1(19)) +
          (occ_func_0_1(9) * occ_func_0_1(21)) +
          (occ_func_0_1(42) * occ_func_0_1(40)) +
          (occ_func_0_1(2) * occ_func_0_1(28)) +
          (occ_func_0_1(11) * occ_func_0_1(34)) +
          (occ_func_0_1(33) * occ_func_0_1(27)) +
          (occ_func_0_1(8) * occ_func_0_1(42)) +
          (occ_func_0_1(5) * occ_func_0_1(41)) +
          (occ_func_0_1(19) * occ_func_0_1(20)) +
          (occ_func_0_1(5) * occ_func_0_1(22)) +
          (occ_func_0_1(8) * occ_func_0_1(24)) +
          (occ_func_0_1(39) * occ_func_0_1(37)) +
          (occ_func_0_1(6) * occ_func_0_1(26)) +
          (occ_func_0_1(7) * occ_func_0_1(27)) +
          (occ_func_0_1(35) * occ_func_0_1(34)) +
          (occ_func_0_1(8) * occ_func_0_1(20)) +
          (occ_func_0_1(5) * occ_func_0_1(19)) +
          (occ_func_0_1(41) * occ_func_0_1(42)) +
          (occ_func_0_1(1) * occ_func_0_1(31)) +
          (occ_func_0_1(12) * occ_func_0_1(38)) +
          (occ_func_0_1(30) * occ_func_0_1(23)) +
          (occ_func_0_1(9) * occ_func_0_1(25)) +
          (occ_func_0_1(4) * occ_func_0_1(23)) +
          (occ_func_0_1(36) * occ_func_0_1(38)) +
          (occ_func_0_1(11) * occ_func_0_1(29)) +
          (occ_func_0_1(2) * occ_func_0_1(22)) +
          (occ_func_0_1(32) * occ_func_0_1(39)) +
          (occ_func_0_1(12) * occ_func_0_1(35)) +
          (occ_func_0_1(1) * occ_func_0_1(28)) +
          (occ_func_0_1(26) * occ_func_0_1(33)) +
          (occ_func_0_1(3) * occ_func_0_1(24)) +
          (occ_func_0_1(10) * occ_func_0_1(29)) +
          (occ_func_0_1(37) * occ_func_0_1(32))) /
         24.0;
}

/**** Basis functions for orbit 3, 2****
#Points: 3
MaxLength: 5.6568542  MinLength: 2.8284271
 0.0000000   0.0000000   0.0000000 A B C
 1.0000000   0.0000000   0.0000000 A B C
 2.0000000   0.0000000   0.0000000 A B C
****/
double test_Clexulator::eval_bfunc_3_2_0() const {
  return ((occ_func_0_0(0) * occ_func_0_0(12) * occ_func_0_0(54)) +
          (occ_func_0_0(0) * occ_func_0_0(8) * occ_func_0_0(50)) +
          (occ_func_0_0(0) * occ_func_0_0(7) * occ_func_0_0(49)) +
          (occ_func_0_0(0) * occ_func_0_0(11) * occ_func_0_0(53)) +
          (occ_func_0_0(0) * occ_func_0_0(9) * occ_func_0_0(51)) +
          (occ_func_0_0(0) * occ_func_0_0(10) * occ_func_0_0(52))) /
         6.0;
}
double test_Clexulator::eval_bfunc_3_2_1() const {
  return (((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(12) *
                occ_func_0_0(54) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(12) *
                occ_func_0_1(54))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(8) *
                occ_func_0_0(50) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(8) *
                occ_func_0_1(50))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(7) *
                occ_func_0_0(49) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(7) *
                occ_func_0_1(49))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(11) *
                occ_func_0_0(53) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(11) *
                occ_func_0_1(53))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(9) *
                occ_func_0_0(51) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(9) *
                occ_func_0_1(51))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(10) *
                occ_func_0_0(52) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(10) *
                occ_func_0_1(52)))) /
         6.0;
}
double test_Clexulator::eval_bfunc_3_2_2() const {
  return ((occ_func_0_0(0) * occ_func_0_1(12) * occ_func_0_0(54)) +
          (occ_func_0_0(0) * occ_func_0_1(8) * occ_func_0_0(50)) +
          (occ_func_0_0(0) * occ_func_0_1(7) * occ_func_0_0(49)) +
          (occ_func_0_0(0) * occ_func_0_1(11) * occ_func_0_0(53)) +
          (occ_func_0_0(0) * occ_func_0_1(9) * occ_func_0_0(51)) +
          (occ_func_0_0(0) * occ_func_0_1(10) * occ_func_0_0(52))) /
         6.0;
}
double test_Clexulator::eval_bfunc_3_2_3() const {
  return (((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(12) *
                occ_func_0_0(54) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(12) *
                occ_func_0_1(54))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(8) *
                occ_func_0_0(50) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(8) *
                occ_func_0_1(50))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(7) *
                occ_func_0_0(49) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(7) *
                occ_func_0_1(49))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(11) *
                occ_func_0_0(53) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(11) *
                occ_func_0_1(53))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(9) *
                occ_func_0_0(51) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(9) *
                occ_func_0_1(51))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(10) *
                occ_func_0_0(52) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(10) *
                occ_func_0_1(52)))) /
         6.0;
}
double test_Clexulator::eval_bfunc_3_2_4() const {
  return ((occ_func_0_1(0) * occ_func_0_0(12) * occ_func_0_1(54)) +
          (occ_func_0_1(0) * occ_func_0_0(8) * occ_func_0_1(50)) +
          (occ_func_0_1(0) * occ_func_0_0(7) * occ_func_0_1(49)) +
          (occ_func_0_1(0) * occ_func_0_0(11) * occ_func_0_1(53)) +
          (occ_func_0_1(0) * occ_func_0_0(9) * occ_func_0_1(51)) +
          (occ_func_0_1(0) * occ_func_0_0(10) * occ_func_0_1(52))) /
         6.0;
}
double test_Clexulator::eval_bfunc_3_2_5() const {
  return ((occ_func_0_1(0) * occ_func_0_1(12) * occ_func_0_1(54)) +
          (occ_func_0_1(0) * occ_func_0_1(8) * occ_func_0_1(50)) +
          (occ_func_0_1(0) * occ_func_0_1(7) * occ_func_0_1(49)) +
          (occ_func_0_1(0) * occ_func_0_1(11) * occ_func_0_1(53)) +
          (occ_func_0_1(0) * occ_func_0_1(9) * occ_func_0_1(51)) +
          (occ_func_0_1(0) * occ_func_0_1(10) * occ_func_0_1(52))) /
         6.0;
}

double test_Clexulator::site_eval_at_0_bfunc_3_2_0() const {
  return ((occ_func_0_0(0) * occ_func_0_0(12) * occ_func_0_0(54)) +
          (occ_func_0_0(1) * occ_func_0_0(0) * occ_func_0_0(12)) +
          (occ_func_0_0(43) * occ_func_0_0(1) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(8) * occ_func_0_0(50)) +
          (occ_func_0_0(5) * occ_func_0_0(0) * occ_func_0_0(8)) +
          (occ_func_0_0(47) * occ_func_0_0(5) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(7) * occ_func_0_0(49)) +
          (occ_func_0_0(6) * occ_func_0_0(0) * occ_func_0_0(7)) +
          (occ_func_0_0(48) * occ_func_0_0(6) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(11) * occ_func_0_0(53)) +
          (occ_func_0_0(2) * occ_func_0_0(0) * occ_func_0_0(11)) +
          (occ_func_0_0(44) * occ_func_0_0(2) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(9) * occ_func_0_0(51)) +
          (occ_func_0_0(4) * occ_func_0_0(0) * occ_func_0_0(9)) +
          (occ_func_0_0(46) * occ_func_0_0(4) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(10) * occ_func_0_0(52)) +
          (occ_func_0_0(3) * occ_func_0_0(0) * occ_func_0_0(10)) +
          (occ_func_0_0(45) * occ_func_0_0(3) * occ_func_0_0(0))) /
         6.0;
}
double test_Clexulator::site_eval_at_0_bfunc_3_2_1() const {
  return (((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(12) *
                occ_func_0_0(54) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(12) *
                occ_func_0_1(54))) +
          ((0.7071067812 * occ_func_0_1(1) * occ_func_0_0(0) *
                occ_func_0_0(12) +
            0.7071067812 * occ_func_0_0(1) * occ_func_0_0(0) *
                occ_func_0_1(12))) +
          ((0.7071067812 * occ_func_0_1(43) * occ_func_0_0(1) *
                occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(43) * occ_func_0_0(1) *
                occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(8) *
                occ_func_0_0(50) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(8) *
                occ_func_0_1(50))) +
          ((0.7071067812 * occ_func_0_1(5) * occ_func_0_0(0) * occ_func_0_0(8) +
            0.7071067812 * occ_func_0_0(5) * occ_func_0_0(0) *
                occ_func_0_1(8))) +
          ((0.7071067812 * occ_func_0_1(47) * occ_func_0_0(5) *
                occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(47) * occ_func_0_0(5) *
                occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(7) *
                occ_func_0_0(49) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(7) *
                occ_func_0_1(49))) +
          ((0.7071067812 * occ_func_0_1(6) * occ_func_0_0(0) * occ_func_0_0(7) +
            0.7071067812 * occ_func_0_0(6) * occ_func_0_0(0) *
                occ_func_0_1(7))) +
          ((0.7071067812 * occ_func_0_1(48) * occ_func_0_0(6) *
                occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(48) * occ_func_0_0(6) *
                occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(11) *
                occ_func_0_0(53) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(11) *
                occ_func_0_1(53))) +
          ((0.7071067812 * occ_func_0_1(2) * occ_func_0_0(0) *
                occ_func_0_0(11) +
            0.7071067812 * occ_func_0_0(2) * occ_func_0_0(0) *
                occ_func_0_1(11))) +
          ((0.7071067812 * occ_func_0_1(44) * occ_func_0_0(2) *
                occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(44) * occ_func_0_0(2) *
                occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(9) *
                occ_func_0_0(51) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(9) *
                occ_func_0_1(51))) +
          ((0.7071067812 * occ_func_0_1(4) * occ_func_0_0(0) * occ_func_0_0(9) +
            0.7071067812 * occ_func_0_0(4) * occ_func_0_0(0) *
                occ_func_0_1(9))) +
          ((0.7071067812 * occ_func_0_1(46) * occ_func_0_0(4) *
                occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(46) * occ_func_0_0(4) *
                occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(10) *
                occ_func_0_0(52) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(10) *
                occ_func_0_1(52))) +
          ((0.7071067812 * occ_func_0_1(3) * occ_func_0_0(0) *
                occ_func_0_0(10) +
            0.7071067812 * occ_func_0_0(3) * occ_func_0_0(0) *
                occ_func_0_1(10))) +
          ((0.7071067812 * occ_func_0_1(45) * occ_func_0_0(3) *
                occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(45) * occ_func_0_0(3) *
                occ_func_0_1(0)))) /
         6.0;
}
double test_Clexulator::site_eval_at_0_bfunc_3_2_2() const {
  return ((occ_func_0_0(0) * occ_func_0_1(12) * occ_func_0_0(54)) +
          (occ_func_0_0(1) * occ_func_0_1(0) * occ_func_0_0(12)) +
          (occ_func_0_0(43) * occ_func_0_1(1) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_1(8) * occ_func_0_0(50)) +
          (occ_func_0_0(5) * occ_func_0_1(0) * occ_func_0_0(8)) +
          (occ_func_0_0(47) * occ_func_0_1(5) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_1(7) * occ_func_0_0(49)) +
          (occ_func_0_0(6) * occ_func_0_1(0) * occ_func_0_0(7)) +
          (occ_func_0_0(48) * occ_func_0_1(6) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_1(11) * occ_func_0_0(53)) +
          (occ_func_0_0(2) * occ_func_0_1(0) * occ_func_0_0(11)) +
          (occ_func_0_0(44) * occ_func_0_1(2) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_1(9) * occ_func_0_0(51)) +
          (occ_func_0_0(4) * occ_func_0_1(0) * occ_func_0_0(9)) +
          (occ_func_0_0(46) * occ_func_0_1(4) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_1(10) * occ_func_0_0(52)) +
          (occ_func_0_0(3) * occ_func_0_1(0) * occ_func_0_0(10)) +
          (occ_func_0_0(45) * occ_func_0_1(3) * occ_func_0_0(0))) /
         6.0;
}
double test_Clexulator::site_eval_at_0_bfunc_3_2_3() const {
  return (((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(12) *
                occ_func_0_0(54) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(12) *
                occ_func_0_1(54))) +
          ((0.7071067812 * occ_func_0_1(1) * occ_func_0_1(0) *
                occ_func_0_0(12) +
            0.7071067812 * occ_func_0_0(1) * occ_func_0_1(0) *
                occ_func_0_1(12))) +
          ((0.7071067812 * occ_func_0_1(43) * occ_func_0_1(1) *
                occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(43) * occ_func_0_1(1) *
                occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(8) *
                occ_func_0_0(50) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(8) *
                occ_func_0_1(50))) +
          ((0.7071067812 * occ_func_0_1(5) * occ_func_0_1(0) * occ_func_0_0(8) +
            0.7071067812 * occ_func_0_0(5) * occ_func_0_1(0) *
                occ_func_0_1(8))) +
          ((0.7071067812 * occ_func_0_1(47) * occ_func_0_1(5) *
                occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(47) * occ_func_0_1(5) *
                occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(7) *
                occ_func_0_0(49) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(7) *
                occ_func_0_1(49))) +
          ((0.7071067812 * occ_func_0_1(6) * occ_func_0_1(0) * occ_func_0_0(7) +
            0.7071067812 * occ_func_0_0(6) * occ_func_0_1(0) *
                occ_func_0_1(7))) +
          ((0.7071067812 * occ_func_0_1(48) * occ_func_0_1(6) *
                occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(48) * occ_func_0_1(6) *
                occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(11) *
                occ_func_0_0(53) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(11) *
                occ_func_0_1(53))) +
          ((0.7071067812 * occ_func_0_1(2) * occ_func_0_1(0) *
                occ_func_0_0(11) +
            0.7071067812 * occ_func_0_0(2) * occ_func_0_1(0) *
                occ_func_0_1(11))) +
          ((0.7071067812 * occ_func_0_1(44) * occ_func_0_1(2) *
                occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(44) * occ_func_0_1(2) *
                occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(9) *
                occ_func_0_0(51) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(9) *
                occ_func_0_1(51))) +
          ((0.7071067812 * occ_func_0_1(4) * occ_func_0_1(0) * occ_func_0_0(9) +
            0.7071067812 * occ_func_0_0(4) * occ_func_0_1(0) *
                occ_func_0_1(9))) +
          ((0.7071067812 * occ_func_0_1(46) * occ_func_0_1(4) *
                occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(46) * occ_func_0_1(4) *
                occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(10) *
                occ_func_0_0(52) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(10) *
                occ_func_0_1(52))) +
          ((0.7071067812 * occ_func_0_1(3) * occ_func_0_1(0) *
                occ_func_0_0(10) +
            0.7071067812 * occ_func_0_0(3) * occ_func_0_1(0) *
                occ_func_0_1(10))) +
          ((0.7071067812 * occ_func_0_1(45) * occ_func_0_1(3) *
                occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(45) * occ_func_0_1(3) *
                occ_func_0_1(0)))) /
         6.0;
}
double test_Clexulator::site_eval_at_0_bfunc_3_2_4() const {
  return ((occ_func_0_1(0) * occ_func_0_0(12) * occ_func_0_1(54)) +
          (occ_func_0_1(1) * occ_func_0_0(0) * occ_func_0_1(12)) +
          (occ_func_0_1(43) * occ_func_0_0(1) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_0(8) * occ_func_0_1(50)) +
          (occ_func_0_1(5) * occ_func_0_0(0) * occ_func_0_1(8)) +
          (occ_func_0_1(47) * occ_func_0_0(5) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_0(7) * occ_func_0_1(49)) +
          (occ_func_0_1(6) * occ_func_0_0(0) * occ_func_0_1(7)) +
          (occ_func_0_1(48) * occ_func_0_0(6) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_0(11) * occ_func_0_1(53)) +
          (occ_func_0_1(2) * occ_func_0_0(0) * occ_func_0_1(11)) +
          (occ_func_0_1(44) * occ_func_0_0(2) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_0(9) * occ_func_0_1(51)) +
          (occ_func_0_1(4) * occ_func_0_0(0) * occ_func_0_1(9)) +
          (occ_func_0_1(46) * occ_func_0_0(4) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_0(10) * occ_func_0_1(52)) +
          (occ_func_0_1(3) * occ_func_0_0(0) * occ_func_0_1(10)) +
          (occ_func_0_1(45) * occ_func_0_0(3) * occ_func_0_1(0))) /
         6.0;
}
double test_Clexulator::site_eval_at_0_bfunc_3_2_5() const {
  return ((occ_func_0_1(0) * occ_func_0_1(12) * occ_func_0_1(54)) +
          (occ_func_0_1(1) * occ_func_0_1(0) * occ_func_0_1(12)) +
          (occ_func_0_1(43) * occ_func_0_1(1) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(8) * occ_func_0_1(50)) +
          (occ_func_0_1(5) * occ_func_0_1(0) * occ_func_0_1(8)) +
          (occ_func_0_1(47) * occ_func_0_1(5) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(7) * occ_func_0_1(49)) +
          (occ_func_0_1(6) * occ_func_0_1(0) * occ_func_0_1(7)) +
          (occ_func_0_1(48) * occ_func_0_1(6) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(11) * occ_func_0_1(53)) +
          (occ_func_0_1(2) * occ_func_0_1(0) * occ_func_0_1(11)) +
          (occ_func_0_1(44) * occ_func_0_1(2) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(9) * occ_func_0_1(51)) +
          (occ_func_0_1(4) * occ_func_0_1(0) * occ_func_0_1(9)) +
          (occ_func_0_1(46) * occ_func_0_1(4) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(10) * occ_func_0_1(52)) +
          (occ_func_0_1(3) * occ_func_0_1(0) * occ_func_0_1(10)) +
          (occ_func_0_1(45) * occ_func_0_1(3) * occ_func_0_1(0))) /
         6.0;
}

double test_Clexulator::delta_site_eval_at_0_bfunc_3_2_0(int occ_i,
                                                         int occ_f) const {
  return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) *
         ((occ_func_0_0(12) * occ_func_0_0(54)) +
          (occ_func_0_0(1) * occ_func_0_0(12)) +
          (occ_func_0_0(43) * occ_func_0_0(1)) +
          (occ_func_0_0(8) * occ_func_0_0(50)) +
          (occ_func_0_0(5) * occ_func_0_0(8)) +
          (occ_func_0_0(47) * occ_func_0_0(5)) +
          (occ_func_0_0(7) * occ_func_0_0(49)) +
          (occ_func_0_0(6) * occ_func_0_0(7)) +
          (occ_func_0_0(48) * occ_func_0_0(6)) +
          (occ_func_0_0(11) * occ_func_0_0(53)) +
          (occ_func_0_0(2) * occ_func_0_0(11)) +
          (occ_func_0_0(44) * occ_func_0_0(2)) +
          (occ_func_0_0(9) * occ_func_0_0(51)) +
          (occ_func_0_0(4) * occ_func_0_0(9)) +
          (occ_func_0_0(46) * occ_func_0_0(4)) +
          (occ_func_0_0(10) * occ_func_0_0(52)) +
          (occ_func_0_0(3) * occ_func_0_0(10)) +
          (occ_func_0_0(45) * occ_func_0_0(3))) /
         6.0;
}
double test_Clexulator::delta_site_eval_at_0_bfunc_3_2_1(int occ_i,
                                                         int occ_f) const {
  return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) *
             ((0.7071067812 * occ_func_0_0(12) * occ_func_0_1(54)) +
              ((0.7071067812 * occ_func_0_1(1) * occ_func_0_0(12) +
                0.7071067812 * occ_func_0_0(1) * occ_func_0_1(12))) +
              (0.7071067812 * occ_func_0_1(43) * occ_func_0_0(1)) +
              (0.7071067812 * occ_func_0_0(8) * occ_func_0_1(50)) +
              ((0.7071067812 * occ_func_0_1(5) * occ_func_0_0(8) +
                0.7071067812 * occ_func_0_0(5) * occ_func_0_1(8))) +
              (0.7071067812 * occ_func_0_1(47) * occ_func_0_0(5)) +
              (0.7071067812 * occ_func_0_0(7) * occ_func_0_1(49)) +
              ((0.7071067812 * occ_func_0_1(6) * occ_func_0_0(7) +
                0.7071067812 * occ_func_0_0(6) * occ_func_0_1(7))) +
              (0.7071067812 * occ_func_0_1(48) * occ_func_0_0(6)) +
              (0.7071067812 * occ_func_0_0(11) * occ_func_0_1(53)) +
              ((0.7071067812 * occ_func_0_1(2) * occ_func_0_0(11) +
                0.7071067812 * occ_func_0_0(2) * occ_func_0_1(11))) +
              (0.7071067812 * occ_func_0_1(44) * occ_func_0_0(2)) +
              (0.7071067812 * occ_func_0_0(9) * occ_func_0_1(51)) +
              ((0.7071067812 * occ_func_0_1(4) * occ_func_0_0(9) +
                0.7071067812 * occ_func_0_0(4) * occ_func_0_1(9))) +
              (0.7071067812 * occ_func_0_1(46) * occ_func_0_0(4)) +
              (0.7071067812 * occ_func_0_0(10) * occ_func_0_1(52)) +
              ((0.7071067812 * occ_func_0_1(3) * occ_func_0_0(10) +
                0.7071067812 * occ_func_0_0(3) * occ_func_0_1(10))) +
              (0.7071067812 * occ_func_0_1(45) * occ_func_0_0(3))) /
             6.0 +
         (m_occ_func_0_1[occ_f] - m_occ_func_0_1[occ_i]) *
             ((0.7071067812 * occ_func_0_0(12) * occ_func_0_0(54)) +
              (0.7071067812 * occ_func_0_0(43) * occ_func_0_0(1)) +
              (0.7071067812 * occ_func_0_0(8) * occ_func_0_0(50)) +
              (0.7071067812 * occ_func_0_0(47) * occ_func_0_0(5)) +
              (0.7071067812 * occ_func_0_0(7) * occ_func_0_0(49)) +
              (0.7071067812 * occ_func_0_0(48) * occ_func_0_0(6)) +
              (0.7071067812 * occ_func_0_0(11) * occ_func_0_0(53)) +
              (0.7071067812 * occ_func_0_0(44) * occ_func_0_0(2)) +
              (0.7071067812 * occ_func_0_0(9) * occ_func_0_0(51)) +
              (0.7071067812 * occ_func_0_0(46) * occ_func_0_0(4)) +
              (0.7071067812 * occ_func_0_0(10) * occ_func_0_0(52)) +
              (0.7071067812 * occ_func_0_0(45) * occ_func_0_0(3))) /
             6.0;
}
double test_Clexulator::delta_site_eval_at_0_bfunc_3_2_2(int occ_i,
                                                         int occ_f) const {
  return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) *
             ((occ_func_0_1(12) * occ_func_0_0(54)) +
              (occ_func_0_0(43) * occ_func_0_1(1)) +
              (occ_func_0_1(8) * occ_func_0_0(50)) +
              (occ_func_0_0(47) * occ_func_0_1(5)) +
              (occ_func_0_1(7) * occ_func_0_0(49)) +
              (occ_func_0_0(48) * occ_func_0_1(6)) +
              (occ_func_0_1(11) * occ_func_0_0(53)) +
              (occ_func_0_0(44) * occ_func_0_1(2)) +
              (occ_func_0_1(9) * occ_func_0_0(51)) +
              (occ_func_0_0(46) * occ_func_0_1(4)) +
              (occ_func_0_1(10) * occ_func_0_0(52)) +
              (occ_func_0_0(45) * occ_func_0_1(3))) /
             6.0 +
         (m_occ_func_0_1[occ_f] - m_occ_func_0_1[occ_i]) *
             ((occ_func_0_0(1) * occ_func_0_0(12)) +
              (occ_func_0_0(5) * occ_func_0_0(8)) +
              (occ_func_0_0(6) * occ_func_0_0(7)) +
              (occ_func_0_0(2) * occ_func_0_0(11)) +
              (occ_func_0_0(4) * occ_func_0_0(9)) +
              (occ_func_0_0(3) * occ_func_0_0(10))) /
             6.0;
}
double test_Clexulator::delta_site_eval_at_0_bfunc_3_2_3(int occ_i,
                                                         int occ_f) const {
  return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) *
             ((0.7071067812 * occ_func_0_1(12) * occ_func_0_1(54)) +
              (0.7071067812 * occ_func_0_1(43) * occ_func_0_1(1)) +
              (0.7071067812 * occ_func_0_1(8) * occ_func_0_1(50)) +
              (0.7071067812 * occ_func_0_1(47) * occ_func_0_1(5)) +
              (0.7071067812 * occ_func_0_1(7) * occ_func_0_1(49)) +
              (0.7071067812 * occ_func_0_1(48) * occ_func_0_1(6)) +
              (0.7071067812 * occ_func_0_1(11) * occ_func_0_1(53)) +
              (0.7071067812 * occ_func_0_1(44) * occ_func_0_1(2)) +
              (0.7071067812 * occ_func_0_1(9) * occ_func_0_1(51)) +
              (0.7071067812 * occ_func_0_1(46) * occ_func_0_1(4)) +
              (0.7071067812 * occ_func_0_1(10) * occ_func_0_1(52)) +
              (0.7071067812 * occ_func_0_1(45) * occ_func_0_1(3))) /
             6.0 +
         (m_occ_func_0_1[occ_f] - m_occ_func_0_1[occ_i]) *
             ((0.7071067812 * occ_func_0_1(12) * occ_func_0_0(54)) +
              ((0.7071067812 * occ_func_0_1(1) * occ_func_0_0(12) +
                0.7071067812 * occ_func_0_0(1) * occ_func_0_1(12))) +
              (0.7071067812 * occ_func_0_0(43) * occ_func_0_1(1)) +
              (0.7071067812 * occ_func_0_1(8) * occ_func_0_0(50)) +
              ((0.7071067812 * occ_func_0_1(5) * occ_func_0_0(8) +
                0.7071067812 * occ_func_0_0(5) * occ_func_0_1(8))) +
              (0.7071067812 * occ_func_0_0(47) * occ_func_0_1(5)) +
              (0.7071067812 * occ_func_0_1(7) * occ_func_0_0(49)) +
              ((0.7071067812 * occ_func_0_1(6) * occ_func_0_0(7) +
                0.7071067812 * occ_func_0_0(6) * occ_func_0_1(7))) +
              (0.7071067812 * occ_func_0_0(48) * occ_func_0_1(6)) +
              (0.7071067812 * occ_func_0_1(11) * occ_func_0_0(53)) +
              ((0.7071067812 * occ_func_0_1(2) * occ_func_0_0(11) +
                0.7071067812 * occ_func_0_0(2) * occ_func_0_1(11))) +
              (0.7071067812 * occ_func_0_0(44) * occ_func_0_1(2)) +
              (0.7071067812 * occ_func_0_1(9) * occ_func_0_0(51)) +
              ((0.7071067812 * occ_func_0_1(4) * occ_func_0_0(9) +
                0.7071067812 * occ_func_0_0(4) * occ_func_0_1(9))) +
              (0.7071067812 * occ_func_0_0(46) * occ_func_0_1(4)) +
              (0.7071067812 * occ_func_0_1(10) * occ_func_0_0(52)) +
              ((0.7071067812 * occ_func_0_1(3) * occ_func_0_0(10) +
                0.7071067812 * occ_func_0_0(3) * occ_func_0_1(10))) +
              (0.7071067812 * occ_func_0_0(45) * occ_func_0_1(3))) /
             6.0;
}
double test_Clexulator::delta_site_eval_at_0_bfunc_3_2_4(int occ_i,
                                                         int occ_f) const {
  return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) *
             ((occ_func_0_1(1) * occ_func_0_1(12)) +
              (occ_func_0_1(5) * occ_func_0_1(8)) +
              (occ_func_0_1(6) * occ_func_0_1(7)) +
              (occ_func_0_1(2) * occ_func_0_1(11)) +
              (occ_func_0_1(4) * occ_func_0_1(9)) +
              (occ_func_0_1(3) * occ_func_0_1(10))) /
             6.0 +
         (m_occ_func_0_1[occ_f] - m_occ_func_0_1[occ_i]) *
             ((occ_func_0_0(12) * occ_func_0_1(54)) +
              (occ_func_0_1(43) * occ_func_0_0(1)) +
              (occ_func_0_0(8) * occ_func_0_1(50)) +
              (occ_func_0_1(47) * occ_func_0_0(5)) +
              (occ_func_0_0(7) * occ_func_0_1(49)) +
              (occ_func_0_1(48) * occ_func_0_0(6)) +
              (occ_func_0_0(11) * occ_func_0_1(53)) +
              (occ_func_0_1(44) * occ_func_0_0(2)) +
              (occ_func_0_0(9) * occ_func_0_1(51)) +
              (occ_func_0_1(46) * occ_func_0_0(4)) +
              (occ_func_0_0(10) * occ_func_0_1(52)) +
              (occ_func_0_1(45) * occ_func_0_0(3))) /
             6.0;
}
double test_Clexulator::delta_site_eval_at_0_bfunc_3_2_5(int occ_i,
                                                         int occ_f) const {
  return (m_occ_func_0_1[occ_f] - m_occ_func_0_1[occ_i]) *
         ((occ_func_0_1(12) * occ_func_0_1(54)) +
          (occ_func_0_1(1) * occ_func_0_1(12)) +
          (occ_func_0_1(43) * occ_func_0_1(1)) +
          (occ_func_0_1(8) * occ_func_0_1(50)) +
          (occ_func_0_1(5) * occ_func_0_1(8)) +
          (occ_func_0_1(47) * occ_func_0_1(5)) +
          (occ_func_0_1(7) * occ_func_0_1(49)) +
          (occ_func_0_1(6) * occ_func_0_1(7)) +
          (occ_func_0_1(48) * occ_func_0_1(6)) +
          (occ_func_0_1(11) * occ_func_0_1(53)) +
          (occ_func_0_1(2) * occ_func_0_1(11)) +
          (occ_func_0_1(44) * occ_func_0_1(2)) +
          (occ_func_0_1(9) * occ_func_0_1(51)) +
          (occ_func_0_1(4) * occ_func_0_1(9)) +
          (occ_func_0_1(46) * occ_func_0_1(4)) +
          (occ_func_0_1(10) * occ_func_0_1(52)) +
          (occ_func_0_1(3) * occ_func_0_1(10)) +
          (occ_func_0_1(45) * occ_func_0_1(3))) /
         6.0;
}

/**** Basis functions for orbit 3, 3****
#Points: 3
MaxLength: 6.9282032  MinLength: 2.8284271
 0.0000000   0.0000000   0.0000000 A B C
 0.0000000   0.0000000   1.0000000 A B C
 1.0000000   1.0000000   1.0000000 A B C
****/
double test_Clexulator::eval_bfunc_3_3_0() const {
  return ((occ_func_0_0(0) * occ_func_0_0(7) * occ_func_0_0(85)) +
          (occ_func_0_0(0) * occ_func_0_0(9) * occ_func_0_0(82)) +
          (occ_func_0_0(0) * occ_func_0_0(5) * occ_func_0_0(81)) +
          (occ_func_0_0(0) * occ_func_0_0(10) * occ_func_0_0(86)) +
          (occ_func_0_0(0) * occ_func_0_0(3) * occ_func_0_0(82)) +
          (occ_func_0_0(0) * occ_func_0_0(12) * occ_func_0_0(86)) +
          (occ_func_0_0(0) * occ_func_0_0(2) * occ_func_0_0(81)) +
          (occ_func_0_0(0) * occ_func_0_0(8) * occ_func_0_0(84)) +
          (occ_func_0_0(0) * occ_func_0_0(11) * occ_func_0_0(84)) +
          (occ_func_0_0(0) * occ_func_0_0(9) * occ_func_0_0(85)) +
          (occ_func_0_0(0) * occ_func_0_0(4) * occ_func_0_0(83)) +
          (occ_func_0_0(0) * occ_func_0_0(12) * occ_func_0_0(85)) +
          (occ_func_0_0(0) * occ_func_0_0(2) * occ_func_0_0(79)) +
          (occ_func_0_0(0) * occ_func_0_0(1) * occ_func_0_0(79)) +
          (occ_func_0_0(0) * occ_func_0_0(5) * occ_func_0_0(83)) +
          (occ_func_0_0(0) * occ_func_0_0(11) * occ_func_0_0(86)) +
          (occ_func_0_0(0) * occ_func_0_0(4) * occ_func_0_0(80)) +
          (occ_func_0_0(0) * occ_func_0_0(8) * occ_func_0_0(82)) +
          (occ_func_0_0(0) * occ_func_0_0(3) * occ_func_0_0(79)) +
          (occ_func_0_0(0) * occ_func_0_0(6) * occ_func_0_0(84)) +
          (occ_func_0_0(0) * occ_func_0_0(7) * occ_func_0_0(81)) +
          (occ_func_0_0(0) * occ_func_0_0(6) * occ_func_0_0(80)) +
          (occ_func_0_0(0) * occ_func_0_0(10) * occ_func_0_0(83)) +
          (occ_func_0_0(0) * occ_func_0_0(1) * occ_func_0_0(80))) /
         24.0;
}
double test_Clexulator::eval_bfunc_3_3_1() const {
  return ((occ_func_0_1(0) * occ_func_0_0(7) * occ_func_0_0(85)) +
          (occ_func_0_1(0) * occ_func_0_0(9) * occ_func_0_0(82)) +
          (occ_func_0_1(0) * occ_func_0_0(5) * occ_func_0_0(81)) +
          (occ_func_0_1(0) * occ_func_0_0(10) * occ_func_0_0(86)) +
          (occ_func_0_1(0) * occ_func_0_0(3) * occ_func_0_0(82)) +
          (occ_func_0_1(0) * occ_func_0_0(12) * occ_func_0_0(86)) +
          (occ_func_0_1(0) * occ_func_0_0(2) * occ_func_0_0(81)) +
          (occ_func_0_1(0) * occ_func_0_0(8) * occ_func_0_0(84)) +
          (occ_func_0_1(0) * occ_func_0_0(11) * occ_func_0_0(84)) +
          (occ_func_0_1(0) * occ_func_0_0(9) * occ_func_0_0(85)) +
          (occ_func_0_1(0) * occ_func_0_0(4) * occ_func_0_0(83)) +
          (occ_func_0_1(0) * occ_func_0_0(12) * occ_func_0_0(85)) +
          (occ_func_0_1(0) * occ_func_0_0(2) * occ_func_0_0(79)) +
          (occ_func_0_1(0) * occ_func_0_0(1) * occ_func_0_0(79)) +
          (occ_func_0_1(0) * occ_func_0_0(5) * occ_func_0_0(83)) +
          (occ_func_0_1(0) * occ_func_0_0(11) * occ_func_0_0(86)) +
          (occ_func_0_1(0) * occ_func_0_0(4) * occ_func_0_0(80)) +
          (occ_func_0_1(0) * occ_func_0_0(8) * occ_func_0_0(82)) +
          (occ_func_0_1(0) * occ_func_0_0(3) * occ_func_0_0(79)) +
          (occ_func_0_1(0) * occ_func_0_0(6) * occ_func_0_0(84)) +
          (occ_func_0_1(0) * occ_func_0_0(7) * occ_func_0_0(81)) +
          (occ_func_0_1(0) * occ_func_0_0(6) * occ_func_0_0(80)) +
          (occ_func_0_1(0) * occ_func_0_0(10) * occ_func_0_0(83)) +
          (occ_func_0_1(0) * occ_func_0_0(1) * occ_func_0_0(80))) /
         24.0;
}
double test_Clexulator::eval_bfunc_3_3_2() const {
  return ((occ_func_0_0(0) * occ_func_0_1(7) * occ_func_0_0(85)) +
          (occ_func_0_0(0) * occ_func_0_1(9) * occ_func_0_0(82)) +
          (occ_func_0_0(0) * occ_func_0_1(5) * occ_func_0_0(81)) +
          (occ_func_0_0(0) * occ_func_0_1(10) * occ_func_0_0(86)) +
          (occ_func_0_0(0) * occ_func_0_1(3) * occ_func_0_0(82)) +
          (occ_func_0_0(0) * occ_func_0_1(12) * occ_func_0_0(86)) +
          (occ_func_0_0(0) * occ_func_0_1(2) * occ_func_0_0(81)) +
          (occ_func_0_0(0) * occ_func_0_1(8) * occ_func_0_0(84)) +
          (occ_func_0_0(0) * occ_func_0_1(11) * occ_func_0_0(84)) +
          (occ_func_0_0(0) * occ_func_0_1(9) * occ_func_0_0(85)) +
          (occ_func_0_0(0) * occ_func_0_1(4) * occ_func_0_0(83)) +
          (occ_func_0_0(0) * occ_func_0_1(12) * occ_func_0_0(85)) +
          (occ_func_0_0(0) * occ_func_0_1(2) * occ_func_0_0(79)) +
          (occ_func_0_0(0) * occ_func_0_1(1) * occ_func_0_0(79)) +
          (occ_func_0_0(0) * occ_func_0_1(5) * occ_func_0_0(83)) +
          (occ_func_0_0(0) * occ_func_0_1(11) * occ_func_0_0(86)) +
          (occ_func_0_0(0) * occ_func_0_1(4) * occ_func_0_0(80)) +
          (occ_func_0_0(0) * occ_func_0_1(8) * occ_func_0_0(82)) +
          (occ_func_0_0(0) * occ_func_0_1(3) * occ_func_0_0(79)) +
          (occ_func_0_0(0) * occ_func_0_1(6) * occ_func_0_0(84)) +
          (occ_func_0_0(0) * occ_func_0_1(7) * occ_func_0_0(81)) +
          (occ_func_0_0(0) * occ_func_0_1(6) * occ_func_0_0(80)) +
          (occ_func_0_0(0) * occ_func_0_1(10) * occ_func_0_0(83)) +
          (occ_func_0_0(0) * occ_func_0_1(1) * occ_func_0_0(80))) /
         24.0;
}
double test_Clexulator::eval_bfunc_3_3_3() const {
  return ((occ_func_0_1(0) * occ_func_0_1(7) * occ_func_0_0(85)) +
          (occ_func_0_1(0) * occ_func_0_1(9) * occ_func_0_0(82)) +
          (occ_func_0_1(0) * occ_func_0_1(5) * occ_func_0_0(81)) +
          (occ_func_0_1(0) * occ_func_0_1(10) * occ_func_0_0(86)) +
          (occ_func_0_1(0) * occ_func_0_1(3) * occ_func_0_0(82)) +
          (occ_func_0_1(0) * occ_func_0_1(12) * occ_func_0_0(86)) +
          (occ_func_0_1(0) * occ_func_0_1(2) * occ_func_0_0(81)) +
          (occ_func_0_1(0) * occ_func_0_1(8) * occ_func_0_0(84)) +
          (occ_func_0_1(0) * occ_func_0_1(11) * occ_func_0_0(84)) +
          (occ_func_0_1(0) * occ_func_0_1(9) * occ_func_0_0(85)) +
          (occ_func_0_1(0) * occ_func_0_1(4) * occ_func_0_0(83)) +
          (occ_func_0_1(0) * occ_func_0_1(12) * occ_func_0_0(85)) +
          (occ_func_0_1(0) * occ_func_0_1(2) * occ_func_0_0(79)) +
          (occ_func_0_1(0) * occ_func_0_1(1) * occ_func_0_0(79)) +
          (occ_func_0_1(0) * occ_func_0_1(5) * occ_func_0_0(83)) +
          (occ_func_0_1(0) * occ_func_0_1(11) * occ_func_0_0(86)) +
          (occ_func_0_1(0) * occ_func_0_1(4) * occ_func_0_0(80)) +
          (occ_func_0_1(0) * occ_func_0_1(8) * occ_func_0_0(82)) +
          (occ_func_0_1(0) * occ_func_0_1(3) * occ_func_0_0(79)) +
          (occ_func_0_1(0) * occ_func_0_1(6) * occ_func_0_0(84)) +
          (occ_func_0_1(0) * occ_func_0_1(7) * occ_func_0_0(81)) +
          (occ_func_0_1(0) * occ_func_0_1(6) * occ_func_0_0(80)) +
          (occ_func_0_1(0) * occ_func_0_1(10) * occ_func_0_0(83)) +
          (occ_func_0_1(0) * occ_func_0_1(1) * occ_func_0_0(80))) /
         24.0;
}
double test_Clexulator::eval_bfunc_3_3_4() const {
  return ((occ_func_0_0(0) * occ_func_0_0(7) * occ_func_0_1(85)) +
          (occ_func_0_0(0) * occ_func_0_0(9) * occ_func_0_1(82)) +
          (occ_func_0_0(0) * occ_func_0_0(5) * occ_func_0_1(81)) +
          (occ_func_0_0(0) * occ_func_0_0(10) * occ_func_0_1(86)) +
          (occ_func_0_0(0) * occ_func_0_0(3) * occ_func_0_1(82)) +
          (occ_func_0_0(0) * occ_func_0_0(12) * occ_func_0_1(86)) +
          (occ_func_0_0(0) * occ_func_0_0(2) * occ_func_0_1(81)) +
          (occ_func_0_0(0) * occ_func_0_0(8) * occ_func_0_1(84)) +
          (occ_func_0_0(0) * occ_func_0_0(11) * occ_func_0_1(84)) +
          (occ_func_0_0(0) * occ_func_0_0(9) * occ_func_0_1(85)) +
          (occ_func_0_0(0) * occ_func_0_0(4) * occ_func_0_1(83)) +
          (occ_func_0_0(0) * occ_func_0_0(12) * occ_func_0_1(85)) +
          (occ_func_0_0(0) * occ_func_0_0(2) * occ_func_0_1(79)) +
          (occ_func_0_0(0) * occ_func_0_0(1) * occ_func_0_1(79)) +
          (occ_func_0_0(0) * occ_func_0_0(5) * occ_func_0_1(83)) +
          (occ_func_0_0(0) * occ_func_0_0(11) * occ_func_0_1(86)) +
          (occ_func_0_0(0) * occ_func_0_0(4) * occ_func_0_1(80)) +
          (occ_func_0_0(0) * occ_func_0_0(8) * occ_func_0_1(82)) +
          (occ_func_0_0(0) * occ_func_0_0(3) * occ_func_0_1(79)) +
          (occ_func_0_0(0) * occ_func_0_0(6) * occ_func_0_1(84)) +
          (occ_func_0_0(0) * occ_func_0_0(7) * occ_func_0_1(81)) +
          (occ_func_0_0(0) * occ_func_0_0(6) * occ_func_0_1(80)) +
          (occ_func_0_0(0) * occ_func_0_0(10) * occ_func_0_1(83)) +
          (occ_func_0_0(0) * occ_func_0_0(1) * occ_func_0_1(80))) /
         24.0;
}
double test_Clexulator::eval_bfunc_3_3_5() const {
  return ((occ_func_0_1(0) * occ_func_0_0(7) * occ_func_0_1(85)) +
          (occ_func_0_1(0) * occ_func_0_0(9) * occ_func_0_1(82)) +
          (occ_func_0_1(0) * occ_func_0_0(5) * occ_func_0_1(81)) +
          (occ_func_0_1(0) * occ_func_0_0(10) * occ_func_0_1(86)) +
          (occ_func_0_1(0) * occ_func_0_0(3) * occ_func_0_1(82)) +
          (occ_func_0_1(0) * occ_func_0_0(12) * occ_func_0_1(86)) +
          (occ_func_0_1(0) * occ_func_0_0(2) * occ_func_0_1(81)) +
          (occ_func_0_1(0) * occ_func_0_0(8) * occ_func_0_1(84)) +
          (occ_func_0_1(0) * occ_func_0_0(11) * occ_func_0_1(84)) +
          (occ_func_0_1(0) * occ_func_0_0(9) * occ_func_0_1(85)) +
          (occ_func_0_1(0) * occ_func_0_0(4) * occ_func_0_1(83)) +
          (occ_func_0_1(0) * occ_func_0_0(12) * occ_func_0_1(85)) +
          (occ_func_0_1(0) * occ_func_0_0(2) * occ_func_0_1(79)) +
          (occ_func_0_1(0) * occ_func_0_0(1) * occ_func_0_1(79)) +
          (occ_func_0_1(0) * occ_func_0_0(5) * occ_func_0_1(83)) +
          (occ_func_0_1(0) * occ_func_0_0(11) * occ_func_0_1(86)) +
          (occ_func_0_1(0) * occ_func_0_0(4) * occ_func_0_1(80)) +
          (occ_func_0_1(0) * occ_func_0_0(8) * occ_func_0_1(82)) +
          (occ_func_0_1(0) * occ_func_0_0(3) * occ_func_0_1(79)) +
          (occ_func_0_1(0) * occ_func_0_0(6) * occ_func_0_1(84)) +
          (occ_func_0_1(0) * occ_func_0_0(7) * occ_func_0_1(81)) +
          (occ_func_0_1(0) * occ_func_0_0(6) * occ_func_0_1(80)) +
          (occ_func_0_1(0) * occ_func_0_0(10) * occ_func_0_1(83)) +
          (occ_func_0_1(0) * occ_func_0_0(1) * occ_func_0_1(80))) /
         24.0;
}
double test_Clexulator::eval_bfunc_3_3_6() const {
  return ((occ_func_0_0(0) * occ_func_0_1(7) * occ_func_0_1(85)) +
          (occ_func_0_0(0) * occ_func_0_1(9) * occ_func_0_1(82)) +
          (occ_func_0_0(0) * occ_func_0_1(5) * occ_func_0_1(81)) +
          (occ_func_0_0(0) * occ_func_0_1(10) * occ_func_0_1(86)) +
          (occ_func_0_0(0) * occ_func_0_1(3) * occ_func_0_1(82)) +
          (occ_func_0_0(0) * occ_func_0_1(12) * occ_func_0_1(86)) +
          (occ_func_0_0(0) * occ_func_0_1(2) * occ_func_0_1(81)) +
          (occ_func_0_0(0) * occ_func_0_1(8) * occ_func_0_1(84)) +
          (occ_func_0_0(0) * occ_func_0_1(11) * occ_func_0_1(84)) +
          (occ_func_0_0(0) * occ_func_0_1(9) * occ_func_0_1(85)) +
          (occ_func_0_0(0) * occ_func_0_1(4) * occ_func_0_1(83)) +
          (occ_func_0_0(0) * occ_func_0_1(12) * occ_func_0_1(85)) +
          (occ_func_0_0(0) * occ_func_0_1(2) * occ_func_0_1(79)) +
          (occ_func_0_0(0) * occ_func_0_1(1) * occ_func_0_1(79)) +
          (occ_func_0_0(0) * occ_func_0_1(5) * occ_func_0_1(83)) +
          (occ_func_0_0(0) * occ_func_0_1(11) * occ_func_0_1(86)) +
          (occ_func_0_0(0) * occ_func_0_1(4) * occ_func_0_1(80)) +
          (occ_func_0_0(0) * occ_func_0_1(8) * occ_func_0_1(82)) +
          (occ_func_0_0(0) * occ_func_0_1(3) * occ_func_0_1(79)) +
          (occ_func_0_0(0) * occ_func_0_1(6) * occ_func_0_1(84)) +
          (occ_func_0_0(0) * occ_func_0_1(7) * occ_func_0_1(81)) +
          (occ_func_0_0(0) * occ_func_0_1(6) * occ_func_0_1(80)) +
          (occ_func_0_0(0) * occ_func_0_1(10) * occ_func_0_1(83)) +
          (occ_func_0_0(0) * occ_func_0_1(1) * occ_func_0_1(80))) /
         24.0;
}
double test_Clexulator::eval_bfunc_3_3_7() const {
  return ((occ_func_0_1(0) * occ_func_0_1(7) * occ_func_0_1(85)) +
          (occ_func_0_1(0) * occ_func_0_1(9) * occ_func_0_1(82)) +
          (occ_func_0_1(0) * occ_func_0_1(5) * occ_func_0_1(81)) +
          (occ_func_0_1(0) * occ_func_0_1(10) * occ_func_0_1(86)) +
          (occ_func_0_1(0) * occ_func_0_1(3) * occ_func_0_1(82)) +
          (occ_func_0_1(0) * occ_func_0_1(12) * occ_func_0_1(86)) +
          (occ_func_0_1(0) * occ_func_0_1(2) * occ_func_0_1(81)) +
          (occ_func_0_1(0) * occ_func_0_1(8) * occ_func_0_1(84)) +
          (occ_func_0_1(0) * occ_func_0_1(11) * occ_func_0_1(84)) +
          (occ_func_0_1(0) * occ_func_0_1(9) * occ_func_0_1(85)) +
          (occ_func_0_1(0) * occ_func_0_1(4) * occ_func_0_1(83)) +
          (occ_func_0_1(0) * occ_func_0_1(12) * occ_func_0_1(85)) +
          (occ_func_0_1(0) * occ_func_0_1(2) * occ_func_0_1(79)) +
          (occ_func_0_1(0) * occ_func_0_1(1) * occ_func_0_1(79)) +
          (occ_func_0_1(0) * occ_func_0_1(5) * occ_func_0_1(83)) +
          (occ_func_0_1(0) * occ_func_0_1(11) * occ_func_0_1(86)) +
          (occ_func_0_1(0) * occ_func_0_1(4) * occ_func_0_1(80)) +
          (occ_func_0_1(0) * occ_func_0_1(8) * occ_func_0_1(82)) +
          (occ_func_0_1(0) * occ_func_0_1(3) * occ_func_0_1(79)) +
          (occ_func_0_1(0) * occ_func_0_1(6) * occ_func_0_1(84)) +
          (occ_func_0_1(0) * occ_func_0_1(7) * occ_func_0_1(81)) +
          (occ_func_0_1(0) * occ_func_0_1(6) * occ_func_0_1(80)) +
          (occ_func_0_1(0) * occ_func_0_1(10) * occ_func_0_1(83)) +
          (occ_func_0_1(0) * occ_func_0_1(1) * occ_func_0_1(80))) /
         24.0;
}

double test_Clexulator::site_eval_at_0_bfunc_3_3_0() const {
  return ((occ_func_0_0(0) * occ_func_0_0(7) * occ_func_0_0(85)) +
          (occ_func_0_0(6) * occ_func_0_0(0) * occ_func_0_0(39)) +
          (occ_func_0_0(80) * occ_func_0_0(22) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(9) * occ_func_0_0(82)) +
          (occ_func_0_0(4) * occ_func_0_0(0) * occ_func_0_0(26)) +
          (occ_func_0_0(83) * occ_func_0_0(35) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(5) * occ_func_0_0(81)) +
          (occ_func_0_0(8) * occ_func_0_0(0) * occ_func_0_0(25)) +
          (occ_func_0_0(84) * occ_func_0_0(36) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(10) * occ_func_0_0(86)) +
          (occ_func_0_0(3) * occ_func_0_0(0) * occ_func_0_0(42)) +
          (occ_func_0_0(79) * occ_func_0_0(19) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(3) * occ_func_0_0(82)) +
          (occ_func_0_0(10) * occ_func_0_0(0) * occ_func_0_0(33)) +
          (occ_func_0_0(83) * occ_func_0_0(28) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(12) * occ_func_0_0(86)) +
          (occ_func_0_0(1) * occ_func_0_0(0) * occ_func_0_0(40)) +
          (occ_func_0_0(79) * occ_func_0_0(21) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(2) * occ_func_0_0(81)) +
          (occ_func_0_0(11) * occ_func_0_0(0) * occ_func_0_0(30)) +
          (occ_func_0_0(84) * occ_func_0_0(31) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(8) * occ_func_0_0(84)) +
          (occ_func_0_0(5) * occ_func_0_0(0) * occ_func_0_0(36)) +
          (occ_func_0_0(81) * occ_func_0_0(25) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(11) * occ_func_0_0(84)) +
          (occ_func_0_0(2) * occ_func_0_0(0) * occ_func_0_0(31)) +
          (occ_func_0_0(81) * occ_func_0_0(30) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(9) * occ_func_0_0(85)) +
          (occ_func_0_0(4) * occ_func_0_0(0) * occ_func_0_0(37)) +
          (occ_func_0_0(80) * occ_func_0_0(24) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(4) * occ_func_0_0(83)) +
          (occ_func_0_0(9) * occ_func_0_0(0) * occ_func_0_0(35)) +
          (occ_func_0_0(82) * occ_func_0_0(26) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(12) * occ_func_0_0(85)) +
          (occ_func_0_0(1) * occ_func_0_0(0) * occ_func_0_0(32)) +
          (occ_func_0_0(80) * occ_func_0_0(29) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(2) * occ_func_0_0(79)) +
          (occ_func_0_0(11) * occ_func_0_0(0) * occ_func_0_0(20)) +
          (occ_func_0_0(86) * occ_func_0_0(41) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(1) * occ_func_0_0(79)) +
          (occ_func_0_0(12) * occ_func_0_0(0) * occ_func_0_0(21)) +
          (occ_func_0_0(86) * occ_func_0_0(40) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(5) * occ_func_0_0(83)) +
          (occ_func_0_0(8) * occ_func_0_0(0) * occ_func_0_0(34)) +
          (occ_func_0_0(82) * occ_func_0_0(27) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(11) * occ_func_0_0(86)) +
          (occ_func_0_0(2) * occ_func_0_0(0) * occ_func_0_0(41)) +
          (occ_func_0_0(79) * occ_func_0_0(20) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(4) * occ_func_0_0(80)) +
          (occ_func_0_0(9) * occ_func_0_0(0) * occ_func_0_0(24)) +
          (occ_func_0_0(85) * occ_func_0_0(37) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(8) * occ_func_0_0(82)) +
          (occ_func_0_0(5) * occ_func_0_0(0) * occ_func_0_0(27)) +
          (occ_func_0_0(83) * occ_func_0_0(34) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(3) * occ_func_0_0(79)) +
          (occ_func_0_0(10) * occ_func_0_0(0) * occ_func_0_0(19)) +
          (occ_func_0_0(86) * occ_func_0_0(42) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(6) * occ_func_0_0(84)) +
          (occ_func_0_0(7) * occ_func_0_0(0) * occ_func_0_0(38)) +
          (occ_func_0_0(81) * occ_func_0_0(23) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(7) * occ_func_0_0(81)) +
          (occ_func_0_0(6) * occ_func_0_0(0) * occ_func_0_0(23)) +
          (occ_func_0_0(84) * occ_func_0_0(38) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(6) * occ_func_0_0(80)) +
          (occ_func_0_0(7) * occ_func_0_0(0) * occ_func_0_0(22)) +
          (occ_func_0_0(85) * occ_func_0_0(39) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(10) * occ_func_0_0(83)) +
          (occ_func_0_0(3) * occ_func_0_0(0) * occ_func_0_0(28)) +
          (occ_func_0_0(82) * occ_func_0_0(33) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(1) * occ_func_0_0(80)) +
          (occ_func_0_0(12) * occ_func_0_0(0) * occ_func_0_0(29)) +
          (occ_func_0_0(85) * occ_func_0_0(32) * occ_func_0_0(0))) /
         24.0;
}
double test_Clexulator::site_eval_at_0_bfunc_3_3_1() const {
  return ((occ_func_0_1(0) * occ_func_0_0(7) * occ_func_0_0(85)) +
          (occ_func_0_1(6) * occ_func_0_0(0) * occ_func_0_0(39)) +
          (occ_func_0_1(80) * occ_func_0_0(22) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_0(9) * occ_func_0_0(82)) +
          (occ_func_0_1(4) * occ_func_0_0(0) * occ_func_0_0(26)) +
          (occ_func_0_1(83) * occ_func_0_0(35) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_0(5) * occ_func_0_0(81)) +
          (occ_func_0_1(8) * occ_func_0_0(0) * occ_func_0_0(25)) +
          (occ_func_0_1(84) * occ_func_0_0(36) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_0(10) * occ_func_0_0(86)) +
          (occ_func_0_1(3) * occ_func_0_0(0) * occ_func_0_0(42)) +
          (occ_func_0_1(79) * occ_func_0_0(19) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_0(3) * occ_func_0_0(82)) +
          (occ_func_0_1(10) * occ_func_0_0(0) * occ_func_0_0(33)) +
          (occ_func_0_1(83) * occ_func_0_0(28) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_0(12) * occ_func_0_0(86)) +
          (occ_func_0_1(1) * occ_func_0_0(0) * occ_func_0_0(40)) +
          (occ_func_0_1(79) * occ_func_0_0(21) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_0(2) * occ_func_0_0(81)) +
          (occ_func_0_1(11) * occ_func_0_0(0) * occ_func_0_0(30)) +
          (occ_func_0_1(84) * occ_func_0_0(31) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_0(8) * occ_func_0_0(84)) +
          (occ_func_0_1(5) * occ_func_0_0(0) * occ_func_0_0(36)) +
          (occ_func_0_1(81) * occ_func_0_0(25) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_0(11) * occ_func_0_0(84)) +
          (occ_func_0_1(2) * occ_func_0_0(0) * occ_func_0_0(31)) +
          (occ_func_0_1(81) * occ_func_0_0(30) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_0(9) * occ_func_0_0(85)) +
          (occ_func_0_1(4) * occ_func_0_0(0) * occ_func_0_0(37)) +
          (occ_func_0_1(80) * occ_func_0_0(24) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_0(4) * occ_func_0_0(83)) +
          (occ_func_0_1(9) * occ_func_0_0(0) * occ_func_0_0(35)) +
          (occ_func_0_1(82) * occ_func_0_0(26) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_0(12) * occ_func_0_0(85)) +
          (occ_func_0_1(1) * occ_func_0_0(0) * occ_func_0_0(32)) +
          (occ_func_0_1(80) * occ_func_0_0(29) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_0(2) * occ_func_0_0(79)) +
          (occ_func_0_1(11) * occ_func_0_0(0) * occ_func_0_0(20)) +
          (occ_func_0_1(86) * occ_func_0_0(41) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_0(1) * occ_func_0_0(79)) +
          (occ_func_0_1(12) * occ_func_0_0(0) * occ_func_0_0(21)) +
          (occ_func_0_1(86) * occ_func_0_0(40) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_0(5) * occ_func_0_0(83)) +
          (occ_func_0_1(8) * occ_func_0_0(0) * occ_func_0_0(34)) +
          (occ_func_0_1(82) * occ_func_0_0(27) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_0(11) * occ_func_0_0(86)) +
          (occ_func_0_1(2) * occ_func_0_0(0) * occ_func_0_0(41)) +
          (occ_func_0_1(79) * occ_func_0_0(20) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_0(4) * occ_func_0_0(80)) +
          (occ_func_0_1(9) * occ_func_0_0(0) * occ_func_0_0(24)) +
          (occ_func_0_1(85) * occ_func_0_0(37) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_0(8) * occ_func_0_0(82)) +
          (occ_func_0_1(5) * occ_func_0_0(0) * occ_func_0_0(27)) +
          (occ_func_0_1(83) * occ_func_0_0(34) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_0(3) * occ_func_0_0(79)) +
          (occ_func_0_1(10) * occ_func_0_0(0) * occ_func_0_0(19)) +
          (occ_func_0_1(86) * occ_func_0_0(42) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_0(6) * occ_func_0_0(84)) +
          (occ_func_0_1(7) * occ_func_0_0(0) * occ_func_0_0(38)) +
          (occ_func_0_1(81) * occ_func_0_0(23) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_0(7) * occ_func_0_0(81)) +
          (occ_func_0_1(6) * occ_func_0_0(0) * occ_func_0_0(23)) +
          (occ_func_0_1(84) * occ_func_0_0(38) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_0(6) * occ_func_0_0(80)) +
          (occ_func_0_1(7) * occ_func_0_0(0) * occ_func_0_0(22)) +
          (occ_func_0_1(85) * occ_func_0_0(39) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_0(10) * occ_func_0_0(83)) +
          (occ_func_0_1(3) * occ_func_0_0(0) * occ_func_0_0(28)) +
          (occ_func_0_1(82) * occ_func_0_0(33) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_0(1) * occ_func_0_0(80)) +
          (occ_func_0_1(12) * occ_func_0_0(0) * occ_func_0_0(29)) +
          (occ_func_0_1(85) * occ_func_0_0(32) * occ_func_0_0(0))) /
         24.0;
}
double test_Clexulator::site_eval_at_0_bfunc_3_3_2() const {
  return ((occ_func_0_0(0) * occ_func_0_1(7) * occ_func_0_0(85)) +
          (occ_func_0_0(6) * occ_func_0_1(0) * occ_func_0_0(39)) +
          (occ_func_0_0(80) * occ_func_0_1(22) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_1(9) * occ_func_0_0(82)) +
          (occ_func_0_0(4) * occ_func_0_1(0) * occ_func_0_0(26)) +
          (occ_func_0_0(83) * occ_func_0_1(35) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_1(5) * occ_func_0_0(81)) +
          (occ_func_0_0(8) * occ_func_0_1(0) * occ_func_0_0(25)) +
          (occ_func_0_0(84) * occ_func_0_1(36) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_1(10) * occ_func_0_0(86)) +
          (occ_func_0_0(3) * occ_func_0_1(0) * occ_func_0_0(42)) +
          (occ_func_0_0(79) * occ_func_0_1(19) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_1(3) * occ_func_0_0(82)) +
          (occ_func_0_0(10) * occ_func_0_1(0) * occ_func_0_0(33)) +
          (occ_func_0_0(83) * occ_func_0_1(28) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_1(12) * occ_func_0_0(86)) +
          (occ_func_0_0(1) * occ_func_0_1(0) * occ_func_0_0(40)) +
          (occ_func_0_0(79) * occ_func_0_1(21) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_1(2) * occ_func_0_0(81)) +
          (occ_func_0_0(11) * occ_func_0_1(0) * occ_func_0_0(30)) +
          (occ_func_0_0(84) * occ_func_0_1(31) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_1(8) * occ_func_0_0(84)) +
          (occ_func_0_0(5) * occ_func_0_1(0) * occ_func_0_0(36)) +
          (occ_func_0_0(81) * occ_func_0_1(25) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_1(11) * occ_func_0_0(84)) +
          (occ_func_0_0(2) * occ_func_0_1(0) * occ_func_0_0(31)) +
          (occ_func_0_0(81) * occ_func_0_1(30) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_1(9) * occ_func_0_0(85)) +
          (occ_func_0_0(4) * occ_func_0_1(0) * occ_func_0_0(37)) +
          (occ_func_0_0(80) * occ_func_0_1(24) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_1(4) * occ_func_0_0(83)) +
          (occ_func_0_0(9) * occ_func_0_1(0) * occ_func_0_0(35)) +
          (occ_func_0_0(82) * occ_func_0_1(26) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_1(12) * occ_func_0_0(85)) +
          (occ_func_0_0(1) * occ_func_0_1(0) * occ_func_0_0(32)) +
          (occ_func_0_0(80) * occ_func_0_1(29) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_1(2) * occ_func_0_0(79)) +
          (occ_func_0_0(11) * occ_func_0_1(0) * occ_func_0_0(20)) +
          (occ_func_0_0(86) * occ_func_0_1(41) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_1(1) * occ_func_0_0(79)) +
          (occ_func_0_0(12) * occ_func_0_1(0) * occ_func_0_0(21)) +
          (occ_func_0_0(86) * occ_func_0_1(40) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_1(5) * occ_func_0_0(83)) +
          (occ_func_0_0(8) * occ_func_0_1(0) * occ_func_0_0(34)) +
          (occ_func_0_0(82) * occ_func_0_1(27) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_1(11) * occ_func_0_0(86)) +
          (occ_func_0_0(2) * occ_func_0_1(0) * occ_func_0_0(41)) +
          (occ_func_0_0(79) * occ_func_0_1(20) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_1(4) * occ_func_0_0(80)) +
          (occ_func_0_0(9) * occ_func_0_1(0) * occ_func_0_0(24)) +
          (occ_func_0_0(85) * occ_func_0_1(37) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_1(8) * occ_func_0_0(82)) +
          (occ_func_0_0(5) * occ_func_0_1(0) * occ_func_0_0(27)) +
          (occ_func_0_0(83) * occ_func_0_1(34) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_1(3) * occ_func_0_0(79)) +
          (occ_func_0_0(10) * occ_func_0_1(0) * occ_func_0_0(19)) +
          (occ_func_0_0(86) * occ_func_0_1(42) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_1(6) * occ_func_0_0(84)) +
          (occ_func_0_0(7) * occ_func_0_1(0) * occ_func_0_0(38)) +
          (occ_func_0_0(81) * occ_func_0_1(23) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_1(7) * occ_func_0_0(81)) +
          (occ_func_0_0(6) * occ_func_0_1(0) * occ_func_0_0(23)) +
          (occ_func_0_0(84) * occ_func_0_1(38) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_1(6) * occ_func_0_0(80)) +
          (occ_func_0_0(7) * occ_func_0_1(0) * occ_func_0_0(22)) +
          (occ_func_0_0(85) * occ_func_0_1(39) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_1(10) * occ_func_0_0(83)) +
          (occ_func_0_0(3) * occ_func_0_1(0) * occ_func_0_0(28)) +
          (occ_func_0_0(82) * occ_func_0_1(33) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_1(1) * occ_func_0_0(80)) +
          (occ_func_0_0(12) * occ_func_0_1(0) * occ_func_0_0(29)) +
          (occ_func_0_0(85) * occ_func_0_1(32) * occ_func_0_0(0))) /
         24.0;
}
double test_Clexulator::site_eval_at_0_bfunc_3_3_3() const {
  return ((occ_func_0_1(0) * occ_func_0_1(7) * occ_func_0_0(85)) +
          (occ_func_0_1(6) * occ_func_0_1(0) * occ_func_0_0(39)) +
          (occ_func_0_1(80) * occ_func_0_1(22) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_1(9) * occ_func_0_0(82)) +
          (occ_func_0_1(4) * occ_func_0_1(0) * occ_func_0_0(26)) +
          (occ_func_0_1(83) * occ_func_0_1(35) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_1(5) * occ_func_0_0(81)) +
          (occ_func_0_1(8) * occ_func_0_1(0) * occ_func_0_0(25)) +
          (occ_func_0_1(84) * occ_func_0_1(36) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_1(10) * occ_func_0_0(86)) +
          (occ_func_0_1(3) * occ_func_0_1(0) * occ_func_0_0(42)) +
          (occ_func_0_1(79) * occ_func_0_1(19) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_1(3) * occ_func_0_0(82)) +
          (occ_func_0_1(10) * occ_func_0_1(0) * occ_func_0_0(33)) +
          (occ_func_0_1(83) * occ_func_0_1(28) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_1(12) * occ_func_0_0(86)) +
          (occ_func_0_1(1) * occ_func_0_1(0) * occ_func_0_0(40)) +
          (occ_func_0_1(79) * occ_func_0_1(21) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_1(2) * occ_func_0_0(81)) +
          (occ_func_0_1(11) * occ_func_0_1(0) * occ_func_0_0(30)) +
          (occ_func_0_1(84) * occ_func_0_1(31) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_1(8) * occ_func_0_0(84)) +
          (occ_func_0_1(5) * occ_func_0_1(0) * occ_func_0_0(36)) +
          (occ_func_0_1(81) * occ_func_0_1(25) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_1(11) * occ_func_0_0(84)) +
          (occ_func_0_1(2) * occ_func_0_1(0) * occ_func_0_0(31)) +
          (occ_func_0_1(81) * occ_func_0_1(30) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_1(9) * occ_func_0_0(85)) +
          (occ_func_0_1(4) * occ_func_0_1(0) * occ_func_0_0(37)) +
          (occ_func_0_1(80) * occ_func_0_1(24) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_1(4) * occ_func_0_0(83)) +
          (occ_func_0_1(9) * occ_func_0_1(0) * occ_func_0_0(35)) +
          (occ_func_0_1(82) * occ_func_0_1(26) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_1(12) * occ_func_0_0(85)) +
          (occ_func_0_1(1) * occ_func_0_1(0) * occ_func_0_0(32)) +
          (occ_func_0_1(80) * occ_func_0_1(29) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_1(2) * occ_func_0_0(79)) +
          (occ_func_0_1(11) * occ_func_0_1(0) * occ_func_0_0(20)) +
          (occ_func_0_1(86) * occ_func_0_1(41) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_1(1) * occ_func_0_0(79)) +
          (occ_func_0_1(12) * occ_func_0_1(0) * occ_func_0_0(21)) +
          (occ_func_0_1(86) * occ_func_0_1(40) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_1(5) * occ_func_0_0(83)) +
          (occ_func_0_1(8) * occ_func_0_1(0) * occ_func_0_0(34)) +
          (occ_func_0_1(82) * occ_func_0_1(27) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_1(11) * occ_func_0_0(86)) +
          (occ_func_0_1(2) * occ_func_0_1(0) * occ_func_0_0(41)) +
          (occ_func_0_1(79) * occ_func_0_1(20) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_1(4) * occ_func_0_0(80)) +
          (occ_func_0_1(9) * occ_func_0_1(0) * occ_func_0_0(24)) +
          (occ_func_0_1(85) * occ_func_0_1(37) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_1(8) * occ_func_0_0(82)) +
          (occ_func_0_1(5) * occ_func_0_1(0) * occ_func_0_0(27)) +
          (occ_func_0_1(83) * occ_func_0_1(34) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_1(3) * occ_func_0_0(79)) +
          (occ_func_0_1(10) * occ_func_0_1(0) * occ_func_0_0(19)) +
          (occ_func_0_1(86) * occ_func_0_1(42) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_1(6) * occ_func_0_0(84)) +
          (occ_func_0_1(7) * occ_func_0_1(0) * occ_func_0_0(38)) +
          (occ_func_0_1(81) * occ_func_0_1(23) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_1(7) * occ_func_0_0(81)) +
          (occ_func_0_1(6) * occ_func_0_1(0) * occ_func_0_0(23)) +
          (occ_func_0_1(84) * occ_func_0_1(38) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_1(6) * occ_func_0_0(80)) +
          (occ_func_0_1(7) * occ_func_0_1(0) * occ_func_0_0(22)) +
          (occ_func_0_1(85) * occ_func_0_1(39) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_1(10) * occ_func_0_0(83)) +
          (occ_func_0_1(3) * occ_func_0_1(0) * occ_func_0_0(28)) +
          (occ_func_0_1(82) * occ_func_0_1(33) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_1(1) * occ_func_0_0(80)) +
          (occ_func_0_1(12) * occ_func_0_1(0) * occ_func_0_0(29)) +
          (occ_func_0_1(85) * occ_func_0_1(32) * occ_func_0_0(0))) /
         24.0;
}
double test_Clexulator::site_eval_at_0_bfunc_3_3_4() const {
  return ((occ_func_0_0(0) * occ_func_0_0(7) * occ_func_0_1(85)) +
          (occ_func_0_0(6) * occ_func_0_0(0) * occ_func_0_1(39)) +
          (occ_func_0_0(80) * occ_func_0_0(22) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_0(9) * occ_func_0_1(82)) +
          (occ_func_0_0(4) * occ_func_0_0(0) * occ_func_0_1(26)) +
          (occ_func_0_0(83) * occ_func_0_0(35) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_0(5) * occ_func_0_1(81)) +
          (occ_func_0_0(8) * occ_func_0_0(0) * occ_func_0_1(25)) +
          (occ_func_0_0(84) * occ_func_0_0(36) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_0(10) * occ_func_0_1(86)) +
          (occ_func_0_0(3) * occ_func_0_0(0) * occ_func_0_1(42)) +
          (occ_func_0_0(79) * occ_func_0_0(19) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_0(3) * occ_func_0_1(82)) +
          (occ_func_0_0(10) * occ_func_0_0(0) * occ_func_0_1(33)) +
          (occ_func_0_0(83) * occ_func_0_0(28) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_0(12) * occ_func_0_1(86)) +
          (occ_func_0_0(1) * occ_func_0_0(0) * occ_func_0_1(40)) +
          (occ_func_0_0(79) * occ_func_0_0(21) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_0(2) * occ_func_0_1(81)) +
          (occ_func_0_0(11) * occ_func_0_0(0) * occ_func_0_1(30)) +
          (occ_func_0_0(84) * occ_func_0_0(31) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_0(8) * occ_func_0_1(84)) +
          (occ_func_0_0(5) * occ_func_0_0(0) * occ_func_0_1(36)) +
          (occ_func_0_0(81) * occ_func_0_0(25) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_0(11) * occ_func_0_1(84)) +
          (occ_func_0_0(2) * occ_func_0_0(0) * occ_func_0_1(31)) +
          (occ_func_0_0(81) * occ_func_0_0(30) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_0(9) * occ_func_0_1(85)) +
          (occ_func_0_0(4) * occ_func_0_0(0) * occ_func_0_1(37)) +
          (occ_func_0_0(80) * occ_func_0_0(24) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_0(4) * occ_func_0_1(83)) +
          (occ_func_0_0(9) * occ_func_0_0(0) * occ_func_0_1(35)) +
          (occ_func_0_0(82) * occ_func_0_0(26) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_0(12) * occ_func_0_1(85)) +
          (occ_func_0_0(1) * occ_func_0_0(0) * occ_func_0_1(32)) +
          (occ_func_0_0(80) * occ_func_0_0(29) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_0(2) * occ_func_0_1(79)) +
          (occ_func_0_0(11) * occ_func_0_0(0) * occ_func_0_1(20)) +
          (occ_func_0_0(86) * occ_func_0_0(41) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_0(1) * occ_func_0_1(79)) +
          (occ_func_0_0(12) * occ_func_0_0(0) * occ_func_0_1(21)) +
          (occ_func_0_0(86) * occ_func_0_0(40) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_0(5) * occ_func_0_1(83)) +
          (occ_func_0_0(8) * occ_func_0_0(0) * occ_func_0_1(34)) +
          (occ_func_0_0(82) * occ_func_0_0(27) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_0(11) * occ_func_0_1(86)) +
          (occ_func_0_0(2) * occ_func_0_0(0) * occ_func_0_1(41)) +
          (occ_func_0_0(79) * occ_func_0_0(20) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_0(4) * occ_func_0_1(80)) +
          (occ_func_0_0(9) * occ_func_0_0(0) * occ_func_0_1(24)) +
          (occ_func_0_0(85) * occ_func_0_0(37) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_0(8) * occ_func_0_1(82)) +
          (occ_func_0_0(5) * occ_func_0_0(0) * occ_func_0_1(27)) +
          (occ_func_0_0(83) * occ_func_0_0(34) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_0(3) * occ_func_0_1(79)) +
          (occ_func_0_0(10) * occ_func_0_0(0) * occ_func_0_1(19)) +
          (occ_func_0_0(86) * occ_func_0_0(42) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_0(6) * occ_func_0_1(84)) +
          (occ_func_0_0(7) * occ_func_0_0(0) * occ_func_0_1(38)) +
          (occ_func_0_0(81) * occ_func_0_0(23) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_0(7) * occ_func_0_1(81)) +
          (occ_func_0_0(6) * occ_func_0_0(0) * occ_func_0_1(23)) +
          (occ_func_0_0(84) * occ_func_0_0(38) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_0(6) * occ_func_0_1(80)) +
          (occ_func_0_0(7) * occ_func_0_0(0) * occ_func_0_1(22)) +
          (occ_func_0_0(85) * occ_func_0_0(39) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_0(10) * occ_func_0_1(83)) +
          (occ_func_0_0(3) * occ_func_0_0(0) * occ_func_0_1(28)) +
          (occ_func_0_0(82) * occ_func_0_0(33) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_0(1) * occ_func_0_1(80)) +
          (occ_func_0_0(12) * occ_func_0_0(0) * occ_func_0_1(29)) +
          (occ_func_0_0(85) * occ_func_0_0(32) * occ_func_0_1(0))) /
         24.0;
}
double test_Clexulator::site_eval_at_0_bfunc_3_3_5() const {
  return ((occ_func_0_1(0) * occ_func_0_0(7) * occ_func_0_1(85)) +
          (occ_func_0_1(6) * occ_func_0_0(0) * occ_func_0_1(39)) +
          (occ_func_0_1(80) * occ_func_0_0(22) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_0(9) * occ_func_0_1(82)) +
          (occ_func_0_1(4) * occ_func_0_0(0) * occ_func_0_1(26)) +
          (occ_func_0_1(83) * occ_func_0_0(35) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_0(5) * occ_func_0_1(81)) +
          (occ_func_0_1(8) * occ_func_0_0(0) * occ_func_0_1(25)) +
          (occ_func_0_1(84) * occ_func_0_0(36) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_0(10) * occ_func_0_1(86)) +
          (occ_func_0_1(3) * occ_func_0_0(0) * occ_func_0_1(42)) +
          (occ_func_0_1(79) * occ_func_0_0(19) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_0(3) * occ_func_0_1(82)) +
          (occ_func_0_1(10) * occ_func_0_0(0) * occ_func_0_1(33)) +
          (occ_func_0_1(83) * occ_func_0_0(28) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_0(12) * occ_func_0_1(86)) +
          (occ_func_0_1(1) * occ_func_0_0(0) * occ_func_0_1(40)) +
          (occ_func_0_1(79) * occ_func_0_0(21) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_0(2) * occ_func_0_1(81)) +
          (occ_func_0_1(11) * occ_func_0_0(0) * occ_func_0_1(30)) +
          (occ_func_0_1(84) * occ_func_0_0(31) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_0(8) * occ_func_0_1(84)) +
          (occ_func_0_1(5) * occ_func_0_0(0) * occ_func_0_1(36)) +
          (occ_func_0_1(81) * occ_func_0_0(25) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_0(11) * occ_func_0_1(84)) +
          (occ_func_0_1(2) * occ_func_0_0(0) * occ_func_0_1(31)) +
          (occ_func_0_1(81) * occ_func_0_0(30) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_0(9) * occ_func_0_1(85)) +
          (occ_func_0_1(4) * occ_func_0_0(0) * occ_func_0_1(37)) +
          (occ_func_0_1(80) * occ_func_0_0(24) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_0(4) * occ_func_0_1(83)) +
          (occ_func_0_1(9) * occ_func_0_0(0) * occ_func_0_1(35)) +
          (occ_func_0_1(82) * occ_func_0_0(26) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_0(12) * occ_func_0_1(85)) +
          (occ_func_0_1(1) * occ_func_0_0(0) * occ_func_0_1(32)) +
          (occ_func_0_1(80) * occ_func_0_0(29) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_0(2) * occ_func_0_1(79)) +
          (occ_func_0_1(11) * occ_func_0_0(0) * occ_func_0_1(20)) +
          (occ_func_0_1(86) * occ_func_0_0(41) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_0(1) * occ_func_0_1(79)) +
          (occ_func_0_1(12) * occ_func_0_0(0) * occ_func_0_1(21)) +
          (occ_func_0_1(86) * occ_func_0_0(40) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_0(5) * occ_func_0_1(83)) +
          (occ_func_0_1(8) * occ_func_0_0(0) * occ_func_0_1(34)) +
          (occ_func_0_1(82) * occ_func_0_0(27) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_0(11) * occ_func_0_1(86)) +
          (occ_func_0_1(2) * occ_func_0_0(0) * occ_func_0_1(41)) +
          (occ_func_0_1(79) * occ_func_0_0(20) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_0(4) * occ_func_0_1(80)) +
          (occ_func_0_1(9) * occ_func_0_0(0) * occ_func_0_1(24)) +
          (occ_func_0_1(85) * occ_func_0_0(37) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_0(8) * occ_func_0_1(82)) +
          (occ_func_0_1(5) * occ_func_0_0(0) * occ_func_0_1(27)) +
          (occ_func_0_1(83) * occ_func_0_0(34) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_0(3) * occ_func_0_1(79)) +
          (occ_func_0_1(10) * occ_func_0_0(0) * occ_func_0_1(19)) +
          (occ_func_0_1(86) * occ_func_0_0(42) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_0(6) * occ_func_0_1(84)) +
          (occ_func_0_1(7) * occ_func_0_0(0) * occ_func_0_1(38)) +
          (occ_func_0_1(81) * occ_func_0_0(23) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_0(7) * occ_func_0_1(81)) +
          (occ_func_0_1(6) * occ_func_0_0(0) * occ_func_0_1(23)) +
          (occ_func_0_1(84) * occ_func_0_0(38) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_0(6) * occ_func_0_1(80)) +
          (occ_func_0_1(7) * occ_func_0_0(0) * occ_func_0_1(22)) +
          (occ_func_0_1(85) * occ_func_0_0(39) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_0(10) * occ_func_0_1(83)) +
          (occ_func_0_1(3) * occ_func_0_0(0) * occ_func_0_1(28)) +
          (occ_func_0_1(82) * occ_func_0_0(33) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_0(1) * occ_func_0_1(80)) +
          (occ_func_0_1(12) * occ_func_0_0(0) * occ_func_0_1(29)) +
          (occ_func_0_1(85) * occ_func_0_0(32) * occ_func_0_1(0))) /
         24.0;
}
double test_Clexulator::site_eval_at_0_bfunc_3_3_6() const {
  return ((occ_func_0_0(0) * occ_func_0_1(7) * occ_func_0_1(85)) +
          (occ_func_0_0(6) * occ_func_0_1(0) * occ_func_0_1(39)) +
          (occ_func_0_0(80) * occ_func_0_1(22) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_1(9) * occ_func_0_1(82)) +
          (occ_func_0_0(4) * occ_func_0_1(0) * occ_func_0_1(26)) +
          (occ_func_0_0(83) * occ_func_0_1(35) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_1(5) * occ_func_0_1(81)) +
          (occ_func_0_0(8) * occ_func_0_1(0) * occ_func_0_1(25)) +
          (occ_func_0_0(84) * occ_func_0_1(36) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_1(10) * occ_func_0_1(86)) +
          (occ_func_0_0(3) * occ_func_0_1(0) * occ_func_0_1(42)) +
          (occ_func_0_0(79) * occ_func_0_1(19) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_1(3) * occ_func_0_1(82)) +
          (occ_func_0_0(10) * occ_func_0_1(0) * occ_func_0_1(33)) +
          (occ_func_0_0(83) * occ_func_0_1(28) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_1(12) * occ_func_0_1(86)) +
          (occ_func_0_0(1) * occ_func_0_1(0) * occ_func_0_1(40)) +
          (occ_func_0_0(79) * occ_func_0_1(21) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_1(2) * occ_func_0_1(81)) +
          (occ_func_0_0(11) * occ_func_0_1(0) * occ_func_0_1(30)) +
          (occ_func_0_0(84) * occ_func_0_1(31) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_1(8) * occ_func_0_1(84)) +
          (occ_func_0_0(5) * occ_func_0_1(0) * occ_func_0_1(36)) +
          (occ_func_0_0(81) * occ_func_0_1(25) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_1(11) * occ_func_0_1(84)) +
          (occ_func_0_0(2) * occ_func_0_1(0) * occ_func_0_1(31)) +
          (occ_func_0_0(81) * occ_func_0_1(30) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_1(9) * occ_func_0_1(85)) +
          (occ_func_0_0(4) * occ_func_0_1(0) * occ_func_0_1(37)) +
          (occ_func_0_0(80) * occ_func_0_1(24) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_1(4) * occ_func_0_1(83)) +
          (occ_func_0_0(9) * occ_func_0_1(0) * occ_func_0_1(35)) +
          (occ_func_0_0(82) * occ_func_0_1(26) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_1(12) * occ_func_0_1(85)) +
          (occ_func_0_0(1) * occ_func_0_1(0) * occ_func_0_1(32)) +
          (occ_func_0_0(80) * occ_func_0_1(29) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_1(2) * occ_func_0_1(79)) +
          (occ_func_0_0(11) * occ_func_0_1(0) * occ_func_0_1(20)) +
          (occ_func_0_0(86) * occ_func_0_1(41) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_1(1) * occ_func_0_1(79)) +
          (occ_func_0_0(12) * occ_func_0_1(0) * occ_func_0_1(21)) +
          (occ_func_0_0(86) * occ_func_0_1(40) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_1(5) * occ_func_0_1(83)) +
          (occ_func_0_0(8) * occ_func_0_1(0) * occ_func_0_1(34)) +
          (occ_func_0_0(82) * occ_func_0_1(27) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_1(11) * occ_func_0_1(86)) +
          (occ_func_0_0(2) * occ_func_0_1(0) * occ_func_0_1(41)) +
          (occ_func_0_0(79) * occ_func_0_1(20) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_1(4) * occ_func_0_1(80)) +
          (occ_func_0_0(9) * occ_func_0_1(0) * occ_func_0_1(24)) +
          (occ_func_0_0(85) * occ_func_0_1(37) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_1(8) * occ_func_0_1(82)) +
          (occ_func_0_0(5) * occ_func_0_1(0) * occ_func_0_1(27)) +
          (occ_func_0_0(83) * occ_func_0_1(34) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_1(3) * occ_func_0_1(79)) +
          (occ_func_0_0(10) * occ_func_0_1(0) * occ_func_0_1(19)) +
          (occ_func_0_0(86) * occ_func_0_1(42) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_1(6) * occ_func_0_1(84)) +
          (occ_func_0_0(7) * occ_func_0_1(0) * occ_func_0_1(38)) +
          (occ_func_0_0(81) * occ_func_0_1(23) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_1(7) * occ_func_0_1(81)) +
          (occ_func_0_0(6) * occ_func_0_1(0) * occ_func_0_1(23)) +
          (occ_func_0_0(84) * occ_func_0_1(38) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_1(6) * occ_func_0_1(80)) +
          (occ_func_0_0(7) * occ_func_0_1(0) * occ_func_0_1(22)) +
          (occ_func_0_0(85) * occ_func_0_1(39) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_1(10) * occ_func_0_1(83)) +
          (occ_func_0_0(3) * occ_func_0_1(0) * occ_func_0_1(28)) +
          (occ_func_0_0(82) * occ_func_0_1(33) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_1(1) * occ_func_0_1(80)) +
          (occ_func_0_0(12) * occ_func_0_1(0) * occ_func_0_1(29)) +
          (occ_func_0_0(85) * occ_func_0_1(32) * occ_func_0_1(0))) /
         24.0;
}
double test_Clexulator::site_eval_at_0_bfunc_3_3_7() const {
  return ((occ_func_0_1(0) * occ_func_0_1(7) * occ_func_0_1(85)) +
          (occ_func_0_1(6) * occ_func_0_1(0) * occ_func_0_1(39)) +
          (occ_func_0_1(80) * occ_func_0_1(22) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(9) * occ_func_0_1(82)) +
          (occ_func_0_1(4) * occ_func_0_1(0) * occ_func_0_1(26)) +
          (occ_func_0_1(83) * occ_func_0_1(35) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(5) * occ_func_0_1(81)) +
          (occ_func_0_1(8) * occ_func_0_1(0) * occ_func_0_1(25)) +
          (occ_func_0_1(84) * occ_func_0_1(36) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(10) * occ_func_0_1(86)) +
          (occ_func_0_1(3) * occ_func_0_1(0) * occ_func_0_1(42)) +
          (occ_func_0_1(79) * occ_func_0_1(19) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(3) * occ_func_0_1(82)) +
          (occ_func_0_1(10) * occ_func_0_1(0) * occ_func_0_1(33)) +
          (occ_func_0_1(83) * occ_func_0_1(28) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(12) * occ_func_0_1(86)) +
          (occ_func_0_1(1) * occ_func_0_1(0) * occ_func_0_1(40)) +
          (occ_func_0_1(79) * occ_func_0_1(21) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(2) * occ_func_0_1(81)) +
          (occ_func_0_1(11) * occ_func_0_1(0) * occ_func_0_1(30)) +
          (occ_func_0_1(84) * occ_func_0_1(31) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(8) * occ_func_0_1(84)) +
          (occ_func_0_1(5) * occ_func_0_1(0) * occ_func_0_1(36)) +
          (occ_func_0_1(81) * occ_func_0_1(25) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(11) * occ_func_0_1(84)) +
          (occ_func_0_1(2) * occ_func_0_1(0) * occ_func_0_1(31)) +
          (occ_func_0_1(81) * occ_func_0_1(30) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(9) * occ_func_0_1(85)) +
          (occ_func_0_1(4) * occ_func_0_1(0) * occ_func_0_1(37)) +
          (occ_func_0_1(80) * occ_func_0_1(24) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(4) * occ_func_0_1(83)) +
          (occ_func_0_1(9) * occ_func_0_1(0) * occ_func_0_1(35)) +
          (occ_func_0_1(82) * occ_func_0_1(26) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(12) * occ_func_0_1(85)) +
          (occ_func_0_1(1) * occ_func_0_1(0) * occ_func_0_1(32)) +
          (occ_func_0_1(80) * occ_func_0_1(29) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(2) * occ_func_0_1(79)) +
          (occ_func_0_1(11) * occ_func_0_1(0) * occ_func_0_1(20)) +
          (occ_func_0_1(86) * occ_func_0_1(41) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(1) * occ_func_0_1(79)) +
          (occ_func_0_1(12) * occ_func_0_1(0) * occ_func_0_1(21)) +
          (occ_func_0_1(86) * occ_func_0_1(40) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(5) * occ_func_0_1(83)) +
          (occ_func_0_1(8) * occ_func_0_1(0) * occ_func_0_1(34)) +
          (occ_func_0_1(82) * occ_func_0_1(27) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(11) * occ_func_0_1(86)) +
          (occ_func_0_1(2) * occ_func_0_1(0) * occ_func_0_1(41)) +
          (occ_func_0_1(79) * occ_func_0_1(20) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(4) * occ_func_0_1(80)) +
          (occ_func_0_1(9) * occ_func_0_1(0) * occ_func_0_1(24)) +
          (occ_func_0_1(85) * occ_func_0_1(37) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(8) * occ_func_0_1(82)) +
          (occ_func_0_1(5) * occ_func_0_1(0) * occ_func_0_1(27)) +
          (occ_func_0_1(83) * occ_func_0_1(34) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(3) * occ_func_0_1(79)) +
          (occ_func_0_1(10) * occ_func_0_1(0) * occ_func_0_1(19)) +
          (occ_func_0_1(86) * occ_func_0_1(42) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(6) * occ_func_0_1(84)) +
          (occ_func_0_1(7) * occ_func_0_1(0) * occ_func_0_1(38)) +
          (occ_func_0_1(81) * occ_func_0_1(23) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(7) * occ_func_0_1(81)) +
          (occ_func_0_1(6) * occ_func_0_1(0) * occ_func_0_1(23)) +
          (occ_func_0_1(84) * occ_func_0_1(38) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(6) * occ_func_0_1(80)) +
          (occ_func_0_1(7) * occ_func_0_1(0) * occ_func_0_1(22)) +
          (occ_func_0_1(85) * occ_func_0_1(39) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(10) * occ_func_0_1(83)) +
          (occ_func_0_1(3) * occ_func_0_1(0) * occ_func_0_1(28)) +
          (occ_func_0_1(82) * occ_func_0_1(33) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(1) * occ_func_0_1(80)) +
          (occ_func_0_1(12) * occ_func_0_1(0) * occ_func_0_1(29)) +
          (occ_func_0_1(85) * occ_func_0_1(32) * occ_func_0_1(0))) /
         24.0;
}

double test_Clexulator::delta_site_eval_at_0_bfunc_3_3_0(int occ_i,
                                                         int occ_f) const {
  return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) *
         ((occ_func_0_0(7) * occ_func_0_0(85)) +
          (occ_func_0_0(6) * occ_func_0_0(39)) +
          (occ_func_0_0(80) * occ_func_0_0(22)) +
          (occ_func_0_0(9) * occ_func_0_0(82)) +
          (occ_func_0_0(4) * occ_func_0_0(26)) +
          (occ_func_0_0(83) * occ_func_0_0(35)) +
          (occ_func_0_0(5) * occ_func_0_0(81)) +
          (occ_func_0_0(8) * occ_func_0_0(25)) +
          (occ_func_0_0(84) * occ_func_0_0(36)) +
          (occ_func_0_0(10) * occ_func_0_0(86)) +
          (occ_func_0_0(3) * occ_func_0_0(42)) +
          (occ_func_0_0(79) * occ_func_0_0(19)) +
          (occ_func_0_0(3) * occ_func_0_0(82)) +
          (occ_func_0_0(10) * occ_func_0_0(33)) +
          (occ_func_0_0(83) * occ_func_0_0(28)) +
          (occ_func_0_0(12) * occ_func_0_0(86)) +
          (occ_func_0_0(1) * occ_func_0_0(40)) +
          (occ_func_0_0(79) * occ_func_0_0(21)) +
          (occ_func_0_0(2) * occ_func_0_0(81)) +
          (occ_func_0_0(11) * occ_func_0_0(30)) +
          (occ_func_0_0(84) * occ_func_0_0(31)) +
          (occ_func_0_0(8) * occ_func_0_0(84)) +
          (occ_func_0_0(5) * occ_func_0_0(36)) +
          (occ_func_0_0(81) * occ_func_0_0(25)) +
          (occ_func_0_0(11) * occ_func_0_0(84)) +
          (occ_func_0_0(2) * occ_func_0_0(31)) +
          (occ_func_0_0(81) * occ_func_0_0(30)) +
          (occ_func_0_0(9) * occ_func_0_0(85)) +
          (occ_func_0_0(4) * occ_func_0_0(37)) +
          (occ_func_0_0(80) * occ_func_0_0(24)) +
          (occ_func_0_0(4) * occ_func_0_0(83)) +
          (occ_func_0_0(9) * occ_func_0_0(35)) +
          (occ_func_0_0(82) * occ_func_0_0(26)) +
          (occ_func_0_0(12) * occ_func_0_0(85)) +
          (occ_func_0_0(1) * occ_func_0_0(32)) +
          (occ_func_0_0(80) * occ_func_0_0(29)) +
          (occ_func_0_0(2) * occ_func_0_0(79)) +
          (occ_func_0_0(11) * occ_func_0_0(20)) +
          (occ_func_0_0(86) * occ_func_0_0(41)) +
          (occ_func_0_0(1) * occ_func_0_0(79)) +
          (occ_func_0_0(12) * occ_func_0_0(21)) +
          (occ_func_0_0(86) * occ_func_0_0(40)) +
          (occ_func_0_0(5) * occ_func_0_0(83)) +
          (occ_func_0_0(8) * occ_func_0_0(34)) +
          (occ_func_0_0(82) * occ_func_0_0(27)) +
          (occ_func_0_0(11) * occ_func_0_0(86)) +
          (occ_func_0_0(2) * occ_func_0_0(41)) +
          (occ_func_0_0(79) * occ_func_0_0(20)) +
          (occ_func_0_0(4) * occ_func_0_0(80)) +
          (occ_func_0_0(9) * occ_func_0_0(24)) +
          (occ_func_0_0(85) * occ_func_0_0(37)) +
          (occ_func_0_0(8) * occ_func_0_0(82)) +
          (occ_func_0_0(5) * occ_func_0_0(27)) +
          (occ_func_0_0(83) * occ_func_0_0(34)) +
          (occ_func_0_0(3) * occ_func_0_0(79)) +
          (occ_func_0_0(10) * occ_func_0_0(19)) +
          (occ_func_0_0(86) * occ_func_0_0(42)) +
          (occ_func_0_0(6) * occ_func_0_0(84)) +
          (occ_func_0_0(7) * occ_func_0_0(38)) +
          (occ_func_0_0(81) * occ_func_0_0(23)) +
          (occ_func_0_0(7) * occ_func_0_0(81)) +
          (occ_func_0_0(6) * occ_func_0_0(23)) +
          (occ_func_0_0(84) * occ_func_0_0(38)) +
          (occ_func_0_0(6) * occ_func_0_0(80)) +
          (occ_func_0_0(7) * occ_func_0_0(22)) +
          (occ_func_0_0(85) * occ_func_0_0(39)) +
          (occ_func_0_0(10) * occ_func_0_0(83)) +
          (occ_func_0_0(3) * occ_func_0_0(28)) +
          (occ_func_0_0(82) * occ_func_0_0(33)) +
          (occ_func_0_0(1) * occ_func_0_0(80)) +
          (occ_func_0_0(12) * occ_func_0_0(29)) +
          (occ_func_0_0(85) * occ_func_0_0(32))) /
         24.0;
}
double test_Clexulator::delta_site_eval_at_0_bfunc_3_3_1(int occ_i,
                                                         int occ_f) const {
  return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) *
             ((occ_func_0_1(6) * occ_func_0_0(39)) +
              (occ_func_0_1(80) * occ_func_0_0(22)) +
              (occ_func_0_1(4) * occ_func_0_0(26)) +
              (occ_func_0_1(83) * occ_func_0_0(35)) +
              (occ_func_0_1(8) * occ_func_0_0(25)) +
              (occ_func_0_1(84) * occ_func_0_0(36)) +
              (occ_func_0_1(3) * occ_func_0_0(42)) +
              (occ_func_0_1(79) * occ_func_0_0(19)) +
              (occ_func_0_1(10) * occ_func_0_0(33)) +
              (occ_func_0_1(83) * occ_func_0_0(28)) +
              (occ_func_0_1(1) * occ_func_0_0(40)) +
              (occ_func_0_1(79) * occ_func_0_0(21)) +
              (occ_func_0_1(11) * occ_func_0_0(30)) +
              (occ_func_0_1(84) * occ_func_0_0(31)) +
              (occ_func_0_1(5) * occ_func_0_0(36)) +
              (occ_func_0_1(81) * occ_func_0_0(25)) +
              (occ_func_0_1(2) * occ_func_0_0(31)) +
              (occ_func_0_1(81) * occ_func_0_0(30)) +
              (occ_func_0_1(4) * occ_func_0_0(37)) +
              (occ_func_0_1(80) * occ_func_0_0(24)) +
              (occ_func_0_1(9) * occ_func_0_0(35)) +
              (occ_func_0_1(82) * occ_func_0_0(26)) +
              (occ_func_0_1(1) * occ_func_0_0(32)) +
              (occ_func_0_1(80) * occ_func_0_0(29)) +
              (occ_func_0_1(11) * occ_func_0_0(20)) +
              (occ_func_0_1(86) * occ_func_0_0(41)) +
              (occ_func_0_1(12) * occ_func_0_0(21)) +
              (occ_func_0_1(86) * occ_func_0_0(40)) +
              (occ_func_0_1(8) * occ_func_0_0(34)) +
              (occ_func_0_1(82) * occ_func_0_0(27)) +
              (occ_func_0_1(2) * occ_func_0_0(41)) +
              (occ_func_0_1(79) * occ_func_0_0(20)) +
              (occ_func_0_1(9) * occ_func_0_0(24)) +
              (occ_func_0_1(85) * occ_func_0_0(37)) +
              (occ_func_0_1(5) * occ_func_0_0(27)) +
              (occ_func_0_1(83) * occ_func_0_0(34)) +
              (occ_func_0_1(10) * occ_func_0_0(19)) +
              (occ_func_0_1(86) * occ_func_0_0(42)) +
              (occ_func_0_1(7) * occ_func_0_0(38)) +
              (occ_func_0_1(81) * occ_func_0_0(23)) +
              (occ_func_0_1(6) * occ_func_0_0(23)) +
              (occ_func_0_1(84) * occ_func_0_0(38)) +
              (occ_func_0_1(7) * occ_func_0_0(22)) +
              (occ_func_0_1(85) * occ_func_0_0(39)) +
              (occ_func_0_1(3) * occ_func_0_0(28)) +
              (occ_func_0_1(82) * occ_func_0_0(33)) +
              (occ_func_0_1(12) * occ_func_0_0(29)) +
              (occ_func_0_1(85) * occ_func_0_0(32))) /
             24.0 +
         (m_occ_func_0_1[occ_f] - m_occ_func_0_1[occ_i]) *
             ((occ_func_0_0(7) * occ_func_0_0(85)) +
              (occ_func_0_0(9) * occ_func_0_0(82)) +
              (occ_func_0_0(5) * occ_func_0_0(81)) +
              (occ_func_0_0(10) * occ_func_0_0(86)) +
              (occ_func_0_0(3) * occ_func_0_0(82)) +
              (occ_func_0_0(12) * occ_func_0_0(86)) +
              (occ_func_0_0(2) * occ_func_0_0(81)) +
              (occ_func_0_0(8) * occ_func_0_0(84)) +
              (occ_func_0_0(11) * occ_func_0_0(84)) +
              (occ_func_0_0(9) * occ_func_0_0(85)) +
              (occ_func_0_0(4) * occ_func_0_0(83)) +
              (occ_func_0_0(12) * occ_func_0_0(85)) +
              (occ_func_0_0(2) * occ_func_0_0(79)) +
              (occ_func_0_0(1) * occ_func_0_0(79)) +
              (occ_func_0_0(5) * occ_func_0_0(83)) +
              (occ_func_0_0(11) * occ_func_0_0(86)) +
              (occ_func_0_0(4) * occ_func_0_0(80)) +
              (occ_func_0_0(8) * occ_func_0_0(82)) +
              (occ_func_0_0(3) * occ_func_0_0(79)) +
              (occ_func_0_0(6) * occ_func_0_0(84)) +
              (occ_func_0_0(7) * occ_func_0_0(81)) +
              (occ_func_0_0(6) * occ_func_0_0(80)) +
              (occ_func_0_0(10) * occ_func_0_0(83)) +
              (occ_func_0_0(1) * occ_func_0_0(80))) /
             24.0;
}
double test_Clexulator::delta_site_eval_at_0_bfunc_3_3_2(int occ_i,
                                                         int occ_f) const {
  return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) *
             ((occ_func_0_1(7) * occ_func_0_0(85)) +
              (occ_func_0_0(80) * occ_func_0_1(22)) +
              (occ_func_0_1(9) * occ_func_0_0(82)) +
              (occ_func_0_0(83) * occ_func_0_1(35)) +
              (occ_func_0_1(5) * occ_func_0_0(81)) +
              (occ_func_0_0(84) * occ_func_0_1(36)) +
              (occ_func_0_1(10) * occ_func_0_0(86)) +
              (occ_func_0_0(79) * occ_func_0_1(19)) +
              (occ_func_0_1(3) * occ_func_0_0(82)) +
              (occ_func_0_0(83) * occ_func_0_1(28)) +
              (occ_func_0_1(12) * occ_func_0_0(86)) +
              (occ_func_0_0(79) * occ_func_0_1(21)) +
              (occ_func_0_1(2) * occ_func_0_0(81)) +
              (occ_func_0_0(84) * occ_func_0_1(31)) +
              (occ_func_0_1(8) * occ_func_0_0(84)) +
              (occ_func_0_0(81) * occ_func_0_1(25)) +
              (occ_func_0_1(11) * occ_func_0_0(84)) +
              (occ_func_0_0(81) * occ_func_0_1(30)) +
              (occ_func_0_1(9) * occ_func_0_0(85)) +
              (occ_func_0_0(80) * occ_func_0_1(24)) +
              (occ_func_0_1(4) * occ_func_0_0(83)) +
              (occ_func_0_0(82) * occ_func_0_1(26)) +
              (occ_func_0_1(12) * occ_func_0_0(85)) +
              (occ_func_0_0(80) * occ_func_0_1(29)) +
              (occ_func_0_1(2) * occ_func_0_0(79)) +
              (occ_func_0_0(86) * occ_func_0_1(41)) +
              (occ_func_0_1(1) * occ_func_0_0(79)) +
              (occ_func_0_0(86) * occ_func_0_1(40)) +
              (occ_func_0_1(5) * occ_func_0_0(83)) +
              (occ_func_0_0(82) * occ_func_0_1(27)) +
              (occ_func_0_1(11) * occ_func_0_0(86)) +
              (occ_func_0_0(79) * occ_func_0_1(20)) +
              (occ_func_0_1(4) * occ_func_0_0(80)) +
              (occ_func_0_0(85) * occ_func_0_1(37)) +
              (occ_func_0_1(8) * occ_func_0_0(82)) +
              (occ_func_0_0(83) * occ_func_0_1(34)) +
              (occ_func_0_1(3) * occ_func_0_0(79)) +
              (occ_func_0_0(86) * occ_func_0_1(42)) +
              (occ_func_0_1(6) * occ_func_0_0(84)) +
              (occ_func_0_0(81) * occ_func_0_1(23)) +
              (occ_func_0_1(7) * occ_func_0_0(81)) +
              (occ_func_0_0(84) * occ_func_0_1(38)) +
              (occ_func_0_1(6) * occ_func_0_0(80)) +
              (occ_func_0_0(85) * occ_func_0_1(39)) +
              (occ_func_0_1(10) * occ_func_0_0(83)) +
              (occ_func_0_0(82) * occ_func_0_1(33)) +
              (occ_func_0_1(1) * occ_func_0_0(80)) +
              (occ_func_0_0(85) * occ_func_0_1(32))) /
             24.0 +
         (m_occ_func_0_1[occ_f] - m_occ_func_0_1[occ_i]) *
             ((occ_func_0_0(6) * occ_func_0_0(39)) +
              (occ_func_0_0(4) * occ_func_0_0(26)) +
              (occ_func_0_0(8) * occ_func_0_0(25)) +
              (occ_func_0_0(3) * occ_func_0_0(42)) +
              (occ_func_0_0(10) * occ_func_0_0(33)) +
              (occ_func_0_0(1) * occ_func_0_0(40)) +
              (occ_func_0_0(11) * occ_func_0_0(30)) +
              (occ_func_0_0(5) * occ_func_0_0(36)) +
              (occ_func_0_0(2) * occ_func_0_0(31)) +
              (occ_func_0_0(4) * occ_func_0_0(37)) +
              (occ_func_0_0(9) * occ_func_0_0(35)) +
              (occ_func_0_0(1) * occ_func_0_0(32)) +
              (occ_func_0_0(11) * occ_func_0_0(20)) +
              (occ_func_0_0(12) * occ_func_0_0(21)) +
              (occ_func_0_0(8) * occ_func_0_0(34)) +
              (occ_func_0_0(2) * occ_func_0_0(41)) +
              (occ_func_0_0(9) * occ_func_0_0(24)) +
              (occ_func_0_0(5) * occ_func_0_0(27)) +
              (occ_func_0_0(10) * occ_func_0_0(19)) +
              (occ_func_0_0(7) * occ_func_0_0(38)) +
              (occ_func_0_0(6) * occ_func_0_0(23)) +
              (occ_func_0_0(7) * occ_func_0_0(22)) +
              (occ_func_0_0(3) * occ_func_0_0(28)) +
              (occ_func_0_0(12) * occ_func_0_0(29))) /
             24.0;
}
double test_Clexulator::delta_site_eval_at_0_bfunc_3_3_3(int occ_i,
                                                         int occ_f) const {
  return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) *
             ((occ_func_0_1(80) * occ_func_0_1(22)) +
              (occ_func_0_1(83) * occ_func_0_1(35)) +
              (occ_func_0_1(84) * occ_func_0_1(36)) +
              (occ_func_0_1(79) * occ_func_0_1(19)) +
              (occ_func_0_1(83) * occ_func_0_1(28)) +
              (occ_func_0_1(79) * occ_func_0_1(21)) +
              (occ_func_0_1(84) * occ_func_0_1(31)) +
              (occ_func_0_1(81) * occ_func_0_1(25)) +
              (occ_func_0_1(81) * occ_func_0_1(30)) +
              (occ_func_0_1(80) * occ_func_0_1(24)) +
              (occ_func_0_1(82) * occ_func_0_1(26)) +
              (occ_func_0_1(80) * occ_func_0_1(29)) +
              (occ_func_0_1(86) * occ_func_0_1(41)) +
              (occ_func_0_1(86) * occ_func_0_1(40)) +
              (occ_func_0_1(82) * occ_func_0_1(27)) +
              (occ_func_0_1(79) * occ_func_0_1(20)) +
              (occ_func_0_1(85) * occ_func_0_1(37)) +
              (occ_func_0_1(83) * occ_func_0_1(34)) +
              (occ_func_0_1(86) * occ_func_0_1(42)) +
              (occ_func_0_1(81) * occ_func_0_1(23)) +
              (occ_func_0_1(84) * occ_func_0_1(38)) +
              (occ_func_0_1(85) * occ_func_0_1(39)) +
              (occ_func_0_1(82) * occ_func_0_1(33)) +
              (occ_func_0_1(85) * occ_func_0_1(32))) /
             24.0 +
         (m_occ_func_0_1[occ_f] - m_occ_func_0_1[occ_i]) *
             ((occ_func_0_1(7) * occ_func_0_0(85)) +
              (occ_func_0_1(6) * occ_func_0_0(39)) +
              (occ_func_0_1(9) * occ_func_0_0(82)) +
              (occ_func_0_1(4) * occ_func_0_0(26)) +
              (occ_func_0_1(5) * occ_func_0_0(81)) +
              (occ_func_0_1(8) * occ_func_0_0(25)) +
              (occ_func_0_1(10) * occ_func_0_0(86)) +
              (occ_func_0_1(3) * occ_func_0_0(42)) +
              (occ_func_0_1(3) * occ_func_0_0(82)) +
              (occ_func_0_1(10) * occ_func_0_0(33)) +
              (occ_func_0_1(12) * occ_func_0_0(86)) +
              (occ_func_0_1(1) * occ_func_0_0(40)) +
              (occ_func_0_1(2) * occ_func_0_0(81)) +
              (occ_func_0_1(11) * occ_func_0_0(30)) +
              (occ_func_0_1(8) * occ_func_0_0(84)) +
              (occ_func_0_1(5) * occ_func_0_0(36)) +
              (occ_func_0_1(11) * occ_func_0_0(84)) +
              (occ_func_0_1(2) * occ_func_0_0(31)) +
              (occ_func_0_1(9) * occ_func_0_0(85)) +
              (occ_func_0_1(4) * occ_func_0_0(37)) +
              (occ_func_0_1(4) * occ_func_0_0(83)) +
              (occ_func_0_1(9) * occ_func_0_0(35)) +
              (occ_func_0_1(12) * occ_func_0_0(85)) +
              (occ_func_0_1(1) * occ_func_0_0(32)) +
              (occ_func_0_1(2) * occ_func_0_0(79)) +
              (occ_func_0_1(11) * occ_func_0_0(20)) +
              (occ_func_0_1(1) * occ_func_0_0(79)) +
              (occ_func_0_1(12) * occ_func_0_0(21)) +
              (occ_func_0_1(5) * occ_func_0_0(83)) +
              (occ_func_0_1(8) * occ_func_0_0(34)) +
              (occ_func_0_1(11) * occ_func_0_0(86)) +
              (occ_func_0_1(2) * occ_func_0_0(41)) +
              (occ_func_0_1(4) * occ_func_0_0(80)) +
              (occ_func_0_1(9) * occ_func_0_0(24)) +
              (occ_func_0_1(8) * occ_func_0_0(82)) +
              (occ_func_0_1(5) * occ_func_0_0(27)) +
              (occ_func_0_1(3) * occ_func_0_0(79)) +
              (occ_func_0_1(10) * occ_func_0_0(19)) +
              (occ_func_0_1(6) * occ_func_0_0(84)) +
              (occ_func_0_1(7) * occ_func_0_0(38)) +
              (occ_func_0_1(7) * occ_func_0_0(81)) +
              (occ_func_0_1(6) * occ_func_0_0(23)) +
              (occ_func_0_1(6) * occ_func_0_0(80)) +
              (occ_func_0_1(7) * occ_func_0_0(22)) +
              (occ_func_0_1(10) * occ_func_0_0(83)) +
              (occ_func_0_1(3) * occ_func_0_0(28)) +
              (occ_func_0_1(1) * occ_func_0_0(80)) +
              (occ_func_0_1(12) * occ_func_0_0(29))) /
             24.0;
}
double test_Clexulator::delta_site_eval_at_0_bfunc_3_3_4(int occ_i,
                                                         int occ_f) const {
  return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) *
             ((occ_func_0_0(7) * occ_func_0_1(85)) +
              (occ_func_0_0(6) * occ_func_0_1(39)) +
              (occ_func_0_0(9) * occ_func_0_1(82)) +
              (occ_func_0_0(4) * occ_func_0_1(26)) +
              (occ_func_0_0(5) * occ_func_0_1(81)) +
              (occ_func_0_0(8) * occ_func_0_1(25)) +
              (occ_func_0_0(10) * occ_func_0_1(86)) +
              (occ_func_0_0(3) * occ_func_0_1(42)) +
              (occ_func_0_0(3) * occ_func_0_1(82)) +
              (occ_func_0_0(10) * occ_func_0_1(33)) +
              (occ_func_0_0(12) * occ_func_0_1(86)) +
              (occ_func_0_0(1) * occ_func_0_1(40)) +
              (occ_func_0_0(2) * occ_func_0_1(81)) +
              (occ_func_0_0(11) * occ_func_0_1(30)) +
              (occ_func_0_0(8) * occ_func_0_1(84)) +
              (occ_func_0_0(5) * occ_func_0_1(36)) +
              (occ_func_0_0(11) * occ_func_0_1(84)) +
              (occ_func_0_0(2) * occ_func_0_1(31)) +
              (occ_func_0_0(9) * occ_func_0_1(85)) +
              (occ_func_0_0(4) * occ_func_0_1(37)) +
              (occ_func_0_0(4) * occ_func_0_1(83)) +
              (occ_func_0_0(9) * occ_func_0_1(35)) +
              (occ_func_0_0(12) * occ_func_0_1(85)) +
              (occ_func_0_0(1) * occ_func_0_1(32)) +
              (occ_func_0_0(2) * occ_func_0_1(79)) +
              (occ_func_0_0(11) * occ_func_0_1(20)) +
              (occ_func_0_0(1) * occ_func_0_1(79)) +
              (occ_func_0_0(12) * occ_func_0_1(21)) +
              (occ_func_0_0(5) * occ_func_0_1(83)) +
              (occ_func_0_0(8) * occ_func_0_1(34)) +
              (occ_func_0_0(11) * occ_func_0_1(86)) +
              (occ_func_0_0(2) * occ_func_0_1(41)) +
              (occ_func_0_0(4) * occ_func_0_1(80)) +
              (occ_func_0_0(9) * occ_func_0_1(24)) +
              (occ_func_0_0(8) * occ_func_0_1(82)) +
              (occ_func_0_0(5) * occ_func_0_1(27)) +
              (occ_func_0_0(3) * occ_func_0_1(79)) +
              (occ_func_0_0(10) * occ_func_0_1(19)) +
              (occ_func_0_0(6) * occ_func_0_1(84)) +
              (occ_func_0_0(7) * occ_func_0_1(38)) +
              (occ_func_0_0(7) * occ_func_0_1(81)) +
              (occ_func_0_0(6) * occ_func_0_1(23)) +
              (occ_func_0_0(6) * occ_func_0_1(80)) +
              (occ_func_0_0(7) * occ_func_0_1(22)) +
              (occ_func_0_0(10) * occ_func_0_1(83)) +
              (occ_func_0_0(3) * occ_func_0_1(28)) +
              (occ_func_0_0(1) * occ_func_0_1(80)) +
              (occ_func_0_0(12) * occ_func_0_1(29))) /
             24.0 +
         (m_occ_func_0_1[occ_f] - m_occ_func_0_1[occ_i]) *
             ((occ_func_0_0(80) * occ_func_0_0(22)) +
              (occ_func_0_0(83) * occ_func_0_0(35)) +
              (occ_func_0_0(84) * occ_func_0_0(36)) +
              (occ_func_0_0(79) * occ_func_0_0(19)) +
              (occ_func_0_0(83) * occ_func_0_0(28)) +
              (occ_func_0_0(79) * occ_func_0_0(21)) +
              (occ_func_0_0(84) * occ_func_0_0(31)) +
              (occ_func_0_0(81) * occ_func_0_0(25)) +
              (occ_func_0_0(81) * occ_func_0_0(30)) +
              (occ_func_0_0(80) * occ_func_0_0(24)) +
              (occ_func_0_0(82) * occ_func_0_0(26)) +
              (occ_func_0_0(80) * occ_func_0_0(29)) +
              (occ_func_0_0(86) * occ_func_0_0(41)) +
              (occ_func_0_0(86) * occ_func_0_0(40)) +
              (occ_func_0_0(82) * occ_func_0_0(27)) +
              (occ_func_0_0(79) * occ_func_0_0(20)) +
              (occ_func_0_0(85) * occ_func_0_0(37)) +
              (occ_func_0_0(83) * occ_func_0_0(34)) +
              (occ_func_0_0(86) * occ_func_0_0(42)) +
              (occ_func_0_0(81) * occ_func_0_0(23)) +
              (occ_func_0_0(84) * occ_func_0_0(38)) +
              (occ_func_0_0(85) * occ_func_0_0(39)) +
              (occ_func_0_0(82) * occ_func_0_0(33)) +
              (occ_func_0_0(85) * occ_func_0_0(32))) /
             24.0;
}
double test_Clexulator::delta_site_eval_at_0_bfunc_3_3_5(int occ_i,
                                                         int occ_f) const {
  return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) *
             ((occ_func_0_1(6) * occ_func_0_1(39)) +
              (occ_func_0_1(4) * occ_func_0_1(26)) +
              (occ_func_0_1(8) * occ_func_0_1(25)) +
              (occ_func_0_1(3) * occ_func_0_1(42)) +
              (occ_func_0_1(10) * occ_func_0_1(33)) +
              (occ_func_0_1(1) * occ_func_0_1(40)) +
              (occ_func_0_1(11) * occ_func_0_1(30)) +
              (occ_func_0_1(5) * occ_func_0_1(36)) +
              (occ_func_0_1(2) * occ_func_0_1(31)) +
              (occ_func_0_1(4) * occ_func_0_1(37)) +
              (occ_func_0_1(9) * occ_func_0_1(35)) +
              (occ_func_0_1(1) * occ_func_0_1(32)) +
              (occ_func_0_1(11) * occ_func_0_1(20)) +
              (occ_func_0_1(12) * occ_func_0_1(21)) +
              (occ_func_0_1(8) * occ_func_0_1(34)) +
              (occ_func_0_1(2) * occ_func_0_1(41)) +
              (occ_func_0_1(9) * occ_func_0_1(24)) +
              (occ_func_0_1(5) * occ_func_0_1(27)) +
              (occ_func_0_1(10) * occ_func_0_1(19)) +
              (occ_func_0_1(7) * occ_func_0_1(38)) +
              (occ_func_0_1(6) * occ_func_0_1(23)) +
              (occ_func_0_1(7) * occ_func_0_1(22)) +
              (occ_func_0_1(3) * occ_func_0_1(28)) +
              (occ_func_0_1(12) * occ_func_0_1(29))) /
             24.0 +
         (m_occ_func_0_1[occ_f] - m_occ_func_0_1[occ_i]) *
             ((occ_func_0_0(7) * occ_func_0_1(85)) +
              (occ_func_0_1(80) * occ_func_0_0(22)) +
              (occ_func_0_0(9) * occ_func_0_1(82)) +
              (occ_func_0_1(83) * occ_func_0_0(35)) +
              (occ_func_0_0(5) * occ_func_0_1(81)) +
              (occ_func_0_1(84) * occ_func_0_0(36)) +
              (occ_func_0_0(10) * occ_func_0_1(86)) +
              (occ_func_0_1(79) * occ_func_0_0(19)) +
              (occ_func_0_0(3) * occ_func_0_1(82)) +
              (occ_func_0_1(83) * occ_func_0_0(28)) +
              (occ_func_0_0(12) * occ_func_0_1(86)) +
              (occ_func_0_1(79) * occ_func_0_0(21)) +
              (occ_func_0_0(2) * occ_func_0_1(81)) +
              (occ_func_0_1(84) * occ_func_0_0(31)) +
              (occ_func_0_0(8) * occ_func_0_1(84)) +
              (occ_func_0_1(81) * occ_func_0_0(25)) +
              (occ_func_0_0(11) * occ_func_0_1(84)) +
              (occ_func_0_1(81) * occ_func_0_0(30)) +
              (occ_func_0_0(9) * occ_func_0_1(85)) +
              (occ_func_0_1(80) * occ_func_0_0(24)) +
              (occ_func_0_0(4) * occ_func_0_1(83)) +
              (occ_func_0_1(82) * occ_func_0_0(26)) +
              (occ_func_0_0(12) * occ_func_0_1(85)) +
              (occ_func_0_1(80) * occ_func_0_0(29)) +
              (occ_func_0_0(2) * occ_func_0_1(79)) +
              (occ_func_0_1(86) * occ_func_0_0(41)) +
              (occ_func_0_0(1) * occ_func_0_1(79)) +
              (occ_func_0_1(86) * occ_func_0_0(40)) +
              (occ_func_0_0(5) * occ_func_0_1(83)) +
              (occ_func_0_1(82) * occ_func_0_0(27)) +
              (occ_func_0_0(11) * occ_func_0_1(86)) +
              (occ_func_0_1(79) * occ_func_0_0(20)) +
              (occ_func_0_0(4) * occ_func_0_1(80)) +
              (occ_func_0_1(85) * occ_func_0_0(37)) +
              (occ_func_0_0(8) * occ_func_0_1(82)) +
              (occ_func_0_1(83) * occ_func_0_0(34)) +
              (occ_func_0_0(3) * occ_func_0_1(79)) +
              (occ_func_0_1(86) * occ_func_0_0(42)) +
              (occ_func_0_0(6) * occ_func_0_1(84)) +
              (occ_func_0_1(81) * occ_func_0_0(23)) +
              (occ_func_0_0(7) * occ_func_0_1(81)) +
              (occ_func_0_1(84) * occ_func_0_0(38)) +
              (occ_func_0_0(6) * occ_func_0_1(80)) +
              (occ_func_0_1(85) * occ_func_0_0(39)) +
              (occ_func_0_0(10) * occ_func_0_1(83)) +
              (occ_func_0_1(82) * occ_func_0_0(33)) +
              (occ_func_0_0(1) * occ_func_0_1(80)) +
              (occ_func_0_1(85) * occ_func_0_0(32))) /
             24.0;
}
double test_Clexulator::delta_site_eval_at_0_bfunc_3_3_6(int occ_i,
                                                         int occ_f) const {
  return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) *
             ((occ_func_0_1(7) * occ_func_0_1(85)) +
              (occ_func_0_1(9) * occ_func_0_1(82)) +
              (occ_func_0_1(5) * occ_func_0_1(81)) +
              (occ_func_0_1(10) * occ_func_0_1(86)) +
              (occ_func_0_1(3) * occ_func_0_1(82)) +
              (occ_func_0_1(12) * occ_func_0_1(86)) +
              (occ_func_0_1(2) * occ_func_0_1(81)) +
              (occ_func_0_1(8) * occ_func_0_1(84)) +
              (occ_func_0_1(11) * occ_func_0_1(84)) +
              (occ_func_0_1(9) * occ_func_0_1(85)) +
              (occ_func_0_1(4) * occ_func_0_1(83)) +
              (occ_func_0_1(12) * occ_func_0_1(85)) +
              (occ_func_0_1(2) * occ_func_0_1(79)) +
              (occ_func_0_1(1) * occ_func_0_1(79)) +
              (occ_func_0_1(5) * occ_func_0_1(83)) +
              (occ_func_0_1(11) * occ_func_0_1(86)) +
              (occ_func_0_1(4) * occ_func_0_1(80)) +
              (occ_func_0_1(8) * occ_func_0_1(82)) +
              (occ_func_0_1(3) * occ_func_0_1(79)) +
              (occ_func_0_1(6) * occ_func_0_1(84)) +
              (occ_func_0_1(7) * occ_func_0_1(81)) +
              (occ_func_0_1(6) * occ_func_0_1(80)) +
              (occ_func_0_1(10) * occ_func_0_1(83)) +
              (occ_func_0_1(1) * occ_func_0_1(80))) /
             24.0 +
         (m_occ_func_0_1[occ_f] - m_occ_func_0_1[occ_i]) *
             ((occ_func_0_0(6) * occ_func_0_1(39)) +
              (occ_func_0_0(80) * occ_func_0_1(22)) +
              (occ_func_0_0(4) * occ_func_0_1(26)) +
              (occ_func_0_0(83) * occ_func_0_1(35)) +
              (occ_func_0_0(8) * occ_func_0_1(25)) +
              (occ_func_0_0(84) * occ_func_0_1(36)) +
              (occ_func_0_0(3) * occ_func_0_1(42)) +
              (occ_func_0_0(79) * occ_func_0_1(19)) +
              (occ_func_0_0(10) * occ_func_0_1(33)) +
              (occ_func_0_0(83) * occ_func_0_1(28)) +
              (occ_func_0_0(1) * occ_func_0_1(40)) +
              (occ_func_0_0(79) * occ_func_0_1(21)) +
              (occ_func_0_0(11) * occ_func_0_1(30)) +
              (occ_func_0_0(84) * occ_func_0_1(31)) +
              (occ_func_0_0(5) * occ_func_0_1(36)) +
              (occ_func_0_0(81) * occ_func_0_1(25)) +
              (occ_func_0_0(2) * occ_func_0_1(31)) +
              (occ_func_0_0(81) * occ_func_0_1(30)) +
              (occ_func_0_0(4) * occ_func_0_1(37)) +
              (occ_func_0_0(80) * occ_func_0_1(24)) +
              (occ_func_0_0(9) * occ_func_0_1(35)) +
              (occ_func_0_0(82) * occ_func_0_1(26)) +
              (occ_func_0_0(1) * occ_func_0_1(32)) +
              (occ_func_0_0(80) * occ_func_0_1(29)) +
              (occ_func_0_0(11) * occ_func_0_1(20)) +
              (occ_func_0_0(86) * occ_func_0_1(41)) +
              (occ_func_0_0(12) * occ_func_0_1(21)) +
              (occ_func_0_0(86) * occ_func_0_1(40)) +
              (occ_func_0_0(8) * occ_func_0_1(34)) +
              (occ_func_0_0(82) * occ_func_0_1(27)) +
              (occ_func_0_0(2) * occ_func_0_1(41)) +
              (occ_func_0_0(79) * occ_func_0_1(20)) +
              (occ_func_0_0(9) * occ_func_0_1(24)) +
              (occ_func_0_0(85) * occ_func_0_1(37)) +
              (occ_func_0_0(5) * occ_func_0_1(27)) +
              (occ_func_0_0(83) * occ_func_0_1(34)) +
              (occ_func_0_0(10) * occ_func_0_1(19)) +
              (occ_func_0_0(86) * occ_func_0_1(42)) +
              (occ_func_0_0(7) * occ_func_0_1(38)) +
              (occ_func_0_0(81) * occ_func_0_1(23)) +
              (occ_func_0_0(6) * occ_func_0_1(23)) +
              (occ_func_0_0(84) * occ_func_0_1(38)) +
              (occ_func_0_0(7) * occ_func_0_1(22)) +
              (occ_func_0_0(85) * occ_func_0_1(39)) +
              (occ_func_0_0(3) * occ_func_0_1(28)) +
              (occ_func_0_0(82) * occ_func_0_1(33)) +
              (occ_func_0_0(12) * occ_func_0_1(29)) +
              (occ_func_0_0(85) * occ_func_0_1(32))) /
             24.0;
}
double test_Clexulator::delta_site_eval_at_0_bfunc_3_3_7(int occ_i,
                                                         int occ_f) const {
  return (m_occ_func_0_1[occ_f] - m_occ_func_0_1[occ_i]) *
         ((occ_func_0_1(7) * occ_func_0_1(85)) +
          (occ_func_0_1(6) * occ_func_0_1(39)) +
          (occ_func_0_1(80) * occ_func_0_1(22)) +
          (occ_func_0_1(9) * occ_func_0_1(82)) +
          (occ_func_0_1(4) * occ_func_0_1(26)) +
          (occ_func_0_1(83) * occ_func_0_1(35)) +
          (occ_func_0_1(5) * occ_func_0_1(81)) +
          (occ_func_0_1(8) * occ_func_0_1(25)) +
          (occ_func_0_1(84) * occ_func_0_1(36)) +
          (occ_func_0_1(10) * occ_func_0_1(86)) +
          (occ_func_0_1(3) * occ_func_0_1(42)) +
          (occ_func_0_1(79) * occ_func_0_1(19)) +
          (occ_func_0_1(3) * occ_func_0_1(82)) +
          (occ_func_0_1(10) * occ_func_0_1(33)) +
          (occ_func_0_1(83) * occ_func_0_1(28)) +
          (occ_func_0_1(12) * occ_func_0_1(86)) +
          (occ_func_0_1(1) * occ_func_0_1(40)) +
          (occ_func_0_1(79) * occ_func_0_1(21)) +
          (occ_func_0_1(2) * occ_func_0_1(81)) +
          (occ_func_0_1(11) * occ_func_0_1(30)) +
          (occ_func_0_1(84) * occ_func_0_1(31)) +
          (occ_func_0_1(8) * occ_func_0_1(84)) +
          (occ_func_0_1(5) * occ_func_0_1(36)) +
          (occ_func_0_1(81) * occ_func_0_1(25)) +
          (occ_func_0_1(11) * occ_func_0_1(84)) +
          (occ_func_0_1(2) * occ_func_0_1(31)) +
          (occ_func_0_1(81) * occ_func_0_1(30)) +
          (occ_func_0_1(9) * occ_func_0_1(85)) +
          (occ_func_0_1(4) * occ_func_0_1(37)) +
          (occ_func_0_1(80) * occ_func_0_1(24)) +
          (occ_func_0_1(4) * occ_func_0_1(83)) +
          (occ_func_0_1(9) * occ_func_0_1(35)) +
          (occ_func_0_1(82) * occ_func_0_1(26)) +
          (occ_func_0_1(12) * occ_func_0_1(85)) +
          (occ_func_0_1(1) * occ_func_0_1(32)) +
          (occ_func_0_1(80) * occ_func_0_1(29)) +
          (occ_func_0_1(2) * occ_func_0_1(79)) +
          (occ_func_0_1(11) * occ_func_0_1(20)) +
          (occ_func_0_1(86) * occ_func_0_1(41)) +
          (occ_func_0_1(1) * occ_func_0_1(79)) +
          (occ_func_0_1(12) * occ_func_0_1(21)) +
          (occ_func_0_1(86) * occ_func_0_1(40)) +
          (occ_func_0_1(5) * occ_func_0_1(83)) +
          (occ_func_0_1(8) * occ_func_0_1(34)) +
          (occ_func_0_1(82) * occ_func_0_1(27)) +
          (occ_func_0_1(11) * occ_func_0_1(86)) +
          (occ_func_0_1(2) * occ_func_0_1(41)) +
          (occ_func_0_1(79) * occ_func_0_1(20)) +
          (occ_func_0_1(4) * occ_func_0_1(80)) +
          (occ_func_0_1(9) * occ_func_0_1(24)) +
          (occ_func_0_1(85) * occ_func_0_1(37)) +
          (occ_func_0_1(8) * occ_func_0_1(82)) +
          (occ_func_0_1(5) * occ_func_0_1(27)) +
          (occ_func_0_1(83) * occ_func_0_1(34)) +
          (occ_func_0_1(3) * occ_func_0_1(79)) +
          (occ_func_0_1(10) * occ_func_0_1(19)) +
          (occ_func_0_1(86) * occ_func_0_1(42)) +
          (occ_func_0_1(6) * occ_func_0_1(84)) +
          (occ_func_0_1(7) * occ_func_0_1(38)) +
          (occ_func_0_1(81) * occ_func_0_1(23)) +
          (occ_func_0_1(7) * occ_func_0_1(81)) +
          (occ_func_0_1(6) * occ_func_0_1(23)) +
          (occ_func_0_1(84) * occ_func_0_1(38)) +
          (occ_func_0_1(6) * occ_func_0_1(80)) +
          (occ_func_0_1(7) * occ_func_0_1(22)) +
          (occ_func_0_1(85) * occ_func_0_1(39)) +
          (occ_func_0_1(10) * occ_func_0_1(83)) +
          (occ_func_0_1(3) * occ_func_0_1(28)) +
          (occ_func_0_1(82) * occ_func_0_1(33)) +
          (occ_func_0_1(1) * occ_func_0_1(80)) +
          (occ_func_0_1(12) * occ_func_0_1(29)) +
          (occ_func_0_1(85) * occ_func_0_1(32))) /
         24.0;
}

/**** Basis functions for orbit 3, 4****
#Points: 3
MaxLength: 8.4852814  MinLength: 2.8284271
 0.0000000   0.0000000   0.0000000 A B C
 2.0000000   0.0000000   0.0000000 A B C
 3.0000000   0.0000000   0.0000000 A B C
****/
double test_Clexulator::eval_bfunc_3_4_0() const {
  return ((occ_func_0_0(0) * occ_func_0_0(54) * occ_func_0_0(175)) +
          (occ_func_0_0(0) * occ_func_0_0(50) * occ_func_0_0(160)) +
          (occ_func_0_0(0) * occ_func_0_0(49) * occ_func_0_0(159)) +
          (occ_func_0_0(0) * occ_func_0_0(53) * occ_func_0_0(174)) +
          (occ_func_0_0(0) * occ_func_0_0(51) * occ_func_0_0(161)) +
          (occ_func_0_0(0) * occ_func_0_0(52) * occ_func_0_0(171)) +
          (occ_func_0_0(0) * occ_func_0_0(47) * occ_func_0_0(157)) +
          (occ_func_0_0(0) * occ_func_0_0(48) * occ_func_0_0(158)) +
          (occ_func_0_0(0) * occ_func_0_0(45) * occ_func_0_0(146)) +
          (occ_func_0_0(0) * occ_func_0_0(44) * occ_func_0_0(143)) +
          (occ_func_0_0(0) * occ_func_0_0(46) * occ_func_0_0(156)) +
          (occ_func_0_0(0) * occ_func_0_0(43) * occ_func_0_0(142))) /
         12.0;
}
double test_Clexulator::eval_bfunc_3_4_1() const {
  return ((occ_func_0_1(0) * occ_func_0_0(54) * occ_func_0_0(175)) +
          (occ_func_0_1(0) * occ_func_0_0(50) * occ_func_0_0(160)) +
          (occ_func_0_1(0) * occ_func_0_0(49) * occ_func_0_0(159)) +
          (occ_func_0_1(0) * occ_func_0_0(53) * occ_func_0_0(174)) +
          (occ_func_0_1(0) * occ_func_0_0(51) * occ_func_0_0(161)) +
          (occ_func_0_1(0) * occ_func_0_0(52) * occ_func_0_0(171)) +
          (occ_func_0_1(0) * occ_func_0_0(47) * occ_func_0_0(157)) +
          (occ_func_0_1(0) * occ_func_0_0(48) * occ_func_0_0(158)) +
          (occ_func_0_1(0) * occ_func_0_0(45) * occ_func_0_0(146)) +
          (occ_func_0_1(0) * occ_func_0_0(44) * occ_func_0_0(143)) +
          (occ_func_0_1(0) * occ_func_0_0(46) * occ_func_0_0(156)) +
          (occ_func_0_1(0) * occ_func_0_0(43) * occ_func_0_0(142))) /
         12.0;
}
double test_Clexulator::eval_bfunc_3_4_2() const {
  return ((occ_func_0_0(0) * occ_func_0_1(54) * occ_func_0_0(175)) +
          (occ_func_0_0(0) * occ_func_0_1(50) * occ_func_0_0(160)) +
          (occ_func_0_0(0) * occ_func_0_1(49) * occ_func_0_0(159)) +
          (occ_func_0_0(0) * occ_func_0_1(53) * occ_func_0_0(174)) +
          (occ_func_0_0(0) * occ_func_0_1(51) * occ_func_0_0(161)) +
          (occ_func_0_0(0) * occ_func_0_1(52) * occ_func_0_0(171)) +
          (occ_func_0_0(0) * occ_func_0_1(47) * occ_func_0_0(157)) +
          (occ_func_0_0(0) * occ_func_0_1(48) * occ_func_0_0(158)) +
          (occ_func_0_0(0) * occ_func_0_1(45) * occ_func_0_0(146)) +
          (occ_func_0_0(0) * occ_func_0_1(44) * occ_func_0_0(143)) +
          (occ_func_0_0(0) * occ_func_0_1(46) * occ_func_0_0(156)) +
          (occ_func_0_0(0) * occ_func_0_1(43) * occ_func_0_0(142))) /
         12.0;
}
double test_Clexulator::eval_bfunc_3_4_3() const {
  return ((occ_func_0_1(0) * occ_func_0_1(54) * occ_func_0_0(175)) +
          (occ_func_0_1(0) * occ_func_0_1(50) * occ_func_0_0(160)) +
          (occ_func_0_1(0) * occ_func_0_1(49) * occ_func_0_0(159)) +
          (occ_func_0_1(0) * occ_func_0_1(53) * occ_func_0_0(174)) +
          (occ_func_0_1(0) * occ_func_0_1(51) * occ_func_0_0(161)) +
          (occ_func_0_1(0) * occ_func_0_1(52) * occ_func_0_0(171)) +
          (occ_func_0_1(0) * occ_func_0_1(47) * occ_func_0_0(157)) +
          (occ_func_0_1(0) * occ_func_0_1(48) * occ_func_0_0(158)) +
          (occ_func_0_1(0) * occ_func_0_1(45) * occ_func_0_0(146)) +
          (occ_func_0_1(0) * occ_func_0_1(44) * occ_func_0_0(143)) +
          (occ_func_0_1(0) * occ_func_0_1(46) * occ_func_0_0(156)) +
          (occ_func_0_1(0) * occ_func_0_1(43) * occ_func_0_0(142))) /
         12.0;
}
double test_Clexulator::eval_bfunc_3_4_4() const {
  return ((occ_func_0_0(0) * occ_func_0_0(54) * occ_func_0_1(175)) +
          (occ_func_0_0(0) * occ_func_0_0(50) * occ_func_0_1(160)) +
          (occ_func_0_0(0) * occ_func_0_0(49) * occ_func_0_1(159)) +
          (occ_func_0_0(0) * occ_func_0_0(53) * occ_func_0_1(174)) +
          (occ_func_0_0(0) * occ_func_0_0(51) * occ_func_0_1(161)) +
          (occ_func_0_0(0) * occ_func_0_0(52) * occ_func_0_1(171)) +
          (occ_func_0_0(0) * occ_func_0_0(47) * occ_func_0_1(157)) +
          (occ_func_0_0(0) * occ_func_0_0(48) * occ_func_0_1(158)) +
          (occ_func_0_0(0) * occ_func_0_0(45) * occ_func_0_1(146)) +
          (occ_func_0_0(0) * occ_func_0_0(44) * occ_func_0_1(143)) +
          (occ_func_0_0(0) * occ_func_0_0(46) * occ_func_0_1(156)) +
          (occ_func_0_0(0) * occ_func_0_0(43) * occ_func_0_1(142))) /
         12.0;
}
double test_Clexulator::eval_bfunc_3_4_5() const {
  return ((occ_func_0_1(0) * occ_func_0_0(54) * occ_func_0_1(175)) +
          (occ_func_0_1(0) * occ_func_0_0(50) * occ_func_0_1(160)) +
          (occ_func_0_1(0) * occ_func_0_0(49) * occ_func_0_1(159)) +
          (occ_func_0_1(0) * occ_func_0_0(53) * occ_func_0_1(174)) +
          (occ_func_0_1(0) * occ_func_0_0(51) * occ_func_0_1(161)) +
          (occ_func_0_1(0) * occ_func_0_0(52) * occ_func_0_1(171)) +
          (occ_func_0_1(0) * occ_func_0_0(47) * occ_func_0_1(157)) +
          (occ_func_0_1(0) * occ_func_0_0(48) * occ_func_0_1(158)) +
          (occ_func_0_1(0) * occ_func_0_0(45) * occ_func_0_1(146)) +
          (occ_func_0_1(0) * occ_func_0_0(44) * occ_func_0_1(143)) +
          (occ_func_0_1(0) * occ_func_0_0(46) * occ_func_0_1(156)) +
          (occ_func_0_1(0) * occ_func_0_0(43) * occ_func_0_1(142))) /
         12.0;
}
double test_Clexulator::eval_bfunc_3_4_6() const {
  return ((occ_func_0_0(0) * occ_func_0_1(54) * occ_func_0_1(175)) +
          (occ_func_0_0(0) * occ_func_0_1(50) * occ_func_0_1(160)) +
          (occ_func_0_0(0) * occ_func_0_1(49) * occ_func_0_1(159)) +
          (occ_func_0_0(0) * occ_func_0_1(53) * occ_func_0_1(174)) +
          (occ_func_0_0(0) * occ_func_0_1(51) * occ_func_0_1(161)) +
          (occ_func_0_0(0) * occ_func_0_1(52) * occ_func_0_1(171)) +
          (occ_func_0_0(0) * occ_func_0_1(47) * occ_func_0_1(157)) +
          (occ_func_0_0(0) * occ_func_0_1(48) * occ_func_0_1(158)) +
          (occ_func_0_0(0) * occ_func_0_1(45) * occ_func_0_1(146)) +
          (occ_func_0_0(0) * occ_func_0_1(44) * occ_func_0_1(143)) +
          (occ_func_0_0(0) * occ_func_0_1(46) * occ_func_0_1(156)) +
          (occ_func_0_0(0) * occ_func_0_1(43) * occ_func_0_1(142))) /
         12.0;
}
double test_Clexulator::eval_bfunc_3_4_7() const {
  return ((occ_func_0_1(0) * occ_func_0_1(54) * occ_func_0_1(175)) +
          (occ_func_0_1(0) * occ_func_0_1(50) * occ_func_0_1(160)) +
          (occ_func_0_1(0) * occ_func_0_1(49) * occ_func_0_1(159)) +
          (occ_func_0_1(0) * occ_func_0_1(53) * occ_func_0_1(174)) +
          (occ_func_0_1(0) * occ_func_0_1(51) * occ_func_0_1(161)) +
          (occ_func_0_1(0) * occ_func_0_1(52) * occ_func_0_1(171)) +
          (occ_func_0_1(0) * occ_func_0_1(47) * occ_func_0_1(157)) +
          (occ_func_0_1(0) * occ_func_0_1(48) * occ_func_0_1(158)) +
          (occ_func_0_1(0) * occ_func_0_1(45) * occ_func_0_1(146)) +
          (occ_func_0_1(0) * occ_func_0_1(44) * occ_func_0_1(143)) +
          (occ_func_0_1(0) * occ_func_0_1(46) * occ_func_0_1(156)) +
          (occ_func_0_1(0) * occ_func_0_1(43) * occ_func_0_1(142))) /
         12.0;
}

double test_Clexulator::site_eval_at_0_bfunc_3_4_0() const {
  return ((occ_func_0_0(0) * occ_func_0_0(54) * occ_func_0_0(175)) +
          (occ_func_0_0(43) * occ_func_0_0(0) * occ_func_0_0(12)) +
          (occ_func_0_0(142) * occ_func_0_0(1) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(50) * occ_func_0_0(160)) +
          (occ_func_0_0(47) * occ_func_0_0(0) * occ_func_0_0(8)) +
          (occ_func_0_0(157) * occ_func_0_0(5) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(49) * occ_func_0_0(159)) +
          (occ_func_0_0(48) * occ_func_0_0(0) * occ_func_0_0(7)) +
          (occ_func_0_0(158) * occ_func_0_0(6) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(53) * occ_func_0_0(174)) +
          (occ_func_0_0(44) * occ_func_0_0(0) * occ_func_0_0(11)) +
          (occ_func_0_0(143) * occ_func_0_0(2) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(51) * occ_func_0_0(161)) +
          (occ_func_0_0(46) * occ_func_0_0(0) * occ_func_0_0(9)) +
          (occ_func_0_0(156) * occ_func_0_0(4) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(52) * occ_func_0_0(171)) +
          (occ_func_0_0(45) * occ_func_0_0(0) * occ_func_0_0(10)) +
          (occ_func_0_0(146) * occ_func_0_0(3) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(47) * occ_func_0_0(157)) +
          (occ_func_0_0(50) * occ_func_0_0(0) * occ_func_0_0(5)) +
          (occ_func_0_0(160) * occ_func_0_0(8) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(48) * occ_func_0_0(158)) +
          (occ_func_0_0(49) * occ_func_0_0(0) * occ_func_0_0(6)) +
          (occ_func_0_0(159) * occ_func_0_0(7) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(45) * occ_func_0_0(146)) +
          (occ_func_0_0(52) * occ_func_0_0(0) * occ_func_0_0(3)) +
          (occ_func_0_0(171) * occ_func_0_0(10) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(44) * occ_func_0_0(143)) +
          (occ_func_0_0(53) * occ_func_0_0(0) * occ_func_0_0(2)) +
          (occ_func_0_0(174) * occ_func_0_0(11) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(46) * occ_func_0_0(156)) +
          (occ_func_0_0(51) * occ_func_0_0(0) * occ_func_0_0(4)) +
          (occ_func_0_0(161) * occ_func_0_0(9) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(43) * occ_func_0_0(142)) +
          (occ_func_0_0(54) * occ_func_0_0(0) * occ_func_0_0(1)) +
          (occ_func_0_0(175) * occ_func_0_0(12) * occ_func_0_0(0))) /
         12.0;
}
double test_Clexulator::site_eval_at_0_bfunc_3_4_1() const {
  return ((occ_func_0_1(0) * occ_func_0_0(54) * occ_func_0_0(175)) +
          (occ_func_0_1(43) * occ_func_0_0(0) * occ_func_0_0(12)) +
          (occ_func_0_1(142) * occ_func_0_0(1) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_0(50) * occ_func_0_0(160)) +
          (occ_func_0_1(47) * occ_func_0_0(0) * occ_func_0_0(8)) +
          (occ_func_0_1(157) * occ_func_0_0(5) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_0(49) * occ_func_0_0(159)) +
          (occ_func_0_1(48) * occ_func_0_0(0) * occ_func_0_0(7)) +
          (occ_func_0_1(158) * occ_func_0_0(6) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_0(53) * occ_func_0_0(174)) +
          (occ_func_0_1(44) * occ_func_0_0(0) * occ_func_0_0(11)) +
          (occ_func_0_1(143) * occ_func_0_0(2) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_0(51) * occ_func_0_0(161)) +
          (occ_func_0_1(46) * occ_func_0_0(0) * occ_func_0_0(9)) +
          (occ_func_0_1(156) * occ_func_0_0(4) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_0(52) * occ_func_0_0(171)) +
          (occ_func_0_1(45) * occ_func_0_0(0) * occ_func_0_0(10)) +
          (occ_func_0_1(146) * occ_func_0_0(3) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_0(47) * occ_func_0_0(157)) +
          (occ_func_0_1(50) * occ_func_0_0(0) * occ_func_0_0(5)) +
          (occ_func_0_1(160) * occ_func_0_0(8) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_0(48) * occ_func_0_0(158)) +
          (occ_func_0_1(49) * occ_func_0_0(0) * occ_func_0_0(6)) +
          (occ_func_0_1(159) * occ_func_0_0(7) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_0(45) * occ_func_0_0(146)) +
          (occ_func_0_1(52) * occ_func_0_0(0) * occ_func_0_0(3)) +
          (occ_func_0_1(171) * occ_func_0_0(10) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_0(44) * occ_func_0_0(143)) +
          (occ_func_0_1(53) * occ_func_0_0(0) * occ_func_0_0(2)) +
          (occ_func_0_1(174) * occ_func_0_0(11) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_0(46) * occ_func_0_0(156)) +
          (occ_func_0_1(51) * occ_func_0_0(0) * occ_func_0_0(4)) +
          (occ_func_0_1(161) * occ_func_0_0(9) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_0(43) * occ_func_0_0(142)) +
          (occ_func_0_1(54) * occ_func_0_0(0) * occ_func_0_0(1)) +
          (occ_func_0_1(175) * occ_func_0_0(12) * occ_func_0_0(0))) /
         12.0;
}
double test_Clexulator::site_eval_at_0_bfunc_3_4_2() const {
  return ((occ_func_0_0(0) * occ_func_0_1(54) * occ_func_0_0(175)) +
          (occ_func_0_0(43) * occ_func_0_1(0) * occ_func_0_0(12)) +
          (occ_func_0_0(142) * occ_func_0_1(1) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_1(50) * occ_func_0_0(160)) +
          (occ_func_0_0(47) * occ_func_0_1(0) * occ_func_0_0(8)) +
          (occ_func_0_0(157) * occ_func_0_1(5) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_1(49) * occ_func_0_0(159)) +
          (occ_func_0_0(48) * occ_func_0_1(0) * occ_func_0_0(7)) +
          (occ_func_0_0(158) * occ_func_0_1(6) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_1(53) * occ_func_0_0(174)) +
          (occ_func_0_0(44) * occ_func_0_1(0) * occ_func_0_0(11)) +
          (occ_func_0_0(143) * occ_func_0_1(2) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_1(51) * occ_func_0_0(161)) +
          (occ_func_0_0(46) * occ_func_0_1(0) * occ_func_0_0(9)) +
          (occ_func_0_0(156) * occ_func_0_1(4) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_1(52) * occ_func_0_0(171)) +
          (occ_func_0_0(45) * occ_func_0_1(0) * occ_func_0_0(10)) +
          (occ_func_0_0(146) * occ_func_0_1(3) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_1(47) * occ_func_0_0(157)) +
          (occ_func_0_0(50) * occ_func_0_1(0) * occ_func_0_0(5)) +
          (occ_func_0_0(160) * occ_func_0_1(8) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_1(48) * occ_func_0_0(158)) +
          (occ_func_0_0(49) * occ_func_0_1(0) * occ_func_0_0(6)) +
          (occ_func_0_0(159) * occ_func_0_1(7) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_1(45) * occ_func_0_0(146)) +
          (occ_func_0_0(52) * occ_func_0_1(0) * occ_func_0_0(3)) +
          (occ_func_0_0(171) * occ_func_0_1(10) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_1(44) * occ_func_0_0(143)) +
          (occ_func_0_0(53) * occ_func_0_1(0) * occ_func_0_0(2)) +
          (occ_func_0_0(174) * occ_func_0_1(11) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_1(46) * occ_func_0_0(156)) +
          (occ_func_0_0(51) * occ_func_0_1(0) * occ_func_0_0(4)) +
          (occ_func_0_0(161) * occ_func_0_1(9) * occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_1(43) * occ_func_0_0(142)) +
          (occ_func_0_0(54) * occ_func_0_1(0) * occ_func_0_0(1)) +
          (occ_func_0_0(175) * occ_func_0_1(12) * occ_func_0_0(0))) /
         12.0;
}
double test_Clexulator::site_eval_at_0_bfunc_3_4_3() const {
  return ((occ_func_0_1(0) * occ_func_0_1(54) * occ_func_0_0(175)) +
          (occ_func_0_1(43) * occ_func_0_1(0) * occ_func_0_0(12)) +
          (occ_func_0_1(142) * occ_func_0_1(1) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_1(50) * occ_func_0_0(160)) +
          (occ_func_0_1(47) * occ_func_0_1(0) * occ_func_0_0(8)) +
          (occ_func_0_1(157) * occ_func_0_1(5) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_1(49) * occ_func_0_0(159)) +
          (occ_func_0_1(48) * occ_func_0_1(0) * occ_func_0_0(7)) +
          (occ_func_0_1(158) * occ_func_0_1(6) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_1(53) * occ_func_0_0(174)) +
          (occ_func_0_1(44) * occ_func_0_1(0) * occ_func_0_0(11)) +
          (occ_func_0_1(143) * occ_func_0_1(2) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_1(51) * occ_func_0_0(161)) +
          (occ_func_0_1(46) * occ_func_0_1(0) * occ_func_0_0(9)) +
          (occ_func_0_1(156) * occ_func_0_1(4) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_1(52) * occ_func_0_0(171)) +
          (occ_func_0_1(45) * occ_func_0_1(0) * occ_func_0_0(10)) +
          (occ_func_0_1(146) * occ_func_0_1(3) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_1(47) * occ_func_0_0(157)) +
          (occ_func_0_1(50) * occ_func_0_1(0) * occ_func_0_0(5)) +
          (occ_func_0_1(160) * occ_func_0_1(8) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_1(48) * occ_func_0_0(158)) +
          (occ_func_0_1(49) * occ_func_0_1(0) * occ_func_0_0(6)) +
          (occ_func_0_1(159) * occ_func_0_1(7) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_1(45) * occ_func_0_0(146)) +
          (occ_func_0_1(52) * occ_func_0_1(0) * occ_func_0_0(3)) +
          (occ_func_0_1(171) * occ_func_0_1(10) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_1(44) * occ_func_0_0(143)) +
          (occ_func_0_1(53) * occ_func_0_1(0) * occ_func_0_0(2)) +
          (occ_func_0_1(174) * occ_func_0_1(11) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_1(46) * occ_func_0_0(156)) +
          (occ_func_0_1(51) * occ_func_0_1(0) * occ_func_0_0(4)) +
          (occ_func_0_1(161) * occ_func_0_1(9) * occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_1(43) * occ_func_0_0(142)) +
          (occ_func_0_1(54) * occ_func_0_1(0) * occ_func_0_0(1)) +
          (occ_func_0_1(175) * occ_func_0_1(12) * occ_func_0_0(0))) /
         12.0;
}
double test_Clexulator::site_eval_at_0_bfunc_3_4_4() const {
  return ((occ_func_0_0(0) * occ_func_0_0(54) * occ_func_0_1(175)) +
          (occ_func_0_0(43) * occ_func_0_0(0) * occ_func_0_1(12)) +
          (occ_func_0_0(142) * occ_func_0_0(1) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_0(50) * occ_func_0_1(160)) +
          (occ_func_0_0(47) * occ_func_0_0(0) * occ_func_0_1(8)) +
          (occ_func_0_0(157) * occ_func_0_0(5) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_0(49) * occ_func_0_1(159)) +
          (occ_func_0_0(48) * occ_func_0_0(0) * occ_func_0_1(7)) +
          (occ_func_0_0(158) * occ_func_0_0(6) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_0(53) * occ_func_0_1(174)) +
          (occ_func_0_0(44) * occ_func_0_0(0) * occ_func_0_1(11)) +
          (occ_func_0_0(143) * occ_func_0_0(2) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_0(51) * occ_func_0_1(161)) +
          (occ_func_0_0(46) * occ_func_0_0(0) * occ_func_0_1(9)) +
          (occ_func_0_0(156) * occ_func_0_0(4) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_0(52) * occ_func_0_1(171)) +
          (occ_func_0_0(45) * occ_func_0_0(0) * occ_func_0_1(10)) +
          (occ_func_0_0(146) * occ_func_0_0(3) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_0(47) * occ_func_0_1(157)) +
          (occ_func_0_0(50) * occ_func_0_0(0) * occ_func_0_1(5)) +
          (occ_func_0_0(160) * occ_func_0_0(8) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_0(48) * occ_func_0_1(158)) +
          (occ_func_0_0(49) * occ_func_0_0(0) * occ_func_0_1(6)) +
          (occ_func_0_0(159) * occ_func_0_0(7) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_0(45) * occ_func_0_1(146)) +
          (occ_func_0_0(52) * occ_func_0_0(0) * occ_func_0_1(3)) +
          (occ_func_0_0(171) * occ_func_0_0(10) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_0(44) * occ_func_0_1(143)) +
          (occ_func_0_0(53) * occ_func_0_0(0) * occ_func_0_1(2)) +
          (occ_func_0_0(174) * occ_func_0_0(11) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_0(46) * occ_func_0_1(156)) +
          (occ_func_0_0(51) * occ_func_0_0(0) * occ_func_0_1(4)) +
          (occ_func_0_0(161) * occ_func_0_0(9) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_0(43) * occ_func_0_1(142)) +
          (occ_func_0_0(54) * occ_func_0_0(0) * occ_func_0_1(1)) +
          (occ_func_0_0(175) * occ_func_0_0(12) * occ_func_0_1(0))) /
         12.0;
}
double test_Clexulator::site_eval_at_0_bfunc_3_4_5() const {
  return ((occ_func_0_1(0) * occ_func_0_0(54) * occ_func_0_1(175)) +
          (occ_func_0_1(43) * occ_func_0_0(0) * occ_func_0_1(12)) +
          (occ_func_0_1(142) * occ_func_0_0(1) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_0(50) * occ_func_0_1(160)) +
          (occ_func_0_1(47) * occ_func_0_0(0) * occ_func_0_1(8)) +
          (occ_func_0_1(157) * occ_func_0_0(5) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_0(49) * occ_func_0_1(159)) +
          (occ_func_0_1(48) * occ_func_0_0(0) * occ_func_0_1(7)) +
          (occ_func_0_1(158) * occ_func_0_0(6) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_0(53) * occ_func_0_1(174)) +
          (occ_func_0_1(44) * occ_func_0_0(0) * occ_func_0_1(11)) +
          (occ_func_0_1(143) * occ_func_0_0(2) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_0(51) * occ_func_0_1(161)) +
          (occ_func_0_1(46) * occ_func_0_0(0) * occ_func_0_1(9)) +
          (occ_func_0_1(156) * occ_func_0_0(4) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_0(52) * occ_func_0_1(171)) +
          (occ_func_0_1(45) * occ_func_0_0(0) * occ_func_0_1(10)) +
          (occ_func_0_1(146) * occ_func_0_0(3) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_0(47) * occ_func_0_1(157)) +
          (occ_func_0_1(50) * occ_func_0_0(0) * occ_func_0_1(5)) +
          (occ_func_0_1(160) * occ_func_0_0(8) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_0(48) * occ_func_0_1(158)) +
          (occ_func_0_1(49) * occ_func_0_0(0) * occ_func_0_1(6)) +
          (occ_func_0_1(159) * occ_func_0_0(7) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_0(45) * occ_func_0_1(146)) +
          (occ_func_0_1(52) * occ_func_0_0(0) * occ_func_0_1(3)) +
          (occ_func_0_1(171) * occ_func_0_0(10) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_0(44) * occ_func_0_1(143)) +
          (occ_func_0_1(53) * occ_func_0_0(0) * occ_func_0_1(2)) +
          (occ_func_0_1(174) * occ_func_0_0(11) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_0(46) * occ_func_0_1(156)) +
          (occ_func_0_1(51) * occ_func_0_0(0) * occ_func_0_1(4)) +
          (occ_func_0_1(161) * occ_func_0_0(9) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_0(43) * occ_func_0_1(142)) +
          (occ_func_0_1(54) * occ_func_0_0(0) * occ_func_0_1(1)) +
          (occ_func_0_1(175) * occ_func_0_0(12) * occ_func_0_1(0))) /
         12.0;
}
double test_Clexulator::site_eval_at_0_bfunc_3_4_6() const {
  return ((occ_func_0_0(0) * occ_func_0_1(54) * occ_func_0_1(175)) +
          (occ_func_0_0(43) * occ_func_0_1(0) * occ_func_0_1(12)) +
          (occ_func_0_0(142) * occ_func_0_1(1) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_1(50) * occ_func_0_1(160)) +
          (occ_func_0_0(47) * occ_func_0_1(0) * occ_func_0_1(8)) +
          (occ_func_0_0(157) * occ_func_0_1(5) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_1(49) * occ_func_0_1(159)) +
          (occ_func_0_0(48) * occ_func_0_1(0) * occ_func_0_1(7)) +
          (occ_func_0_0(158) * occ_func_0_1(6) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_1(53) * occ_func_0_1(174)) +
          (occ_func_0_0(44) * occ_func_0_1(0) * occ_func_0_1(11)) +
          (occ_func_0_0(143) * occ_func_0_1(2) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_1(51) * occ_func_0_1(161)) +
          (occ_func_0_0(46) * occ_func_0_1(0) * occ_func_0_1(9)) +
          (occ_func_0_0(156) * occ_func_0_1(4) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_1(52) * occ_func_0_1(171)) +
          (occ_func_0_0(45) * occ_func_0_1(0) * occ_func_0_1(10)) +
          (occ_func_0_0(146) * occ_func_0_1(3) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_1(47) * occ_func_0_1(157)) +
          (occ_func_0_0(50) * occ_func_0_1(0) * occ_func_0_1(5)) +
          (occ_func_0_0(160) * occ_func_0_1(8) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_1(48) * occ_func_0_1(158)) +
          (occ_func_0_0(49) * occ_func_0_1(0) * occ_func_0_1(6)) +
          (occ_func_0_0(159) * occ_func_0_1(7) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_1(45) * occ_func_0_1(146)) +
          (occ_func_0_0(52) * occ_func_0_1(0) * occ_func_0_1(3)) +
          (occ_func_0_0(171) * occ_func_0_1(10) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_1(44) * occ_func_0_1(143)) +
          (occ_func_0_0(53) * occ_func_0_1(0) * occ_func_0_1(2)) +
          (occ_func_0_0(174) * occ_func_0_1(11) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_1(46) * occ_func_0_1(156)) +
          (occ_func_0_0(51) * occ_func_0_1(0) * occ_func_0_1(4)) +
          (occ_func_0_0(161) * occ_func_0_1(9) * occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_1(43) * occ_func_0_1(142)) +
          (occ_func_0_0(54) * occ_func_0_1(0) * occ_func_0_1(1)) +
          (occ_func_0_0(175) * occ_func_0_1(12) * occ_func_0_1(0))) /
         12.0;
}
double test_Clexulator::site_eval_at_0_bfunc_3_4_7() const {
  return ((occ_func_0_1(0) * occ_func_0_1(54) * occ_func_0_1(175)) +
          (occ_func_0_1(43) * occ_func_0_1(0) * occ_func_0_1(12)) +
          (occ_func_0_1(142) * occ_func_0_1(1) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(50) * occ_func_0_1(160)) +
          (occ_func_0_1(47) * occ_func_0_1(0) * occ_func_0_1(8)) +
          (occ_func_0_1(157) * occ_func_0_1(5) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(49) * occ_func_0_1(159)) +
          (occ_func_0_1(48) * occ_func_0_1(0) * occ_func_0_1(7)) +
          (occ_func_0_1(158) * occ_func_0_1(6) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(53) * occ_func_0_1(174)) +
          (occ_func_0_1(44) * occ_func_0_1(0) * occ_func_0_1(11)) +
          (occ_func_0_1(143) * occ_func_0_1(2) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(51) * occ_func_0_1(161)) +
          (occ_func_0_1(46) * occ_func_0_1(0) * occ_func_0_1(9)) +
          (occ_func_0_1(156) * occ_func_0_1(4) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(52) * occ_func_0_1(171)) +
          (occ_func_0_1(45) * occ_func_0_1(0) * occ_func_0_1(10)) +
          (occ_func_0_1(146) * occ_func_0_1(3) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(47) * occ_func_0_1(157)) +
          (occ_func_0_1(50) * occ_func_0_1(0) * occ_func_0_1(5)) +
          (occ_func_0_1(160) * occ_func_0_1(8) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(48) * occ_func_0_1(158)) +
          (occ_func_0_1(49) * occ_func_0_1(0) * occ_func_0_1(6)) +
          (occ_func_0_1(159) * occ_func_0_1(7) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(45) * occ_func_0_1(146)) +
          (occ_func_0_1(52) * occ_func_0_1(0) * occ_func_0_1(3)) +
          (occ_func_0_1(171) * occ_func_0_1(10) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(44) * occ_func_0_1(143)) +
          (occ_func_0_1(53) * occ_func_0_1(0) * occ_func_0_1(2)) +
          (occ_func_0_1(174) * occ_func_0_1(11) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(46) * occ_func_0_1(156)) +
          (occ_func_0_1(51) * occ_func_0_1(0) * occ_func_0_1(4)) +
          (occ_func_0_1(161) * occ_func_0_1(9) * occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(43) * occ_func_0_1(142)) +
          (occ_func_0_1(54) * occ_func_0_1(0) * occ_func_0_1(1)) +
          (occ_func_0_1(175) * occ_func_0_1(12) * occ_func_0_1(0))) /
         12.0;
}

double test_Clexulator::delta_site_eval_at_0_bfunc_3_4_0(int occ_i,
                                                         int occ_f) const {
  return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) *
         ((occ_func_0_0(54) * occ_func_0_0(175)) +
          (occ_func_0_0(43) * occ_func_0_0(12)) +
          (occ_func_0_0(142) * occ_func_0_0(1)) +
          (occ_func_0_0(50) * occ_func_0_0(160)) +
          (occ_func_0_0(47) * occ_func_0_0(8)) +
          (occ_func_0_0(157) * occ_func_0_0(5)) +
          (occ_func_0_0(49) * occ_func_0_0(159)) +
          (occ_func_0_0(48) * occ_func_0_0(7)) +
          (occ_func_0_0(158) * occ_func_0_0(6)) +
          (occ_func_0_0(53) * occ_func_0_0(174)) +
          (occ_func_0_0(44) * occ_func_0_0(11)) +
          (occ_func_0_0(143) * occ_func_0_0(2)) +
          (occ_func_0_0(51) * occ_func_0_0(161)) +
          (occ_func_0_0(46) * occ_func_0_0(9)) +
          (occ_func_0_0(156) * occ_func_0_0(4)) +
          (occ_func_0_0(52) * occ_func_0_0(171)) +
          (occ_func_0_0(45) * occ_func_0_0(10)) +
          (occ_func_0_0(146) * occ_func_0_0(3)) +
          (occ_func_0_0(47) * occ_func_0_0(157)) +
          (occ_func_0_0(50) * occ_func_0_0(5)) +
          (occ_func_0_0(160) * occ_func_0_0(8)) +
          (occ_func_0_0(48) * occ_func_0_0(158)) +
          (occ_func_0_0(49) * occ_func_0_0(6)) +
          (occ_func_0_0(159) * occ_func_0_0(7)) +
          (occ_func_0_0(45) * occ_func_0_0(146)) +
          (occ_func_0_0(52) * occ_func_0_0(3)) +
          (occ_func_0_0(171) * occ_func_0_0(10)) +
          (occ_func_0_0(44) * occ_func_0_0(143)) +
          (occ_func_0_0(53) * occ_func_0_0(2)) +
          (occ_func_0_0(174) * occ_func_0_0(11)) +
          (occ_func_0_0(46) * occ_func_0_0(156)) +
          (occ_func_0_0(51) * occ_func_0_0(4)) +
          (occ_func_0_0(161) * occ_func_0_0(9)) +
          (occ_func_0_0(43) * occ_func_0_0(142)) +
          (occ_func_0_0(54) * occ_func_0_0(1)) +
          (occ_func_0_0(175) * occ_func_0_0(12))) /
         12.0;
}
double test_Clexulator::delta_site_eval_at_0_bfunc_3_4_1(int occ_i,
                                                         int occ_f) const {
  return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) *
             ((occ_func_0_1(43) * occ_func_0_0(12)) +
              (occ_func_0_1(142) * occ_func_0_0(1)) +
              (occ_func_0_1(47) * occ_func_0_0(8)) +
              (occ_func_0_1(157) * occ_func_0_0(5)) +
              (occ_func_0_1(48) * occ_func_0_0(7)) +
              (occ_func_0_1(158) * occ_func_0_0(6)) +
              (occ_func_0_1(44) * occ_func_0_0(11)) +
              (occ_func_0_1(143) * occ_func_0_0(2)) +
              (occ_func_0_1(46) * occ_func_0_0(9)) +
              (occ_func_0_1(156) * occ_func_0_0(4)) +
              (occ_func_0_1(45) * occ_func_0_0(10)) +
              (occ_func_0_1(146) * occ_func_0_0(3)) +
              (occ_func_0_1(50) * occ_func_0_0(5)) +
              (occ_func_0_1(160) * occ_func_0_0(8)) +
              (occ_func_0_1(49) * occ_func_0_0(6)) +
              (occ_func_0_1(159) * occ_func_0_0(7)) +
              (occ_func_0_1(52) * occ_func_0_0(3)) +
              (occ_func_0_1(171) * occ_func_0_0(10)) +
              (occ_func_0_1(53) * occ_func_0_0(2)) +
              (occ_func_0_1(174) * occ_func_0_0(11)) +
              (occ_func_0_1(51) * occ_func_0_0(4)) +
              (occ_func_0_1(161) * occ_func_0_0(9)) +
              (occ_func_0_1(54) * occ_func_0_0(1)) +
              (occ_func_0_1(175) * occ_func_0_0(12))) /
             12.0 +
         (m_occ_func_0_1[occ_f] - m_occ_func_0_1[occ_i]) *
             ((occ_func_0_0(54) * occ_func_0_0(175)) +
              (occ_func_0_0(50) * occ_func_0_0(160)) +
              (occ_func_0_0(49) * occ_func_0_0(159)) +
              (occ_func_0_0(53) * occ_func_0_0(174)) +
              (occ_func_0_0(51) * occ_func_0_0(161)) +
              (occ_func_0_0(52) * occ_func_0_0(171)) +
              (occ_func_0_0(47) * occ_func_0_0(157)) +
              (occ_func_0_0(48) * occ_func_0_0(158)) +
              (occ_func_0_0(45) * occ_func_0_0(146)) +
              (occ_func_0_0(44) * occ_func_0_0(143)) +
              (occ_func_0_0(46) * occ_func_0_0(156)) +
              (occ_func_0_0(43) * occ_func_0_0(142))) /
             12.0;
}
double test_Clexulator::delta_site_eval_at_0_bfunc_3_4_2(int occ_i,
                                                         int occ_f) const {
  return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) *
             ((occ_func_0_1(54) * occ_func_0_0(175)) +
              (occ_func_0_0(142) * occ_func_0_1(1)) +
              (occ_func_0_1(50) * occ_func_0_0(160)) +
              (occ_func_0_0(157) * occ_func_0_1(5)) +
              (occ_func_0_1(49) * occ_func_0_0(159)) +
              (occ_func_0_0(158) * occ_func_0_1(6)) +
              (occ_func_0_1(53) * occ_func_0_0(174)) +
              (occ_func_0_0(143) * occ_func_0_1(2)) +
              (occ_func_0_1(51) * occ_func_0_0(161)) +
              (occ_func_0_0(156) * occ_func_0_1(4)) +
              (occ_func_0_1(52) * occ_func_0_0(171)) +
              (occ_func_0_0(146) * occ_func_0_1(3)) +
              (occ_func_0_1(47) * occ_func_0_0(157)) +
              (occ_func_0_0(160) * occ_func_0_1(8)) +
              (occ_func_0_1(48) * occ_func_0_0(158)) +
              (occ_func_0_0(159) * occ_func_0_1(7)) +
              (occ_func_0_1(45) * occ_func_0_0(146)) +
              (occ_func_0_0(171) * occ_func_0_1(10)) +
              (occ_func_0_1(44) * occ_func_0_0(143)) +
              (occ_func_0_0(174) * occ_func_0_1(11)) +
              (occ_func_0_1(46) * occ_func_0_0(156)) +
              (occ_func_0_0(161) * occ_func_0_1(9)) +
              (occ_func_0_1(43) * occ_func_0_0(142)) +
              (occ_func_0_0(175) * occ_func_0_1(12))) /
             12.0 +
         (m_occ_func_0_1[occ_f] - m_occ_func_0_1[occ_i]) *
             ((occ_func_0_0(43) * occ_func_0_0(12)) +
              (occ_func_0_0(47) * occ_func_0_0(8)) +
              (occ_func_0_0(48) * occ_func_0_0(7)) +
              (occ_func_0_0(44) * occ_func_0_0(11)) +
              (occ_func_0_0(46) * occ_func_0_0(9)) +
              (occ_func_0_0(45) * occ_func_0_0(10)) +
              (occ_func_0_0(50) * occ_func_0_0(5)) +
              (occ_func_0_0(49) * occ_func_0_0(6)) +
              (occ_func_0_0(52) * occ_func_0_0(3)) +
              (occ_func_0_0(53) * occ_func_0_0(2)) +
              (occ_func_0_0(51) * occ_func_0_0(4)) +
              (occ_func_0_0(54) * occ_func_0_0(1))) /
             12.0;
}
double test_Clexulator::delta_site_eval_at_0_bfunc_3_4_3(int occ_i,
                                                         int occ_f) const {
  return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) *
             ((occ_func_0_1(142) * occ_func_0_1(1)) +
              (occ_func_0_1(157) * occ_func_0_1(5)) +
              (occ_func_0_1(158) * occ_func_0_1(6)) +
              (occ_func_0_1(143) * occ_func_0_1(2)) +
              (occ_func_0_1(156) * occ_func_0_1(4)) +
              (occ_func_0_1(146) * occ_func_0_1(3)) +
              (occ_func_0_1(160) * occ_func_0_1(8)) +
              (occ_func_0_1(159) * occ_func_0_1(7)) +
              (occ_func_0_1(171) * occ_func_0_1(10)) +
              (occ_func_0_1(174) * occ_func_0_1(11)) +
              (occ_func_0_1(161) * occ_func_0_1(9)) +
              (occ_func_0_1(175) * occ_func_0_1(12))) /
             12.0 +
         (m_occ_func_0_1[occ_f] - m_occ_func_0_1[occ_i]) *
             ((occ_func_0_1(54) * occ_func_0_0(175)) +
              (occ_func_0_1(43) * occ_func_0_0(12)) +
              (occ_func_0_1(50) * occ_func_0_0(160)) +
              (occ_func_0_1(47) * occ_func_0_0(8)) +
              (occ_func_0_1(49) * occ_func_0_0(159)) +
              (occ_func_0_1(48) * occ_func_0_0(7)) +
              (occ_func_0_1(53) * occ_func_0_0(174)) +
              (occ_func_0_1(44) * occ_func_0_0(11)) +
              (occ_func_0_1(51) * occ_func_0_0(161)) +
              (occ_func_0_1(46) * occ_func_0_0(9)) +
              (occ_func_0_1(52) * occ_func_0_0(171)) +
              (occ_func_0_1(45) * occ_func_0_0(10)) +
              (occ_func_0_1(47) * occ_func_0_0(157)) +
              (occ_func_0_1(50) * occ_func_0_0(5)) +
              (occ_func_0_1(48) * occ_func_0_0(158)) +
              (occ_func_0_1(49) * occ_func_0_0(6)) +
              (occ_func_0_1(45) * occ_func_0_0(146)) +
              (occ_func_0_1(52) * occ_func_0_0(3)) +
              (occ_func_0_1(44) * occ_func_0_0(143)) +
              (occ_func_0_1(53) * occ_func_0_0(2)) +
              (occ_func_0_1(46) * occ_func_0_0(156)) +
              (occ_func_0_1(51) * occ_func_0_0(4)) +
              (occ_func_0_1(43) * occ_func_0_0(142)) +
              (occ_func_0_1(54) * occ_func_0_0(1))) /
             12.0;
}
double test_Clexulator::delta_site_eval_at_0_bfunc_3_4_4(int occ_i,
                                                         int occ_f) const {
  return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) *
             ((occ_func_0_0(54) * occ_func_0_1(175)) +
              (occ_func_0_0(43) * occ_func_0_1(12)) +
              (occ_func_0_0(50) * occ_func_0_1(160)) +
              (occ_func_0_0(47) * occ_func_0_1(8)) +
              (occ_func_0_0(49) * occ_func_0_1(159)) +
              (occ_func_0_0(48) * occ_func_0_1(7)) +
              (occ_func_0_0(53) * occ_func_0_1(174)) +
              (occ_func_0_0(44) * occ_func_0_1(11)) +
              (occ_func_0_0(51) * occ_func_0_1(161)) +
              (occ_func_0_0(46) * occ_func_0_1(9)) +
              (occ_func_0_0(52) * occ_func_0_1(171)) +
              (occ_func_0_0(45) * occ_func_0_1(10)) +
              (occ_func_0_0(47) * occ_func_0_1(157)) +
              (occ_func_0_0(50) * occ_func_0_1(5)) +
              (occ_func_0_0(48) * occ_func_0_1(158)) +
              (occ_func_0_0(49) * occ_func_0_1(6)) +
              (occ_func_0_0(45) * occ_func_0_1(146)) +
              (occ_func_0_0(52) * occ_func_0_1(3)) +
              (occ_func_0_0(44) * occ_func_0_1(143)) +
              (occ_func_0_0(53) * occ_func_0_1(2)) +
              (occ_func_0_0(46) * occ_func_0_1(156)) +
              (occ_func_0_0(51) * occ_func_0_1(4)) +
              (occ_func_0_0(43) * occ_func_0_1(142)) +
              (occ_func_0_0(54) * occ_func_0_1(1))) /
             12.0 +
         (m_occ_func_0_1[occ_f] - m_occ_func_0_1[occ_i]) *
             ((occ_func_0_0(142) * occ_func_0_0(1)) +
              (occ_func_0_0(157) * occ_func_0_0(5)) +
              (occ_func_0_0(158) * occ_func_0_0(6)) +
              (occ_func_0_0(143) * occ_func_0_0(2)) +
              (occ_func_0_0(156) * occ_func_0_0(4)) +
              (occ_func_0_0(146) * occ_func_0_0(3)) +
              (occ_func_0_0(160) * occ_func_0_0(8)) +
              (occ_func_0_0(159) * occ_func_0_0(7)) +
              (occ_func_0_0(171) * occ_func_0_0(10)) +
              (occ_func_0_0(174) * occ_func_0_0(11)) +
              (occ_func_0_0(161) * occ_func_0_0(9)) +
              (occ_func_0_0(175) * occ_func_0_0(12))) /
             12.0;
}
double test_Clexulator::delta_site_eval_at_0_bfunc_3_4_5(int occ_i,
                                                         int occ_f) const {
  return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) *
             ((occ_func_0_1(43) * occ_func_0_1(12)) +
              (occ_func_0_1(47) * occ_func_0_1(8)) +
              (occ_func_0_1(48) * occ_func_0_1(7)) +
              (occ_func_0_1(44) * occ_func_0_1(11)) +
              (occ_func_0_1(46) * occ_func_0_1(9)) +
              (occ_func_0_1(45) * occ_func_0_1(10)) +
              (occ_func_0_1(50) * occ_func_0_1(5)) +
              (occ_func_0_1(49) * occ_func_0_1(6)) +
              (occ_func_0_1(52) * occ_func_0_1(3)) +
              (occ_func_0_1(53) * occ_func_0_1(2)) +
              (occ_func_0_1(51) * occ_func_0_1(4)) +
              (occ_func_0_1(54) * occ_func_0_1(1))) /
             12.0 +
         (m_occ_func_0_1[occ_f] - m_occ_func_0_1[occ_i]) *
             ((occ_func_0_0(54) * occ_func_0_1(175)) +
              (occ_func_0_1(142) * occ_func_0_0(1)) +
              (occ_func_0_0(50) * occ_func_0_1(160)) +
              (occ_func_0_1(157) * occ_func_0_0(5)) +
              (occ_func_0_0(49) * occ_func_0_1(159)) +
              (occ_func_0_1(158) * occ_func_0_0(6)) +
              (occ_func_0_0(53) * occ_func_0_1(174)) +
              (occ_func_0_1(143) * occ_func_0_0(2)) +
              (occ_func_0_0(51) * occ_func_0_1(161)) +
              (occ_func_0_1(156) * occ_func_0_0(4)) +
              (occ_func_0_0(52) * occ_func_0_1(171)) +
              (occ_func_0_1(146) * occ_func_0_0(3)) +
              (occ_func_0_0(47) * occ_func_0_1(157)) +
              (occ_func_0_1(160) * occ_func_0_0(8)) +
              (occ_func_0_0(48) * occ_func_0_1(158)) +
              (occ_func_0_1(159) * occ_func_0_0(7)) +
              (occ_func_0_0(45) * occ_func_0_1(146)) +
              (occ_func_0_1(171) * occ_func_0_0(10)) +
              (occ_func_0_0(44) * occ_func_0_1(143)) +
              (occ_func_0_1(174) * occ_func_0_0(11)) +
              (occ_func_0_0(46) * occ_func_0_1(156)) +
              (occ_func_0_1(161) * occ_func_0_0(9)) +
              (occ_func_0_0(43) * occ_func_0_1(142)) +
              (occ_func_0_1(175) * occ_func_0_0(12))) /
             12.0;
}
double test_Clexulator::delta_site_eval_at_0_bfunc_3_4_6(int occ_i,
                                                         int occ_f) const {
  return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) *
             ((occ_func_0_1(54) * occ_func_0_1(175)) +
              (occ_func_0_1(50) * occ_func_0_1(160)) +
              (occ_func_0_1(49) * occ_func_0_1(159)) +
              (occ_func_0_1(53) * occ_func_0_1(174)) +
              (occ_func_0_1(51) * occ_func_0_1(161)) +
              (occ_func_0_1(52) * occ_func_0_1(171)) +
              (occ_func_0_1(47) * occ_func_0_1(157)) +
              (occ_func_0_1(48) * occ_func_0_1(158)) +
              (occ_func_0_1(45) * occ_func_0_1(146)) +
              (occ_func_0_1(44) * occ_func_0_1(143)) +
              (occ_func_0_1(46) * occ_func_0_1(156)) +
              (occ_func_0_1(43) * occ_func_0_1(142))) /
             12.0 +
         (m_occ_func_0_1[occ_f] - m_occ_func_0_1[occ_i]) *
             ((occ_func_0_0(43) * occ_func_0_1(12)) +
              (occ_func_0_0(142) * occ_func_0_1(1)) +
              (occ_func_0_0(47) * occ_func_0_1(8)) +
              (occ_func_0_0(157) * occ_func_0_1(5)) +
              (occ_func_0_0(48) * occ_func_0_1(7)) +
              (occ_func_0_0(158) * occ_func_0_1(6)) +
              (occ_func_0_0(44) * occ_func_0_1(11)) +
              (occ_func_0_0(143) * occ_func_0_1(2)) +
              (occ_func_0_0(46) * occ_func_0_1(9)) +
              (occ_func_0_0(156) * occ_func_0_1(4)) +
              (occ_func_0_0(45) * occ_func_0_1(10)) +
              (occ_func_0_0(146) * occ_func_0_1(3)) +
              (occ_func_0_0(50) * occ_func_0_1(5)) +
              (occ_func_0_0(160) * occ_func_0_1(8)) +
              (occ_func_0_0(49) * occ_func_0_1(6)) +
              (occ_func_0_0(159) * occ_func_0_1(7)) +
              (occ_func_0_0(52) * occ_func_0_1(3)) +
              (occ_func_0_0(171) * occ_func_0_1(10)) +
              (occ_func_0_0(53) * occ_func_0_1(2)) +
              (occ_func_0_0(174) * occ_func_0_1(11)) +
              (occ_func_0_0(51) * occ_func_0_1(4)) +
              (occ_func_0_0(161) * occ_func_0_1(9)) +
              (occ_func_0_0(54) * occ_func_0_1(1)) +
              (occ_func_0_0(175) * occ_func_0_1(12))) /
             12.0;
}
double test_Clexulator::delta_site_eval_at_0_bfunc_3_4_7(int occ_i,
                                                         int occ_f) const {
  return (m_occ_func_0_1[occ_f] - m_occ_func_0_1[occ_i]) *
         ((occ_func_0_1(54) * occ_func_0_1(175)) +
          (occ_func_0_1(43) * occ_func_0_1(12)) +
          (occ_func_0_1(142) * occ_func_0_1(1)) +
          (occ_func_0_1(50) * occ_func_0_1(160)) +
          (occ_func_0_1(47) * occ_func_0_1(8)) +
          (occ_func_0_1(157) * occ_func_0_1(5)) +
          (occ_func_0_1(49) * occ_func_0_1(159)) +
          (occ_func_0_1(48) * occ_func_0_1(7)) +
          (occ_func_0_1(158) * occ_func_0_1(6)) +
          (occ_func_0_1(53) * occ_func_0_1(174)) +
          (occ_func_0_1(44) * occ_func_0_1(11)) +
          (occ_func_0_1(143) * occ_func_0_1(2)) +
          (occ_func_0_1(51) * occ_func_0_1(161)) +
          (occ_func_0_1(46) * occ_func_0_1(9)) +
          (occ_func_0_1(156) * occ_func_0_1(4)) +
          (occ_func_0_1(52) * occ_func_0_1(171)) +
          (occ_func_0_1(45) * occ_func_0_1(10)) +
          (occ_func_0_1(146) * occ_func_0_1(3)) +
          (occ_func_0_1(47) * occ_func_0_1(157)) +
          (occ_func_0_1(50) * occ_func_0_1(5)) +
          (occ_func_0_1(160) * occ_func_0_1(8)) +
          (occ_func_0_1(48) * occ_func_0_1(158)) +
          (occ_func_0_1(49) * occ_func_0_1(6)) +
          (occ_func_0_1(159) * occ_func_0_1(7)) +
          (occ_func_0_1(45) * occ_func_0_1(146)) +
          (occ_func_0_1(52) * occ_func_0_1(3)) +
          (occ_func_0_1(171) * occ_func_0_1(10)) +
          (occ_func_0_1(44) * occ_func_0_1(143)) +
          (occ_func_0_1(53) * occ_func_0_1(2)) +
          (occ_func_0_1(174) * occ_func_0_1(11)) +
          (occ_func_0_1(46) * occ_func_0_1(156)) +
          (occ_func_0_1(51) * occ_func_0_1(4)) +
          (occ_func_0_1(161) * occ_func_0_1(9)) +
          (occ_func_0_1(43) * occ_func_0_1(142)) +
          (occ_func_0_1(54) * occ_func_0_1(1)) +
          (occ_func_0_1(175) * occ_func_0_1(12))) /
         12.0;
}

/**** Basis functions for orbit 4, 0****
#Points: 4
MaxLength: 6.9282032  MinLength: 2.8284271
 0.0000000   0.0000000   0.0000000 A B C
 0.0000000   1.0000000   0.0000000 A B C
 0.0000000   0.0000000   1.0000000 A B C
 1.0000000   1.0000000   1.0000000 A B C
****/
double test_Clexulator::eval_bfunc_4_0_0() const {
  return ((occ_func_0_0(0) * occ_func_0_0(9) * occ_func_0_0(7) *
           occ_func_0_0(85)) +
          (occ_func_0_0(0) * occ_func_0_0(3) * occ_func_0_0(9) *
           occ_func_0_0(82)) +
          (occ_func_0_0(0) * occ_func_0_0(2) * occ_func_0_0(5) *
           occ_func_0_0(81)) +
          (occ_func_0_0(0) * occ_func_0_0(12) * occ_func_0_0(10) *
           occ_func_0_0(86)) +
          (occ_func_0_0(0) * occ_func_0_0(8) * occ_func_0_0(3) *
           occ_func_0_0(82)) +
          (occ_func_0_0(0) * occ_func_0_0(11) * occ_func_0_0(12) *
           occ_func_0_0(86)) +
          (occ_func_0_0(0) * occ_func_0_0(7) * occ_func_0_0(2) *
           occ_func_0_0(81)) +
          (occ_func_0_0(0) * occ_func_0_0(6) * occ_func_0_0(8) *
           occ_func_0_0(84)) +
          (occ_func_0_0(0) * occ_func_0_0(8) * occ_func_0_0(11) *
           occ_func_0_0(84)) +
          (occ_func_0_0(0) * occ_func_0_0(12) * occ_func_0_0(9) *
           occ_func_0_0(85)) +
          (occ_func_0_0(0) * occ_func_0_0(5) * occ_func_0_0(4) *
           occ_func_0_0(83)) +
          (occ_func_0_0(0) * occ_func_0_0(7) * occ_func_0_0(12) *
           occ_func_0_0(85)) +
          (occ_func_0_0(0) * occ_func_0_0(1) * occ_func_0_0(2) *
           occ_func_0_0(79)) +
          (occ_func_0_0(0) * occ_func_0_0(3) * occ_func_0_0(1) *
           occ_func_0_0(79)) +
          (occ_func_0_0(0) * occ_func_0_0(10) * occ_func_0_0(5) *
           occ_func_0_0(83)) +
          (occ_func_0_0(0) * occ_func_0_0(10) * occ_func_0_0(11) *
           occ_func_0_0(86)) +
          (occ_func_0_0(0) * occ_func_0_0(6) * occ_func_0_0(4) *
           occ_func_0_0(80)) +
          (occ_func_0_0(0) * occ_func_0_0(9) * occ_func_0_0(8) *
           occ_func_0_0(82)) +
          (occ_func_0_0(0) * occ_func_0_0(2) * occ_func_0_0(3) *
           occ_func_0_0(79)) +
          (occ_func_0_0(0) * occ_func_0_0(11) * occ_func_0_0(6) *
           occ_func_0_0(84)) +
          (occ_func_0_0(0) * occ_func_0_0(5) * occ_func_0_0(7) *
           occ_func_0_0(81)) +
          (occ_func_0_0(0) * occ_func_0_0(1) * occ_func_0_0(6) *
           occ_func_0_0(80)) +
          (occ_func_0_0(0) * occ_func_0_0(4) * occ_func_0_0(10) *
           occ_func_0_0(83)) +
          (occ_func_0_0(0) * occ_func_0_0(4) * occ_func_0_0(1) *
           occ_func_0_0(80))) /
         24.0;
}
double test_Clexulator::eval_bfunc_4_0_1() const {
  return ((occ_func_0_1(0) * occ_func_0_0(9) * occ_func_0_0(7) *
           occ_func_0_0(85)) +
          (occ_func_0_1(0) * occ_func_0_0(3) * occ_func_0_0(9) *
           occ_func_0_0(82)) +
          (occ_func_0_1(0) * occ_func_0_0(2) * occ_func_0_0(5) *
           occ_func_0_0(81)) +
          (occ_func_0_1(0) * occ_func_0_0(12) * occ_func_0_0(10) *
           occ_func_0_0(86)) +
          (occ_func_0_1(0) * occ_func_0_0(8) * occ_func_0_0(3) *
           occ_func_0_0(82)) +
          (occ_func_0_1(0) * occ_func_0_0(11) * occ_func_0_0(12) *
           occ_func_0_0(86)) +
          (occ_func_0_1(0) * occ_func_0_0(7) * occ_func_0_0(2) *
           occ_func_0_0(81)) +
          (occ_func_0_1(0) * occ_func_0_0(6) * occ_func_0_0(8) *
           occ_func_0_0(84)) +
          (occ_func_0_1(0) * occ_func_0_0(8) * occ_func_0_0(11) *
           occ_func_0_0(84)) +
          (occ_func_0_1(0) * occ_func_0_0(12) * occ_func_0_0(9) *
           occ_func_0_0(85)) +
          (occ_func_0_1(0) * occ_func_0_0(5) * occ_func_0_0(4) *
           occ_func_0_0(83)) +
          (occ_func_0_1(0) * occ_func_0_0(7) * occ_func_0_0(12) *
           occ_func_0_0(85)) +
          (occ_func_0_1(0) * occ_func_0_0(1) * occ_func_0_0(2) *
           occ_func_0_0(79)) +
          (occ_func_0_1(0) * occ_func_0_0(3) * occ_func_0_0(1) *
           occ_func_0_0(79)) +
          (occ_func_0_1(0) * occ_func_0_0(10) * occ_func_0_0(5) *
           occ_func_0_0(83)) +
          (occ_func_0_1(0) * occ_func_0_0(10) * occ_func_0_0(11) *
           occ_func_0_0(86)) +
          (occ_func_0_1(0) * occ_func_0_0(6) * occ_func_0_0(4) *
           occ_func_0_0(80)) +
          (occ_func_0_1(0) * occ_func_0_0(9) * occ_func_0_0(8) *
           occ_func_0_0(82)) +
          (occ_func_0_1(0) * occ_func_0_0(2) * occ_func_0_0(3) *
           occ_func_0_0(79)) +
          (occ_func_0_1(0) * occ_func_0_0(11) * occ_func_0_0(6) *
           occ_func_0_0(84)) +
          (occ_func_0_1(0) * occ_func_0_0(5) * occ_func_0_0(7) *
           occ_func_0_0(81)) +
          (occ_func_0_1(0) * occ_func_0_0(1) * occ_func_0_0(6) *
           occ_func_0_0(80)) +
          (occ_func_0_1(0) * occ_func_0_0(4) * occ_func_0_0(10) *
           occ_func_0_0(83)) +
          (occ_func_0_1(0) * occ_func_0_0(4) * occ_func_0_0(1) *
           occ_func_0_0(80))) /
         24.0;
}
double test_Clexulator::eval_bfunc_4_0_2() const {
  return (((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(9) * occ_func_0_0(7) *
                occ_func_0_0(85) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(9) * occ_func_0_1(7) *
                occ_func_0_0(85))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(3) * occ_func_0_0(9) *
                occ_func_0_0(82) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(3) * occ_func_0_1(9) *
                occ_func_0_0(82))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(2) * occ_func_0_0(5) *
                occ_func_0_0(81) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(2) * occ_func_0_1(5) *
                occ_func_0_0(81))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(12) *
                occ_func_0_0(10) * occ_func_0_0(86) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(12) *
                occ_func_0_1(10) * occ_func_0_0(86))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(8) * occ_func_0_0(3) *
                occ_func_0_0(82) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(8) * occ_func_0_1(3) *
                occ_func_0_0(82))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(11) *
                occ_func_0_0(12) * occ_func_0_0(86) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(11) *
                occ_func_0_1(12) * occ_func_0_0(86))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(7) * occ_func_0_0(2) *
                occ_func_0_0(81) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(7) * occ_func_0_1(2) *
                occ_func_0_0(81))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(6) * occ_func_0_0(8) *
                occ_func_0_0(84) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(6) * occ_func_0_1(8) *
                occ_func_0_0(84))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(8) *
                occ_func_0_0(11) * occ_func_0_0(84) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(8) *
                occ_func_0_1(11) * occ_func_0_0(84))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(12) *
                occ_func_0_0(9) * occ_func_0_0(85) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(12) *
                occ_func_0_1(9) * occ_func_0_0(85))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(5) * occ_func_0_0(4) *
                occ_func_0_0(83) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(5) * occ_func_0_1(4) *
                occ_func_0_0(83))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(7) *
                occ_func_0_0(12) * occ_func_0_0(85) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(7) *
                occ_func_0_1(12) * occ_func_0_0(85))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(1) * occ_func_0_0(2) *
                occ_func_0_0(79) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(1) * occ_func_0_1(2) *
                occ_func_0_0(79))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(3) * occ_func_0_0(1) *
                occ_func_0_0(79) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(3) * occ_func_0_1(1) *
                occ_func_0_0(79))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(10) *
                occ_func_0_0(5) * occ_func_0_0(83) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(10) *
                occ_func_0_1(5) * occ_func_0_0(83))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(10) *
                occ_func_0_0(11) * occ_func_0_0(86) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(10) *
                occ_func_0_1(11) * occ_func_0_0(86))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(6) * occ_func_0_0(4) *
                occ_func_0_0(80) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(6) * occ_func_0_1(4) *
                occ_func_0_0(80))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(9) * occ_func_0_0(8) *
                occ_func_0_0(82) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(9) * occ_func_0_1(8) *
                occ_func_0_0(82))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(2) * occ_func_0_0(3) *
                occ_func_0_0(79) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(2) * occ_func_0_1(3) *
                occ_func_0_0(79))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(11) *
                occ_func_0_0(6) * occ_func_0_0(84) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(11) *
                occ_func_0_1(6) * occ_func_0_0(84))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(5) * occ_func_0_0(7) *
                occ_func_0_0(81) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(5) * occ_func_0_1(7) *
                occ_func_0_0(81))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(1) * occ_func_0_0(6) *
                occ_func_0_0(80) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(1) * occ_func_0_1(6) *
                occ_func_0_0(80))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(4) *
                occ_func_0_0(10) * occ_func_0_0(83) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(4) *
                occ_func_0_1(10) * occ_func_0_0(83))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(4) * occ_func_0_0(1) *
                occ_func_0_0(80) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(4) * occ_func_0_1(1) *
                occ_func_0_0(80)))) /
         24.0;
}
double test_Clexulator::eval_bfunc_4_0_3() const {
  return (((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(9) * occ_func_0_0(7) *
                occ_func_0_0(85) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(9) * occ_func_0_1(7) *
                occ_func_0_0(85))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(3) * occ_func_0_0(9) *
                occ_func_0_0(82) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(3) * occ_func_0_1(9) *
                occ_func_0_0(82))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(2) * occ_func_0_0(5) *
                occ_func_0_0(81) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(2) * occ_func_0_1(5) *
                occ_func_0_0(81))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(12) *
                occ_func_0_0(10) * occ_func_0_0(86) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(12) *
                occ_func_0_1(10) * occ_func_0_0(86))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(8) * occ_func_0_0(3) *
                occ_func_0_0(82) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(8) * occ_func_0_1(3) *
                occ_func_0_0(82))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(11) *
                occ_func_0_0(12) * occ_func_0_0(86) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(11) *
                occ_func_0_1(12) * occ_func_0_0(86))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(7) * occ_func_0_0(2) *
                occ_func_0_0(81) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(7) * occ_func_0_1(2) *
                occ_func_0_0(81))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(6) * occ_func_0_0(8) *
                occ_func_0_0(84) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(6) * occ_func_0_1(8) *
                occ_func_0_0(84))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(8) *
                occ_func_0_0(11) * occ_func_0_0(84) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(8) *
                occ_func_0_1(11) * occ_func_0_0(84))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(12) *
                occ_func_0_0(9) * occ_func_0_0(85) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(12) *
                occ_func_0_1(9) * occ_func_0_0(85))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(5) * occ_func_0_0(4) *
                occ_func_0_0(83) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(5) * occ_func_0_1(4) *
                occ_func_0_0(83))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(7) *
                occ_func_0_0(12) * occ_func_0_0(85) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(7) *
                occ_func_0_1(12) * occ_func_0_0(85))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(1) * occ_func_0_0(2) *
                occ_func_0_0(79) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(1) * occ_func_0_1(2) *
                occ_func_0_0(79))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(3) * occ_func_0_0(1) *
                occ_func_0_0(79) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(3) * occ_func_0_1(1) *
                occ_func_0_0(79))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(10) *
                occ_func_0_0(5) * occ_func_0_0(83) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(10) *
                occ_func_0_1(5) * occ_func_0_0(83))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(10) *
                occ_func_0_0(11) * occ_func_0_0(86) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(10) *
                occ_func_0_1(11) * occ_func_0_0(86))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(6) * occ_func_0_0(4) *
                occ_func_0_0(80) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(6) * occ_func_0_1(4) *
                occ_func_0_0(80))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(9) * occ_func_0_0(8) *
                occ_func_0_0(82) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(9) * occ_func_0_1(8) *
                occ_func_0_0(82))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(2) * occ_func_0_0(3) *
                occ_func_0_0(79) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(2) * occ_func_0_1(3) *
                occ_func_0_0(79))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(11) *
                occ_func_0_0(6) * occ_func_0_0(84) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(11) *
                occ_func_0_1(6) * occ_func_0_0(84))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(5) * occ_func_0_0(7) *
                occ_func_0_0(81) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(5) * occ_func_0_1(7) *
                occ_func_0_0(81))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(1) * occ_func_0_0(6) *
                occ_func_0_0(80) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(1) * occ_func_0_1(6) *
                occ_func_0_0(80))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(4) *
                occ_func_0_0(10) * occ_func_0_0(83) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(4) *
                occ_func_0_1(10) * occ_func_0_0(83))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(4) * occ_func_0_0(1) *
                occ_func_0_0(80) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(4) * occ_func_0_1(1) *
                occ_func_0_0(80)))) /
         24.0;
}
double test_Clexulator::eval_bfunc_4_0_4() const {
  return ((occ_func_0_0(0) * occ_func_0_1(9) * occ_func_0_1(7) *
           occ_func_0_0(85)) +
          (occ_func_0_0(0) * occ_func_0_1(3) * occ_func_0_1(9) *
           occ_func_0_0(82)) +
          (occ_func_0_0(0) * occ_func_0_1(2) * occ_func_0_1(5) *
           occ_func_0_0(81)) +
          (occ_func_0_0(0) * occ_func_0_1(12) * occ_func_0_1(10) *
           occ_func_0_0(86)) +
          (occ_func_0_0(0) * occ_func_0_1(8) * occ_func_0_1(3) *
           occ_func_0_0(82)) +
          (occ_func_0_0(0) * occ_func_0_1(11) * occ_func_0_1(12) *
           occ_func_0_0(86)) +
          (occ_func_0_0(0) * occ_func_0_1(7) * occ_func_0_1(2) *
           occ_func_0_0(81)) +
          (occ_func_0_0(0) * occ_func_0_1(6) * occ_func_0_1(8) *
           occ_func_0_0(84)) +
          (occ_func_0_0(0) * occ_func_0_1(8) * occ_func_0_1(11) *
           occ_func_0_0(84)) +
          (occ_func_0_0(0) * occ_func_0_1(12) * occ_func_0_1(9) *
           occ_func_0_0(85)) +
          (occ_func_0_0(0) * occ_func_0_1(5) * occ_func_0_1(4) *
           occ_func_0_0(83)) +
          (occ_func_0_0(0) * occ_func_0_1(7) * occ_func_0_1(12) *
           occ_func_0_0(85)) +
          (occ_func_0_0(0) * occ_func_0_1(1) * occ_func_0_1(2) *
           occ_func_0_0(79)) +
          (occ_func_0_0(0) * occ_func_0_1(3) * occ_func_0_1(1) *
           occ_func_0_0(79)) +
          (occ_func_0_0(0) * occ_func_0_1(10) * occ_func_0_1(5) *
           occ_func_0_0(83)) +
          (occ_func_0_0(0) * occ_func_0_1(10) * occ_func_0_1(11) *
           occ_func_0_0(86)) +
          (occ_func_0_0(0) * occ_func_0_1(6) * occ_func_0_1(4) *
           occ_func_0_0(80)) +
          (occ_func_0_0(0) * occ_func_0_1(9) * occ_func_0_1(8) *
           occ_func_0_0(82)) +
          (occ_func_0_0(0) * occ_func_0_1(2) * occ_func_0_1(3) *
           occ_func_0_0(79)) +
          (occ_func_0_0(0) * occ_func_0_1(11) * occ_func_0_1(6) *
           occ_func_0_0(84)) +
          (occ_func_0_0(0) * occ_func_0_1(5) * occ_func_0_1(7) *
           occ_func_0_0(81)) +
          (occ_func_0_0(0) * occ_func_0_1(1) * occ_func_0_1(6) *
           occ_func_0_0(80)) +
          (occ_func_0_0(0) * occ_func_0_1(4) * occ_func_0_1(10) *
           occ_func_0_0(83)) +
          (occ_func_0_0(0) * occ_func_0_1(4) * occ_func_0_1(1) *
           occ_func_0_0(80))) /
         24.0;
}
double test_Clexulator::eval_bfunc_4_0_5() const {
  return ((occ_func_0_1(0) * occ_func_0_1(9) * occ_func_0_1(7) *
           occ_func_0_0(85)) +
          (occ_func_0_1(0) * occ_func_0_1(3) * occ_func_0_1(9) *
           occ_func_0_0(82)) +
          (occ_func_0_1(0) * occ_func_0_1(2) * occ_func_0_1(5) *
           occ_func_0_0(81)) +
          (occ_func_0_1(0) * occ_func_0_1(12) * occ_func_0_1(10) *
           occ_func_0_0(86)) +
          (occ_func_0_1(0) * occ_func_0_1(8) * occ_func_0_1(3) *
           occ_func_0_0(82)) +
          (occ_func_0_1(0) * occ_func_0_1(11) * occ_func_0_1(12) *
           occ_func_0_0(86)) +
          (occ_func_0_1(0) * occ_func_0_1(7) * occ_func_0_1(2) *
           occ_func_0_0(81)) +
          (occ_func_0_1(0) * occ_func_0_1(6) * occ_func_0_1(8) *
           occ_func_0_0(84)) +
          (occ_func_0_1(0) * occ_func_0_1(8) * occ_func_0_1(11) *
           occ_func_0_0(84)) +
          (occ_func_0_1(0) * occ_func_0_1(12) * occ_func_0_1(9) *
           occ_func_0_0(85)) +
          (occ_func_0_1(0) * occ_func_0_1(5) * occ_func_0_1(4) *
           occ_func_0_0(83)) +
          (occ_func_0_1(0) * occ_func_0_1(7) * occ_func_0_1(12) *
           occ_func_0_0(85)) +
          (occ_func_0_1(0) * occ_func_0_1(1) * occ_func_0_1(2) *
           occ_func_0_0(79)) +
          (occ_func_0_1(0) * occ_func_0_1(3) * occ_func_0_1(1) *
           occ_func_0_0(79)) +
          (occ_func_0_1(0) * occ_func_0_1(10) * occ_func_0_1(5) *
           occ_func_0_0(83)) +
          (occ_func_0_1(0) * occ_func_0_1(10) * occ_func_0_1(11) *
           occ_func_0_0(86)) +
          (occ_func_0_1(0) * occ_func_0_1(6) * occ_func_0_1(4) *
           occ_func_0_0(80)) +
          (occ_func_0_1(0) * occ_func_0_1(9) * occ_func_0_1(8) *
           occ_func_0_0(82)) +
          (occ_func_0_1(0) * occ_func_0_1(2) * occ_func_0_1(3) *
           occ_func_0_0(79)) +
          (occ_func_0_1(0) * occ_func_0_1(11) * occ_func_0_1(6) *
           occ_func_0_0(84)) +
          (occ_func_0_1(0) * occ_func_0_1(5) * occ_func_0_1(7) *
           occ_func_0_0(81)) +
          (occ_func_0_1(0) * occ_func_0_1(1) * occ_func_0_1(6) *
           occ_func_0_0(80)) +
          (occ_func_0_1(0) * occ_func_0_1(4) * occ_func_0_1(10) *
           occ_func_0_0(83)) +
          (occ_func_0_1(0) * occ_func_0_1(4) * occ_func_0_1(1) *
           occ_func_0_0(80))) /
         24.0;
}
double test_Clexulator::eval_bfunc_4_0_6() const {
  return ((occ_func_0_0(0) * occ_func_0_0(9) * occ_func_0_0(7) *
           occ_func_0_1(85)) +
          (occ_func_0_0(0) * occ_func_0_0(3) * occ_func_0_0(9) *
           occ_func_0_1(82)) +
          (occ_func_0_0(0) * occ_func_0_0(2) * occ_func_0_0(5) *
           occ_func_0_1(81)) +
          (occ_func_0_0(0) * occ_func_0_0(12) * occ_func_0_0(10) *
           occ_func_0_1(86)) +
          (occ_func_0_0(0) * occ_func_0_0(8) * occ_func_0_0(3) *
           occ_func_0_1(82)) +
          (occ_func_0_0(0) * occ_func_0_0(11) * occ_func_0_0(12) *
           occ_func_0_1(86)) +
          (occ_func_0_0(0) * occ_func_0_0(7) * occ_func_0_0(2) *
           occ_func_0_1(81)) +
          (occ_func_0_0(0) * occ_func_0_0(6) * occ_func_0_0(8) *
           occ_func_0_1(84)) +
          (occ_func_0_0(0) * occ_func_0_0(8) * occ_func_0_0(11) *
           occ_func_0_1(84)) +
          (occ_func_0_0(0) * occ_func_0_0(12) * occ_func_0_0(9) *
           occ_func_0_1(85)) +
          (occ_func_0_0(0) * occ_func_0_0(5) * occ_func_0_0(4) *
           occ_func_0_1(83)) +
          (occ_func_0_0(0) * occ_func_0_0(7) * occ_func_0_0(12) *
           occ_func_0_1(85)) +
          (occ_func_0_0(0) * occ_func_0_0(1) * occ_func_0_0(2) *
           occ_func_0_1(79)) +
          (occ_func_0_0(0) * occ_func_0_0(3) * occ_func_0_0(1) *
           occ_func_0_1(79)) +
          (occ_func_0_0(0) * occ_func_0_0(10) * occ_func_0_0(5) *
           occ_func_0_1(83)) +
          (occ_func_0_0(0) * occ_func_0_0(10) * occ_func_0_0(11) *
           occ_func_0_1(86)) +
          (occ_func_0_0(0) * occ_func_0_0(6) * occ_func_0_0(4) *
           occ_func_0_1(80)) +
          (occ_func_0_0(0) * occ_func_0_0(9) * occ_func_0_0(8) *
           occ_func_0_1(82)) +
          (occ_func_0_0(0) * occ_func_0_0(2) * occ_func_0_0(3) *
           occ_func_0_1(79)) +
          (occ_func_0_0(0) * occ_func_0_0(11) * occ_func_0_0(6) *
           occ_func_0_1(84)) +
          (occ_func_0_0(0) * occ_func_0_0(5) * occ_func_0_0(7) *
           occ_func_0_1(81)) +
          (occ_func_0_0(0) * occ_func_0_0(1) * occ_func_0_0(6) *
           occ_func_0_1(80)) +
          (occ_func_0_0(0) * occ_func_0_0(4) * occ_func_0_0(10) *
           occ_func_0_1(83)) +
          (occ_func_0_0(0) * occ_func_0_0(4) * occ_func_0_0(1) *
           occ_func_0_1(80))) /
         24.0;
}
double test_Clexulator::eval_bfunc_4_0_7() const {
  return ((occ_func_0_1(0) * occ_func_0_0(9) * occ_func_0_0(7) *
           occ_func_0_1(85)) +
          (occ_func_0_1(0) * occ_func_0_0(3) * occ_func_0_0(9) *
           occ_func_0_1(82)) +
          (occ_func_0_1(0) * occ_func_0_0(2) * occ_func_0_0(5) *
           occ_func_0_1(81)) +
          (occ_func_0_1(0) * occ_func_0_0(12) * occ_func_0_0(10) *
           occ_func_0_1(86)) +
          (occ_func_0_1(0) * occ_func_0_0(8) * occ_func_0_0(3) *
           occ_func_0_1(82)) +
          (occ_func_0_1(0) * occ_func_0_0(11) * occ_func_0_0(12) *
           occ_func_0_1(86)) +
          (occ_func_0_1(0) * occ_func_0_0(7) * occ_func_0_0(2) *
           occ_func_0_1(81)) +
          (occ_func_0_1(0) * occ_func_0_0(6) * occ_func_0_0(8) *
           occ_func_0_1(84)) +
          (occ_func_0_1(0) * occ_func_0_0(8) * occ_func_0_0(11) *
           occ_func_0_1(84)) +
          (occ_func_0_1(0) * occ_func_0_0(12) * occ_func_0_0(9) *
           occ_func_0_1(85)) +
          (occ_func_0_1(0) * occ_func_0_0(5) * occ_func_0_0(4) *
           occ_func_0_1(83)) +
          (occ_func_0_1(0) * occ_func_0_0(7) * occ_func_0_0(12) *
           occ_func_0_1(85)) +
          (occ_func_0_1(0) * occ_func_0_0(1) * occ_func_0_0(2) *
           occ_func_0_1(79)) +
          (occ_func_0_1(0) * occ_func_0_0(3) * occ_func_0_0(1) *
           occ_func_0_1(79)) +
          (occ_func_0_1(0) * occ_func_0_0(10) * occ_func_0_0(5) *
           occ_func_0_1(83)) +
          (occ_func_0_1(0) * occ_func_0_0(10) * occ_func_0_0(11) *
           occ_func_0_1(86)) +
          (occ_func_0_1(0) * occ_func_0_0(6) * occ_func_0_0(4) *
           occ_func_0_1(80)) +
          (occ_func_0_1(0) * occ_func_0_0(9) * occ_func_0_0(8) *
           occ_func_0_1(82)) +
          (occ_func_0_1(0) * occ_func_0_0(2) * occ_func_0_0(3) *
           occ_func_0_1(79)) +
          (occ_func_0_1(0) * occ_func_0_0(11) * occ_func_0_0(6) *
           occ_func_0_1(84)) +
          (occ_func_0_1(0) * occ_func_0_0(5) * occ_func_0_0(7) *
           occ_func_0_1(81)) +
          (occ_func_0_1(0) * occ_func_0_0(1) * occ_func_0_0(6) *
           occ_func_0_1(80)) +
          (occ_func_0_1(0) * occ_func_0_0(4) * occ_func_0_0(10) *
           occ_func_0_1(83)) +
          (occ_func_0_1(0) * occ_func_0_0(4) * occ_func_0_0(1) *
           occ_func_0_1(80))) /
         24.0;
}
double test_Clexulator::eval_bfunc_4_0_8() const {
  return (((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(9) * occ_func_0_0(7) *
                occ_func_0_1(85) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(9) * occ_func_0_1(7) *
                occ_func_0_1(85))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(3) * occ_func_0_0(9) *
                occ_func_0_1(82) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(3) * occ_func_0_1(9) *
                occ_func_0_1(82))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(2) * occ_func_0_0(5) *
                occ_func_0_1(81) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(2) * occ_func_0_1(5) *
                occ_func_0_1(81))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(12) *
                occ_func_0_0(10) * occ_func_0_1(86) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(12) *
                occ_func_0_1(10) * occ_func_0_1(86))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(8) * occ_func_0_0(3) *
                occ_func_0_1(82) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(8) * occ_func_0_1(3) *
                occ_func_0_1(82))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(11) *
                occ_func_0_0(12) * occ_func_0_1(86) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(11) *
                occ_func_0_1(12) * occ_func_0_1(86))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(7) * occ_func_0_0(2) *
                occ_func_0_1(81) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(7) * occ_func_0_1(2) *
                occ_func_0_1(81))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(6) * occ_func_0_0(8) *
                occ_func_0_1(84) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(6) * occ_func_0_1(8) *
                occ_func_0_1(84))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(8) *
                occ_func_0_0(11) * occ_func_0_1(84) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(8) *
                occ_func_0_1(11) * occ_func_0_1(84))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(12) *
                occ_func_0_0(9) * occ_func_0_1(85) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(12) *
                occ_func_0_1(9) * occ_func_0_1(85))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(5) * occ_func_0_0(4) *
                occ_func_0_1(83) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(5) * occ_func_0_1(4) *
                occ_func_0_1(83))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(7) *
                occ_func_0_0(12) * occ_func_0_1(85) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(7) *
                occ_func_0_1(12) * occ_func_0_1(85))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(1) * occ_func_0_0(2) *
                occ_func_0_1(79) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(1) * occ_func_0_1(2) *
                occ_func_0_1(79))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(3) * occ_func_0_0(1) *
                occ_func_0_1(79) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(3) * occ_func_0_1(1) *
                occ_func_0_1(79))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(10) *
                occ_func_0_0(5) * occ_func_0_1(83) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(10) *
                occ_func_0_1(5) * occ_func_0_1(83))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(10) *
                occ_func_0_0(11) * occ_func_0_1(86) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(10) *
                occ_func_0_1(11) * occ_func_0_1(86))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(6) * occ_func_0_0(4) *
                occ_func_0_1(80) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(6) * occ_func_0_1(4) *
                occ_func_0_1(80))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(9) * occ_func_0_0(8) *
                occ_func_0_1(82) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(9) * occ_func_0_1(8) *
                occ_func_0_1(82))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(2) * occ_func_0_0(3) *
                occ_func_0_1(79) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(2) * occ_func_0_1(3) *
                occ_func_0_1(79))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(11) *
                occ_func_0_0(6) * occ_func_0_1(84) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(11) *
                occ_func_0_1(6) * occ_func_0_1(84))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(5) * occ_func_0_0(7) *
                occ_func_0_1(81) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(5) * occ_func_0_1(7) *
                occ_func_0_1(81))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(1) * occ_func_0_0(6) *
                occ_func_0_1(80) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(1) * occ_func_0_1(6) *
                occ_func_0_1(80))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(4) *
                occ_func_0_0(10) * occ_func_0_1(83) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(4) *
                occ_func_0_1(10) * occ_func_0_1(83))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(4) * occ_func_0_0(1) *
                occ_func_0_1(80) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(4) * occ_func_0_1(1) *
                occ_func_0_1(80)))) /
         24.0;
}
double test_Clexulator::eval_bfunc_4_0_9() const {
  return (((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(9) * occ_func_0_0(7) *
                occ_func_0_1(85) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(9) * occ_func_0_1(7) *
                occ_func_0_1(85))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(3) * occ_func_0_0(9) *
                occ_func_0_1(82) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(3) * occ_func_0_1(9) *
                occ_func_0_1(82))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(2) * occ_func_0_0(5) *
                occ_func_0_1(81) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(2) * occ_func_0_1(5) *
                occ_func_0_1(81))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(12) *
                occ_func_0_0(10) * occ_func_0_1(86) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(12) *
                occ_func_0_1(10) * occ_func_0_1(86))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(8) * occ_func_0_0(3) *
                occ_func_0_1(82) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(8) * occ_func_0_1(3) *
                occ_func_0_1(82))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(11) *
                occ_func_0_0(12) * occ_func_0_1(86) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(11) *
                occ_func_0_1(12) * occ_func_0_1(86))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(7) * occ_func_0_0(2) *
                occ_func_0_1(81) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(7) * occ_func_0_1(2) *
                occ_func_0_1(81))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(6) * occ_func_0_0(8) *
                occ_func_0_1(84) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(6) * occ_func_0_1(8) *
                occ_func_0_1(84))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(8) *
                occ_func_0_0(11) * occ_func_0_1(84) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(8) *
                occ_func_0_1(11) * occ_func_0_1(84))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(12) *
                occ_func_0_0(9) * occ_func_0_1(85) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(12) *
                occ_func_0_1(9) * occ_func_0_1(85))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(5) * occ_func_0_0(4) *
                occ_func_0_1(83) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(5) * occ_func_0_1(4) *
                occ_func_0_1(83))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(7) *
                occ_func_0_0(12) * occ_func_0_1(85) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(7) *
                occ_func_0_1(12) * occ_func_0_1(85))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(1) * occ_func_0_0(2) *
                occ_func_0_1(79) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(1) * occ_func_0_1(2) *
                occ_func_0_1(79))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(3) * occ_func_0_0(1) *
                occ_func_0_1(79) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(3) * occ_func_0_1(1) *
                occ_func_0_1(79))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(10) *
                occ_func_0_0(5) * occ_func_0_1(83) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(10) *
                occ_func_0_1(5) * occ_func_0_1(83))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(10) *
                occ_func_0_0(11) * occ_func_0_1(86) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(10) *
                occ_func_0_1(11) * occ_func_0_1(86))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(6) * occ_func_0_0(4) *
                occ_func_0_1(80) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(6) * occ_func_0_1(4) *
                occ_func_0_1(80))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(9) * occ_func_0_0(8) *
                occ_func_0_1(82) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(9) * occ_func_0_1(8) *
                occ_func_0_1(82))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(2) * occ_func_0_0(3) *
                occ_func_0_1(79) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(2) * occ_func_0_1(3) *
                occ_func_0_1(79))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(11) *
                occ_func_0_0(6) * occ_func_0_1(84) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(11) *
                occ_func_0_1(6) * occ_func_0_1(84))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(5) * occ_func_0_0(7) *
                occ_func_0_1(81) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(5) * occ_func_0_1(7) *
                occ_func_0_1(81))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(1) * occ_func_0_0(6) *
                occ_func_0_1(80) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(1) * occ_func_0_1(6) *
                occ_func_0_1(80))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(4) *
                occ_func_0_0(10) * occ_func_0_1(83) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(4) *
                occ_func_0_1(10) * occ_func_0_1(83))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(4) * occ_func_0_0(1) *
                occ_func_0_1(80) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(4) * occ_func_0_1(1) *
                occ_func_0_1(80)))) /
         24.0;
}
double test_Clexulator::eval_bfunc_4_0_10() const {
  return ((occ_func_0_0(0) * occ_func_0_1(9) * occ_func_0_1(7) *
           occ_func_0_1(85)) +
          (occ_func_0_0(0) * occ_func_0_1(3) * occ_func_0_1(9) *
           occ_func_0_1(82)) +
          (occ_func_0_0(0) * occ_func_0_1(2) * occ_func_0_1(5) *
           occ_func_0_1(81)) +
          (occ_func_0_0(0) * occ_func_0_1(12) * occ_func_0_1(10) *
           occ_func_0_1(86)) +
          (occ_func_0_0(0) * occ_func_0_1(8) * occ_func_0_1(3) *
           occ_func_0_1(82)) +
          (occ_func_0_0(0) * occ_func_0_1(11) * occ_func_0_1(12) *
           occ_func_0_1(86)) +
          (occ_func_0_0(0) * occ_func_0_1(7) * occ_func_0_1(2) *
           occ_func_0_1(81)) +
          (occ_func_0_0(0) * occ_func_0_1(6) * occ_func_0_1(8) *
           occ_func_0_1(84)) +
          (occ_func_0_0(0) * occ_func_0_1(8) * occ_func_0_1(11) *
           occ_func_0_1(84)) +
          (occ_func_0_0(0) * occ_func_0_1(12) * occ_func_0_1(9) *
           occ_func_0_1(85)) +
          (occ_func_0_0(0) * occ_func_0_1(5) * occ_func_0_1(4) *
           occ_func_0_1(83)) +
          (occ_func_0_0(0) * occ_func_0_1(7) * occ_func_0_1(12) *
           occ_func_0_1(85)) +
          (occ_func_0_0(0) * occ_func_0_1(1) * occ_func_0_1(2) *
           occ_func_0_1(79)) +
          (occ_func_0_0(0) * occ_func_0_1(3) * occ_func_0_1(1) *
           occ_func_0_1(79)) +
          (occ_func_0_0(0) * occ_func_0_1(10) * occ_func_0_1(5) *
           occ_func_0_1(83)) +
          (occ_func_0_0(0) * occ_func_0_1(10) * occ_func_0_1(11) *
           occ_func_0_1(86)) +
          (occ_func_0_0(0) * occ_func_0_1(6) * occ_func_0_1(4) *
           occ_func_0_1(80)) +
          (occ_func_0_0(0) * occ_func_0_1(9) * occ_func_0_1(8) *
           occ_func_0_1(82)) +
          (occ_func_0_0(0) * occ_func_0_1(2) * occ_func_0_1(3) *
           occ_func_0_1(79)) +
          (occ_func_0_0(0) * occ_func_0_1(11) * occ_func_0_1(6) *
           occ_func_0_1(84)) +
          (occ_func_0_0(0) * occ_func_0_1(5) * occ_func_0_1(7) *
           occ_func_0_1(81)) +
          (occ_func_0_0(0) * occ_func_0_1(1) * occ_func_0_1(6) *
           occ_func_0_1(80)) +
          (occ_func_0_0(0) * occ_func_0_1(4) * occ_func_0_1(10) *
           occ_func_0_1(83)) +
          (occ_func_0_0(0) * occ_func_0_1(4) * occ_func_0_1(1) *
           occ_func_0_1(80))) /
         24.0;
}
double test_Clexulator::eval_bfunc_4_0_11() const {
  return ((occ_func_0_1(0) * occ_func_0_1(9) * occ_func_0_1(7) *
           occ_func_0_1(85)) +
          (occ_func_0_1(0) * occ_func_0_1(3) * occ_func_0_1(9) *
           occ_func_0_1(82)) +
          (occ_func_0_1(0) * occ_func_0_1(2) * occ_func_0_1(5) *
           occ_func_0_1(81)) +
          (occ_func_0_1(0) * occ_func_0_1(12) * occ_func_0_1(10) *
           occ_func_0_1(86)) +
          (occ_func_0_1(0) * occ_func_0_1(8) * occ_func_0_1(3) *
           occ_func_0_1(82)) +
          (occ_func_0_1(0) * occ_func_0_1(11) * occ_func_0_1(12) *
           occ_func_0_1(86)) +
          (occ_func_0_1(0) * occ_func_0_1(7) * occ_func_0_1(2) *
           occ_func_0_1(81)) +
          (occ_func_0_1(0) * occ_func_0_1(6) * occ_func_0_1(8) *
           occ_func_0_1(84)) +
          (occ_func_0_1(0) * occ_func_0_1(8) * occ_func_0_1(11) *
           occ_func_0_1(84)) +
          (occ_func_0_1(0) * occ_func_0_1(12) * occ_func_0_1(9) *
           occ_func_0_1(85)) +
          (occ_func_0_1(0) * occ_func_0_1(5) * occ_func_0_1(4) *
           occ_func_0_1(83)) +
          (occ_func_0_1(0) * occ_func_0_1(7) * occ_func_0_1(12) *
           occ_func_0_1(85)) +
          (occ_func_0_1(0) * occ_func_0_1(1) * occ_func_0_1(2) *
           occ_func_0_1(79)) +
          (occ_func_0_1(0) * occ_func_0_1(3) * occ_func_0_1(1) *
           occ_func_0_1(79)) +
          (occ_func_0_1(0) * occ_func_0_1(10) * occ_func_0_1(5) *
           occ_func_0_1(83)) +
          (occ_func_0_1(0) * occ_func_0_1(10) * occ_func_0_1(11) *
           occ_func_0_1(86)) +
          (occ_func_0_1(0) * occ_func_0_1(6) * occ_func_0_1(4) *
           occ_func_0_1(80)) +
          (occ_func_0_1(0) * occ_func_0_1(9) * occ_func_0_1(8) *
           occ_func_0_1(82)) +
          (occ_func_0_1(0) * occ_func_0_1(2) * occ_func_0_1(3) *
           occ_func_0_1(79)) +
          (occ_func_0_1(0) * occ_func_0_1(11) * occ_func_0_1(6) *
           occ_func_0_1(84)) +
          (occ_func_0_1(0) * occ_func_0_1(5) * occ_func_0_1(7) *
           occ_func_0_1(81)) +
          (occ_func_0_1(0) * occ_func_0_1(1) * occ_func_0_1(6) *
           occ_func_0_1(80)) +
          (occ_func_0_1(0) * occ_func_0_1(4) * occ_func_0_1(10) *
           occ_func_0_1(83)) +
          (occ_func_0_1(0) * occ_func_0_1(4) * occ_func_0_1(1) *
           occ_func_0_1(80))) /
         24.0;
}

double test_Clexulator::site_eval_at_0_bfunc_4_0_0() const {
  return ((occ_func_0_0(0) * occ_func_0_0(9) * occ_func_0_0(7) *
           occ_func_0_0(85)) +
          (occ_func_0_0(4) * occ_func_0_0(0) * occ_func_0_0(5) *
           occ_func_0_0(37)) +
          (occ_func_0_0(6) * occ_func_0_0(8) * occ_func_0_0(0) *
           occ_func_0_0(39)) +
          (occ_func_0_0(80) * occ_func_0_0(24) * occ_func_0_0(22) *
           occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(3) * occ_func_0_0(9) *
           occ_func_0_0(82)) +
          (occ_func_0_0(10) * occ_func_0_0(0) * occ_func_0_0(12) *
           occ_func_0_0(33)) +
          (occ_func_0_0(4) * occ_func_0_0(1) * occ_func_0_0(0) *
           occ_func_0_0(26)) +
          (occ_func_0_0(83) * occ_func_0_0(28) * occ_func_0_0(35) *
           occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(2) * occ_func_0_0(5) *
           occ_func_0_0(81)) +
          (occ_func_0_0(11) * occ_func_0_0(0) * occ_func_0_0(10) *
           occ_func_0_0(30)) +
          (occ_func_0_0(8) * occ_func_0_0(3) * occ_func_0_0(0) *
           occ_func_0_0(25)) +
          (occ_func_0_0(84) * occ_func_0_0(31) * occ_func_0_0(36) *
           occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(12) * occ_func_0_0(10) *
           occ_func_0_0(86)) +
          (occ_func_0_0(1) * occ_func_0_0(0) * occ_func_0_0(4) *
           occ_func_0_0(40)) +
          (occ_func_0_0(3) * occ_func_0_0(9) * occ_func_0_0(0) *
           occ_func_0_0(42)) +
          (occ_func_0_0(79) * occ_func_0_0(21) * occ_func_0_0(19) *
           occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(8) * occ_func_0_0(3) *
           occ_func_0_0(82)) +
          (occ_func_0_0(5) * occ_func_0_0(0) * occ_func_0_0(2) *
           occ_func_0_0(27)) +
          (occ_func_0_0(10) * occ_func_0_0(11) * occ_func_0_0(0) *
           occ_func_0_0(33)) +
          (occ_func_0_0(83) * occ_func_0_0(34) * occ_func_0_0(28) *
           occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(11) * occ_func_0_0(12) *
           occ_func_0_0(86)) +
          (occ_func_0_0(2) * occ_func_0_0(0) * occ_func_0_0(7) *
           occ_func_0_0(41)) +
          (occ_func_0_0(1) * occ_func_0_0(6) * occ_func_0_0(0) *
           occ_func_0_0(40)) +
          (occ_func_0_0(79) * occ_func_0_0(20) * occ_func_0_0(21) *
           occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(7) * occ_func_0_0(2) *
           occ_func_0_0(81)) +
          (occ_func_0_0(6) * occ_func_0_0(0) * occ_func_0_0(1) *
           occ_func_0_0(23)) +
          (occ_func_0_0(11) * occ_func_0_0(12) * occ_func_0_0(0) *
           occ_func_0_0(30)) +
          (occ_func_0_0(84) * occ_func_0_0(38) * occ_func_0_0(31) *
           occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(6) * occ_func_0_0(8) *
           occ_func_0_0(84)) +
          (occ_func_0_0(7) * occ_func_0_0(0) * occ_func_0_0(9) *
           occ_func_0_0(38)) +
          (occ_func_0_0(5) * occ_func_0_0(4) * occ_func_0_0(0) *
           occ_func_0_0(36)) +
          (occ_func_0_0(81) * occ_func_0_0(23) * occ_func_0_0(25) *
           occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(8) * occ_func_0_0(11) *
           occ_func_0_0(84)) +
          (occ_func_0_0(5) * occ_func_0_0(0) * occ_func_0_0(10) *
           occ_func_0_0(36)) +
          (occ_func_0_0(2) * occ_func_0_0(3) * occ_func_0_0(0) *
           occ_func_0_0(31)) +
          (occ_func_0_0(81) * occ_func_0_0(25) * occ_func_0_0(30) *
           occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(12) * occ_func_0_0(9) *
           occ_func_0_0(85)) +
          (occ_func_0_0(1) * occ_func_0_0(0) * occ_func_0_0(3) *
           occ_func_0_0(32)) +
          (occ_func_0_0(4) * occ_func_0_0(10) * occ_func_0_0(0) *
           occ_func_0_0(37)) +
          (occ_func_0_0(80) * occ_func_0_0(29) * occ_func_0_0(24) *
           occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(5) * occ_func_0_0(4) *
           occ_func_0_0(83)) +
          (occ_func_0_0(8) * occ_func_0_0(0) * occ_func_0_0(6) *
           occ_func_0_0(34)) +
          (occ_func_0_0(9) * occ_func_0_0(7) * occ_func_0_0(0) *
           occ_func_0_0(35)) +
          (occ_func_0_0(82) * occ_func_0_0(27) * occ_func_0_0(26) *
           occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(7) * occ_func_0_0(12) *
           occ_func_0_0(85)) +
          (occ_func_0_0(6) * occ_func_0_0(0) * occ_func_0_0(11) *
           occ_func_0_0(39)) +
          (occ_func_0_0(1) * occ_func_0_0(2) * occ_func_0_0(0) *
           occ_func_0_0(32)) +
          (occ_func_0_0(80) * occ_func_0_0(22) * occ_func_0_0(29) *
           occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(1) * occ_func_0_0(2) *
           occ_func_0_0(79)) +
          (occ_func_0_0(12) * occ_func_0_0(0) * occ_func_0_0(7) *
           occ_func_0_0(21)) +
          (occ_func_0_0(11) * occ_func_0_0(6) * occ_func_0_0(0) *
           occ_func_0_0(20)) +
          (occ_func_0_0(86) * occ_func_0_0(40) * occ_func_0_0(41) *
           occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(3) * occ_func_0_0(1) *
           occ_func_0_0(79)) +
          (occ_func_0_0(10) * occ_func_0_0(0) * occ_func_0_0(4) *
           occ_func_0_0(19)) +
          (occ_func_0_0(12) * occ_func_0_0(9) * occ_func_0_0(0) *
           occ_func_0_0(21)) +
          (occ_func_0_0(86) * occ_func_0_0(42) * occ_func_0_0(40) *
           occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(10) * occ_func_0_0(5) *
           occ_func_0_0(83)) +
          (occ_func_0_0(3) * occ_func_0_0(0) * occ_func_0_0(2) *
           occ_func_0_0(28)) +
          (occ_func_0_0(8) * occ_func_0_0(11) * occ_func_0_0(0) *
           occ_func_0_0(34)) +
          (occ_func_0_0(82) * occ_func_0_0(33) * occ_func_0_0(27) *
           occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(10) * occ_func_0_0(11) *
           occ_func_0_0(86)) +
          (occ_func_0_0(3) * occ_func_0_0(0) * occ_func_0_0(8) *
           occ_func_0_0(42)) +
          (occ_func_0_0(2) * occ_func_0_0(5) * occ_func_0_0(0) *
           occ_func_0_0(41)) +
          (occ_func_0_0(79) * occ_func_0_0(19) * occ_func_0_0(20) *
           occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(6) * occ_func_0_0(4) *
           occ_func_0_0(80)) +
          (occ_func_0_0(7) * occ_func_0_0(0) * occ_func_0_0(5) *
           occ_func_0_0(22)) +
          (occ_func_0_0(9) * occ_func_0_0(8) * occ_func_0_0(0) *
           occ_func_0_0(24)) +
          (occ_func_0_0(85) * occ_func_0_0(39) * occ_func_0_0(37) *
           occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(9) * occ_func_0_0(8) *
           occ_func_0_0(82)) +
          (occ_func_0_0(4) * occ_func_0_0(0) * occ_func_0_0(6) *
           occ_func_0_0(26)) +
          (occ_func_0_0(5) * occ_func_0_0(7) * occ_func_0_0(0) *
           occ_func_0_0(27)) +
          (occ_func_0_0(83) * occ_func_0_0(35) * occ_func_0_0(34) *
           occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(2) * occ_func_0_0(3) *
           occ_func_0_0(79)) +
          (occ_func_0_0(11) * occ_func_0_0(0) * occ_func_0_0(8) *
           occ_func_0_0(20)) +
          (occ_func_0_0(10) * occ_func_0_0(5) * occ_func_0_0(0) *
           occ_func_0_0(19)) +
          (occ_func_0_0(86) * occ_func_0_0(41) * occ_func_0_0(42) *
           occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(11) * occ_func_0_0(6) *
           occ_func_0_0(84)) +
          (occ_func_0_0(2) * occ_func_0_0(0) * occ_func_0_0(1) *
           occ_func_0_0(31)) +
          (occ_func_0_0(7) * occ_func_0_0(12) * occ_func_0_0(0) *
           occ_func_0_0(38)) +
          (occ_func_0_0(81) * occ_func_0_0(30) * occ_func_0_0(23) *
           occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(5) * occ_func_0_0(7) *
           occ_func_0_0(81)) +
          (occ_func_0_0(8) * occ_func_0_0(0) * occ_func_0_0(9) *
           occ_func_0_0(25)) +
          (occ_func_0_0(6) * occ_func_0_0(4) * occ_func_0_0(0) *
           occ_func_0_0(23)) +
          (occ_func_0_0(84) * occ_func_0_0(36) * occ_func_0_0(38) *
           occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(1) * occ_func_0_0(6) *
           occ_func_0_0(80)) +
          (occ_func_0_0(12) * occ_func_0_0(0) * occ_func_0_0(11) *
           occ_func_0_0(29)) +
          (occ_func_0_0(7) * occ_func_0_0(2) * occ_func_0_0(0) *
           occ_func_0_0(22)) +
          (occ_func_0_0(85) * occ_func_0_0(32) * occ_func_0_0(39) *
           occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(4) * occ_func_0_0(10) *
           occ_func_0_0(83)) +
          (occ_func_0_0(9) * occ_func_0_0(0) * occ_func_0_0(12) *
           occ_func_0_0(35)) +
          (occ_func_0_0(3) * occ_func_0_0(1) * occ_func_0_0(0) *
           occ_func_0_0(28)) +
          (occ_func_0_0(82) * occ_func_0_0(26) * occ_func_0_0(33) *
           occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(4) * occ_func_0_0(1) *
           occ_func_0_0(80)) +
          (occ_func_0_0(9) * occ_func_0_0(0) * occ_func_0_0(3) *
           occ_func_0_0(24)) +
          (occ_func_0_0(12) * occ_func_0_0(10) * occ_func_0_0(0) *
           occ_func_0_0(29)) +
          (occ_func_0_0(85) * occ_func_0_0(37) * occ_func_0_0(32) *
           occ_func_0_0(0))) /
         24.0;
}
double test_Clexulator::site_eval_at_0_bfunc_4_0_1() const {
  return ((occ_func_0_1(0) * occ_func_0_0(9) * occ_func_0_0(7) *
           occ_func_0_0(85)) +
          (occ_func_0_1(4) * occ_func_0_0(0) * occ_func_0_0(5) *
           occ_func_0_0(37)) +
          (occ_func_0_1(6) * occ_func_0_0(8) * occ_func_0_0(0) *
           occ_func_0_0(39)) +
          (occ_func_0_1(80) * occ_func_0_0(24) * occ_func_0_0(22) *
           occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_0(3) * occ_func_0_0(9) *
           occ_func_0_0(82)) +
          (occ_func_0_1(10) * occ_func_0_0(0) * occ_func_0_0(12) *
           occ_func_0_0(33)) +
          (occ_func_0_1(4) * occ_func_0_0(1) * occ_func_0_0(0) *
           occ_func_0_0(26)) +
          (occ_func_0_1(83) * occ_func_0_0(28) * occ_func_0_0(35) *
           occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_0(2) * occ_func_0_0(5) *
           occ_func_0_0(81)) +
          (occ_func_0_1(11) * occ_func_0_0(0) * occ_func_0_0(10) *
           occ_func_0_0(30)) +
          (occ_func_0_1(8) * occ_func_0_0(3) * occ_func_0_0(0) *
           occ_func_0_0(25)) +
          (occ_func_0_1(84) * occ_func_0_0(31) * occ_func_0_0(36) *
           occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_0(12) * occ_func_0_0(10) *
           occ_func_0_0(86)) +
          (occ_func_0_1(1) * occ_func_0_0(0) * occ_func_0_0(4) *
           occ_func_0_0(40)) +
          (occ_func_0_1(3) * occ_func_0_0(9) * occ_func_0_0(0) *
           occ_func_0_0(42)) +
          (occ_func_0_1(79) * occ_func_0_0(21) * occ_func_0_0(19) *
           occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_0(8) * occ_func_0_0(3) *
           occ_func_0_0(82)) +
          (occ_func_0_1(5) * occ_func_0_0(0) * occ_func_0_0(2) *
           occ_func_0_0(27)) +
          (occ_func_0_1(10) * occ_func_0_0(11) * occ_func_0_0(0) *
           occ_func_0_0(33)) +
          (occ_func_0_1(83) * occ_func_0_0(34) * occ_func_0_0(28) *
           occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_0(11) * occ_func_0_0(12) *
           occ_func_0_0(86)) +
          (occ_func_0_1(2) * occ_func_0_0(0) * occ_func_0_0(7) *
           occ_func_0_0(41)) +
          (occ_func_0_1(1) * occ_func_0_0(6) * occ_func_0_0(0) *
           occ_func_0_0(40)) +
          (occ_func_0_1(79) * occ_func_0_0(20) * occ_func_0_0(21) *
           occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_0(7) * occ_func_0_0(2) *
           occ_func_0_0(81)) +
          (occ_func_0_1(6) * occ_func_0_0(0) * occ_func_0_0(1) *
           occ_func_0_0(23)) +
          (occ_func_0_1(11) * occ_func_0_0(12) * occ_func_0_0(0) *
           occ_func_0_0(30)) +
          (occ_func_0_1(84) * occ_func_0_0(38) * occ_func_0_0(31) *
           occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_0(6) * occ_func_0_0(8) *
           occ_func_0_0(84)) +
          (occ_func_0_1(7) * occ_func_0_0(0) * occ_func_0_0(9) *
           occ_func_0_0(38)) +
          (occ_func_0_1(5) * occ_func_0_0(4) * occ_func_0_0(0) *
           occ_func_0_0(36)) +
          (occ_func_0_1(81) * occ_func_0_0(23) * occ_func_0_0(25) *
           occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_0(8) * occ_func_0_0(11) *
           occ_func_0_0(84)) +
          (occ_func_0_1(5) * occ_func_0_0(0) * occ_func_0_0(10) *
           occ_func_0_0(36)) +
          (occ_func_0_1(2) * occ_func_0_0(3) * occ_func_0_0(0) *
           occ_func_0_0(31)) +
          (occ_func_0_1(81) * occ_func_0_0(25) * occ_func_0_0(30) *
           occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_0(12) * occ_func_0_0(9) *
           occ_func_0_0(85)) +
          (occ_func_0_1(1) * occ_func_0_0(0) * occ_func_0_0(3) *
           occ_func_0_0(32)) +
          (occ_func_0_1(4) * occ_func_0_0(10) * occ_func_0_0(0) *
           occ_func_0_0(37)) +
          (occ_func_0_1(80) * occ_func_0_0(29) * occ_func_0_0(24) *
           occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_0(5) * occ_func_0_0(4) *
           occ_func_0_0(83)) +
          (occ_func_0_1(8) * occ_func_0_0(0) * occ_func_0_0(6) *
           occ_func_0_0(34)) +
          (occ_func_0_1(9) * occ_func_0_0(7) * occ_func_0_0(0) *
           occ_func_0_0(35)) +
          (occ_func_0_1(82) * occ_func_0_0(27) * occ_func_0_0(26) *
           occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_0(7) * occ_func_0_0(12) *
           occ_func_0_0(85)) +
          (occ_func_0_1(6) * occ_func_0_0(0) * occ_func_0_0(11) *
           occ_func_0_0(39)) +
          (occ_func_0_1(1) * occ_func_0_0(2) * occ_func_0_0(0) *
           occ_func_0_0(32)) +
          (occ_func_0_1(80) * occ_func_0_0(22) * occ_func_0_0(29) *
           occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_0(1) * occ_func_0_0(2) *
           occ_func_0_0(79)) +
          (occ_func_0_1(12) * occ_func_0_0(0) * occ_func_0_0(7) *
           occ_func_0_0(21)) +
          (occ_func_0_1(11) * occ_func_0_0(6) * occ_func_0_0(0) *
           occ_func_0_0(20)) +
          (occ_func_0_1(86) * occ_func_0_0(40) * occ_func_0_0(41) *
           occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_0(3) * occ_func_0_0(1) *
           occ_func_0_0(79)) +
          (occ_func_0_1(10) * occ_func_0_0(0) * occ_func_0_0(4) *
           occ_func_0_0(19)) +
          (occ_func_0_1(12) * occ_func_0_0(9) * occ_func_0_0(0) *
           occ_func_0_0(21)) +
          (occ_func_0_1(86) * occ_func_0_0(42) * occ_func_0_0(40) *
           occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_0(10) * occ_func_0_0(5) *
           occ_func_0_0(83)) +
          (occ_func_0_1(3) * occ_func_0_0(0) * occ_func_0_0(2) *
           occ_func_0_0(28)) +
          (occ_func_0_1(8) * occ_func_0_0(11) * occ_func_0_0(0) *
           occ_func_0_0(34)) +
          (occ_func_0_1(82) * occ_func_0_0(33) * occ_func_0_0(27) *
           occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_0(10) * occ_func_0_0(11) *
           occ_func_0_0(86)) +
          (occ_func_0_1(3) * occ_func_0_0(0) * occ_func_0_0(8) *
           occ_func_0_0(42)) +
          (occ_func_0_1(2) * occ_func_0_0(5) * occ_func_0_0(0) *
           occ_func_0_0(41)) +
          (occ_func_0_1(79) * occ_func_0_0(19) * occ_func_0_0(20) *
           occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_0(6) * occ_func_0_0(4) *
           occ_func_0_0(80)) +
          (occ_func_0_1(7) * occ_func_0_0(0) * occ_func_0_0(5) *
           occ_func_0_0(22)) +
          (occ_func_0_1(9) * occ_func_0_0(8) * occ_func_0_0(0) *
           occ_func_0_0(24)) +
          (occ_func_0_1(85) * occ_func_0_0(39) * occ_func_0_0(37) *
           occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_0(9) * occ_func_0_0(8) *
           occ_func_0_0(82)) +
          (occ_func_0_1(4) * occ_func_0_0(0) * occ_func_0_0(6) *
           occ_func_0_0(26)) +
          (occ_func_0_1(5) * occ_func_0_0(7) * occ_func_0_0(0) *
           occ_func_0_0(27)) +
          (occ_func_0_1(83) * occ_func_0_0(35) * occ_func_0_0(34) *
           occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_0(2) * occ_func_0_0(3) *
           occ_func_0_0(79)) +
          (occ_func_0_1(11) * occ_func_0_0(0) * occ_func_0_0(8) *
           occ_func_0_0(20)) +
          (occ_func_0_1(10) * occ_func_0_0(5) * occ_func_0_0(0) *
           occ_func_0_0(19)) +
          (occ_func_0_1(86) * occ_func_0_0(41) * occ_func_0_0(42) *
           occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_0(11) * occ_func_0_0(6) *
           occ_func_0_0(84)) +
          (occ_func_0_1(2) * occ_func_0_0(0) * occ_func_0_0(1) *
           occ_func_0_0(31)) +
          (occ_func_0_1(7) * occ_func_0_0(12) * occ_func_0_0(0) *
           occ_func_0_0(38)) +
          (occ_func_0_1(81) * occ_func_0_0(30) * occ_func_0_0(23) *
           occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_0(5) * occ_func_0_0(7) *
           occ_func_0_0(81)) +
          (occ_func_0_1(8) * occ_func_0_0(0) * occ_func_0_0(9) *
           occ_func_0_0(25)) +
          (occ_func_0_1(6) * occ_func_0_0(4) * occ_func_0_0(0) *
           occ_func_0_0(23)) +
          (occ_func_0_1(84) * occ_func_0_0(36) * occ_func_0_0(38) *
           occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_0(1) * occ_func_0_0(6) *
           occ_func_0_0(80)) +
          (occ_func_0_1(12) * occ_func_0_0(0) * occ_func_0_0(11) *
           occ_func_0_0(29)) +
          (occ_func_0_1(7) * occ_func_0_0(2) * occ_func_0_0(0) *
           occ_func_0_0(22)) +
          (occ_func_0_1(85) * occ_func_0_0(32) * occ_func_0_0(39) *
           occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_0(4) * occ_func_0_0(10) *
           occ_func_0_0(83)) +
          (occ_func_0_1(9) * occ_func_0_0(0) * occ_func_0_0(12) *
           occ_func_0_0(35)) +
          (occ_func_0_1(3) * occ_func_0_0(1) * occ_func_0_0(0) *
           occ_func_0_0(28)) +
          (occ_func_0_1(82) * occ_func_0_0(26) * occ_func_0_0(33) *
           occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_0(4) * occ_func_0_0(1) *
           occ_func_0_0(80)) +
          (occ_func_0_1(9) * occ_func_0_0(0) * occ_func_0_0(3) *
           occ_func_0_0(24)) +
          (occ_func_0_1(12) * occ_func_0_0(10) * occ_func_0_0(0) *
           occ_func_0_0(29)) +
          (occ_func_0_1(85) * occ_func_0_0(37) * occ_func_0_0(32) *
           occ_func_0_0(0))) /
         24.0;
}
double test_Clexulator::site_eval_at_0_bfunc_4_0_2() const {
  return (((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(9) * occ_func_0_0(7) *
                occ_func_0_0(85) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(9) * occ_func_0_1(7) *
                occ_func_0_0(85))) +
          ((0.7071067812 * occ_func_0_0(4) * occ_func_0_1(0) * occ_func_0_0(5) *
                occ_func_0_0(37) +
            0.7071067812 * occ_func_0_0(4) * occ_func_0_0(0) * occ_func_0_1(5) *
                occ_func_0_0(37))) +
          ((0.7071067812 * occ_func_0_0(6) * occ_func_0_1(8) * occ_func_0_0(0) *
                occ_func_0_0(39) +
            0.7071067812 * occ_func_0_0(6) * occ_func_0_0(8) * occ_func_0_1(0) *
                occ_func_0_0(39))) +
          ((0.7071067812 * occ_func_0_0(80) * occ_func_0_1(24) *
                occ_func_0_0(22) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(80) * occ_func_0_0(24) *
                occ_func_0_1(22) * occ_func_0_0(0))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(3) * occ_func_0_0(9) *
                occ_func_0_0(82) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(3) * occ_func_0_1(9) *
                occ_func_0_0(82))) +
          ((0.7071067812 * occ_func_0_0(10) * occ_func_0_1(0) *
                occ_func_0_0(12) * occ_func_0_0(33) +
            0.7071067812 * occ_func_0_0(10) * occ_func_0_0(0) *
                occ_func_0_1(12) * occ_func_0_0(33))) +
          ((0.7071067812 * occ_func_0_0(4) * occ_func_0_1(1) * occ_func_0_0(0) *
                occ_func_0_0(26) +
            0.7071067812 * occ_func_0_0(4) * occ_func_0_0(1) * occ_func_0_1(0) *
                occ_func_0_0(26))) +
          ((0.7071067812 * occ_func_0_0(83) * occ_func_0_1(28) *
                occ_func_0_0(35) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(83) * occ_func_0_0(28) *
                occ_func_0_1(35) * occ_func_0_0(0))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(2) * occ_func_0_0(5) *
                occ_func_0_0(81) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(2) * occ_func_0_1(5) *
                occ_func_0_0(81))) +
          ((0.7071067812 * occ_func_0_0(11) * occ_func_0_1(0) *
                occ_func_0_0(10) * occ_func_0_0(30) +
            0.7071067812 * occ_func_0_0(11) * occ_func_0_0(0) *
                occ_func_0_1(10) * occ_func_0_0(30))) +
          ((0.7071067812 * occ_func_0_0(8) * occ_func_0_1(3) * occ_func_0_0(0) *
                occ_func_0_0(25) +
            0.7071067812 * occ_func_0_0(8) * occ_func_0_0(3) * occ_func_0_1(0) *
                occ_func_0_0(25))) +
          ((0.7071067812 * occ_func_0_0(84) * occ_func_0_1(31) *
                occ_func_0_0(36) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(84) * occ_func_0_0(31) *
                occ_func_0_1(36) * occ_func_0_0(0))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(12) *
                occ_func_0_0(10) * occ_func_0_0(86) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(12) *
                occ_func_0_1(10) * occ_func_0_0(86))) +
          ((0.7071067812 * occ_func_0_0(1) * occ_func_0_1(0) * occ_func_0_0(4) *
                occ_func_0_0(40) +
            0.7071067812 * occ_func_0_0(1) * occ_func_0_0(0) * occ_func_0_1(4) *
                occ_func_0_0(40))) +
          ((0.7071067812 * occ_func_0_0(3) * occ_func_0_1(9) * occ_func_0_0(0) *
                occ_func_0_0(42) +
            0.7071067812 * occ_func_0_0(3) * occ_func_0_0(9) * occ_func_0_1(0) *
                occ_func_0_0(42))) +
          ((0.7071067812 * occ_func_0_0(79) * occ_func_0_1(21) *
                occ_func_0_0(19) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(79) * occ_func_0_0(21) *
                occ_func_0_1(19) * occ_func_0_0(0))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(8) * occ_func_0_0(3) *
                occ_func_0_0(82) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(8) * occ_func_0_1(3) *
                occ_func_0_0(82))) +
          ((0.7071067812 * occ_func_0_0(5) * occ_func_0_1(0) * occ_func_0_0(2) *
                occ_func_0_0(27) +
            0.7071067812 * occ_func_0_0(5) * occ_func_0_0(0) * occ_func_0_1(2) *
                occ_func_0_0(27))) +
          ((0.7071067812 * occ_func_0_0(10) * occ_func_0_1(11) *
                occ_func_0_0(0) * occ_func_0_0(33) +
            0.7071067812 * occ_func_0_0(10) * occ_func_0_0(11) *
                occ_func_0_1(0) * occ_func_0_0(33))) +
          ((0.7071067812 * occ_func_0_0(83) * occ_func_0_1(34) *
                occ_func_0_0(28) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(83) * occ_func_0_0(34) *
                occ_func_0_1(28) * occ_func_0_0(0))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(11) *
                occ_func_0_0(12) * occ_func_0_0(86) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(11) *
                occ_func_0_1(12) * occ_func_0_0(86))) +
          ((0.7071067812 * occ_func_0_0(2) * occ_func_0_1(0) * occ_func_0_0(7) *
                occ_func_0_0(41) +
            0.7071067812 * occ_func_0_0(2) * occ_func_0_0(0) * occ_func_0_1(7) *
                occ_func_0_0(41))) +
          ((0.7071067812 * occ_func_0_0(1) * occ_func_0_1(6) * occ_func_0_0(0) *
                occ_func_0_0(40) +
            0.7071067812 * occ_func_0_0(1) * occ_func_0_0(6) * occ_func_0_1(0) *
                occ_func_0_0(40))) +
          ((0.7071067812 * occ_func_0_0(79) * occ_func_0_1(20) *
                occ_func_0_0(21) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(79) * occ_func_0_0(20) *
                occ_func_0_1(21) * occ_func_0_0(0))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(7) * occ_func_0_0(2) *
                occ_func_0_0(81) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(7) * occ_func_0_1(2) *
                occ_func_0_0(81))) +
          ((0.7071067812 * occ_func_0_0(6) * occ_func_0_1(0) * occ_func_0_0(1) *
                occ_func_0_0(23) +
            0.7071067812 * occ_func_0_0(6) * occ_func_0_0(0) * occ_func_0_1(1) *
                occ_func_0_0(23))) +
          ((0.7071067812 * occ_func_0_0(11) * occ_func_0_1(12) *
                occ_func_0_0(0) * occ_func_0_0(30) +
            0.7071067812 * occ_func_0_0(11) * occ_func_0_0(12) *
                occ_func_0_1(0) * occ_func_0_0(30))) +
          ((0.7071067812 * occ_func_0_0(84) * occ_func_0_1(38) *
                occ_func_0_0(31) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(84) * occ_func_0_0(38) *
                occ_func_0_1(31) * occ_func_0_0(0))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(6) * occ_func_0_0(8) *
                occ_func_0_0(84) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(6) * occ_func_0_1(8) *
                occ_func_0_0(84))) +
          ((0.7071067812 * occ_func_0_0(7) * occ_func_0_1(0) * occ_func_0_0(9) *
                occ_func_0_0(38) +
            0.7071067812 * occ_func_0_0(7) * occ_func_0_0(0) * occ_func_0_1(9) *
                occ_func_0_0(38))) +
          ((0.7071067812 * occ_func_0_0(5) * occ_func_0_1(4) * occ_func_0_0(0) *
                occ_func_0_0(36) +
            0.7071067812 * occ_func_0_0(5) * occ_func_0_0(4) * occ_func_0_1(0) *
                occ_func_0_0(36))) +
          ((0.7071067812 * occ_func_0_0(81) * occ_func_0_1(23) *
                occ_func_0_0(25) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(81) * occ_func_0_0(23) *
                occ_func_0_1(25) * occ_func_0_0(0))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(8) *
                occ_func_0_0(11) * occ_func_0_0(84) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(8) *
                occ_func_0_1(11) * occ_func_0_0(84))) +
          ((0.7071067812 * occ_func_0_0(5) * occ_func_0_1(0) *
                occ_func_0_0(10) * occ_func_0_0(36) +
            0.7071067812 * occ_func_0_0(5) * occ_func_0_0(0) *
                occ_func_0_1(10) * occ_func_0_0(36))) +
          ((0.7071067812 * occ_func_0_0(2) * occ_func_0_1(3) * occ_func_0_0(0) *
                occ_func_0_0(31) +
            0.7071067812 * occ_func_0_0(2) * occ_func_0_0(3) * occ_func_0_1(0) *
                occ_func_0_0(31))) +
          ((0.7071067812 * occ_func_0_0(81) * occ_func_0_1(25) *
                occ_func_0_0(30) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(81) * occ_func_0_0(25) *
                occ_func_0_1(30) * occ_func_0_0(0))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(12) *
                occ_func_0_0(9) * occ_func_0_0(85) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(12) *
                occ_func_0_1(9) * occ_func_0_0(85))) +
          ((0.7071067812 * occ_func_0_0(1) * occ_func_0_1(0) * occ_func_0_0(3) *
                occ_func_0_0(32) +
            0.7071067812 * occ_func_0_0(1) * occ_func_0_0(0) * occ_func_0_1(3) *
                occ_func_0_0(32))) +
          ((0.7071067812 * occ_func_0_0(4) * occ_func_0_1(10) *
                occ_func_0_0(0) * occ_func_0_0(37) +
            0.7071067812 * occ_func_0_0(4) * occ_func_0_0(10) *
                occ_func_0_1(0) * occ_func_0_0(37))) +
          ((0.7071067812 * occ_func_0_0(80) * occ_func_0_1(29) *
                occ_func_0_0(24) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(80) * occ_func_0_0(29) *
                occ_func_0_1(24) * occ_func_0_0(0))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(5) * occ_func_0_0(4) *
                occ_func_0_0(83) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(5) * occ_func_0_1(4) *
                occ_func_0_0(83))) +
          ((0.7071067812 * occ_func_0_0(8) * occ_func_0_1(0) * occ_func_0_0(6) *
                occ_func_0_0(34) +
            0.7071067812 * occ_func_0_0(8) * occ_func_0_0(0) * occ_func_0_1(6) *
                occ_func_0_0(34))) +
          ((0.7071067812 * occ_func_0_0(9) * occ_func_0_1(7) * occ_func_0_0(0) *
                occ_func_0_0(35) +
            0.7071067812 * occ_func_0_0(9) * occ_func_0_0(7) * occ_func_0_1(0) *
                occ_func_0_0(35))) +
          ((0.7071067812 * occ_func_0_0(82) * occ_func_0_1(27) *
                occ_func_0_0(26) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(82) * occ_func_0_0(27) *
                occ_func_0_1(26) * occ_func_0_0(0))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(7) *
                occ_func_0_0(12) * occ_func_0_0(85) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(7) *
                occ_func_0_1(12) * occ_func_0_0(85))) +
          ((0.7071067812 * occ_func_0_0(6) * occ_func_0_1(0) *
                occ_func_0_0(11) * occ_func_0_0(39) +
            0.7071067812 * occ_func_0_0(6) * occ_func_0_0(0) *
                occ_func_0_1(11) * occ_func_0_0(39))) +
          ((0.7071067812 * occ_func_0_0(1) * occ_func_0_1(2) * occ_func_0_0(0) *
                occ_func_0_0(32) +
            0.7071067812 * occ_func_0_0(1) * occ_func_0_0(2) * occ_func_0_1(0) *
                occ_func_0_0(32))) +
          ((0.7071067812 * occ_func_0_0(80) * occ_func_0_1(22) *
                occ_func_0_0(29) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(80) * occ_func_0_0(22) *
                occ_func_0_1(29) * occ_func_0_0(0))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(1) * occ_func_0_0(2) *
                occ_func_0_0(79) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(1) * occ_func_0_1(2) *
                occ_func_0_0(79))) +
          ((0.7071067812 * occ_func_0_0(12) * occ_func_0_1(0) *
                occ_func_0_0(7) * occ_func_0_0(21) +
            0.7071067812 * occ_func_0_0(12) * occ_func_0_0(0) *
                occ_func_0_1(7) * occ_func_0_0(21))) +
          ((0.7071067812 * occ_func_0_0(11) * occ_func_0_1(6) *
                occ_func_0_0(0) * occ_func_0_0(20) +
            0.7071067812 * occ_func_0_0(11) * occ_func_0_0(6) *
                occ_func_0_1(0) * occ_func_0_0(20))) +
          ((0.7071067812 * occ_func_0_0(86) * occ_func_0_1(40) *
                occ_func_0_0(41) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(86) * occ_func_0_0(40) *
                occ_func_0_1(41) * occ_func_0_0(0))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(3) * occ_func_0_0(1) *
                occ_func_0_0(79) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(3) * occ_func_0_1(1) *
                occ_func_0_0(79))) +
          ((0.7071067812 * occ_func_0_0(10) * occ_func_0_1(0) *
                occ_func_0_0(4) * occ_func_0_0(19) +
            0.7071067812 * occ_func_0_0(10) * occ_func_0_0(0) *
                occ_func_0_1(4) * occ_func_0_0(19))) +
          ((0.7071067812 * occ_func_0_0(12) * occ_func_0_1(9) *
                occ_func_0_0(0) * occ_func_0_0(21) +
            0.7071067812 * occ_func_0_0(12) * occ_func_0_0(9) *
                occ_func_0_1(0) * occ_func_0_0(21))) +
          ((0.7071067812 * occ_func_0_0(86) * occ_func_0_1(42) *
                occ_func_0_0(40) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(86) * occ_func_0_0(42) *
                occ_func_0_1(40) * occ_func_0_0(0))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(10) *
                occ_func_0_0(5) * occ_func_0_0(83) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(10) *
                occ_func_0_1(5) * occ_func_0_0(83))) +
          ((0.7071067812 * occ_func_0_0(3) * occ_func_0_1(0) * occ_func_0_0(2) *
                occ_func_0_0(28) +
            0.7071067812 * occ_func_0_0(3) * occ_func_0_0(0) * occ_func_0_1(2) *
                occ_func_0_0(28))) +
          ((0.7071067812 * occ_func_0_0(8) * occ_func_0_1(11) *
                occ_func_0_0(0) * occ_func_0_0(34) +
            0.7071067812 * occ_func_0_0(8) * occ_func_0_0(11) *
                occ_func_0_1(0) * occ_func_0_0(34))) +
          ((0.7071067812 * occ_func_0_0(82) * occ_func_0_1(33) *
                occ_func_0_0(27) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(82) * occ_func_0_0(33) *
                occ_func_0_1(27) * occ_func_0_0(0))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(10) *
                occ_func_0_0(11) * occ_func_0_0(86) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(10) *
                occ_func_0_1(11) * occ_func_0_0(86))) +
          ((0.7071067812 * occ_func_0_0(3) * occ_func_0_1(0) * occ_func_0_0(8) *
                occ_func_0_0(42) +
            0.7071067812 * occ_func_0_0(3) * occ_func_0_0(0) * occ_func_0_1(8) *
                occ_func_0_0(42))) +
          ((0.7071067812 * occ_func_0_0(2) * occ_func_0_1(5) * occ_func_0_0(0) *
                occ_func_0_0(41) +
            0.7071067812 * occ_func_0_0(2) * occ_func_0_0(5) * occ_func_0_1(0) *
                occ_func_0_0(41))) +
          ((0.7071067812 * occ_func_0_0(79) * occ_func_0_1(19) *
                occ_func_0_0(20) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(79) * occ_func_0_0(19) *
                occ_func_0_1(20) * occ_func_0_0(0))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(6) * occ_func_0_0(4) *
                occ_func_0_0(80) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(6) * occ_func_0_1(4) *
                occ_func_0_0(80))) +
          ((0.7071067812 * occ_func_0_0(7) * occ_func_0_1(0) * occ_func_0_0(5) *
                occ_func_0_0(22) +
            0.7071067812 * occ_func_0_0(7) * occ_func_0_0(0) * occ_func_0_1(5) *
                occ_func_0_0(22))) +
          ((0.7071067812 * occ_func_0_0(9) * occ_func_0_1(8) * occ_func_0_0(0) *
                occ_func_0_0(24) +
            0.7071067812 * occ_func_0_0(9) * occ_func_0_0(8) * occ_func_0_1(0) *
                occ_func_0_0(24))) +
          ((0.7071067812 * occ_func_0_0(85) * occ_func_0_1(39) *
                occ_func_0_0(37) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(85) * occ_func_0_0(39) *
                occ_func_0_1(37) * occ_func_0_0(0))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(9) * occ_func_0_0(8) *
                occ_func_0_0(82) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(9) * occ_func_0_1(8) *
                occ_func_0_0(82))) +
          ((0.7071067812 * occ_func_0_0(4) * occ_func_0_1(0) * occ_func_0_0(6) *
                occ_func_0_0(26) +
            0.7071067812 * occ_func_0_0(4) * occ_func_0_0(0) * occ_func_0_1(6) *
                occ_func_0_0(26))) +
          ((0.7071067812 * occ_func_0_0(5) * occ_func_0_1(7) * occ_func_0_0(0) *
                occ_func_0_0(27) +
            0.7071067812 * occ_func_0_0(5) * occ_func_0_0(7) * occ_func_0_1(0) *
                occ_func_0_0(27))) +
          ((0.7071067812 * occ_func_0_0(83) * occ_func_0_1(35) *
                occ_func_0_0(34) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(83) * occ_func_0_0(35) *
                occ_func_0_1(34) * occ_func_0_0(0))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(2) * occ_func_0_0(3) *
                occ_func_0_0(79) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(2) * occ_func_0_1(3) *
                occ_func_0_0(79))) +
          ((0.7071067812 * occ_func_0_0(11) * occ_func_0_1(0) *
                occ_func_0_0(8) * occ_func_0_0(20) +
            0.7071067812 * occ_func_0_0(11) * occ_func_0_0(0) *
                occ_func_0_1(8) * occ_func_0_0(20))) +
          ((0.7071067812 * occ_func_0_0(10) * occ_func_0_1(5) *
                occ_func_0_0(0) * occ_func_0_0(19) +
            0.7071067812 * occ_func_0_0(10) * occ_func_0_0(5) *
                occ_func_0_1(0) * occ_func_0_0(19))) +
          ((0.7071067812 * occ_func_0_0(86) * occ_func_0_1(41) *
                occ_func_0_0(42) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(86) * occ_func_0_0(41) *
                occ_func_0_1(42) * occ_func_0_0(0))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(11) *
                occ_func_0_0(6) * occ_func_0_0(84) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(11) *
                occ_func_0_1(6) * occ_func_0_0(84))) +
          ((0.7071067812 * occ_func_0_0(2) * occ_func_0_1(0) * occ_func_0_0(1) *
                occ_func_0_0(31) +
            0.7071067812 * occ_func_0_0(2) * occ_func_0_0(0) * occ_func_0_1(1) *
                occ_func_0_0(31))) +
          ((0.7071067812 * occ_func_0_0(7) * occ_func_0_1(12) *
                occ_func_0_0(0) * occ_func_0_0(38) +
            0.7071067812 * occ_func_0_0(7) * occ_func_0_0(12) *
                occ_func_0_1(0) * occ_func_0_0(38))) +
          ((0.7071067812 * occ_func_0_0(81) * occ_func_0_1(30) *
                occ_func_0_0(23) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(81) * occ_func_0_0(30) *
                occ_func_0_1(23) * occ_func_0_0(0))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(5) * occ_func_0_0(7) *
                occ_func_0_0(81) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(5) * occ_func_0_1(7) *
                occ_func_0_0(81))) +
          ((0.7071067812 * occ_func_0_0(8) * occ_func_0_1(0) * occ_func_0_0(9) *
                occ_func_0_0(25) +
            0.7071067812 * occ_func_0_0(8) * occ_func_0_0(0) * occ_func_0_1(9) *
                occ_func_0_0(25))) +
          ((0.7071067812 * occ_func_0_0(6) * occ_func_0_1(4) * occ_func_0_0(0) *
                occ_func_0_0(23) +
            0.7071067812 * occ_func_0_0(6) * occ_func_0_0(4) * occ_func_0_1(0) *
                occ_func_0_0(23))) +
          ((0.7071067812 * occ_func_0_0(84) * occ_func_0_1(36) *
                occ_func_0_0(38) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(84) * occ_func_0_0(36) *
                occ_func_0_1(38) * occ_func_0_0(0))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(1) * occ_func_0_0(6) *
                occ_func_0_0(80) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(1) * occ_func_0_1(6) *
                occ_func_0_0(80))) +
          ((0.7071067812 * occ_func_0_0(12) * occ_func_0_1(0) *
                occ_func_0_0(11) * occ_func_0_0(29) +
            0.7071067812 * occ_func_0_0(12) * occ_func_0_0(0) *
                occ_func_0_1(11) * occ_func_0_0(29))) +
          ((0.7071067812 * occ_func_0_0(7) * occ_func_0_1(2) * occ_func_0_0(0) *
                occ_func_0_0(22) +
            0.7071067812 * occ_func_0_0(7) * occ_func_0_0(2) * occ_func_0_1(0) *
                occ_func_0_0(22))) +
          ((0.7071067812 * occ_func_0_0(85) * occ_func_0_1(32) *
                occ_func_0_0(39) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(85) * occ_func_0_0(32) *
                occ_func_0_1(39) * occ_func_0_0(0))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(4) *
                occ_func_0_0(10) * occ_func_0_0(83) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(4) *
                occ_func_0_1(10) * occ_func_0_0(83))) +
          ((0.7071067812 * occ_func_0_0(9) * occ_func_0_1(0) *
                occ_func_0_0(12) * occ_func_0_0(35) +
            0.7071067812 * occ_func_0_0(9) * occ_func_0_0(0) *
                occ_func_0_1(12) * occ_func_0_0(35))) +
          ((0.7071067812 * occ_func_0_0(3) * occ_func_0_1(1) * occ_func_0_0(0) *
                occ_func_0_0(28) +
            0.7071067812 * occ_func_0_0(3) * occ_func_0_0(1) * occ_func_0_1(0) *
                occ_func_0_0(28))) +
          ((0.7071067812 * occ_func_0_0(82) * occ_func_0_1(26) *
                occ_func_0_0(33) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(82) * occ_func_0_0(26) *
                occ_func_0_1(33) * occ_func_0_0(0))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(4) * occ_func_0_0(1) *
                occ_func_0_0(80) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(4) * occ_func_0_1(1) *
                occ_func_0_0(80))) +
          ((0.7071067812 * occ_func_0_0(9) * occ_func_0_1(0) * occ_func_0_0(3) *
                occ_func_0_0(24) +
            0.7071067812 * occ_func_0_0(9) * occ_func_0_0(0) * occ_func_0_1(3) *
                occ_func_0_0(24))) +
          ((0.7071067812 * occ_func_0_0(12) * occ_func_0_1(10) *
                occ_func_0_0(0) * occ_func_0_0(29) +
            0.7071067812 * occ_func_0_0(12) * occ_func_0_0(10) *
                occ_func_0_1(0) * occ_func_0_0(29))) +
          ((0.7071067812 * occ_func_0_0(85) * occ_func_0_1(37) *
                occ_func_0_0(32) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(85) * occ_func_0_0(37) *
                occ_func_0_1(32) * occ_func_0_0(0)))) /
         24.0;
}
double test_Clexulator::site_eval_at_0_bfunc_4_0_3() const {
  return (((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(9) * occ_func_0_0(7) *
                occ_func_0_0(85) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(9) * occ_func_0_1(7) *
                occ_func_0_0(85))) +
          ((0.7071067812 * occ_func_0_1(4) * occ_func_0_1(0) * occ_func_0_0(5) *
                occ_func_0_0(37) +
            0.7071067812 * occ_func_0_1(4) * occ_func_0_0(0) * occ_func_0_1(5) *
                occ_func_0_0(37))) +
          ((0.7071067812 * occ_func_0_1(6) * occ_func_0_1(8) * occ_func_0_0(0) *
                occ_func_0_0(39) +
            0.7071067812 * occ_func_0_1(6) * occ_func_0_0(8) * occ_func_0_1(0) *
                occ_func_0_0(39))) +
          ((0.7071067812 * occ_func_0_1(80) * occ_func_0_1(24) *
                occ_func_0_0(22) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_1(80) * occ_func_0_0(24) *
                occ_func_0_1(22) * occ_func_0_0(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(3) * occ_func_0_0(9) *
                occ_func_0_0(82) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(3) * occ_func_0_1(9) *
                occ_func_0_0(82))) +
          ((0.7071067812 * occ_func_0_1(10) * occ_func_0_1(0) *
                occ_func_0_0(12) * occ_func_0_0(33) +
            0.7071067812 * occ_func_0_1(10) * occ_func_0_0(0) *
                occ_func_0_1(12) * occ_func_0_0(33))) +
          ((0.7071067812 * occ_func_0_1(4) * occ_func_0_1(1) * occ_func_0_0(0) *
                occ_func_0_0(26) +
            0.7071067812 * occ_func_0_1(4) * occ_func_0_0(1) * occ_func_0_1(0) *
                occ_func_0_0(26))) +
          ((0.7071067812 * occ_func_0_1(83) * occ_func_0_1(28) *
                occ_func_0_0(35) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_1(83) * occ_func_0_0(28) *
                occ_func_0_1(35) * occ_func_0_0(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(2) * occ_func_0_0(5) *
                occ_func_0_0(81) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(2) * occ_func_0_1(5) *
                occ_func_0_0(81))) +
          ((0.7071067812 * occ_func_0_1(11) * occ_func_0_1(0) *
                occ_func_0_0(10) * occ_func_0_0(30) +
            0.7071067812 * occ_func_0_1(11) * occ_func_0_0(0) *
                occ_func_0_1(10) * occ_func_0_0(30))) +
          ((0.7071067812 * occ_func_0_1(8) * occ_func_0_1(3) * occ_func_0_0(0) *
                occ_func_0_0(25) +
            0.7071067812 * occ_func_0_1(8) * occ_func_0_0(3) * occ_func_0_1(0) *
                occ_func_0_0(25))) +
          ((0.7071067812 * occ_func_0_1(84) * occ_func_0_1(31) *
                occ_func_0_0(36) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_1(84) * occ_func_0_0(31) *
                occ_func_0_1(36) * occ_func_0_0(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(12) *
                occ_func_0_0(10) * occ_func_0_0(86) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(12) *
                occ_func_0_1(10) * occ_func_0_0(86))) +
          ((0.7071067812 * occ_func_0_1(1) * occ_func_0_1(0) * occ_func_0_0(4) *
                occ_func_0_0(40) +
            0.7071067812 * occ_func_0_1(1) * occ_func_0_0(0) * occ_func_0_1(4) *
                occ_func_0_0(40))) +
          ((0.7071067812 * occ_func_0_1(3) * occ_func_0_1(9) * occ_func_0_0(0) *
                occ_func_0_0(42) +
            0.7071067812 * occ_func_0_1(3) * occ_func_0_0(9) * occ_func_0_1(0) *
                occ_func_0_0(42))) +
          ((0.7071067812 * occ_func_0_1(79) * occ_func_0_1(21) *
                occ_func_0_0(19) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_1(79) * occ_func_0_0(21) *
                occ_func_0_1(19) * occ_func_0_0(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(8) * occ_func_0_0(3) *
                occ_func_0_0(82) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(8) * occ_func_0_1(3) *
                occ_func_0_0(82))) +
          ((0.7071067812 * occ_func_0_1(5) * occ_func_0_1(0) * occ_func_0_0(2) *
                occ_func_0_0(27) +
            0.7071067812 * occ_func_0_1(5) * occ_func_0_0(0) * occ_func_0_1(2) *
                occ_func_0_0(27))) +
          ((0.7071067812 * occ_func_0_1(10) * occ_func_0_1(11) *
                occ_func_0_0(0) * occ_func_0_0(33) +
            0.7071067812 * occ_func_0_1(10) * occ_func_0_0(11) *
                occ_func_0_1(0) * occ_func_0_0(33))) +
          ((0.7071067812 * occ_func_0_1(83) * occ_func_0_1(34) *
                occ_func_0_0(28) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_1(83) * occ_func_0_0(34) *
                occ_func_0_1(28) * occ_func_0_0(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(11) *
                occ_func_0_0(12) * occ_func_0_0(86) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(11) *
                occ_func_0_1(12) * occ_func_0_0(86))) +
          ((0.7071067812 * occ_func_0_1(2) * occ_func_0_1(0) * occ_func_0_0(7) *
                occ_func_0_0(41) +
            0.7071067812 * occ_func_0_1(2) * occ_func_0_0(0) * occ_func_0_1(7) *
                occ_func_0_0(41))) +
          ((0.7071067812 * occ_func_0_1(1) * occ_func_0_1(6) * occ_func_0_0(0) *
                occ_func_0_0(40) +
            0.7071067812 * occ_func_0_1(1) * occ_func_0_0(6) * occ_func_0_1(0) *
                occ_func_0_0(40))) +
          ((0.7071067812 * occ_func_0_1(79) * occ_func_0_1(20) *
                occ_func_0_0(21) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_1(79) * occ_func_0_0(20) *
                occ_func_0_1(21) * occ_func_0_0(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(7) * occ_func_0_0(2) *
                occ_func_0_0(81) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(7) * occ_func_0_1(2) *
                occ_func_0_0(81))) +
          ((0.7071067812 * occ_func_0_1(6) * occ_func_0_1(0) * occ_func_0_0(1) *
                occ_func_0_0(23) +
            0.7071067812 * occ_func_0_1(6) * occ_func_0_0(0) * occ_func_0_1(1) *
                occ_func_0_0(23))) +
          ((0.7071067812 * occ_func_0_1(11) * occ_func_0_1(12) *
                occ_func_0_0(0) * occ_func_0_0(30) +
            0.7071067812 * occ_func_0_1(11) * occ_func_0_0(12) *
                occ_func_0_1(0) * occ_func_0_0(30))) +
          ((0.7071067812 * occ_func_0_1(84) * occ_func_0_1(38) *
                occ_func_0_0(31) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_1(84) * occ_func_0_0(38) *
                occ_func_0_1(31) * occ_func_0_0(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(6) * occ_func_0_0(8) *
                occ_func_0_0(84) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(6) * occ_func_0_1(8) *
                occ_func_0_0(84))) +
          ((0.7071067812 * occ_func_0_1(7) * occ_func_0_1(0) * occ_func_0_0(9) *
                occ_func_0_0(38) +
            0.7071067812 * occ_func_0_1(7) * occ_func_0_0(0) * occ_func_0_1(9) *
                occ_func_0_0(38))) +
          ((0.7071067812 * occ_func_0_1(5) * occ_func_0_1(4) * occ_func_0_0(0) *
                occ_func_0_0(36) +
            0.7071067812 * occ_func_0_1(5) * occ_func_0_0(4) * occ_func_0_1(0) *
                occ_func_0_0(36))) +
          ((0.7071067812 * occ_func_0_1(81) * occ_func_0_1(23) *
                occ_func_0_0(25) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_1(81) * occ_func_0_0(23) *
                occ_func_0_1(25) * occ_func_0_0(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(8) *
                occ_func_0_0(11) * occ_func_0_0(84) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(8) *
                occ_func_0_1(11) * occ_func_0_0(84))) +
          ((0.7071067812 * occ_func_0_1(5) * occ_func_0_1(0) *
                occ_func_0_0(10) * occ_func_0_0(36) +
            0.7071067812 * occ_func_0_1(5) * occ_func_0_0(0) *
                occ_func_0_1(10) * occ_func_0_0(36))) +
          ((0.7071067812 * occ_func_0_1(2) * occ_func_0_1(3) * occ_func_0_0(0) *
                occ_func_0_0(31) +
            0.7071067812 * occ_func_0_1(2) * occ_func_0_0(3) * occ_func_0_1(0) *
                occ_func_0_0(31))) +
          ((0.7071067812 * occ_func_0_1(81) * occ_func_0_1(25) *
                occ_func_0_0(30) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_1(81) * occ_func_0_0(25) *
                occ_func_0_1(30) * occ_func_0_0(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(12) *
                occ_func_0_0(9) * occ_func_0_0(85) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(12) *
                occ_func_0_1(9) * occ_func_0_0(85))) +
          ((0.7071067812 * occ_func_0_1(1) * occ_func_0_1(0) * occ_func_0_0(3) *
                occ_func_0_0(32) +
            0.7071067812 * occ_func_0_1(1) * occ_func_0_0(0) * occ_func_0_1(3) *
                occ_func_0_0(32))) +
          ((0.7071067812 * occ_func_0_1(4) * occ_func_0_1(10) *
                occ_func_0_0(0) * occ_func_0_0(37) +
            0.7071067812 * occ_func_0_1(4) * occ_func_0_0(10) *
                occ_func_0_1(0) * occ_func_0_0(37))) +
          ((0.7071067812 * occ_func_0_1(80) * occ_func_0_1(29) *
                occ_func_0_0(24) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_1(80) * occ_func_0_0(29) *
                occ_func_0_1(24) * occ_func_0_0(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(5) * occ_func_0_0(4) *
                occ_func_0_0(83) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(5) * occ_func_0_1(4) *
                occ_func_0_0(83))) +
          ((0.7071067812 * occ_func_0_1(8) * occ_func_0_1(0) * occ_func_0_0(6) *
                occ_func_0_0(34) +
            0.7071067812 * occ_func_0_1(8) * occ_func_0_0(0) * occ_func_0_1(6) *
                occ_func_0_0(34))) +
          ((0.7071067812 * occ_func_0_1(9) * occ_func_0_1(7) * occ_func_0_0(0) *
                occ_func_0_0(35) +
            0.7071067812 * occ_func_0_1(9) * occ_func_0_0(7) * occ_func_0_1(0) *
                occ_func_0_0(35))) +
          ((0.7071067812 * occ_func_0_1(82) * occ_func_0_1(27) *
                occ_func_0_0(26) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_1(82) * occ_func_0_0(27) *
                occ_func_0_1(26) * occ_func_0_0(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(7) *
                occ_func_0_0(12) * occ_func_0_0(85) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(7) *
                occ_func_0_1(12) * occ_func_0_0(85))) +
          ((0.7071067812 * occ_func_0_1(6) * occ_func_0_1(0) *
                occ_func_0_0(11) * occ_func_0_0(39) +
            0.7071067812 * occ_func_0_1(6) * occ_func_0_0(0) *
                occ_func_0_1(11) * occ_func_0_0(39))) +
          ((0.7071067812 * occ_func_0_1(1) * occ_func_0_1(2) * occ_func_0_0(0) *
                occ_func_0_0(32) +
            0.7071067812 * occ_func_0_1(1) * occ_func_0_0(2) * occ_func_0_1(0) *
                occ_func_0_0(32))) +
          ((0.7071067812 * occ_func_0_1(80) * occ_func_0_1(22) *
                occ_func_0_0(29) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_1(80) * occ_func_0_0(22) *
                occ_func_0_1(29) * occ_func_0_0(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(1) * occ_func_0_0(2) *
                occ_func_0_0(79) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(1) * occ_func_0_1(2) *
                occ_func_0_0(79))) +
          ((0.7071067812 * occ_func_0_1(12) * occ_func_0_1(0) *
                occ_func_0_0(7) * occ_func_0_0(21) +
            0.7071067812 * occ_func_0_1(12) * occ_func_0_0(0) *
                occ_func_0_1(7) * occ_func_0_0(21))) +
          ((0.7071067812 * occ_func_0_1(11) * occ_func_0_1(6) *
                occ_func_0_0(0) * occ_func_0_0(20) +
            0.7071067812 * occ_func_0_1(11) * occ_func_0_0(6) *
                occ_func_0_1(0) * occ_func_0_0(20))) +
          ((0.7071067812 * occ_func_0_1(86) * occ_func_0_1(40) *
                occ_func_0_0(41) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_1(86) * occ_func_0_0(40) *
                occ_func_0_1(41) * occ_func_0_0(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(3) * occ_func_0_0(1) *
                occ_func_0_0(79) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(3) * occ_func_0_1(1) *
                occ_func_0_0(79))) +
          ((0.7071067812 * occ_func_0_1(10) * occ_func_0_1(0) *
                occ_func_0_0(4) * occ_func_0_0(19) +
            0.7071067812 * occ_func_0_1(10) * occ_func_0_0(0) *
                occ_func_0_1(4) * occ_func_0_0(19))) +
          ((0.7071067812 * occ_func_0_1(12) * occ_func_0_1(9) *
                occ_func_0_0(0) * occ_func_0_0(21) +
            0.7071067812 * occ_func_0_1(12) * occ_func_0_0(9) *
                occ_func_0_1(0) * occ_func_0_0(21))) +
          ((0.7071067812 * occ_func_0_1(86) * occ_func_0_1(42) *
                occ_func_0_0(40) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_1(86) * occ_func_0_0(42) *
                occ_func_0_1(40) * occ_func_0_0(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(10) *
                occ_func_0_0(5) * occ_func_0_0(83) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(10) *
                occ_func_0_1(5) * occ_func_0_0(83))) +
          ((0.7071067812 * occ_func_0_1(3) * occ_func_0_1(0) * occ_func_0_0(2) *
                occ_func_0_0(28) +
            0.7071067812 * occ_func_0_1(3) * occ_func_0_0(0) * occ_func_0_1(2) *
                occ_func_0_0(28))) +
          ((0.7071067812 * occ_func_0_1(8) * occ_func_0_1(11) *
                occ_func_0_0(0) * occ_func_0_0(34) +
            0.7071067812 * occ_func_0_1(8) * occ_func_0_0(11) *
                occ_func_0_1(0) * occ_func_0_0(34))) +
          ((0.7071067812 * occ_func_0_1(82) * occ_func_0_1(33) *
                occ_func_0_0(27) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_1(82) * occ_func_0_0(33) *
                occ_func_0_1(27) * occ_func_0_0(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(10) *
                occ_func_0_0(11) * occ_func_0_0(86) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(10) *
                occ_func_0_1(11) * occ_func_0_0(86))) +
          ((0.7071067812 * occ_func_0_1(3) * occ_func_0_1(0) * occ_func_0_0(8) *
                occ_func_0_0(42) +
            0.7071067812 * occ_func_0_1(3) * occ_func_0_0(0) * occ_func_0_1(8) *
                occ_func_0_0(42))) +
          ((0.7071067812 * occ_func_0_1(2) * occ_func_0_1(5) * occ_func_0_0(0) *
                occ_func_0_0(41) +
            0.7071067812 * occ_func_0_1(2) * occ_func_0_0(5) * occ_func_0_1(0) *
                occ_func_0_0(41))) +
          ((0.7071067812 * occ_func_0_1(79) * occ_func_0_1(19) *
                occ_func_0_0(20) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_1(79) * occ_func_0_0(19) *
                occ_func_0_1(20) * occ_func_0_0(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(6) * occ_func_0_0(4) *
                occ_func_0_0(80) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(6) * occ_func_0_1(4) *
                occ_func_0_0(80))) +
          ((0.7071067812 * occ_func_0_1(7) * occ_func_0_1(0) * occ_func_0_0(5) *
                occ_func_0_0(22) +
            0.7071067812 * occ_func_0_1(7) * occ_func_0_0(0) * occ_func_0_1(5) *
                occ_func_0_0(22))) +
          ((0.7071067812 * occ_func_0_1(9) * occ_func_0_1(8) * occ_func_0_0(0) *
                occ_func_0_0(24) +
            0.7071067812 * occ_func_0_1(9) * occ_func_0_0(8) * occ_func_0_1(0) *
                occ_func_0_0(24))) +
          ((0.7071067812 * occ_func_0_1(85) * occ_func_0_1(39) *
                occ_func_0_0(37) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_1(85) * occ_func_0_0(39) *
                occ_func_0_1(37) * occ_func_0_0(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(9) * occ_func_0_0(8) *
                occ_func_0_0(82) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(9) * occ_func_0_1(8) *
                occ_func_0_0(82))) +
          ((0.7071067812 * occ_func_0_1(4) * occ_func_0_1(0) * occ_func_0_0(6) *
                occ_func_0_0(26) +
            0.7071067812 * occ_func_0_1(4) * occ_func_0_0(0) * occ_func_0_1(6) *
                occ_func_0_0(26))) +
          ((0.7071067812 * occ_func_0_1(5) * occ_func_0_1(7) * occ_func_0_0(0) *
                occ_func_0_0(27) +
            0.7071067812 * occ_func_0_1(5) * occ_func_0_0(7) * occ_func_0_1(0) *
                occ_func_0_0(27))) +
          ((0.7071067812 * occ_func_0_1(83) * occ_func_0_1(35) *
                occ_func_0_0(34) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_1(83) * occ_func_0_0(35) *
                occ_func_0_1(34) * occ_func_0_0(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(2) * occ_func_0_0(3) *
                occ_func_0_0(79) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(2) * occ_func_0_1(3) *
                occ_func_0_0(79))) +
          ((0.7071067812 * occ_func_0_1(11) * occ_func_0_1(0) *
                occ_func_0_0(8) * occ_func_0_0(20) +
            0.7071067812 * occ_func_0_1(11) * occ_func_0_0(0) *
                occ_func_0_1(8) * occ_func_0_0(20))) +
          ((0.7071067812 * occ_func_0_1(10) * occ_func_0_1(5) *
                occ_func_0_0(0) * occ_func_0_0(19) +
            0.7071067812 * occ_func_0_1(10) * occ_func_0_0(5) *
                occ_func_0_1(0) * occ_func_0_0(19))) +
          ((0.7071067812 * occ_func_0_1(86) * occ_func_0_1(41) *
                occ_func_0_0(42) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_1(86) * occ_func_0_0(41) *
                occ_func_0_1(42) * occ_func_0_0(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(11) *
                occ_func_0_0(6) * occ_func_0_0(84) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(11) *
                occ_func_0_1(6) * occ_func_0_0(84))) +
          ((0.7071067812 * occ_func_0_1(2) * occ_func_0_1(0) * occ_func_0_0(1) *
                occ_func_0_0(31) +
            0.7071067812 * occ_func_0_1(2) * occ_func_0_0(0) * occ_func_0_1(1) *
                occ_func_0_0(31))) +
          ((0.7071067812 * occ_func_0_1(7) * occ_func_0_1(12) *
                occ_func_0_0(0) * occ_func_0_0(38) +
            0.7071067812 * occ_func_0_1(7) * occ_func_0_0(12) *
                occ_func_0_1(0) * occ_func_0_0(38))) +
          ((0.7071067812 * occ_func_0_1(81) * occ_func_0_1(30) *
                occ_func_0_0(23) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_1(81) * occ_func_0_0(30) *
                occ_func_0_1(23) * occ_func_0_0(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(5) * occ_func_0_0(7) *
                occ_func_0_0(81) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(5) * occ_func_0_1(7) *
                occ_func_0_0(81))) +
          ((0.7071067812 * occ_func_0_1(8) * occ_func_0_1(0) * occ_func_0_0(9) *
                occ_func_0_0(25) +
            0.7071067812 * occ_func_0_1(8) * occ_func_0_0(0) * occ_func_0_1(9) *
                occ_func_0_0(25))) +
          ((0.7071067812 * occ_func_0_1(6) * occ_func_0_1(4) * occ_func_0_0(0) *
                occ_func_0_0(23) +
            0.7071067812 * occ_func_0_1(6) * occ_func_0_0(4) * occ_func_0_1(0) *
                occ_func_0_0(23))) +
          ((0.7071067812 * occ_func_0_1(84) * occ_func_0_1(36) *
                occ_func_0_0(38) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_1(84) * occ_func_0_0(36) *
                occ_func_0_1(38) * occ_func_0_0(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(1) * occ_func_0_0(6) *
                occ_func_0_0(80) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(1) * occ_func_0_1(6) *
                occ_func_0_0(80))) +
          ((0.7071067812 * occ_func_0_1(12) * occ_func_0_1(0) *
                occ_func_0_0(11) * occ_func_0_0(29) +
            0.7071067812 * occ_func_0_1(12) * occ_func_0_0(0) *
                occ_func_0_1(11) * occ_func_0_0(29))) +
          ((0.7071067812 * occ_func_0_1(7) * occ_func_0_1(2) * occ_func_0_0(0) *
                occ_func_0_0(22) +
            0.7071067812 * occ_func_0_1(7) * occ_func_0_0(2) * occ_func_0_1(0) *
                occ_func_0_0(22))) +
          ((0.7071067812 * occ_func_0_1(85) * occ_func_0_1(32) *
                occ_func_0_0(39) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_1(85) * occ_func_0_0(32) *
                occ_func_0_1(39) * occ_func_0_0(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(4) *
                occ_func_0_0(10) * occ_func_0_0(83) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(4) *
                occ_func_0_1(10) * occ_func_0_0(83))) +
          ((0.7071067812 * occ_func_0_1(9) * occ_func_0_1(0) *
                occ_func_0_0(12) * occ_func_0_0(35) +
            0.7071067812 * occ_func_0_1(9) * occ_func_0_0(0) *
                occ_func_0_1(12) * occ_func_0_0(35))) +
          ((0.7071067812 * occ_func_0_1(3) * occ_func_0_1(1) * occ_func_0_0(0) *
                occ_func_0_0(28) +
            0.7071067812 * occ_func_0_1(3) * occ_func_0_0(1) * occ_func_0_1(0) *
                occ_func_0_0(28))) +
          ((0.7071067812 * occ_func_0_1(82) * occ_func_0_1(26) *
                occ_func_0_0(33) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_1(82) * occ_func_0_0(26) *
                occ_func_0_1(33) * occ_func_0_0(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(4) * occ_func_0_0(1) *
                occ_func_0_0(80) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(4) * occ_func_0_1(1) *
                occ_func_0_0(80))) +
          ((0.7071067812 * occ_func_0_1(9) * occ_func_0_1(0) * occ_func_0_0(3) *
                occ_func_0_0(24) +
            0.7071067812 * occ_func_0_1(9) * occ_func_0_0(0) * occ_func_0_1(3) *
                occ_func_0_0(24))) +
          ((0.7071067812 * occ_func_0_1(12) * occ_func_0_1(10) *
                occ_func_0_0(0) * occ_func_0_0(29) +
            0.7071067812 * occ_func_0_1(12) * occ_func_0_0(10) *
                occ_func_0_1(0) * occ_func_0_0(29))) +
          ((0.7071067812 * occ_func_0_1(85) * occ_func_0_1(37) *
                occ_func_0_0(32) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_1(85) * occ_func_0_0(37) *
                occ_func_0_1(32) * occ_func_0_0(0)))) /
         24.0;
}
double test_Clexulator::site_eval_at_0_bfunc_4_0_4() const {
  return ((occ_func_0_0(0) * occ_func_0_1(9) * occ_func_0_1(7) *
           occ_func_0_0(85)) +
          (occ_func_0_0(4) * occ_func_0_1(0) * occ_func_0_1(5) *
           occ_func_0_0(37)) +
          (occ_func_0_0(6) * occ_func_0_1(8) * occ_func_0_1(0) *
           occ_func_0_0(39)) +
          (occ_func_0_0(80) * occ_func_0_1(24) * occ_func_0_1(22) *
           occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_1(3) * occ_func_0_1(9) *
           occ_func_0_0(82)) +
          (occ_func_0_0(10) * occ_func_0_1(0) * occ_func_0_1(12) *
           occ_func_0_0(33)) +
          (occ_func_0_0(4) * occ_func_0_1(1) * occ_func_0_1(0) *
           occ_func_0_0(26)) +
          (occ_func_0_0(83) * occ_func_0_1(28) * occ_func_0_1(35) *
           occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_1(2) * occ_func_0_1(5) *
           occ_func_0_0(81)) +
          (occ_func_0_0(11) * occ_func_0_1(0) * occ_func_0_1(10) *
           occ_func_0_0(30)) +
          (occ_func_0_0(8) * occ_func_0_1(3) * occ_func_0_1(0) *
           occ_func_0_0(25)) +
          (occ_func_0_0(84) * occ_func_0_1(31) * occ_func_0_1(36) *
           occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_1(12) * occ_func_0_1(10) *
           occ_func_0_0(86)) +
          (occ_func_0_0(1) * occ_func_0_1(0) * occ_func_0_1(4) *
           occ_func_0_0(40)) +
          (occ_func_0_0(3) * occ_func_0_1(9) * occ_func_0_1(0) *
           occ_func_0_0(42)) +
          (occ_func_0_0(79) * occ_func_0_1(21) * occ_func_0_1(19) *
           occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_1(8) * occ_func_0_1(3) *
           occ_func_0_0(82)) +
          (occ_func_0_0(5) * occ_func_0_1(0) * occ_func_0_1(2) *
           occ_func_0_0(27)) +
          (occ_func_0_0(10) * occ_func_0_1(11) * occ_func_0_1(0) *
           occ_func_0_0(33)) +
          (occ_func_0_0(83) * occ_func_0_1(34) * occ_func_0_1(28) *
           occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_1(11) * occ_func_0_1(12) *
           occ_func_0_0(86)) +
          (occ_func_0_0(2) * occ_func_0_1(0) * occ_func_0_1(7) *
           occ_func_0_0(41)) +
          (occ_func_0_0(1) * occ_func_0_1(6) * occ_func_0_1(0) *
           occ_func_0_0(40)) +
          (occ_func_0_0(79) * occ_func_0_1(20) * occ_func_0_1(21) *
           occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_1(7) * occ_func_0_1(2) *
           occ_func_0_0(81)) +
          (occ_func_0_0(6) * occ_func_0_1(0) * occ_func_0_1(1) *
           occ_func_0_0(23)) +
          (occ_func_0_0(11) * occ_func_0_1(12) * occ_func_0_1(0) *
           occ_func_0_0(30)) +
          (occ_func_0_0(84) * occ_func_0_1(38) * occ_func_0_1(31) *
           occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_1(6) * occ_func_0_1(8) *
           occ_func_0_0(84)) +
          (occ_func_0_0(7) * occ_func_0_1(0) * occ_func_0_1(9) *
           occ_func_0_0(38)) +
          (occ_func_0_0(5) * occ_func_0_1(4) * occ_func_0_1(0) *
           occ_func_0_0(36)) +
          (occ_func_0_0(81) * occ_func_0_1(23) * occ_func_0_1(25) *
           occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_1(8) * occ_func_0_1(11) *
           occ_func_0_0(84)) +
          (occ_func_0_0(5) * occ_func_0_1(0) * occ_func_0_1(10) *
           occ_func_0_0(36)) +
          (occ_func_0_0(2) * occ_func_0_1(3) * occ_func_0_1(0) *
           occ_func_0_0(31)) +
          (occ_func_0_0(81) * occ_func_0_1(25) * occ_func_0_1(30) *
           occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_1(12) * occ_func_0_1(9) *
           occ_func_0_0(85)) +
          (occ_func_0_0(1) * occ_func_0_1(0) * occ_func_0_1(3) *
           occ_func_0_0(32)) +
          (occ_func_0_0(4) * occ_func_0_1(10) * occ_func_0_1(0) *
           occ_func_0_0(37)) +
          (occ_func_0_0(80) * occ_func_0_1(29) * occ_func_0_1(24) *
           occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_1(5) * occ_func_0_1(4) *
           occ_func_0_0(83)) +
          (occ_func_0_0(8) * occ_func_0_1(0) * occ_func_0_1(6) *
           occ_func_0_0(34)) +
          (occ_func_0_0(9) * occ_func_0_1(7) * occ_func_0_1(0) *
           occ_func_0_0(35)) +
          (occ_func_0_0(82) * occ_func_0_1(27) * occ_func_0_1(26) *
           occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_1(7) * occ_func_0_1(12) *
           occ_func_0_0(85)) +
          (occ_func_0_0(6) * occ_func_0_1(0) * occ_func_0_1(11) *
           occ_func_0_0(39)) +
          (occ_func_0_0(1) * occ_func_0_1(2) * occ_func_0_1(0) *
           occ_func_0_0(32)) +
          (occ_func_0_0(80) * occ_func_0_1(22) * occ_func_0_1(29) *
           occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_1(1) * occ_func_0_1(2) *
           occ_func_0_0(79)) +
          (occ_func_0_0(12) * occ_func_0_1(0) * occ_func_0_1(7) *
           occ_func_0_0(21)) +
          (occ_func_0_0(11) * occ_func_0_1(6) * occ_func_0_1(0) *
           occ_func_0_0(20)) +
          (occ_func_0_0(86) * occ_func_0_1(40) * occ_func_0_1(41) *
           occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_1(3) * occ_func_0_1(1) *
           occ_func_0_0(79)) +
          (occ_func_0_0(10) * occ_func_0_1(0) * occ_func_0_1(4) *
           occ_func_0_0(19)) +
          (occ_func_0_0(12) * occ_func_0_1(9) * occ_func_0_1(0) *
           occ_func_0_0(21)) +
          (occ_func_0_0(86) * occ_func_0_1(42) * occ_func_0_1(40) *
           occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_1(10) * occ_func_0_1(5) *
           occ_func_0_0(83)) +
          (occ_func_0_0(3) * occ_func_0_1(0) * occ_func_0_1(2) *
           occ_func_0_0(28)) +
          (occ_func_0_0(8) * occ_func_0_1(11) * occ_func_0_1(0) *
           occ_func_0_0(34)) +
          (occ_func_0_0(82) * occ_func_0_1(33) * occ_func_0_1(27) *
           occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_1(10) * occ_func_0_1(11) *
           occ_func_0_0(86)) +
          (occ_func_0_0(3) * occ_func_0_1(0) * occ_func_0_1(8) *
           occ_func_0_0(42)) +
          (occ_func_0_0(2) * occ_func_0_1(5) * occ_func_0_1(0) *
           occ_func_0_0(41)) +
          (occ_func_0_0(79) * occ_func_0_1(19) * occ_func_0_1(20) *
           occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_1(6) * occ_func_0_1(4) *
           occ_func_0_0(80)) +
          (occ_func_0_0(7) * occ_func_0_1(0) * occ_func_0_1(5) *
           occ_func_0_0(22)) +
          (occ_func_0_0(9) * occ_func_0_1(8) * occ_func_0_1(0) *
           occ_func_0_0(24)) +
          (occ_func_0_0(85) * occ_func_0_1(39) * occ_func_0_1(37) *
           occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_1(9) * occ_func_0_1(8) *
           occ_func_0_0(82)) +
          (occ_func_0_0(4) * occ_func_0_1(0) * occ_func_0_1(6) *
           occ_func_0_0(26)) +
          (occ_func_0_0(5) * occ_func_0_1(7) * occ_func_0_1(0) *
           occ_func_0_0(27)) +
          (occ_func_0_0(83) * occ_func_0_1(35) * occ_func_0_1(34) *
           occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_1(2) * occ_func_0_1(3) *
           occ_func_0_0(79)) +
          (occ_func_0_0(11) * occ_func_0_1(0) * occ_func_0_1(8) *
           occ_func_0_0(20)) +
          (occ_func_0_0(10) * occ_func_0_1(5) * occ_func_0_1(0) *
           occ_func_0_0(19)) +
          (occ_func_0_0(86) * occ_func_0_1(41) * occ_func_0_1(42) *
           occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_1(11) * occ_func_0_1(6) *
           occ_func_0_0(84)) +
          (occ_func_0_0(2) * occ_func_0_1(0) * occ_func_0_1(1) *
           occ_func_0_0(31)) +
          (occ_func_0_0(7) * occ_func_0_1(12) * occ_func_0_1(0) *
           occ_func_0_0(38)) +
          (occ_func_0_0(81) * occ_func_0_1(30) * occ_func_0_1(23) *
           occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_1(5) * occ_func_0_1(7) *
           occ_func_0_0(81)) +
          (occ_func_0_0(8) * occ_func_0_1(0) * occ_func_0_1(9) *
           occ_func_0_0(25)) +
          (occ_func_0_0(6) * occ_func_0_1(4) * occ_func_0_1(0) *
           occ_func_0_0(23)) +
          (occ_func_0_0(84) * occ_func_0_1(36) * occ_func_0_1(38) *
           occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_1(1) * occ_func_0_1(6) *
           occ_func_0_0(80)) +
          (occ_func_0_0(12) * occ_func_0_1(0) * occ_func_0_1(11) *
           occ_func_0_0(29)) +
          (occ_func_0_0(7) * occ_func_0_1(2) * occ_func_0_1(0) *
           occ_func_0_0(22)) +
          (occ_func_0_0(85) * occ_func_0_1(32) * occ_func_0_1(39) *
           occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_1(4) * occ_func_0_1(10) *
           occ_func_0_0(83)) +
          (occ_func_0_0(9) * occ_func_0_1(0) * occ_func_0_1(12) *
           occ_func_0_0(35)) +
          (occ_func_0_0(3) * occ_func_0_1(1) * occ_func_0_1(0) *
           occ_func_0_0(28)) +
          (occ_func_0_0(82) * occ_func_0_1(26) * occ_func_0_1(33) *
           occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_1(4) * occ_func_0_1(1) *
           occ_func_0_0(80)) +
          (occ_func_0_0(9) * occ_func_0_1(0) * occ_func_0_1(3) *
           occ_func_0_0(24)) +
          (occ_func_0_0(12) * occ_func_0_1(10) * occ_func_0_1(0) *
           occ_func_0_0(29)) +
          (occ_func_0_0(85) * occ_func_0_1(37) * occ_func_0_1(32) *
           occ_func_0_0(0))) /
         24.0;
}
double test_Clexulator::site_eval_at_0_bfunc_4_0_5() const {
  return ((occ_func_0_1(0) * occ_func_0_1(9) * occ_func_0_1(7) *
           occ_func_0_0(85)) +
          (occ_func_0_1(4) * occ_func_0_1(0) * occ_func_0_1(5) *
           occ_func_0_0(37)) +
          (occ_func_0_1(6) * occ_func_0_1(8) * occ_func_0_1(0) *
           occ_func_0_0(39)) +
          (occ_func_0_1(80) * occ_func_0_1(24) * occ_func_0_1(22) *
           occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_1(3) * occ_func_0_1(9) *
           occ_func_0_0(82)) +
          (occ_func_0_1(10) * occ_func_0_1(0) * occ_func_0_1(12) *
           occ_func_0_0(33)) +
          (occ_func_0_1(4) * occ_func_0_1(1) * occ_func_0_1(0) *
           occ_func_0_0(26)) +
          (occ_func_0_1(83) * occ_func_0_1(28) * occ_func_0_1(35) *
           occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_1(2) * occ_func_0_1(5) *
           occ_func_0_0(81)) +
          (occ_func_0_1(11) * occ_func_0_1(0) * occ_func_0_1(10) *
           occ_func_0_0(30)) +
          (occ_func_0_1(8) * occ_func_0_1(3) * occ_func_0_1(0) *
           occ_func_0_0(25)) +
          (occ_func_0_1(84) * occ_func_0_1(31) * occ_func_0_1(36) *
           occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_1(12) * occ_func_0_1(10) *
           occ_func_0_0(86)) +
          (occ_func_0_1(1) * occ_func_0_1(0) * occ_func_0_1(4) *
           occ_func_0_0(40)) +
          (occ_func_0_1(3) * occ_func_0_1(9) * occ_func_0_1(0) *
           occ_func_0_0(42)) +
          (occ_func_0_1(79) * occ_func_0_1(21) * occ_func_0_1(19) *
           occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_1(8) * occ_func_0_1(3) *
           occ_func_0_0(82)) +
          (occ_func_0_1(5) * occ_func_0_1(0) * occ_func_0_1(2) *
           occ_func_0_0(27)) +
          (occ_func_0_1(10) * occ_func_0_1(11) * occ_func_0_1(0) *
           occ_func_0_0(33)) +
          (occ_func_0_1(83) * occ_func_0_1(34) * occ_func_0_1(28) *
           occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_1(11) * occ_func_0_1(12) *
           occ_func_0_0(86)) +
          (occ_func_0_1(2) * occ_func_0_1(0) * occ_func_0_1(7) *
           occ_func_0_0(41)) +
          (occ_func_0_1(1) * occ_func_0_1(6) * occ_func_0_1(0) *
           occ_func_0_0(40)) +
          (occ_func_0_1(79) * occ_func_0_1(20) * occ_func_0_1(21) *
           occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_1(7) * occ_func_0_1(2) *
           occ_func_0_0(81)) +
          (occ_func_0_1(6) * occ_func_0_1(0) * occ_func_0_1(1) *
           occ_func_0_0(23)) +
          (occ_func_0_1(11) * occ_func_0_1(12) * occ_func_0_1(0) *
           occ_func_0_0(30)) +
          (occ_func_0_1(84) * occ_func_0_1(38) * occ_func_0_1(31) *
           occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_1(6) * occ_func_0_1(8) *
           occ_func_0_0(84)) +
          (occ_func_0_1(7) * occ_func_0_1(0) * occ_func_0_1(9) *
           occ_func_0_0(38)) +
          (occ_func_0_1(5) * occ_func_0_1(4) * occ_func_0_1(0) *
           occ_func_0_0(36)) +
          (occ_func_0_1(81) * occ_func_0_1(23) * occ_func_0_1(25) *
           occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_1(8) * occ_func_0_1(11) *
           occ_func_0_0(84)) +
          (occ_func_0_1(5) * occ_func_0_1(0) * occ_func_0_1(10) *
           occ_func_0_0(36)) +
          (occ_func_0_1(2) * occ_func_0_1(3) * occ_func_0_1(0) *
           occ_func_0_0(31)) +
          (occ_func_0_1(81) * occ_func_0_1(25) * occ_func_0_1(30) *
           occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_1(12) * occ_func_0_1(9) *
           occ_func_0_0(85)) +
          (occ_func_0_1(1) * occ_func_0_1(0) * occ_func_0_1(3) *
           occ_func_0_0(32)) +
          (occ_func_0_1(4) * occ_func_0_1(10) * occ_func_0_1(0) *
           occ_func_0_0(37)) +
          (occ_func_0_1(80) * occ_func_0_1(29) * occ_func_0_1(24) *
           occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_1(5) * occ_func_0_1(4) *
           occ_func_0_0(83)) +
          (occ_func_0_1(8) * occ_func_0_1(0) * occ_func_0_1(6) *
           occ_func_0_0(34)) +
          (occ_func_0_1(9) * occ_func_0_1(7) * occ_func_0_1(0) *
           occ_func_0_0(35)) +
          (occ_func_0_1(82) * occ_func_0_1(27) * occ_func_0_1(26) *
           occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_1(7) * occ_func_0_1(12) *
           occ_func_0_0(85)) +
          (occ_func_0_1(6) * occ_func_0_1(0) * occ_func_0_1(11) *
           occ_func_0_0(39)) +
          (occ_func_0_1(1) * occ_func_0_1(2) * occ_func_0_1(0) *
           occ_func_0_0(32)) +
          (occ_func_0_1(80) * occ_func_0_1(22) * occ_func_0_1(29) *
           occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_1(1) * occ_func_0_1(2) *
           occ_func_0_0(79)) +
          (occ_func_0_1(12) * occ_func_0_1(0) * occ_func_0_1(7) *
           occ_func_0_0(21)) +
          (occ_func_0_1(11) * occ_func_0_1(6) * occ_func_0_1(0) *
           occ_func_0_0(20)) +
          (occ_func_0_1(86) * occ_func_0_1(40) * occ_func_0_1(41) *
           occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_1(3) * occ_func_0_1(1) *
           occ_func_0_0(79)) +
          (occ_func_0_1(10) * occ_func_0_1(0) * occ_func_0_1(4) *
           occ_func_0_0(19)) +
          (occ_func_0_1(12) * occ_func_0_1(9) * occ_func_0_1(0) *
           occ_func_0_0(21)) +
          (occ_func_0_1(86) * occ_func_0_1(42) * occ_func_0_1(40) *
           occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_1(10) * occ_func_0_1(5) *
           occ_func_0_0(83)) +
          (occ_func_0_1(3) * occ_func_0_1(0) * occ_func_0_1(2) *
           occ_func_0_0(28)) +
          (occ_func_0_1(8) * occ_func_0_1(11) * occ_func_0_1(0) *
           occ_func_0_0(34)) +
          (occ_func_0_1(82) * occ_func_0_1(33) * occ_func_0_1(27) *
           occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_1(10) * occ_func_0_1(11) *
           occ_func_0_0(86)) +
          (occ_func_0_1(3) * occ_func_0_1(0) * occ_func_0_1(8) *
           occ_func_0_0(42)) +
          (occ_func_0_1(2) * occ_func_0_1(5) * occ_func_0_1(0) *
           occ_func_0_0(41)) +
          (occ_func_0_1(79) * occ_func_0_1(19) * occ_func_0_1(20) *
           occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_1(6) * occ_func_0_1(4) *
           occ_func_0_0(80)) +
          (occ_func_0_1(7) * occ_func_0_1(0) * occ_func_0_1(5) *
           occ_func_0_0(22)) +
          (occ_func_0_1(9) * occ_func_0_1(8) * occ_func_0_1(0) *
           occ_func_0_0(24)) +
          (occ_func_0_1(85) * occ_func_0_1(39) * occ_func_0_1(37) *
           occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_1(9) * occ_func_0_1(8) *
           occ_func_0_0(82)) +
          (occ_func_0_1(4) * occ_func_0_1(0) * occ_func_0_1(6) *
           occ_func_0_0(26)) +
          (occ_func_0_1(5) * occ_func_0_1(7) * occ_func_0_1(0) *
           occ_func_0_0(27)) +
          (occ_func_0_1(83) * occ_func_0_1(35) * occ_func_0_1(34) *
           occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_1(2) * occ_func_0_1(3) *
           occ_func_0_0(79)) +
          (occ_func_0_1(11) * occ_func_0_1(0) * occ_func_0_1(8) *
           occ_func_0_0(20)) +
          (occ_func_0_1(10) * occ_func_0_1(5) * occ_func_0_1(0) *
           occ_func_0_0(19)) +
          (occ_func_0_1(86) * occ_func_0_1(41) * occ_func_0_1(42) *
           occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_1(11) * occ_func_0_1(6) *
           occ_func_0_0(84)) +
          (occ_func_0_1(2) * occ_func_0_1(0) * occ_func_0_1(1) *
           occ_func_0_0(31)) +
          (occ_func_0_1(7) * occ_func_0_1(12) * occ_func_0_1(0) *
           occ_func_0_0(38)) +
          (occ_func_0_1(81) * occ_func_0_1(30) * occ_func_0_1(23) *
           occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_1(5) * occ_func_0_1(7) *
           occ_func_0_0(81)) +
          (occ_func_0_1(8) * occ_func_0_1(0) * occ_func_0_1(9) *
           occ_func_0_0(25)) +
          (occ_func_0_1(6) * occ_func_0_1(4) * occ_func_0_1(0) *
           occ_func_0_0(23)) +
          (occ_func_0_1(84) * occ_func_0_1(36) * occ_func_0_1(38) *
           occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_1(1) * occ_func_0_1(6) *
           occ_func_0_0(80)) +
          (occ_func_0_1(12) * occ_func_0_1(0) * occ_func_0_1(11) *
           occ_func_0_0(29)) +
          (occ_func_0_1(7) * occ_func_0_1(2) * occ_func_0_1(0) *
           occ_func_0_0(22)) +
          (occ_func_0_1(85) * occ_func_0_1(32) * occ_func_0_1(39) *
           occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_1(4) * occ_func_0_1(10) *
           occ_func_0_0(83)) +
          (occ_func_0_1(9) * occ_func_0_1(0) * occ_func_0_1(12) *
           occ_func_0_0(35)) +
          (occ_func_0_1(3) * occ_func_0_1(1) * occ_func_0_1(0) *
           occ_func_0_0(28)) +
          (occ_func_0_1(82) * occ_func_0_1(26) * occ_func_0_1(33) *
           occ_func_0_0(0)) +
          (occ_func_0_1(0) * occ_func_0_1(4) * occ_func_0_1(1) *
           occ_func_0_0(80)) +
          (occ_func_0_1(9) * occ_func_0_1(0) * occ_func_0_1(3) *
           occ_func_0_0(24)) +
          (occ_func_0_1(12) * occ_func_0_1(10) * occ_func_0_1(0) *
           occ_func_0_0(29)) +
          (occ_func_0_1(85) * occ_func_0_1(37) * occ_func_0_1(32) *
           occ_func_0_0(0))) /
         24.0;
}
double test_Clexulator::site_eval_at_0_bfunc_4_0_6() const {
  return ((occ_func_0_0(0) * occ_func_0_0(9) * occ_func_0_0(7) *
           occ_func_0_1(85)) +
          (occ_func_0_0(4) * occ_func_0_0(0) * occ_func_0_0(5) *
           occ_func_0_1(37)) +
          (occ_func_0_0(6) * occ_func_0_0(8) * occ_func_0_0(0) *
           occ_func_0_1(39)) +
          (occ_func_0_0(80) * occ_func_0_0(24) * occ_func_0_0(22) *
           occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_0(3) * occ_func_0_0(9) *
           occ_func_0_1(82)) +
          (occ_func_0_0(10) * occ_func_0_0(0) * occ_func_0_0(12) *
           occ_func_0_1(33)) +
          (occ_func_0_0(4) * occ_func_0_0(1) * occ_func_0_0(0) *
           occ_func_0_1(26)) +
          (occ_func_0_0(83) * occ_func_0_0(28) * occ_func_0_0(35) *
           occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_0(2) * occ_func_0_0(5) *
           occ_func_0_1(81)) +
          (occ_func_0_0(11) * occ_func_0_0(0) * occ_func_0_0(10) *
           occ_func_0_1(30)) +
          (occ_func_0_0(8) * occ_func_0_0(3) * occ_func_0_0(0) *
           occ_func_0_1(25)) +
          (occ_func_0_0(84) * occ_func_0_0(31) * occ_func_0_0(36) *
           occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_0(12) * occ_func_0_0(10) *
           occ_func_0_1(86)) +
          (occ_func_0_0(1) * occ_func_0_0(0) * occ_func_0_0(4) *
           occ_func_0_1(40)) +
          (occ_func_0_0(3) * occ_func_0_0(9) * occ_func_0_0(0) *
           occ_func_0_1(42)) +
          (occ_func_0_0(79) * occ_func_0_0(21) * occ_func_0_0(19) *
           occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_0(8) * occ_func_0_0(3) *
           occ_func_0_1(82)) +
          (occ_func_0_0(5) * occ_func_0_0(0) * occ_func_0_0(2) *
           occ_func_0_1(27)) +
          (occ_func_0_0(10) * occ_func_0_0(11) * occ_func_0_0(0) *
           occ_func_0_1(33)) +
          (occ_func_0_0(83) * occ_func_0_0(34) * occ_func_0_0(28) *
           occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_0(11) * occ_func_0_0(12) *
           occ_func_0_1(86)) +
          (occ_func_0_0(2) * occ_func_0_0(0) * occ_func_0_0(7) *
           occ_func_0_1(41)) +
          (occ_func_0_0(1) * occ_func_0_0(6) * occ_func_0_0(0) *
           occ_func_0_1(40)) +
          (occ_func_0_0(79) * occ_func_0_0(20) * occ_func_0_0(21) *
           occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_0(7) * occ_func_0_0(2) *
           occ_func_0_1(81)) +
          (occ_func_0_0(6) * occ_func_0_0(0) * occ_func_0_0(1) *
           occ_func_0_1(23)) +
          (occ_func_0_0(11) * occ_func_0_0(12) * occ_func_0_0(0) *
           occ_func_0_1(30)) +
          (occ_func_0_0(84) * occ_func_0_0(38) * occ_func_0_0(31) *
           occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_0(6) * occ_func_0_0(8) *
           occ_func_0_1(84)) +
          (occ_func_0_0(7) * occ_func_0_0(0) * occ_func_0_0(9) *
           occ_func_0_1(38)) +
          (occ_func_0_0(5) * occ_func_0_0(4) * occ_func_0_0(0) *
           occ_func_0_1(36)) +
          (occ_func_0_0(81) * occ_func_0_0(23) * occ_func_0_0(25) *
           occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_0(8) * occ_func_0_0(11) *
           occ_func_0_1(84)) +
          (occ_func_0_0(5) * occ_func_0_0(0) * occ_func_0_0(10) *
           occ_func_0_1(36)) +
          (occ_func_0_0(2) * occ_func_0_0(3) * occ_func_0_0(0) *
           occ_func_0_1(31)) +
          (occ_func_0_0(81) * occ_func_0_0(25) * occ_func_0_0(30) *
           occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_0(12) * occ_func_0_0(9) *
           occ_func_0_1(85)) +
          (occ_func_0_0(1) * occ_func_0_0(0) * occ_func_0_0(3) *
           occ_func_0_1(32)) +
          (occ_func_0_0(4) * occ_func_0_0(10) * occ_func_0_0(0) *
           occ_func_0_1(37)) +
          (occ_func_0_0(80) * occ_func_0_0(29) * occ_func_0_0(24) *
           occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_0(5) * occ_func_0_0(4) *
           occ_func_0_1(83)) +
          (occ_func_0_0(8) * occ_func_0_0(0) * occ_func_0_0(6) *
           occ_func_0_1(34)) +
          (occ_func_0_0(9) * occ_func_0_0(7) * occ_func_0_0(0) *
           occ_func_0_1(35)) +
          (occ_func_0_0(82) * occ_func_0_0(27) * occ_func_0_0(26) *
           occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_0(7) * occ_func_0_0(12) *
           occ_func_0_1(85)) +
          (occ_func_0_0(6) * occ_func_0_0(0) * occ_func_0_0(11) *
           occ_func_0_1(39)) +
          (occ_func_0_0(1) * occ_func_0_0(2) * occ_func_0_0(0) *
           occ_func_0_1(32)) +
          (occ_func_0_0(80) * occ_func_0_0(22) * occ_func_0_0(29) *
           occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_0(1) * occ_func_0_0(2) *
           occ_func_0_1(79)) +
          (occ_func_0_0(12) * occ_func_0_0(0) * occ_func_0_0(7) *
           occ_func_0_1(21)) +
          (occ_func_0_0(11) * occ_func_0_0(6) * occ_func_0_0(0) *
           occ_func_0_1(20)) +
          (occ_func_0_0(86) * occ_func_0_0(40) * occ_func_0_0(41) *
           occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_0(3) * occ_func_0_0(1) *
           occ_func_0_1(79)) +
          (occ_func_0_0(10) * occ_func_0_0(0) * occ_func_0_0(4) *
           occ_func_0_1(19)) +
          (occ_func_0_0(12) * occ_func_0_0(9) * occ_func_0_0(0) *
           occ_func_0_1(21)) +
          (occ_func_0_0(86) * occ_func_0_0(42) * occ_func_0_0(40) *
           occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_0(10) * occ_func_0_0(5) *
           occ_func_0_1(83)) +
          (occ_func_0_0(3) * occ_func_0_0(0) * occ_func_0_0(2) *
           occ_func_0_1(28)) +
          (occ_func_0_0(8) * occ_func_0_0(11) * occ_func_0_0(0) *
           occ_func_0_1(34)) +
          (occ_func_0_0(82) * occ_func_0_0(33) * occ_func_0_0(27) *
           occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_0(10) * occ_func_0_0(11) *
           occ_func_0_1(86)) +
          (occ_func_0_0(3) * occ_func_0_0(0) * occ_func_0_0(8) *
           occ_func_0_1(42)) +
          (occ_func_0_0(2) * occ_func_0_0(5) * occ_func_0_0(0) *
           occ_func_0_1(41)) +
          (occ_func_0_0(79) * occ_func_0_0(19) * occ_func_0_0(20) *
           occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_0(6) * occ_func_0_0(4) *
           occ_func_0_1(80)) +
          (occ_func_0_0(7) * occ_func_0_0(0) * occ_func_0_0(5) *
           occ_func_0_1(22)) +
          (occ_func_0_0(9) * occ_func_0_0(8) * occ_func_0_0(0) *
           occ_func_0_1(24)) +
          (occ_func_0_0(85) * occ_func_0_0(39) * occ_func_0_0(37) *
           occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_0(9) * occ_func_0_0(8) *
           occ_func_0_1(82)) +
          (occ_func_0_0(4) * occ_func_0_0(0) * occ_func_0_0(6) *
           occ_func_0_1(26)) +
          (occ_func_0_0(5) * occ_func_0_0(7) * occ_func_0_0(0) *
           occ_func_0_1(27)) +
          (occ_func_0_0(83) * occ_func_0_0(35) * occ_func_0_0(34) *
           occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_0(2) * occ_func_0_0(3) *
           occ_func_0_1(79)) +
          (occ_func_0_0(11) * occ_func_0_0(0) * occ_func_0_0(8) *
           occ_func_0_1(20)) +
          (occ_func_0_0(10) * occ_func_0_0(5) * occ_func_0_0(0) *
           occ_func_0_1(19)) +
          (occ_func_0_0(86) * occ_func_0_0(41) * occ_func_0_0(42) *
           occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_0(11) * occ_func_0_0(6) *
           occ_func_0_1(84)) +
          (occ_func_0_0(2) * occ_func_0_0(0) * occ_func_0_0(1) *
           occ_func_0_1(31)) +
          (occ_func_0_0(7) * occ_func_0_0(12) * occ_func_0_0(0) *
           occ_func_0_1(38)) +
          (occ_func_0_0(81) * occ_func_0_0(30) * occ_func_0_0(23) *
           occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_0(5) * occ_func_0_0(7) *
           occ_func_0_1(81)) +
          (occ_func_0_0(8) * occ_func_0_0(0) * occ_func_0_0(9) *
           occ_func_0_1(25)) +
          (occ_func_0_0(6) * occ_func_0_0(4) * occ_func_0_0(0) *
           occ_func_0_1(23)) +
          (occ_func_0_0(84) * occ_func_0_0(36) * occ_func_0_0(38) *
           occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_0(1) * occ_func_0_0(6) *
           occ_func_0_1(80)) +
          (occ_func_0_0(12) * occ_func_0_0(0) * occ_func_0_0(11) *
           occ_func_0_1(29)) +
          (occ_func_0_0(7) * occ_func_0_0(2) * occ_func_0_0(0) *
           occ_func_0_1(22)) +
          (occ_func_0_0(85) * occ_func_0_0(32) * occ_func_0_0(39) *
           occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_0(4) * occ_func_0_0(10) *
           occ_func_0_1(83)) +
          (occ_func_0_0(9) * occ_func_0_0(0) * occ_func_0_0(12) *
           occ_func_0_1(35)) +
          (occ_func_0_0(3) * occ_func_0_0(1) * occ_func_0_0(0) *
           occ_func_0_1(28)) +
          (occ_func_0_0(82) * occ_func_0_0(26) * occ_func_0_0(33) *
           occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_0(4) * occ_func_0_0(1) *
           occ_func_0_1(80)) +
          (occ_func_0_0(9) * occ_func_0_0(0) * occ_func_0_0(3) *
           occ_func_0_1(24)) +
          (occ_func_0_0(12) * occ_func_0_0(10) * occ_func_0_0(0) *
           occ_func_0_1(29)) +
          (occ_func_0_0(85) * occ_func_0_0(37) * occ_func_0_0(32) *
           occ_func_0_1(0))) /
         24.0;
}
double test_Clexulator::site_eval_at_0_bfunc_4_0_7() const {
  return ((occ_func_0_1(0) * occ_func_0_0(9) * occ_func_0_0(7) *
           occ_func_0_1(85)) +
          (occ_func_0_1(4) * occ_func_0_0(0) * occ_func_0_0(5) *
           occ_func_0_1(37)) +
          (occ_func_0_1(6) * occ_func_0_0(8) * occ_func_0_0(0) *
           occ_func_0_1(39)) +
          (occ_func_0_1(80) * occ_func_0_0(24) * occ_func_0_0(22) *
           occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_0(3) * occ_func_0_0(9) *
           occ_func_0_1(82)) +
          (occ_func_0_1(10) * occ_func_0_0(0) * occ_func_0_0(12) *
           occ_func_0_1(33)) +
          (occ_func_0_1(4) * occ_func_0_0(1) * occ_func_0_0(0) *
           occ_func_0_1(26)) +
          (occ_func_0_1(83) * occ_func_0_0(28) * occ_func_0_0(35) *
           occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_0(2) * occ_func_0_0(5) *
           occ_func_0_1(81)) +
          (occ_func_0_1(11) * occ_func_0_0(0) * occ_func_0_0(10) *
           occ_func_0_1(30)) +
          (occ_func_0_1(8) * occ_func_0_0(3) * occ_func_0_0(0) *
           occ_func_0_1(25)) +
          (occ_func_0_1(84) * occ_func_0_0(31) * occ_func_0_0(36) *
           occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_0(12) * occ_func_0_0(10) *
           occ_func_0_1(86)) +
          (occ_func_0_1(1) * occ_func_0_0(0) * occ_func_0_0(4) *
           occ_func_0_1(40)) +
          (occ_func_0_1(3) * occ_func_0_0(9) * occ_func_0_0(0) *
           occ_func_0_1(42)) +
          (occ_func_0_1(79) * occ_func_0_0(21) * occ_func_0_0(19) *
           occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_0(8) * occ_func_0_0(3) *
           occ_func_0_1(82)) +
          (occ_func_0_1(5) * occ_func_0_0(0) * occ_func_0_0(2) *
           occ_func_0_1(27)) +
          (occ_func_0_1(10) * occ_func_0_0(11) * occ_func_0_0(0) *
           occ_func_0_1(33)) +
          (occ_func_0_1(83) * occ_func_0_0(34) * occ_func_0_0(28) *
           occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_0(11) * occ_func_0_0(12) *
           occ_func_0_1(86)) +
          (occ_func_0_1(2) * occ_func_0_0(0) * occ_func_0_0(7) *
           occ_func_0_1(41)) +
          (occ_func_0_1(1) * occ_func_0_0(6) * occ_func_0_0(0) *
           occ_func_0_1(40)) +
          (occ_func_0_1(79) * occ_func_0_0(20) * occ_func_0_0(21) *
           occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_0(7) * occ_func_0_0(2) *
           occ_func_0_1(81)) +
          (occ_func_0_1(6) * occ_func_0_0(0) * occ_func_0_0(1) *
           occ_func_0_1(23)) +
          (occ_func_0_1(11) * occ_func_0_0(12) * occ_func_0_0(0) *
           occ_func_0_1(30)) +
          (occ_func_0_1(84) * occ_func_0_0(38) * occ_func_0_0(31) *
           occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_0(6) * occ_func_0_0(8) *
           occ_func_0_1(84)) +
          (occ_func_0_1(7) * occ_func_0_0(0) * occ_func_0_0(9) *
           occ_func_0_1(38)) +
          (occ_func_0_1(5) * occ_func_0_0(4) * occ_func_0_0(0) *
           occ_func_0_1(36)) +
          (occ_func_0_1(81) * occ_func_0_0(23) * occ_func_0_0(25) *
           occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_0(8) * occ_func_0_0(11) *
           occ_func_0_1(84)) +
          (occ_func_0_1(5) * occ_func_0_0(0) * occ_func_0_0(10) *
           occ_func_0_1(36)) +
          (occ_func_0_1(2) * occ_func_0_0(3) * occ_func_0_0(0) *
           occ_func_0_1(31)) +
          (occ_func_0_1(81) * occ_func_0_0(25) * occ_func_0_0(30) *
           occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_0(12) * occ_func_0_0(9) *
           occ_func_0_1(85)) +
          (occ_func_0_1(1) * occ_func_0_0(0) * occ_func_0_0(3) *
           occ_func_0_1(32)) +
          (occ_func_0_1(4) * occ_func_0_0(10) * occ_func_0_0(0) *
           occ_func_0_1(37)) +
          (occ_func_0_1(80) * occ_func_0_0(29) * occ_func_0_0(24) *
           occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_0(5) * occ_func_0_0(4) *
           occ_func_0_1(83)) +
          (occ_func_0_1(8) * occ_func_0_0(0) * occ_func_0_0(6) *
           occ_func_0_1(34)) +
          (occ_func_0_1(9) * occ_func_0_0(7) * occ_func_0_0(0) *
           occ_func_0_1(35)) +
          (occ_func_0_1(82) * occ_func_0_0(27) * occ_func_0_0(26) *
           occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_0(7) * occ_func_0_0(12) *
           occ_func_0_1(85)) +
          (occ_func_0_1(6) * occ_func_0_0(0) * occ_func_0_0(11) *
           occ_func_0_1(39)) +
          (occ_func_0_1(1) * occ_func_0_0(2) * occ_func_0_0(0) *
           occ_func_0_1(32)) +
          (occ_func_0_1(80) * occ_func_0_0(22) * occ_func_0_0(29) *
           occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_0(1) * occ_func_0_0(2) *
           occ_func_0_1(79)) +
          (occ_func_0_1(12) * occ_func_0_0(0) * occ_func_0_0(7) *
           occ_func_0_1(21)) +
          (occ_func_0_1(11) * occ_func_0_0(6) * occ_func_0_0(0) *
           occ_func_0_1(20)) +
          (occ_func_0_1(86) * occ_func_0_0(40) * occ_func_0_0(41) *
           occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_0(3) * occ_func_0_0(1) *
           occ_func_0_1(79)) +
          (occ_func_0_1(10) * occ_func_0_0(0) * occ_func_0_0(4) *
           occ_func_0_1(19)) +
          (occ_func_0_1(12) * occ_func_0_0(9) * occ_func_0_0(0) *
           occ_func_0_1(21)) +
          (occ_func_0_1(86) * occ_func_0_0(42) * occ_func_0_0(40) *
           occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_0(10) * occ_func_0_0(5) *
           occ_func_0_1(83)) +
          (occ_func_0_1(3) * occ_func_0_0(0) * occ_func_0_0(2) *
           occ_func_0_1(28)) +
          (occ_func_0_1(8) * occ_func_0_0(11) * occ_func_0_0(0) *
           occ_func_0_1(34)) +
          (occ_func_0_1(82) * occ_func_0_0(33) * occ_func_0_0(27) *
           occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_0(10) * occ_func_0_0(11) *
           occ_func_0_1(86)) +
          (occ_func_0_1(3) * occ_func_0_0(0) * occ_func_0_0(8) *
           occ_func_0_1(42)) +
          (occ_func_0_1(2) * occ_func_0_0(5) * occ_func_0_0(0) *
           occ_func_0_1(41)) +
          (occ_func_0_1(79) * occ_func_0_0(19) * occ_func_0_0(20) *
           occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_0(6) * occ_func_0_0(4) *
           occ_func_0_1(80)) +
          (occ_func_0_1(7) * occ_func_0_0(0) * occ_func_0_0(5) *
           occ_func_0_1(22)) +
          (occ_func_0_1(9) * occ_func_0_0(8) * occ_func_0_0(0) *
           occ_func_0_1(24)) +
          (occ_func_0_1(85) * occ_func_0_0(39) * occ_func_0_0(37) *
           occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_0(9) * occ_func_0_0(8) *
           occ_func_0_1(82)) +
          (occ_func_0_1(4) * occ_func_0_0(0) * occ_func_0_0(6) *
           occ_func_0_1(26)) +
          (occ_func_0_1(5) * occ_func_0_0(7) * occ_func_0_0(0) *
           occ_func_0_1(27)) +
          (occ_func_0_1(83) * occ_func_0_0(35) * occ_func_0_0(34) *
           occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_0(2) * occ_func_0_0(3) *
           occ_func_0_1(79)) +
          (occ_func_0_1(11) * occ_func_0_0(0) * occ_func_0_0(8) *
           occ_func_0_1(20)) +
          (occ_func_0_1(10) * occ_func_0_0(5) * occ_func_0_0(0) *
           occ_func_0_1(19)) +
          (occ_func_0_1(86) * occ_func_0_0(41) * occ_func_0_0(42) *
           occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_0(11) * occ_func_0_0(6) *
           occ_func_0_1(84)) +
          (occ_func_0_1(2) * occ_func_0_0(0) * occ_func_0_0(1) *
           occ_func_0_1(31)) +
          (occ_func_0_1(7) * occ_func_0_0(12) * occ_func_0_0(0) *
           occ_func_0_1(38)) +
          (occ_func_0_1(81) * occ_func_0_0(30) * occ_func_0_0(23) *
           occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_0(5) * occ_func_0_0(7) *
           occ_func_0_1(81)) +
          (occ_func_0_1(8) * occ_func_0_0(0) * occ_func_0_0(9) *
           occ_func_0_1(25)) +
          (occ_func_0_1(6) * occ_func_0_0(4) * occ_func_0_0(0) *
           occ_func_0_1(23)) +
          (occ_func_0_1(84) * occ_func_0_0(36) * occ_func_0_0(38) *
           occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_0(1) * occ_func_0_0(6) *
           occ_func_0_1(80)) +
          (occ_func_0_1(12) * occ_func_0_0(0) * occ_func_0_0(11) *
           occ_func_0_1(29)) +
          (occ_func_0_1(7) * occ_func_0_0(2) * occ_func_0_0(0) *
           occ_func_0_1(22)) +
          (occ_func_0_1(85) * occ_func_0_0(32) * occ_func_0_0(39) *
           occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_0(4) * occ_func_0_0(10) *
           occ_func_0_1(83)) +
          (occ_func_0_1(9) * occ_func_0_0(0) * occ_func_0_0(12) *
           occ_func_0_1(35)) +
          (occ_func_0_1(3) * occ_func_0_0(1) * occ_func_0_0(0) *
           occ_func_0_1(28)) +
          (occ_func_0_1(82) * occ_func_0_0(26) * occ_func_0_0(33) *
           occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_0(4) * occ_func_0_0(1) *
           occ_func_0_1(80)) +
          (occ_func_0_1(9) * occ_func_0_0(0) * occ_func_0_0(3) *
           occ_func_0_1(24)) +
          (occ_func_0_1(12) * occ_func_0_0(10) * occ_func_0_0(0) *
           occ_func_0_1(29)) +
          (occ_func_0_1(85) * occ_func_0_0(37) * occ_func_0_0(32) *
           occ_func_0_1(0))) /
         24.0;
}
double test_Clexulator::site_eval_at_0_bfunc_4_0_8() const {
  return (((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(9) * occ_func_0_0(7) *
                occ_func_0_1(85) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(9) * occ_func_0_1(7) *
                occ_func_0_1(85))) +
          ((0.7071067812 * occ_func_0_0(4) * occ_func_0_1(0) * occ_func_0_0(5) *
                occ_func_0_1(37) +
            0.7071067812 * occ_func_0_0(4) * occ_func_0_0(0) * occ_func_0_1(5) *
                occ_func_0_1(37))) +
          ((0.7071067812 * occ_func_0_0(6) * occ_func_0_1(8) * occ_func_0_0(0) *
                occ_func_0_1(39) +
            0.7071067812 * occ_func_0_0(6) * occ_func_0_0(8) * occ_func_0_1(0) *
                occ_func_0_1(39))) +
          ((0.7071067812 * occ_func_0_0(80) * occ_func_0_1(24) *
                occ_func_0_0(22) * occ_func_0_1(0) +
            0.7071067812 * occ_func_0_0(80) * occ_func_0_0(24) *
                occ_func_0_1(22) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(3) * occ_func_0_0(9) *
                occ_func_0_1(82) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(3) * occ_func_0_1(9) *
                occ_func_0_1(82))) +
          ((0.7071067812 * occ_func_0_0(10) * occ_func_0_1(0) *
                occ_func_0_0(12) * occ_func_0_1(33) +
            0.7071067812 * occ_func_0_0(10) * occ_func_0_0(0) *
                occ_func_0_1(12) * occ_func_0_1(33))) +
          ((0.7071067812 * occ_func_0_0(4) * occ_func_0_1(1) * occ_func_0_0(0) *
                occ_func_0_1(26) +
            0.7071067812 * occ_func_0_0(4) * occ_func_0_0(1) * occ_func_0_1(0) *
                occ_func_0_1(26))) +
          ((0.7071067812 * occ_func_0_0(83) * occ_func_0_1(28) *
                occ_func_0_0(35) * occ_func_0_1(0) +
            0.7071067812 * occ_func_0_0(83) * occ_func_0_0(28) *
                occ_func_0_1(35) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(2) * occ_func_0_0(5) *
                occ_func_0_1(81) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(2) * occ_func_0_1(5) *
                occ_func_0_1(81))) +
          ((0.7071067812 * occ_func_0_0(11) * occ_func_0_1(0) *
                occ_func_0_0(10) * occ_func_0_1(30) +
            0.7071067812 * occ_func_0_0(11) * occ_func_0_0(0) *
                occ_func_0_1(10) * occ_func_0_1(30))) +
          ((0.7071067812 * occ_func_0_0(8) * occ_func_0_1(3) * occ_func_0_0(0) *
                occ_func_0_1(25) +
            0.7071067812 * occ_func_0_0(8) * occ_func_0_0(3) * occ_func_0_1(0) *
                occ_func_0_1(25))) +
          ((0.7071067812 * occ_func_0_0(84) * occ_func_0_1(31) *
                occ_func_0_0(36) * occ_func_0_1(0) +
            0.7071067812 * occ_func_0_0(84) * occ_func_0_0(31) *
                occ_func_0_1(36) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(12) *
                occ_func_0_0(10) * occ_func_0_1(86) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(12) *
                occ_func_0_1(10) * occ_func_0_1(86))) +
          ((0.7071067812 * occ_func_0_0(1) * occ_func_0_1(0) * occ_func_0_0(4) *
                occ_func_0_1(40) +
            0.7071067812 * occ_func_0_0(1) * occ_func_0_0(0) * occ_func_0_1(4) *
                occ_func_0_1(40))) +
          ((0.7071067812 * occ_func_0_0(3) * occ_func_0_1(9) * occ_func_0_0(0) *
                occ_func_0_1(42) +
            0.7071067812 * occ_func_0_0(3) * occ_func_0_0(9) * occ_func_0_1(0) *
                occ_func_0_1(42))) +
          ((0.7071067812 * occ_func_0_0(79) * occ_func_0_1(21) *
                occ_func_0_0(19) * occ_func_0_1(0) +
            0.7071067812 * occ_func_0_0(79) * occ_func_0_0(21) *
                occ_func_0_1(19) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(8) * occ_func_0_0(3) *
                occ_func_0_1(82) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(8) * occ_func_0_1(3) *
                occ_func_0_1(82))) +
          ((0.7071067812 * occ_func_0_0(5) * occ_func_0_1(0) * occ_func_0_0(2) *
                occ_func_0_1(27) +
            0.7071067812 * occ_func_0_0(5) * occ_func_0_0(0) * occ_func_0_1(2) *
                occ_func_0_1(27))) +
          ((0.7071067812 * occ_func_0_0(10) * occ_func_0_1(11) *
                occ_func_0_0(0) * occ_func_0_1(33) +
            0.7071067812 * occ_func_0_0(10) * occ_func_0_0(11) *
                occ_func_0_1(0) * occ_func_0_1(33))) +
          ((0.7071067812 * occ_func_0_0(83) * occ_func_0_1(34) *
                occ_func_0_0(28) * occ_func_0_1(0) +
            0.7071067812 * occ_func_0_0(83) * occ_func_0_0(34) *
                occ_func_0_1(28) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(11) *
                occ_func_0_0(12) * occ_func_0_1(86) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(11) *
                occ_func_0_1(12) * occ_func_0_1(86))) +
          ((0.7071067812 * occ_func_0_0(2) * occ_func_0_1(0) * occ_func_0_0(7) *
                occ_func_0_1(41) +
            0.7071067812 * occ_func_0_0(2) * occ_func_0_0(0) * occ_func_0_1(7) *
                occ_func_0_1(41))) +
          ((0.7071067812 * occ_func_0_0(1) * occ_func_0_1(6) * occ_func_0_0(0) *
                occ_func_0_1(40) +
            0.7071067812 * occ_func_0_0(1) * occ_func_0_0(6) * occ_func_0_1(0) *
                occ_func_0_1(40))) +
          ((0.7071067812 * occ_func_0_0(79) * occ_func_0_1(20) *
                occ_func_0_0(21) * occ_func_0_1(0) +
            0.7071067812 * occ_func_0_0(79) * occ_func_0_0(20) *
                occ_func_0_1(21) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(7) * occ_func_0_0(2) *
                occ_func_0_1(81) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(7) * occ_func_0_1(2) *
                occ_func_0_1(81))) +
          ((0.7071067812 * occ_func_0_0(6) * occ_func_0_1(0) * occ_func_0_0(1) *
                occ_func_0_1(23) +
            0.7071067812 * occ_func_0_0(6) * occ_func_0_0(0) * occ_func_0_1(1) *
                occ_func_0_1(23))) +
          ((0.7071067812 * occ_func_0_0(11) * occ_func_0_1(12) *
                occ_func_0_0(0) * occ_func_0_1(30) +
            0.7071067812 * occ_func_0_0(11) * occ_func_0_0(12) *
                occ_func_0_1(0) * occ_func_0_1(30))) +
          ((0.7071067812 * occ_func_0_0(84) * occ_func_0_1(38) *
                occ_func_0_0(31) * occ_func_0_1(0) +
            0.7071067812 * occ_func_0_0(84) * occ_func_0_0(38) *
                occ_func_0_1(31) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(6) * occ_func_0_0(8) *
                occ_func_0_1(84) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(6) * occ_func_0_1(8) *
                occ_func_0_1(84))) +
          ((0.7071067812 * occ_func_0_0(7) * occ_func_0_1(0) * occ_func_0_0(9) *
                occ_func_0_1(38) +
            0.7071067812 * occ_func_0_0(7) * occ_func_0_0(0) * occ_func_0_1(9) *
                occ_func_0_1(38))) +
          ((0.7071067812 * occ_func_0_0(5) * occ_func_0_1(4) * occ_func_0_0(0) *
                occ_func_0_1(36) +
            0.7071067812 * occ_func_0_0(5) * occ_func_0_0(4) * occ_func_0_1(0) *
                occ_func_0_1(36))) +
          ((0.7071067812 * occ_func_0_0(81) * occ_func_0_1(23) *
                occ_func_0_0(25) * occ_func_0_1(0) +
            0.7071067812 * occ_func_0_0(81) * occ_func_0_0(23) *
                occ_func_0_1(25) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(8) *
                occ_func_0_0(11) * occ_func_0_1(84) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(8) *
                occ_func_0_1(11) * occ_func_0_1(84))) +
          ((0.7071067812 * occ_func_0_0(5) * occ_func_0_1(0) *
                occ_func_0_0(10) * occ_func_0_1(36) +
            0.7071067812 * occ_func_0_0(5) * occ_func_0_0(0) *
                occ_func_0_1(10) * occ_func_0_1(36))) +
          ((0.7071067812 * occ_func_0_0(2) * occ_func_0_1(3) * occ_func_0_0(0) *
                occ_func_0_1(31) +
            0.7071067812 * occ_func_0_0(2) * occ_func_0_0(3) * occ_func_0_1(0) *
                occ_func_0_1(31))) +
          ((0.7071067812 * occ_func_0_0(81) * occ_func_0_1(25) *
                occ_func_0_0(30) * occ_func_0_1(0) +
            0.7071067812 * occ_func_0_0(81) * occ_func_0_0(25) *
                occ_func_0_1(30) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(12) *
                occ_func_0_0(9) * occ_func_0_1(85) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(12) *
                occ_func_0_1(9) * occ_func_0_1(85))) +
          ((0.7071067812 * occ_func_0_0(1) * occ_func_0_1(0) * occ_func_0_0(3) *
                occ_func_0_1(32) +
            0.7071067812 * occ_func_0_0(1) * occ_func_0_0(0) * occ_func_0_1(3) *
                occ_func_0_1(32))) +
          ((0.7071067812 * occ_func_0_0(4) * occ_func_0_1(10) *
                occ_func_0_0(0) * occ_func_0_1(37) +
            0.7071067812 * occ_func_0_0(4) * occ_func_0_0(10) *
                occ_func_0_1(0) * occ_func_0_1(37))) +
          ((0.7071067812 * occ_func_0_0(80) * occ_func_0_1(29) *
                occ_func_0_0(24) * occ_func_0_1(0) +
            0.7071067812 * occ_func_0_0(80) * occ_func_0_0(29) *
                occ_func_0_1(24) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(5) * occ_func_0_0(4) *
                occ_func_0_1(83) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(5) * occ_func_0_1(4) *
                occ_func_0_1(83))) +
          ((0.7071067812 * occ_func_0_0(8) * occ_func_0_1(0) * occ_func_0_0(6) *
                occ_func_0_1(34) +
            0.7071067812 * occ_func_0_0(8) * occ_func_0_0(0) * occ_func_0_1(6) *
                occ_func_0_1(34))) +
          ((0.7071067812 * occ_func_0_0(9) * occ_func_0_1(7) * occ_func_0_0(0) *
                occ_func_0_1(35) +
            0.7071067812 * occ_func_0_0(9) * occ_func_0_0(7) * occ_func_0_1(0) *
                occ_func_0_1(35))) +
          ((0.7071067812 * occ_func_0_0(82) * occ_func_0_1(27) *
                occ_func_0_0(26) * occ_func_0_1(0) +
            0.7071067812 * occ_func_0_0(82) * occ_func_0_0(27) *
                occ_func_0_1(26) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(7) *
                occ_func_0_0(12) * occ_func_0_1(85) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(7) *
                occ_func_0_1(12) * occ_func_0_1(85))) +
          ((0.7071067812 * occ_func_0_0(6) * occ_func_0_1(0) *
                occ_func_0_0(11) * occ_func_0_1(39) +
            0.7071067812 * occ_func_0_0(6) * occ_func_0_0(0) *
                occ_func_0_1(11) * occ_func_0_1(39))) +
          ((0.7071067812 * occ_func_0_0(1) * occ_func_0_1(2) * occ_func_0_0(0) *
                occ_func_0_1(32) +
            0.7071067812 * occ_func_0_0(1) * occ_func_0_0(2) * occ_func_0_1(0) *
                occ_func_0_1(32))) +
          ((0.7071067812 * occ_func_0_0(80) * occ_func_0_1(22) *
                occ_func_0_0(29) * occ_func_0_1(0) +
            0.7071067812 * occ_func_0_0(80) * occ_func_0_0(22) *
                occ_func_0_1(29) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(1) * occ_func_0_0(2) *
                occ_func_0_1(79) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(1) * occ_func_0_1(2) *
                occ_func_0_1(79))) +
          ((0.7071067812 * occ_func_0_0(12) * occ_func_0_1(0) *
                occ_func_0_0(7) * occ_func_0_1(21) +
            0.7071067812 * occ_func_0_0(12) * occ_func_0_0(0) *
                occ_func_0_1(7) * occ_func_0_1(21))) +
          ((0.7071067812 * occ_func_0_0(11) * occ_func_0_1(6) *
                occ_func_0_0(0) * occ_func_0_1(20) +
            0.7071067812 * occ_func_0_0(11) * occ_func_0_0(6) *
                occ_func_0_1(0) * occ_func_0_1(20))) +
          ((0.7071067812 * occ_func_0_0(86) * occ_func_0_1(40) *
                occ_func_0_0(41) * occ_func_0_1(0) +
            0.7071067812 * occ_func_0_0(86) * occ_func_0_0(40) *
                occ_func_0_1(41) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(3) * occ_func_0_0(1) *
                occ_func_0_1(79) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(3) * occ_func_0_1(1) *
                occ_func_0_1(79))) +
          ((0.7071067812 * occ_func_0_0(10) * occ_func_0_1(0) *
                occ_func_0_0(4) * occ_func_0_1(19) +
            0.7071067812 * occ_func_0_0(10) * occ_func_0_0(0) *
                occ_func_0_1(4) * occ_func_0_1(19))) +
          ((0.7071067812 * occ_func_0_0(12) * occ_func_0_1(9) *
                occ_func_0_0(0) * occ_func_0_1(21) +
            0.7071067812 * occ_func_0_0(12) * occ_func_0_0(9) *
                occ_func_0_1(0) * occ_func_0_1(21))) +
          ((0.7071067812 * occ_func_0_0(86) * occ_func_0_1(42) *
                occ_func_0_0(40) * occ_func_0_1(0) +
            0.7071067812 * occ_func_0_0(86) * occ_func_0_0(42) *
                occ_func_0_1(40) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(10) *
                occ_func_0_0(5) * occ_func_0_1(83) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(10) *
                occ_func_0_1(5) * occ_func_0_1(83))) +
          ((0.7071067812 * occ_func_0_0(3) * occ_func_0_1(0) * occ_func_0_0(2) *
                occ_func_0_1(28) +
            0.7071067812 * occ_func_0_0(3) * occ_func_0_0(0) * occ_func_0_1(2) *
                occ_func_0_1(28))) +
          ((0.7071067812 * occ_func_0_0(8) * occ_func_0_1(11) *
                occ_func_0_0(0) * occ_func_0_1(34) +
            0.7071067812 * occ_func_0_0(8) * occ_func_0_0(11) *
                occ_func_0_1(0) * occ_func_0_1(34))) +
          ((0.7071067812 * occ_func_0_0(82) * occ_func_0_1(33) *
                occ_func_0_0(27) * occ_func_0_1(0) +
            0.7071067812 * occ_func_0_0(82) * occ_func_0_0(33) *
                occ_func_0_1(27) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(10) *
                occ_func_0_0(11) * occ_func_0_1(86) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(10) *
                occ_func_0_1(11) * occ_func_0_1(86))) +
          ((0.7071067812 * occ_func_0_0(3) * occ_func_0_1(0) * occ_func_0_0(8) *
                occ_func_0_1(42) +
            0.7071067812 * occ_func_0_0(3) * occ_func_0_0(0) * occ_func_0_1(8) *
                occ_func_0_1(42))) +
          ((0.7071067812 * occ_func_0_0(2) * occ_func_0_1(5) * occ_func_0_0(0) *
                occ_func_0_1(41) +
            0.7071067812 * occ_func_0_0(2) * occ_func_0_0(5) * occ_func_0_1(0) *
                occ_func_0_1(41))) +
          ((0.7071067812 * occ_func_0_0(79) * occ_func_0_1(19) *
                occ_func_0_0(20) * occ_func_0_1(0) +
            0.7071067812 * occ_func_0_0(79) * occ_func_0_0(19) *
                occ_func_0_1(20) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(6) * occ_func_0_0(4) *
                occ_func_0_1(80) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(6) * occ_func_0_1(4) *
                occ_func_0_1(80))) +
          ((0.7071067812 * occ_func_0_0(7) * occ_func_0_1(0) * occ_func_0_0(5) *
                occ_func_0_1(22) +
            0.7071067812 * occ_func_0_0(7) * occ_func_0_0(0) * occ_func_0_1(5) *
                occ_func_0_1(22))) +
          ((0.7071067812 * occ_func_0_0(9) * occ_func_0_1(8) * occ_func_0_0(0) *
                occ_func_0_1(24) +
            0.7071067812 * occ_func_0_0(9) * occ_func_0_0(8) * occ_func_0_1(0) *
                occ_func_0_1(24))) +
          ((0.7071067812 * occ_func_0_0(85) * occ_func_0_1(39) *
                occ_func_0_0(37) * occ_func_0_1(0) +
            0.7071067812 * occ_func_0_0(85) * occ_func_0_0(39) *
                occ_func_0_1(37) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(9) * occ_func_0_0(8) *
                occ_func_0_1(82) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(9) * occ_func_0_1(8) *
                occ_func_0_1(82))) +
          ((0.7071067812 * occ_func_0_0(4) * occ_func_0_1(0) * occ_func_0_0(6) *
                occ_func_0_1(26) +
            0.7071067812 * occ_func_0_0(4) * occ_func_0_0(0) * occ_func_0_1(6) *
                occ_func_0_1(26))) +
          ((0.7071067812 * occ_func_0_0(5) * occ_func_0_1(7) * occ_func_0_0(0) *
                occ_func_0_1(27) +
            0.7071067812 * occ_func_0_0(5) * occ_func_0_0(7) * occ_func_0_1(0) *
                occ_func_0_1(27))) +
          ((0.7071067812 * occ_func_0_0(83) * occ_func_0_1(35) *
                occ_func_0_0(34) * occ_func_0_1(0) +
            0.7071067812 * occ_func_0_0(83) * occ_func_0_0(35) *
                occ_func_0_1(34) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(2) * occ_func_0_0(3) *
                occ_func_0_1(79) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(2) * occ_func_0_1(3) *
                occ_func_0_1(79))) +
          ((0.7071067812 * occ_func_0_0(11) * occ_func_0_1(0) *
                occ_func_0_0(8) * occ_func_0_1(20) +
            0.7071067812 * occ_func_0_0(11) * occ_func_0_0(0) *
                occ_func_0_1(8) * occ_func_0_1(20))) +
          ((0.7071067812 * occ_func_0_0(10) * occ_func_0_1(5) *
                occ_func_0_0(0) * occ_func_0_1(19) +
            0.7071067812 * occ_func_0_0(10) * occ_func_0_0(5) *
                occ_func_0_1(0) * occ_func_0_1(19))) +
          ((0.7071067812 * occ_func_0_0(86) * occ_func_0_1(41) *
                occ_func_0_0(42) * occ_func_0_1(0) +
            0.7071067812 * occ_func_0_0(86) * occ_func_0_0(41) *
                occ_func_0_1(42) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(11) *
                occ_func_0_0(6) * occ_func_0_1(84) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(11) *
                occ_func_0_1(6) * occ_func_0_1(84))) +
          ((0.7071067812 * occ_func_0_0(2) * occ_func_0_1(0) * occ_func_0_0(1) *
                occ_func_0_1(31) +
            0.7071067812 * occ_func_0_0(2) * occ_func_0_0(0) * occ_func_0_1(1) *
                occ_func_0_1(31))) +
          ((0.7071067812 * occ_func_0_0(7) * occ_func_0_1(12) *
                occ_func_0_0(0) * occ_func_0_1(38) +
            0.7071067812 * occ_func_0_0(7) * occ_func_0_0(12) *
                occ_func_0_1(0) * occ_func_0_1(38))) +
          ((0.7071067812 * occ_func_0_0(81) * occ_func_0_1(30) *
                occ_func_0_0(23) * occ_func_0_1(0) +
            0.7071067812 * occ_func_0_0(81) * occ_func_0_0(30) *
                occ_func_0_1(23) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(5) * occ_func_0_0(7) *
                occ_func_0_1(81) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(5) * occ_func_0_1(7) *
                occ_func_0_1(81))) +
          ((0.7071067812 * occ_func_0_0(8) * occ_func_0_1(0) * occ_func_0_0(9) *
                occ_func_0_1(25) +
            0.7071067812 * occ_func_0_0(8) * occ_func_0_0(0) * occ_func_0_1(9) *
                occ_func_0_1(25))) +
          ((0.7071067812 * occ_func_0_0(6) * occ_func_0_1(4) * occ_func_0_0(0) *
                occ_func_0_1(23) +
            0.7071067812 * occ_func_0_0(6) * occ_func_0_0(4) * occ_func_0_1(0) *
                occ_func_0_1(23))) +
          ((0.7071067812 * occ_func_0_0(84) * occ_func_0_1(36) *
                occ_func_0_0(38) * occ_func_0_1(0) +
            0.7071067812 * occ_func_0_0(84) * occ_func_0_0(36) *
                occ_func_0_1(38) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(1) * occ_func_0_0(6) *
                occ_func_0_1(80) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(1) * occ_func_0_1(6) *
                occ_func_0_1(80))) +
          ((0.7071067812 * occ_func_0_0(12) * occ_func_0_1(0) *
                occ_func_0_0(11) * occ_func_0_1(29) +
            0.7071067812 * occ_func_0_0(12) * occ_func_0_0(0) *
                occ_func_0_1(11) * occ_func_0_1(29))) +
          ((0.7071067812 * occ_func_0_0(7) * occ_func_0_1(2) * occ_func_0_0(0) *
                occ_func_0_1(22) +
            0.7071067812 * occ_func_0_0(7) * occ_func_0_0(2) * occ_func_0_1(0) *
                occ_func_0_1(22))) +
          ((0.7071067812 * occ_func_0_0(85) * occ_func_0_1(32) *
                occ_func_0_0(39) * occ_func_0_1(0) +
            0.7071067812 * occ_func_0_0(85) * occ_func_0_0(32) *
                occ_func_0_1(39) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(4) *
                occ_func_0_0(10) * occ_func_0_1(83) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(4) *
                occ_func_0_1(10) * occ_func_0_1(83))) +
          ((0.7071067812 * occ_func_0_0(9) * occ_func_0_1(0) *
                occ_func_0_0(12) * occ_func_0_1(35) +
            0.7071067812 * occ_func_0_0(9) * occ_func_0_0(0) *
                occ_func_0_1(12) * occ_func_0_1(35))) +
          ((0.7071067812 * occ_func_0_0(3) * occ_func_0_1(1) * occ_func_0_0(0) *
                occ_func_0_1(28) +
            0.7071067812 * occ_func_0_0(3) * occ_func_0_0(1) * occ_func_0_1(0) *
                occ_func_0_1(28))) +
          ((0.7071067812 * occ_func_0_0(82) * occ_func_0_1(26) *
                occ_func_0_0(33) * occ_func_0_1(0) +
            0.7071067812 * occ_func_0_0(82) * occ_func_0_0(26) *
                occ_func_0_1(33) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(4) * occ_func_0_0(1) *
                occ_func_0_1(80) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(4) * occ_func_0_1(1) *
                occ_func_0_1(80))) +
          ((0.7071067812 * occ_func_0_0(9) * occ_func_0_1(0) * occ_func_0_0(3) *
                occ_func_0_1(24) +
            0.7071067812 * occ_func_0_0(9) * occ_func_0_0(0) * occ_func_0_1(3) *
                occ_func_0_1(24))) +
          ((0.7071067812 * occ_func_0_0(12) * occ_func_0_1(10) *
                occ_func_0_0(0) * occ_func_0_1(29) +
            0.7071067812 * occ_func_0_0(12) * occ_func_0_0(10) *
                occ_func_0_1(0) * occ_func_0_1(29))) +
          ((0.7071067812 * occ_func_0_0(85) * occ_func_0_1(37) *
                occ_func_0_0(32) * occ_func_0_1(0) +
            0.7071067812 * occ_func_0_0(85) * occ_func_0_0(37) *
                occ_func_0_1(32) * occ_func_0_1(0)))) /
         24.0;
}
double test_Clexulator::site_eval_at_0_bfunc_4_0_9() const {
  return (((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(9) * occ_func_0_0(7) *
                occ_func_0_1(85) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(9) * occ_func_0_1(7) *
                occ_func_0_1(85))) +
          ((0.7071067812 * occ_func_0_1(4) * occ_func_0_1(0) * occ_func_0_0(5) *
                occ_func_0_1(37) +
            0.7071067812 * occ_func_0_1(4) * occ_func_0_0(0) * occ_func_0_1(5) *
                occ_func_0_1(37))) +
          ((0.7071067812 * occ_func_0_1(6) * occ_func_0_1(8) * occ_func_0_0(0) *
                occ_func_0_1(39) +
            0.7071067812 * occ_func_0_1(6) * occ_func_0_0(8) * occ_func_0_1(0) *
                occ_func_0_1(39))) +
          ((0.7071067812 * occ_func_0_1(80) * occ_func_0_1(24) *
                occ_func_0_0(22) * occ_func_0_1(0) +
            0.7071067812 * occ_func_0_1(80) * occ_func_0_0(24) *
                occ_func_0_1(22) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(3) * occ_func_0_0(9) *
                occ_func_0_1(82) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(3) * occ_func_0_1(9) *
                occ_func_0_1(82))) +
          ((0.7071067812 * occ_func_0_1(10) * occ_func_0_1(0) *
                occ_func_0_0(12) * occ_func_0_1(33) +
            0.7071067812 * occ_func_0_1(10) * occ_func_0_0(0) *
                occ_func_0_1(12) * occ_func_0_1(33))) +
          ((0.7071067812 * occ_func_0_1(4) * occ_func_0_1(1) * occ_func_0_0(0) *
                occ_func_0_1(26) +
            0.7071067812 * occ_func_0_1(4) * occ_func_0_0(1) * occ_func_0_1(0) *
                occ_func_0_1(26))) +
          ((0.7071067812 * occ_func_0_1(83) * occ_func_0_1(28) *
                occ_func_0_0(35) * occ_func_0_1(0) +
            0.7071067812 * occ_func_0_1(83) * occ_func_0_0(28) *
                occ_func_0_1(35) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(2) * occ_func_0_0(5) *
                occ_func_0_1(81) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(2) * occ_func_0_1(5) *
                occ_func_0_1(81))) +
          ((0.7071067812 * occ_func_0_1(11) * occ_func_0_1(0) *
                occ_func_0_0(10) * occ_func_0_1(30) +
            0.7071067812 * occ_func_0_1(11) * occ_func_0_0(0) *
                occ_func_0_1(10) * occ_func_0_1(30))) +
          ((0.7071067812 * occ_func_0_1(8) * occ_func_0_1(3) * occ_func_0_0(0) *
                occ_func_0_1(25) +
            0.7071067812 * occ_func_0_1(8) * occ_func_0_0(3) * occ_func_0_1(0) *
                occ_func_0_1(25))) +
          ((0.7071067812 * occ_func_0_1(84) * occ_func_0_1(31) *
                occ_func_0_0(36) * occ_func_0_1(0) +
            0.7071067812 * occ_func_0_1(84) * occ_func_0_0(31) *
                occ_func_0_1(36) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(12) *
                occ_func_0_0(10) * occ_func_0_1(86) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(12) *
                occ_func_0_1(10) * occ_func_0_1(86))) +
          ((0.7071067812 * occ_func_0_1(1) * occ_func_0_1(0) * occ_func_0_0(4) *
                occ_func_0_1(40) +
            0.7071067812 * occ_func_0_1(1) * occ_func_0_0(0) * occ_func_0_1(4) *
                occ_func_0_1(40))) +
          ((0.7071067812 * occ_func_0_1(3) * occ_func_0_1(9) * occ_func_0_0(0) *
                occ_func_0_1(42) +
            0.7071067812 * occ_func_0_1(3) * occ_func_0_0(9) * occ_func_0_1(0) *
                occ_func_0_1(42))) +
          ((0.7071067812 * occ_func_0_1(79) * occ_func_0_1(21) *
                occ_func_0_0(19) * occ_func_0_1(0) +
            0.7071067812 * occ_func_0_1(79) * occ_func_0_0(21) *
                occ_func_0_1(19) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(8) * occ_func_0_0(3) *
                occ_func_0_1(82) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(8) * occ_func_0_1(3) *
                occ_func_0_1(82))) +
          ((0.7071067812 * occ_func_0_1(5) * occ_func_0_1(0) * occ_func_0_0(2) *
                occ_func_0_1(27) +
            0.7071067812 * occ_func_0_1(5) * occ_func_0_0(0) * occ_func_0_1(2) *
                occ_func_0_1(27))) +
          ((0.7071067812 * occ_func_0_1(10) * occ_func_0_1(11) *
                occ_func_0_0(0) * occ_func_0_1(33) +
            0.7071067812 * occ_func_0_1(10) * occ_func_0_0(11) *
                occ_func_0_1(0) * occ_func_0_1(33))) +
          ((0.7071067812 * occ_func_0_1(83) * occ_func_0_1(34) *
                occ_func_0_0(28) * occ_func_0_1(0) +
            0.7071067812 * occ_func_0_1(83) * occ_func_0_0(34) *
                occ_func_0_1(28) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(11) *
                occ_func_0_0(12) * occ_func_0_1(86) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(11) *
                occ_func_0_1(12) * occ_func_0_1(86))) +
          ((0.7071067812 * occ_func_0_1(2) * occ_func_0_1(0) * occ_func_0_0(7) *
                occ_func_0_1(41) +
            0.7071067812 * occ_func_0_1(2) * occ_func_0_0(0) * occ_func_0_1(7) *
                occ_func_0_1(41))) +
          ((0.7071067812 * occ_func_0_1(1) * occ_func_0_1(6) * occ_func_0_0(0) *
                occ_func_0_1(40) +
            0.7071067812 * occ_func_0_1(1) * occ_func_0_0(6) * occ_func_0_1(0) *
                occ_func_0_1(40))) +
          ((0.7071067812 * occ_func_0_1(79) * occ_func_0_1(20) *
                occ_func_0_0(21) * occ_func_0_1(0) +
            0.7071067812 * occ_func_0_1(79) * occ_func_0_0(20) *
                occ_func_0_1(21) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(7) * occ_func_0_0(2) *
                occ_func_0_1(81) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(7) * occ_func_0_1(2) *
                occ_func_0_1(81))) +
          ((0.7071067812 * occ_func_0_1(6) * occ_func_0_1(0) * occ_func_0_0(1) *
                occ_func_0_1(23) +
            0.7071067812 * occ_func_0_1(6) * occ_func_0_0(0) * occ_func_0_1(1) *
                occ_func_0_1(23))) +
          ((0.7071067812 * occ_func_0_1(11) * occ_func_0_1(12) *
                occ_func_0_0(0) * occ_func_0_1(30) +
            0.7071067812 * occ_func_0_1(11) * occ_func_0_0(12) *
                occ_func_0_1(0) * occ_func_0_1(30))) +
          ((0.7071067812 * occ_func_0_1(84) * occ_func_0_1(38) *
                occ_func_0_0(31) * occ_func_0_1(0) +
            0.7071067812 * occ_func_0_1(84) * occ_func_0_0(38) *
                occ_func_0_1(31) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(6) * occ_func_0_0(8) *
                occ_func_0_1(84) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(6) * occ_func_0_1(8) *
                occ_func_0_1(84))) +
          ((0.7071067812 * occ_func_0_1(7) * occ_func_0_1(0) * occ_func_0_0(9) *
                occ_func_0_1(38) +
            0.7071067812 * occ_func_0_1(7) * occ_func_0_0(0) * occ_func_0_1(9) *
                occ_func_0_1(38))) +
          ((0.7071067812 * occ_func_0_1(5) * occ_func_0_1(4) * occ_func_0_0(0) *
                occ_func_0_1(36) +
            0.7071067812 * occ_func_0_1(5) * occ_func_0_0(4) * occ_func_0_1(0) *
                occ_func_0_1(36))) +
          ((0.7071067812 * occ_func_0_1(81) * occ_func_0_1(23) *
                occ_func_0_0(25) * occ_func_0_1(0) +
            0.7071067812 * occ_func_0_1(81) * occ_func_0_0(23) *
                occ_func_0_1(25) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(8) *
                occ_func_0_0(11) * occ_func_0_1(84) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(8) *
                occ_func_0_1(11) * occ_func_0_1(84))) +
          ((0.7071067812 * occ_func_0_1(5) * occ_func_0_1(0) *
                occ_func_0_0(10) * occ_func_0_1(36) +
            0.7071067812 * occ_func_0_1(5) * occ_func_0_0(0) *
                occ_func_0_1(10) * occ_func_0_1(36))) +
          ((0.7071067812 * occ_func_0_1(2) * occ_func_0_1(3) * occ_func_0_0(0) *
                occ_func_0_1(31) +
            0.7071067812 * occ_func_0_1(2) * occ_func_0_0(3) * occ_func_0_1(0) *
                occ_func_0_1(31))) +
          ((0.7071067812 * occ_func_0_1(81) * occ_func_0_1(25) *
                occ_func_0_0(30) * occ_func_0_1(0) +
            0.7071067812 * occ_func_0_1(81) * occ_func_0_0(25) *
                occ_func_0_1(30) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(12) *
                occ_func_0_0(9) * occ_func_0_1(85) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(12) *
                occ_func_0_1(9) * occ_func_0_1(85))) +
          ((0.7071067812 * occ_func_0_1(1) * occ_func_0_1(0) * occ_func_0_0(3) *
                occ_func_0_1(32) +
            0.7071067812 * occ_func_0_1(1) * occ_func_0_0(0) * occ_func_0_1(3) *
                occ_func_0_1(32))) +
          ((0.7071067812 * occ_func_0_1(4) * occ_func_0_1(10) *
                occ_func_0_0(0) * occ_func_0_1(37) +
            0.7071067812 * occ_func_0_1(4) * occ_func_0_0(10) *
                occ_func_0_1(0) * occ_func_0_1(37))) +
          ((0.7071067812 * occ_func_0_1(80) * occ_func_0_1(29) *
                occ_func_0_0(24) * occ_func_0_1(0) +
            0.7071067812 * occ_func_0_1(80) * occ_func_0_0(29) *
                occ_func_0_1(24) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(5) * occ_func_0_0(4) *
                occ_func_0_1(83) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(5) * occ_func_0_1(4) *
                occ_func_0_1(83))) +
          ((0.7071067812 * occ_func_0_1(8) * occ_func_0_1(0) * occ_func_0_0(6) *
                occ_func_0_1(34) +
            0.7071067812 * occ_func_0_1(8) * occ_func_0_0(0) * occ_func_0_1(6) *
                occ_func_0_1(34))) +
          ((0.7071067812 * occ_func_0_1(9) * occ_func_0_1(7) * occ_func_0_0(0) *
                occ_func_0_1(35) +
            0.7071067812 * occ_func_0_1(9) * occ_func_0_0(7) * occ_func_0_1(0) *
                occ_func_0_1(35))) +
          ((0.7071067812 * occ_func_0_1(82) * occ_func_0_1(27) *
                occ_func_0_0(26) * occ_func_0_1(0) +
            0.7071067812 * occ_func_0_1(82) * occ_func_0_0(27) *
                occ_func_0_1(26) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(7) *
                occ_func_0_0(12) * occ_func_0_1(85) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(7) *
                occ_func_0_1(12) * occ_func_0_1(85))) +
          ((0.7071067812 * occ_func_0_1(6) * occ_func_0_1(0) *
                occ_func_0_0(11) * occ_func_0_1(39) +
            0.7071067812 * occ_func_0_1(6) * occ_func_0_0(0) *
                occ_func_0_1(11) * occ_func_0_1(39))) +
          ((0.7071067812 * occ_func_0_1(1) * occ_func_0_1(2) * occ_func_0_0(0) *
                occ_func_0_1(32) +
            0.7071067812 * occ_func_0_1(1) * occ_func_0_0(2) * occ_func_0_1(0) *
                occ_func_0_1(32))) +
          ((0.7071067812 * occ_func_0_1(80) * occ_func_0_1(22) *
                occ_func_0_0(29) * occ_func_0_1(0) +
            0.7071067812 * occ_func_0_1(80) * occ_func_0_0(22) *
                occ_func_0_1(29) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(1) * occ_func_0_0(2) *
                occ_func_0_1(79) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(1) * occ_func_0_1(2) *
                occ_func_0_1(79))) +
          ((0.7071067812 * occ_func_0_1(12) * occ_func_0_1(0) *
                occ_func_0_0(7) * occ_func_0_1(21) +
            0.7071067812 * occ_func_0_1(12) * occ_func_0_0(0) *
                occ_func_0_1(7) * occ_func_0_1(21))) +
          ((0.7071067812 * occ_func_0_1(11) * occ_func_0_1(6) *
                occ_func_0_0(0) * occ_func_0_1(20) +
            0.7071067812 * occ_func_0_1(11) * occ_func_0_0(6) *
                occ_func_0_1(0) * occ_func_0_1(20))) +
          ((0.7071067812 * occ_func_0_1(86) * occ_func_0_1(40) *
                occ_func_0_0(41) * occ_func_0_1(0) +
            0.7071067812 * occ_func_0_1(86) * occ_func_0_0(40) *
                occ_func_0_1(41) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(3) * occ_func_0_0(1) *
                occ_func_0_1(79) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(3) * occ_func_0_1(1) *
                occ_func_0_1(79))) +
          ((0.7071067812 * occ_func_0_1(10) * occ_func_0_1(0) *
                occ_func_0_0(4) * occ_func_0_1(19) +
            0.7071067812 * occ_func_0_1(10) * occ_func_0_0(0) *
                occ_func_0_1(4) * occ_func_0_1(19))) +
          ((0.7071067812 * occ_func_0_1(12) * occ_func_0_1(9) *
                occ_func_0_0(0) * occ_func_0_1(21) +
            0.7071067812 * occ_func_0_1(12) * occ_func_0_0(9) *
                occ_func_0_1(0) * occ_func_0_1(21))) +
          ((0.7071067812 * occ_func_0_1(86) * occ_func_0_1(42) *
                occ_func_0_0(40) * occ_func_0_1(0) +
            0.7071067812 * occ_func_0_1(86) * occ_func_0_0(42) *
                occ_func_0_1(40) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(10) *
                occ_func_0_0(5) * occ_func_0_1(83) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(10) *
                occ_func_0_1(5) * occ_func_0_1(83))) +
          ((0.7071067812 * occ_func_0_1(3) * occ_func_0_1(0) * occ_func_0_0(2) *
                occ_func_0_1(28) +
            0.7071067812 * occ_func_0_1(3) * occ_func_0_0(0) * occ_func_0_1(2) *
                occ_func_0_1(28))) +
          ((0.7071067812 * occ_func_0_1(8) * occ_func_0_1(11) *
                occ_func_0_0(0) * occ_func_0_1(34) +
            0.7071067812 * occ_func_0_1(8) * occ_func_0_0(11) *
                occ_func_0_1(0) * occ_func_0_1(34))) +
          ((0.7071067812 * occ_func_0_1(82) * occ_func_0_1(33) *
                occ_func_0_0(27) * occ_func_0_1(0) +
            0.7071067812 * occ_func_0_1(82) * occ_func_0_0(33) *
                occ_func_0_1(27) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(10) *
                occ_func_0_0(11) * occ_func_0_1(86) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(10) *
                occ_func_0_1(11) * occ_func_0_1(86))) +
          ((0.7071067812 * occ_func_0_1(3) * occ_func_0_1(0) * occ_func_0_0(8) *
                occ_func_0_1(42) +
            0.7071067812 * occ_func_0_1(3) * occ_func_0_0(0) * occ_func_0_1(8) *
                occ_func_0_1(42))) +
          ((0.7071067812 * occ_func_0_1(2) * occ_func_0_1(5) * occ_func_0_0(0) *
                occ_func_0_1(41) +
            0.7071067812 * occ_func_0_1(2) * occ_func_0_0(5) * occ_func_0_1(0) *
                occ_func_0_1(41))) +
          ((0.7071067812 * occ_func_0_1(79) * occ_func_0_1(19) *
                occ_func_0_0(20) * occ_func_0_1(0) +
            0.7071067812 * occ_func_0_1(79) * occ_func_0_0(19) *
                occ_func_0_1(20) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(6) * occ_func_0_0(4) *
                occ_func_0_1(80) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(6) * occ_func_0_1(4) *
                occ_func_0_1(80))) +
          ((0.7071067812 * occ_func_0_1(7) * occ_func_0_1(0) * occ_func_0_0(5) *
                occ_func_0_1(22) +
            0.7071067812 * occ_func_0_1(7) * occ_func_0_0(0) * occ_func_0_1(5) *
                occ_func_0_1(22))) +
          ((0.7071067812 * occ_func_0_1(9) * occ_func_0_1(8) * occ_func_0_0(0) *
                occ_func_0_1(24) +
            0.7071067812 * occ_func_0_1(9) * occ_func_0_0(8) * occ_func_0_1(0) *
                occ_func_0_1(24))) +
          ((0.7071067812 * occ_func_0_1(85) * occ_func_0_1(39) *
                occ_func_0_0(37) * occ_func_0_1(0) +
            0.7071067812 * occ_func_0_1(85) * occ_func_0_0(39) *
                occ_func_0_1(37) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(9) * occ_func_0_0(8) *
                occ_func_0_1(82) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(9) * occ_func_0_1(8) *
                occ_func_0_1(82))) +
          ((0.7071067812 * occ_func_0_1(4) * occ_func_0_1(0) * occ_func_0_0(6) *
                occ_func_0_1(26) +
            0.7071067812 * occ_func_0_1(4) * occ_func_0_0(0) * occ_func_0_1(6) *
                occ_func_0_1(26))) +
          ((0.7071067812 * occ_func_0_1(5) * occ_func_0_1(7) * occ_func_0_0(0) *
                occ_func_0_1(27) +
            0.7071067812 * occ_func_0_1(5) * occ_func_0_0(7) * occ_func_0_1(0) *
                occ_func_0_1(27))) +
          ((0.7071067812 * occ_func_0_1(83) * occ_func_0_1(35) *
                occ_func_0_0(34) * occ_func_0_1(0) +
            0.7071067812 * occ_func_0_1(83) * occ_func_0_0(35) *
                occ_func_0_1(34) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(2) * occ_func_0_0(3) *
                occ_func_0_1(79) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(2) * occ_func_0_1(3) *
                occ_func_0_1(79))) +
          ((0.7071067812 * occ_func_0_1(11) * occ_func_0_1(0) *
                occ_func_0_0(8) * occ_func_0_1(20) +
            0.7071067812 * occ_func_0_1(11) * occ_func_0_0(0) *
                occ_func_0_1(8) * occ_func_0_1(20))) +
          ((0.7071067812 * occ_func_0_1(10) * occ_func_0_1(5) *
                occ_func_0_0(0) * occ_func_0_1(19) +
            0.7071067812 * occ_func_0_1(10) * occ_func_0_0(5) *
                occ_func_0_1(0) * occ_func_0_1(19))) +
          ((0.7071067812 * occ_func_0_1(86) * occ_func_0_1(41) *
                occ_func_0_0(42) * occ_func_0_1(0) +
            0.7071067812 * occ_func_0_1(86) * occ_func_0_0(41) *
                occ_func_0_1(42) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(11) *
                occ_func_0_0(6) * occ_func_0_1(84) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(11) *
                occ_func_0_1(6) * occ_func_0_1(84))) +
          ((0.7071067812 * occ_func_0_1(2) * occ_func_0_1(0) * occ_func_0_0(1) *
                occ_func_0_1(31) +
            0.7071067812 * occ_func_0_1(2) * occ_func_0_0(0) * occ_func_0_1(1) *
                occ_func_0_1(31))) +
          ((0.7071067812 * occ_func_0_1(7) * occ_func_0_1(12) *
                occ_func_0_0(0) * occ_func_0_1(38) +
            0.7071067812 * occ_func_0_1(7) * occ_func_0_0(12) *
                occ_func_0_1(0) * occ_func_0_1(38))) +
          ((0.7071067812 * occ_func_0_1(81) * occ_func_0_1(30) *
                occ_func_0_0(23) * occ_func_0_1(0) +
            0.7071067812 * occ_func_0_1(81) * occ_func_0_0(30) *
                occ_func_0_1(23) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(5) * occ_func_0_0(7) *
                occ_func_0_1(81) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(5) * occ_func_0_1(7) *
                occ_func_0_1(81))) +
          ((0.7071067812 * occ_func_0_1(8) * occ_func_0_1(0) * occ_func_0_0(9) *
                occ_func_0_1(25) +
            0.7071067812 * occ_func_0_1(8) * occ_func_0_0(0) * occ_func_0_1(9) *
                occ_func_0_1(25))) +
          ((0.7071067812 * occ_func_0_1(6) * occ_func_0_1(4) * occ_func_0_0(0) *
                occ_func_0_1(23) +
            0.7071067812 * occ_func_0_1(6) * occ_func_0_0(4) * occ_func_0_1(0) *
                occ_func_0_1(23))) +
          ((0.7071067812 * occ_func_0_1(84) * occ_func_0_1(36) *
                occ_func_0_0(38) * occ_func_0_1(0) +
            0.7071067812 * occ_func_0_1(84) * occ_func_0_0(36) *
                occ_func_0_1(38) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(1) * occ_func_0_0(6) *
                occ_func_0_1(80) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(1) * occ_func_0_1(6) *
                occ_func_0_1(80))) +
          ((0.7071067812 * occ_func_0_1(12) * occ_func_0_1(0) *
                occ_func_0_0(11) * occ_func_0_1(29) +
            0.7071067812 * occ_func_0_1(12) * occ_func_0_0(0) *
                occ_func_0_1(11) * occ_func_0_1(29))) +
          ((0.7071067812 * occ_func_0_1(7) * occ_func_0_1(2) * occ_func_0_0(0) *
                occ_func_0_1(22) +
            0.7071067812 * occ_func_0_1(7) * occ_func_0_0(2) * occ_func_0_1(0) *
                occ_func_0_1(22))) +
          ((0.7071067812 * occ_func_0_1(85) * occ_func_0_1(32) *
                occ_func_0_0(39) * occ_func_0_1(0) +
            0.7071067812 * occ_func_0_1(85) * occ_func_0_0(32) *
                occ_func_0_1(39) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(4) *
                occ_func_0_0(10) * occ_func_0_1(83) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(4) *
                occ_func_0_1(10) * occ_func_0_1(83))) +
          ((0.7071067812 * occ_func_0_1(9) * occ_func_0_1(0) *
                occ_func_0_0(12) * occ_func_0_1(35) +
            0.7071067812 * occ_func_0_1(9) * occ_func_0_0(0) *
                occ_func_0_1(12) * occ_func_0_1(35))) +
          ((0.7071067812 * occ_func_0_1(3) * occ_func_0_1(1) * occ_func_0_0(0) *
                occ_func_0_1(28) +
            0.7071067812 * occ_func_0_1(3) * occ_func_0_0(1) * occ_func_0_1(0) *
                occ_func_0_1(28))) +
          ((0.7071067812 * occ_func_0_1(82) * occ_func_0_1(26) *
                occ_func_0_0(33) * occ_func_0_1(0) +
            0.7071067812 * occ_func_0_1(82) * occ_func_0_0(26) *
                occ_func_0_1(33) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(4) * occ_func_0_0(1) *
                occ_func_0_1(80) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(4) * occ_func_0_1(1) *
                occ_func_0_1(80))) +
          ((0.7071067812 * occ_func_0_1(9) * occ_func_0_1(0) * occ_func_0_0(3) *
                occ_func_0_1(24) +
            0.7071067812 * occ_func_0_1(9) * occ_func_0_0(0) * occ_func_0_1(3) *
                occ_func_0_1(24))) +
          ((0.7071067812 * occ_func_0_1(12) * occ_func_0_1(10) *
                occ_func_0_0(0) * occ_func_0_1(29) +
            0.7071067812 * occ_func_0_1(12) * occ_func_0_0(10) *
                occ_func_0_1(0) * occ_func_0_1(29))) +
          ((0.7071067812 * occ_func_0_1(85) * occ_func_0_1(37) *
                occ_func_0_0(32) * occ_func_0_1(0) +
            0.7071067812 * occ_func_0_1(85) * occ_func_0_0(37) *
                occ_func_0_1(32) * occ_func_0_1(0)))) /
         24.0;
}
double test_Clexulator::site_eval_at_0_bfunc_4_0_10() const {
  return ((occ_func_0_0(0) * occ_func_0_1(9) * occ_func_0_1(7) *
           occ_func_0_1(85)) +
          (occ_func_0_0(4) * occ_func_0_1(0) * occ_func_0_1(5) *
           occ_func_0_1(37)) +
          (occ_func_0_0(6) * occ_func_0_1(8) * occ_func_0_1(0) *
           occ_func_0_1(39)) +
          (occ_func_0_0(80) * occ_func_0_1(24) * occ_func_0_1(22) *
           occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_1(3) * occ_func_0_1(9) *
           occ_func_0_1(82)) +
          (occ_func_0_0(10) * occ_func_0_1(0) * occ_func_0_1(12) *
           occ_func_0_1(33)) +
          (occ_func_0_0(4) * occ_func_0_1(1) * occ_func_0_1(0) *
           occ_func_0_1(26)) +
          (occ_func_0_0(83) * occ_func_0_1(28) * occ_func_0_1(35) *
           occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_1(2) * occ_func_0_1(5) *
           occ_func_0_1(81)) +
          (occ_func_0_0(11) * occ_func_0_1(0) * occ_func_0_1(10) *
           occ_func_0_1(30)) +
          (occ_func_0_0(8) * occ_func_0_1(3) * occ_func_0_1(0) *
           occ_func_0_1(25)) +
          (occ_func_0_0(84) * occ_func_0_1(31) * occ_func_0_1(36) *
           occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_1(12) * occ_func_0_1(10) *
           occ_func_0_1(86)) +
          (occ_func_0_0(1) * occ_func_0_1(0) * occ_func_0_1(4) *
           occ_func_0_1(40)) +
          (occ_func_0_0(3) * occ_func_0_1(9) * occ_func_0_1(0) *
           occ_func_0_1(42)) +
          (occ_func_0_0(79) * occ_func_0_1(21) * occ_func_0_1(19) *
           occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_1(8) * occ_func_0_1(3) *
           occ_func_0_1(82)) +
          (occ_func_0_0(5) * occ_func_0_1(0) * occ_func_0_1(2) *
           occ_func_0_1(27)) +
          (occ_func_0_0(10) * occ_func_0_1(11) * occ_func_0_1(0) *
           occ_func_0_1(33)) +
          (occ_func_0_0(83) * occ_func_0_1(34) * occ_func_0_1(28) *
           occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_1(11) * occ_func_0_1(12) *
           occ_func_0_1(86)) +
          (occ_func_0_0(2) * occ_func_0_1(0) * occ_func_0_1(7) *
           occ_func_0_1(41)) +
          (occ_func_0_0(1) * occ_func_0_1(6) * occ_func_0_1(0) *
           occ_func_0_1(40)) +
          (occ_func_0_0(79) * occ_func_0_1(20) * occ_func_0_1(21) *
           occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_1(7) * occ_func_0_1(2) *
           occ_func_0_1(81)) +
          (occ_func_0_0(6) * occ_func_0_1(0) * occ_func_0_1(1) *
           occ_func_0_1(23)) +
          (occ_func_0_0(11) * occ_func_0_1(12) * occ_func_0_1(0) *
           occ_func_0_1(30)) +
          (occ_func_0_0(84) * occ_func_0_1(38) * occ_func_0_1(31) *
           occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_1(6) * occ_func_0_1(8) *
           occ_func_0_1(84)) +
          (occ_func_0_0(7) * occ_func_0_1(0) * occ_func_0_1(9) *
           occ_func_0_1(38)) +
          (occ_func_0_0(5) * occ_func_0_1(4) * occ_func_0_1(0) *
           occ_func_0_1(36)) +
          (occ_func_0_0(81) * occ_func_0_1(23) * occ_func_0_1(25) *
           occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_1(8) * occ_func_0_1(11) *
           occ_func_0_1(84)) +
          (occ_func_0_0(5) * occ_func_0_1(0) * occ_func_0_1(10) *
           occ_func_0_1(36)) +
          (occ_func_0_0(2) * occ_func_0_1(3) * occ_func_0_1(0) *
           occ_func_0_1(31)) +
          (occ_func_0_0(81) * occ_func_0_1(25) * occ_func_0_1(30) *
           occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_1(12) * occ_func_0_1(9) *
           occ_func_0_1(85)) +
          (occ_func_0_0(1) * occ_func_0_1(0) * occ_func_0_1(3) *
           occ_func_0_1(32)) +
          (occ_func_0_0(4) * occ_func_0_1(10) * occ_func_0_1(0) *
           occ_func_0_1(37)) +
          (occ_func_0_0(80) * occ_func_0_1(29) * occ_func_0_1(24) *
           occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_1(5) * occ_func_0_1(4) *
           occ_func_0_1(83)) +
          (occ_func_0_0(8) * occ_func_0_1(0) * occ_func_0_1(6) *
           occ_func_0_1(34)) +
          (occ_func_0_0(9) * occ_func_0_1(7) * occ_func_0_1(0) *
           occ_func_0_1(35)) +
          (occ_func_0_0(82) * occ_func_0_1(27) * occ_func_0_1(26) *
           occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_1(7) * occ_func_0_1(12) *
           occ_func_0_1(85)) +
          (occ_func_0_0(6) * occ_func_0_1(0) * occ_func_0_1(11) *
           occ_func_0_1(39)) +
          (occ_func_0_0(1) * occ_func_0_1(2) * occ_func_0_1(0) *
           occ_func_0_1(32)) +
          (occ_func_0_0(80) * occ_func_0_1(22) * occ_func_0_1(29) *
           occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_1(1) * occ_func_0_1(2) *
           occ_func_0_1(79)) +
          (occ_func_0_0(12) * occ_func_0_1(0) * occ_func_0_1(7) *
           occ_func_0_1(21)) +
          (occ_func_0_0(11) * occ_func_0_1(6) * occ_func_0_1(0) *
           occ_func_0_1(20)) +
          (occ_func_0_0(86) * occ_func_0_1(40) * occ_func_0_1(41) *
           occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_1(3) * occ_func_0_1(1) *
           occ_func_0_1(79)) +
          (occ_func_0_0(10) * occ_func_0_1(0) * occ_func_0_1(4) *
           occ_func_0_1(19)) +
          (occ_func_0_0(12) * occ_func_0_1(9) * occ_func_0_1(0) *
           occ_func_0_1(21)) +
          (occ_func_0_0(86) * occ_func_0_1(42) * occ_func_0_1(40) *
           occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_1(10) * occ_func_0_1(5) *
           occ_func_0_1(83)) +
          (occ_func_0_0(3) * occ_func_0_1(0) * occ_func_0_1(2) *
           occ_func_0_1(28)) +
          (occ_func_0_0(8) * occ_func_0_1(11) * occ_func_0_1(0) *
           occ_func_0_1(34)) +
          (occ_func_0_0(82) * occ_func_0_1(33) * occ_func_0_1(27) *
           occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_1(10) * occ_func_0_1(11) *
           occ_func_0_1(86)) +
          (occ_func_0_0(3) * occ_func_0_1(0) * occ_func_0_1(8) *
           occ_func_0_1(42)) +
          (occ_func_0_0(2) * occ_func_0_1(5) * occ_func_0_1(0) *
           occ_func_0_1(41)) +
          (occ_func_0_0(79) * occ_func_0_1(19) * occ_func_0_1(20) *
           occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_1(6) * occ_func_0_1(4) *
           occ_func_0_1(80)) +
          (occ_func_0_0(7) * occ_func_0_1(0) * occ_func_0_1(5) *
           occ_func_0_1(22)) +
          (occ_func_0_0(9) * occ_func_0_1(8) * occ_func_0_1(0) *
           occ_func_0_1(24)) +
          (occ_func_0_0(85) * occ_func_0_1(39) * occ_func_0_1(37) *
           occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_1(9) * occ_func_0_1(8) *
           occ_func_0_1(82)) +
          (occ_func_0_0(4) * occ_func_0_1(0) * occ_func_0_1(6) *
           occ_func_0_1(26)) +
          (occ_func_0_0(5) * occ_func_0_1(7) * occ_func_0_1(0) *
           occ_func_0_1(27)) +
          (occ_func_0_0(83) * occ_func_0_1(35) * occ_func_0_1(34) *
           occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_1(2) * occ_func_0_1(3) *
           occ_func_0_1(79)) +
          (occ_func_0_0(11) * occ_func_0_1(0) * occ_func_0_1(8) *
           occ_func_0_1(20)) +
          (occ_func_0_0(10) * occ_func_0_1(5) * occ_func_0_1(0) *
           occ_func_0_1(19)) +
          (occ_func_0_0(86) * occ_func_0_1(41) * occ_func_0_1(42) *
           occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_1(11) * occ_func_0_1(6) *
           occ_func_0_1(84)) +
          (occ_func_0_0(2) * occ_func_0_1(0) * occ_func_0_1(1) *
           occ_func_0_1(31)) +
          (occ_func_0_0(7) * occ_func_0_1(12) * occ_func_0_1(0) *
           occ_func_0_1(38)) +
          (occ_func_0_0(81) * occ_func_0_1(30) * occ_func_0_1(23) *
           occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_1(5) * occ_func_0_1(7) *
           occ_func_0_1(81)) +
          (occ_func_0_0(8) * occ_func_0_1(0) * occ_func_0_1(9) *
           occ_func_0_1(25)) +
          (occ_func_0_0(6) * occ_func_0_1(4) * occ_func_0_1(0) *
           occ_func_0_1(23)) +
          (occ_func_0_0(84) * occ_func_0_1(36) * occ_func_0_1(38) *
           occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_1(1) * occ_func_0_1(6) *
           occ_func_0_1(80)) +
          (occ_func_0_0(12) * occ_func_0_1(0) * occ_func_0_1(11) *
           occ_func_0_1(29)) +
          (occ_func_0_0(7) * occ_func_0_1(2) * occ_func_0_1(0) *
           occ_func_0_1(22)) +
          (occ_func_0_0(85) * occ_func_0_1(32) * occ_func_0_1(39) *
           occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_1(4) * occ_func_0_1(10) *
           occ_func_0_1(83)) +
          (occ_func_0_0(9) * occ_func_0_1(0) * occ_func_0_1(12) *
           occ_func_0_1(35)) +
          (occ_func_0_0(3) * occ_func_0_1(1) * occ_func_0_1(0) *
           occ_func_0_1(28)) +
          (occ_func_0_0(82) * occ_func_0_1(26) * occ_func_0_1(33) *
           occ_func_0_1(0)) +
          (occ_func_0_0(0) * occ_func_0_1(4) * occ_func_0_1(1) *
           occ_func_0_1(80)) +
          (occ_func_0_0(9) * occ_func_0_1(0) * occ_func_0_1(3) *
           occ_func_0_1(24)) +
          (occ_func_0_0(12) * occ_func_0_1(10) * occ_func_0_1(0) *
           occ_func_0_1(29)) +
          (occ_func_0_0(85) * occ_func_0_1(37) * occ_func_0_1(32) *
           occ_func_0_1(0))) /
         24.0;
}
double test_Clexulator::site_eval_at_0_bfunc_4_0_11() const {
  return ((occ_func_0_1(0) * occ_func_0_1(9) * occ_func_0_1(7) *
           occ_func_0_1(85)) +
          (occ_func_0_1(4) * occ_func_0_1(0) * occ_func_0_1(5) *
           occ_func_0_1(37)) +
          (occ_func_0_1(6) * occ_func_0_1(8) * occ_func_0_1(0) *
           occ_func_0_1(39)) +
          (occ_func_0_1(80) * occ_func_0_1(24) * occ_func_0_1(22) *
           occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(3) * occ_func_0_1(9) *
           occ_func_0_1(82)) +
          (occ_func_0_1(10) * occ_func_0_1(0) * occ_func_0_1(12) *
           occ_func_0_1(33)) +
          (occ_func_0_1(4) * occ_func_0_1(1) * occ_func_0_1(0) *
           occ_func_0_1(26)) +
          (occ_func_0_1(83) * occ_func_0_1(28) * occ_func_0_1(35) *
           occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(2) * occ_func_0_1(5) *
           occ_func_0_1(81)) +
          (occ_func_0_1(11) * occ_func_0_1(0) * occ_func_0_1(10) *
           occ_func_0_1(30)) +
          (occ_func_0_1(8) * occ_func_0_1(3) * occ_func_0_1(0) *
           occ_func_0_1(25)) +
          (occ_func_0_1(84) * occ_func_0_1(31) * occ_func_0_1(36) *
           occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(12) * occ_func_0_1(10) *
           occ_func_0_1(86)) +
          (occ_func_0_1(1) * occ_func_0_1(0) * occ_func_0_1(4) *
           occ_func_0_1(40)) +
          (occ_func_0_1(3) * occ_func_0_1(9) * occ_func_0_1(0) *
           occ_func_0_1(42)) +
          (occ_func_0_1(79) * occ_func_0_1(21) * occ_func_0_1(19) *
           occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(8) * occ_func_0_1(3) *
           occ_func_0_1(82)) +
          (occ_func_0_1(5) * occ_func_0_1(0) * occ_func_0_1(2) *
           occ_func_0_1(27)) +
          (occ_func_0_1(10) * occ_func_0_1(11) * occ_func_0_1(0) *
           occ_func_0_1(33)) +
          (occ_func_0_1(83) * occ_func_0_1(34) * occ_func_0_1(28) *
           occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(11) * occ_func_0_1(12) *
           occ_func_0_1(86)) +
          (occ_func_0_1(2) * occ_func_0_1(0) * occ_func_0_1(7) *
           occ_func_0_1(41)) +
          (occ_func_0_1(1) * occ_func_0_1(6) * occ_func_0_1(0) *
           occ_func_0_1(40)) +
          (occ_func_0_1(79) * occ_func_0_1(20) * occ_func_0_1(21) *
           occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(7) * occ_func_0_1(2) *
           occ_func_0_1(81)) +
          (occ_func_0_1(6) * occ_func_0_1(0) * occ_func_0_1(1) *
           occ_func_0_1(23)) +
          (occ_func_0_1(11) * occ_func_0_1(12) * occ_func_0_1(0) *
           occ_func_0_1(30)) +
          (occ_func_0_1(84) * occ_func_0_1(38) * occ_func_0_1(31) *
           occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(6) * occ_func_0_1(8) *
           occ_func_0_1(84)) +
          (occ_func_0_1(7) * occ_func_0_1(0) * occ_func_0_1(9) *
           occ_func_0_1(38)) +
          (occ_func_0_1(5) * occ_func_0_1(4) * occ_func_0_1(0) *
           occ_func_0_1(36)) +
          (occ_func_0_1(81) * occ_func_0_1(23) * occ_func_0_1(25) *
           occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(8) * occ_func_0_1(11) *
           occ_func_0_1(84)) +
          (occ_func_0_1(5) * occ_func_0_1(0) * occ_func_0_1(10) *
           occ_func_0_1(36)) +
          (occ_func_0_1(2) * occ_func_0_1(3) * occ_func_0_1(0) *
           occ_func_0_1(31)) +
          (occ_func_0_1(81) * occ_func_0_1(25) * occ_func_0_1(30) *
           occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(12) * occ_func_0_1(9) *
           occ_func_0_1(85)) +
          (occ_func_0_1(1) * occ_func_0_1(0) * occ_func_0_1(3) *
           occ_func_0_1(32)) +
          (occ_func_0_1(4) * occ_func_0_1(10) * occ_func_0_1(0) *
           occ_func_0_1(37)) +
          (occ_func_0_1(80) * occ_func_0_1(29) * occ_func_0_1(24) *
           occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(5) * occ_func_0_1(4) *
           occ_func_0_1(83)) +
          (occ_func_0_1(8) * occ_func_0_1(0) * occ_func_0_1(6) *
           occ_func_0_1(34)) +
          (occ_func_0_1(9) * occ_func_0_1(7) * occ_func_0_1(0) *
           occ_func_0_1(35)) +
          (occ_func_0_1(82) * occ_func_0_1(27) * occ_func_0_1(26) *
           occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(7) * occ_func_0_1(12) *
           occ_func_0_1(85)) +
          (occ_func_0_1(6) * occ_func_0_1(0) * occ_func_0_1(11) *
           occ_func_0_1(39)) +
          (occ_func_0_1(1) * occ_func_0_1(2) * occ_func_0_1(0) *
           occ_func_0_1(32)) +
          (occ_func_0_1(80) * occ_func_0_1(22) * occ_func_0_1(29) *
           occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(1) * occ_func_0_1(2) *
           occ_func_0_1(79)) +
          (occ_func_0_1(12) * occ_func_0_1(0) * occ_func_0_1(7) *
           occ_func_0_1(21)) +
          (occ_func_0_1(11) * occ_func_0_1(6) * occ_func_0_1(0) *
           occ_func_0_1(20)) +
          (occ_func_0_1(86) * occ_func_0_1(40) * occ_func_0_1(41) *
           occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(3) * occ_func_0_1(1) *
           occ_func_0_1(79)) +
          (occ_func_0_1(10) * occ_func_0_1(0) * occ_func_0_1(4) *
           occ_func_0_1(19)) +
          (occ_func_0_1(12) * occ_func_0_1(9) * occ_func_0_1(0) *
           occ_func_0_1(21)) +
          (occ_func_0_1(86) * occ_func_0_1(42) * occ_func_0_1(40) *
           occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(10) * occ_func_0_1(5) *
           occ_func_0_1(83)) +
          (occ_func_0_1(3) * occ_func_0_1(0) * occ_func_0_1(2) *
           occ_func_0_1(28)) +
          (occ_func_0_1(8) * occ_func_0_1(11) * occ_func_0_1(0) *
           occ_func_0_1(34)) +
          (occ_func_0_1(82) * occ_func_0_1(33) * occ_func_0_1(27) *
           occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(10) * occ_func_0_1(11) *
           occ_func_0_1(86)) +
          (occ_func_0_1(3) * occ_func_0_1(0) * occ_func_0_1(8) *
           occ_func_0_1(42)) +
          (occ_func_0_1(2) * occ_func_0_1(5) * occ_func_0_1(0) *
           occ_func_0_1(41)) +
          (occ_func_0_1(79) * occ_func_0_1(19) * occ_func_0_1(20) *
           occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(6) * occ_func_0_1(4) *
           occ_func_0_1(80)) +
          (occ_func_0_1(7) * occ_func_0_1(0) * occ_func_0_1(5) *
           occ_func_0_1(22)) +
          (occ_func_0_1(9) * occ_func_0_1(8) * occ_func_0_1(0) *
           occ_func_0_1(24)) +
          (occ_func_0_1(85) * occ_func_0_1(39) * occ_func_0_1(37) *
           occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(9) * occ_func_0_1(8) *
           occ_func_0_1(82)) +
          (occ_func_0_1(4) * occ_func_0_1(0) * occ_func_0_1(6) *
           occ_func_0_1(26)) +
          (occ_func_0_1(5) * occ_func_0_1(7) * occ_func_0_1(0) *
           occ_func_0_1(27)) +
          (occ_func_0_1(83) * occ_func_0_1(35) * occ_func_0_1(34) *
           occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(2) * occ_func_0_1(3) *
           occ_func_0_1(79)) +
          (occ_func_0_1(11) * occ_func_0_1(0) * occ_func_0_1(8) *
           occ_func_0_1(20)) +
          (occ_func_0_1(10) * occ_func_0_1(5) * occ_func_0_1(0) *
           occ_func_0_1(19)) +
          (occ_func_0_1(86) * occ_func_0_1(41) * occ_func_0_1(42) *
           occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(11) * occ_func_0_1(6) *
           occ_func_0_1(84)) +
          (occ_func_0_1(2) * occ_func_0_1(0) * occ_func_0_1(1) *
           occ_func_0_1(31)) +
          (occ_func_0_1(7) * occ_func_0_1(12) * occ_func_0_1(0) *
           occ_func_0_1(38)) +
          (occ_func_0_1(81) * occ_func_0_1(30) * occ_func_0_1(23) *
           occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(5) * occ_func_0_1(7) *
           occ_func_0_1(81)) +
          (occ_func_0_1(8) * occ_func_0_1(0) * occ_func_0_1(9) *
           occ_func_0_1(25)) +
          (occ_func_0_1(6) * occ_func_0_1(4) * occ_func_0_1(0) *
           occ_func_0_1(23)) +
          (occ_func_0_1(84) * occ_func_0_1(36) * occ_func_0_1(38) *
           occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(1) * occ_func_0_1(6) *
           occ_func_0_1(80)) +
          (occ_func_0_1(12) * occ_func_0_1(0) * occ_func_0_1(11) *
           occ_func_0_1(29)) +
          (occ_func_0_1(7) * occ_func_0_1(2) * occ_func_0_1(0) *
           occ_func_0_1(22)) +
          (occ_func_0_1(85) * occ_func_0_1(32) * occ_func_0_1(39) *
           occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(4) * occ_func_0_1(10) *
           occ_func_0_1(83)) +
          (occ_func_0_1(9) * occ_func_0_1(0) * occ_func_0_1(12) *
           occ_func_0_1(35)) +
          (occ_func_0_1(3) * occ_func_0_1(1) * occ_func_0_1(0) *
           occ_func_0_1(28)) +
          (occ_func_0_1(82) * occ_func_0_1(26) * occ_func_0_1(33) *
           occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(4) * occ_func_0_1(1) *
           occ_func_0_1(80)) +
          (occ_func_0_1(9) * occ_func_0_1(0) * occ_func_0_1(3) *
           occ_func_0_1(24)) +
          (occ_func_0_1(12) * occ_func_0_1(10) * occ_func_0_1(0) *
           occ_func_0_1(29)) +
          (occ_func_0_1(85) * occ_func_0_1(37) * occ_func_0_1(32) *
           occ_func_0_1(0))) /
         24.0;
}

double test_Clexulator::delta_site_eval_at_0_bfunc_4_0_0(int occ_i,
                                                         int occ_f) const {
  return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) *
         ((occ_func_0_0(9) * occ_func_0_0(7) * occ_func_0_0(85)) +
          (occ_func_0_0(4) * occ_func_0_0(5) * occ_func_0_0(37)) +
          (occ_func_0_0(6) * occ_func_0_0(8) * occ_func_0_0(39)) +
          (occ_func_0_0(80) * occ_func_0_0(24) * occ_func_0_0(22)) +
          (occ_func_0_0(3) * occ_func_0_0(9) * occ_func_0_0(82)) +
          (occ_func_0_0(10) * occ_func_0_0(12) * occ_func_0_0(33)) +
          (occ_func_0_0(4) * occ_func_0_0(1) * occ_func_0_0(26)) +
          (occ_func_0_0(83) * occ_func_0_0(28) * occ_func_0_0(35)) +
          (occ_func_0_0(2) * occ_func_0_0(5) * occ_func_0_0(81)) +
          (occ_func_0_0(11) * occ_func_0_0(10) * occ_func_0_0(30)) +
          (occ_func_0_0(8) * occ_func_0_0(3) * occ_func_0_0(25)) +
          (occ_func_0_0(84) * occ_func_0_0(31) * occ_func_0_0(36)) +
          (occ_func_0_0(12) * occ_func_0_0(10) * occ_func_0_0(86)) +
          (occ_func_0_0(1) * occ_func_0_0(4) * occ_func_0_0(40)) +
          (occ_func_0_0(3) * occ_func_0_0(9) * occ_func_0_0(42)) +
          (occ_func_0_0(79) * occ_func_0_0(21) * occ_func_0_0(19)) +
          (occ_func_0_0(8) * occ_func_0_0(3) * occ_func_0_0(82)) +
          (occ_func_0_0(5) * occ_func_0_0(2) * occ_func_0_0(27)) +
          (occ_func_0_0(10) * occ_func_0_0(11) * occ_func_0_0(33)) +
          (occ_func_0_0(83) * occ_func_0_0(34) * occ_func_0_0(28)) +
          (occ_func_0_0(11) * occ_func_0_0(12) * occ_func_0_0(86)) +
          (occ_func_0_0(2) * occ_func_0_0(7) * occ_func_0_0(41)) +
          (occ_func_0_0(1) * occ_func_0_0(6) * occ_func_0_0(40)) +
          (occ_func_0_0(79) * occ_func_0_0(20) * occ_func_0_0(21)) +
          (occ_func_0_0(7) * occ_func_0_0(2) * occ_func_0_0(81)) +
          (occ_func_0_0(6) * occ_func_0_0(1) * occ_func_0_0(23)) +
          (occ_func_0_0(11) * occ_func_0_0(12) * occ_func_0_0(30)) +
          (occ_func_0_0(84) * occ_func_0_0(38) * occ_func_0_0(31)) +
          (occ_func_0_0(6) * occ_func_0_0(8) * occ_func_0_0(84)) +
          (occ_func_0_0(7) * occ_func_0_0(9) * occ_func_0_0(38)) +
          (occ_func_0_0(5) * occ_func_0_0(4) * occ_func_0_0(36)) +
          (occ_func_0_0(81) * occ_func_0_0(23) * occ_func_0_0(25)) +
          (occ_func_0_0(8) * occ_func_0_0(11) * occ_func_0_0(84)) +
          (occ_func_0_0(5) * occ_func_0_0(10) * occ_func_0_0(36)) +
          (occ_func_0_0(2) * occ_func_0_0(3) * occ_func_0_0(31)) +
          (occ_func_0_0(81) * occ_func_0_0(25) * occ_func_0_0(30)) +
          (occ_func_0_0(12) * occ_func_0_0(9) * occ_func_0_0(85)) +
          (occ_func_0_0(1) * occ_func_0_0(3) * occ_func_0_0(32)) +
          (occ_func_0_0(4) * occ_func_0_0(10) * occ_func_0_0(37)) +
          (occ_func_0_0(80) * occ_func_0_0(29) * occ_func_0_0(24)) +
          (occ_func_0_0(5) * occ_func_0_0(4) * occ_func_0_0(83)) +
          (occ_func_0_0(8) * occ_func_0_0(6) * occ_func_0_0(34)) +
          (occ_func_0_0(9) * occ_func_0_0(7) * occ_func_0_0(35)) +
          (occ_func_0_0(82) * occ_func_0_0(27) * occ_func_0_0(26)) +
          (occ_func_0_0(7) * occ_func_0_0(12) * occ_func_0_0(85)) +
          (occ_func_0_0(6) * occ_func_0_0(11) * occ_func_0_0(39)) +
          (occ_func_0_0(1) * occ_func_0_0(2) * occ_func_0_0(32)) +
          (occ_func_0_0(80) * occ_func_0_0(22) * occ_func_0_0(29)) +
          (occ_func_0_0(1) * occ_func_0_0(2) * occ_func_0_0(79)) +
          (occ_func_0_0(12) * occ_func_0_0(7) * occ_func_0_0(21)) +
          (occ_func_0_0(11) * occ_func_0_0(6) * occ_func_0_0(20)) +
          (occ_func_0_0(86) * occ_func_0_0(40) * occ_func_0_0(41)) +
          (occ_func_0_0(3) * occ_func_0_0(1) * occ_func_0_0(79)) +
          (occ_func_0_0(10) * occ_func_0_0(4) * occ_func_0_0(19)) +
          (occ_func_0_0(12) * occ_func_0_0(9) * occ_func_0_0(21)) +
          (occ_func_0_0(86) * occ_func_0_0(42) * occ_func_0_0(40)) +
          (occ_func_0_0(10) * occ_func_0_0(5) * occ_func_0_0(83)) +
          (occ_func_0_0(3) * occ_func_0_0(2) * occ_func_0_0(28)) +
          (occ_func_0_0(8) * occ_func_0_0(11) * occ_func_0_0(34)) +
          (occ_func_0_0(82) * occ_func_0_0(33) * occ_func_0_0(27)) +
          (occ_func_0_0(10) * occ_func_0_0(11) * occ_func_0_0(86)) +
          (occ_func_0_0(3) * occ_func_0_0(8) * occ_func_0_0(42)) +
          (occ_func_0_0(2) * occ_func_0_0(5) * occ_func_0_0(41)) +
          (occ_func_0_0(79) * occ_func_0_0(19) * occ_func_0_0(20)) +
          (occ_func_0_0(6) * occ_func_0_0(4) * occ_func_0_0(80)) +
          (occ_func_0_0(7) * occ_func_0_0(5) * occ_func_0_0(22)) +
          (occ_func_0_0(9) * occ_func_0_0(8) * occ_func_0_0(24)) +
          (occ_func_0_0(85) * occ_func_0_0(39) * occ_func_0_0(37)) +
          (occ_func_0_0(9) * occ_func_0_0(8) * occ_func_0_0(82)) +
          (occ_func_0_0(4) * occ_func_0_0(6) * occ_func_0_0(26)) +
          (occ_func_0_0(5) * occ_func_0_0(7) * occ_func_0_0(27)) +
          (occ_func_0_0(83) * occ_func_0_0(35) * occ_func_0_0(34)) +
          (occ_func_0_0(2) * occ_func_0_0(3) * occ_func_0_0(79)) +
          (occ_func_0_0(11) * occ_func_0_0(8) * occ_func_0_0(20)) +
          (occ_func_0_0(10) * occ_func_0_0(5) * occ_func_0_0(19)) +
          (occ_func_0_0(86) * occ_func_0_0(41) * occ_func_0_0(42)) +
          (occ_func_0_0(11) * occ_func_0_0(6) * occ_func_0_0(84)) +
          (occ_func_0_0(2) * occ_func_0_0(1) * occ_func_0_0(31)) +
          (occ_func_0_0(7) * occ_func_0_0(12) * occ_func_0_0(38)) +
          (occ_func_0_0(81) * occ_func_0_0(30) * occ_func_0_0(23)) +
          (occ_func_0_0(5) * occ_func_0_0(7) * occ_func_0_0(81)) +
          (occ_func_0_0(8) * occ_func_0_0(9) * occ_func_0_0(25)) +
          (occ_func_0_0(6) * occ_func_0_0(4) * occ_func_0_0(23)) +
          (occ_func_0_0(84) * occ_func_0_0(36) * occ_func_0_0(38)) +
          (occ_func_0_0(1) * occ_func_0_0(6) * occ_func_0_0(80)) +
          (occ_func_0_0(12) * occ_func_0_0(11) * occ_func_0_0(29)) +
          (occ_func_0_0(7) * occ_func_0_0(2) * occ_func_0_0(22)) +
          (occ_func_0_0(85) * occ_func_0_0(32) * occ_func_0_0(39)) +
          (occ_func_0_0(4) * occ_func_0_0(10) * occ_func_0_0(83)) +
          (occ_func_0_0(9) * occ_func_0_0(12) * occ_func_0_0(35)) +
          (occ_func_0_0(3) * occ_func_0_0(1) * occ_func_0_0(28)) +
          (occ_func_0_0(82) * occ_func_0_0(26) * occ_func_0_0(33)) +
          (occ_func_0_0(4) * occ_func_0_0(1) * occ_func_0_0(80)) +
          (occ_func_0_0(9) * occ_func_0_0(3) * occ_func_0_0(24)) +
          (occ_func_0_0(12) * occ_func_0_0(10) * occ_func_0_0(29)) +
          (occ_func_0_0(85) * occ_func_0_0(37) * occ_func_0_0(32))) /
         24.0;
}
double test_Clexulator::delta_site_eval_at_0_bfunc_4_0_1(int occ_i,
                                                         int occ_f) const {
  return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) *
             ((occ_func_0_1(4) * occ_func_0_0(5) * occ_func_0_0(37)) +
              (occ_func_0_1(6) * occ_func_0_0(8) * occ_func_0_0(39)) +
              (occ_func_0_1(80) * occ_func_0_0(24) * occ_func_0_0(22)) +
              (occ_func_0_1(10) * occ_func_0_0(12) * occ_func_0_0(33)) +
              (occ_func_0_1(4) * occ_func_0_0(1) * occ_func_0_0(26)) +
              (occ_func_0_1(83) * occ_func_0_0(28) * occ_func_0_0(35)) +
              (occ_func_0_1(11) * occ_func_0_0(10) * occ_func_0_0(30)) +
              (occ_func_0_1(8) * occ_func_0_0(3) * occ_func_0_0(25)) +
              (occ_func_0_1(84) * occ_func_0_0(31) * occ_func_0_0(36)) +
              (occ_func_0_1(1) * occ_func_0_0(4) * occ_func_0_0(40)) +
              (occ_func_0_1(3) * occ_func_0_0(9) * occ_func_0_0(42)) +
              (occ_func_0_1(79) * occ_func_0_0(21) * occ_func_0_0(19)) +
              (occ_func_0_1(5) * occ_func_0_0(2) * occ_func_0_0(27)) +
              (occ_func_0_1(10) * occ_func_0_0(11) * occ_func_0_0(33)) +
              (occ_func_0_1(83) * occ_func_0_0(34) * occ_func_0_0(28)) +
              (occ_func_0_1(2) * occ_func_0_0(7) * occ_func_0_0(41)) +
              (occ_func_0_1(1) * occ_func_0_0(6) * occ_func_0_0(40)) +
              (occ_func_0_1(79) * occ_func_0_0(20) * occ_func_0_0(21)) +
              (occ_func_0_1(6) * occ_func_0_0(1) * occ_func_0_0(23)) +
              (occ_func_0_1(11) * occ_func_0_0(12) * occ_func_0_0(30)) +
              (occ_func_0_1(84) * occ_func_0_0(38) * occ_func_0_0(31)) +
              (occ_func_0_1(7) * occ_func_0_0(9) * occ_func_0_0(38)) +
              (occ_func_0_1(5) * occ_func_0_0(4) * occ_func_0_0(36)) +
              (occ_func_0_1(81) * occ_func_0_0(23) * occ_func_0_0(25)) +
              (occ_func_0_1(5) * occ_func_0_0(10) * occ_func_0_0(36)) +
              (occ_func_0_1(2) * occ_func_0_0(3) * occ_func_0_0(31)) +
              (occ_func_0_1(81) * occ_func_0_0(25) * occ_func_0_0(30)) +
              (occ_func_0_1(1) * occ_func_0_0(3) * occ_func_0_0(32)) +
              (occ_func_0_1(4) * occ_func_0_0(10) * occ_func_0_0(37)) +
              (occ_func_0_1(80) * occ_func_0_0(29) * occ_func_0_0(24)) +
              (occ_func_0_1(8) * occ_func_0_0(6) * occ_func_0_0(34)) +
              (occ_func_0_1(9) * occ_func_0_0(7) * occ_func_0_0(35)) +
              (occ_func_0_1(82) * occ_func_0_0(27) * occ_func_0_0(26)) +
              (occ_func_0_1(6) * occ_func_0_0(11) * occ_func_0_0(39)) +
              (occ_func_0_1(1) * occ_func_0_0(2) * occ_func_0_0(32)) +
              (occ_func_0_1(80) * occ_func_0_0(22) * occ_func_0_0(29)) +
              (occ_func_0_1(12) * occ_func_0_0(7) * occ_func_0_0(21)) +
              (occ_func_0_1(11) * occ_func_0_0(6) * occ_func_0_0(20)) +
              (occ_func_0_1(86) * occ_func_0_0(40) * occ_func_0_0(41)) +
              (occ_func_0_1(10) * occ_func_0_0(4) * occ_func_0_0(19)) +
              (occ_func_0_1(12) * occ_func_0_0(9) * occ_func_0_0(21)) +
              (occ_func_0_1(86) * occ_func_0_0(42) * occ_func_0_0(40)) +
              (occ_func_0_1(3) * occ_func_0_0(2) * occ_func_0_0(28)) +
              (occ_func_0_1(8) * occ_func_0_0(11) * occ_func_0_0(34)) +
              (occ_func_0_1(82) * occ_func_0_0(33) * occ_func_0_0(27)) +
              (occ_func_0_1(3) * occ_func_0_0(8) * occ_func_0_0(42)) +
              (occ_func_0_1(2) * occ_func_0_0(5) * occ_func_0_0(41)) +
              (occ_func_0_1(79) * occ_func_0_0(19) * occ_func_0_0(20)) +
              (occ_func_0_1(7) * occ_func_0_0(5) * occ_func_0_0(22)) +
              (occ_func_0_1(9) * occ_func_0_0(8) * occ_func_0_0(24)) +
              (occ_func_0_1(85) * occ_func_0_0(39) * occ_func_0_0(37)) +
              (occ_func_0_1(4) * occ_func_0_0(6) * occ_func_0_0(26)) +
              (occ_func_0_1(5) * occ_func_0_0(7) * occ_func_0_0(27)) +
              (occ_func_0_1(83) * occ_func_0_0(35) * occ_func_0_0(34)) +
              (occ_func_0_1(11) * occ_func_0_0(8) * occ_func_0_0(20)) +
              (occ_func_0_1(10) * occ_func_0_0(5) * occ_func_0_0(19)) +
              (occ_func_0_1(86) * occ_func_0_0(41) * occ_func_0_0(42)) +
              (occ_func_0_1(2) * occ_func_0_0(1) * occ_func_0_0(31)) +
              (occ_func_0_1(7) * occ_func_0_0(12) * occ_func_0_0(38)) +
              (occ_func_0_1(81) * occ_func_0_0(30) * occ_func_0_0(23)) +
              (occ_func_0_1(8) * occ_func_0_0(9) * occ_func_0_0(25)) +
              (occ_func_0_1(6) * occ_func_0_0(4) * occ_func_0_0(23)) +
              (occ_func_0_1(84) * occ_func_0_0(36) * occ_func_0_0(38)) +
              (occ_func_0_1(12) * occ_func_0_0(11) * occ_func_0_0(29)) +
              (occ_func_0_1(7) * occ_func_0_0(2) * occ_func_0_0(22)) +
              (occ_func_0_1(85) * occ_func_0_0(32) * occ_func_0_0(39)) +
              (occ_func_0_1(9) * occ_func_0_0(12) * occ_func_0_0(35)) +
              (occ_func_0_1(3) * occ_func_0_0(1) * occ_func_0_0(28)) +
              (occ_func_0_1(82) * occ_func_0_0(26) * occ_func_0_0(33)) +
              (occ_func_0_1(9) * occ_func_0_0(3) * occ_func_0_0(24)) +
              (occ_func_0_1(12) * occ_func_0_0(10) * occ_func_0_0(29)) +
              (occ_func_0_1(85) * occ_func_0_0(37) * occ_func_0_0(32))) /
             24.0 +
         (m_occ_func_0_1[occ_f] - m_occ_func_0_1[occ_i]) *
             ((occ_func_0_0(9) * occ_func_0_0(7) * occ_func_0_0(85)) +
              (occ_func_0_0(3) * occ_func_0_0(9) * occ_func_0_0(82)) +
              (occ_func_0_0(2) * occ_func_0_0(5) * occ_func_0_0(81)) +
              (occ_func_0_0(12) * occ_func_0_0(10) * occ_func_0_0(86)) +
              (occ_func_0_0(8) * occ_func_0_0(3) * occ_func_0_0(82)) +
              (occ_func_0_0(11) * occ_func_0_0(12) * occ_func_0_0(86)) +
              (occ_func_0_0(7) * occ_func_0_0(2) * occ_func_0_0(81)) +
              (occ_func_0_0(6) * occ_func_0_0(8) * occ_func_0_0(84)) +
              (occ_func_0_0(8) * occ_func_0_0(11) * occ_func_0_0(84)) +
              (occ_func_0_0(12) * occ_func_0_0(9) * occ_func_0_0(85)) +
              (occ_func_0_0(5) * occ_func_0_0(4) * occ_func_0_0(83)) +
              (occ_func_0_0(7) * occ_func_0_0(12) * occ_func_0_0(85)) +
              (occ_func_0_0(1) * occ_func_0_0(2) * occ_func_0_0(79)) +
              (occ_func_0_0(3) * occ_func_0_0(1) * occ_func_0_0(79)) +
              (occ_func_0_0(10) * occ_func_0_0(5) * occ_func_0_0(83)) +
              (occ_func_0_0(10) * occ_func_0_0(11) * occ_func_0_0(86)) +
              (occ_func_0_0(6) * occ_func_0_0(4) * occ_func_0_0(80)) +
              (occ_func_0_0(9) * occ_func_0_0(8) * occ_func_0_0(82)) +
              (occ_func_0_0(2) * occ_func_0_0(3) * occ_func_0_0(79)) +
              (occ_func_0_0(11) * occ_func_0_0(6) * occ_func_0_0(84)) +
              (occ_func_0_0(5) * occ_func_0_0(7) * occ_func_0_0(81)) +
              (occ_func_0_0(1) * occ_func_0_0(6) * occ_func_0_0(80)) +
              (occ_func_0_0(4) * occ_func_0_0(10) * occ_func_0_0(83)) +
              (occ_func_0_0(4) * occ_func_0_0(1) * occ_func_0_0(80))) /
             24.0;
}
double test_Clexulator::delta_site_eval_at_0_bfunc_4_0_2(int occ_i,
                                                         int occ_f) const {
  return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) *
             (((0.7071067812 * occ_func_0_1(9) * occ_func_0_0(7) *
                    occ_func_0_0(85) +
                0.7071067812 * occ_func_0_0(9) * occ_func_0_1(7) *
                    occ_func_0_0(85))) +
              (0.7071067812 * occ_func_0_0(4) * occ_func_0_1(5) *
               occ_func_0_0(37)) +
              (0.7071067812 * occ_func_0_0(6) * occ_func_0_1(8) *
               occ_func_0_0(39)) +
              ((0.7071067812 * occ_func_0_0(80) * occ_func_0_1(24) *
                    occ_func_0_0(22) +
                0.7071067812 * occ_func_0_0(80) * occ_func_0_0(24) *
                    occ_func_0_1(22))) +
              ((0.7071067812 * occ_func_0_1(3) * occ_func_0_0(9) *
                    occ_func_0_0(82) +
                0.7071067812 * occ_func_0_0(3) * occ_func_0_1(9) *
                    occ_func_0_0(82))) +
              (0.7071067812 * occ_func_0_0(10) * occ_func_0_1(12) *
               occ_func_0_0(33)) +
              (0.7071067812 * occ_func_0_0(4) * occ_func_0_1(1) *
               occ_func_0_0(26)) +
              ((0.7071067812 * occ_func_0_0(83) * occ_func_0_1(28) *
                    occ_func_0_0(35) +
                0.7071067812 * occ_func_0_0(83) * occ_func_0_0(28) *
                    occ_func_0_1(35))) +
              ((0.7071067812 * occ_func_0_1(2) * occ_func_0_0(5) *
                    occ_func_0_0(81) +
                0.7071067812 * occ_func_0_0(2) * occ_func_0_1(5) *
                    occ_func_0_0(81))) +
              (0.7071067812 * occ_func_0_0(11) * occ_func_0_1(10) *
               occ_func_0_0(30)) +
              (0.7071067812 * occ_func_0_0(8) * occ_func_0_1(3) *
               occ_func_0_0(25)) +
              ((0.7071067812 * occ_func_0_0(84) * occ_func_0_1(31) *
                    occ_func_0_0(36) +
                0.7071067812 * occ_func_0_0(84) * occ_func_0_0(31) *
                    occ_func_0_1(36))) +
              ((0.7071067812 * occ_func_0_1(12) * occ_func_0_0(10) *
                    occ_func_0_0(86) +
                0.7071067812 * occ_func_0_0(12) * occ_func_0_1(10) *
                    occ_func_0_0(86))) +
              (0.7071067812 * occ_func_0_0(1) * occ_func_0_1(4) *
               occ_func_0_0(40)) +
              (0.7071067812 * occ_func_0_0(3) * occ_func_0_1(9) *
               occ_func_0_0(42)) +
              ((0.7071067812 * occ_func_0_0(79) * occ_func_0_1(21) *
                    occ_func_0_0(19) +
                0.7071067812 * occ_func_0_0(79) * occ_func_0_0(21) *
                    occ_func_0_1(19))) +
              ((0.7071067812 * occ_func_0_1(8) * occ_func_0_0(3) *
                    occ_func_0_0(82) +
                0.7071067812 * occ_func_0_0(8) * occ_func_0_1(3) *
                    occ_func_0_0(82))) +
              (0.7071067812 * occ_func_0_0(5) * occ_func_0_1(2) *
               occ_func_0_0(27)) +
              (0.7071067812 * occ_func_0_0(10) * occ_func_0_1(11) *
               occ_func_0_0(33)) +
              ((0.7071067812 * occ_func_0_0(83) * occ_func_0_1(34) *
                    occ_func_0_0(28) +
                0.7071067812 * occ_func_0_0(83) * occ_func_0_0(34) *
                    occ_func_0_1(28))) +
              ((0.7071067812 * occ_func_0_1(11) * occ_func_0_0(12) *
                    occ_func_0_0(86) +
                0.7071067812 * occ_func_0_0(11) * occ_func_0_1(12) *
                    occ_func_0_0(86))) +
              (0.7071067812 * occ_func_0_0(2) * occ_func_0_1(7) *
               occ_func_0_0(41)) +
              (0.7071067812 * occ_func_0_0(1) * occ_func_0_1(6) *
               occ_func_0_0(40)) +
              ((0.7071067812 * occ_func_0_0(79) * occ_func_0_1(20) *
                    occ_func_0_0(21) +
                0.7071067812 * occ_func_0_0(79) * occ_func_0_0(20) *
                    occ_func_0_1(21))) +
              ((0.7071067812 * occ_func_0_1(7) * occ_func_0_0(2) *
                    occ_func_0_0(81) +
                0.7071067812 * occ_func_0_0(7) * occ_func_0_1(2) *
                    occ_func_0_0(81))) +
              (0.7071067812 * occ_func_0_0(6) * occ_func_0_1(1) *
               occ_func_0_0(23)) +
              (0.7071067812 * occ_func_0_0(11) * occ_func_0_1(12) *
               occ_func_0_0(30)) +
              ((0.7071067812 * occ_func_0_0(84) * occ_func_0_1(38) *
                    occ_func_0_0(31) +
                0.7071067812 * occ_func_0_0(84) * occ_func_0_0(38) *
                    occ_func_0_1(31))) +
              ((0.7071067812 * occ_func_0_1(6) * occ_func_0_0(8) *
                    occ_func_0_0(84) +
                0.7071067812 * occ_func_0_0(6) * occ_func_0_1(8) *
                    occ_func_0_0(84))) +
              (0.7071067812 * occ_func_0_0(7) * occ_func_0_1(9) *
               occ_func_0_0(38)) +
              (0.7071067812 * occ_func_0_0(5) * occ_func_0_1(4) *
               occ_func_0_0(36)) +
              ((0.7071067812 * occ_func_0_0(81) * occ_func_0_1(23) *
                    occ_func_0_0(25) +
                0.7071067812 * occ_func_0_0(81) * occ_func_0_0(23) *
                    occ_func_0_1(25))) +
              ((0.7071067812 * occ_func_0_1(8) * occ_func_0_0(11) *
                    occ_func_0_0(84) +
                0.7071067812 * occ_func_0_0(8) * occ_func_0_1(11) *
                    occ_func_0_0(84))) +
              (0.7071067812 * occ_func_0_0(5) * occ_func_0_1(10) *
               occ_func_0_0(36)) +
              (0.7071067812 * occ_func_0_0(2) * occ_func_0_1(3) *
               occ_func_0_0(31)) +
              ((0.7071067812 * occ_func_0_0(81) * occ_func_0_1(25) *
                    occ_func_0_0(30) +
                0.7071067812 * occ_func_0_0(81) * occ_func_0_0(25) *
                    occ_func_0_1(30))) +
              ((0.7071067812 * occ_func_0_1(12) * occ_func_0_0(9) *
                    occ_func_0_0(85) +
                0.7071067812 * occ_func_0_0(12) * occ_func_0_1(9) *
                    occ_func_0_0(85))) +
              (0.7071067812 * occ_func_0_0(1) * occ_func_0_1(3) *
               occ_func_0_0(32)) +
              (0.7071067812 * occ_func_0_0(4) * occ_func_0_1(10) *
               occ_func_0_0(37)) +
              ((0.7071067812 * occ_func_0_0(80) * occ_func_0_1(29) *
                    occ_func_0_0(24) +
                0.7071067812 * occ_func_0_0(80) * occ_func_0_0(29) *
                    occ_func_0_1(24))) +
              ((0.7071067812 * occ_func_0_1(5) * occ_func_0_0(4) *
                    occ_func_0_0(83) +
                0.7071067812 * occ_func_0_0(5) * occ_func_0_1(4) *
                    occ_func_0_0(83))) +
              (0.7071067812 * occ_func_0_0(8) * occ_func_0_1(6) *
               occ_func_0_0(34)) +
              (0.7071067812 * occ_func_0_0(9) * occ_func_0_1(7) *
               occ_func_0_0(35)) +
              ((0.7071067812 * occ_func_0_0(82) * occ_func_0_1(27) *
                    occ_func_0_0(26) +
                0.7071067812 * occ_func_0_0(82) * occ_func_0_0(27) *
                    occ_func_0_1(26))) +
              ((0.7071067812 * occ_func_0_1(7) * occ_func_0_0(12) *
                    occ_func_0_0(85) +
                0.7071067812 * occ_func_0_0(7) * occ_func_0_1(12) *
                    occ_func_0_0(85))) +
              (0.7071067812 * occ_func_0_0(6) * occ_func_0_1(11) *
               occ_func_0_0(39)) +
              (0.7071067812 * occ_func_0_0(1) * occ_func_0_1(2) *
               occ_func_0_0(32)) +
              ((0.7071067812 * occ_func_0_0(80) * occ_func_0_1(22) *
                    occ_func_0_0(29) +
                0.7071067812 * occ_func_0_0(80) * occ_func_0_0(22) *
                    occ_func_0_1(29))) +
              ((0.7071067812 * occ_func_0_1(1) * occ_func_0_0(2) *
                    occ_func_0_0(79) +
                0.7071067812 * occ_func_0_0(1) * occ_func_0_1(2) *
                    occ_func_0_0(79))) +
              (0.7071067812 * occ_func_0_0(12) * occ_func_0_1(7) *
               occ_func_0_0(21)) +
              (0.7071067812 * occ_func_0_0(11) * occ_func_0_1(6) *
               occ_func_0_0(20)) +
              ((0.7071067812 * occ_func_0_0(86) * occ_func_0_1(40) *
                    occ_func_0_0(41) +
                0.7071067812 * occ_func_0_0(86) * occ_func_0_0(40) *
                    occ_func_0_1(41))) +
              ((0.7071067812 * occ_func_0_1(3) * occ_func_0_0(1) *
                    occ_func_0_0(79) +
                0.7071067812 * occ_func_0_0(3) * occ_func_0_1(1) *
                    occ_func_0_0(79))) +
              (0.7071067812 * occ_func_0_0(10) * occ_func_0_1(4) *
               occ_func_0_0(19)) +
              (0.7071067812 * occ_func_0_0(12) * occ_func_0_1(9) *
               occ_func_0_0(21)) +
              ((0.7071067812 * occ_func_0_0(86) * occ_func_0_1(42) *
                    occ_func_0_0(40) +
                0.7071067812 * occ_func_0_0(86) * occ_func_0_0(42) *
                    occ_func_0_1(40))) +
              ((0.7071067812 * occ_func_0_1(10) * occ_func_0_0(5) *
                    occ_func_0_0(83) +
                0.7071067812 * occ_func_0_0(10) * occ_func_0_1(5) *
                    occ_func_0_0(83))) +
              (0.7071067812 * occ_func_0_0(3) * occ_func_0_1(2) *
               occ_func_0_0(28)) +
              (0.7071067812 * occ_func_0_0(8) * occ_func_0_1(11) *
               occ_func_0_0(34)) +
              ((0.7071067812 * occ_func_0_0(82) * occ_func_0_1(33) *
                    occ_func_0_0(27) +
                0.7071067812 * occ_func_0_0(82) * occ_func_0_0(33) *
                    occ_func_0_1(27))) +
              ((0.7071067812 * occ_func_0_1(10) * occ_func_0_0(11) *
                    occ_func_0_0(86) +
                0.7071067812 * occ_func_0_0(10) * occ_func_0_1(11) *
                    occ_func_0_0(86))) +
              (0.7071067812 * occ_func_0_0(3) * occ_func_0_1(8) *
               occ_func_0_0(42)) +
              (0.7071067812 * occ_func_0_0(2) * occ_func_0_1(5) *
               occ_func_0_0(41)) +
              ((0.7071067812 * occ_func_0_0(79) * occ_func_0_1(19) *
                    occ_func_0_0(20) +
                0.7071067812 * occ_func_0_0(79) * occ_func_0_0(19) *
                    occ_func_0_1(20))) +
              ((0.7071067812 * occ_func_0_1(6) * occ_func_0_0(4) *
                    occ_func_0_0(80) +
                0.7071067812 * occ_func_0_0(6) * occ_func_0_1(4) *
                    occ_func_0_0(80))) +
              (0.7071067812 * occ_func_0_0(7) * occ_func_0_1(5) *
               occ_func_0_0(22)) +
              (0.7071067812 * occ_func_0_0(9) * occ_func_0_1(8) *
               occ_func_0_0(24)) +
              ((0.7071067812 * occ_func_0_0(85) * occ_func_0_1(39) *
                    occ_func_0_0(37) +
                0.7071067812 * occ_func_0_0(85) * occ_func_0_0(39) *
                    occ_func_0_1(37))) +
              ((0.7071067812 * occ_func_0_1(9) * occ_func_0_0(8) *
                    occ_func_0_0(82) +
                0.7071067812 * occ_func_0_0(9) * occ_func_0_1(8) *
                    occ_func_0_0(82))) +
              (0.7071067812 * occ_func_0_0(4) * occ_func_0_1(6) *
               occ_func_0_0(26)) +
              (0.7071067812 * occ_func_0_0(5) * occ_func_0_1(7) *
               occ_func_0_0(27)) +
              ((0.7071067812 * occ_func_0_0(83) * occ_func_0_1(35) *
                    occ_func_0_0(34) +
                0.7071067812 * occ_func_0_0(83) * occ_func_0_0(35) *
                    occ_func_0_1(34))) +
              ((0.7071067812 * occ_func_0_1(2) * occ_func_0_0(3) *
                    occ_func_0_0(79) +
                0.7071067812 * occ_func_0_0(2) * occ_func_0_1(3) *
                    occ_func_0_0(79))) +
              (0.7071067812 * occ_func_0_0(11) * occ_func_0_1(8) *
               occ_func_0_0(20)) +
              (0.7071067812 * occ_func_0_0(10) * occ_func_0_1(5) *
               occ_func_0_0(19)) +
              ((0.7071067812 * occ_func_0_0(86) * occ_func_0_1(41) *
                    occ_func_0_0(42) +
                0.7071067812 * occ_func_0_0(86) * occ_func_0_0(41) *
                    occ_func_0_1(42))) +
              ((0.7071067812 * occ_func_0_1(11) * occ_func_0_0(6) *
                    occ_func_0_0(84) +
                0.7071067812 * occ_func_0_0(11) * occ_func_0_1(6) *
                    occ_func_0_0(84))) +
              (0.7071067812 * occ_func_0_0(2) * occ_func_0_1(1) *
               occ_func_0_0(31)) +
              (0.7071067812 * occ_func_0_0(7) * occ_func_0_1(12) *
               occ_func_0_0(38)) +
              ((0.7071067812 * occ_func_0_0(81) * occ_func_0_1(30) *
                    occ_func_0_0(23) +
                0.7071067812 * occ_func_0_0(81) * occ_func_0_0(30) *
                    occ_func_0_1(23))) +
              ((0.7071067812 * occ_func_0_1(5) * occ_func_0_0(7) *
                    occ_func_0_0(81) +
                0.7071067812 * occ_func_0_0(5) * occ_func_0_1(7) *
                    occ_func_0_0(81))) +
              (0.7071067812 * occ_func_0_0(8) * occ_func_0_1(9) *
               occ_func_0_0(25)) +
              (0.7071067812 * occ_func_0_0(6) * occ_func_0_1(4) *
               occ_func_0_0(23)) +
              ((0.7071067812 * occ_func_0_0(84) * occ_func_0_1(36) *
                    occ_func_0_0(38) +
                0.7071067812 * occ_func_0_0(84) * occ_func_0_0(36) *
                    occ_func_0_1(38))) +
              ((0.7071067812 * occ_func_0_1(1) * occ_func_0_0(6) *
                    occ_func_0_0(80) +
                0.7071067812 * occ_func_0_0(1) * occ_func_0_1(6) *
                    occ_func_0_0(80))) +
              (0.7071067812 * occ_func_0_0(12) * occ_func_0_1(11) *
               occ_func_0_0(29)) +
              (0.7071067812 * occ_func_0_0(7) * occ_func_0_1(2) *
               occ_func_0_0(22)) +
              ((0.7071067812 * occ_func_0_0(85) * occ_func_0_1(32) *
                    occ_func_0_0(39) +
                0.7071067812 * occ_func_0_0(85) * occ_func_0_0(32) *
                    occ_func_0_1(39))) +
              ((0.7071067812 * occ_func_0_1(4) * occ_func_0_0(10) *
                    occ_func_0_0(83) +
                0.7071067812 * occ_func_0_0(4) * occ_func_0_1(10) *
                    occ_func_0_0(83))) +
              (0.7071067812 * occ_func_0_0(9) * occ_func_0_1(12) *
               occ_func_0_0(35)) +
              (0.7071067812 * occ_func_0_0(3) * occ_func_0_1(1) *
               occ_func_0_0(28)) +
              ((0.7071067812 * occ_func_0_0(82) * occ_func_0_1(26) *
                    occ_func_0_0(33) +
                0.7071067812 * occ_func_0_0(82) * occ_func_0_0(26) *
                    occ_func_0_1(33))) +
              ((0.7071067812 * occ_func_0_1(4) * occ_func_0_0(1) *
                    occ_func_0_0(80) +
                0.7071067812 * occ_func_0_0(4) * occ_func_0_1(1) *
                    occ_func_0_0(80))) +
              (0.7071067812 * occ_func_0_0(9) * occ_func_0_1(3) *
               occ_func_0_0(24)) +
              (0.7071067812 * occ_func_0_0(12) * occ_func_0_1(10) *
               occ_func_0_0(29)) +
              ((0.7071067812 * occ_func_0_0(85) * occ_func_0_1(37) *
                    occ_func_0_0(32) +
                0.7071067812 * occ_func_0_0(85) * occ_func_0_0(37) *
                    occ_func_0_1(32)))) /
             24.0 +
         (m_occ_func_0_1[occ_f] - m_occ_func_0_1[occ_i]) *
             ((0.7071067812 * occ_func_0_0(4) * occ_func_0_0(5) *
               occ_func_0_0(37)) +
              (0.7071067812 * occ_func_0_0(6) * occ_func_0_0(8) *
               occ_func_0_0(39)) +
              (0.7071067812 * occ_func_0_0(10) * occ_func_0_0(12) *
               occ_func_0_0(33)) +
              (0.7071067812 * occ_func_0_0(4) * occ_func_0_0(1) *
               occ_func_0_0(26)) +
              (0.7071067812 * occ_func_0_0(11) * occ_func_0_0(10) *
               occ_func_0_0(30)) +
              (0.7071067812 * occ_func_0_0(8) * occ_func_0_0(3) *
               occ_func_0_0(25)) +
              (0.7071067812 * occ_func_0_0(1) * occ_func_0_0(4) *
               occ_func_0_0(40)) +
              (0.7071067812 * occ_func_0_0(3) * occ_func_0_0(9) *
               occ_func_0_0(42)) +
              (0.7071067812 * occ_func_0_0(5) * occ_func_0_0(2) *
               occ_func_0_0(27)) +
              (0.7071067812 * occ_func_0_0(10) * occ_func_0_0(11) *
               occ_func_0_0(33)) +
              (0.7071067812 * occ_func_0_0(2) * occ_func_0_0(7) *
               occ_func_0_0(41)) +
              (0.7071067812 * occ_func_0_0(1) * occ_func_0_0(6) *
               occ_func_0_0(40)) +
              (0.7071067812 * occ_func_0_0(6) * occ_func_0_0(1) *
               occ_func_0_0(23)) +
              (0.7071067812 * occ_func_0_0(11) * occ_func_0_0(12) *
               occ_func_0_0(30)) +
              (0.7071067812 * occ_func_0_0(7) * occ_func_0_0(9) *
               occ_func_0_0(38)) +
              (0.7071067812 * occ_func_0_0(5) * occ_func_0_0(4) *
               occ_func_0_0(36)) +
              (0.7071067812 * occ_func_0_0(5) * occ_func_0_0(10) *
               occ_func_0_0(36)) +
              (0.7071067812 * occ_func_0_0(2) * occ_func_0_0(3) *
               occ_func_0_0(31)) +
              (0.7071067812 * occ_func_0_0(1) * occ_func_0_0(3) *
               occ_func_0_0(32)) +
              (0.7071067812 * occ_func_0_0(4) * occ_func_0_0(10) *
               occ_func_0_0(37)) +
              (0.7071067812 * occ_func_0_0(8) * occ_func_0_0(6) *
               occ_func_0_0(34)) +
              (0.7071067812 * occ_func_0_0(9) * occ_func_0_0(7) *
               occ_func_0_0(35)) +
              (0.7071067812 * occ_func_0_0(6) * occ_func_0_0(11) *
               occ_func_0_0(39)) +
              (0.7071067812 * occ_func_0_0(1) * occ_func_0_0(2) *
               occ_func_0_0(32)) +
              (0.7071067812 * occ_func_0_0(12) * occ_func_0_0(7) *
               occ_func_0_0(21)) +
              (0.7071067812 * occ_func_0_0(11) * occ_func_0_0(6) *
               occ_func_0_0(20)) +
              (0.7071067812 * occ_func_0_0(10) * occ_func_0_0(4) *
               occ_func_0_0(19)) +
              (0.7071067812 * occ_func_0_0(12) * occ_func_0_0(9) *
               occ_func_0_0(21)) +
              (0.7071067812 * occ_func_0_0(3) * occ_func_0_0(2) *
               occ_func_0_0(28)) +
              (0.7071067812 * occ_func_0_0(8) * occ_func_0_0(11) *
               occ_func_0_0(34)) +
              (0.7071067812 * occ_func_0_0(3) * occ_func_0_0(8) *
               occ_func_0_0(42)) +
              (0.7071067812 * occ_func_0_0(2) * occ_func_0_0(5) *
               occ_func_0_0(41)) +
              (0.7071067812 * occ_func_0_0(7) * occ_func_0_0(5) *
               occ_func_0_0(22)) +
              (0.7071067812 * occ_func_0_0(9) * occ_func_0_0(8) *
               occ_func_0_0(24)) +
              (0.7071067812 * occ_func_0_0(4) * occ_func_0_0(6) *
               occ_func_0_0(26)) +
              (0.7071067812 * occ_func_0_0(5) * occ_func_0_0(7) *
               occ_func_0_0(27)) +
              (0.7071067812 * occ_func_0_0(11) * occ_func_0_0(8) *
               occ_func_0_0(20)) +
              (0.7071067812 * occ_func_0_0(10) * occ_func_0_0(5) *
               occ_func_0_0(19)) +
              (0.7071067812 * occ_func_0_0(2) * occ_func_0_0(1) *
               occ_func_0_0(31)) +
              (0.7071067812 * occ_func_0_0(7) * occ_func_0_0(12) *
               occ_func_0_0(38)) +
              (0.7071067812 * occ_func_0_0(8) * occ_func_0_0(9) *
               occ_func_0_0(25)) +
              (0.7071067812 * occ_func_0_0(6) * occ_func_0_0(4) *
               occ_func_0_0(23)) +
              (0.7071067812 * occ_func_0_0(12) * occ_func_0_0(11) *
               occ_func_0_0(29)) +
              (0.7071067812 * occ_func_0_0(7) * occ_func_0_0(2) *
               occ_func_0_0(22)) +
              (0.7071067812 * occ_func_0_0(9) * occ_func_0_0(12) *
               occ_func_0_0(35)) +
              (0.7071067812 * occ_func_0_0(3) * occ_func_0_0(1) *
               occ_func_0_0(28)) +
              (0.7071067812 * occ_func_0_0(9) * occ_func_0_0(3) *
               occ_func_0_0(24)) +
              (0.7071067812 * occ_func_0_0(12) * occ_func_0_0(10) *
               occ_func_0_0(29))) /
             24.0;
}
double test_Clexulator::delta_site_eval_at_0_bfunc_4_0_3(int occ_i,
                                                         int occ_f) const {
  return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) *
             ((0.7071067812 * occ_func_0_1(4) * occ_func_0_1(5) *
               occ_func_0_0(37)) +
              (0.7071067812 * occ_func_0_1(6) * occ_func_0_1(8) *
               occ_func_0_0(39)) +
              ((0.7071067812 * occ_func_0_1(80) * occ_func_0_1(24) *
                    occ_func_0_0(22) +
                0.7071067812 * occ_func_0_1(80) * occ_func_0_0(24) *
                    occ_func_0_1(22))) +
              (0.7071067812 * occ_func_0_1(10) * occ_func_0_1(12) *
               occ_func_0_0(33)) +
              (0.7071067812 * occ_func_0_1(4) * occ_func_0_1(1) *
               occ_func_0_0(26)) +
              ((0.7071067812 * occ_func_0_1(83) * occ_func_0_1(28) *
                    occ_func_0_0(35) +
                0.7071067812 * occ_func_0_1(83) * occ_func_0_0(28) *
                    occ_func_0_1(35))) +
              (0.7071067812 * occ_func_0_1(11) * occ_func_0_1(10) *
               occ_func_0_0(30)) +
              (0.7071067812 * occ_func_0_1(8) * occ_func_0_1(3) *
               occ_func_0_0(25)) +
              ((0.7071067812 * occ_func_0_1(84) * occ_func_0_1(31) *
                    occ_func_0_0(36) +
                0.7071067812 * occ_func_0_1(84) * occ_func_0_0(31) *
                    occ_func_0_1(36))) +
              (0.7071067812 * occ_func_0_1(1) * occ_func_0_1(4) *
               occ_func_0_0(40)) +
              (0.7071067812 * occ_func_0_1(3) * occ_func_0_1(9) *
               occ_func_0_0(42)) +
              ((0.7071067812 * occ_func_0_1(79) * occ_func_0_1(21) *
                    occ_func_0_0(19) +
                0.7071067812 * occ_func_0_1(79) * occ_func_0_0(21) *
                    occ_func_0_1(19))) +
              (0.7071067812 * occ_func_0_1(5) * occ_func_0_1(2) *
               occ_func_0_0(27)) +
              (0.7071067812 * occ_func_0_1(10) * occ_func_0_1(11) *
               occ_func_0_0(33)) +
              ((0.7071067812 * occ_func_0_1(83) * occ_func_0_1(34) *
                    occ_func_0_0(28) +
                0.7071067812 * occ_func_0_1(83) * occ_func_0_0(34) *
                    occ_func_0_1(28))) +
              (0.7071067812 * occ_func_0_1(2) * occ_func_0_1(7) *
               occ_func_0_0(41)) +
              (0.7071067812 * occ_func_0_1(1) * occ_func_0_1(6) *
               occ_func_0_0(40)) +
              ((0.7071067812 * occ_func_0_1(79) * occ_func_0_1(20) *
                    occ_func_0_0(21) +
                0.7071067812 * occ_func_0_1(79) * occ_func_0_0(20) *
                    occ_func_0_1(21))) +
              (0.7071067812 * occ_func_0_1(6) * occ_func_0_1(1) *
               occ_func_0_0(23)) +
              (0.7071067812 * occ_func_0_1(11) * occ_func_0_1(12) *
               occ_func_0_0(30)) +
              ((0.7071067812 * occ_func_0_1(84) * occ_func_0_1(38) *
                    occ_func_0_0(31) +
                0.7071067812 * occ_func_0_1(84) * occ_func_0_0(38) *
                    occ_func_0_1(31))) +
              (0.7071067812 * occ_func_0_1(7) * occ_func_0_1(9) *
               occ_func_0_0(38)) +
              (0.7071067812 * occ_func_0_1(5) * occ_func_0_1(4) *
               occ_func_0_0(36)) +
              ((0.7071067812 * occ_func_0_1(81) * occ_func_0_1(23) *
                    occ_func_0_0(25) +
                0.7071067812 * occ_func_0_1(81) * occ_func_0_0(23) *
                    occ_func_0_1(25))) +
              (0.7071067812 * occ_func_0_1(5) * occ_func_0_1(10) *
               occ_func_0_0(36)) +
              (0.7071067812 * occ_func_0_1(2) * occ_func_0_1(3) *
               occ_func_0_0(31)) +
              ((0.7071067812 * occ_func_0_1(81) * occ_func_0_1(25) *
                    occ_func_0_0(30) +
                0.7071067812 * occ_func_0_1(81) * occ_func_0_0(25) *
                    occ_func_0_1(30))) +
              (0.7071067812 * occ_func_0_1(1) * occ_func_0_1(3) *
               occ_func_0_0(32)) +
              (0.7071067812 * occ_func_0_1(4) * occ_func_0_1(10) *
               occ_func_0_0(37)) +
              ((0.7071067812 * occ_func_0_1(80) * occ_func_0_1(29) *
                    occ_func_0_0(24) +
                0.7071067812 * occ_func_0_1(80) * occ_func_0_0(29) *
                    occ_func_0_1(24))) +
              (0.7071067812 * occ_func_0_1(8) * occ_func_0_1(6) *
               occ_func_0_0(34)) +
              (0.7071067812 * occ_func_0_1(9) * occ_func_0_1(7) *
               occ_func_0_0(35)) +
              ((0.7071067812 * occ_func_0_1(82) * occ_func_0_1(27) *
                    occ_func_0_0(26) +
                0.7071067812 * occ_func_0_1(82) * occ_func_0_0(27) *
                    occ_func_0_1(26))) +
              (0.7071067812 * occ_func_0_1(6) * occ_func_0_1(11) *
               occ_func_0_0(39)) +
              (0.7071067812 * occ_func_0_1(1) * occ_func_0_1(2) *
               occ_func_0_0(32)) +
              ((0.7071067812 * occ_func_0_1(80) * occ_func_0_1(22) *
                    occ_func_0_0(29) +
                0.7071067812 * occ_func_0_1(80) * occ_func_0_0(22) *
                    occ_func_0_1(29))) +
              (0.7071067812 * occ_func_0_1(12) * occ_func_0_1(7) *
               occ_func_0_0(21)) +
              (0.7071067812 * occ_func_0_1(11) * occ_func_0_1(6) *
               occ_func_0_0(20)) +
              ((0.7071067812 * occ_func_0_1(86) * occ_func_0_1(40) *
                    occ_func_0_0(41) +
                0.7071067812 * occ_func_0_1(86) * occ_func_0_0(40) *
                    occ_func_0_1(41))) +
              (0.7071067812 * occ_func_0_1(10) * occ_func_0_1(4) *
               occ_func_0_0(19)) +
              (0.7071067812 * occ_func_0_1(12) * occ_func_0_1(9) *
               occ_func_0_0(21)) +
              ((0.7071067812 * occ_func_0_1(86) * occ_func_0_1(42) *
                    occ_func_0_0(40) +
                0.7071067812 * occ_func_0_1(86) * occ_func_0_0(42) *
                    occ_func_0_1(40))) +
              (0.7071067812 * occ_func_0_1(3) * occ_func_0_1(2) *
               occ_func_0_0(28)) +
              (0.7071067812 * occ_func_0_1(8) * occ_func_0_1(11) *
               occ_func_0_0(34)) +
              ((0.7071067812 * occ_func_0_1(82) * occ_func_0_1(33) *
                    occ_func_0_0(27) +
                0.7071067812 * occ_func_0_1(82) * occ_func_0_0(33) *
                    occ_func_0_1(27))) +
              (0.7071067812 * occ_func_0_1(3) * occ_func_0_1(8) *
               occ_func_0_0(42)) +
              (0.7071067812 * occ_func_0_1(2) * occ_func_0_1(5) *
               occ_func_0_0(41)) +
              ((0.7071067812 * occ_func_0_1(79) * occ_func_0_1(19) *
                    occ_func_0_0(20) +
                0.7071067812 * occ_func_0_1(79) * occ_func_0_0(19) *
                    occ_func_0_1(20))) +
              (0.7071067812 * occ_func_0_1(7) * occ_func_0_1(5) *
               occ_func_0_0(22)) +
              (0.7071067812 * occ_func_0_1(9) * occ_func_0_1(8) *
               occ_func_0_0(24)) +
              ((0.7071067812 * occ_func_0_1(85) * occ_func_0_1(39) *
                    occ_func_0_0(37) +
                0.7071067812 * occ_func_0_1(85) * occ_func_0_0(39) *
                    occ_func_0_1(37))) +
              (0.7071067812 * occ_func_0_1(4) * occ_func_0_1(6) *
               occ_func_0_0(26)) +
              (0.7071067812 * occ_func_0_1(5) * occ_func_0_1(7) *
               occ_func_0_0(27)) +
              ((0.7071067812 * occ_func_0_1(83) * occ_func_0_1(35) *
                    occ_func_0_0(34) +
                0.7071067812 * occ_func_0_1(83) * occ_func_0_0(35) *
                    occ_func_0_1(34))) +
              (0.7071067812 * occ_func_0_1(11) * occ_func_0_1(8) *
               occ_func_0_0(20)) +
              (0.7071067812 * occ_func_0_1(10) * occ_func_0_1(5) *
               occ_func_0_0(19)) +
              ((0.7071067812 * occ_func_0_1(86) * occ_func_0_1(41) *
                    occ_func_0_0(42) +
                0.7071067812 * occ_func_0_1(86) * occ_func_0_0(41) *
                    occ_func_0_1(42))) +
              (0.7071067812 * occ_func_0_1(2) * occ_func_0_1(1) *
               occ_func_0_0(31)) +
              (0.7071067812 * occ_func_0_1(7) * occ_func_0_1(12) *
               occ_func_0_0(38)) +
              ((0.7071067812 * occ_func_0_1(81) * occ_func_0_1(30) *
                    occ_func_0_0(23) +
                0.7071067812 * occ_func_0_1(81) * occ_func_0_0(30) *
                    occ_func_0_1(23))) +
              (0.7071067812 * occ_func_0_1(8) * occ_func_0_1(9) *
               occ_func_0_0(25)) +
              (0.7071067812 * occ_func_0_1(6) * occ_func_0_1(4) *
               occ_func_0_0(23)) +
              ((0.7071067812 * occ_func_0_1(84) * occ_func_0_1(36) *
                    occ_func_0_0(38) +
                0.7071067812 * occ_func_0_1(84) * occ_func_0_0(36) *
                    occ_func_0_1(38))) +
              (0.7071067812 * occ_func_0_1(12) * occ_func_0_1(11) *
               occ_func_0_0(29)) +
              (0.7071067812 * occ_func_0_1(7) * occ_func_0_1(2) *
               occ_func_0_0(22)) +
              ((0.7071067812 * occ_func_0_1(85) * occ_func_0_1(32) *
                    occ_func_0_0(39) +
                0.7071067812 * occ_func_0_1(85) * occ_func_0_0(32) *
                    occ_func_0_1(39))) +
              (0.7071067812 * occ_func_0_1(9) * occ_func_0_1(12) *
               occ_func_0_0(35)) +
              (0.7071067812 * occ_func_0_1(3) * occ_func_0_1(1) *
               occ_func_0_0(28)) +
              ((0.7071067812 * occ_func_0_1(82) * occ_func_0_1(26) *
                    occ_func_0_0(33) +
                0.7071067812 * occ_func_0_1(82) * occ_func_0_0(26) *
                    occ_func_0_1(33))) +
              (0.7071067812 * occ_func_0_1(9) * occ_func_0_1(3) *
               occ_func_0_0(24)) +
              (0.7071067812 * occ_func_0_1(12) * occ_func_0_1(10) *
               occ_func_0_0(29)) +
              ((0.7071067812 * occ_func_0_1(85) * occ_func_0_1(37) *
                    occ_func_0_0(32) +
                0.7071067812 * occ_func_0_1(85) * occ_func_0_0(37) *
                    occ_func_0_1(32)))) /
             24.0 +
         (m_occ_func_0_1[occ_f] - m_occ_func_0_1[occ_i]) *
             (((0.7071067812 * occ_func_0_1(9) * occ_func_0_0(7) *
                    occ_func_0_0(85) +
                0.7071067812 * occ_func_0_0(9) * occ_func_0_1(7) *
                    occ_func_0_0(85))) +
              (0.7071067812 * occ_func_0_1(4) * occ_func_0_0(5) *
               occ_func_0_0(37)) +
              (0.7071067812 * occ_func_0_1(6) * occ_func_0_0(8) *
               occ_func_0_0(39)) +
              ((0.7071067812 * occ_func_0_1(3) * occ_func_0_0(9) *
                    occ_func_0_0(82) +
                0.7071067812 * occ_func_0_0(3) * occ_func_0_1(9) *
                    occ_func_0_0(82))) +
              (0.7071067812 * occ_func_0_1(10) * occ_func_0_0(12) *
               occ_func_0_0(33)) +
              (0.7071067812 * occ_func_0_1(4) * occ_func_0_0(1) *
               occ_func_0_0(26)) +
              ((0.7071067812 * occ_func_0_1(2) * occ_func_0_0(5) *
                    occ_func_0_0(81) +
                0.7071067812 * occ_func_0_0(2) * occ_func_0_1(5) *
                    occ_func_0_0(81))) +
              (0.7071067812 * occ_func_0_1(11) * occ_func_0_0(10) *
               occ_func_0_0(30)) +
              (0.7071067812 * occ_func_0_1(8) * occ_func_0_0(3) *
               occ_func_0_0(25)) +
              ((0.7071067812 * occ_func_0_1(12) * occ_func_0_0(10) *
                    occ_func_0_0(86) +
                0.7071067812 * occ_func_0_0(12) * occ_func_0_1(10) *
                    occ_func_0_0(86))) +
              (0.7071067812 * occ_func_0_1(1) * occ_func_0_0(4) *
               occ_func_0_0(40)) +
              (0.7071067812 * occ_func_0_1(3) * occ_func_0_0(9) *
               occ_func_0_0(42)) +
              ((0.7071067812 * occ_func_0_1(8) * occ_func_0_0(3) *
                    occ_func_0_0(82) +
                0.7071067812 * occ_func_0_0(8) * occ_func_0_1(3) *
                    occ_func_0_0(82))) +
              (0.7071067812 * occ_func_0_1(5) * occ_func_0_0(2) *
               occ_func_0_0(27)) +
              (0.7071067812 * occ_func_0_1(10) * occ_func_0_0(11) *
               occ_func_0_0(33)) +
              ((0.7071067812 * occ_func_0_1(11) * occ_func_0_0(12) *
                    occ_func_0_0(86) +
                0.7071067812 * occ_func_0_0(11) * occ_func_0_1(12) *
                    occ_func_0_0(86))) +
              (0.7071067812 * occ_func_0_1(2) * occ_func_0_0(7) *
               occ_func_0_0(41)) +
              (0.7071067812 * occ_func_0_1(1) * occ_func_0_0(6) *
               occ_func_0_0(40)) +
              ((0.7071067812 * occ_func_0_1(7) * occ_func_0_0(2) *
                    occ_func_0_0(81) +
                0.7071067812 * occ_func_0_0(7) * occ_func_0_1(2) *
                    occ_func_0_0(81))) +
              (0.7071067812 * occ_func_0_1(6) * occ_func_0_0(1) *
               occ_func_0_0(23)) +
              (0.7071067812 * occ_func_0_1(11) * occ_func_0_0(12) *
               occ_func_0_0(30)) +
              ((0.7071067812 * occ_func_0_1(6) * occ_func_0_0(8) *
                    occ_func_0_0(84) +
                0.7071067812 * occ_func_0_0(6) * occ_func_0_1(8) *
                    occ_func_0_0(84))) +
              (0.7071067812 * occ_func_0_1(7) * occ_func_0_0(9) *
               occ_func_0_0(38)) +
              (0.7071067812 * occ_func_0_1(5) * occ_func_0_0(4) *
               occ_func_0_0(36)) +
              ((0.7071067812 * occ_func_0_1(8) * occ_func_0_0(11) *
                    occ_func_0_0(84) +
                0.7071067812 * occ_func_0_0(8) * occ_func_0_1(11) *
                    occ_func_0_0(84))) +
              (0.7071067812 * occ_func_0_1(5) * occ_func_0_0(10) *
               occ_func_0_0(36)) +
              (0.7071067812 * occ_func_0_1(2) * occ_func_0_0(3) *
               occ_func_0_0(31)) +
              ((0.7071067812 * occ_func_0_1(12) * occ_func_0_0(9) *
                    occ_func_0_0(85) +
                0.7071067812 * occ_func_0_0(12) * occ_func_0_1(9) *
                    occ_func_0_0(85))) +
              (0.7071067812 * occ_func_0_1(1) * occ_func_0_0(3) *
               occ_func_0_0(32)) +
              (0.7071067812 * occ_func_0_1(4) * occ_func_0_0(10) *
               occ_func_0_0(37)) +
              ((0.7071067812 * occ_func_0_1(5) * occ_func_0_0(4) *
                    occ_func_0_0(83) +
                0.7071067812 * occ_func_0_0(5) * occ_func_0_1(4) *
                    occ_func_0_0(83))) +
              (0.7071067812 * occ_func_0_1(8) * occ_func_0_0(6) *
               occ_func_0_0(34)) +
              (0.7071067812 * occ_func_0_1(9) * occ_func_0_0(7) *
               occ_func_0_0(35)) +
              ((0.7071067812 * occ_func_0_1(7) * occ_func_0_0(12) *
                    occ_func_0_0(85) +
                0.7071067812 * occ_func_0_0(7) * occ_func_0_1(12) *
                    occ_func_0_0(85))) +
              (0.7071067812 * occ_func_0_1(6) * occ_func_0_0(11) *
               occ_func_0_0(39)) +
              (0.7071067812 * occ_func_0_1(1) * occ_func_0_0(2) *
               occ_func_0_0(32)) +
              ((0.7071067812 * occ_func_0_1(1) * occ_func_0_0(2) *
                    occ_func_0_0(79) +
                0.7071067812 * occ_func_0_0(1) * occ_func_0_1(2) *
                    occ_func_0_0(79))) +
              (0.7071067812 * occ_func_0_1(12) * occ_func_0_0(7) *
               occ_func_0_0(21)) +
              (0.7071067812 * occ_func_0_1(11) * occ_func_0_0(6) *
               occ_func_0_0(20)) +
              ((0.7071067812 * occ_func_0_1(3) * occ_func_0_0(1) *
                    occ_func_0_0(79) +
                0.7071067812 * occ_func_0_0(3) * occ_func_0_1(1) *
                    occ_func_0_0(79))) +
              (0.7071067812 * occ_func_0_1(10) * occ_func_0_0(4) *
               occ_func_0_0(19)) +
              (0.7071067812 * occ_func_0_1(12) * occ_func_0_0(9) *
               occ_func_0_0(21)) +
              ((0.7071067812 * occ_func_0_1(10) * occ_func_0_0(5) *
                    occ_func_0_0(83) +
                0.7071067812 * occ_func_0_0(10) * occ_func_0_1(5) *
                    occ_func_0_0(83))) +
              (0.7071067812 * occ_func_0_1(3) * occ_func_0_0(2) *
               occ_func_0_0(28)) +
              (0.7071067812 * occ_func_0_1(8) * occ_func_0_0(11) *
               occ_func_0_0(34)) +
              ((0.7071067812 * occ_func_0_1(10) * occ_func_0_0(11) *
                    occ_func_0_0(86) +
                0.7071067812 * occ_func_0_0(10) * occ_func_0_1(11) *
                    occ_func_0_0(86))) +
              (0.7071067812 * occ_func_0_1(3) * occ_func_0_0(8) *
               occ_func_0_0(42)) +
              (0.7071067812 * occ_func_0_1(2) * occ_func_0_0(5) *
               occ_func_0_0(41)) +
              ((0.7071067812 * occ_func_0_1(6) * occ_func_0_0(4) *
                    occ_func_0_0(80) +
                0.7071067812 * occ_func_0_0(6) * occ_func_0_1(4) *
                    occ_func_0_0(80))) +
              (0.7071067812 * occ_func_0_1(7) * occ_func_0_0(5) *
               occ_func_0_0(22)) +
              (0.7071067812 * occ_func_0_1(9) * occ_func_0_0(8) *
               occ_func_0_0(24)) +
              ((0.7071067812 * occ_func_0_1(9) * occ_func_0_0(8) *
                    occ_func_0_0(82) +
                0.7071067812 * occ_func_0_0(9) * occ_func_0_1(8) *
                    occ_func_0_0(82))) +
              (0.7071067812 * occ_func_0_1(4) * occ_func_0_0(6) *
               occ_func_0_0(26)) +
              (0.7071067812 * occ_func_0_1(5) * occ_func_0_0(7) *
               occ_func_0_0(27)) +
              ((0.7071067812 * occ_func_0_1(2) * occ_func_0_0(3) *
                    occ_func_0_0(79) +
                0.7071067812 * occ_func_0_0(2) * occ_func_0_1(3) *
                    occ_func_0_0(79))) +
              (0.7071067812 * occ_func_0_1(11) * occ_func_0_0(8) *
               occ_func_0_0(20)) +
              (0.7071067812 * occ_func_0_1(10) * occ_func_0_0(5) *
               occ_func_0_0(19)) +
              ((0.7071067812 * occ_func_0_1(11) * occ_func_0_0(6) *
                    occ_func_0_0(84) +
                0.7071067812 * occ_func_0_0(11) * occ_func_0_1(6) *
                    occ_func_0_0(84))) +
              (0.7071067812 * occ_func_0_1(2) * occ_func_0_0(1) *
               occ_func_0_0(31)) +
              (0.7071067812 * occ_func_0_1(7) * occ_func_0_0(12) *
               occ_func_0_0(38)) +
              ((0.7071067812 * occ_func_0_1(5) * occ_func_0_0(7) *
                    occ_func_0_0(81) +
                0.7071067812 * occ_func_0_0(5) * occ_func_0_1(7) *
                    occ_func_0_0(81))) +
              (0.7071067812 * occ_func_0_1(8) * occ_func_0_0(9) *
               occ_func_0_0(25)) +
              (0.7071067812 * occ_func_0_1(6) * occ_func_0_0(4) *
               occ_func_0_0(23)) +
              ((0.7071067812 * occ_func_0_1(1) * occ_func_0_0(6) *
                    occ_func_0_0(80) +
                0.7071067812 * occ_func_0_0(1) * occ_func_0_1(6) *
                    occ_func_0_0(80))) +
              (0.7071067812 * occ_func_0_1(12) * occ_func_0_0(11) *
               occ_func_0_0(29)) +
              (0.7071067812 * occ_func_0_1(7) * occ_func_0_0(2) *
               occ_func_0_0(22)) +
              ((0.7071067812 * occ_func_0_1(4) * occ_func_0_0(10) *
                    occ_func_0_0(83) +
                0.7071067812 * occ_func_0_0(4) * occ_func_0_1(10) *
                    occ_func_0_0(83))) +
              (0.7071067812 * occ_func_0_1(9) * occ_func_0_0(12) *
               occ_func_0_0(35)) +
              (0.7071067812 * occ_func_0_1(3) * occ_func_0_0(1) *
               occ_func_0_0(28)) +
              ((0.7071067812 * occ_func_0_1(4) * occ_func_0_0(1) *
                    occ_func_0_0(80) +
                0.7071067812 * occ_func_0_0(4) * occ_func_0_1(1) *
                    occ_func_0_0(80))) +
              (0.7071067812 * occ_func_0_1(9) * occ_func_0_0(3) *
               occ_func_0_0(24)) +
              (0.7071067812 * occ_func_0_1(12) * occ_func_0_0(10) *
               occ_func_0_0(29))) /
             24.0;
}
double test_Clexulator::delta_site_eval_at_0_bfunc_4_0_4(int occ_i,
                                                         int occ_f) const {
  return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) *
             ((occ_func_0_1(9) * occ_func_0_1(7) * occ_func_0_0(85)) +
              (occ_func_0_0(80) * occ_func_0_1(24) * occ_func_0_1(22)) +
              (occ_func_0_1(3) * occ_func_0_1(9) * occ_func_0_0(82)) +
              (occ_func_0_0(83) * occ_func_0_1(28) * occ_func_0_1(35)) +
              (occ_func_0_1(2) * occ_func_0_1(5) * occ_func_0_0(81)) +
              (occ_func_0_0(84) * occ_func_0_1(31) * occ_func_0_1(36)) +
              (occ_func_0_1(12) * occ_func_0_1(10) * occ_func_0_0(86)) +
              (occ_func_0_0(79) * occ_func_0_1(21) * occ_func_0_1(19)) +
              (occ_func_0_1(8) * occ_func_0_1(3) * occ_func_0_0(82)) +
              (occ_func_0_0(83) * occ_func_0_1(34) * occ_func_0_1(28)) +
              (occ_func_0_1(11) * occ_func_0_1(12) * occ_func_0_0(86)) +
              (occ_func_0_0(79) * occ_func_0_1(20) * occ_func_0_1(21)) +
              (occ_func_0_1(7) * occ_func_0_1(2) * occ_func_0_0(81)) +
              (occ_func_0_0(84) * occ_func_0_1(38) * occ_func_0_1(31)) +
              (occ_func_0_1(6) * occ_func_0_1(8) * occ_func_0_0(84)) +
              (occ_func_0_0(81) * occ_func_0_1(23) * occ_func_0_1(25)) +
              (occ_func_0_1(8) * occ_func_0_1(11) * occ_func_0_0(84)) +
              (occ_func_0_0(81) * occ_func_0_1(25) * occ_func_0_1(30)) +
              (occ_func_0_1(12) * occ_func_0_1(9) * occ_func_0_0(85)) +
              (occ_func_0_0(80) * occ_func_0_1(29) * occ_func_0_1(24)) +
              (occ_func_0_1(5) * occ_func_0_1(4) * occ_func_0_0(83)) +
              (occ_func_0_0(82) * occ_func_0_1(27) * occ_func_0_1(26)) +
              (occ_func_0_1(7) * occ_func_0_1(12) * occ_func_0_0(85)) +
              (occ_func_0_0(80) * occ_func_0_1(22) * occ_func_0_1(29)) +
              (occ_func_0_1(1) * occ_func_0_1(2) * occ_func_0_0(79)) +
              (occ_func_0_0(86) * occ_func_0_1(40) * occ_func_0_1(41)) +
              (occ_func_0_1(3) * occ_func_0_1(1) * occ_func_0_0(79)) +
              (occ_func_0_0(86) * occ_func_0_1(42) * occ_func_0_1(40)) +
              (occ_func_0_1(10) * occ_func_0_1(5) * occ_func_0_0(83)) +
              (occ_func_0_0(82) * occ_func_0_1(33) * occ_func_0_1(27)) +
              (occ_func_0_1(10) * occ_func_0_1(11) * occ_func_0_0(86)) +
              (occ_func_0_0(79) * occ_func_0_1(19) * occ_func_0_1(20)) +
              (occ_func_0_1(6) * occ_func_0_1(4) * occ_func_0_0(80)) +
              (occ_func_0_0(85) * occ_func_0_1(39) * occ_func_0_1(37)) +
              (occ_func_0_1(9) * occ_func_0_1(8) * occ_func_0_0(82)) +
              (occ_func_0_0(83) * occ_func_0_1(35) * occ_func_0_1(34)) +
              (occ_func_0_1(2) * occ_func_0_1(3) * occ_func_0_0(79)) +
              (occ_func_0_0(86) * occ_func_0_1(41) * occ_func_0_1(42)) +
              (occ_func_0_1(11) * occ_func_0_1(6) * occ_func_0_0(84)) +
              (occ_func_0_0(81) * occ_func_0_1(30) * occ_func_0_1(23)) +
              (occ_func_0_1(5) * occ_func_0_1(7) * occ_func_0_0(81)) +
              (occ_func_0_0(84) * occ_func_0_1(36) * occ_func_0_1(38)) +
              (occ_func_0_1(1) * occ_func_0_1(6) * occ_func_0_0(80)) +
              (occ_func_0_0(85) * occ_func_0_1(32) * occ_func_0_1(39)) +
              (occ_func_0_1(4) * occ_func_0_1(10) * occ_func_0_0(83)) +
              (occ_func_0_0(82) * occ_func_0_1(26) * occ_func_0_1(33)) +
              (occ_func_0_1(4) * occ_func_0_1(1) * occ_func_0_0(80)) +
              (occ_func_0_0(85) * occ_func_0_1(37) * occ_func_0_1(32))) /
             24.0 +
         (m_occ_func_0_1[occ_f] - m_occ_func_0_1[occ_i]) *
             ((occ_func_0_0(4) * occ_func_0_1(5) * occ_func_0_0(37)) +
              (occ_func_0_0(6) * occ_func_0_1(8) * occ_func_0_0(39)) +
              (occ_func_0_0(10) * occ_func_0_1(12) * occ_func_0_0(33)) +
              (occ_func_0_0(4) * occ_func_0_1(1) * occ_func_0_0(26)) +
              (occ_func_0_0(11) * occ_func_0_1(10) * occ_func_0_0(30)) +
              (occ_func_0_0(8) * occ_func_0_1(3) * occ_func_0_0(25)) +
              (occ_func_0_0(1) * occ_func_0_1(4) * occ_func_0_0(40)) +
              (occ_func_0_0(3) * occ_func_0_1(9) * occ_func_0_0(42)) +
              (occ_func_0_0(5) * occ_func_0_1(2) * occ_func_0_0(27)) +
              (occ_func_0_0(10) * occ_func_0_1(11) * occ_func_0_0(33)) +
              (occ_func_0_0(2) * occ_func_0_1(7) * occ_func_0_0(41)) +
              (occ_func_0_0(1) * occ_func_0_1(6) * occ_func_0_0(40)) +
              (occ_func_0_0(6) * occ_func_0_1(1) * occ_func_0_0(23)) +
              (occ_func_0_0(11) * occ_func_0_1(12) * occ_func_0_0(30)) +
              (occ_func_0_0(7) * occ_func_0_1(9) * occ_func_0_0(38)) +
              (occ_func_0_0(5) * occ_func_0_1(4) * occ_func_0_0(36)) +
              (occ_func_0_0(5) * occ_func_0_1(10) * occ_func_0_0(36)) +
              (occ_func_0_0(2) * occ_func_0_1(3) * occ_func_0_0(31)) +
              (occ_func_0_0(1) * occ_func_0_1(3) * occ_func_0_0(32)) +
              (occ_func_0_0(4) * occ_func_0_1(10) * occ_func_0_0(37)) +
              (occ_func_0_0(8) * occ_func_0_1(6) * occ_func_0_0(34)) +
              (occ_func_0_0(9) * occ_func_0_1(7) * occ_func_0_0(35)) +
              (occ_func_0_0(6) * occ_func_0_1(11) * occ_func_0_0(39)) +
              (occ_func_0_0(1) * occ_func_0_1(2) * occ_func_0_0(32)) +
              (occ_func_0_0(12) * occ_func_0_1(7) * occ_func_0_0(21)) +
              (occ_func_0_0(11) * occ_func_0_1(6) * occ_func_0_0(20)) +
              (occ_func_0_0(10) * occ_func_0_1(4) * occ_func_0_0(19)) +
              (occ_func_0_0(12) * occ_func_0_1(9) * occ_func_0_0(21)) +
              (occ_func_0_0(3) * occ_func_0_1(2) * occ_func_0_0(28)) +
              (occ_func_0_0(8) * occ_func_0_1(11) * occ_func_0_0(34)) +
              (occ_func_0_0(3) * occ_func_0_1(8) * occ_func_0_0(42)) +
              (occ_func_0_0(2) * occ_func_0_1(5) * occ_func_0_0(41)) +
              (occ_func_0_0(7) * occ_func_0_1(5) * occ_func_0_0(22)) +
              (occ_func_0_0(9) * occ_func_0_1(8) * occ_func_0_0(24)) +
              (occ_func_0_0(4) * occ_func_0_1(6) * occ_func_0_0(26)) +
              (occ_func_0_0(5) * occ_func_0_1(7) * occ_func_0_0(27)) +
              (occ_func_0_0(11) * occ_func_0_1(8) * occ_func_0_0(20)) +
              (occ_func_0_0(10) * occ_func_0_1(5) * occ_func_0_0(19)) +
              (occ_func_0_0(2) * occ_func_0_1(1) * occ_func_0_0(31)) +
              (occ_func_0_0(7) * occ_func_0_1(12) * occ_func_0_0(38)) +
              (occ_func_0_0(8) * occ_func_0_1(9) * occ_func_0_0(25)) +
              (occ_func_0_0(6) * occ_func_0_1(4) * occ_func_0_0(23)) +
              (occ_func_0_0(12) * occ_func_0_1(11) * occ_func_0_0(29)) +
              (occ_func_0_0(7) * occ_func_0_1(2) * occ_func_0_0(22)) +
              (occ_func_0_0(9) * occ_func_0_1(12) * occ_func_0_0(35)) +
              (occ_func_0_0(3) * occ_func_0_1(1) * occ_func_0_0(28)) +
              (occ_func_0_0(9) * occ_func_0_1(3) * occ_func_0_0(24)) +
              (occ_func_0_0(12) * occ_func_0_1(10) * occ_func_0_0(29))) /
             24.0;
}
double test_Clexulator::delta_site_eval_at_0_bfunc_4_0_5(int occ_i,
                                                         int occ_f) const {
  return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) *
             ((occ_func_0_1(80) * occ_func_0_1(24) * occ_func_0_1(22)) +
              (occ_func_0_1(83) * occ_func_0_1(28) * occ_func_0_1(35)) +
              (occ_func_0_1(84) * occ_func_0_1(31) * occ_func_0_1(36)) +
              (occ_func_0_1(79) * occ_func_0_1(21) * occ_func_0_1(19)) +
              (occ_func_0_1(83) * occ_func_0_1(34) * occ_func_0_1(28)) +
              (occ_func_0_1(79) * occ_func_0_1(20) * occ_func_0_1(21)) +
              (occ_func_0_1(84) * occ_func_0_1(38) * occ_func_0_1(31)) +
              (occ_func_0_1(81) * occ_func_0_1(23) * occ_func_0_1(25)) +
              (occ_func_0_1(81) * occ_func_0_1(25) * occ_func_0_1(30)) +
              (occ_func_0_1(80) * occ_func_0_1(29) * occ_func_0_1(24)) +
              (occ_func_0_1(82) * occ_func_0_1(27) * occ_func_0_1(26)) +
              (occ_func_0_1(80) * occ_func_0_1(22) * occ_func_0_1(29)) +
              (occ_func_0_1(86) * occ_func_0_1(40) * occ_func_0_1(41)) +
              (occ_func_0_1(86) * occ_func_0_1(42) * occ_func_0_1(40)) +
              (occ_func_0_1(82) * occ_func_0_1(33) * occ_func_0_1(27)) +
              (occ_func_0_1(79) * occ_func_0_1(19) * occ_func_0_1(20)) +
              (occ_func_0_1(85) * occ_func_0_1(39) * occ_func_0_1(37)) +
              (occ_func_0_1(83) * occ_func_0_1(35) * occ_func_0_1(34)) +
              (occ_func_0_1(86) * occ_func_0_1(41) * occ_func_0_1(42)) +
              (occ_func_0_1(81) * occ_func_0_1(30) * occ_func_0_1(23)) +
              (occ_func_0_1(84) * occ_func_0_1(36) * occ_func_0_1(38)) +
              (occ_func_0_1(85) * occ_func_0_1(32) * occ_func_0_1(39)) +
              (occ_func_0_1(82) * occ_func_0_1(26) * occ_func_0_1(33)) +
              (occ_func_0_1(85) * occ_func_0_1(37) * occ_func_0_1(32))) /
             24.0 +
         (m_occ_func_0_1[occ_f] - m_occ_func_0_1[occ_i]) *
             ((occ_func_0_1(9) * occ_func_0_1(7) * occ_func_0_0(85)) +
              (occ_func_0_1(4) * occ_func_0_1(5) * occ_func_0_0(37)) +
              (occ_func_0_1(6) * occ_func_0_1(8) * occ_func_0_0(39)) +
              (occ_func_0_1(3) * occ_func_0_1(9) * occ_func_0_0(82)) +
              (occ_func_0_1(10) * occ_func_0_1(12) * occ_func_0_0(33)) +
              (occ_func_0_1(4) * occ_func_0_1(1) * occ_func_0_0(26)) +
              (occ_func_0_1(2) * occ_func_0_1(5) * occ_func_0_0(81)) +
              (occ_func_0_1(11) * occ_func_0_1(10) * occ_func_0_0(30)) +
              (occ_func_0_1(8) * occ_func_0_1(3) * occ_func_0_0(25)) +
              (occ_func_0_1(12) * occ_func_0_1(10) * occ_func_0_0(86)) +
              (occ_func_0_1(1) * occ_func_0_1(4) * occ_func_0_0(40)) +
              (occ_func_0_1(3) * occ_func_0_1(9) * occ_func_0_0(42)) +
              (occ_func_0_1(8) * occ_func_0_1(3) * occ_func_0_0(82)) +
              (occ_func_0_1(5) * occ_func_0_1(2) * occ_func_0_0(27)) +
              (occ_func_0_1(10) * occ_func_0_1(11) * occ_func_0_0(33)) +
              (occ_func_0_1(11) * occ_func_0_1(12) * occ_func_0_0(86)) +
              (occ_func_0_1(2) * occ_func_0_1(7) * occ_func_0_0(41)) +
              (occ_func_0_1(1) * occ_func_0_1(6) * occ_func_0_0(40)) +
              (occ_func_0_1(7) * occ_func_0_1(2) * occ_func_0_0(81)) +
              (occ_func_0_1(6) * occ_func_0_1(1) * occ_func_0_0(23)) +
              (occ_func_0_1(11) * occ_func_0_1(12) * occ_func_0_0(30)) +
              (occ_func_0_1(6) * occ_func_0_1(8) * occ_func_0_0(84)) +
              (occ_func_0_1(7) * occ_func_0_1(9) * occ_func_0_0(38)) +
              (occ_func_0_1(5) * occ_func_0_1(4) * occ_func_0_0(36)) +
              (occ_func_0_1(8) * occ_func_0_1(11) * occ_func_0_0(84)) +
              (occ_func_0_1(5) * occ_func_0_1(10) * occ_func_0_0(36)) +
              (occ_func_0_1(2) * occ_func_0_1(3) * occ_func_0_0(31)) +
              (occ_func_0_1(12) * occ_func_0_1(9) * occ_func_0_0(85)) +
              (occ_func_0_1(1) * occ_func_0_1(3) * occ_func_0_0(32)) +
              (occ_func_0_1(4) * occ_func_0_1(10) * occ_func_0_0(37)) +
              (occ_func_0_1(5) * occ_func_0_1(4) * occ_func_0_0(83)) +
              (occ_func_0_1(8) * occ_func_0_1(6) * occ_func_0_0(34)) +
              (occ_func_0_1(9) * occ_func_0_1(7) * occ_func_0_0(35)) +
              (occ_func_0_1(7) * occ_func_0_1(12) * occ_func_0_0(85)) +
              (occ_func_0_1(6) * occ_func_0_1(11) * occ_func_0_0(39)) +
              (occ_func_0_1(1) * occ_func_0_1(2) * occ_func_0_0(32)) +
              (occ_func_0_1(1) * occ_func_0_1(2) * occ_func_0_0(79)) +
              (occ_func_0_1(12) * occ_func_0_1(7) * occ_func_0_0(21)) +
              (occ_func_0_1(11) * occ_func_0_1(6) * occ_func_0_0(20)) +
              (occ_func_0_1(3) * occ_func_0_1(1) * occ_func_0_0(79)) +
              (occ_func_0_1(10) * occ_func_0_1(4) * occ_func_0_0(19)) +
              (occ_func_0_1(12) * occ_func_0_1(9) * occ_func_0_0(21)) +
              (occ_func_0_1(10) * occ_func_0_1(5) * occ_func_0_0(83)) +
              (occ_func_0_1(3) * occ_func_0_1(2) * occ_func_0_0(28)) +
              (occ_func_0_1(8) * occ_func_0_1(11) * occ_func_0_0(34)) +
              (occ_func_0_1(10) * occ_func_0_1(11) * occ_func_0_0(86)) +
              (occ_func_0_1(3) * occ_func_0_1(8) * occ_func_0_0(42)) +
              (occ_func_0_1(2) * occ_func_0_1(5) * occ_func_0_0(41)) +
              (occ_func_0_1(6) * occ_func_0_1(4) * occ_func_0_0(80)) +
              (occ_func_0_1(7) * occ_func_0_1(5) * occ_func_0_0(22)) +
              (occ_func_0_1(9) * occ_func_0_1(8) * occ_func_0_0(24)) +
              (occ_func_0_1(9) * occ_func_0_1(8) * occ_func_0_0(82)) +
              (occ_func_0_1(4) * occ_func_0_1(6) * occ_func_0_0(26)) +
              (occ_func_0_1(5) * occ_func_0_1(7) * occ_func_0_0(27)) +
              (occ_func_0_1(2) * occ_func_0_1(3) * occ_func_0_0(79)) +
              (occ_func_0_1(11) * occ_func_0_1(8) * occ_func_0_0(20)) +
              (occ_func_0_1(10) * occ_func_0_1(5) * occ_func_0_0(19)) +
              (occ_func_0_1(11) * occ_func_0_1(6) * occ_func_0_0(84)) +
              (occ_func_0_1(2) * occ_func_0_1(1) * occ_func_0_0(31)) +
              (occ_func_0_1(7) * occ_func_0_1(12) * occ_func_0_0(38)) +
              (occ_func_0_1(5) * occ_func_0_1(7) * occ_func_0_0(81)) +
              (occ_func_0_1(8) * occ_func_0_1(9) * occ_func_0_0(25)) +
              (occ_func_0_1(6) * occ_func_0_1(4) * occ_func_0_0(23)) +
              (occ_func_0_1(1) * occ_func_0_1(6) * occ_func_0_0(80)) +
              (occ_func_0_1(12) * occ_func_0_1(11) * occ_func_0_0(29)) +
              (occ_func_0_1(7) * occ_func_0_1(2) * occ_func_0_0(22)) +
              (occ_func_0_1(4) * occ_func_0_1(10) * occ_func_0_0(83)) +
              (occ_func_0_1(9) * occ_func_0_1(12) * occ_func_0_0(35)) +
              (occ_func_0_1(3) * occ_func_0_1(1) * occ_func_0_0(28)) +
              (occ_func_0_1(4) * occ_func_0_1(1) * occ_func_0_0(80)) +
              (occ_func_0_1(9) * occ_func_0_1(3) * occ_func_0_0(24)) +
              (occ_func_0_1(12) * occ_func_0_1(10) * occ_func_0_0(29))) /
             24.0;
}
double test_Clexulator::delta_site_eval_at_0_bfunc_4_0_6(int occ_i,
                                                         int occ_f) const {
  return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) *
             ((occ_func_0_0(9) * occ_func_0_0(7) * occ_func_0_1(85)) +
              (occ_func_0_0(4) * occ_func_0_0(5) * occ_func_0_1(37)) +
              (occ_func_0_0(6) * occ_func_0_0(8) * occ_func_0_1(39)) +
              (occ_func_0_0(3) * occ_func_0_0(9) * occ_func_0_1(82)) +
              (occ_func_0_0(10) * occ_func_0_0(12) * occ_func_0_1(33)) +
              (occ_func_0_0(4) * occ_func_0_0(1) * occ_func_0_1(26)) +
              (occ_func_0_0(2) * occ_func_0_0(5) * occ_func_0_1(81)) +
              (occ_func_0_0(11) * occ_func_0_0(10) * occ_func_0_1(30)) +
              (occ_func_0_0(8) * occ_func_0_0(3) * occ_func_0_1(25)) +
              (occ_func_0_0(12) * occ_func_0_0(10) * occ_func_0_1(86)) +
              (occ_func_0_0(1) * occ_func_0_0(4) * occ_func_0_1(40)) +
              (occ_func_0_0(3) * occ_func_0_0(9) * occ_func_0_1(42)) +
              (occ_func_0_0(8) * occ_func_0_0(3) * occ_func_0_1(82)) +
              (occ_func_0_0(5) * occ_func_0_0(2) * occ_func_0_1(27)) +
              (occ_func_0_0(10) * occ_func_0_0(11) * occ_func_0_1(33)) +
              (occ_func_0_0(11) * occ_func_0_0(12) * occ_func_0_1(86)) +
              (occ_func_0_0(2) * occ_func_0_0(7) * occ_func_0_1(41)) +
              (occ_func_0_0(1) * occ_func_0_0(6) * occ_func_0_1(40)) +
              (occ_func_0_0(7) * occ_func_0_0(2) * occ_func_0_1(81)) +
              (occ_func_0_0(6) * occ_func_0_0(1) * occ_func_0_1(23)) +
              (occ_func_0_0(11) * occ_func_0_0(12) * occ_func_0_1(30)) +
              (occ_func_0_0(6) * occ_func_0_0(8) * occ_func_0_1(84)) +
              (occ_func_0_0(7) * occ_func_0_0(9) * occ_func_0_1(38)) +
              (occ_func_0_0(5) * occ_func_0_0(4) * occ_func_0_1(36)) +
              (occ_func_0_0(8) * occ_func_0_0(11) * occ_func_0_1(84)) +
              (occ_func_0_0(5) * occ_func_0_0(10) * occ_func_0_1(36)) +
              (occ_func_0_0(2) * occ_func_0_0(3) * occ_func_0_1(31)) +
              (occ_func_0_0(12) * occ_func_0_0(9) * occ_func_0_1(85)) +
              (occ_func_0_0(1) * occ_func_0_0(3) * occ_func_0_1(32)) +
              (occ_func_0_0(4) * occ_func_0_0(10) * occ_func_0_1(37)) +
              (occ_func_0_0(5) * occ_func_0_0(4) * occ_func_0_1(83)) +
              (occ_func_0_0(8) * occ_func_0_0(6) * occ_func_0_1(34)) +
              (occ_func_0_0(9) * occ_func_0_0(7) * occ_func_0_1(35)) +
              (occ_func_0_0(7) * occ_func_0_0(12) * occ_func_0_1(85)) +
              (occ_func_0_0(6) * occ_func_0_0(11) * occ_func_0_1(39)) +
              (occ_func_0_0(1) * occ_func_0_0(2) * occ_func_0_1(32)) +
              (occ_func_0_0(1) * occ_func_0_0(2) * occ_func_0_1(79)) +
              (occ_func_0_0(12) * occ_func_0_0(7) * occ_func_0_1(21)) +
              (occ_func_0_0(11) * occ_func_0_0(6) * occ_func_0_1(20)) +
              (occ_func_0_0(3) * occ_func_0_0(1) * occ_func_0_1(79)) +
              (occ_func_0_0(10) * occ_func_0_0(4) * occ_func_0_1(19)) +
              (occ_func_0_0(12) * occ_func_0_0(9) * occ_func_0_1(21)) +
              (occ_func_0_0(10) * occ_func_0_0(5) * occ_func_0_1(83)) +
              (occ_func_0_0(3) * occ_func_0_0(2) * occ_func_0_1(28)) +
              (occ_func_0_0(8) * occ_func_0_0(11) * occ_func_0_1(34)) +
              (occ_func_0_0(10) * occ_func_0_0(11) * occ_func_0_1(86)) +
              (occ_func_0_0(3) * occ_func_0_0(8) * occ_func_0_1(42)) +
              (occ_func_0_0(2) * occ_func_0_0(5) * occ_func_0_1(41)) +
              (occ_func_0_0(6) * occ_func_0_0(4) * occ_func_0_1(80)) +
              (occ_func_0_0(7) * occ_func_0_0(5) * occ_func_0_1(22)) +
              (occ_func_0_0(9) * occ_func_0_0(8) * occ_func_0_1(24)) +
              (occ_func_0_0(9) * occ_func_0_0(8) * occ_func_0_1(82)) +
              (occ_func_0_0(4) * occ_func_0_0(6) * occ_func_0_1(26)) +
              (occ_func_0_0(5) * occ_func_0_0(7) * occ_func_0_1(27)) +
              (occ_func_0_0(2) * occ_func_0_0(3) * occ_func_0_1(79)) +
              (occ_func_0_0(11) * occ_func_0_0(8) * occ_func_0_1(20)) +
              (occ_func_0_0(10) * occ_func_0_0(5) * occ_func_0_1(19)) +
              (occ_func_0_0(11) * occ_func_0_0(6) * occ_func_0_1(84)) +
              (occ_func_0_0(2) * occ_func_0_0(1) * occ_func_0_1(31)) +
              (occ_func_0_0(7) * occ_func_0_0(12) * occ_func_0_1(38)) +
              (occ_func_0_0(5) * occ_func_0_0(7) * occ_func_0_1(81)) +
              (occ_func_0_0(8) * occ_func_0_0(9) * occ_func_0_1(25)) +
              (occ_func_0_0(6) * occ_func_0_0(4) * occ_func_0_1(23)) +
              (occ_func_0_0(1) * occ_func_0_0(6) * occ_func_0_1(80)) +
              (occ_func_0_0(12) * occ_func_0_0(11) * occ_func_0_1(29)) +
              (occ_func_0_0(7) * occ_func_0_0(2) * occ_func_0_1(22)) +
              (occ_func_0_0(4) * occ_func_0_0(10) * occ_func_0_1(83)) +
              (occ_func_0_0(9) * occ_func_0_0(12) * occ_func_0_1(35)) +
              (occ_func_0_0(3) * occ_func_0_0(1) * occ_func_0_1(28)) +
              (occ_func_0_0(4) * occ_func_0_0(1) * occ_func_0_1(80)) +
              (occ_func_0_0(9) * occ_func_0_0(3) * occ_func_0_1(24)) +
              (occ_func_0_0(12) * occ_func_0_0(10) * occ_func_0_1(29))) /
             24.0 +
         (m_occ_func_0_1[occ_f] - m_occ_func_0_1[occ_i]) *
             ((occ_func_0_0(80) * occ_func_0_0(24) * occ_func_0_0(22)) +
              (occ_func_0_0(83) * occ_func_0_0(28) * occ_func_0_0(35)) +
              (occ_func_0_0(84) * occ_func_0_0(31) * occ_func_0_0(36)) +
              (occ_func_0_0(79) * occ_func_0_0(21) * occ_func_0_0(19)) +
              (occ_func_0_0(83) * occ_func_0_0(34) * occ_func_0_0(28)) +
              (occ_func_0_0(79) * occ_func_0_0(20) * occ_func_0_0(21)) +
              (occ_func_0_0(84) * occ_func_0_0(38) * occ_func_0_0(31)) +
              (occ_func_0_0(81) * occ_func_0_0(23) * occ_func_0_0(25)) +
              (occ_func_0_0(81) * occ_func_0_0(25) * occ_func_0_0(30)) +
              (occ_func_0_0(80) * occ_func_0_0(29) * occ_func_0_0(24)) +
              (occ_func_0_0(82) * occ_func_0_0(27) * occ_func_0_0(26)) +
              (occ_func_0_0(80) * occ_func_0_0(22) * occ_func_0_0(29)) +
              (occ_func_0_0(86) * occ_func_0_0(40) * occ_func_0_0(41)) +
              (occ_func_0_0(86) * occ_func_0_0(42) * occ_func_0_0(40)) +
              (occ_func_0_0(82) * occ_func_0_0(33) * occ_func_0_0(27)) +
              (occ_func_0_0(79) * occ_func_0_0(19) * occ_func_0_0(20)) +
              (occ_func_0_0(85) * occ_func_0_0(39) * occ_func_0_0(37)) +
              (occ_func_0_0(83) * occ_func_0_0(35) * occ_func_0_0(34)) +
              (occ_func_0_0(86) * occ_func_0_0(41) * occ_func_0_0(42)) +
              (occ_func_0_0(81) * occ_func_0_0(30) * occ_func_0_0(23)) +
              (occ_func_0_0(84) * occ_func_0_0(36) * occ_func_0_0(38)) +
              (occ_func_0_0(85) * occ_func_0_0(32) * occ_func_0_0(39)) +
              (occ_func_0_0(82) * occ_func_0_0(26) * occ_func_0_0(33)) +
              (occ_func_0_0(85) * occ_func_0_0(37) * occ_func_0_0(32))) /
             24.0;
}
double test_Clexulator::delta_site_eval_at_0_bfunc_4_0_7(int occ_i,
                                                         int occ_f) const {
  return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) *
             ((occ_func_0_1(4) * occ_func_0_0(5) * occ_func_0_1(37)) +
              (occ_func_0_1(6) * occ_func_0_0(8) * occ_func_0_1(39)) +
              (occ_func_0_1(10) * occ_func_0_0(12) * occ_func_0_1(33)) +
              (occ_func_0_1(4) * occ_func_0_0(1) * occ_func_0_1(26)) +
              (occ_func_0_1(11) * occ_func_0_0(10) * occ_func_0_1(30)) +
              (occ_func_0_1(8) * occ_func_0_0(3) * occ_func_0_1(25)) +
              (occ_func_0_1(1) * occ_func_0_0(4) * occ_func_0_1(40)) +
              (occ_func_0_1(3) * occ_func_0_0(9) * occ_func_0_1(42)) +
              (occ_func_0_1(5) * occ_func_0_0(2) * occ_func_0_1(27)) +
              (occ_func_0_1(10) * occ_func_0_0(11) * occ_func_0_1(33)) +
              (occ_func_0_1(2) * occ_func_0_0(7) * occ_func_0_1(41)) +
              (occ_func_0_1(1) * occ_func_0_0(6) * occ_func_0_1(40)) +
              (occ_func_0_1(6) * occ_func_0_0(1) * occ_func_0_1(23)) +
              (occ_func_0_1(11) * occ_func_0_0(12) * occ_func_0_1(30)) +
              (occ_func_0_1(7) * occ_func_0_0(9) * occ_func_0_1(38)) +
              (occ_func_0_1(5) * occ_func_0_0(4) * occ_func_0_1(36)) +
              (occ_func_0_1(5) * occ_func_0_0(10) * occ_func_0_1(36)) +
              (occ_func_0_1(2) * occ_func_0_0(3) * occ_func_0_1(31)) +
              (occ_func_0_1(1) * occ_func_0_0(3) * occ_func_0_1(32)) +
              (occ_func_0_1(4) * occ_func_0_0(10) * occ_func_0_1(37)) +
              (occ_func_0_1(8) * occ_func_0_0(6) * occ_func_0_1(34)) +
              (occ_func_0_1(9) * occ_func_0_0(7) * occ_func_0_1(35)) +
              (occ_func_0_1(6) * occ_func_0_0(11) * occ_func_0_1(39)) +
              (occ_func_0_1(1) * occ_func_0_0(2) * occ_func_0_1(32)) +
              (occ_func_0_1(12) * occ_func_0_0(7) * occ_func_0_1(21)) +
              (occ_func_0_1(11) * occ_func_0_0(6) * occ_func_0_1(20)) +
              (occ_func_0_1(10) * occ_func_0_0(4) * occ_func_0_1(19)) +
              (occ_func_0_1(12) * occ_func_0_0(9) * occ_func_0_1(21)) +
              (occ_func_0_1(3) * occ_func_0_0(2) * occ_func_0_1(28)) +
              (occ_func_0_1(8) * occ_func_0_0(11) * occ_func_0_1(34)) +
              (occ_func_0_1(3) * occ_func_0_0(8) * occ_func_0_1(42)) +
              (occ_func_0_1(2) * occ_func_0_0(5) * occ_func_0_1(41)) +
              (occ_func_0_1(7) * occ_func_0_0(5) * occ_func_0_1(22)) +
              (occ_func_0_1(9) * occ_func_0_0(8) * occ_func_0_1(24)) +
              (occ_func_0_1(4) * occ_func_0_0(6) * occ_func_0_1(26)) +
              (occ_func_0_1(5) * occ_func_0_0(7) * occ_func_0_1(27)) +
              (occ_func_0_1(11) * occ_func_0_0(8) * occ_func_0_1(20)) +
              (occ_func_0_1(10) * occ_func_0_0(5) * occ_func_0_1(19)) +
              (occ_func_0_1(2) * occ_func_0_0(1) * occ_func_0_1(31)) +
              (occ_func_0_1(7) * occ_func_0_0(12) * occ_func_0_1(38)) +
              (occ_func_0_1(8) * occ_func_0_0(9) * occ_func_0_1(25)) +
              (occ_func_0_1(6) * occ_func_0_0(4) * occ_func_0_1(23)) +
              (occ_func_0_1(12) * occ_func_0_0(11) * occ_func_0_1(29)) +
              (occ_func_0_1(7) * occ_func_0_0(2) * occ_func_0_1(22)) +
              (occ_func_0_1(9) * occ_func_0_0(12) * occ_func_0_1(35)) +
              (occ_func_0_1(3) * occ_func_0_0(1) * occ_func_0_1(28)) +
              (occ_func_0_1(9) * occ_func_0_0(3) * occ_func_0_1(24)) +
              (occ_func_0_1(12) * occ_func_0_0(10) * occ_func_0_1(29))) /
             24.0 +
         (m_occ_func_0_1[occ_f] - m_occ_func_0_1[occ_i]) *
             ((occ_func_0_0(9) * occ_func_0_0(7) * occ_func_0_1(85)) +
              (occ_func_0_1(80) * occ_func_0_0(24) * occ_func_0_0(22)) +
              (occ_func_0_0(3) * occ_func_0_0(9) * occ_func_0_1(82)) +
              (occ_func_0_1(83) * occ_func_0_0(28) * occ_func_0_0(35)) +
              (occ_func_0_0(2) * occ_func_0_0(5) * occ_func_0_1(81)) +
              (occ_func_0_1(84) * occ_func_0_0(31) * occ_func_0_0(36)) +
              (occ_func_0_0(12) * occ_func_0_0(10) * occ_func_0_1(86)) +
              (occ_func_0_1(79) * occ_func_0_0(21) * occ_func_0_0(19)) +
              (occ_func_0_0(8) * occ_func_0_0(3) * occ_func_0_1(82)) +
              (occ_func_0_1(83) * occ_func_0_0(34) * occ_func_0_0(28)) +
              (occ_func_0_0(11) * occ_func_0_0(12) * occ_func_0_1(86)) +
              (occ_func_0_1(79) * occ_func_0_0(20) * occ_func_0_0(21)) +
              (occ_func_0_0(7) * occ_func_0_0(2) * occ_func_0_1(81)) +
              (occ_func_0_1(84) * occ_func_0_0(38) * occ_func_0_0(31)) +
              (occ_func_0_0(6) * occ_func_0_0(8) * occ_func_0_1(84)) +
              (occ_func_0_1(81) * occ_func_0_0(23) * occ_func_0_0(25)) +
              (occ_func_0_0(8) * occ_func_0_0(11) * occ_func_0_1(84)) +
              (occ_func_0_1(81) * occ_func_0_0(25) * occ_func_0_0(30)) +
              (occ_func_0_0(12) * occ_func_0_0(9) * occ_func_0_1(85)) +
              (occ_func_0_1(80) * occ_func_0_0(29) * occ_func_0_0(24)) +
              (occ_func_0_0(5) * occ_func_0_0(4) * occ_func_0_1(83)) +
              (occ_func_0_1(82) * occ_func_0_0(27) * occ_func_0_0(26)) +
              (occ_func_0_0(7) * occ_func_0_0(12) * occ_func_0_1(85)) +
              (occ_func_0_1(80) * occ_func_0_0(22) * occ_func_0_0(29)) +
              (occ_func_0_0(1) * occ_func_0_0(2) * occ_func_0_1(79)) +
              (occ_func_0_1(86) * occ_func_0_0(40) * occ_func_0_0(41)) +
              (occ_func_0_0(3) * occ_func_0_0(1) * occ_func_0_1(79)) +
              (occ_func_0_1(86) * occ_func_0_0(42) * occ_func_0_0(40)) +
              (occ_func_0_0(10) * occ_func_0_0(5) * occ_func_0_1(83)) +
              (occ_func_0_1(82) * occ_func_0_0(33) * occ_func_0_0(27)) +
              (occ_func_0_0(10) * occ_func_0_0(11) * occ_func_0_1(86)) +
              (occ_func_0_1(79) * occ_func_0_0(19) * occ_func_0_0(20)) +
              (occ_func_0_0(6) * occ_func_0_0(4) * occ_func_0_1(80)) +
              (occ_func_0_1(85) * occ_func_0_0(39) * occ_func_0_0(37)) +
              (occ_func_0_0(9) * occ_func_0_0(8) * occ_func_0_1(82)) +
              (occ_func_0_1(83) * occ_func_0_0(35) * occ_func_0_0(34)) +
              (occ_func_0_0(2) * occ_func_0_0(3) * occ_func_0_1(79)) +
              (occ_func_0_1(86) * occ_func_0_0(41) * occ_func_0_0(42)) +
              (occ_func_0_0(11) * occ_func_0_0(6) * occ_func_0_1(84)) +
              (occ_func_0_1(81) * occ_func_0_0(30) * occ_func_0_0(23)) +
              (occ_func_0_0(5) * occ_func_0_0(7) * occ_func_0_1(81)) +
              (occ_func_0_1(84) * occ_func_0_0(36) * occ_func_0_0(38)) +
              (occ_func_0_0(1) * occ_func_0_0(6) * occ_func_0_1(80)) +
              (occ_func_0_1(85) * occ_func_0_0(32) * occ_func_0_0(39)) +
              (occ_func_0_0(4) * occ_func_0_0(10) * occ_func_0_1(83)) +
              (occ_func_0_1(82) * occ_func_0_0(26) * occ_func_0_0(33)) +
              (occ_func_0_0(4) * occ_func_0_0(1) * occ_func_0_1(80)) +
              (occ_func_0_1(85) * occ_func_0_0(37) * occ_func_0_0(32))) /
             24.0;
}
double test_Clexulator::delta_site_eval_at_0_bfunc_4_0_8(int occ_i,
                                                         int occ_f) const {
  return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) *
             (((0.7071067812 * occ_func_0_1(9) * occ_func_0_0(7) *
                    occ_func_0_1(85) +
                0.7071067812 * occ_func_0_0(9) * occ_func_0_1(7) *
                    occ_func_0_1(85))) +
              (0.7071067812 * occ_func_0_0(4) * occ_func_0_1(5) *
               occ_func_0_1(37)) +
              (0.7071067812 * occ_func_0_0(6) * occ_func_0_1(8) *
               occ_func_0_1(39)) +
              ((0.7071067812 * occ_func_0_1(3) * occ_func_0_0(9) *
                    occ_func_0_1(82) +
                0.7071067812 * occ_func_0_0(3) * occ_func_0_1(9) *
                    occ_func_0_1(82))) +
              (0.7071067812 * occ_func_0_0(10) * occ_func_0_1(12) *
               occ_func_0_1(33)) +
              (0.7071067812 * occ_func_0_0(4) * occ_func_0_1(1) *
               occ_func_0_1(26)) +
              ((0.7071067812 * occ_func_0_1(2) * occ_func_0_0(5) *
                    occ_func_0_1(81) +
                0.7071067812 * occ_func_0_0(2) * occ_func_0_1(5) *
                    occ_func_0_1(81))) +
              (0.7071067812 * occ_func_0_0(11) * occ_func_0_1(10) *
               occ_func_0_1(30)) +
              (0.7071067812 * occ_func_0_0(8) * occ_func_0_1(3) *
               occ_func_0_1(25)) +
              ((0.7071067812 * occ_func_0_1(12) * occ_func_0_0(10) *
                    occ_func_0_1(86) +
                0.7071067812 * occ_func_0_0(12) * occ_func_0_1(10) *
                    occ_func_0_1(86))) +
              (0.7071067812 * occ_func_0_0(1) * occ_func_0_1(4) *
               occ_func_0_1(40)) +
              (0.7071067812 * occ_func_0_0(3) * occ_func_0_1(9) *
               occ_func_0_1(42)) +
              ((0.7071067812 * occ_func_0_1(8) * occ_func_0_0(3) *
                    occ_func_0_1(82) +
                0.7071067812 * occ_func_0_0(8) * occ_func_0_1(3) *
                    occ_func_0_1(82))) +
              (0.7071067812 * occ_func_0_0(5) * occ_func_0_1(2) *
               occ_func_0_1(27)) +
              (0.7071067812 * occ_func_0_0(10) * occ_func_0_1(11) *
               occ_func_0_1(33)) +
              ((0.7071067812 * occ_func_0_1(11) * occ_func_0_0(12) *
                    occ_func_0_1(86) +
                0.7071067812 * occ_func_0_0(11) * occ_func_0_1(12) *
                    occ_func_0_1(86))) +
              (0.7071067812 * occ_func_0_0(2) * occ_func_0_1(7) *
               occ_func_0_1(41)) +
              (0.7071067812 * occ_func_0_0(1) * occ_func_0_1(6) *
               occ_func_0_1(40)) +
              ((0.7071067812 * occ_func_0_1(7) * occ_func_0_0(2) *
                    occ_func_0_1(81) +
                0.7071067812 * occ_func_0_0(7) * occ_func_0_1(2) *
                    occ_func_0_1(81))) +
              (0.7071067812 * occ_func_0_0(6) * occ_func_0_1(1) *
               occ_func_0_1(23)) +
              (0.7071067812 * occ_func_0_0(11) * occ_func_0_1(12) *
               occ_func_0_1(30)) +
              ((0.7071067812 * occ_func_0_1(6) * occ_func_0_0(8) *
                    occ_func_0_1(84) +
                0.7071067812 * occ_func_0_0(6) * occ_func_0_1(8) *
                    occ_func_0_1(84))) +
              (0.7071067812 * occ_func_0_0(7) * occ_func_0_1(9) *
               occ_func_0_1(38)) +
              (0.7071067812 * occ_func_0_0(5) * occ_func_0_1(4) *
               occ_func_0_1(36)) +
              ((0.7071067812 * occ_func_0_1(8) * occ_func_0_0(11) *
                    occ_func_0_1(84) +
                0.7071067812 * occ_func_0_0(8) * occ_func_0_1(11) *
                    occ_func_0_1(84))) +
              (0.7071067812 * occ_func_0_0(5) * occ_func_0_1(10) *
               occ_func_0_1(36)) +
              (0.7071067812 * occ_func_0_0(2) * occ_func_0_1(3) *
               occ_func_0_1(31)) +
              ((0.7071067812 * occ_func_0_1(12) * occ_func_0_0(9) *
                    occ_func_0_1(85) +
                0.7071067812 * occ_func_0_0(12) * occ_func_0_1(9) *
                    occ_func_0_1(85))) +
              (0.7071067812 * occ_func_0_0(1) * occ_func_0_1(3) *
               occ_func_0_1(32)) +
              (0.7071067812 * occ_func_0_0(4) * occ_func_0_1(10) *
               occ_func_0_1(37)) +
              ((0.7071067812 * occ_func_0_1(5) * occ_func_0_0(4) *
                    occ_func_0_1(83) +
                0.7071067812 * occ_func_0_0(5) * occ_func_0_1(4) *
                    occ_func_0_1(83))) +
              (0.7071067812 * occ_func_0_0(8) * occ_func_0_1(6) *
               occ_func_0_1(34)) +
              (0.7071067812 * occ_func_0_0(9) * occ_func_0_1(7) *
               occ_func_0_1(35)) +
              ((0.7071067812 * occ_func_0_1(7) * occ_func_0_0(12) *
                    occ_func_0_1(85) +
                0.7071067812 * occ_func_0_0(7) * occ_func_0_1(12) *
                    occ_func_0_1(85))) +
              (0.7071067812 * occ_func_0_0(6) * occ_func_0_1(11) *
               occ_func_0_1(39)) +
              (0.7071067812 * occ_func_0_0(1) * occ_func_0_1(2) *
               occ_func_0_1(32)) +
              ((0.7071067812 * occ_func_0_1(1) * occ_func_0_0(2) *
                    occ_func_0_1(79) +
                0.7071067812 * occ_func_0_0(1) * occ_func_0_1(2) *
                    occ_func_0_1(79))) +
              (0.7071067812 * occ_func_0_0(12) * occ_func_0_1(7) *
               occ_func_0_1(21)) +
              (0.7071067812 * occ_func_0_0(11) * occ_func_0_1(6) *
               occ_func_0_1(20)) +
              ((0.7071067812 * occ_func_0_1(3) * occ_func_0_0(1) *
                    occ_func_0_1(79) +
                0.7071067812 * occ_func_0_0(3) * occ_func_0_1(1) *
                    occ_func_0_1(79))) +
              (0.7071067812 * occ_func_0_0(10) * occ_func_0_1(4) *
               occ_func_0_1(19)) +
              (0.7071067812 * occ_func_0_0(12) * occ_func_0_1(9) *
               occ_func_0_1(21)) +
              ((0.7071067812 * occ_func_0_1(10) * occ_func_0_0(5) *
                    occ_func_0_1(83) +
                0.7071067812 * occ_func_0_0(10) * occ_func_0_1(5) *
                    occ_func_0_1(83))) +
              (0.7071067812 * occ_func_0_0(3) * occ_func_0_1(2) *
               occ_func_0_1(28)) +
              (0.7071067812 * occ_func_0_0(8) * occ_func_0_1(11) *
               occ_func_0_1(34)) +
              ((0.7071067812 * occ_func_0_1(10) * occ_func_0_0(11) *
                    occ_func_0_1(86) +
                0.7071067812 * occ_func_0_0(10) * occ_func_0_1(11) *
                    occ_func_0_1(86))) +
              (0.7071067812 * occ_func_0_0(3) * occ_func_0_1(8) *
               occ_func_0_1(42)) +
              (0.7071067812 * occ_func_0_0(2) * occ_func_0_1(5) *
               occ_func_0_1(41)) +
              ((0.7071067812 * occ_func_0_1(6) * occ_func_0_0(4) *
                    occ_func_0_1(80) +
                0.7071067812 * occ_func_0_0(6) * occ_func_0_1(4) *
                    occ_func_0_1(80))) +
              (0.7071067812 * occ_func_0_0(7) * occ_func_0_1(5) *
               occ_func_0_1(22)) +
              (0.7071067812 * occ_func_0_0(9) * occ_func_0_1(8) *
               occ_func_0_1(24)) +
              ((0.7071067812 * occ_func_0_1(9) * occ_func_0_0(8) *
                    occ_func_0_1(82) +
                0.7071067812 * occ_func_0_0(9) * occ_func_0_1(8) *
                    occ_func_0_1(82))) +
              (0.7071067812 * occ_func_0_0(4) * occ_func_0_1(6) *
               occ_func_0_1(26)) +
              (0.7071067812 * occ_func_0_0(5) * occ_func_0_1(7) *
               occ_func_0_1(27)) +
              ((0.7071067812 * occ_func_0_1(2) * occ_func_0_0(3) *
                    occ_func_0_1(79) +
                0.7071067812 * occ_func_0_0(2) * occ_func_0_1(3) *
                    occ_func_0_1(79))) +
              (0.7071067812 * occ_func_0_0(11) * occ_func_0_1(8) *
               occ_func_0_1(20)) +
              (0.7071067812 * occ_func_0_0(10) * occ_func_0_1(5) *
               occ_func_0_1(19)) +
              ((0.7071067812 * occ_func_0_1(11) * occ_func_0_0(6) *
                    occ_func_0_1(84) +
                0.7071067812 * occ_func_0_0(11) * occ_func_0_1(6) *
                    occ_func_0_1(84))) +
              (0.7071067812 * occ_func_0_0(2) * occ_func_0_1(1) *
               occ_func_0_1(31)) +
              (0.7071067812 * occ_func_0_0(7) * occ_func_0_1(12) *
               occ_func_0_1(38)) +
              ((0.7071067812 * occ_func_0_1(5) * occ_func_0_0(7) *
                    occ_func_0_1(81) +
                0.7071067812 * occ_func_0_0(5) * occ_func_0_1(7) *
                    occ_func_0_1(81))) +
              (0.7071067812 * occ_func_0_0(8) * occ_func_0_1(9) *
               occ_func_0_1(25)) +
              (0.7071067812 * occ_func_0_0(6) * occ_func_0_1(4) *
               occ_func_0_1(23)) +
              ((0.7071067812 * occ_func_0_1(1) * occ_func_0_0(6) *
                    occ_func_0_1(80) +
                0.7071067812 * occ_func_0_0(1) * occ_func_0_1(6) *
                    occ_func_0_1(80))) +
              (0.7071067812 * occ_func_0_0(12) * occ_func_0_1(11) *
               occ_func_0_1(29)) +
              (0.7071067812 * occ_func_0_0(7) * occ_func_0_1(2) *
               occ_func_0_1(22)) +
              ((0.7071067812 * occ_func_0_1(4) * occ_func_0_0(10) *
                    occ_func_0_1(83) +
                0.7071067812 * occ_func_0_0(4) * occ_func_0_1(10) *
                    occ_func_0_1(83))) +
              (0.7071067812 * occ_func_0_0(9) * occ_func_0_1(12) *
               occ_func_0_1(35)) +
              (0.7071067812 * occ_func_0_0(3) * occ_func_0_1(1) *
               occ_func_0_1(28)) +
              ((0.7071067812 * occ_func_0_1(4) * occ_func_0_0(1) *
                    occ_func_0_1(80) +
                0.7071067812 * occ_func_0_0(4) * occ_func_0_1(1) *
                    occ_func_0_1(80))) +
              (0.7071067812 * occ_func_0_0(9) * occ_func_0_1(3) *
               occ_func_0_1(24)) +
              (0.7071067812 * occ_func_0_0(12) * occ_func_0_1(10) *
               occ_func_0_1(29))) /
             24.0 +
         (m_occ_func_0_1[occ_f] - m_occ_func_0_1[occ_i]) *
             ((0.7071067812 * occ_func_0_0(4) * occ_func_0_0(5) *
               occ_func_0_1(37)) +
              (0.7071067812 * occ_func_0_0(6) * occ_func_0_0(8) *
               occ_func_0_1(39)) +
              ((0.7071067812 * occ_func_0_0(80) * occ_func_0_1(24) *
                    occ_func_0_0(22) +
                0.7071067812 * occ_func_0_0(80) * occ_func_0_0(24) *
                    occ_func_0_1(22))) +
              (0.7071067812 * occ_func_0_0(10) * occ_func_0_0(12) *
               occ_func_0_1(33)) +
              (0.7071067812 * occ_func_0_0(4) * occ_func_0_0(1) *
               occ_func_0_1(26)) +
              ((0.7071067812 * occ_func_0_0(83) * occ_func_0_1(28) *
                    occ_func_0_0(35) +
                0.7071067812 * occ_func_0_0(83) * occ_func_0_0(28) *
                    occ_func_0_1(35))) +
              (0.7071067812 * occ_func_0_0(11) * occ_func_0_0(10) *
               occ_func_0_1(30)) +
              (0.7071067812 * occ_func_0_0(8) * occ_func_0_0(3) *
               occ_func_0_1(25)) +
              ((0.7071067812 * occ_func_0_0(84) * occ_func_0_1(31) *
                    occ_func_0_0(36) +
                0.7071067812 * occ_func_0_0(84) * occ_func_0_0(31) *
                    occ_func_0_1(36))) +
              (0.7071067812 * occ_func_0_0(1) * occ_func_0_0(4) *
               occ_func_0_1(40)) +
              (0.7071067812 * occ_func_0_0(3) * occ_func_0_0(9) *
               occ_func_0_1(42)) +
              ((0.7071067812 * occ_func_0_0(79) * occ_func_0_1(21) *
                    occ_func_0_0(19) +
                0.7071067812 * occ_func_0_0(79) * occ_func_0_0(21) *
                    occ_func_0_1(19))) +
              (0.7071067812 * occ_func_0_0(5) * occ_func_0_0(2) *
               occ_func_0_1(27)) +
              (0.7071067812 * occ_func_0_0(10) * occ_func_0_0(11) *
               occ_func_0_1(33)) +
              ((0.7071067812 * occ_func_0_0(83) * occ_func_0_1(34) *
                    occ_func_0_0(28) +
                0.7071067812 * occ_func_0_0(83) * occ_func_0_0(34) *
                    occ_func_0_1(28))) +
              (0.7071067812 * occ_func_0_0(2) * occ_func_0_0(7) *
               occ_func_0_1(41)) +
              (0.7071067812 * occ_func_0_0(1) * occ_func_0_0(6) *
               occ_func_0_1(40)) +
              ((0.7071067812 * occ_func_0_0(79) * occ_func_0_1(20) *
                    occ_func_0_0(21) +
                0.7071067812 * occ_func_0_0(79) * occ_func_0_0(20) *
                    occ_func_0_1(21))) +
              (0.7071067812 * occ_func_0_0(6) * occ_func_0_0(1) *
               occ_func_0_1(23)) +
              (0.7071067812 * occ_func_0_0(11) * occ_func_0_0(12) *
               occ_func_0_1(30)) +
              ((0.7071067812 * occ_func_0_0(84) * occ_func_0_1(38) *
                    occ_func_0_0(31) +
                0.7071067812 * occ_func_0_0(84) * occ_func_0_0(38) *
                    occ_func_0_1(31))) +
              (0.7071067812 * occ_func_0_0(7) * occ_func_0_0(9) *
               occ_func_0_1(38)) +
              (0.7071067812 * occ_func_0_0(5) * occ_func_0_0(4) *
               occ_func_0_1(36)) +
              ((0.7071067812 * occ_func_0_0(81) * occ_func_0_1(23) *
                    occ_func_0_0(25) +
                0.7071067812 * occ_func_0_0(81) * occ_func_0_0(23) *
                    occ_func_0_1(25))) +
              (0.7071067812 * occ_func_0_0(5) * occ_func_0_0(10) *
               occ_func_0_1(36)) +
              (0.7071067812 * occ_func_0_0(2) * occ_func_0_0(3) *
               occ_func_0_1(31)) +
              ((0.7071067812 * occ_func_0_0(81) * occ_func_0_1(25) *
                    occ_func_0_0(30) +
                0.7071067812 * occ_func_0_0(81) * occ_func_0_0(25) *
                    occ_func_0_1(30))) +
              (0.7071067812 * occ_func_0_0(1) * occ_func_0_0(3) *
               occ_func_0_1(32)) +
              (0.7071067812 * occ_func_0_0(4) * occ_func_0_0(10) *
               occ_func_0_1(37)) +
              ((0.7071067812 * occ_func_0_0(80) * occ_func_0_1(29) *
                    occ_func_0_0(24) +
                0.7071067812 * occ_func_0_0(80) * occ_func_0_0(29) *
                    occ_func_0_1(24))) +
              (0.7071067812 * occ_func_0_0(8) * occ_func_0_0(6) *
               occ_func_0_1(34)) +
              (0.7071067812 * occ_func_0_0(9) * occ_func_0_0(7) *
               occ_func_0_1(35)) +
              ((0.7071067812 * occ_func_0_0(82) * occ_func_0_1(27) *
                    occ_func_0_0(26) +
                0.7071067812 * occ_func_0_0(82) * occ_func_0_0(27) *
                    occ_func_0_1(26))) +
              (0.7071067812 * occ_func_0_0(6) * occ_func_0_0(11) *
               occ_func_0_1(39)) +
              (0.7071067812 * occ_func_0_0(1) * occ_func_0_0(2) *
               occ_func_0_1(32)) +
              ((0.7071067812 * occ_func_0_0(80) * occ_func_0_1(22) *
                    occ_func_0_0(29) +
                0.7071067812 * occ_func_0_0(80) * occ_func_0_0(22) *
                    occ_func_0_1(29))) +
              (0.7071067812 * occ_func_0_0(12) * occ_func_0_0(7) *
               occ_func_0_1(21)) +
              (0.7071067812 * occ_func_0_0(11) * occ_func_0_0(6) *
               occ_func_0_1(20)) +
              ((0.7071067812 * occ_func_0_0(86) * occ_func_0_1(40) *
                    occ_func_0_0(41) +
                0.7071067812 * occ_func_0_0(86) * occ_func_0_0(40) *
                    occ_func_0_1(41))) +
              (0.7071067812 * occ_func_0_0(10) * occ_func_0_0(4) *
               occ_func_0_1(19)) +
              (0.7071067812 * occ_func_0_0(12) * occ_func_0_0(9) *
               occ_func_0_1(21)) +
              ((0.7071067812 * occ_func_0_0(86) * occ_func_0_1(42) *
                    occ_func_0_0(40) +
                0.7071067812 * occ_func_0_0(86) * occ_func_0_0(42) *
                    occ_func_0_1(40))) +
              (0.7071067812 * occ_func_0_0(3) * occ_func_0_0(2) *
               occ_func_0_1(28)) +
              (0.7071067812 * occ_func_0_0(8) * occ_func_0_0(11) *
               occ_func_0_1(34)) +
              ((0.7071067812 * occ_func_0_0(82) * occ_func_0_1(33) *
                    occ_func_0_0(27) +
                0.7071067812 * occ_func_0_0(82) * occ_func_0_0(33) *
                    occ_func_0_1(27))) +
              (0.7071067812 * occ_func_0_0(3) * occ_func_0_0(8) *
               occ_func_0_1(42)) +
              (0.7071067812 * occ_func_0_0(2) * occ_func_0_0(5) *
               occ_func_0_1(41)) +
              ((0.7071067812 * occ_func_0_0(79) * occ_func_0_1(19) *
                    occ_func_0_0(20) +
                0.7071067812 * occ_func_0_0(79) * occ_func_0_0(19) *
                    occ_func_0_1(20))) +
              (0.7071067812 * occ_func_0_0(7) * occ_func_0_0(5) *
               occ_func_0_1(22)) +
              (0.7071067812 * occ_func_0_0(9) * occ_func_0_0(8) *
               occ_func_0_1(24)) +
              ((0.7071067812 * occ_func_0_0(85) * occ_func_0_1(39) *
                    occ_func_0_0(37) +
                0.7071067812 * occ_func_0_0(85) * occ_func_0_0(39) *
                    occ_func_0_1(37))) +
              (0.7071067812 * occ_func_0_0(4) * occ_func_0_0(6) *
               occ_func_0_1(26)) +
              (0.7071067812 * occ_func_0_0(5) * occ_func_0_0(7) *
               occ_func_0_1(27)) +
              ((0.7071067812 * occ_func_0_0(83) * occ_func_0_1(35) *
                    occ_func_0_0(34) +
                0.7071067812 * occ_func_0_0(83) * occ_func_0_0(35) *
                    occ_func_0_1(34))) +
              (0.7071067812 * occ_func_0_0(11) * occ_func_0_0(8) *
               occ_func_0_1(20)) +
              (0.7071067812 * occ_func_0_0(10) * occ_func_0_0(5) *
               occ_func_0_1(19)) +
              ((0.7071067812 * occ_func_0_0(86) * occ_func_0_1(41) *
                    occ_func_0_0(42) +
                0.7071067812 * occ_func_0_0(86) * occ_func_0_0(41) *
                    occ_func_0_1(42))) +
              (0.7071067812 * occ_func_0_0(2) * occ_func_0_0(1) *
               occ_func_0_1(31)) +
              (0.7071067812 * occ_func_0_0(7) * occ_func_0_0(12) *
               occ_func_0_1(38)) +
              ((0.7071067812 * occ_func_0_0(81) * occ_func_0_1(30) *
                    occ_func_0_0(23) +
                0.7071067812 * occ_func_0_0(81) * occ_func_0_0(30) *
                    occ_func_0_1(23))) +
              (0.7071067812 * occ_func_0_0(8) * occ_func_0_0(9) *
               occ_func_0_1(25)) +
              (0.7071067812 * occ_func_0_0(6) * occ_func_0_0(4) *
               occ_func_0_1(23)) +
              ((0.7071067812 * occ_func_0_0(84) * occ_func_0_1(36) *
                    occ_func_0_0(38) +
                0.7071067812 * occ_func_0_0(84) * occ_func_0_0(36) *
                    occ_func_0_1(38))) +
              (0.7071067812 * occ_func_0_0(12) * occ_func_0_0(11) *
               occ_func_0_1(29)) +
              (0.7071067812 * occ_func_0_0(7) * occ_func_0_0(2) *
               occ_func_0_1(22)) +
              ((0.7071067812 * occ_func_0_0(85) * occ_func_0_1(32) *
                    occ_func_0_0(39) +
                0.7071067812 * occ_func_0_0(85) * occ_func_0_0(32) *
                    occ_func_0_1(39))) +
              (0.7071067812 * occ_func_0_0(9) * occ_func_0_0(12) *
               occ_func_0_1(35)) +
              (0.7071067812 * occ_func_0_0(3) * occ_func_0_0(1) *
               occ_func_0_1(28)) +
              ((0.7071067812 * occ_func_0_0(82) * occ_func_0_1(26) *
                    occ_func_0_0(33) +
                0.7071067812 * occ_func_0_0(82) * occ_func_0_0(26) *
                    occ_func_0_1(33))) +
              (0.7071067812 * occ_func_0_0(9) * occ_func_0_0(3) *
               occ_func_0_1(24)) +
              (0.7071067812 * occ_func_0_0(12) * occ_func_0_0(10) *
               occ_func_0_1(29)) +
              ((0.7071067812 * occ_func_0_0(85) * occ_func_0_1(37) *
                    occ_func_0_0(32) +
                0.7071067812 * occ_func_0_0(85) * occ_func_0_0(37) *
                    occ_func_0_1(32)))) /
             24.0;
}
double test_Clexulator::delta_site_eval_at_0_bfunc_4_0_9(int occ_i,
                                                         int occ_f) const {
  return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) *
             ((0.7071067812 * occ_func_0_1(4) * occ_func_0_1(5) *
               occ_func_0_1(37)) +
              (0.7071067812 * occ_func_0_1(6) * occ_func_0_1(8) *
               occ_func_0_1(39)) +
              (0.7071067812 * occ_func_0_1(10) * occ_func_0_1(12) *
               occ_func_0_1(33)) +
              (0.7071067812 * occ_func_0_1(4) * occ_func_0_1(1) *
               occ_func_0_1(26)) +
              (0.7071067812 * occ_func_0_1(11) * occ_func_0_1(10) *
               occ_func_0_1(30)) +
              (0.7071067812 * occ_func_0_1(8) * occ_func_0_1(3) *
               occ_func_0_1(25)) +
              (0.7071067812 * occ_func_0_1(1) * occ_func_0_1(4) *
               occ_func_0_1(40)) +
              (0.7071067812 * occ_func_0_1(3) * occ_func_0_1(9) *
               occ_func_0_1(42)) +
              (0.7071067812 * occ_func_0_1(5) * occ_func_0_1(2) *
               occ_func_0_1(27)) +
              (0.7071067812 * occ_func_0_1(10) * occ_func_0_1(11) *
               occ_func_0_1(33)) +
              (0.7071067812 * occ_func_0_1(2) * occ_func_0_1(7) *
               occ_func_0_1(41)) +
              (0.7071067812 * occ_func_0_1(1) * occ_func_0_1(6) *
               occ_func_0_1(40)) +
              (0.7071067812 * occ_func_0_1(6) * occ_func_0_1(1) *
               occ_func_0_1(23)) +
              (0.7071067812 * occ_func_0_1(11) * occ_func_0_1(12) *
               occ_func_0_1(30)) +
              (0.7071067812 * occ_func_0_1(7) * occ_func_0_1(9) *
               occ_func_0_1(38)) +
              (0.7071067812 * occ_func_0_1(5) * occ_func_0_1(4) *
               occ_func_0_1(36)) +
              (0.7071067812 * occ_func_0_1(5) * occ_func_0_1(10) *
               occ_func_0_1(36)) +
              (0.7071067812 * occ_func_0_1(2) * occ_func_0_1(3) *
               occ_func_0_1(31)) +
              (0.7071067812 * occ_func_0_1(1) * occ_func_0_1(3) *
               occ_func_0_1(32)) +
              (0.7071067812 * occ_func_0_1(4) * occ_func_0_1(10) *
               occ_func_0_1(37)) +
              (0.7071067812 * occ_func_0_1(8) * occ_func_0_1(6) *
               occ_func_0_1(34)) +
              (0.7071067812 * occ_func_0_1(9) * occ_func_0_1(7) *
               occ_func_0_1(35)) +
              (0.7071067812 * occ_func_0_1(6) * occ_func_0_1(11) *
               occ_func_0_1(39)) +
              (0.7071067812 * occ_func_0_1(1) * occ_func_0_1(2) *
               occ_func_0_1(32)) +
              (0.7071067812 * occ_func_0_1(12) * occ_func_0_1(7) *
               occ_func_0_1(21)) +
              (0.7071067812 * occ_func_0_1(11) * occ_func_0_1(6) *
               occ_func_0_1(20)) +
              (0.7071067812 * occ_func_0_1(10) * occ_func_0_1(4) *
               occ_func_0_1(19)) +
              (0.7071067812 * occ_func_0_1(12) * occ_func_0_1(9) *
               occ_func_0_1(21)) +
              (0.7071067812 * occ_func_0_1(3) * occ_func_0_1(2) *
               occ_func_0_1(28)) +
              (0.7071067812 * occ_func_0_1(8) * occ_func_0_1(11) *
               occ_func_0_1(34)) +
              (0.7071067812 * occ_func_0_1(3) * occ_func_0_1(8) *
               occ_func_0_1(42)) +
              (0.7071067812 * occ_func_0_1(2) * occ_func_0_1(5) *
               occ_func_0_1(41)) +
              (0.7071067812 * occ_func_0_1(7) * occ_func_0_1(5) *
               occ_func_0_1(22)) +
              (0.7071067812 * occ_func_0_1(9) * occ_func_0_1(8) *
               occ_func_0_1(24)) +
              (0.7071067812 * occ_func_0_1(4) * occ_func_0_1(6) *
               occ_func_0_1(26)) +
              (0.7071067812 * occ_func_0_1(5) * occ_func_0_1(7) *
               occ_func_0_1(27)) +
              (0.7071067812 * occ_func_0_1(11) * occ_func_0_1(8) *
               occ_func_0_1(20)) +
              (0.7071067812 * occ_func_0_1(10) * occ_func_0_1(5) *
               occ_func_0_1(19)) +
              (0.7071067812 * occ_func_0_1(2) * occ_func_0_1(1) *
               occ_func_0_1(31)) +
              (0.7071067812 * occ_func_0_1(7) * occ_func_0_1(12) *
               occ_func_0_1(38)) +
              (0.7071067812 * occ_func_0_1(8) * occ_func_0_1(9) *
               occ_func_0_1(25)) +
              (0.7071067812 * occ_func_0_1(6) * occ_func_0_1(4) *
               occ_func_0_1(23)) +
              (0.7071067812 * occ_func_0_1(12) * occ_func_0_1(11) *
               occ_func_0_1(29)) +
              (0.7071067812 * occ_func_0_1(7) * occ_func_0_1(2) *
               occ_func_0_1(22)) +
              (0.7071067812 * occ_func_0_1(9) * occ_func_0_1(12) *
               occ_func_0_1(35)) +
              (0.7071067812 * occ_func_0_1(3) * occ_func_0_1(1) *
               occ_func_0_1(28)) +
              (0.7071067812 * occ_func_0_1(9) * occ_func_0_1(3) *
               occ_func_0_1(24)) +
              (0.7071067812 * occ_func_0_1(12) * occ_func_0_1(10) *
               occ_func_0_1(29))) /
             24.0 +
         (m_occ_func_0_1[occ_f] - m_occ_func_0_1[occ_i]) *
             (((0.7071067812 * occ_func_0_1(9) * occ_func_0_0(7) *
                    occ_func_0_1(85) +
                0.7071067812 * occ_func_0_0(9) * occ_func_0_1(7) *
                    occ_func_0_1(85))) +
              (0.7071067812 * occ_func_0_1(4) * occ_func_0_0(5) *
               occ_func_0_1(37)) +
              (0.7071067812 * occ_func_0_1(6) * occ_func_0_0(8) *
               occ_func_0_1(39)) +
              ((0.7071067812 * occ_func_0_1(80) * occ_func_0_1(24) *
                    occ_func_0_0(22) +
                0.7071067812 * occ_func_0_1(80) * occ_func_0_0(24) *
                    occ_func_0_1(22))) +
              ((0.7071067812 * occ_func_0_1(3) * occ_func_0_0(9) *
                    occ_func_0_1(82) +
                0.7071067812 * occ_func_0_0(3) * occ_func_0_1(9) *
                    occ_func_0_1(82))) +
              (0.7071067812 * occ_func_0_1(10) * occ_func_0_0(12) *
               occ_func_0_1(33)) +
              (0.7071067812 * occ_func_0_1(4) * occ_func_0_0(1) *
               occ_func_0_1(26)) +
              ((0.7071067812 * occ_func_0_1(83) * occ_func_0_1(28) *
                    occ_func_0_0(35) +
                0.7071067812 * occ_func_0_1(83) * occ_func_0_0(28) *
                    occ_func_0_1(35))) +
              ((0.7071067812 * occ_func_0_1(2) * occ_func_0_0(5) *
                    occ_func_0_1(81) +
                0.7071067812 * occ_func_0_0(2) * occ_func_0_1(5) *
                    occ_func_0_1(81))) +
              (0.7071067812 * occ_func_0_1(11) * occ_func_0_0(10) *
               occ_func_0_1(30)) +
              (0.7071067812 * occ_func_0_1(8) * occ_func_0_0(3) *
               occ_func_0_1(25)) +
              ((0.7071067812 * occ_func_0_1(84) * occ_func_0_1(31) *
                    occ_func_0_0(36) +
                0.7071067812 * occ_func_0_1(84) * occ_func_0_0(31) *
                    occ_func_0_1(36))) +
              ((0.7071067812 * occ_func_0_1(12) * occ_func_0_0(10) *
                    occ_func_0_1(86) +
                0.7071067812 * occ_func_0_0(12) * occ_func_0_1(10) *
                    occ_func_0_1(86))) +
              (0.7071067812 * occ_func_0_1(1) * occ_func_0_0(4) *
               occ_func_0_1(40)) +
              (0.7071067812 * occ_func_0_1(3) * occ_func_0_0(9) *
               occ_func_0_1(42)) +
              ((0.7071067812 * occ_func_0_1(79) * occ_func_0_1(21) *
                    occ_func_0_0(19) +
                0.7071067812 * occ_func_0_1(79) * occ_func_0_0(21) *
                    occ_func_0_1(19))) +
              ((0.7071067812 * occ_func_0_1(8) * occ_func_0_0(3) *
                    occ_func_0_1(82) +
                0.7071067812 * occ_func_0_0(8) * occ_func_0_1(3) *
                    occ_func_0_1(82))) +
              (0.7071067812 * occ_func_0_1(5) * occ_func_0_0(2) *
               occ_func_0_1(27)) +
              (0.7071067812 * occ_func_0_1(10) * occ_func_0_0(11) *
               occ_func_0_1(33)) +
              ((0.7071067812 * occ_func_0_1(83) * occ_func_0_1(34) *
                    occ_func_0_0(28) +
                0.7071067812 * occ_func_0_1(83) * occ_func_0_0(34) *
                    occ_func_0_1(28))) +
              ((0.7071067812 * occ_func_0_1(11) * occ_func_0_0(12) *
                    occ_func_0_1(86) +
                0.7071067812 * occ_func_0_0(11) * occ_func_0_1(12) *
                    occ_func_0_1(86))) +
              (0.7071067812 * occ_func_0_1(2) * occ_func_0_0(7) *
               occ_func_0_1(41)) +
              (0.7071067812 * occ_func_0_1(1) * occ_func_0_0(6) *
               occ_func_0_1(40)) +
              ((0.7071067812 * occ_func_0_1(79) * occ_func_0_1(20) *
                    occ_func_0_0(21) +
                0.7071067812 * occ_func_0_1(79) * occ_func_0_0(20) *
                    occ_func_0_1(21))) +
              ((0.7071067812 * occ_func_0_1(7) * occ_func_0_0(2) *
                    occ_func_0_1(81) +
                0.7071067812 * occ_func_0_0(7) * occ_func_0_1(2) *
                    occ_func_0_1(81))) +
              (0.7071067812 * occ_func_0_1(6) * occ_func_0_0(1) *
               occ_func_0_1(23)) +
              (0.7071067812 * occ_func_0_1(11) * occ_func_0_0(12) *
               occ_func_0_1(30)) +
              ((0.7071067812 * occ_func_0_1(84) * occ_func_0_1(38) *
                    occ_func_0_0(31) +
                0.7071067812 * occ_func_0_1(84) * occ_func_0_0(38) *
                    occ_func_0_1(31))) +
              ((0.7071067812 * occ_func_0_1(6) * occ_func_0_0(8) *
                    occ_func_0_1(84) +
                0.7071067812 * occ_func_0_0(6) * occ_func_0_1(8) *
                    occ_func_0_1(84))) +
              (0.7071067812 * occ_func_0_1(7) * occ_func_0_0(9) *
               occ_func_0_1(38)) +
              (0.7071067812 * occ_func_0_1(5) * occ_func_0_0(4) *
               occ_func_0_1(36)) +
              ((0.7071067812 * occ_func_0_1(81) * occ_func_0_1(23) *
                    occ_func_0_0(25) +
                0.7071067812 * occ_func_0_1(81) * occ_func_0_0(23) *
                    occ_func_0_1(25))) +
              ((0.7071067812 * occ_func_0_1(8) * occ_func_0_0(11) *
                    occ_func_0_1(84) +
                0.7071067812 * occ_func_0_0(8) * occ_func_0_1(11) *
                    occ_func_0_1(84))) +
              (0.7071067812 * occ_func_0_1(5) * occ_func_0_0(10) *
               occ_func_0_1(36)) +
              (0.7071067812 * occ_func_0_1(2) * occ_func_0_0(3) *
               occ_func_0_1(31)) +
              ((0.7071067812 * occ_func_0_1(81) * occ_func_0_1(25) *
                    occ_func_0_0(30) +
                0.7071067812 * occ_func_0_1(81) * occ_func_0_0(25) *
                    occ_func_0_1(30))) +
              ((0.7071067812 * occ_func_0_1(12) * occ_func_0_0(9) *
                    occ_func_0_1(85) +
                0.7071067812 * occ_func_0_0(12) * occ_func_0_1(9) *
                    occ_func_0_1(85))) +
              (0.7071067812 * occ_func_0_1(1) * occ_func_0_0(3) *
               occ_func_0_1(32)) +
              (0.7071067812 * occ_func_0_1(4) * occ_func_0_0(10) *
               occ_func_0_1(37)) +
              ((0.7071067812 * occ_func_0_1(80) * occ_func_0_1(29) *
                    occ_func_0_0(24) +
                0.7071067812 * occ_func_0_1(80) * occ_func_0_0(29) *
                    occ_func_0_1(24))) +
              ((0.7071067812 * occ_func_0_1(5) * occ_func_0_0(4) *
                    occ_func_0_1(83) +
                0.7071067812 * occ_func_0_0(5) * occ_func_0_1(4) *
                    occ_func_0_1(83))) +
              (0.7071067812 * occ_func_0_1(8) * occ_func_0_0(6) *
               occ_func_0_1(34)) +
              (0.7071067812 * occ_func_0_1(9) * occ_func_0_0(7) *
               occ_func_0_1(35)) +
              ((0.7071067812 * occ_func_0_1(82) * occ_func_0_1(27) *
                    occ_func_0_0(26) +
                0.7071067812 * occ_func_0_1(82) * occ_func_0_0(27) *
                    occ_func_0_1(26))) +
              ((0.7071067812 * occ_func_0_1(7) * occ_func_0_0(12) *
                    occ_func_0_1(85) +
                0.7071067812 * occ_func_0_0(7) * occ_func_0_1(12) *
                    occ_func_0_1(85))) +
              (0.7071067812 * occ_func_0_1(6) * occ_func_0_0(11) *
               occ_func_0_1(39)) +
              (0.7071067812 * occ_func_0_1(1) * occ_func_0_0(2) *
               occ_func_0_1(32)) +
              ((0.7071067812 * occ_func_0_1(80) * occ_func_0_1(22) *
                    occ_func_0_0(29) +
                0.7071067812 * occ_func_0_1(80) * occ_func_0_0(22) *
                    occ_func_0_1(29))) +
              ((0.7071067812 * occ_func_0_1(1) * occ_func_0_0(2) *
                    occ_func_0_1(79) +
                0.7071067812 * occ_func_0_0(1) * occ_func_0_1(2) *
                    occ_func_0_1(79))) +
              (0.7071067812 * occ_func_0_1(12) * occ_func_0_0(7) *
               occ_func_0_1(21)) +
              (0.7071067812 * occ_func_0_1(11) * occ_func_0_0(6) *
               occ_func_0_1(20)) +
              ((0.7071067812 * occ_func_0_1(86) * occ_func_0_1(40) *
                    occ_func_0_0(41) +
                0.7071067812 * occ_func_0_1(86) * occ_func_0_0(40) *
                    occ_func_0_1(41))) +
              ((0.7071067812 * occ_func_0_1(3) * occ_func_0_0(1) *
                    occ_func_0_1(79) +
                0.7071067812 * occ_func_0_0(3) * occ_func_0_1(1) *
                    occ_func_0_1(79))) +
              (0.7071067812 * occ_func_0_1(10) * occ_func_0_0(4) *
               occ_func_0_1(19)) +
              (0.7071067812 * occ_func_0_1(12) * occ_func_0_0(9) *
               occ_func_0_1(21)) +
              ((0.7071067812 * occ_func_0_1(86) * occ_func_0_1(42) *
                    occ_func_0_0(40) +
                0.7071067812 * occ_func_0_1(86) * occ_func_0_0(42) *
                    occ_func_0_1(40))) +
              ((0.7071067812 * occ_func_0_1(10) * occ_func_0_0(5) *
                    occ_func_0_1(83) +
                0.7071067812 * occ_func_0_0(10) * occ_func_0_1(5) *
                    occ_func_0_1(83))) +
              (0.7071067812 * occ_func_0_1(3) * occ_func_0_0(2) *
               occ_func_0_1(28)) +
              (0.7071067812 * occ_func_0_1(8) * occ_func_0_0(11) *
               occ_func_0_1(34)) +
              ((0.7071067812 * occ_func_0_1(82) * occ_func_0_1(33) *
                    occ_func_0_0(27) +
                0.7071067812 * occ_func_0_1(82) * occ_func_0_0(33) *
                    occ_func_0_1(27))) +
              ((0.7071067812 * occ_func_0_1(10) * occ_func_0_0(11) *
                    occ_func_0_1(86) +
                0.7071067812 * occ_func_0_0(10) * occ_func_0_1(11) *
                    occ_func_0_1(86))) +
              (0.7071067812 * occ_func_0_1(3) * occ_func_0_0(8) *
               occ_func_0_1(42)) +
              (0.7071067812 * occ_func_0_1(2) * occ_func_0_0(5) *
               occ_func_0_1(41)) +
              ((0.7071067812 * occ_func_0_1(79) * occ_func_0_1(19) *
                    occ_func_0_0(20) +
                0.7071067812 * occ_func_0_1(79) * occ_func_0_0(19) *
                    occ_func_0_1(20))) +
              ((0.7071067812 * occ_func_0_1(6) * occ_func_0_0(4) *
                    occ_func_0_1(80) +
                0.7071067812 * occ_func_0_0(6) * occ_func_0_1(4) *
                    occ_func_0_1(80))) +
              (0.7071067812 * occ_func_0_1(7) * occ_func_0_0(5) *
               occ_func_0_1(22)) +
              (0.7071067812 * occ_func_0_1(9) * occ_func_0_0(8) *
               occ_func_0_1(24)) +
              ((0.7071067812 * occ_func_0_1(85) * occ_func_0_1(39) *
                    occ_func_0_0(37) +
                0.7071067812 * occ_func_0_1(85) * occ_func_0_0(39) *
                    occ_func_0_1(37))) +
              ((0.7071067812 * occ_func_0_1(9) * occ_func_0_0(8) *
                    occ_func_0_1(82) +
                0.7071067812 * occ_func_0_0(9) * occ_func_0_1(8) *
                    occ_func_0_1(82))) +
              (0.7071067812 * occ_func_0_1(4) * occ_func_0_0(6) *
               occ_func_0_1(26)) +
              (0.7071067812 * occ_func_0_1(5) * occ_func_0_0(7) *
               occ_func_0_1(27)) +
              ((0.7071067812 * occ_func_0_1(83) * occ_func_0_1(35) *
                    occ_func_0_0(34) +
                0.7071067812 * occ_func_0_1(83) * occ_func_0_0(35) *
                    occ_func_0_1(34))) +
              ((0.7071067812 * occ_func_0_1(2) * occ_func_0_0(3) *
                    occ_func_0_1(79) +
                0.7071067812 * occ_func_0_0(2) * occ_func_0_1(3) *
                    occ_func_0_1(79))) +
              (0.7071067812 * occ_func_0_1(11) * occ_func_0_0(8) *
               occ_func_0_1(20)) +
              (0.7071067812 * occ_func_0_1(10) * occ_func_0_0(5) *
               occ_func_0_1(19)) +
              ((0.7071067812 * occ_func_0_1(86) * occ_func_0_1(41) *
                    occ_func_0_0(42) +
                0.7071067812 * occ_func_0_1(86) * occ_func_0_0(41) *
                    occ_func_0_1(42))) +
              ((0.7071067812 * occ_func_0_1(11) * occ_func_0_0(6) *
                    occ_func_0_1(84) +
                0.7071067812 * occ_func_0_0(11) * occ_func_0_1(6) *
                    occ_func_0_1(84))) +
              (0.7071067812 * occ_func_0_1(2) * occ_func_0_0(1) *
               occ_func_0_1(31)) +
              (0.7071067812 * occ_func_0_1(7) * occ_func_0_0(12) *
               occ_func_0_1(38)) +
              ((0.7071067812 * occ_func_0_1(81) * occ_func_0_1(30) *
                    occ_func_0_0(23) +
                0.7071067812 * occ_func_0_1(81) * occ_func_0_0(30) *
                    occ_func_0_1(23))) +
              ((0.7071067812 * occ_func_0_1(5) * occ_func_0_0(7) *
                    occ_func_0_1(81) +
                0.7071067812 * occ_func_0_0(5) * occ_func_0_1(7) *
                    occ_func_0_1(81))) +
              (0.7071067812 * occ_func_0_1(8) * occ_func_0_0(9) *
               occ_func_0_1(25)) +
              (0.7071067812 * occ_func_0_1(6) * occ_func_0_0(4) *
               occ_func_0_1(23)) +
              ((0.7071067812 * occ_func_0_1(84) * occ_func_0_1(36) *
                    occ_func_0_0(38) +
                0.7071067812 * occ_func_0_1(84) * occ_func_0_0(36) *
                    occ_func_0_1(38))) +
              ((0.7071067812 * occ_func_0_1(1) * occ_func_0_0(6) *
                    occ_func_0_1(80) +
                0.7071067812 * occ_func_0_0(1) * occ_func_0_1(6) *
                    occ_func_0_1(80))) +
              (0.7071067812 * occ_func_0_1(12) * occ_func_0_0(11) *
               occ_func_0_1(29)) +
              (0.7071067812 * occ_func_0_1(7) * occ_func_0_0(2) *
               occ_func_0_1(22)) +
              ((0.7071067812 * occ_func_0_1(85) * occ_func_0_1(32) *
                    occ_func_0_0(39) +
                0.7071067812 * occ_func_0_1(85) * occ_func_0_0(32) *
                    occ_func_0_1(39))) +
              ((0.7071067812 * occ_func_0_1(4) * occ_func_0_0(10) *
                    occ_func_0_1(83) +
                0.7071067812 * occ_func_0_0(4) * occ_func_0_1(10) *
                    occ_func_0_1(83))) +
              (0.7071067812 * occ_func_0_1(9) * occ_func_0_0(12) *
               occ_func_0_1(35)) +
              (0.7071067812 * occ_func_0_1(3) * occ_func_0_0(1) *
               occ_func_0_1(28)) +
              ((0.7071067812 * occ_func_0_1(82) * occ_func_0_1(26) *
                    occ_func_0_0(33) +
                0.7071067812 * occ_func_0_1(82) * occ_func_0_0(26) *
                    occ_func_0_1(33))) +
              ((0.7071067812 * occ_func_0_1(4) * occ_func_0_0(1) *
                    occ_func_0_1(80) +
                0.7071067812 * occ_func_0_0(4) * occ_func_0_1(1) *
                    occ_func_0_1(80))) +
              (0.7071067812 * occ_func_0_1(9) * occ_func_0_0(3) *
               occ_func_0_1(24)) +
              (0.7071067812 * occ_func_0_1(12) * occ_func_0_0(10) *
               occ_func_0_1(29)) +
              ((0.7071067812 * occ_func_0_1(85) * occ_func_0_1(37) *
                    occ_func_0_0(32) +
                0.7071067812 * occ_func_0_1(85) * occ_func_0_0(37) *
                    occ_func_0_1(32)))) /
             24.0;
}
double test_Clexulator::delta_site_eval_at_0_bfunc_4_0_10(int occ_i,
                                                          int occ_f) const {
  return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) *
             ((occ_func_0_1(9) * occ_func_0_1(7) * occ_func_0_1(85)) +
              (occ_func_0_1(3) * occ_func_0_1(9) * occ_func_0_1(82)) +
              (occ_func_0_1(2) * occ_func_0_1(5) * occ_func_0_1(81)) +
              (occ_func_0_1(12) * occ_func_0_1(10) * occ_func_0_1(86)) +
              (occ_func_0_1(8) * occ_func_0_1(3) * occ_func_0_1(82)) +
              (occ_func_0_1(11) * occ_func_0_1(12) * occ_func_0_1(86)) +
              (occ_func_0_1(7) * occ_func_0_1(2) * occ_func_0_1(81)) +
              (occ_func_0_1(6) * occ_func_0_1(8) * occ_func_0_1(84)) +
              (occ_func_0_1(8) * occ_func_0_1(11) * occ_func_0_1(84)) +
              (occ_func_0_1(12) * occ_func_0_1(9) * occ_func_0_1(85)) +
              (occ_func_0_1(5) * occ_func_0_1(4) * occ_func_0_1(83)) +
              (occ_func_0_1(7) * occ_func_0_1(12) * occ_func_0_1(85)) +
              (occ_func_0_1(1) * occ_func_0_1(2) * occ_func_0_1(79)) +
              (occ_func_0_1(3) * occ_func_0_1(1) * occ_func_0_1(79)) +
              (occ_func_0_1(10) * occ_func_0_1(5) * occ_func_0_1(83)) +
              (occ_func_0_1(10) * occ_func_0_1(11) * occ_func_0_1(86)) +
              (occ_func_0_1(6) * occ_func_0_1(4) * occ_func_0_1(80)) +
              (occ_func_0_1(9) * occ_func_0_1(8) * occ_func_0_1(82)) +
              (occ_func_0_1(2) * occ_func_0_1(3) * occ_func_0_1(79)) +
              (occ_func_0_1(11) * occ_func_0_1(6) * occ_func_0_1(84)) +
              (occ_func_0_1(5) * occ_func_0_1(7) * occ_func_0_1(81)) +
              (occ_func_0_1(1) * occ_func_0_1(6) * occ_func_0_1(80)) +
              (occ_func_0_1(4) * occ_func_0_1(10) * occ_func_0_1(83)) +
              (occ_func_0_1(4) * occ_func_0_1(1) * occ_func_0_1(80))) /
             24.0 +
         (m_occ_func_0_1[occ_f] - m_occ_func_0_1[occ_i]) *
             ((occ_func_0_0(4) * occ_func_0_1(5) * occ_func_0_1(37)) +
              (occ_func_0_0(6) * occ_func_0_1(8) * occ_func_0_1(39)) +
              (occ_func_0_0(80) * occ_func_0_1(24) * occ_func_0_1(22)) +
              (occ_func_0_0(10) * occ_func_0_1(12) * occ_func_0_1(33)) +
              (occ_func_0_0(4) * occ_func_0_1(1) * occ_func_0_1(26)) +
              (occ_func_0_0(83) * occ_func_0_1(28) * occ_func_0_1(35)) +
              (occ_func_0_0(11) * occ_func_0_1(10) * occ_func_0_1(30)) +
              (occ_func_0_0(8) * occ_func_0_1(3) * occ_func_0_1(25)) +
              (occ_func_0_0(84) * occ_func_0_1(31) * occ_func_0_1(36)) +
              (occ_func_0_0(1) * occ_func_0_1(4) * occ_func_0_1(40)) +
              (occ_func_0_0(3) * occ_func_0_1(9) * occ_func_0_1(42)) +
              (occ_func_0_0(79) * occ_func_0_1(21) * occ_func_0_1(19)) +
              (occ_func_0_0(5) * occ_func_0_1(2) * occ_func_0_1(27)) +
              (occ_func_0_0(10) * occ_func_0_1(11) * occ_func_0_1(33)) +
              (occ_func_0_0(83) * occ_func_0_1(34) * occ_func_0_1(28)) +
              (occ_func_0_0(2) * occ_func_0_1(7) * occ_func_0_1(41)) +
              (occ_func_0_0(1) * occ_func_0_1(6) * occ_func_0_1(40)) +
              (occ_func_0_0(79) * occ_func_0_1(20) * occ_func_0_1(21)) +
              (occ_func_0_0(6) * occ_func_0_1(1) * occ_func_0_1(23)) +
              (occ_func_0_0(11) * occ_func_0_1(12) * occ_func_0_1(30)) +
              (occ_func_0_0(84) * occ_func_0_1(38) * occ_func_0_1(31)) +
              (occ_func_0_0(7) * occ_func_0_1(9) * occ_func_0_1(38)) +
              (occ_func_0_0(5) * occ_func_0_1(4) * occ_func_0_1(36)) +
              (occ_func_0_0(81) * occ_func_0_1(23) * occ_func_0_1(25)) +
              (occ_func_0_0(5) * occ_func_0_1(10) * occ_func_0_1(36)) +
              (occ_func_0_0(2) * occ_func_0_1(3) * occ_func_0_1(31)) +
              (occ_func_0_0(81) * occ_func_0_1(25) * occ_func_0_1(30)) +
              (occ_func_0_0(1) * occ_func_0_1(3) * occ_func_0_1(32)) +
              (occ_func_0_0(4) * occ_func_0_1(10) * occ_func_0_1(37)) +
              (occ_func_0_0(80) * occ_func_0_1(29) * occ_func_0_1(24)) +
              (occ_func_0_0(8) * occ_func_0_1(6) * occ_func_0_1(34)) +
              (occ_func_0_0(9) * occ_func_0_1(7) * occ_func_0_1(35)) +
              (occ_func_0_0(82) * occ_func_0_1(27) * occ_func_0_1(26)) +
              (occ_func_0_0(6) * occ_func_0_1(11) * occ_func_0_1(39)) +
              (occ_func_0_0(1) * occ_func_0_1(2) * occ_func_0_1(32)) +
              (occ_func_0_0(80) * occ_func_0_1(22) * occ_func_0_1(29)) +
              (occ_func_0_0(12) * occ_func_0_1(7) * occ_func_0_1(21)) +
              (occ_func_0_0(11) * occ_func_0_1(6) * occ_func_0_1(20)) +
              (occ_func_0_0(86) * occ_func_0_1(40) * occ_func_0_1(41)) +
              (occ_func_0_0(10) * occ_func_0_1(4) * occ_func_0_1(19)) +
              (occ_func_0_0(12) * occ_func_0_1(9) * occ_func_0_1(21)) +
              (occ_func_0_0(86) * occ_func_0_1(42) * occ_func_0_1(40)) +
              (occ_func_0_0(3) * occ_func_0_1(2) * occ_func_0_1(28)) +
              (occ_func_0_0(8) * occ_func_0_1(11) * occ_func_0_1(34)) +
              (occ_func_0_0(82) * occ_func_0_1(33) * occ_func_0_1(27)) +
              (occ_func_0_0(3) * occ_func_0_1(8) * occ_func_0_1(42)) +
              (occ_func_0_0(2) * occ_func_0_1(5) * occ_func_0_1(41)) +
              (occ_func_0_0(79) * occ_func_0_1(19) * occ_func_0_1(20)) +
              (occ_func_0_0(7) * occ_func_0_1(5) * occ_func_0_1(22)) +
              (occ_func_0_0(9) * occ_func_0_1(8) * occ_func_0_1(24)) +
              (occ_func_0_0(85) * occ_func_0_1(39) * occ_func_0_1(37)) +
              (occ_func_0_0(4) * occ_func_0_1(6) * occ_func_0_1(26)) +
              (occ_func_0_0(5) * occ_func_0_1(7) * occ_func_0_1(27)) +
              (occ_func_0_0(83) * occ_func_0_1(35) * occ_func_0_1(34)) +
              (occ_func_0_0(11) * occ_func_0_1(8) * occ_func_0_1(20)) +
              (occ_func_0_0(10) * occ_func_0_1(5) * occ_func_0_1(19)) +
              (occ_func_0_0(86) * occ_func_0_1(41) * occ_func_0_1(42)) +
              (occ_func_0_0(2) * occ_func_0_1(1) * occ_func_0_1(31)) +
              (occ_func_0_0(7) * occ_func_0_1(12) * occ_func_0_1(38)) +
              (occ_func_0_0(81) * occ_func_0_1(30) * occ_func_0_1(23)) +
              (occ_func_0_0(8) * occ_func_0_1(9) * occ_func_0_1(25)) +
              (occ_func_0_0(6) * occ_func_0_1(4) * occ_func_0_1(23)) +
              (occ_func_0_0(84) * occ_func_0_1(36) * occ_func_0_1(38)) +
              (occ_func_0_0(12) * occ_func_0_1(11) * occ_func_0_1(29)) +
              (occ_func_0_0(7) * occ_func_0_1(2) * occ_func_0_1(22)) +
              (occ_func_0_0(85) * occ_func_0_1(32) * occ_func_0_1(39)) +
              (occ_func_0_0(9) * occ_func_0_1(12) * occ_func_0_1(35)) +
              (occ_func_0_0(3) * occ_func_0_1(1) * occ_func_0_1(28)) +
              (occ_func_0_0(82) * occ_func_0_1(26) * occ_func_0_1(33)) +
              (occ_func_0_0(9) * occ_func_0_1(3) * occ_func_0_1(24)) +
              (occ_func_0_0(12) * occ_func_0_1(10) * occ_func_0_1(29)) +
              (occ_func_0_0(85) * occ_func_0_1(37) * occ_func_0_1(32))) /
             24.0;
}
double test_Clexulator::delta_site_eval_at_0_bfunc_4_0_11(int occ_i,
                                                          int occ_f) const {
  return (m_occ_func_0_1[occ_f] - m_occ_func_0_1[occ_i]) *
         ((occ_func_0_1(9) * occ_func_0_1(7) * occ_func_0_1(85)) +
          (occ_func_0_1(4) * occ_func_0_1(5) * occ_func_0_1(37)) +
          (occ_func_0_1(6) * occ_func_0_1(8) * occ_func_0_1(39)) +
          (occ_func_0_1(80) * occ_func_0_1(24) * occ_func_0_1(22)) +
          (occ_func_0_1(3) * occ_func_0_1(9) * occ_func_0_1(82)) +
          (occ_func_0_1(10) * occ_func_0_1(12) * occ_func_0_1(33)) +
          (occ_func_0_1(4) * occ_func_0_1(1) * occ_func_0_1(26)) +
          (occ_func_0_1(83) * occ_func_0_1(28) * occ_func_0_1(35)) +
          (occ_func_0_1(2) * occ_func_0_1(5) * occ_func_0_1(81)) +
          (occ_func_0_1(11) * occ_func_0_1(10) * occ_func_0_1(30)) +
          (occ_func_0_1(8) * occ_func_0_1(3) * occ_func_0_1(25)) +
          (occ_func_0_1(84) * occ_func_0_1(31) * occ_func_0_1(36)) +
          (occ_func_0_1(12) * occ_func_0_1(10) * occ_func_0_1(86)) +
          (occ_func_0_1(1) * occ_func_0_1(4) * occ_func_0_1(40)) +
          (occ_func_0_1(3) * occ_func_0_1(9) * occ_func_0_1(42)) +
          (occ_func_0_1(79) * occ_func_0_1(21) * occ_func_0_1(19)) +
          (occ_func_0_1(8) * occ_func_0_1(3) * occ_func_0_1(82)) +
          (occ_func_0_1(5) * occ_func_0_1(2) * occ_func_0_1(27)) +
          (occ_func_0_1(10) * occ_func_0_1(11) * occ_func_0_1(33)) +
          (occ_func_0_1(83) * occ_func_0_1(34) * occ_func_0_1(28)) +
          (occ_func_0_1(11) * occ_func_0_1(12) * occ_func_0_1(86)) +
          (occ_func_0_1(2) * occ_func_0_1(7) * occ_func_0_1(41)) +
          (occ_func_0_1(1) * occ_func_0_1(6) * occ_func_0_1(40)) +
          (occ_func_0_1(79) * occ_func_0_1(20) * occ_func_0_1(21)) +
          (occ_func_0_1(7) * occ_func_0_1(2) * occ_func_0_1(81)) +
          (occ_func_0_1(6) * occ_func_0_1(1) * occ_func_0_1(23)) +
          (occ_func_0_1(11) * occ_func_0_1(12) * occ_func_0_1(30)) +
          (occ_func_0_1(84) * occ_func_0_1(38) * occ_func_0_1(31)) +
          (occ_func_0_1(6) * occ_func_0_1(8) * occ_func_0_1(84)) +
          (occ_func_0_1(7) * occ_func_0_1(9) * occ_func_0_1(38)) +
          (occ_func_0_1(5) * occ_func_0_1(4) * occ_func_0_1(36)) +
          (occ_func_0_1(81) * occ_func_0_1(23) * occ_func_0_1(25)) +
          (occ_func_0_1(8) * occ_func_0_1(11) * occ_func_0_1(84)) +
          (occ_func_0_1(5) * occ_func_0_1(10) * occ_func_0_1(36)) +
          (occ_func_0_1(2) * occ_func_0_1(3) * occ_func_0_1(31)) +
          (occ_func_0_1(81) * occ_func_0_1(25) * occ_func_0_1(30)) +
          (occ_func_0_1(12) * occ_func_0_1(9) * occ_func_0_1(85)) +
          (occ_func_0_1(1) * occ_func_0_1(3) * occ_func_0_1(32)) +
          (occ_func_0_1(4) * occ_func_0_1(10) * occ_func_0_1(37)) +
          (occ_func_0_1(80) * occ_func_0_1(29) * occ_func_0_1(24)) +
          (occ_func_0_1(5) * occ_func_0_1(4) * occ_func_0_1(83)) +
          (occ_func_0_1(8) * occ_func_0_1(6) * occ_func_0_1(34)) +
          (occ_func_0_1(9) * occ_func_0_1(7) * occ_func_0_1(35)) +
          (occ_func_0_1(82) * occ_func_0_1(27) * occ_func_0_1(26)) +
          (occ_func_0_1(7) * occ_func_0_1(12) * occ_func_0_1(85)) +
          (occ_func_0_1(6) * occ_func_0_1(11) * occ_func_0_1(39)) +
          (occ_func_0_1(1) * occ_func_0_1(2) * occ_func_0_1(32)) +
          (occ_func_0_1(80) * occ_func_0_1(22) * occ_func_0_1(29)) +
          (occ_func_0_1(1) * occ_func_0_1(2) * occ_func_0_1(79)) +
          (occ_func_0_1(12) * occ_func_0_1(7) * occ_func_0_1(21)) +
          (occ_func_0_1(11) * occ_func_0_1(6) * occ_func_0_1(20)) +
          (occ_func_0_1(86) * occ_func_0_1(40) * occ_func_0_1(41)) +
          (occ_func_0_1(3) * occ_func_0_1(1) * occ_func_0_1(79)) +
          (occ_func_0_1(10) * occ_func_0_1(4) * occ_func_0_1(19)) +
          (occ_func_0_1(12) * occ_func_0_1(9) * occ_func_0_1(21)) +
          (occ_func_0_1(86) * occ_func_0_1(42) * occ_func_0_1(40)) +
          (occ_func_0_1(10) * occ_func_0_1(5) * occ_func_0_1(83)) +
          (occ_func_0_1(3) * occ_func_0_1(2) * occ_func_0_1(28)) +
          (occ_func_0_1(8) * occ_func_0_1(11) * occ_func_0_1(34)) +
          (occ_func_0_1(82) * occ_func_0_1(33) * occ_func_0_1(27)) +
          (occ_func_0_1(10) * occ_func_0_1(11) * occ_func_0_1(86)) +
          (occ_func_0_1(3) * occ_func_0_1(8) * occ_func_0_1(42)) +
          (occ_func_0_1(2) * occ_func_0_1(5) * occ_func_0_1(41)) +
          (occ_func_0_1(79) * occ_func_0_1(19) * occ_func_0_1(20)) +
          (occ_func_0_1(6) * occ_func_0_1(4) * occ_func_0_1(80)) +
          (occ_func_0_1(7) * occ_func_0_1(5) * occ_func_0_1(22)) +
          (occ_func_0_1(9) * occ_func_0_1(8) * occ_func_0_1(24)) +
          (occ_func_0_1(85) * occ_func_0_1(39) * occ_func_0_1(37)) +
          (occ_func_0_1(9) * occ_func_0_1(8) * occ_func_0_1(82)) +
          (occ_func_0_1(4) * occ_func_0_1(6) * occ_func_0_1(26)) +
          (occ_func_0_1(5) * occ_func_0_1(7) * occ_func_0_1(27)) +
          (occ_func_0_1(83) * occ_func_0_1(35) * occ_func_0_1(34)) +
          (occ_func_0_1(2) * occ_func_0_1(3) * occ_func_0_1(79)) +
          (occ_func_0_1(11) * occ_func_0_1(8) * occ_func_0_1(20)) +
          (occ_func_0_1(10) * occ_func_0_1(5) * occ_func_0_1(19)) +
          (occ_func_0_1(86) * occ_func_0_1(41) * occ_func_0_1(42)) +
          (occ_func_0_1(11) * occ_func_0_1(6) * occ_func_0_1(84)) +
          (occ_func_0_1(2) * occ_func_0_1(1) * occ_func_0_1(31)) +
          (occ_func_0_1(7) * occ_func_0_1(12) * occ_func_0_1(38)) +
          (occ_func_0_1(81) * occ_func_0_1(30) * occ_func_0_1(23)) +
          (occ_func_0_1(5) * occ_func_0_1(7) * occ_func_0_1(81)) +
          (occ_func_0_1(8) * occ_func_0_1(9) * occ_func_0_1(25)) +
          (occ_func_0_1(6) * occ_func_0_1(4) * occ_func_0_1(23)) +
          (occ_func_0_1(84) * occ_func_0_1(36) * occ_func_0_1(38)) +
          (occ_func_0_1(1) * occ_func_0_1(6) * occ_func_0_1(80)) +
          (occ_func_0_1(12) * occ_func_0_1(11) * occ_func_0_1(29)) +
          (occ_func_0_1(7) * occ_func_0_1(2) * occ_func_0_1(22)) +
          (occ_func_0_1(85) * occ_func_0_1(32) * occ_func_0_1(39)) +
          (occ_func_0_1(4) * occ_func_0_1(10) * occ_func_0_1(83)) +
          (occ_func_0_1(9) * occ_func_0_1(12) * occ_func_0_1(35)) +
          (occ_func_0_1(3) * occ_func_0_1(1) * occ_func_0_1(28)) +
          (occ_func_0_1(82) * occ_func_0_1(26) * occ_func_0_1(33)) +
          (occ_func_0_1(4) * occ_func_0_1(1) * occ_func_0_1(80)) +
          (occ_func_0_1(9) * occ_func_0_1(3) * occ_func_0_1(24)) +
          (occ_func_0_1(12) * occ_func_0_1(10) * occ_func_0_1(29)) +
          (occ_func_0_1(85) * occ_func_0_1(37) * occ_func_0_1(32))) /
         24.0;
}

/**** Basis functions for orbit 4, 1****
#Points: 4
MaxLength: 8.4852814  MinLength: 2.8284271
 0.0000000   0.0000000   0.0000000 A B C
 1.0000000   0.0000000   0.0000000 A B C
 2.0000000   0.0000000   0.0000000 A B C
 3.0000000   0.0000000   0.0000000 A B C
****/
double test_Clexulator::eval_bfunc_4_1_0() const {
  return ((occ_func_0_0(0) * occ_func_0_0(12) * occ_func_0_0(54) *
           occ_func_0_0(175)) +
          (occ_func_0_0(0) * occ_func_0_0(8) * occ_func_0_0(50) *
           occ_func_0_0(160)) +
          (occ_func_0_0(0) * occ_func_0_0(7) * occ_func_0_0(49) *
           occ_func_0_0(159)) +
          (occ_func_0_0(0) * occ_func_0_0(11) * occ_func_0_0(53) *
           occ_func_0_0(174)) +
          (occ_func_0_0(0) * occ_func_0_0(9) * occ_func_0_0(51) *
           occ_func_0_0(161)) +
          (occ_func_0_0(0) * occ_func_0_0(10) * occ_func_0_0(52) *
           occ_func_0_0(171))) /
         6.0;
}
double test_Clexulator::eval_bfunc_4_1_1() const {
  return (((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(12) *
                occ_func_0_0(54) * occ_func_0_0(175) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(12) *
                occ_func_0_0(54) * occ_func_0_1(175))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(8) *
                occ_func_0_0(50) * occ_func_0_0(160) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(8) *
                occ_func_0_0(50) * occ_func_0_1(160))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(7) *
                occ_func_0_0(49) * occ_func_0_0(159) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(7) *
                occ_func_0_0(49) * occ_func_0_1(159))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(11) *
                occ_func_0_0(53) * occ_func_0_0(174) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(11) *
                occ_func_0_0(53) * occ_func_0_1(174))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(9) *
                occ_func_0_0(51) * occ_func_0_0(161) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(9) *
                occ_func_0_0(51) * occ_func_0_1(161))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(10) *
                occ_func_0_0(52) * occ_func_0_0(171) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(10) *
                occ_func_0_0(52) * occ_func_0_1(171)))) /
         6.0;
}
double test_Clexulator::eval_bfunc_4_1_2() const {
  return (((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(12) *
                occ_func_0_0(54) * occ_func_0_0(175) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(12) *
                occ_func_0_1(54) * occ_func_0_0(175))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(8) *
                occ_func_0_0(50) * occ_func_0_0(160) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(8) *
                occ_func_0_1(50) * occ_func_0_0(160))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(7) *
                occ_func_0_0(49) * occ_func_0_0(159) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(7) *
                occ_func_0_1(49) * occ_func_0_0(159))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(11) *
                occ_func_0_0(53) * occ_func_0_0(174) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(11) *
                occ_func_0_1(53) * occ_func_0_0(174))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(9) *
                occ_func_0_0(51) * occ_func_0_0(161) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(9) *
                occ_func_0_1(51) * occ_func_0_0(161))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(10) *
                occ_func_0_0(52) * occ_func_0_0(171) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(10) *
                occ_func_0_1(52) * occ_func_0_0(171)))) /
         6.0;
}
double test_Clexulator::eval_bfunc_4_1_3() const {
  return (((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(12) *
                occ_func_0_0(54) * occ_func_0_0(175) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(12) *
                occ_func_0_1(54) * occ_func_0_1(175))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(8) *
                occ_func_0_0(50) * occ_func_0_0(160) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(8) *
                occ_func_0_1(50) * occ_func_0_1(160))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(7) *
                occ_func_0_0(49) * occ_func_0_0(159) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(7) *
                occ_func_0_1(49) * occ_func_0_1(159))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(11) *
                occ_func_0_0(53) * occ_func_0_0(174) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(11) *
                occ_func_0_1(53) * occ_func_0_1(174))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(9) *
                occ_func_0_0(51) * occ_func_0_0(161) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(9) *
                occ_func_0_1(51) * occ_func_0_1(161))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(10) *
                occ_func_0_0(52) * occ_func_0_0(171) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(10) *
                occ_func_0_1(52) * occ_func_0_1(171)))) /
         6.0;
}
double test_Clexulator::eval_bfunc_4_1_4() const {
  return (((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(12) *
                occ_func_0_1(54) * occ_func_0_0(175) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(12) *
                occ_func_0_0(54) * occ_func_0_1(175))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(8) *
                occ_func_0_1(50) * occ_func_0_0(160) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(8) *
                occ_func_0_0(50) * occ_func_0_1(160))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(7) *
                occ_func_0_1(49) * occ_func_0_0(159) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(7) *
                occ_func_0_0(49) * occ_func_0_1(159))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(11) *
                occ_func_0_1(53) * occ_func_0_0(174) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(11) *
                occ_func_0_0(53) * occ_func_0_1(174))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(9) *
                occ_func_0_1(51) * occ_func_0_0(161) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(9) *
                occ_func_0_0(51) * occ_func_0_1(161))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(10) *
                occ_func_0_1(52) * occ_func_0_0(171) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(10) *
                occ_func_0_0(52) * occ_func_0_1(171)))) /
         6.0;
}
double test_Clexulator::eval_bfunc_4_1_5() const {
  return ((occ_func_0_0(0) * occ_func_0_1(12) * occ_func_0_1(54) *
           occ_func_0_0(175)) +
          (occ_func_0_0(0) * occ_func_0_1(8) * occ_func_0_1(50) *
           occ_func_0_0(160)) +
          (occ_func_0_0(0) * occ_func_0_1(7) * occ_func_0_1(49) *
           occ_func_0_0(159)) +
          (occ_func_0_0(0) * occ_func_0_1(11) * occ_func_0_1(53) *
           occ_func_0_0(174)) +
          (occ_func_0_0(0) * occ_func_0_1(9) * occ_func_0_1(51) *
           occ_func_0_0(161)) +
          (occ_func_0_0(0) * occ_func_0_1(10) * occ_func_0_1(52) *
           occ_func_0_0(171))) /
         6.0;
}
double test_Clexulator::eval_bfunc_4_1_6() const {
  return (((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(12) *
                occ_func_0_1(54) * occ_func_0_0(175) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(12) *
                occ_func_0_1(54) * occ_func_0_1(175))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(8) *
                occ_func_0_1(50) * occ_func_0_0(160) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(8) *
                occ_func_0_1(50) * occ_func_0_1(160))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(7) *
                occ_func_0_1(49) * occ_func_0_0(159) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(7) *
                occ_func_0_1(49) * occ_func_0_1(159))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(11) *
                occ_func_0_1(53) * occ_func_0_0(174) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(11) *
                occ_func_0_1(53) * occ_func_0_1(174))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(9) *
                occ_func_0_1(51) * occ_func_0_0(161) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(9) *
                occ_func_0_1(51) * occ_func_0_1(161))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(10) *
                occ_func_0_1(52) * occ_func_0_0(171) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(10) *
                occ_func_0_1(52) * occ_func_0_1(171)))) /
         6.0;
}
double test_Clexulator::eval_bfunc_4_1_7() const {
  return ((occ_func_0_1(0) * occ_func_0_0(12) * occ_func_0_0(54) *
           occ_func_0_1(175)) +
          (occ_func_0_1(0) * occ_func_0_0(8) * occ_func_0_0(50) *
           occ_func_0_1(160)) +
          (occ_func_0_1(0) * occ_func_0_0(7) * occ_func_0_0(49) *
           occ_func_0_1(159)) +
          (occ_func_0_1(0) * occ_func_0_0(11) * occ_func_0_0(53) *
           occ_func_0_1(174)) +
          (occ_func_0_1(0) * occ_func_0_0(9) * occ_func_0_0(51) *
           occ_func_0_1(161)) +
          (occ_func_0_1(0) * occ_func_0_0(10) * occ_func_0_0(52) *
           occ_func_0_1(171))) /
         6.0;
}
double test_Clexulator::eval_bfunc_4_1_8() const {
  return (((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(12) *
                occ_func_0_0(54) * occ_func_0_1(175) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(12) *
                occ_func_0_1(54) * occ_func_0_1(175))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(8) *
                occ_func_0_0(50) * occ_func_0_1(160) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(8) *
                occ_func_0_1(50) * occ_func_0_1(160))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(7) *
                occ_func_0_0(49) * occ_func_0_1(159) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(7) *
                occ_func_0_1(49) * occ_func_0_1(159))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(11) *
                occ_func_0_0(53) * occ_func_0_1(174) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(11) *
                occ_func_0_1(53) * occ_func_0_1(174))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(9) *
                occ_func_0_0(51) * occ_func_0_1(161) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(9) *
                occ_func_0_1(51) * occ_func_0_1(161))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(10) *
                occ_func_0_0(52) * occ_func_0_1(171) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(10) *
                occ_func_0_1(52) * occ_func_0_1(171)))) /
         6.0;
}
double test_Clexulator::eval_bfunc_4_1_9() const {
  return ((occ_func_0_1(0) * occ_func_0_1(12) * occ_func_0_1(54) *
           occ_func_0_1(175)) +
          (occ_func_0_1(0) * occ_func_0_1(8) * occ_func_0_1(50) *
           occ_func_0_1(160)) +
          (occ_func_0_1(0) * occ_func_0_1(7) * occ_func_0_1(49) *
           occ_func_0_1(159)) +
          (occ_func_0_1(0) * occ_func_0_1(11) * occ_func_0_1(53) *
           occ_func_0_1(174)) +
          (occ_func_0_1(0) * occ_func_0_1(9) * occ_func_0_1(51) *
           occ_func_0_1(161)) +
          (occ_func_0_1(0) * occ_func_0_1(10) * occ_func_0_1(52) *
           occ_func_0_1(171))) /
         6.0;
}

double test_Clexulator::site_eval_at_0_bfunc_4_1_0() const {
  return ((occ_func_0_0(0) * occ_func_0_0(12) * occ_func_0_0(54) *
           occ_func_0_0(175)) +
          (occ_func_0_0(1) * occ_func_0_0(0) * occ_func_0_0(12) *
           occ_func_0_0(54)) +
          (occ_func_0_0(43) * occ_func_0_0(1) * occ_func_0_0(0) *
           occ_func_0_0(12)) +
          (occ_func_0_0(142) * occ_func_0_0(43) * occ_func_0_0(1) *
           occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(8) * occ_func_0_0(50) *
           occ_func_0_0(160)) +
          (occ_func_0_0(5) * occ_func_0_0(0) * occ_func_0_0(8) *
           occ_func_0_0(50)) +
          (occ_func_0_0(47) * occ_func_0_0(5) * occ_func_0_0(0) *
           occ_func_0_0(8)) +
          (occ_func_0_0(157) * occ_func_0_0(47) * occ_func_0_0(5) *
           occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(7) * occ_func_0_0(49) *
           occ_func_0_0(159)) +
          (occ_func_0_0(6) * occ_func_0_0(0) * occ_func_0_0(7) *
           occ_func_0_0(49)) +
          (occ_func_0_0(48) * occ_func_0_0(6) * occ_func_0_0(0) *
           occ_func_0_0(7)) +
          (occ_func_0_0(158) * occ_func_0_0(48) * occ_func_0_0(6) *
           occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(11) * occ_func_0_0(53) *
           occ_func_0_0(174)) +
          (occ_func_0_0(2) * occ_func_0_0(0) * occ_func_0_0(11) *
           occ_func_0_0(53)) +
          (occ_func_0_0(44) * occ_func_0_0(2) * occ_func_0_0(0) *
           occ_func_0_0(11)) +
          (occ_func_0_0(143) * occ_func_0_0(44) * occ_func_0_0(2) *
           occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(9) * occ_func_0_0(51) *
           occ_func_0_0(161)) +
          (occ_func_0_0(4) * occ_func_0_0(0) * occ_func_0_0(9) *
           occ_func_0_0(51)) +
          (occ_func_0_0(46) * occ_func_0_0(4) * occ_func_0_0(0) *
           occ_func_0_0(9)) +
          (occ_func_0_0(156) * occ_func_0_0(46) * occ_func_0_0(4) *
           occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_0(10) * occ_func_0_0(52) *
           occ_func_0_0(171)) +
          (occ_func_0_0(3) * occ_func_0_0(0) * occ_func_0_0(10) *
           occ_func_0_0(52)) +
          (occ_func_0_0(45) * occ_func_0_0(3) * occ_func_0_0(0) *
           occ_func_0_0(10)) +
          (occ_func_0_0(146) * occ_func_0_0(45) * occ_func_0_0(3) *
           occ_func_0_0(0))) /
         6.0;
}
double test_Clexulator::site_eval_at_0_bfunc_4_1_1() const {
  return (((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(12) *
                occ_func_0_0(54) * occ_func_0_0(175) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(12) *
                occ_func_0_0(54) * occ_func_0_1(175))) +
          ((0.7071067812 * occ_func_0_1(1) * occ_func_0_0(0) *
                occ_func_0_0(12) * occ_func_0_0(54) +
            0.7071067812 * occ_func_0_0(1) * occ_func_0_0(0) *
                occ_func_0_0(12) * occ_func_0_1(54))) +
          ((0.7071067812 * occ_func_0_1(43) * occ_func_0_0(1) *
                occ_func_0_0(0) * occ_func_0_0(12) +
            0.7071067812 * occ_func_0_0(43) * occ_func_0_0(1) *
                occ_func_0_0(0) * occ_func_0_1(12))) +
          ((0.7071067812 * occ_func_0_1(142) * occ_func_0_0(43) *
                occ_func_0_0(1) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(142) * occ_func_0_0(43) *
                occ_func_0_0(1) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(8) *
                occ_func_0_0(50) * occ_func_0_0(160) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(8) *
                occ_func_0_0(50) * occ_func_0_1(160))) +
          ((0.7071067812 * occ_func_0_1(5) * occ_func_0_0(0) * occ_func_0_0(8) *
                occ_func_0_0(50) +
            0.7071067812 * occ_func_0_0(5) * occ_func_0_0(0) * occ_func_0_0(8) *
                occ_func_0_1(50))) +
          ((0.7071067812 * occ_func_0_1(47) * occ_func_0_0(5) *
                occ_func_0_0(0) * occ_func_0_0(8) +
            0.7071067812 * occ_func_0_0(47) * occ_func_0_0(5) *
                occ_func_0_0(0) * occ_func_0_1(8))) +
          ((0.7071067812 * occ_func_0_1(157) * occ_func_0_0(47) *
                occ_func_0_0(5) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(157) * occ_func_0_0(47) *
                occ_func_0_0(5) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(7) *
                occ_func_0_0(49) * occ_func_0_0(159) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(7) *
                occ_func_0_0(49) * occ_func_0_1(159))) +
          ((0.7071067812 * occ_func_0_1(6) * occ_func_0_0(0) * occ_func_0_0(7) *
                occ_func_0_0(49) +
            0.7071067812 * occ_func_0_0(6) * occ_func_0_0(0) * occ_func_0_0(7) *
                occ_func_0_1(49))) +
          ((0.7071067812 * occ_func_0_1(48) * occ_func_0_0(6) *
                occ_func_0_0(0) * occ_func_0_0(7) +
            0.7071067812 * occ_func_0_0(48) * occ_func_0_0(6) *
                occ_func_0_0(0) * occ_func_0_1(7))) +
          ((0.7071067812 * occ_func_0_1(158) * occ_func_0_0(48) *
                occ_func_0_0(6) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(158) * occ_func_0_0(48) *
                occ_func_0_0(6) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(11) *
                occ_func_0_0(53) * occ_func_0_0(174) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(11) *
                occ_func_0_0(53) * occ_func_0_1(174))) +
          ((0.7071067812 * occ_func_0_1(2) * occ_func_0_0(0) *
                occ_func_0_0(11) * occ_func_0_0(53) +
            0.7071067812 * occ_func_0_0(2) * occ_func_0_0(0) *
                occ_func_0_0(11) * occ_func_0_1(53))) +
          ((0.7071067812 * occ_func_0_1(44) * occ_func_0_0(2) *
                occ_func_0_0(0) * occ_func_0_0(11) +
            0.7071067812 * occ_func_0_0(44) * occ_func_0_0(2) *
                occ_func_0_0(0) * occ_func_0_1(11))) +
          ((0.7071067812 * occ_func_0_1(143) * occ_func_0_0(44) *
                occ_func_0_0(2) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(143) * occ_func_0_0(44) *
                occ_func_0_0(2) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(9) *
                occ_func_0_0(51) * occ_func_0_0(161) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(9) *
                occ_func_0_0(51) * occ_func_0_1(161))) +
          ((0.7071067812 * occ_func_0_1(4) * occ_func_0_0(0) * occ_func_0_0(9) *
                occ_func_0_0(51) +
            0.7071067812 * occ_func_0_0(4) * occ_func_0_0(0) * occ_func_0_0(9) *
                occ_func_0_1(51))) +
          ((0.7071067812 * occ_func_0_1(46) * occ_func_0_0(4) *
                occ_func_0_0(0) * occ_func_0_0(9) +
            0.7071067812 * occ_func_0_0(46) * occ_func_0_0(4) *
                occ_func_0_0(0) * occ_func_0_1(9))) +
          ((0.7071067812 * occ_func_0_1(156) * occ_func_0_0(46) *
                occ_func_0_0(4) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(156) * occ_func_0_0(46) *
                occ_func_0_0(4) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(10) *
                occ_func_0_0(52) * occ_func_0_0(171) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(10) *
                occ_func_0_0(52) * occ_func_0_1(171))) +
          ((0.7071067812 * occ_func_0_1(3) * occ_func_0_0(0) *
                occ_func_0_0(10) * occ_func_0_0(52) +
            0.7071067812 * occ_func_0_0(3) * occ_func_0_0(0) *
                occ_func_0_0(10) * occ_func_0_1(52))) +
          ((0.7071067812 * occ_func_0_1(45) * occ_func_0_0(3) *
                occ_func_0_0(0) * occ_func_0_0(10) +
            0.7071067812 * occ_func_0_0(45) * occ_func_0_0(3) *
                occ_func_0_0(0) * occ_func_0_1(10))) +
          ((0.7071067812 * occ_func_0_1(146) * occ_func_0_0(45) *
                occ_func_0_0(3) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(146) * occ_func_0_0(45) *
                occ_func_0_0(3) * occ_func_0_1(0)))) /
         6.0;
}
double test_Clexulator::site_eval_at_0_bfunc_4_1_2() const {
  return (((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(12) *
                occ_func_0_0(54) * occ_func_0_0(175) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(12) *
                occ_func_0_1(54) * occ_func_0_0(175))) +
          ((0.7071067812 * occ_func_0_0(1) * occ_func_0_1(0) *
                occ_func_0_0(12) * occ_func_0_0(54) +
            0.7071067812 * occ_func_0_0(1) * occ_func_0_0(0) *
                occ_func_0_1(12) * occ_func_0_0(54))) +
          ((0.7071067812 * occ_func_0_0(43) * occ_func_0_1(1) *
                occ_func_0_0(0) * occ_func_0_0(12) +
            0.7071067812 * occ_func_0_0(43) * occ_func_0_0(1) *
                occ_func_0_1(0) * occ_func_0_0(12))) +
          ((0.7071067812 * occ_func_0_0(142) * occ_func_0_1(43) *
                occ_func_0_0(1) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(142) * occ_func_0_0(43) *
                occ_func_0_1(1) * occ_func_0_0(0))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(8) *
                occ_func_0_0(50) * occ_func_0_0(160) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(8) *
                occ_func_0_1(50) * occ_func_0_0(160))) +
          ((0.7071067812 * occ_func_0_0(5) * occ_func_0_1(0) * occ_func_0_0(8) *
                occ_func_0_0(50) +
            0.7071067812 * occ_func_0_0(5) * occ_func_0_0(0) * occ_func_0_1(8) *
                occ_func_0_0(50))) +
          ((0.7071067812 * occ_func_0_0(47) * occ_func_0_1(5) *
                occ_func_0_0(0) * occ_func_0_0(8) +
            0.7071067812 * occ_func_0_0(47) * occ_func_0_0(5) *
                occ_func_0_1(0) * occ_func_0_0(8))) +
          ((0.7071067812 * occ_func_0_0(157) * occ_func_0_1(47) *
                occ_func_0_0(5) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(157) * occ_func_0_0(47) *
                occ_func_0_1(5) * occ_func_0_0(0))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(7) *
                occ_func_0_0(49) * occ_func_0_0(159) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(7) *
                occ_func_0_1(49) * occ_func_0_0(159))) +
          ((0.7071067812 * occ_func_0_0(6) * occ_func_0_1(0) * occ_func_0_0(7) *
                occ_func_0_0(49) +
            0.7071067812 * occ_func_0_0(6) * occ_func_0_0(0) * occ_func_0_1(7) *
                occ_func_0_0(49))) +
          ((0.7071067812 * occ_func_0_0(48) * occ_func_0_1(6) *
                occ_func_0_0(0) * occ_func_0_0(7) +
            0.7071067812 * occ_func_0_0(48) * occ_func_0_0(6) *
                occ_func_0_1(0) * occ_func_0_0(7))) +
          ((0.7071067812 * occ_func_0_0(158) * occ_func_0_1(48) *
                occ_func_0_0(6) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(158) * occ_func_0_0(48) *
                occ_func_0_1(6) * occ_func_0_0(0))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(11) *
                occ_func_0_0(53) * occ_func_0_0(174) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(11) *
                occ_func_0_1(53) * occ_func_0_0(174))) +
          ((0.7071067812 * occ_func_0_0(2) * occ_func_0_1(0) *
                occ_func_0_0(11) * occ_func_0_0(53) +
            0.7071067812 * occ_func_0_0(2) * occ_func_0_0(0) *
                occ_func_0_1(11) * occ_func_0_0(53))) +
          ((0.7071067812 * occ_func_0_0(44) * occ_func_0_1(2) *
                occ_func_0_0(0) * occ_func_0_0(11) +
            0.7071067812 * occ_func_0_0(44) * occ_func_0_0(2) *
                occ_func_0_1(0) * occ_func_0_0(11))) +
          ((0.7071067812 * occ_func_0_0(143) * occ_func_0_1(44) *
                occ_func_0_0(2) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(143) * occ_func_0_0(44) *
                occ_func_0_1(2) * occ_func_0_0(0))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(9) *
                occ_func_0_0(51) * occ_func_0_0(161) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(9) *
                occ_func_0_1(51) * occ_func_0_0(161))) +
          ((0.7071067812 * occ_func_0_0(4) * occ_func_0_1(0) * occ_func_0_0(9) *
                occ_func_0_0(51) +
            0.7071067812 * occ_func_0_0(4) * occ_func_0_0(0) * occ_func_0_1(9) *
                occ_func_0_0(51))) +
          ((0.7071067812 * occ_func_0_0(46) * occ_func_0_1(4) *
                occ_func_0_0(0) * occ_func_0_0(9) +
            0.7071067812 * occ_func_0_0(46) * occ_func_0_0(4) *
                occ_func_0_1(0) * occ_func_0_0(9))) +
          ((0.7071067812 * occ_func_0_0(156) * occ_func_0_1(46) *
                occ_func_0_0(4) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(156) * occ_func_0_0(46) *
                occ_func_0_1(4) * occ_func_0_0(0))) +
          ((0.7071067812 * occ_func_0_0(0) * occ_func_0_1(10) *
                occ_func_0_0(52) * occ_func_0_0(171) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(10) *
                occ_func_0_1(52) * occ_func_0_0(171))) +
          ((0.7071067812 * occ_func_0_0(3) * occ_func_0_1(0) *
                occ_func_0_0(10) * occ_func_0_0(52) +
            0.7071067812 * occ_func_0_0(3) * occ_func_0_0(0) *
                occ_func_0_1(10) * occ_func_0_0(52))) +
          ((0.7071067812 * occ_func_0_0(45) * occ_func_0_1(3) *
                occ_func_0_0(0) * occ_func_0_0(10) +
            0.7071067812 * occ_func_0_0(45) * occ_func_0_0(3) *
                occ_func_0_1(0) * occ_func_0_0(10))) +
          ((0.7071067812 * occ_func_0_0(146) * occ_func_0_1(45) *
                occ_func_0_0(3) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(146) * occ_func_0_0(45) *
                occ_func_0_1(3) * occ_func_0_0(0)))) /
         6.0;
}
double test_Clexulator::site_eval_at_0_bfunc_4_1_3() const {
  return (((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(12) *
                occ_func_0_0(54) * occ_func_0_0(175) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(12) *
                occ_func_0_1(54) * occ_func_0_1(175))) +
          ((0.7071067812 * occ_func_0_1(1) * occ_func_0_1(0) *
                occ_func_0_0(12) * occ_func_0_0(54) +
            0.7071067812 * occ_func_0_0(1) * occ_func_0_0(0) *
                occ_func_0_1(12) * occ_func_0_1(54))) +
          ((0.7071067812 * occ_func_0_1(43) * occ_func_0_1(1) *
                occ_func_0_0(0) * occ_func_0_0(12) +
            0.7071067812 * occ_func_0_0(43) * occ_func_0_0(1) *
                occ_func_0_1(0) * occ_func_0_1(12))) +
          ((0.7071067812 * occ_func_0_1(142) * occ_func_0_1(43) *
                occ_func_0_0(1) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(142) * occ_func_0_0(43) *
                occ_func_0_1(1) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(8) *
                occ_func_0_0(50) * occ_func_0_0(160) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(8) *
                occ_func_0_1(50) * occ_func_0_1(160))) +
          ((0.7071067812 * occ_func_0_1(5) * occ_func_0_1(0) * occ_func_0_0(8) *
                occ_func_0_0(50) +
            0.7071067812 * occ_func_0_0(5) * occ_func_0_0(0) * occ_func_0_1(8) *
                occ_func_0_1(50))) +
          ((0.7071067812 * occ_func_0_1(47) * occ_func_0_1(5) *
                occ_func_0_0(0) * occ_func_0_0(8) +
            0.7071067812 * occ_func_0_0(47) * occ_func_0_0(5) *
                occ_func_0_1(0) * occ_func_0_1(8))) +
          ((0.7071067812 * occ_func_0_1(157) * occ_func_0_1(47) *
                occ_func_0_0(5) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(157) * occ_func_0_0(47) *
                occ_func_0_1(5) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(7) *
                occ_func_0_0(49) * occ_func_0_0(159) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(7) *
                occ_func_0_1(49) * occ_func_0_1(159))) +
          ((0.7071067812 * occ_func_0_1(6) * occ_func_0_1(0) * occ_func_0_0(7) *
                occ_func_0_0(49) +
            0.7071067812 * occ_func_0_0(6) * occ_func_0_0(0) * occ_func_0_1(7) *
                occ_func_0_1(49))) +
          ((0.7071067812 * occ_func_0_1(48) * occ_func_0_1(6) *
                occ_func_0_0(0) * occ_func_0_0(7) +
            0.7071067812 * occ_func_0_0(48) * occ_func_0_0(6) *
                occ_func_0_1(0) * occ_func_0_1(7))) +
          ((0.7071067812 * occ_func_0_1(158) * occ_func_0_1(48) *
                occ_func_0_0(6) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(158) * occ_func_0_0(48) *
                occ_func_0_1(6) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(11) *
                occ_func_0_0(53) * occ_func_0_0(174) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(11) *
                occ_func_0_1(53) * occ_func_0_1(174))) +
          ((0.7071067812 * occ_func_0_1(2) * occ_func_0_1(0) *
                occ_func_0_0(11) * occ_func_0_0(53) +
            0.7071067812 * occ_func_0_0(2) * occ_func_0_0(0) *
                occ_func_0_1(11) * occ_func_0_1(53))) +
          ((0.7071067812 * occ_func_0_1(44) * occ_func_0_1(2) *
                occ_func_0_0(0) * occ_func_0_0(11) +
            0.7071067812 * occ_func_0_0(44) * occ_func_0_0(2) *
                occ_func_0_1(0) * occ_func_0_1(11))) +
          ((0.7071067812 * occ_func_0_1(143) * occ_func_0_1(44) *
                occ_func_0_0(2) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(143) * occ_func_0_0(44) *
                occ_func_0_1(2) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(9) *
                occ_func_0_0(51) * occ_func_0_0(161) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(9) *
                occ_func_0_1(51) * occ_func_0_1(161))) +
          ((0.7071067812 * occ_func_0_1(4) * occ_func_0_1(0) * occ_func_0_0(9) *
                occ_func_0_0(51) +
            0.7071067812 * occ_func_0_0(4) * occ_func_0_0(0) * occ_func_0_1(9) *
                occ_func_0_1(51))) +
          ((0.7071067812 * occ_func_0_1(46) * occ_func_0_1(4) *
                occ_func_0_0(0) * occ_func_0_0(9) +
            0.7071067812 * occ_func_0_0(46) * occ_func_0_0(4) *
                occ_func_0_1(0) * occ_func_0_1(9))) +
          ((0.7071067812 * occ_func_0_1(156) * occ_func_0_1(46) *
                occ_func_0_0(4) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(156) * occ_func_0_0(46) *
                occ_func_0_1(4) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(10) *
                occ_func_0_0(52) * occ_func_0_0(171) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_0(10) *
                occ_func_0_1(52) * occ_func_0_1(171))) +
          ((0.7071067812 * occ_func_0_1(3) * occ_func_0_1(0) *
                occ_func_0_0(10) * occ_func_0_0(52) +
            0.7071067812 * occ_func_0_0(3) * occ_func_0_0(0) *
                occ_func_0_1(10) * occ_func_0_1(52))) +
          ((0.7071067812 * occ_func_0_1(45) * occ_func_0_1(3) *
                occ_func_0_0(0) * occ_func_0_0(10) +
            0.7071067812 * occ_func_0_0(45) * occ_func_0_0(3) *
                occ_func_0_1(0) * occ_func_0_1(10))) +
          ((0.7071067812 * occ_func_0_1(146) * occ_func_0_1(45) *
                occ_func_0_0(3) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(146) * occ_func_0_0(45) *
                occ_func_0_1(3) * occ_func_0_1(0)))) /
         6.0;
}
double test_Clexulator::site_eval_at_0_bfunc_4_1_4() const {
  return (((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(12) *
                occ_func_0_1(54) * occ_func_0_0(175) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(12) *
                occ_func_0_0(54) * occ_func_0_1(175))) +
          ((0.7071067812 * occ_func_0_1(1) * occ_func_0_0(0) *
                occ_func_0_1(12) * occ_func_0_0(54) +
            0.7071067812 * occ_func_0_0(1) * occ_func_0_1(0) *
                occ_func_0_0(12) * occ_func_0_1(54))) +
          ((0.7071067812 * occ_func_0_1(43) * occ_func_0_0(1) *
                occ_func_0_1(0) * occ_func_0_0(12) +
            0.7071067812 * occ_func_0_0(43) * occ_func_0_1(1) *
                occ_func_0_0(0) * occ_func_0_1(12))) +
          ((0.7071067812 * occ_func_0_1(142) * occ_func_0_0(43) *
                occ_func_0_1(1) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(142) * occ_func_0_1(43) *
                occ_func_0_0(1) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(8) *
                occ_func_0_1(50) * occ_func_0_0(160) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(8) *
                occ_func_0_0(50) * occ_func_0_1(160))) +
          ((0.7071067812 * occ_func_0_1(5) * occ_func_0_0(0) * occ_func_0_1(8) *
                occ_func_0_0(50) +
            0.7071067812 * occ_func_0_0(5) * occ_func_0_1(0) * occ_func_0_0(8) *
                occ_func_0_1(50))) +
          ((0.7071067812 * occ_func_0_1(47) * occ_func_0_0(5) *
                occ_func_0_1(0) * occ_func_0_0(8) +
            0.7071067812 * occ_func_0_0(47) * occ_func_0_1(5) *
                occ_func_0_0(0) * occ_func_0_1(8))) +
          ((0.7071067812 * occ_func_0_1(157) * occ_func_0_0(47) *
                occ_func_0_1(5) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(157) * occ_func_0_1(47) *
                occ_func_0_0(5) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(7) *
                occ_func_0_1(49) * occ_func_0_0(159) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(7) *
                occ_func_0_0(49) * occ_func_0_1(159))) +
          ((0.7071067812 * occ_func_0_1(6) * occ_func_0_0(0) * occ_func_0_1(7) *
                occ_func_0_0(49) +
            0.7071067812 * occ_func_0_0(6) * occ_func_0_1(0) * occ_func_0_0(7) *
                occ_func_0_1(49))) +
          ((0.7071067812 * occ_func_0_1(48) * occ_func_0_0(6) *
                occ_func_0_1(0) * occ_func_0_0(7) +
            0.7071067812 * occ_func_0_0(48) * occ_func_0_1(6) *
                occ_func_0_0(0) * occ_func_0_1(7))) +
          ((0.7071067812 * occ_func_0_1(158) * occ_func_0_0(48) *
                occ_func_0_1(6) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(158) * occ_func_0_1(48) *
                occ_func_0_0(6) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(11) *
                occ_func_0_1(53) * occ_func_0_0(174) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(11) *
                occ_func_0_0(53) * occ_func_0_1(174))) +
          ((0.7071067812 * occ_func_0_1(2) * occ_func_0_0(0) *
                occ_func_0_1(11) * occ_func_0_0(53) +
            0.7071067812 * occ_func_0_0(2) * occ_func_0_1(0) *
                occ_func_0_0(11) * occ_func_0_1(53))) +
          ((0.7071067812 * occ_func_0_1(44) * occ_func_0_0(2) *
                occ_func_0_1(0) * occ_func_0_0(11) +
            0.7071067812 * occ_func_0_0(44) * occ_func_0_1(2) *
                occ_func_0_0(0) * occ_func_0_1(11))) +
          ((0.7071067812 * occ_func_0_1(143) * occ_func_0_0(44) *
                occ_func_0_1(2) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(143) * occ_func_0_1(44) *
                occ_func_0_0(2) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(9) *
                occ_func_0_1(51) * occ_func_0_0(161) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(9) *
                occ_func_0_0(51) * occ_func_0_1(161))) +
          ((0.7071067812 * occ_func_0_1(4) * occ_func_0_0(0) * occ_func_0_1(9) *
                occ_func_0_0(51) +
            0.7071067812 * occ_func_0_0(4) * occ_func_0_1(0) * occ_func_0_0(9) *
                occ_func_0_1(51))) +
          ((0.7071067812 * occ_func_0_1(46) * occ_func_0_0(4) *
                occ_func_0_1(0) * occ_func_0_0(9) +
            0.7071067812 * occ_func_0_0(46) * occ_func_0_1(4) *
                occ_func_0_0(0) * occ_func_0_1(9))) +
          ((0.7071067812 * occ_func_0_1(156) * occ_func_0_0(46) *
                occ_func_0_1(4) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(156) * occ_func_0_1(46) *
                occ_func_0_0(4) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_0(10) *
                occ_func_0_1(52) * occ_func_0_0(171) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(10) *
                occ_func_0_0(52) * occ_func_0_1(171))) +
          ((0.7071067812 * occ_func_0_1(3) * occ_func_0_0(0) *
                occ_func_0_1(10) * occ_func_0_0(52) +
            0.7071067812 * occ_func_0_0(3) * occ_func_0_1(0) *
                occ_func_0_0(10) * occ_func_0_1(52))) +
          ((0.7071067812 * occ_func_0_1(45) * occ_func_0_0(3) *
                occ_func_0_1(0) * occ_func_0_0(10) +
            0.7071067812 * occ_func_0_0(45) * occ_func_0_1(3) *
                occ_func_0_0(0) * occ_func_0_1(10))) +
          ((0.7071067812 * occ_func_0_1(146) * occ_func_0_0(45) *
                occ_func_0_1(3) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(146) * occ_func_0_1(45) *
                occ_func_0_0(3) * occ_func_0_1(0)))) /
         6.0;
}
double test_Clexulator::site_eval_at_0_bfunc_4_1_5() const {
  return ((occ_func_0_0(0) * occ_func_0_1(12) * occ_func_0_1(54) *
           occ_func_0_0(175)) +
          (occ_func_0_0(1) * occ_func_0_1(0) * occ_func_0_1(12) *
           occ_func_0_0(54)) +
          (occ_func_0_0(43) * occ_func_0_1(1) * occ_func_0_1(0) *
           occ_func_0_0(12)) +
          (occ_func_0_0(142) * occ_func_0_1(43) * occ_func_0_1(1) *
           occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_1(8) * occ_func_0_1(50) *
           occ_func_0_0(160)) +
          (occ_func_0_0(5) * occ_func_0_1(0) * occ_func_0_1(8) *
           occ_func_0_0(50)) +
          (occ_func_0_0(47) * occ_func_0_1(5) * occ_func_0_1(0) *
           occ_func_0_0(8)) +
          (occ_func_0_0(157) * occ_func_0_1(47) * occ_func_0_1(5) *
           occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_1(7) * occ_func_0_1(49) *
           occ_func_0_0(159)) +
          (occ_func_0_0(6) * occ_func_0_1(0) * occ_func_0_1(7) *
           occ_func_0_0(49)) +
          (occ_func_0_0(48) * occ_func_0_1(6) * occ_func_0_1(0) *
           occ_func_0_0(7)) +
          (occ_func_0_0(158) * occ_func_0_1(48) * occ_func_0_1(6) *
           occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_1(11) * occ_func_0_1(53) *
           occ_func_0_0(174)) +
          (occ_func_0_0(2) * occ_func_0_1(0) * occ_func_0_1(11) *
           occ_func_0_0(53)) +
          (occ_func_0_0(44) * occ_func_0_1(2) * occ_func_0_1(0) *
           occ_func_0_0(11)) +
          (occ_func_0_0(143) * occ_func_0_1(44) * occ_func_0_1(2) *
           occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_1(9) * occ_func_0_1(51) *
           occ_func_0_0(161)) +
          (occ_func_0_0(4) * occ_func_0_1(0) * occ_func_0_1(9) *
           occ_func_0_0(51)) +
          (occ_func_0_0(46) * occ_func_0_1(4) * occ_func_0_1(0) *
           occ_func_0_0(9)) +
          (occ_func_0_0(156) * occ_func_0_1(46) * occ_func_0_1(4) *
           occ_func_0_0(0)) +
          (occ_func_0_0(0) * occ_func_0_1(10) * occ_func_0_1(52) *
           occ_func_0_0(171)) +
          (occ_func_0_0(3) * occ_func_0_1(0) * occ_func_0_1(10) *
           occ_func_0_0(52)) +
          (occ_func_0_0(45) * occ_func_0_1(3) * occ_func_0_1(0) *
           occ_func_0_0(10)) +
          (occ_func_0_0(146) * occ_func_0_1(45) * occ_func_0_1(3) *
           occ_func_0_0(0))) /
         6.0;
}
double test_Clexulator::site_eval_at_0_bfunc_4_1_6() const {
  return (((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(12) *
                occ_func_0_1(54) * occ_func_0_0(175) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(12) *
                occ_func_0_1(54) * occ_func_0_1(175))) +
          ((0.7071067812 * occ_func_0_1(1) * occ_func_0_1(0) *
                occ_func_0_1(12) * occ_func_0_0(54) +
            0.7071067812 * occ_func_0_0(1) * occ_func_0_1(0) *
                occ_func_0_1(12) * occ_func_0_1(54))) +
          ((0.7071067812 * occ_func_0_1(43) * occ_func_0_1(1) *
                occ_func_0_1(0) * occ_func_0_0(12) +
            0.7071067812 * occ_func_0_0(43) * occ_func_0_1(1) *
                occ_func_0_1(0) * occ_func_0_1(12))) +
          ((0.7071067812 * occ_func_0_1(142) * occ_func_0_1(43) *
                occ_func_0_1(1) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(142) * occ_func_0_1(43) *
                occ_func_0_1(1) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(8) *
                occ_func_0_1(50) * occ_func_0_0(160) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(8) *
                occ_func_0_1(50) * occ_func_0_1(160))) +
          ((0.7071067812 * occ_func_0_1(5) * occ_func_0_1(0) * occ_func_0_1(8) *
                occ_func_0_0(50) +
            0.7071067812 * occ_func_0_0(5) * occ_func_0_1(0) * occ_func_0_1(8) *
                occ_func_0_1(50))) +
          ((0.7071067812 * occ_func_0_1(47) * occ_func_0_1(5) *
                occ_func_0_1(0) * occ_func_0_0(8) +
            0.7071067812 * occ_func_0_0(47) * occ_func_0_1(5) *
                occ_func_0_1(0) * occ_func_0_1(8))) +
          ((0.7071067812 * occ_func_0_1(157) * occ_func_0_1(47) *
                occ_func_0_1(5) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(157) * occ_func_0_1(47) *
                occ_func_0_1(5) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(7) *
                occ_func_0_1(49) * occ_func_0_0(159) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(7) *
                occ_func_0_1(49) * occ_func_0_1(159))) +
          ((0.7071067812 * occ_func_0_1(6) * occ_func_0_1(0) * occ_func_0_1(7) *
                occ_func_0_0(49) +
            0.7071067812 * occ_func_0_0(6) * occ_func_0_1(0) * occ_func_0_1(7) *
                occ_func_0_1(49))) +
          ((0.7071067812 * occ_func_0_1(48) * occ_func_0_1(6) *
                occ_func_0_1(0) * occ_func_0_0(7) +
            0.7071067812 * occ_func_0_0(48) * occ_func_0_1(6) *
                occ_func_0_1(0) * occ_func_0_1(7))) +
          ((0.7071067812 * occ_func_0_1(158) * occ_func_0_1(48) *
                occ_func_0_1(6) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(158) * occ_func_0_1(48) *
                occ_func_0_1(6) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(11) *
                occ_func_0_1(53) * occ_func_0_0(174) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(11) *
                occ_func_0_1(53) * occ_func_0_1(174))) +
          ((0.7071067812 * occ_func_0_1(2) * occ_func_0_1(0) *
                occ_func_0_1(11) * occ_func_0_0(53) +
            0.7071067812 * occ_func_0_0(2) * occ_func_0_1(0) *
                occ_func_0_1(11) * occ_func_0_1(53))) +
          ((0.7071067812 * occ_func_0_1(44) * occ_func_0_1(2) *
                occ_func_0_1(0) * occ_func_0_0(11) +
            0.7071067812 * occ_func_0_0(44) * occ_func_0_1(2) *
                occ_func_0_1(0) * occ_func_0_1(11))) +
          ((0.7071067812 * occ_func_0_1(143) * occ_func_0_1(44) *
                occ_func_0_1(2) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(143) * occ_func_0_1(44) *
                occ_func_0_1(2) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(9) *
                occ_func_0_1(51) * occ_func_0_0(161) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(9) *
                occ_func_0_1(51) * occ_func_0_1(161))) +
          ((0.7071067812 * occ_func_0_1(4) * occ_func_0_1(0) * occ_func_0_1(9) *
                occ_func_0_0(51) +
            0.7071067812 * occ_func_0_0(4) * occ_func_0_1(0) * occ_func_0_1(9) *
                occ_func_0_1(51))) +
          ((0.7071067812 * occ_func_0_1(46) * occ_func_0_1(4) *
                occ_func_0_1(0) * occ_func_0_0(9) +
            0.7071067812 * occ_func_0_0(46) * occ_func_0_1(4) *
                occ_func_0_1(0) * occ_func_0_1(9))) +
          ((0.7071067812 * occ_func_0_1(156) * occ_func_0_1(46) *
                occ_func_0_1(4) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(156) * occ_func_0_1(46) *
                occ_func_0_1(4) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(10) *
                occ_func_0_1(52) * occ_func_0_0(171) +
            0.7071067812 * occ_func_0_0(0) * occ_func_0_1(10) *
                occ_func_0_1(52) * occ_func_0_1(171))) +
          ((0.7071067812 * occ_func_0_1(3) * occ_func_0_1(0) *
                occ_func_0_1(10) * occ_func_0_0(52) +
            0.7071067812 * occ_func_0_0(3) * occ_func_0_1(0) *
                occ_func_0_1(10) * occ_func_0_1(52))) +
          ((0.7071067812 * occ_func_0_1(45) * occ_func_0_1(3) *
                occ_func_0_1(0) * occ_func_0_0(10) +
            0.7071067812 * occ_func_0_0(45) * occ_func_0_1(3) *
                occ_func_0_1(0) * occ_func_0_1(10))) +
          ((0.7071067812 * occ_func_0_1(146) * occ_func_0_1(45) *
                occ_func_0_1(3) * occ_func_0_0(0) +
            0.7071067812 * occ_func_0_0(146) * occ_func_0_1(45) *
                occ_func_0_1(3) * occ_func_0_1(0)))) /
         6.0;
}
double test_Clexulator::site_eval_at_0_bfunc_4_1_7() const {
  return ((occ_func_0_1(0) * occ_func_0_0(12) * occ_func_0_0(54) *
           occ_func_0_1(175)) +
          (occ_func_0_1(1) * occ_func_0_0(0) * occ_func_0_0(12) *
           occ_func_0_1(54)) +
          (occ_func_0_1(43) * occ_func_0_0(1) * occ_func_0_0(0) *
           occ_func_0_1(12)) +
          (occ_func_0_1(142) * occ_func_0_0(43) * occ_func_0_0(1) *
           occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_0(8) * occ_func_0_0(50) *
           occ_func_0_1(160)) +
          (occ_func_0_1(5) * occ_func_0_0(0) * occ_func_0_0(8) *
           occ_func_0_1(50)) +
          (occ_func_0_1(47) * occ_func_0_0(5) * occ_func_0_0(0) *
           occ_func_0_1(8)) +
          (occ_func_0_1(157) * occ_func_0_0(47) * occ_func_0_0(5) *
           occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_0(7) * occ_func_0_0(49) *
           occ_func_0_1(159)) +
          (occ_func_0_1(6) * occ_func_0_0(0) * occ_func_0_0(7) *
           occ_func_0_1(49)) +
          (occ_func_0_1(48) * occ_func_0_0(6) * occ_func_0_0(0) *
           occ_func_0_1(7)) +
          (occ_func_0_1(158) * occ_func_0_0(48) * occ_func_0_0(6) *
           occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_0(11) * occ_func_0_0(53) *
           occ_func_0_1(174)) +
          (occ_func_0_1(2) * occ_func_0_0(0) * occ_func_0_0(11) *
           occ_func_0_1(53)) +
          (occ_func_0_1(44) * occ_func_0_0(2) * occ_func_0_0(0) *
           occ_func_0_1(11)) +
          (occ_func_0_1(143) * occ_func_0_0(44) * occ_func_0_0(2) *
           occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_0(9) * occ_func_0_0(51) *
           occ_func_0_1(161)) +
          (occ_func_0_1(4) * occ_func_0_0(0) * occ_func_0_0(9) *
           occ_func_0_1(51)) +
          (occ_func_0_1(46) * occ_func_0_0(4) * occ_func_0_0(0) *
           occ_func_0_1(9)) +
          (occ_func_0_1(156) * occ_func_0_0(46) * occ_func_0_0(4) *
           occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_0(10) * occ_func_0_0(52) *
           occ_func_0_1(171)) +
          (occ_func_0_1(3) * occ_func_0_0(0) * occ_func_0_0(10) *
           occ_func_0_1(52)) +
          (occ_func_0_1(45) * occ_func_0_0(3) * occ_func_0_0(0) *
           occ_func_0_1(10)) +
          (occ_func_0_1(146) * occ_func_0_0(45) * occ_func_0_0(3) *
           occ_func_0_1(0))) /
         6.0;
}
double test_Clexulator::site_eval_at_0_bfunc_4_1_8() const {
  return (((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(12) *
                occ_func_0_0(54) * occ_func_0_1(175) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(12) *
                occ_func_0_1(54) * occ_func_0_1(175))) +
          ((0.7071067812 * occ_func_0_1(1) * occ_func_0_1(0) *
                occ_func_0_0(12) * occ_func_0_1(54) +
            0.7071067812 * occ_func_0_1(1) * occ_func_0_0(0) *
                occ_func_0_1(12) * occ_func_0_1(54))) +
          ((0.7071067812 * occ_func_0_1(43) * occ_func_0_1(1) *
                occ_func_0_0(0) * occ_func_0_1(12) +
            0.7071067812 * occ_func_0_1(43) * occ_func_0_0(1) *
                occ_func_0_1(0) * occ_func_0_1(12))) +
          ((0.7071067812 * occ_func_0_1(142) * occ_func_0_1(43) *
                occ_func_0_0(1) * occ_func_0_1(0) +
            0.7071067812 * occ_func_0_1(142) * occ_func_0_0(43) *
                occ_func_0_1(1) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(8) *
                occ_func_0_0(50) * occ_func_0_1(160) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(8) *
                occ_func_0_1(50) * occ_func_0_1(160))) +
          ((0.7071067812 * occ_func_0_1(5) * occ_func_0_1(0) * occ_func_0_0(8) *
                occ_func_0_1(50) +
            0.7071067812 * occ_func_0_1(5) * occ_func_0_0(0) * occ_func_0_1(8) *
                occ_func_0_1(50))) +
          ((0.7071067812 * occ_func_0_1(47) * occ_func_0_1(5) *
                occ_func_0_0(0) * occ_func_0_1(8) +
            0.7071067812 * occ_func_0_1(47) * occ_func_0_0(5) *
                occ_func_0_1(0) * occ_func_0_1(8))) +
          ((0.7071067812 * occ_func_0_1(157) * occ_func_0_1(47) *
                occ_func_0_0(5) * occ_func_0_1(0) +
            0.7071067812 * occ_func_0_1(157) * occ_func_0_0(47) *
                occ_func_0_1(5) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(7) *
                occ_func_0_0(49) * occ_func_0_1(159) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(7) *
                occ_func_0_1(49) * occ_func_0_1(159))) +
          ((0.7071067812 * occ_func_0_1(6) * occ_func_0_1(0) * occ_func_0_0(7) *
                occ_func_0_1(49) +
            0.7071067812 * occ_func_0_1(6) * occ_func_0_0(0) * occ_func_0_1(7) *
                occ_func_0_1(49))) +
          ((0.7071067812 * occ_func_0_1(48) * occ_func_0_1(6) *
                occ_func_0_0(0) * occ_func_0_1(7) +
            0.7071067812 * occ_func_0_1(48) * occ_func_0_0(6) *
                occ_func_0_1(0) * occ_func_0_1(7))) +
          ((0.7071067812 * occ_func_0_1(158) * occ_func_0_1(48) *
                occ_func_0_0(6) * occ_func_0_1(0) +
            0.7071067812 * occ_func_0_1(158) * occ_func_0_0(48) *
                occ_func_0_1(6) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(11) *
                occ_func_0_0(53) * occ_func_0_1(174) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(11) *
                occ_func_0_1(53) * occ_func_0_1(174))) +
          ((0.7071067812 * occ_func_0_1(2) * occ_func_0_1(0) *
                occ_func_0_0(11) * occ_func_0_1(53) +
            0.7071067812 * occ_func_0_1(2) * occ_func_0_0(0) *
                occ_func_0_1(11) * occ_func_0_1(53))) +
          ((0.7071067812 * occ_func_0_1(44) * occ_func_0_1(2) *
                occ_func_0_0(0) * occ_func_0_1(11) +
            0.7071067812 * occ_func_0_1(44) * occ_func_0_0(2) *
                occ_func_0_1(0) * occ_func_0_1(11))) +
          ((0.7071067812 * occ_func_0_1(143) * occ_func_0_1(44) *
                occ_func_0_0(2) * occ_func_0_1(0) +
            0.7071067812 * occ_func_0_1(143) * occ_func_0_0(44) *
                occ_func_0_1(2) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(9) *
                occ_func_0_0(51) * occ_func_0_1(161) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(9) *
                occ_func_0_1(51) * occ_func_0_1(161))) +
          ((0.7071067812 * occ_func_0_1(4) * occ_func_0_1(0) * occ_func_0_0(9) *
                occ_func_0_1(51) +
            0.7071067812 * occ_func_0_1(4) * occ_func_0_0(0) * occ_func_0_1(9) *
                occ_func_0_1(51))) +
          ((0.7071067812 * occ_func_0_1(46) * occ_func_0_1(4) *
                occ_func_0_0(0) * occ_func_0_1(9) +
            0.7071067812 * occ_func_0_1(46) * occ_func_0_0(4) *
                occ_func_0_1(0) * occ_func_0_1(9))) +
          ((0.7071067812 * occ_func_0_1(156) * occ_func_0_1(46) *
                occ_func_0_0(4) * occ_func_0_1(0) +
            0.7071067812 * occ_func_0_1(156) * occ_func_0_0(46) *
                occ_func_0_1(4) * occ_func_0_1(0))) +
          ((0.7071067812 * occ_func_0_1(0) * occ_func_0_1(10) *
                occ_func_0_0(52) * occ_func_0_1(171) +
            0.7071067812 * occ_func_0_1(0) * occ_func_0_0(10) *
                occ_func_0_1(52) * occ_func_0_1(171))) +
          ((0.7071067812 * occ_func_0_1(3) * occ_func_0_1(0) *
                occ_func_0_0(10) * occ_func_0_1(52) +
            0.7071067812 * occ_func_0_1(3) * occ_func_0_0(0) *
                occ_func_0_1(10) * occ_func_0_1(52))) +
          ((0.7071067812 * occ_func_0_1(45) * occ_func_0_1(3) *
                occ_func_0_0(0) * occ_func_0_1(10) +
            0.7071067812 * occ_func_0_1(45) * occ_func_0_0(3) *
                occ_func_0_1(0) * occ_func_0_1(10))) +
          ((0.7071067812 * occ_func_0_1(146) * occ_func_0_1(45) *
                occ_func_0_0(3) * occ_func_0_1(0) +
            0.7071067812 * occ_func_0_1(146) * occ_func_0_0(45) *
                occ_func_0_1(3) * occ_func_0_1(0)))) /
         6.0;
}
double test_Clexulator::site_eval_at_0_bfunc_4_1_9() const {
  return ((occ_func_0_1(0) * occ_func_0_1(12) * occ_func_0_1(54) *
           occ_func_0_1(175)) +
          (occ_func_0_1(1) * occ_func_0_1(0) * occ_func_0_1(12) *
           occ_func_0_1(54)) +
          (occ_func_0_1(43) * occ_func_0_1(1) * occ_func_0_1(0) *
           occ_func_0_1(12)) +
          (occ_func_0_1(142) * occ_func_0_1(43) * occ_func_0_1(1) *
           occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(8) * occ_func_0_1(50) *
           occ_func_0_1(160)) +
          (occ_func_0_1(5) * occ_func_0_1(0) * occ_func_0_1(8) *
           occ_func_0_1(50)) +
          (occ_func_0_1(47) * occ_func_0_1(5) * occ_func_0_1(0) *
           occ_func_0_1(8)) +
          (occ_func_0_1(157) * occ_func_0_1(47) * occ_func_0_1(5) *
           occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(7) * occ_func_0_1(49) *
           occ_func_0_1(159)) +
          (occ_func_0_1(6) * occ_func_0_1(0) * occ_func_0_1(7) *
           occ_func_0_1(49)) +
          (occ_func_0_1(48) * occ_func_0_1(6) * occ_func_0_1(0) *
           occ_func_0_1(7)) +
          (occ_func_0_1(158) * occ_func_0_1(48) * occ_func_0_1(6) *
           occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(11) * occ_func_0_1(53) *
           occ_func_0_1(174)) +
          (occ_func_0_1(2) * occ_func_0_1(0) * occ_func_0_1(11) *
           occ_func_0_1(53)) +
          (occ_func_0_1(44) * occ_func_0_1(2) * occ_func_0_1(0) *
           occ_func_0_1(11)) +
          (occ_func_0_1(143) * occ_func_0_1(44) * occ_func_0_1(2) *
           occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(9) * occ_func_0_1(51) *
           occ_func_0_1(161)) +
          (occ_func_0_1(4) * occ_func_0_1(0) * occ_func_0_1(9) *
           occ_func_0_1(51)) +
          (occ_func_0_1(46) * occ_func_0_1(4) * occ_func_0_1(0) *
           occ_func_0_1(9)) +
          (occ_func_0_1(156) * occ_func_0_1(46) * occ_func_0_1(4) *
           occ_func_0_1(0)) +
          (occ_func_0_1(0) * occ_func_0_1(10) * occ_func_0_1(52) *
           occ_func_0_1(171)) +
          (occ_func_0_1(3) * occ_func_0_1(0) * occ_func_0_1(10) *
           occ_func_0_1(52)) +
          (occ_func_0_1(45) * occ_func_0_1(3) * occ_func_0_1(0) *
           occ_func_0_1(10)) +
          (occ_func_0_1(146) * occ_func_0_1(45) * occ_func_0_1(3) *
           occ_func_0_1(0))) /
         6.0;
}

double test_Clexulator::delta_site_eval_at_0_bfunc_4_1_0(int occ_i,
                                                         int occ_f) const {
  return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) *
         ((occ_func_0_0(12) * occ_func_0_0(54) * occ_func_0_0(175)) +
          (occ_func_0_0(1) * occ_func_0_0(12) * occ_func_0_0(54)) +
          (occ_func_0_0(43) * occ_func_0_0(1) * occ_func_0_0(12)) +
          (occ_func_0_0(142) * occ_func_0_0(43) * occ_func_0_0(1)) +
          (occ_func_0_0(8) * occ_func_0_0(50) * occ_func_0_0(160)) +
          (occ_func_0_0(5) * occ_func_0_0(8) * occ_func_0_0(50)) +
          (occ_func_0_0(47) * occ_func_0_0(5) * occ_func_0_0(8)) +
          (occ_func_0_0(157) * occ_func_0_0(47) * occ_func_0_0(5)) +
          (occ_func_0_0(7) * occ_func_0_0(49) * occ_func_0_0(159)) +
          (occ_func_0_0(6) * occ_func_0_0(7) * occ_func_0_0(49)) +
          (occ_func_0_0(48) * occ_func_0_0(6) * occ_func_0_0(7)) +
          (occ_func_0_0(158) * occ_func_0_0(48) * occ_func_0_0(6)) +
          (occ_func_0_0(11) * occ_func_0_0(53) * occ_func_0_0(174)) +
          (occ_func_0_0(2) * occ_func_0_0(11) * occ_func_0_0(53)) +
          (occ_func_0_0(44) * occ_func_0_0(2) * occ_func_0_0(11)) +
          (occ_func_0_0(143) * occ_func_0_0(44) * occ_func_0_0(2)) +
          (occ_func_0_0(9) * occ_func_0_0(51) * occ_func_0_0(161)) +
          (occ_func_0_0(4) * occ_func_0_0(9) * occ_func_0_0(51)) +
          (occ_func_0_0(46) * occ_func_0_0(4) * occ_func_0_0(9)) +
          (occ_func_0_0(156) * occ_func_0_0(46) * occ_func_0_0(4)) +
          (occ_func_0_0(10) * occ_func_0_0(52) * occ_func_0_0(171)) +
          (occ_func_0_0(3) * occ_func_0_0(10) * occ_func_0_0(52)) +
          (occ_func_0_0(45) * occ_func_0_0(3) * occ_func_0_0(10)) +
          (occ_func_0_0(146) * occ_func_0_0(45) * occ_func_0_0(3))) /
         6.0;
}
double test_Clexulator::delta_site_eval_at_0_bfunc_4_1_1(int occ_i,
                                                         int occ_f) const {
  return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) *
             ((0.7071067812 * occ_func_0_0(12) * occ_func_0_0(54) *
               occ_func_0_1(175)) +
              ((0.7071067812 * occ_func_0_1(1) * occ_func_0_0(12) *
                    occ_func_0_0(54) +
                0.7071067812 * occ_func_0_0(1) * occ_func_0_0(12) *
                    occ_func_0_1(54))) +
              ((0.7071067812 * occ_func_0_1(43) * occ_func_0_0(1) *
                    occ_func_0_0(12) +
                0.7071067812 * occ_func_0_0(43) * occ_func_0_0(1) *
                    occ_func_0_1(12))) +
              (0.7071067812 * occ_func_0_1(142) * occ_func_0_0(43) *
               occ_func_0_0(1)) +
              (0.7071067812 * occ_func_0_0(8) * occ_func_0_0(50) *
               occ_func_0_1(160)) +
              ((0.7071067812 * occ_func_0_1(5) * occ_func_0_0(8) *
                    occ_func_0_0(50) +
                0.7071067812 * occ_func_0_0(5) * occ_func_0_0(8) *
                    occ_func_0_1(50))) +
              ((0.7071067812 * occ_func_0_1(47) * occ_func_0_0(5) *
                    occ_func_0_0(8) +
                0.7071067812 * occ_func_0_0(47) * occ_func_0_0(5) *
                    occ_func_0_1(8))) +
              (0.7071067812 * occ_func_0_1(157) * occ_func_0_0(47) *
               occ_func_0_0(5)) +
              (0.7071067812 * occ_func_0_0(7) * occ_func_0_0(49) *
               occ_func_0_1(159)) +
              ((0.7071067812 * occ_func_0_1(6) * occ_func_0_0(7) *
                    occ_func_0_0(49) +
                0.7071067812 * occ_func_0_0(6) * occ_func_0_0(7) *
                    occ_func_0_1(49))) +
              ((0.7071067812 * occ_func_0_1(48) * occ_func_0_0(6) *
                    occ_func_0_0(7) +
                0.7071067812 * occ_func_0_0(48) * occ_func_0_0(6) *
                    occ_func_0_1(7))) +
              (0.7071067812 * occ_func_0_1(158) * occ_func_0_0(48) *
               occ_func_0_0(6)) +
              (0.7071067812 * occ_func_0_0(11) * occ_func_0_0(53) *
               occ_func_0_1(174)) +
              ((0.7071067812 * occ_func_0_1(2) * occ_func_0_0(11) *
                    occ_func_0_0(53) +
                0.7071067812 * occ_func_0_0(2) * occ_func_0_0(11) *
                    occ_func_0_1(53))) +
              ((0.7071067812 * occ_func_0_1(44) * occ_func_0_0(2) *
                    occ_func_0_0(11) +
                0.7071067812 * occ_func_0_0(44) * occ_func_0_0(2) *
                    occ_func_0_1(11))) +
              (0.7071067812 * occ_func_0_1(143) * occ_func_0_0(44) *
               occ_func_0_0(2)) +
              (0.7071067812 * occ_func_0_0(9) * occ_func_0_0(51) *
               occ_func_0_1(161)) +
              ((0.7071067812 * occ_func_0_1(4) * occ_func_0_0(9) *
                    occ_func_0_0(51) +
                0.7071067812 * occ_func_0_0(4) * occ_func_0_0(9) *
                    occ_func_0_1(51))) +
              ((0.7071067812 * occ_func_0_1(46) * occ_func_0_0(4) *
                    occ_func_0_0(9) +
                0.7071067812 * occ_func_0_0(46) * occ_func_0_0(4) *
                    occ_func_0_1(9))) +
              (0.7071067812 * occ_func_0_1(156) * occ_func_0_0(46) *
               occ_func_0_0(4)) +
              (0.7071067812 * occ_func_0_0(10) * occ_func_0_0(52) *
               occ_func_0_1(171)) +
              ((0.7071067812 * occ_func_0_1(3) * occ_func_0_0(10) *
                    occ_func_0_0(52) +
                0.7071067812 * occ_func_0_0(3) * occ_func_0_0(10) *
                    occ_func_0_1(52))) +
              ((0.7071067812 * occ_func_0_1(45) * occ_func_0_0(3) *
                    occ_func_0_0(10) +
                0.7071067812 * occ_func_0_0(45) * occ_func_0_0(3) *
                    occ_func_0_1(10))) +
              (0.7071067812 * occ_func_0_1(146) * occ_func_0_0(45) *
               occ_func_0_0(3))) /
             6.0 +
         (m_occ_func_0_1[occ_f] - m_occ_func_0_1[occ_i]) *
             ((0.7071067812 * occ_func_0_0(12) * occ_func_0_0(54) *
               occ_func_0_0(175)) +
              (0.7071067812 * occ_func_0_0(142) * occ_func_0_0(43) *
               occ_func_0_0(1)) +
              (0.7071067812 * occ_func_0_0(8) * occ_func_0_0(50) *
               occ_func_0_0(160)) +
              (0.7071067812 * occ_func_0_0(157) * occ_func_0_0(47) *
               occ_func_0_0(5)) +
              (0.7071067812 * occ_func_0_0(7) * occ_func_0_0(49) *
               occ_func_0_0(159)) +
              (0.7071067812 * occ_func_0_0(158) * occ_func_0_0(48) *
               occ_func_0_0(6)) +
              (0.7071067812 * occ_func_0_0(11) * occ_func_0_0(53) *
               occ_func_0_0(174)) +
              (0.7071067812 * occ_func_0_0(143) * occ_func_0_0(44) *
               occ_func_0_0(2)) +
              (0.7071067812 * occ_func_0_0(9) * occ_func_0_0(51) *
               occ_func_0_0(161)) +
              (0.7071067812 * occ_func_0_0(156) * occ_func_0_0(46) *
               occ_func_0_0(4)) +
              (0.7071067812 * occ_func_0_0(10) * occ_func_0_0(52) *
               occ_func_0_0(171)) +
              (0.7071067812 * occ_func_0_0(146) * occ_func_0_0(45) *
               occ_func_0_0(3))) /
             6.0;
}
double test_Clexulator::delta_site_eval_at_0_bfunc_4_1_2(int occ_i,
                                                         int occ_f) const {
  return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) *
             (((0.7071067812 * occ_func_0_1(12) * occ_func_0_0(54) *
                    occ_func_0_0(175) +
                0.7071067812 * occ_func_0_0(12) * occ_func_0_1(54) *
                    occ_func_0_0(175))) +
              (0.7071067812 * occ_func_0_0(1) * occ_func_0_1(12) *
               occ_func_0_0(54)) +
              (0.7071067812 * occ_func_0_0(43) * occ_func_0_1(1) *
               occ_func_0_0(12)) +
              ((0.7071067812 * occ_func_0_0(142) * occ_func_0_1(43) *
                    occ_func_0_0(1) +
                0.7071067812 * occ_func_0_0(142) * occ_func_0_0(43) *
                    occ_func_0_1(1))) +
              ((0.7071067812 * occ_func_0_1(8) * occ_func_0_0(50) *
                    occ_func_0_0(160) +
                0.7071067812 * occ_func_0_0(8) * occ_func_0_1(50) *
                    occ_func_0_0(160))) +
              (0.7071067812 * occ_func_0_0(5) * occ_func_0_1(8) *
               occ_func_0_0(50)) +
              (0.7071067812 * occ_func_0_0(47) * occ_func_0_1(5) *
               occ_func_0_0(8)) +
              ((0.7071067812 * occ_func_0_0(157) * occ_func_0_1(47) *
                    occ_func_0_0(5) +
                0.7071067812 * occ_func_0_0(157) * occ_func_0_0(47) *
                    occ_func_0_1(5))) +
              ((0.7071067812 * occ_func_0_1(7) * occ_func_0_0(49) *
                    occ_func_0_0(159) +
                0.7071067812 * occ_func_0_0(7) * occ_func_0_1(49) *
                    occ_func_0_0(159))) +
              (0.7071067812 * occ_func_0_0(6) * occ_func_0_1(7) *
               occ_func_0_0(49)) +
              (0.7071067812 * occ_func_0_0(48) * occ_func_0_1(6) *
               occ_func_0_0(7)) +
              ((0.7071067812 * occ_func_0_0(158) * occ_func_0_1(48) *
                    occ_func_0_0(6) +
                0.7071067812 * occ_func_0_0(158) * occ_func_0_0(48) *
                    occ_func_0_1(6))) +
              ((0.7071067812 * occ_func_0_1(11) * occ_func_0_0(53) *
                    occ_func_0_0(174) +
                0.7071067812 * occ_func_0_0(11) * occ_func_0_1(53) *
                    occ_func_0_0(174))) +
              (0.7071067812 * occ_func_0_0(2) * occ_func_0_1(11) *
               occ_func_0_0(53)) +
              (0.7071067812 * occ_func_0_0(44) * occ_func_0_1(2) *
               occ_func_0_0(11)) +
              ((0.7071067812 * occ_func_0_0(143) * occ_func_0_1(44) *
                    occ_func_0_0(2) +
                0.7071067812 * occ_func_0_0(143) * occ_func_0_0(44) *
                    occ_func_0_1(2))) +
              ((0.7071067812 * occ_func_0_1(9) * occ_func_0_0(51) *
                    occ_func_0_0(161) +
                0.7071067812 * occ_func_0_0(9) * occ_func_0_1(51) *
                    occ_func_0_0(161))) +
              (0.7071067812 * occ_func_0_0(4) * occ_func_0_1(9) *
               occ_func_0_0(51)) +
              (0.7071067812 * occ_func_0_0(46) * occ_func_0_1(4) *
               occ_func_0_0(9)) +
              ((0.7071067812 * occ_func_0_0(156) * occ_func_0_1(46) *
                    occ_func_0_0(4) +
                0.7071067812 * occ_func_0_0(156) * occ_func_0_0(46) *
                    occ_func_0_1(4))) +
              ((0.7071067812 * occ_func_0_1(10) * occ_func_0_0(52) *
                    occ_func_0_0(171) +
                0.7071067812 * occ_func_0_0(10) * occ_func_0_1(52) *
                    occ_func_0_0(171))) +
              (0.7071067812 * occ_func_0_0(3) * occ_func_0_1(10) *
               occ_func_0_0(52)) +
              (0.7071067812 * occ_func_0_0(45) * occ_func_0_1(3) *
               occ_func_0_0(10)) +
              ((0.7071067812 * occ_func_0_0(146) * occ_func_0_1(45) *
                    occ_func_0_0(3) +
                0.7071067812 * occ_func_0_0(146) * occ_func_0_0(45) *
                    occ_func_0_1(3)))) /
             6.0 +
         (m_occ_func_0_1[occ_f] - m_occ_func_0_1[occ_i]) *
             ((0.7071067812 * occ_func_0_0(1) * occ_func_0_0(12) *
               occ_func_0_0(54)) +
              (0.7071067812 * occ_func_0_0(43) * occ_func_0_0(1) *
               occ_func_0_0(12)) +
              (0.7071067812 * occ_func_0_0(5) * occ_func_0_0(8) *
               occ_func_0_0(50)) +
              (0.7071067812 * occ_func_0_0(47) * occ_func_0_0(5) *
               occ_func_0_0(8)) +
              (0.7071067812 * occ_func_0_0(6) * occ_func_0_0(7) *
               occ_func_0_0(49)) +
              (0.7071067812 * occ_func_0_0(48) * occ_func_0_0(6) *
               occ_func_0_0(7)) +
              (0.7071067812 * occ_func_0_0(2) * occ_func_0_0(11) *
               occ_func_0_0(53)) +
              (0.7071067812 * occ_func_0_0(44) * occ_func_0_0(2) *
               occ_func_0_0(11)) +
              (0.7071067812 * occ_func_0_0(4) * occ_func_0_0(9) *
               occ_func_0_0(51)) +
              (0.7071067812 * occ_func_0_0(46) * occ_func_0_0(4) *
               occ_func_0_0(9)) +
              (0.7071067812 * occ_func_0_0(3) * occ_func_0_0(10) *
               occ_func_0_0(52)) +
              (0.7071067812 * occ_func_0_0(45) * occ_func_0_0(3) *
               occ_func_0_0(10))) /
             6.0;
}
double test_Clexulator::delta_site_eval_at_0_bfunc_4_1_3(int occ_i,
                                                         int occ_f) const {
  return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) *
             ((0.7071067812 * occ_func_0_0(12) * occ_func_0_1(54) *
               occ_func_0_1(175)) +
              (0.7071067812 * occ_func_0_0(1) * occ_func_0_1(12) *
               occ_func_0_1(54)) +
              (0.7071067812 * occ_func_0_1(43) * occ_func_0_1(1) *
               occ_func_0_0(12)) +
              (0.7071067812 * occ_func_0_1(142) * occ_func_0_1(43) *
               occ_func_0_0(1)) +
              (0.7071067812 * occ_func_0_0(8) * occ_func_0_1(50) *
               occ_func_0_1(160)) +
              (0.7071067812 * occ_func_0_0(5) * occ_func_0_1(8) *
               occ_func_0_1(50)) +
              (0.7071067812 * occ_func_0_1(47) * occ_func_0_1(5) *
               occ_func_0_0(8)) +
              (0.7071067812 * occ_func_0_1(157) * occ_func_0_1(47) *
               occ_func_0_0(5)) +
              (0.7071067812 * occ_func_0_0(7) * occ_func_0_1(49) *
               occ_func_0_1(159)) +
              (0.7071067812 * occ_func_0_0(6) * occ_func_0_1(7) *
               occ_func_0_1(49)) +
              (0.7071067812 * occ_func_0_1(48) * occ_func_0_1(6) *
               occ_func_0_0(7)) +
              (0.7071067812 * occ_func_0_1(158) * occ_func_0_1(48) *
               occ_func_0_0(6)) +
              (0.7071067812 * occ_func_0_0(11) * occ_func_0_1(53) *
               occ_func_0_1(174)) +
              (0.7071067812 * occ_func_0_0(2) * occ_func_0_1(11) *
               occ_func_0_1(53)) +
              (0.7071067812 * occ_func_0_1(44) * occ_func_0_1(2) *
               occ_func_0_0(11)) +
              (0.7071067812 * occ_func_0_1(143) * occ_func_0_1(44) *
               occ_func_0_0(2)) +
              (0.7071067812 * occ_func_0_0(9) * occ_func_0_1(51) *
               occ_func_0_1(161)) +
              (0.7071067812 * occ_func_0_0(4) * occ_func_0_1(9) *
               occ_func_0_1(51)) +
              (0.7071067812 * occ_func_0_1(46) * occ_func_0_1(4) *
               occ_func_0_0(9)) +
              (0.7071067812 * occ_func_0_1(156) * occ_func_0_1(46) *
               occ_func_0_0(4)) +
              (0.7071067812 * occ_func_0_0(10) * occ_func_0_1(52) *
               occ_func_0_1(171)) +
              (0.7071067812 * occ_func_0_0(3) * occ_func_0_1(10) *
               occ_func_0_1(52)) +
              (0.7071067812 * occ_func_0_1(45) * occ_func_0_1(3) *
               occ_func_0_0(10)) +
              (0.7071067812 * occ_func_0_1(146) * occ_func_0_1(45) *
               occ_func_0_0(3))) /
             6.0 +
         (m_occ_func_0_1[occ_f] - m_occ_func_0_1[occ_i]) *
             ((0.7071067812 * occ_func_0_1(12) * occ_func_0_0(54) *
               occ_func_0_0(175)) +
              (0.7071067812 * occ_func_0_1(1) * occ_func_0_0(12) *
               occ_func_0_0(54)) +
              (0.7071067812 * occ_func_0_0(43) * occ_func_0_0(1) *
               occ_func_0_1(12)) +
              (0.7071067812 * occ_func_0_0(142) * occ_func_0_0(43) *
               occ_func_0_1(1)) +
              (0.7071067812 * occ_func_0_1(8) * occ_func_0_0(50) *
               occ_func_0_0(160)) +
              (0.7071067812 * occ_func_0_1(5) * occ_func_0_0(8) *
               occ_func_0_0(50)) +
              (0.7071067812 * occ_func_0_0(47) * occ_func_0_0(5) *
               occ_func_0_1(8)) +
              (0.7071067812 * occ_func_0_0(157) * occ_func_0_0(47) *
               occ_func_0_1(5)) +
              (0.7071067812 * occ_func_0_1(7) * occ_func_0_0(49) *
               occ_func_0_0(159)) +
              (0.7071067812 * occ_func_0_1(6) * occ_func_0_0(7) *
               occ_func_0_0(49)) +
              (0.7071067812 * occ_func_0_0(48) * occ_func_0_0(6) *
               occ_func_0_1(7)) +
              (0.7071067812 * occ_func_0_0(158) * occ_func_0_0(48) *
               occ_func_0_1(6)) +
              (0.7071067812 * occ_func_0_1(11) * occ_func_0_0(53) *
               occ_func_0_0(174)) +
              (0.7071067812 * occ_func_0_1(2) * occ_func_0_0(11) *
               occ_func_0_0(53)) +
              (0.7071067812 * occ_func_0_0(44) * occ_func_0_0(2) *
               occ_func_0_1(11)) +
              (0.7071067812 * occ_func_0_0(143) * occ_func_0_0(44) *
               occ_func_0_1(2)) +
              (0.7071067812 * occ_func_0_1(9) * occ_func_0_0(51) *
               occ_func_0_0(161)) +
              (0.7071067812 * occ_func_0_1(4) * occ_func_0_0(9) *
               occ_func_0_0(51)) +
              (0.7071067812 * occ_func_0_0(46) * occ_func_0_0(4) *
               occ_func_0_1(9)) +
              (0.7071067812 * occ_func_0_0(156) * occ_func_0_0(46) *
               occ_func_0_1(4)) +
              (0.7071067812 * occ_func_0_1(10) * occ_func_0_0(52) *
               occ_func_0_0(171)) +
              (0.7071067812 * occ_func_0_1(3) * occ_func_0_0(10) *
               occ_func_0_0(52)) +
              (0.7071067812 * occ_func_0_0(45) * occ_func_0_0(3) *
               occ_func_0_1(10)) +
              (0.7071067812 * occ_func_0_0(146) * occ_func_0_0(45) *
               occ_func_0_1(3))) /
             6.0;
}
double test_Clexulator::delta_site_eval_at_0_bfunc_4_1_4(int occ_i,
                                                         int occ_f) const {
  return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) *
             ((0.7071067812 * occ_func_0_1(12) * occ_func_0_0(54) *
               occ_func_0_1(175)) +
              (0.7071067812 * occ_func_0_1(1) * occ_func_0_1(12) *
               occ_func_0_0(54)) +
              (0.7071067812 * occ_func_0_0(43) * occ_func_0_1(1) *
               occ_func_0_1(12)) +
              (0.7071067812 * occ_func_0_1(142) * occ_func_0_0(43) *
               occ_func_0_1(1)) +
              (0.7071067812 * occ_func_0_1(8) * occ_func_0_0(50) *
               occ_func_0_1(160)) +
              (0.7071067812 * occ_func_0_1(5) * occ_func_0_1(8) *
               occ_func_0_0(50)) +
              (0.7071067812 * occ_func_0_0(47) * occ_func_0_1(5) *
               occ_func_0_1(8)) +
              (0.7071067812 * occ_func_0_1(157) * occ_func_0_0(47) *
               occ_func_0_1(5)) +
              (0.7071067812 * occ_func_0_1(7) * occ_func_0_0(49) *
               occ_func_0_1(159)) +
              (0.7071067812 * occ_func_0_1(6) * occ_func_0_1(7) *
               occ_func_0_0(49)) +
              (0.7071067812 * occ_func_0_0(48) * occ_func_0_1(6) *
               occ_func_0_1(7)) +
              (0.7071067812 * occ_func_0_1(158) * occ_func_0_0(48) *
               occ_func_0_1(6)) +
              (0.7071067812 * occ_func_0_1(11) * occ_func_0_0(53) *
               occ_func_0_1(174)) +
              (0.7071067812 * occ_func_0_1(2) * occ_func_0_1(11) *
               occ_func_0_0(53)) +
              (0.7071067812 * occ_func_0_0(44) * occ_func_0_1(2) *
               occ_func_0_1(11)) +
              (0.7071067812 * occ_func_0_1(143) * occ_func_0_0(44) *
               occ_func_0_1(2)) +
              (0.7071067812 * occ_func_0_1(9) * occ_func_0_0(51) *
               occ_func_0_1(161)) +
              (0.7071067812 * occ_func_0_1(4) * occ_func_0_1(9) *
               occ_func_0_0(51)) +
              (0.7071067812 * occ_func_0_0(46) * occ_func_0_1(4) *
               occ_func_0_1(9)) +
              (0.7071067812 * occ_func_0_1(156) * occ_func_0_0(46) *
               occ_func_0_1(4)) +
              (0.7071067812 * occ_func_0_1(10) * occ_func_0_0(52) *
               occ_func_0_1(171)) +
              (0.7071067812 * occ_func_0_1(3) * occ_func_0_1(10) *
               occ_func_0_0(52)) +
              (0.7071067812 * occ_func_0_0(45) * occ_func_0_1(3) *
               occ_func_0_1(10)) +
              (0.7071067812 * occ_func_0_1(146) * occ_func_0_0(45) *
               occ_func_0_1(3))) /
             6.0 +
         (m_occ_func_0_1[occ_f] - m_occ_func_0_1[occ_i]) *
             ((0.7071067812 * occ_func_0_0(12) * occ_func_0_1(54) *
               occ_func_0_0(175)) +
              (0.7071067812 * occ_func_0_0(1) * occ_func_0_0(12) *
               occ_func_0_1(54)) +
              (0.7071067812 * occ_func_0_1(43) * occ_func_0_0(1) *
               occ_func_0_0(12)) +
              (0.7071067812 * occ_func_0_0(142) * occ_func_0_1(43) *
               occ_func_0_0(1)) +
              (0.7071067812 * occ_func_0_0(8) * occ_func_0_1(50) *
               occ_func_0_0(160)) +
              (0.7071067812 * occ_func_0_0(5) * occ_func_0_0(8) *
               occ_func_0_1(50)) +
              (0.7071067812 * occ_func_0_1(47) * occ_func_0_0(5) *
               occ_func_0_0(8)) +
              (0.7071067812 * occ_func_0_0(157) * occ_func_0_1(47) *
               occ_func_0_0(5)) +
              (0.7071067812 * occ_func_0_0(7) * occ_func_0_1(49) *
               occ_func_0_0(159)) +
              (0.7071067812 * occ_func_0_0(6) * occ_func_0_0(7) *
               occ_func_0_1(49)) +
              (0.7071067812 * occ_func_0_1(48) * occ_func_0_0(6) *
               occ_func_0_0(7)) +
              (0.7071067812 * occ_func_0_0(158) * occ_func_0_1(48) *
               occ_func_0_0(6)) +
              (0.7071067812 * occ_func_0_0(11) * occ_func_0_1(53) *
               occ_func_0_0(174)) +
              (0.7071067812 * occ_func_0_0(2) * occ_func_0_0(11) *
               occ_func_0_1(53)) +
              (0.7071067812 * occ_func_0_1(44) * occ_func_0_0(2) *
               occ_func_0_0(11)) +
              (0.7071067812 * occ_func_0_0(143) * occ_func_0_1(44) *
               occ_func_0_0(2)) +
              (0.7071067812 * occ_func_0_0(9) * occ_func_0_1(51) *
               occ_func_0_0(161)) +
              (0.7071067812 * occ_func_0_0(4) * occ_func_0_0(9) *
               occ_func_0_1(51)) +
              (0.7071067812 * occ_func_0_1(46) * occ_func_0_0(4) *
               occ_func_0_0(9)) +
              (0.7071067812 * occ_func_0_0(156) * occ_func_0_1(46) *
               occ_func_0_0(4)) +
              (0.7071067812 * occ_func_0_0(10) * occ_func_0_1(52) *
               occ_func_0_0(171)) +
              (0.7071067812 * occ_func_0_0(3) * occ_func_0_0(10) *
               occ_func_0_1(52)) +
              (0.7071067812 * occ_func_0_1(45) * occ_func_0_0(3) *
               occ_func_0_0(10)) +
              (0.7071067812 * occ_func_0_0(146) * occ_func_0_1(45) *
               occ_func_0_0(3))) /
             6.0;
}
double test_Clexulator::delta_site_eval_at_0_bfunc_4_1_5(int occ_i,
                                                         int occ_f) const {
  return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) *
             ((occ_func_0_1(12) * occ_func_0_1(54) * occ_func_0_0(175)) +
              (occ_func_0_0(142) * occ_func_0_1(43) * occ_func_0_1(1)) +
              (occ_func_0_1(8) * occ_func_0_1(50) * occ_func_0_0(160)) +
              (occ_func_0_0(157) * occ_func_0_1(47) * occ_func_0_1(5)) +
              (occ_func_0_1(7) * occ_func_0_1(49) * occ_func_0_0(159)) +
              (occ_func_0_0(158) * occ_func_0_1(48) * occ_func_0_1(6)) +
              (occ_func_0_1(11) * occ_func_0_1(53) * occ_func_0_0(174)) +
              (occ_func_0_0(143) * occ_func_0_1(44) * occ_func_0_1(2)) +
              (occ_func_0_1(9) * occ_func_0_1(51) * occ_func_0_0(161)) +
              (occ_func_0_0(156) * occ_func_0_1(46) * occ_func_0_1(4)) +
              (occ_func_0_1(10) * occ_func_0_1(52) * occ_func_0_0(171)) +
              (occ_func_0_0(146) * occ_func_0_1(45) * occ_func_0_1(3))) /
             6.0 +
         (m_occ_func_0_1[occ_f] - m_occ_func_0_1[occ_i]) *
             ((occ_func_0_0(1) * occ_func_0_1(12) * occ_func_0_0(54)) +
              (occ_func_0_0(43) * occ_func_0_1(1) * occ_func_0_0(12)) +
              (occ_func_0_0(5) * occ_func_0_1(8) * occ_func_0_0(50)) +
              (occ_func_0_0(47) * occ_func_0_1(5) * occ_func_0_0(8)) +
              (occ_func_0_0(6) * occ_func_0_1(7) * occ_func_0_0(49)) +
              (occ_func_0_0(48) * occ_func_0_1(6) * occ_func_0_0(7)) +
              (occ_func_0_0(2) * occ_func_0_1(11) * occ_func_0_0(53)) +
              (occ_func_0_0(44) * occ_func_0_1(2) * occ_func_0_0(11)) +
              (occ_func_0_0(4) * occ_func_0_1(9) * occ_func_0_0(51)) +
              (occ_func_0_0(46) * occ_func_0_1(4) * occ_func_0_0(9)) +
              (occ_func_0_0(3) * occ_func_0_1(10) * occ_func_0_0(52)) +
              (occ_func_0_0(45) * occ_func_0_1(3) * occ_func_0_0(10))) /
             6.0;
}
double test_Clexulator::delta_site_eval_at_0_bfunc_4_1_6(int occ_i,
                                                         int occ_f) const {
  return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) *
             ((0.7071067812 * occ_func_0_1(12) * occ_func_0_1(54) *
               occ_func_0_1(175)) +
              (0.7071067812 * occ_func_0_1(142) * occ_func_0_1(43) *
               occ_func_0_1(1)) +
              (0.7071067812 * occ_func_0_1(8) * occ_func_0_1(50) *
               occ_func_0_1(160)) +
              (0.7071067812 * occ_func_0_1(157) * occ_func_0_1(47) *
               occ_func_0_1(5)) +
              (0.7071067812 * occ_func_0_1(7) * occ_func_0_1(49) *
               occ_func_0_1(159)) +
              (0.7071067812 * occ_func_0_1(158) * occ_func_0_1(48) *
               occ_func_0_1(6)) +
              (0.7071067812 * occ_func_0_1(11) * occ_func_0_1(53) *
               occ_func_0_1(174)) +
              (0.7071067812 * occ_func_0_1(143) * occ_func_0_1(44) *
               occ_func_0_1(2)) +
              (0.7071067812 * occ_func_0_1(9) * occ_func_0_1(51) *
               occ_func_0_1(161)) +
              (0.7071067812 * occ_func_0_1(156) * occ_func_0_1(46) *
               occ_func_0_1(4)) +
              (0.7071067812 * occ_func_0_1(10) * occ_func_0_1(52) *
               occ_func_0_1(171)) +
              (0.7071067812 * occ_func_0_1(146) * occ_func_0_1(45) *
               occ_func_0_1(3))) /
             6.0 +
         (m_occ_func_0_1[occ_f] - m_occ_func_0_1[occ_i]) *
             ((0.7071067812 * occ_func_0_1(12) * occ_func_0_1(54) *
               occ_func_0_0(175)) +
              ((0.7071067812 * occ_func_0_1(1) * occ_func_0_1(12) *
                    occ_func_0_0(54) +
                0.7071067812 * occ_func_0_0(1) * occ_func_0_1(12) *
                    occ_func_0_1(54))) +
              ((0.7071067812 * occ_func_0_1(43) * occ_func_0_1(1) *
                    occ_func_0_0(12) +
                0.7071067812 * occ_func_0_0(43) * occ_func_0_1(1) *
                    occ_func_0_1(12))) +
              (0.7071067812 * occ_func_0_0(142) * occ_func_0_1(43) *
               occ_func_0_1(1)) +
              (0.7071067812 * occ_func_0_1(8) * occ_func_0_1(50) *
               occ_func_0_0(160)) +
              ((0.7071067812 * occ_func_0_1(5) * occ_func_0_1(8) *
                    occ_func_0_0(50) +
                0.7071067812 * occ_func_0_0(5) * occ_func_0_1(8) *
                    occ_func_0_1(50))) +
              ((0.7071067812 * occ_func_0_1(47) * occ_func_0_1(5) *
                    occ_func_0_0(8) +
                0.7071067812 * occ_func_0_0(47) * occ_func_0_1(5) *
                    occ_func_0_1(8))) +
              (0.7071067812 * occ_func_0_0(157) * occ_func_0_1(47) *
               occ_func_0_1(5)) +
              (0.7071067812 * occ_func_0_1(7) * occ_func_0_1(49) *
               occ_func_0_0(159)) +
              ((0.7071067812 * occ_func_0_1(6) * occ_func_0_1(7) *
                    occ_func_0_0(49) +
                0.7071067812 * occ_func_0_0(6) * occ_func_0_1(7) *
                    occ_func_0_1(49))) +
              ((0.7071067812 * occ_func_0_1(48) * occ_func_0_1(6) *
                    occ_func_0_0(7) +
                0.7071067812 * occ_func_0_0(48) * occ_func_0_1(6) *
                    occ_func_0_1(7))) +
              (0.7071067812 * occ_func_0_0(158) * occ_func_0_1(48) *
               occ_func_0_1(6)) +
              (0.7071067812 * occ_func_0_1(11) * occ_func_0_1(53) *
               occ_func_0_0(174)) +
              ((0.7071067812 * occ_func_0_1(2) * occ_func_0_1(11) *
                    occ_func_0_0(53) +
                0.7071067812 * occ_func_0_0(2) * occ_func_0_1(11) *
                    occ_func_0_1(53))) +
              ((0.7071067812 * occ_func_0_1(44) * occ_func_0_1(2) *
                    occ_func_0_0(11) +
                0.7071067812 * occ_func_0_0(44) * occ_func_0_1(2) *
                    occ_func_0_1(11))) +
              (0.7071067812 * occ_func_0_0(143) * occ_func_0_1(44) *
               occ_func_0_1(2)) +
              (0.7071067812 * occ_func_0_1(9) * occ_func_0_1(51) *
               occ_func_0_0(161)) +
              ((0.7071067812 * occ_func_0_1(4) * occ_func_0_1(9) *
                    occ_func_0_0(51) +
                0.7071067812 * occ_func_0_0(4) * occ_func_0_1(9) *
                    occ_func_0_1(51))) +
              ((0.7071067812 * occ_func_0_1(46) * occ_func_0_1(4) *
                    occ_func_0_0(9) +
                0.7071067812 * occ_func_0_0(46) * occ_func_0_1(4) *
                    occ_func_0_1(9))) +
              (0.7071067812 * occ_func_0_0(156) * occ_func_0_1(46) *
               occ_func_0_1(4)) +
              (0.7071067812 * occ_func_0_1(10) * occ_func_0_1(52) *
               occ_func_0_0(171)) +
              ((0.7071067812 * occ_func_0_1(3) * occ_func_0_1(10) *
                    occ_func_0_0(52) +
                0.7071067812 * occ_func_0_0(3) * occ_func_0_1(10) *
                    occ_func_0_1(52))) +
              ((0.7071067812 * occ_func_0_1(45) * occ_func_0_1(3) *
                    occ_func_0_0(10) +
                0.7071067812 * occ_func_0_0(45) * occ_func_0_1(3) *
                    occ_func_0_1(10))) +
              (0.7071067812 * occ_func_0_0(146) * occ_func_0_1(45) *
               occ_func_0_1(3))) /
             6.0;
}
double test_Clexulator::delta_site_eval_at_0_bfunc_4_1_7(int occ_i,
                                                         int occ_f) const {
  return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) *
             ((occ_func_0_1(1) * occ_func_0_0(12) * occ_func_0_1(54)) +
              (occ_func_0_1(43) * occ_func_0_0(1) * occ_func_0_1(12)) +
              (occ_func_0_1(5) * occ_func_0_0(8) * occ_func_0_1(50)) +
              (occ_func_0_1(47) * occ_func_0_0(5) * occ_func_0_1(8)) +
              (occ_func_0_1(6) * occ_func_0_0(7) * occ_func_0_1(49)) +
              (occ_func_0_1(48) * occ_func_0_0(6) * occ_func_0_1(7)) +
              (occ_func_0_1(2) * occ_func_0_0(11) * occ_func_0_1(53)) +
              (occ_func_0_1(44) * occ_func_0_0(2) * occ_func_0_1(11)) +
              (occ_func_0_1(4) * occ_func_0_0(9) * occ_func_0_1(51)) +
              (occ_func_0_1(46) * occ_func_0_0(4) * occ_func_0_1(9)) +
              (occ_func_0_1(3) * occ_func_0_0(10) * occ_func_0_1(52)) +
              (occ_func_0_1(45) * occ_func_0_0(3) * occ_func_0_1(10))) /
             6.0 +
         (m_occ_func_0_1[occ_f] - m_occ_func_0_1[occ_i]) *
             ((occ_func_0_0(12) * occ_func_0_0(54) * occ_func_0_1(175)) +
              (occ_func_0_1(142) * occ_func_0_0(43) * occ_func_0_0(1)) +
              (occ_func_0_0(8) * occ_func_0_0(50) * occ_func_0_1(160)) +
              (occ_func_0_1(157) * occ_func_0_0(47) * occ_func_0_0(5)) +
              (occ_func_0_0(7) * occ_func_0_0(49) * occ_func_0_1(159)) +
              (occ_func_0_1(158) * occ_func_0_0(48) * occ_func_0_0(6)) +
              (occ_func_0_0(11) * occ_func_0_0(53) * occ_func_0_1(174)) +
              (occ_func_0_1(143) * occ_func_0_0(44) * occ_func_0_0(2)) +
              (occ_func_0_0(9) * occ_func_0_0(51) * occ_func_0_1(161)) +
              (occ_func_0_1(156) * occ_func_0_0(46) * occ_func_0_0(4)) +
              (occ_func_0_0(10) * occ_func_0_0(52) * occ_func_0_1(171)) +
              (occ_func_0_1(146) * occ_func_0_0(45) * occ_func_0_0(3))) /
             6.0;
}
double test_Clexulator::delta_site_eval_at_0_bfunc_4_1_8(int occ_i,
                                                         int occ_f) const {
  return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) *
             ((0.7071067812 * occ_func_0_1(1) * occ_func_0_1(12) *
               occ_func_0_1(54)) +
              (0.7071067812 * occ_func_0_1(43) * occ_func_0_1(1) *
               occ_func_0_1(12)) +
              (0.7071067812 * occ_func_0_1(5) * occ_func_0_1(8) *
               occ_func_0_1(50)) +
              (0.7071067812 * occ_func_0_1(47) * occ_func_0_1(5) *
               occ_func_0_1(8)) +
              (0.7071067812 * occ_func_0_1(6) * occ_func_0_1(7) *
               occ_func_0_1(49)) +
              (0.7071067812 * occ_func_0_1(48) * occ_func_0_1(6) *
               occ_func_0_1(7)) +
              (0.7071067812 * occ_func_0_1(2) * occ_func_0_1(11) *
               occ_func_0_1(53)) +
              (0.7071067812 * occ_func_0_1(44) * occ_func_0_1(2) *
               occ_func_0_1(11)) +
              (0.7071067812 * occ_func_0_1(4) * occ_func_0_1(9) *
               occ_func_0_1(51)) +
              (0.7071067812 * occ_func_0_1(46) * occ_func_0_1(4) *
               occ_func_0_1(9)) +
              (0.7071067812 * occ_func_0_1(3) * occ_func_0_1(10) *
               occ_func_0_1(52)) +
              (0.7071067812 * occ_func_0_1(45) * occ_func_0_1(3) *
               occ_func_0_1(10))) /
             6.0 +
         (m_occ_func_0_1[occ_f] - m_occ_func_0_1[occ_i]) *
             (((0.7071067812 * occ_func_0_1(12) * occ_func_0_0(54) *
                    occ_func_0_1(175) +
                0.7071067812 * occ_func_0_0(12) * occ_func_0_1(54) *
                    occ_func_0_1(175))) +
              (0.7071067812 * occ_func_0_1(1) * occ_func_0_0(12) *
               occ_func_0_1(54)) +
              (0.7071067812 * occ_func_0_1(43) * occ_func_0_0(1) *
               occ_func_0_1(12)) +
              ((0.7071067812 * occ_func_0_1(142) * occ_func_0_1(43) *
                    occ_func_0_0(1) +
                0.7071067812 * occ_func_0_1(142) * occ_func_0_0(43) *
                    occ_func_0_1(1))) +
              ((0.7071067812 * occ_func_0_1(8) * occ_func_0_0(50) *
                    occ_func_0_1(160) +
                0.7071067812 * occ_func_0_0(8) * occ_func_0_1(50) *
                    occ_func_0_1(160))) +
              (0.7071067812 * occ_func_0_1(5) * occ_func_0_0(8) *
               occ_func_0_1(50)) +
              (0.7071067812 * occ_func_0_1(47) * occ_func_0_0(5) *
               occ_func_0_1(8)) +
              ((0.7071067812 * occ_func_0_1(157) * occ_func_0_1(47) *
                    occ_func_0_0(5) +
                0.7071067812 * occ_func_0_1(157) * occ_func_0_0(47) *
                    occ_func_0_1(5))) +
              ((0.7071067812 * occ_func_0_1(7) * occ_func_0_0(49) *
                    occ_func_0_1(159) +
                0.7071067812 * occ_func_0_0(7) * occ_func_0_1(49) *
                    occ_func_0_1(159))) +
              (0.7071067812 * occ_func_0_1(6) * occ_func_0_0(7) *
               occ_func_0_1(49)) +
              (0.7071067812 * occ_func_0_1(48) * occ_func_0_0(6) *
               occ_func_0_1(7)) +
              ((0.7071067812 * occ_func_0_1(158) * occ_func_0_1(48) *
                    occ_func_0_0(6) +
                0.7071067812 * occ_func_0_1(158) * occ_func_0_0(48) *
                    occ_func_0_1(6))) +
              ((0.7071067812 * occ_func_0_1(11) * occ_func_0_0(53) *
                    occ_func_0_1(174) +
                0.7071067812 * occ_func_0_0(11) * occ_func_0_1(53) *
                    occ_func_0_1(174))) +
              (0.7071067812 * occ_func_0_1(2) * occ_func_0_0(11) *
               occ_func_0_1(53)) +
              (0.7071067812 * occ_func_0_1(44) * occ_func_0_0(2) *
               occ_func_0_1(11)) +
              ((0.7071067812 * occ_func_0_1(143) * occ_func_0_1(44) *
                    occ_func_0_0(2) +
                0.7071067812 * occ_func_0_1(143) * occ_func_0_0(44) *
                    occ_func_0_1(2))) +
              ((0.7071067812 * occ_func_0_1(9) * occ_func_0_0(51) *
                    occ_func_0_1(161) +
                0.7071067812 * occ_func_0_0(9) * occ_func_0_1(51) *
                    occ_func_0_1(161))) +
              (0.7071067812 * occ_func_0_1(4) * occ_func_0_0(9) *
               occ_func_0_1(51)) +
              (0.7071067812 * occ_func_0_1(46) * occ_func_0_0(4) *
               occ_func_0_1(9)) +
              ((0.7071067812 * occ_func_0_1(156) * occ_func_0_1(46) *
                    occ_func_0_0(4) +
                0.7071067812 * occ_func_0_1(156) * occ_func_0_0(46) *
                    occ_func_0_1(4))) +
              ((0.7071067812 * occ_func_0_1(10) * occ_func_0_0(52) *
                    occ_func_0_1(171) +
                0.7071067812 * occ_func_0_0(10) * occ_func_0_1(52) *
                    occ_func_0_1(171))) +
              (0.7071067812 * occ_func_0_1(3) * occ_func_0_0(10) *
               occ_func_0_1(52)) +
              (0.7071067812 * occ_func_0_1(45) * occ_func_0_0(3) *
               occ_func_0_1(10)) +
              ((0.7071067812 * occ_func_0_1(146) * occ_func_0_1(45) *
                    occ_func_0_0(3) +
                0.7071067812 * occ_func_0_1(146) * occ_func_0_0(45) *
                    occ_func_0_1(3)))) /
             6.0;
}
double test_Clexulator::delta_site_eval_at_0_bfunc_4_1_9(int occ_i,
                                                         int occ_f) const {
  return (m_occ_func_0_1[occ_f] - m_occ_func_0_1[occ_i]) *
         ((occ_func_0_1(12) * occ_func_0_1(54) * occ_func_0_1(175)) +
          (occ_func_0_1(1) * occ_func_0_1(12) * occ_func_0_1(54)) +
          (occ_func_0_1(43) * occ_func_0_1(1) * occ_func_0_1(12)) +
          (occ_func_0_1(142) * occ_func_0_1(43) * occ_func_0_1(1)) +
          (occ_func_0_1(8) * occ_func_0_1(50) * occ_func_0_1(160)) +
          (occ_func_0_1(5) * occ_func_0_1(8) * occ_func_0_1(50)) +
          (occ_func_0_1(47) * occ_func_0_1(5) * occ_func_0_1(8)) +
          (occ_func_0_1(157) * occ_func_0_1(47) * occ_func_0_1(5)) +
          (occ_func_0_1(7) * occ_func_0_1(49) * occ_func_0_1(159)) +
          (occ_func_0_1(6) * occ_func_0_1(7) * occ_func_0_1(49)) +
          (occ_func_0_1(48) * occ_func_0_1(6) * occ_func_0_1(7)) +
          (occ_func_0_1(158) * occ_func_0_1(48) * occ_func_0_1(6)) +
          (occ_func_0_1(11) * occ_func_0_1(53) * occ_func_0_1(174)) +
          (occ_func_0_1(2) * occ_func_0_1(11) * occ_func_0_1(53)) +
          (occ_func_0_1(44) * occ_func_0_1(2) * occ_func_0_1(11)) +
          (occ_func_0_1(143) * occ_func_0_1(44) * occ_func_0_1(2)) +
          (occ_func_0_1(9) * occ_func_0_1(51) * occ_func_0_1(161)) +
          (occ_func_0_1(4) * occ_func_0_1(9) * occ_func_0_1(51)) +
          (occ_func_0_1(46) * occ_func_0_1(4) * occ_func_0_1(9)) +
          (occ_func_0_1(156) * occ_func_0_1(46) * occ_func_0_1(4)) +
          (occ_func_0_1(10) * occ_func_0_1(52) * occ_func_0_1(171)) +
          (occ_func_0_1(3) * occ_func_0_1(10) * occ_func_0_1(52)) +
          (occ_func_0_1(45) * occ_func_0_1(3) * occ_func_0_1(10)) +
          (occ_func_0_1(146) * occ_func_0_1(45) * occ_func_0_1(3))) /
         6.0;
}

}  // namespace CASM

extern "C" {
/// \brief Returns a clexulator::BaseClexulator* owning a test_Clexulator
CASM::clexulator::BaseClexulator *make_test_Clexulator() {
  return new CASM::test_Clexulator();
}
}
