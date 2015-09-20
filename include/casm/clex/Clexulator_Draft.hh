#ifndef MYCLEX1_CLEXULATOR_HH
#define MYCLEX1_CLEXULATOR_HH


namespace CASM {
  class Configuration;
  // PROJECTNAME: MYCLEX1

  namespace MYCLEX1 {

    /****  PRIM  ***
     cubic ZrH2
     4.821393
       0.000    0.500    0.500
       0.500    0.000    0.500
       0.500    0.500    0.000
       Zr   H
       1   2
       Selective Dynamics
       Direct
       0.0000000000000000  0.0000000000000000  0.0000000000000000 Zr Ni Ti
       0.2500000000000000  0.2500000000000000  0.2500000000000000 H Va
       0.7500000000000000  0.7500000000000000  0.7500000000000000 H Va
    **/

    class Clexulator {
    private:
      // number of sites in the supercell, and number of sites in the neighborhood
      int m_site_N, m_nlist_N;

      // Pointer to memory block that holds current nlist
      const int *m_nlist_ptr;

      // Define const lookup table for discrete variables.
      // they get initialized in constructor
      // occfunc_B_M -> B is basis index and M is basis function index
      const double[3] m_occfunc_0_0, m_occfunc_0_1;
      const double[2] m_occfunc_1_0;
      const double[2] m_occfunc_2_0;

      // Pointer to memory block that hold occupation DoF values
      const int *m_occ_ptr;

      // Pointers to memory blocks that hold occupation DoF values
      const double *m_xdisp_ptr;
      const double *m_ydisp_ptr;
      const double *m_zdisp_ptr;

      // Pointer to memory block that holds cluster variables
      // These should probably get calculated externally in a separate Clexulator
      const double *m_c_var_ptr;

    public:
      int nlist_size() const {
        return m_nlist_N;
      };
      int supercell_size() const {
        return m_site_N;
      };
      void set_supercell_size(int _site_N) {
        m_site_N = _site_N
      };

      // populate internal nlist pointer for
      void set_nlist(const int *_nlist_ptr) const;

      void calc_point_corr(double *corr_out);

      void set_occ(const int *_occ_ptr) const;
      void set_xdisp(const double *_occ_ptr) const;
      void set_ydisp(const double *_occ_ptr) const;
      void set_zdisp(const double *_occ_ptr) const;
      void set_cvar(const int *_occ_ptr) const;



      double eval_0_bfunc_2_4_0_0();
      double eval_0_bfunc_2_4_0_1();
      double eval_0_bfunc_2_4_0_2();
      double eval_0_bfunc_2_4_0_3();

      double eval_0_docc_func_0_0_bfunc_2_4_0_3();
      double eval_0_docc_func_0_1_bfunc_2_4_0_3();

    };

    const double ClexClust::sitefunc_0_0 = { -1.0, 0.0, 1.0};
    const double ClexClust::sitefunc_0_1 = {1.0, 0.0, 1.0};
    const double ClexClust::sitefunc_1_0 = { -1.0, 1.0};
    const double ClexClust::sitefunc_2_0 = { -1.0, 1.0};


    /** Branch 3 **/

    /** 5 of 6 Orbits **  Orbit: 3 0  Points: 3  Mult: 12  MinLength: 4.3301270  MaxLength: 5.0000000 **/


    /**     0 of 12 Equivalent Clusters in Orbit 5
                       7.5000000   7.5000000   7.5000000 H
                       2.5000000   7.5000000   7.5000000 H
                       5.0000000   5.0000000  10.0000000 Zr
    **/

    // ClustClust_W_X_Y -> W is number of sites, X is linear orbit index, Y is equiv index
    class ClexClust_3_5_0 : public ClexClust {
    private:
      int[3] nlist_ind;
      int[3] b_ind;
      double[6] clustvar1;
    public:
      void preprocess(Configuration &config, int site_ind);
      double eval(Configuration &config, int site_ind);
      double dp1_1_eval(Configuration &config, int site_ind);
      double dp1_1_eval(Configuration &config, int site_ind);
      double dp1_1_eval(Configuration &config, int site_ind);

    }
#endif
