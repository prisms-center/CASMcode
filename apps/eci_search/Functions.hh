/*
 *  Functions.hh
 */

#ifndef Functions_HH
#define Functions_HH

bool TEST = false;
bool NEW = false;
int PTHREADS = -1;
int MTHREADS = -1;
double HULLTOL = 1.0e-14;

#include <string>
#include "casm/external/Eigen/Dense"
#include "old_Array.hh"
#include "BP_GVec.hh"

/// Class declarations
class Correlation;
class Energy;
class EnergySet;
class ECI;
class ECISet;
class GeneticAlgorithm;

/// Function declarations
void calc_eci(std::string energy_filename, std::string eci_in_filename, std::string corr_in_filename, BP::BP_Vec<ECISet> &population, double hulltol = 1.0e1 - 4);
void calc_cs_eci(std::string energy_filename, std::string eci_in_filename, std::string corr_in_filename, const BP::BP_Vec<double> &mu, int alg, double hulltol = 1.0e1 - 4, unsigned long int max_step = std::numeric_limits<unsigned long int>::max());
void calc_all_eci(int N, std::string energy_filename, std::string eci_in_filename, std::string corr_in_filename);
void calc_directmin_eci(int Nrand, int Nmin, int Nmax, std::string energy_filename, std::string eci_in_filename, std::string corr_in_filename, BP::BP_Vec<ECISet> &population);
void calc_dfsmin_eci(int Nrand, int Nstop, int Nmin, int Nmax, std::string energy_filename, std::string eci_in_filename, std::string corr_in_filename, BP::BP_Vec<ECISet> &population);
void calc_ga_eci(int Npopulation, int Nmin, int Nmax, int Nchildren, int Nmutations, std::string energy_filename, std::string eci_in_filename, std::string corr_in_filename);
void calc_ga_dir_eci(int Npopulation, int Nmin, int Nmax, int Nchildren, int Nmutations, std::string energy_filename, std::string eci_in_filename, std::string corr_in_filename);
void calc_ga_dfs_eci(int Npopulation, int Nmin, int Nmax, int Nchildren, int Nmutations, int Nstop, std::string energy_filename, std::string eci_in_filename, std::string corr_in_filename);

ECISet direct_min(int Nmin, int Nmax, const ECISet &eci_min_A, Correlation &corr, EnergySet &nrg_set, bool print_steps);
ECISet dfs_min(long int Nstop, int Nmin, int Nmax, const ECISet &eci_start, Correlation &corr, EnergySet &nrg_set, bool print_steps);

ECISet calc_FPC_eci(double mu, double prec_shrink, double prec_ECI, const ECISet &eci_set, const Correlation &corr, const EnergySet &nrg_set, bool print_steps, unsigned long int max_step);
ECISet calc_BI_eci(double mu, double prec_shrink, double prec_ECI, const ECISet &eci_set, const Correlation &corr, const EnergySet &nrg_set, bool print_steps, unsigned long int max_step);
void FPC(const Eigen::MatrixXd &M1, const Eigen::VectorXd &V1, Eigen::VectorXd &ECI, double mu, double tau, double prec_shrink, double prec_ECI, bool print_steps, unsigned long int max_step);
void BI(const Eigen::MatrixXd &Cn, const Eigen::VectorXd &En, const Eigen::MatrixXd &M1, Eigen::VectorXd &ECI, double mu, double tau, double prec_shrink, double prec_ECI, bool print_steps, unsigned long int max_step);
double shrink(Eigen::VectorXd &ECI, const Eigen::VectorXd &G, double mu, double tau);

bool is_bitstring(string s);
BP::BP_Vec<ECISet> read_bit_strings_file(ECISet eci_in, string bit_strings_filename, std::string format);
double max_cv(BP::BP_GVec<ECISet> &population, BP::BP_GVec_Member<ECISet> *&worst_parent);
double min_cv(BP::BP_GVec<ECISet> &population, BP::BP_GVec_Member<ECISet> *&worst_parent);
int get_best_index(const BP::BP_Vec<ECISet> &population);
int get_best_index(const BP::BP_GVec<ECISet> &population);
BP::BP_GVec_Member<ECISet> *get_best_member(BP::BP_GVec<ECISet> &population);

std::string unique(std::string filename);
std::string json_ext(std::string filename);
std::string rm_json_ext(std::string filename);
std::string get_format_from_ext(std::string filename);

#endif // Functions_HH
