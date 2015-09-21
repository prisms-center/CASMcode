/*
 *  ECISet.hh
 */

#ifndef ECISet_HH
#define ECISet_HH

#include "jsonParser.hh"
#include "EnergySet.hh"
#include "BP_Vec.hh"
#include "BP_Comb.hh"
#include "casm/external/Eigen/Dense"

double UK = 1e20;

class ECI {
public:
  int orbit;
  int bfunc;
  int label;  // linear index
  int weight; // weight (0=off, 1=on)
  int mult;   // number of equivalent clusters
  int size;   // number of sites in cluster
  double length;
  double value;
  BP::BP_Vec<int> hierarchy;


  ECI() {}

  ECI(const CASM::jsonParser &json);

  ECI(int _label, int _weight, int _mult, int _size, double _length, BP::BP_Vec<int> _hierarchy):
    orbit(-1), bfunc(-1), label(_label), weight(_weight), mult(_mult),
    size(_size), length(_length), value(0), hierarchy(_hierarchy) {
  }

  friend std::ostream &operator<<(std::ostream &outstream, const ECI &e) {
    outstream << std::setprecision(6) << e.label << " \t" <<  e.weight << " \t" << e.mult << " \t" << e.size << " \t" << e.length << " \t" << e.hierarchy.size();
    for(int i = 0; i < e.hierarchy.size(); i++)
      outstream << " \t" << e.hierarchy[i] ;
    return outstream;
  }

  void clear() {
    orbit = bfunc = -1;
    label = weight = mult = size = length = 0;
    hierarchy.clear();
  }

};

CASM::jsonParser &to_json(const ECI &eci, CASM::jsonParser &json);
void from_json(ECI &eci, const CASM::jsonParser &json);


class ECISetState {
public:

  std::string bit_string;
  int Nclust;
  double cv;
  double rms;

  ECISetState() {
    cv = UK;
    rms = UK;
  }

  ECISetState(const std::string &_bit_string, const int &_Nclust, const double &_cv, const double &_rms):
    bit_string(_bit_string), Nclust(_Nclust), cv(_cv), rms(_rms) {
  }

  bool operator<(const ECISetState &RHS) const {
    return (this->cv < RHS.cv);
  }

  bool operator<=(const ECISetState &RHS) const {
    return (this->cv <= RHS.cv);
  }

  bool operator>(const ECISetState &RHS) const {
    return (this->cv > RHS.cv);
  }

  bool operator>=(const ECISetState &RHS) const {
    return (this->cv >= RHS.cv);
  }

  bool operator==(const ECISetState &RHS) const {
    return (bit_string == RHS.bit_string);
  }

};


class ECISet : private BP::BP_Vec<ECI> {
private:

  const EnergySet *nrg;
  const Correlation *corr;

  mutable int Nclust_on;
  mutable bool Nclust_set;
  BP::BP_Vec<int> fix_index;        // 0,1, or -1 for unfixed
  BP::BP_Vec<int> fix_pos_list;		// list of indices of 'fixed' clusters
  BP::BP_Vec<bool> fix_val_list;		// value 'fixed' clusters are fixed to

  int Nmin;
  int Nmax;

  int Nstruct;
  double cv;
  double rms;

  bool singular;

  /// detected input file format: "text" or "json"
  std::string m_format;

public:

  ECISet();

  ECISet(const std::string &eci_in_filename);

  ECISet(const std::string &eci_in_filename, int _Nmin, int _Nmax);

  std::string format() const;

  void read(const std::string &eci_in_filename);

  void set_Nmin(int _Nmin);

  int get_Nmin() const;

  void set_Nmax(int _Nmax);

  int get_Nmax() const;

  bool toggle_allowed(int i) const;

  int get_Nclust_on() const;

  void set_values(const Eigen::VectorXd &ECI);

  void set_values_and_weights(const Eigen::MatrixXd &ECI);

  void clear_weights();

  void set_clust_off(int i);

  void set_clust_on(int i);

  void toggle_clust(int i);

  ECISet &operator=(const BP::BP_Comb &comb);

  ECISetState get_state() const;

  void set_state(const ECISetState &state);

  std::string get_bit_string() const;

  void set_bit_string(const std::string &s);

  unsigned long int size() const;

  int get_weight(unsigned long int i) const;

  int get_size(unsigned long int i) const;

  double get_value(unsigned long int i) const;

  double get_cv() const;

  double get_rms() const;

  int get_Nstruct() const;

  bool get_singular() const;

  void set_cv(double i1);

  void set_rms(double i1);

  void set_Nstruct(int _Nstruct);

  void set_data(const EnergySet &_nrg, const Correlation &_corr);

  bool fix_ok();

  void fix(int i1, bool i2);

  void unfix(int i1);

  void unfix_all();

  void set_fix();

  void set_Rmax_fix(double Rmax);

  void write_ECI(std::ostream &sout);

  void write_ECIin(std::string filename, std::string format);

  void write_ECIin(BP::BP_Vec<int> &index_list, std::string filename, std::string format);

  void write_ECIout(std::string filename, EnergySet &nrg_set, std::string format);

  void write_ECIout(BP::BP_Vec<int> &index_list, std::string filename, EnergySet &nrg_set, std::string format);

  bool operator==(const ECISet &E) const;

  void randomize(int Nclust, MTRand &mtrand);

  void fit();

  void fit(const Correlation &_corr, const EnergySet &nrg_set, bool &_singular);

  void check_cv(const Correlation &_corr, const EnergySet &_nrg_set, bool &_singular);

  static void *fit_threaded(void *arg);

  void set_correlation_matrix(Eigen::MatrixXd &A, const Correlation &_corr, const EnergySet &nrg_set) const;

  CASM::jsonParser &to_json(CASM::jsonParser &json) const;
  void from_json(const CASM::jsonParser &json);

private:

  bool check_if_singular(const Eigen::VectorXd &S) const;

  BP::BP_Vec<double> set_cv_a(const Eigen::MatrixXd &C, const EnergySet &nrg_set) const;


};

CASM::jsonParser &to_json(const ECISet &eciset, CASM::jsonParser &json);
void from_json(ECISet &eciset, const CASM::jsonParser &json);

#endif // ECISet_HH
