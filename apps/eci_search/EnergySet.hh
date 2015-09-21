/*
 *  EnergySet.hh
 */

#ifndef EnergySet_HH
#define EnergySet_HH

#include <iostream>
#include <string>
#include "jsonParser.hh"
#include "BP_Vec.hh"
#include "BP_Geo.hh"
#include "BP_Plot.hh"
#include "casm/external/Eigen/Dense"

class Energy {
public:

  double Ef;
  double weight;
  BP::BP_Vec<double> concentration;
  double dist_from_hull;
  bool dist_from_hull_found;
  std::string name;

  Energy(BP::BP_Vec<std::string> s_list);

  Energy(const CASM::jsonParser &json);

  Energy(double i0, double i1, BP::BP_Vec<double> i2, double i3, std::string i4);

  friend std::ostream &operator<<(std::ostream &outstream, const Energy &e) {
    outstream << std::setw(15) << std::setprecision(9) << e.Ef << " " << std::setw(15) << std::setprecision(9) <<  e.weight << " " ;

    for(int i = 0; i < e.concentration.size(); i++)
      outstream << std::setw(15) << std::setprecision(9) << e.concentration[i] << " " ;

    if(e.dist_from_hull_found)
      outstream << std::setw(15) << std::setprecision(9) << e.dist_from_hull << "    "  ;
    else
      outstream << std::setw(15) << std::setprecision(9) << "-" << "    "  ;

    outstream << std::setw(15) << std::setprecision(9) << std::left << e.name << " " << std::right;


    return outstream;
  }

};

CASM::jsonParser &to_json(const Energy &nrg, CASM::jsonParser &json);
void from_json(Energy &nrg, const CASM::jsonParser &json);

class EnergySet : private BP::BP_Vec< Energy> {
private:

  //BP::BP_Vec< Energy> val;
  int Nstruct_on;
  bool Nstruct_set;

  BP::Geo	hull;
  bool hull_found;				// whether or not the hull has been calculated
  BP::BP_Vec<int> hull_indices;		// indices into (*this) of the ground_states
  BP::BP_Vec<int> non_hull_indices;	// indices into (*this) of the non-ground_states

  bool E_vec_ready;
  Eigen::VectorXd E_vec;

  std::string name;

  /// detected input file format: "text" or "json"
  std::string m_format;

public:

  /// Constructors:

  EnergySet();

  EnergySet(BP::BP_Vec<Energy> energy_list);

  EnergySet(std::string energy_filename);

  /// Function declarations:

  std::string format() const;

  void set(BP::BP_Vec<Energy> energy_list);

  static double calc_clex(const Correlation &corr, const ECISet &eci, unsigned long int i);

  void calc_clex(const Correlation &corr, const ECISet &eci);

  void calc_nrg_diff(BP::BP_Vec<int> &index_list, const Correlation &corr, const ECISet &eci);

  bool E_vec_is_ready() const;

  void set_E_vec();

  const Eigen::VectorXd &get_E_vec() const;

  BP::BP_Vec< Energy> get(BP::BP_Vec<int> indices) const;

  int get_Nstruct_on() const;

  void set_Nstruct_on();

  unsigned long int size() const;

  double get_weight(unsigned long int i) const;

  void set_weight(unsigned long int i, double w);

  double get_Ef(unsigned long int i) const;

  double get_concentration(unsigned long int i, unsigned long int j) const;

  unsigned long int get_concentration_size() const;

  double get_dist_from_hull(unsigned long int i) const;

  std::string get_name() const;

  double find_energy_min();

  bool concentration_is_less(unsigned long int i1, unsigned long int i2) const;

  void calc_hull(bool calc_dist_from_hull, double hulltol = 1.0e-14);

  bool is_hull_found() const;

  BP::BP_Vec<Energy> get_hull_points();

  BP::BP_Vec<int> get_hull_indices();

  BP::BP_Vec<Energy> get_non_hull_points();

  BP::BP_Vec<int> get_non_hull_indices();

  void write_hull(std::string s, std::string format);

  void write_hull(std::ostream &sout);

  void write(std::string s, std::string format);

  void write(std::ostream &sout);

  void plot(BP::BP_Plot &plot, std::string label, std::string s, int line_style, double line_width, int point_style, double point_size, bool point_face, int component1, int component2);

  void write_below_hull(std::string s, BP::BP_Vec<Energy> energy_list, std::string format);

  void weight_nrg(double A, double B, double kT);

  void weight_ERef(double A, double B, double kT, double ERef);

  double calc_rms(EnergySet &that, bool include_weights, double dist_from_hull_cutoff, unsigned long int &count);

  CASM::jsonParser &to_json(CASM::jsonParser &json) const;
  void from_json(const CASM::jsonParser &json);

};

CASM::jsonParser &to_json(const EnergySet &nrgset, CASM::jsonParser &json);
void from_json(EnergySet &nrgset, const CASM::jsonParser &json);


#endif // EnergySet_HH
