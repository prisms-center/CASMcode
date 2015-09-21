/*
 *  eci_search_classes.hh
 *
 *
 *  Created by Brian Puchala on 4/17/2012.
 *  All rights reserved.
 *
 */

#ifndef eci_search_classes_HH
#define eci_search_classes_HH

bool TEST = false;

#include <sstream>

/// Class declarations
class Correlation;
class Energy;
class EnergySet;
class ECI;
class ECISet;
class GeneticAlgorithm;

/// Function declarations
void calc_eci(string energy_filename, string eci_in_filename, string corr_in_filename);
void calc_cs_eci(string energy_filename, string eci_in_filename, string corr_in_filename, const BP_Vec<double> &mu, int alg);
void calc_all_eci(int N, string energy_filename, string eci_in_filename, string corr_in_filename);
void calc_directmin_eci(int Nrand, int Nmin, int Nmax, string energy_filename, string eci_in_filename, string corr_in_filename, BP_Vec<ECISet> &population);
void calc_dfsmin_eci(int Nrand, int Nstop, int Nmin, int Nmax, string energy_filename, string eci_in_filename, string corr_in_filename, BP_Vec<ECISet> &population);
void calc_ga_eci(int Npopulation, int Nmin, int Nmax, int Nchildren, int Nmutations, string energy_filename, string eci_in_filename, string corr_in_filename);
void calc_ga_dir_eci(int Npopulation, int Nmin, int Nmax, int Nchildren, int Nmutations, string energy_filename, string eci_in_filename, string corr_in_filename);
void calc_ga_dfs_eci(int Npopulation, int Nmin, int Nmax, int Nchildren, int Nmutations, int Nstop, string energy_filename, string eci_in_filename, string corr_in_filename);

ECISet direct_min(int Nmin, int Nmax, const ECISet &eci_min_A, Correlation &corr, EnergySet &nrg_set, bool print_steps);
ECISet dfs_min(long int Nstop, int Nmin, int Nmax, const ECISet &eci_start, Correlation &corr, EnergySet &nrg_set, bool print_steps);

void set_correlation_matrix(QXu::Array &A, Correlation &corr, EnergySet &nrg_set, ECISet &eci_in);
void set_energy_vector(QXu::Array &E_vec, EnergySet &nrg_set);
bool check_if_singular(QXu::Array &W);
void fit_Array(Correlation &corr, EnergySet &nrg_set, ECISet &eci_in, bool &singular);

void set_correlation_matrix(MatrixXd &A, const Correlation &corr, const EnergySet &nrg_set, const ECISet &eci_in);
bool check_if_singular(const VectorXd &S);
void fit(Correlation &corr, EnergySet &nrg_set, ECISet &eci_in, bool &singular);

ECISet calc_FPC_eci(double mu, double prec_shrink, double prec_ECI, const ECISet &eci_set, const Correlation &corr, const EnergySet &nrg_set, bool print_steps);
ECISet calc_BI_eci(double mu, double prec_shrink, double prec_ECI, const ECISet &eci_set, const Correlation &corr, const EnergySet &nrg_set, bool print_steps);
void FPC(const MatrixXd &M1, const VectorXd &V1, VectorXd &ECI, double mu, double tau, double prec_shrink, double prec_ECI, bool print_steps);
void BI(const MatrixXd &Cn, const VectorXd &En, const MatrixXd &M1, VectorXd &ECI, double mu, double tau, double prec_shrink, double prec_ECI, bool print_steps);
double shrink(VectorXd &ECI, const VectorXd &G, double mu, double tau);

bool is_bitstring(string s);
BP_Vec<ECISet> read_bit_strings_file(ECISet eci_in, string bit_strings_filename);
double max_cv(BP_GVec<ECISet> &population, BP_GVec_Member<ECISet> *&worst_parent);
double min_cv(BP_GVec<ECISet> &population, BP_GVec_Member<ECISet> *&worst_parent);
int get_best_index(const BP_Vec<ECISet> &population);
int get_best_index(const BP_GVec<ECISet> &population);
BP_GVec_Member<ECISet> *get_best_member(BP_GVec<ECISet> &population);


class Correlation : public BP_Vec< BP_Vec< double> > {
private:
  //BP_Vec< BP_Vec<double> > val;		//val[ configuration][ cluster] = correlation

  bool cv_a_ready;
  BP_Vec<double> cv_a;		// used to calculate cv score

public:
  Correlation() {
    cv_a_ready = false;
  };

  Correlation(string corr_in_filename) {
    cv_a_ready = false;

    BP_Parse file(corr_in_filename);

    int Nclust;
    int Ncon;
    string s1;
    BP_Vec<double> list;

    s1 = file.getline();
    Nclust = next_int(s1);

    s1 = file.getline();
    Ncon = next_int(s1);

    //cout << "Reading " << corr_in_filename << endl;
    //cout << "  Nclust: " << Nclust << " Ncon: " << Ncon << endl;

    //"clusters" label
    s1 = file.getline();

    // data
    do {
      list = file.getline_double();
      if(list.size() != 0)
        add(list);
    }
    while(file.eof() == false);

    //cout << "val.size: " << val.size() << endl;
    if(size() != Ncon) {
      cout << "Error reading '" << corr_in_filename << "': stated #configurations == " << Ncon << ", but found #configurations == " << size() << endl;
      exit(1);
    }

    for(int i = 0; i < size(); i++) {
      //cout << "val[i].size(): " << val[i].size() << endl;
      if((*this)[i].size() != Nclust) {
        cout << "Error: the eci.in file stated #clusters == " << Nclust << ", but reading '" << corr_in_filename << "' found #clusters == " << (*this)[i].size() << " for configuration " << i << endl;
        exit(1);
      }
    }



  };

  void cluster_subset(BP_Vec<int> &index_list) {
    unsigned long int i, j;

    /*cout << "index_list: " << index_list << endl;

    cout << "corr: " << endl;
    for( i=0; i<size(); i++)
    {
    	for(j=0;j<(*this)[i].size(); j++)
    	{
    		cout << "   " <<  (*this)[i][ j];
    	}
    	cout << endl;
    }
    cout << endl;
    */


    BP_Vec< BP_Vec< double> > tmp = (*this);

    /*cout << "tmp corr: " << endl;
    for( i=0; i<tmp.size(); i++)
    {
    	for(j=0;j<tmp[i].size(); j++)
    	{
    		cout << "   " <<  tmp[i][ j];
    	}
    	cout << endl;
    }
    cout << endl;

    cout << "tmp corr subset: " << endl;
    for( i=0; i<tmp.size(); i++)
    {
    	for(j=0; j<index_list.size(); j++)
    	{
    		cout << "   " <<  tmp[i][ index_list[j] ];
    	}
    	cout << endl;
    }
    cout << endl;
    */

    for(i = 0; i < (*this).size(); i++)
      (*this)[i].clear();

    for(i = 0; i < tmp.size(); i++) {
      for(j = 0; j < index_list.size(); j++) {
        (*this)[i].add(tmp[i][ index_list[j]]);
      }
    }

  };

  void write(string filename) {
    unsigned long int i, j;
    BP_Write file(filename);
    file.newfile();

    file << (*this)[0].size() << " # number of clusters" << endl;
    file << (*this).size() << " # number of configurations" << endl;
    file << "clusters" << endl;
    for(i = 0; i < size(); i++) {
      for(j = 0; j < (*this)[i].size(); j++) {
        file << "   " << (*this)[i][j] ;
      }
      file << "\n";
    }
  };

  bool cv_a_is_ready() {
    return cv_a_ready;
  }

  double get_cv_a(unsigned long int i) {
    return cv_a[i];
  };

  void cv_a_reset() {
    cv_a_ready = false;
  }

  void write_covar(const ECISet &eci_in, const EnergySet &nrg_set, string fname) const;

  void set_cv_a(const MatrixXd &C, const EnergySet &nrg_set);
};

class Energy {
public:

  double Ef;
  double weight;
  BP_Vec<double> concentration;
  double dist_from_hull;
  bool dist_from_hull_found;
  string name;

  Energy(BP_Vec<string> s_list) {
    if(s_list.size() < 5) {
      cout << "Error: Energy constructor: s_list.size() = " << s_list.size() << ", expected size >= 5" << endl;
      exit(1);
    }


    Ef = BP::stod(s_list[0]);
    weight =  BP::stod(s_list[1]);
    for(int i = 2; i < s_list.size() - 2; i++)
      concentration.add(BP::stod(s_list[i]));
    dist_from_hull =  BP::stod(s_list[s_list.size() - 2]);
    dist_from_hull_found = true;
    name = s_list[s_list.size() - 1];

  };

  Energy(double i0, double i1, BP_Vec<double> i2, double i3, string i4) {
    Ef = i0;
    weight = i1;
    concentration = i2;
    dist_from_hull = i3;
    dist_from_hull_found = true;
    name = i4;
  };

  friend ostream &operator<<(ostream &outstream, const Energy &e) {
    outstream << setw(15) << setprecision(9) << e.Ef << " " << setw(15) << setprecision(9) <<  e.weight << " " ;

    for(int i = 0; i < e.concentration.size(); i++)
      outstream << setw(15) << setprecision(9) << e.concentration[i] << " " ;

    if(e.dist_from_hull_found)
      outstream << setw(15) << setprecision(9) << e.dist_from_hull << "    "  ;
    else
      outstream << setw(15) << setprecision(9) << "-" << "    "  ;

    outstream << setw(15) << setprecision(9) << left << e.name << " " << right;


    return outstream;
  };

};

class EnergySet : private BP_Vec< Energy> {
private:

  //BP_Vec< Energy> val;
  int Nstruct_on;
  bool Nstruct_set;

  Geo	hull;
  bool hull_found;				// whether or not the hull has been calculated
  BP_Vec<int> hull_indices;		// indices into (*this) of the ground_states
  BP_Vec<int> non_hull_indices;	// indices into (*this) of the non-ground_states

  bool E_vec_ready;
  VectorXd E_vec;

public:

  /// Constructors:

  EnergySet() {
    E_vec_ready = false;
    hull_found = 0;
    Nstruct_set = 0;
    set_Nstruct_on();

  };

  EnergySet(BP_Vec<Energy> energy_list) {
    E_vec_ready = false;
    hull_found = 0;
    Nstruct_set = 0;

    BP_Vec< Energy>::operator=(energy_list);

    set_Nstruct_on();

  };

  EnergySet(string energy_filename) {
    E_vec_ready = false;
    hull_found = 0;
    Nstruct_set = 0;
    BP_Parse file(energy_filename);

    // header line
    file.getline();

    BP_Vec<string> entry;
    // data
    do {
      entry = file.getline_string();
      if(entry.size() >= 5)
        add(entry);
    }
    while(file.eof() == false);

    set_Nstruct_on();
    //cout << "energy_set.size(): " << val.size() << endl;
    //for( int i=0; i<val.size(); i++)
    //	cout << val[i] << endl;
  };

  /// Functions declarations:

  void set(BP_Vec<Energy> energy_list) {
    E_vec_ready = false;
    hull_found = 0;
    Nstruct_set = 0;

    BP_Vec< Energy>::operator=(energy_list);

    set_Nstruct_on();
  };

  void calc_clex(const Correlation &corr, const ECISet &eci);
  void calc_nrg_diff(BP_Vec<int> &index_list, const Correlation &corr, const ECISet &eci);

  /// Functions:

  bool E_vec_is_ready() {
    return E_vec_ready;
  };

  void set_E_vec() {
    //cout << "begin set_E_vec()" << endl;
    // set E_vec to be the weighted energy vector

    E_vec.resize(get_Nstruct_on());

    //cout << "get_Nstruct_on(): " << get_Nstruct_on() << endl;
    //cout << "E_vec.size(): " << E_vec.size() << endl;

    int i, in_i = 0;
    for(i = 0; i < size(); i++)
      if(get_weight(i) != 0) {
        E_vec(in_i) = get_weight(i) * get_Ef(i);
        in_i++;
      }

    E_vec_ready = true;

    //cout << "finish set_E_vec()" << endl;
  };

  const VectorXd &get_E_vec() {
    return E_vec;
  };

  BP_Vec< Energy> get(BP_Vec<int> indices) const {
    // return a BP_Vec< Energy> for the input indices
    //    if any indices are bad, return an empty BP_Vec

    BP_Vec< Energy> energy_list;
    unsigned long int i;
    bool no_errors = true;
    for(i = 0; i < indices.size(); i++) {
      if(indices[i] >= 0 && indices[i] < size())
        energy_list.add((*this)[ indices[i]]);
      else {
        no_errors = false;
        break;
      }
    }

    if(!no_errors)
      energy_list.clear();

    return energy_list;
  }

  int get_Nstruct_on() const {
    return Nstruct_on;
  };

  void set_Nstruct_on() {
    Nstruct_on = 0;
    for(int i = 0; i < size(); i++) {
      if((*this)[i].weight != 0)
        Nstruct_on++;
    }
    Nstruct_set = 1;
  };

  unsigned long int size() const {
    return BP_Vec<Energy>::size();
  };

  double get_weight(unsigned long int i) const {
    return (*this)[i].weight;
  };
  double get_Ef(unsigned long int i) const {
    return (*this)[i].Ef;
  };
  double get_concentration(unsigned long int i, unsigned long int j) const {
    return (*this)[i].concentration[j];
  };
  unsigned long int get_concentration_size() const {
    return (*this)[0].concentration.size();
  };
  double get_dist_from_hull(unsigned long int i) const {
    return (*this)[i].dist_from_hull;
  };

  double find_energy_min() {
    double minVal = get_Ef(0);
    for(int i = 0; i < BP_Vec<Energy>::size(); i++) {
      if(get_Ef(i) <= minVal) {
        minVal = get_Ef(i);
      }
    }
    //        std:cout<<minVal<<"\n";
    return minVal;
  }

  bool concentration_is_less(unsigned long int i1, unsigned long int i2) const {

    // return true if (*this)[i1].x < (*this)[i2].x
    //   for ties, return false

    unsigned long int N = (*this)[i1].concentration.size();
    unsigned long int k = 0;

    bool cont = true;
    while(cont == true) {
      if((*this)[i1].concentration[k] == (*this)[i2].concentration[k]) {
        if(k == (*this)[i1].concentration.size() - 1) {
          // tie
          return false;
        }

        k++;
      }
      else {
        cont = false;
      }
    }


    if((*this)[i1].concentration[k] < (*this)[i2].concentration[k]) {
      return true;
    }

    return false;

  };

  void calc_hull(bool calc_dist_from_hull) {
    // create matrix of points -------------------------
    unsigned long int i, j;
    unsigned long int Nrows = 1 + (*this)[0].concentration.size();	// dimensions are nrg + conc's
    unsigned long int Ncols = size();
    MatrixXd m(Nrows, Ncols);

    //cout << "Nrows: " << Nrows << endl;
    //cout << "Ncols: " << Ncols << endl;


    for(i = 0; i < size(); i++) {
      m(0, i) = get_Ef(i);

      for(j = 0; j < Nrows - 1; j++) {
        m(j + 1, i) = get_concentration(i, j);
      }


    }

    //cout << "m: " << endl;
    //for( i=0; i<m.cols(); i++)
    //{
    //for( j=0; j<m.rows(); j++)
    //	cout << m(j,i) << " " ;
    //cout << endl;
    //}
    //cout << endl;

    // get hull	----------------------------------------
    hull.reset_points(m);
    hull.set_verbosity(0);
    VectorXd bottom_vector = VectorXd::Zero(Nrows);
    bottom_vector(0) = -1;


    //cout << "EnergySet::calc_hull calculating convex hull..." << endl;
    hull_found = hull.calc_CH();
    //cout << "  EnergySet::calc_hull hull calculated" << endl;

    if(hull_found) {
      //hull.write_equivalent_points(cout);
      hull.CH_bottom(bottom_vector);
      hull_indices = hull.CH_verts_indices();

      // calc dist_to_hull ------------------------------
      if(calc_dist_from_hull) {
        for(i = 0; i < size(); i++) {
          (*this)[i].dist_from_hull = hull.CH_dist_to_hull(i);
          (*this)[i].dist_from_hull_found = true;
        }
      }

      // get non_hull_indices --------------------------

      non_hull_indices.clear();
      for(i = 0; i < size(); i++)
        if(!hull_indices.find_first(i, j)) {
          // check for extra hull points... (those that are co-planar with other hull points)
          if(abs((*this)[i].dist_from_hull) < hull.get_Geo_tol()) {
            cout << "*** eci_search has added co-planar points to the hull. Please ignore any warnings about co-planar points. ***" << endl;
            hull_indices.add(i);
          }
          else {
            non_hull_indices.add(i);
          }

        }


      // sort hull_indices by conc -----------------------

      if(hull_indices.size() > 1)
        for(i = 0; i < hull_indices.size() - 1; i++) {
          for(j = hull_indices.size() - 1; j > i; j--) {
            if(concentration_is_less(hull_indices[j], hull_indices[j - 1]))
              hull_indices.swap(j, j - 1);
          }
        }


      // sort non_hull_indices by conc -----------------------

      if(non_hull_indices.size() > 1)
        for(i = 0; i < non_hull_indices.size() - 1; i++) {
          for(j = non_hull_indices.size() - 1; j > i; j--) {
            if(concentration_is_less(non_hull_indices[j], non_hull_indices[j - 1]))
              non_hull_indices.swap(j, j - 1);
          }
        }



      //cout << "print hull points: " << endl;
      //for( i=0; i<hull_indices.size(); i++)
      //	cout << " point: " << (*this)[hull_indices[i]] << endl;

      //cout << endl;
      //cout << "print non-hull points: " << endl;
      //for( i=0; i<non_hull_indices.size(); i++)
      //	cout << " point: " << (*this)[non_hull_indices[i]] << endl;

      //cout << endl;
      //cout << "print non-hull points on hull: " << endl;
      //for( i=0; i<non_hull_indices.size(); i++)
      //	if( abs( (*this)[non_hull_indices[i]].dist_from_hull) < hull.get_Geo_tol())
      //	{
      //		cout << "warning, non-hull point on hull" << endl;
      //		cout << " point: " << (*this)[i] << endl;
      //	}


    }
    else {
      if(calc_dist_from_hull) {
        for(i = 0; i < size(); i++)
          (*this)[i].dist_from_hull_found = false;
      }
    }
  };

  bool is_hull_found() const {
    return hull_found;
  };

  BP_Vec<Energy> get_hull_points() {
    return get(get_hull_indices());
  };

  BP_Vec<int> get_hull_indices() {
    if(!hull_found) {
      calc_hull(false);
    }

    return hull_indices;
  };

  BP_Vec<Energy> get_non_hull_points() {
    return get(get_non_hull_indices());
  };

  BP_Vec<int> get_non_hull_indices() {
    if(!hull_found) {
      calc_hull(false);
    }

    return non_hull_indices;
  };

  void write_hull(string s) {
    BP_Write file(s);
    file.newfile();
    write_hull(file.get_ostream());

  };

  void write_hull(ostream &sout) {
    if(!hull_found) {
      calc_hull(false);
      if(!hull_found)
        return;
    }

    unsigned long int i, j, k, ii;
    BP_Vec< BP_Vec< int> > equiv_points = hull.get_equivalent_points();

    sout << setw(16) << setprecision(9) << "#E_form ";
    sout << setw(16) << setprecision(9) << "weight ";
    sout << setw(16) << setprecision(9) << "concs ";
    sout << setw(16) << setprecision(9) << "dist_from_hull ";
    sout << setw(16) << setprecision(9) << left << "   name " << right;
    sout << "\n";

    // write hull ---------------------------------------
    for(i = 0; i < hull_indices.size(); i++) {
      sout << (*this)[ hull_indices[i]] ;

      for(j = 0; j < equiv_points.size(); j++) {
        if(equiv_points[j].find_first(hull_indices[i], k)) {
          for(ii = 0; ii < equiv_points[j].size(); ii++)
            if(ii != k) {
              sout << (*this)[ equiv_points[j][ii]].name << " ";
            }
        }
      }

      sout << endl;
    }

  };

  void write(string s) {
    BP_Write file(s);
    file.newfile();
    write(file.get_ostream());

  };

  void write(ostream &sout) {
    sout << setw(16) << setprecision(9) << "#E_form ";
    sout << setw(16) << setprecision(9) << "weight ";
    sout << setw(16) << setprecision(9) << "concs ";
    sout << setw(16) << setprecision(9) << "dist_from_hull ";
    sout << setw(16) << setprecision(9) << left << "   name " << right;
    sout << "\n";

    // write hull ---------------------------------------
    for(unsigned long int i = 0; i < size(); i++) {
      sout << (*this)[ i] << endl;
    }
  };

  void plot(BP_Plot &plot, string label, string s, int line_style, double line_width, int point_style, double point_size, bool point_face, int component1, int component2) {
    BP_Vec<double> x_list;
    BP_Vec<double> y_list;
    unsigned long int i, j;

    if(hull_found) {
      //for(i=0; i<hull_indices.size(); i++)
      //{
      //	conc.add( get_concentration(hull_indices[i],component));
      //	nrg.add( get_Ef( hull_indices[i]));
      //}
      //
      //plot[0].add_line( conc, nrg, label);
      //plot[0].last().set_lineprops(s, line_style, line_width);			// s color, line, width 1

      BP_Vec<int> nborlist;
      VectorXd v, nbor_v;

      for(i = 0; i < hull.CH_verts_size(); i++) {
        v = hull.CH_verts_pos(i);
        nborlist = hull.CH_verts_nborverts(i);


        for(j = 0; j < nborlist.size(); j++) {
          x_list.clear();
          y_list.clear();

          if(nborlist[j] > i) {

            nbor_v = hull.CH_verts_pos(nborlist[j]);

            //x_list.add( v(1+component));
            //y_list.add( v(0));

            //x_list.add( nbor_v(1+component));
            //y_list.add( nbor_v(0));

            if(component1 == -1) {
              y_list.add(v(0));
              y_list.add(nbor_v(0));
            }
            else {
              y_list.add(v(1 + component1));
              y_list.add(nbor_v(1 + component1));
            }

            if(component2 == -1) {
              x_list.add(v(0));
              x_list.add(nbor_v(0));
            }
            else {
              x_list.add(v(1 + component2));
              x_list.add(nbor_v(1 + component2));
            }

            plot[0].add_line(x_list, y_list, label);
            plot[0].last().set_lineprops(s, line_style, line_width);
          }
        }
      }
    }

    y_list.clear();
    x_list.clear();

    if(component1 == -1) {
      for(i = 0; i < size(); i++) {
        y_list.add(get_Ef(i));
      }
    }
    else {
      for(i = 0; i < size(); i++) {
        y_list.add(get_concentration(i, component1));
      }
    }

    if(component2 == -1) {
      for(i = 0; i < size(); i++) {
        x_list.add(get_Ef(i));
      }
    }
    else {
      for(i = 0; i < size(); i++) {
        x_list.add(get_concentration(i, component2));
      }
    }



    plot[0].add_points(x_list, y_list, label);
    plot[0].last().set_pointprops(s, point_style, point_size, 0);		// s color, circles, diameter 5pts, edge 0pts
    plot[0].last().set_pointface(point_face);

  };

  void write_below_hull(string s, BP_Vec<Energy> energy_list) {

    if(!hull_found) {
      calc_hull(false);
      if(!hull_found)
        return;
    }

    VectorXd v;
    unsigned long int i, j;

    for(i = 0; i < energy_list.size(); i++) {
      v.resize(1 + energy_list[0].concentration.size());
      v(0) = energy_list[i].Ef;
      for(j = 0; j < energy_list[i].concentration.size(); j++) {
        v(j + 1) =  energy_list[i].concentration[j];
      }

      energy_list[i].dist_from_hull = hull.CH_dist_to_hull(v);
      energy_list[i].dist_from_hull_found = true;
    }

    BP_Write below_hull(s);
    below_hull.newfile();

    below_hull << setw(16) << setprecision(9) << "#E_form ";
    below_hull << setw(16) << setprecision(9) << "weight ";
    below_hull << setw(16) << setprecision(9) << "concs ";
    below_hull << setw(16) << setprecision(9) << "dist_from_hull ";
    below_hull << setw(16) << setprecision(9) << left << "   name " << right;
    below_hull << "\n";

    for(i = 0; i < energy_list.size(); i++)
      if(energy_list[i].dist_from_hull < 0.0)
        below_hull << energy_list[i] << endl;

  };

  void weight_nrg(double A, double B, double kT) {
    E_vec_ready = false;

    for(int i = 0; i < size(); i++)
      (*this)[i].weight = A * exp(-(*this)[i].dist_from_hull / kT) + B;

  };

  void weight_ERef(double A, double B, double kT, double ERef) {
    E_vec_ready = false;

    for(int i = 0; i < size(); i++)
      (*this)[i].weight = (((*this)[i].Ef <= ERef) ? 1.0 : (A * exp(-((*this)[i].Ef - ERef) / kT) + B));

  };

  double calc_rms(EnergySet &that, bool include_weights, double dist_from_hull_cutoff, unsigned long int &count) {
    // calculate the rms between *this and 'that'
    //	  - typically expect that *this is the CLEX energies, and 'that' is the DFT energies
    //   assumes weights of *this and 'that' are the same
    //   if 'dist_from_hull_cutoff' is > 0, then only includes those structures whose
    //   energies in 'that' are below the cutoff

    unsigned long int i, j, k;
    count = 0;
    double rms = 0.0;

    for(i = 0; i < size(); i++)
      if(get_weight(i) != 0) {
        if(include_weights) {
          if(dist_from_hull_cutoff <= 0 || that.get_dist_from_hull(i) < dist_from_hull_cutoff) {
            rms += sqr(get_Ef(i) * get_weight(i) - that.get_Ef(i) * that.get_weight(i));
            count++;
          }
        }
        else {
          if(dist_from_hull_cutoff <= 0 || that.get_dist_from_hull(i) < dist_from_hull_cutoff) {
            rms += sqr(get_Ef(i) - that.get_Ef(i));
            count++;
          }
        }
      }

    //cout << "count: " << count << endl;

    rms /= count;
    return sqrt(rms);

  }

};


class ECI {
public:
  int label;
  int weight;
  int mult;
  int size;
  double length;
  double value;
  BP_Vec<int> hierarchy;


  ECI() {
  };

  ECI(int i1, int i2, int i3, int i4, double i5, BP_Vec<int> i6) {
    label = i1;
    weight = i2;
    mult = i3;
    size = i4;
    length = i5;
    hierarchy = i6;
    value = 0;
  };

  friend ostream &operator<<(ostream &outstream, const ECI &e) {
    outstream << setprecision(6) << e.label << " \t" <<  e.weight << " \t" << e.mult << " \t" << e.size << " \t" << e.length << " \t" << e.hierarchy.size();
    for(int i = 0; i < e.hierarchy.size(); i++)
      outstream << " \t" << e.hierarchy[i] ;
    return outstream;
  };

  void clear() {
    label = weight = mult = size = 0;
    length = 0;
    hierarchy.clear();
  };

};


class ECISet : private BP_Vec<ECI> {
private:

  //BP_Vec<ECI> eci_list;
  mutable int Nclust_on;
  mutable bool Nclust_set;
  BP_Vec<int> fix_pos_list;		// list of indices of 'fixed' clusters
  BP_Vec<bool> fix_val_list;		// value 'fixed' clusters are fixed to

  double cv;
  double rms;

public:

  ECISet() {
    cv = -1;
    rms = -1;
  };

  ECISet(string eci_in_filename) {
    Nclust_set = 0;
    cv = -1;
    rms = -1;

    BP_Parse file(eci_in_filename);
    ECI eci;
    // header line
    file.getline();

    BP_Vec<string> entry;
    // data
    do {
      entry = file.getline_string();
      if(entry.size() != 0) {
        eci.clear();
        eci.label = BP::stoi(entry[0]);
        if(entry[1] == "FixOn") {
          eci.weight = 1;
          fix((*this).size(), 1);
        }
        else if(entry[1] == "FixOff") {
          eci.weight = 0;
          fix((*this).size(), 0);
        }
        else
          eci.weight = BP::stoi(entry[1]);
        eci.mult = BP::stoi(entry[2]);
        eci.size = BP::stoi(entry[3]);
        eci.length = BP::stod(entry[4]);
        for(int i = 6; i < entry.size(); i++)
          eci.hierarchy.add(BP::stoi(entry[i]));
        add(eci);
      }
    }
    while(file.eof() == false);

    get_Nclust_on();

    //cout << "eci_set.size(): " << eci_list.size() << endl;
    //for( int i=0; i<eci_list.size(); i++)
    //	cout << eci_list[i] << endl;
  };

  //ECISet& operator=(const ECISet &eci_set)
  //{
  //	cout << " begin ECISet::operator=ECISet" << endl;
  //	(*this) = eci_set;
  //
  //	return (*this);
  //};

  int get_Nclust_on() const {
    if(Nclust_set) {
      return Nclust_on;
    }
    else {
      Nclust_on = 0;
      for(int i = 0; i < size(); i++) {
        if((*this)[i].weight != 0)
          Nclust_on++;
      }
      Nclust_set = 1;
      return Nclust_on;
    }
  };

  void set_values(QXu::Array &ECI) {
    int in_i = 0;
    for(int i = 0; i < size(); i++)
      if((*this)[i].weight != 0) {
        (*this)[i].value = ECI.zelem(in_i, 0);
        in_i++;
      }
  };

  void set_values(const VectorXd &ECI) {
    int in_i = 0;
    for(int i = 0; i < size(); i++)
      if((*this)[i].weight != 0) {
        (*this)[i].value = ECI(in_i);
        in_i++;
      }
  };

  void set_values_and_weights(const MatrixXd &ECI) {
    clear_weights();

    unsigned long int i;
    double prec = 1e-8;

    for(i = 0; i < size(); i++) {
      if(fabs(ECI(i)) > prec) {
        set_clust_on(i);
        (*this)[i].value = ECI(i);
      }
      else {
        set_clust_off(i);
      }
    }
  };

  void clear_weights() {
    for(int i = 0; i < size(); i++)
      (*this)[i].weight = 0;
    Nclust_on = 0;
    Nclust_set = 1;
  };

  void set_clust_off(int i) {
    if((*this)[i].weight != 0) {
      Nclust_on--;
      (*this)[i].weight = 0;
    }
  };

  void set_clust_on(int i) {
    if((*this)[i].weight == 0) {
      Nclust_on++;
      (*this)[i].weight = 1;
    }
  };

  void toggle_clust(int i) {
    if((*this)[i].weight == 0) {
      Nclust_on++;
      (*this)[i].weight = 1;
    }
    else {
      Nclust_on--;
      (*this)[i].weight = 0;
    }
  };

  ECISet &operator=(const BP_Comb &comb) {
    if(comb.size() == size()) {
      for(int i = 0; i < comb.size(); i++)
        (*this)[i].weight = comb[i];
    }
    Nclust_set = 0;

    return (*this);
  };

  string get_bit_string() const {
    stringstream ss;
    for(int i = 0; i < size(); i++)
      ss << (*this)[i].weight;
    return ss.str();
  }

  void set_bit_string(string s) {
    if(s.size() != size())
      return;

    Nclust_on = 0;
    Nclust_set = 1;
    for(int i = 0; i < s.size(); i++)
      if(s[i] == '0') {
        (*this)[i].weight = 0;
      }
      else {
        (*this)[i].weight = 1;
        Nclust_on++;
      }

  }

  unsigned long int size() const {
    return BP_Vec<ECI>::size();
  };

  int get_weight(unsigned long int i) const {
    return (*this)[i].weight;
  };
  int get_size(unsigned long int i) const {
    return (*this)[i].size;
  };
  double get_value(unsigned long int i) const {
    return (*this)[i].value;
  };
  double get_cv() const {
    return cv;
  };
  double get_rms() const {
    return rms;
  };
  void set_cv(double i1) {
    cv = i1;
  };
  void set_rms(double i1) {
    rms = i1;
  };

  bool fix_ok() {
    for(int i = 0; i < fix_pos_list.size(); i++) {
      if((*this)[ fix_pos_list[i]].weight != fix_val_list[i])
        return false;
    }
    return true;
  };

  void fix(int i1, bool i2) {
    fix_pos_list.add(i1);
    fix_val_list.add(i2);
  };

  void unfix(int i1) {
    for(int i = 0; i < fix_pos_list.size(); i++) {
      if(fix_pos_list[i] == i1) {
        fix_pos_list.remove(i);
        fix_val_list.remove(i);
      }
    }
  };

  void unfix_all() {
    fix_pos_list.erase();
    fix_val_list.erase();
  };

  void set_fix() {
    for(int i = 0; i < fix_pos_list.size(); i++) {
      (*this)[ fix_pos_list[i]].weight = fix_val_list[i];
    }
    Nclust_set = 0;
  };

  void set_Rmax_fix(double Rmax) {
    for(int i = 0; i < size(); i++) {
      if((*this)[ i].length > Rmax)
        fix(i, 0);
    }
  };

  void write_ECI(ostream &sout) {
    for(int i = 0; i < size(); i++) {
      if(get_weight(i) != 0) {
        sout << "i: " << i << " value: " << get_value(i) << endl;

      }

    }

  };

  void write_ECIin(string filename) {
    BP_Write file(filename);
    file.newfile();

    file << "label" << "     " << "weight" << "     " << "mult" << "     " << "size" << "     " << "length" << "     " << "heirarchy" << endl;
    for(int i = 0; i < size(); i++)
      file << (*this)[i] << endl;
  };

  void write_ECIin(BP_Vec<int> &index_list, string filename) {
    BP_Write file(filename);
    file.newfile();

    file << "label" << "     " << "weight" << "     " << "mult" << "     " << "size" << "     " << "length" << "     " << "heirarchy" << endl;
    for(int i = 0; i < index_list.size(); i++)
      file << (*this)[index_list[i] ] << endl;
  };

  void write_ECIout(string filename, EnergySet &nrg_set) {
    BP_Write file(filename);
    file.newfile();

    file << "The number of clusters in fit is " << get_Nclust_on() << endl;
    file << "The number of structure in fit is " << nrg_set.get_Nstruct_on() << endl;
    file << "The rms is " << setprecision(9) << get_rms() << endl;
    file << "Number of structures in WCV calculation is Stru_wcv=" << nrg_set.get_Nstruct_on() << endl;
    file << "The WCV score is " << get_cv() << endl;
    file << "The rms based on the " << "N/A" << " structures having vasp E but not in the fit is: " << "N/A" << endl;
    file << setiosflags(ios::right) << setw(15) << "ECIs" << setiosflags(ios::right) << setw(15) << "ECI/mult" << setiosflags(ios::right) << setw(15) << "Cluster#" << endl;
    for(int i = 0; i < size(); i++) {
      if(get_weight(i) != 0) {
        file << setiosflags(ios::right) << setiosflags(ios::fixed) << setw(15) << setprecision(9) << (*this)[i].value ;
        file << setiosflags(ios::right) << setiosflags(ios::fixed) << setw(15) << setprecision(9) << (*this)[i].value / (*this)[i].mult ;
        file << setiosflags(ios::right) << setiosflags(ios::fixed) << setw(15) << setprecision(9) << (*this)[i].label << endl;
      }
    }

  };

  void write_ECIout(BP_Vec<int> &index_list, string filename, EnergySet &nrg_set) {
    BP_Write file(filename);
    file.newfile();

    file << "The number of clusters in fit is " << get_Nclust_on() << endl;
    file << "The number of structure in fit is " << nrg_set.get_Nstruct_on() << endl;
    file << "The rms is " << setprecision(9) << get_rms() << endl;
    file << "Number of structures in WCV calculation is Stru_wcv=" << nrg_set.get_Nstruct_on() << endl;
    file << "The WCV score is " << get_cv() << endl;
    file << "The rms based on the " << "N/A" << " structures having vasp E but not in the fit is: " << "N/A" << endl;
    file << setiosflags(ios::right) << setw(15) << "ECIs" << setiosflags(ios::right) << setw(15) << "ECI/mult" << setiosflags(ios::right) << setw(15) << "Cluster#" << endl;
    for(int i = 0; i < size(); i++) {
      if(get_weight(i) != 0) {
        file << setiosflags(ios::right) << setiosflags(ios::fixed) << setw(15) << setprecision(9) << (*this)[i].value ;
        file << setiosflags(ios::right) << setiosflags(ios::fixed) << setw(15) << setprecision(9) << (*this)[i].value / (*this)[i].mult ;
        file << setiosflags(ios::right) << setiosflags(ios::fixed) << setw(15) << setprecision(9) << index_list[i] << endl;
      }
    }

  };

  bool operator==(const ECISet &E) {
    if(size() != E.size())
      return false;

    for(int i = 0; i < size(); i++)
      if((*this)[i].weight != E.get_weight(i))
        return false;

    return true;
  };

  void randomize(int Nclust, MTRand &mtrand) {
    unsigned long int i;
    clear_weights();
    set_fix();
    do {
      i = mtrand.randInt(size() - 1);
      set_clust_on(i);
      if(!fix_ok())
        set_clust_off(i);

    }
    while(get_Nclust_on() != Nclust);

  };

};

class GeneticAlgorithm {
private:
  int Nmin;
  int Nmax;
  int Nchildren;
  //int Ngenerations;
  int Nmutations;
  BP_GVec< ECISet> population;

public:

  GeneticAlgorithm() {

  };

  GeneticAlgorithm(BP_Vec< ECISet> &wannabeparents, int min, int max, int children, int mutations) {
    set_population(wannabeparents);
    Nmin = min;
    Nmax = max;
    Nchildren = children;
    Nmutations = mutations;
  };

  ECISet &operator[](unsigned long int i) {
    return population[i];
  };
  unsigned long int size() {
    return population.size();
  };



  ECISet mate(const ECISet &Mom, const ECISet &Dad, MTRand &mtrand) {
    //cout << "begin GeneticAlgorithm:mate()" << endl;
    unsigned long int j;
    // mate
    //cout << "    mate" << endl;
    ECISet child = Mom;
    for(int i = 0; i < Mom.size(); i++) {
      if(Mom.get_weight(i) != Dad.get_weight(i)) {
        if(mtrand.randExc() < 0.5)
          child.set_clust_off(i);
        else
          child.set_clust_on(i);
      }
    }

    // mutate (random with probability 1/Nmutations)
    //cout << "    mutate" << endl;
    for(int i = 0; i < child.size(); i++) {
      if(mtrand.randExc() < (1.0 * Nmutations) / child.size()) {
        child.toggle_clust(i);
        if(!child.fix_ok())
          child.toggle_clust(i);
      }
    }

    // enforce Nmin limit
    //cout << "    enforce Nmin" << endl;
    while(child.get_Nclust_on() < Nmin) {
      j = mtrand.randInt(child.size() - 1);
      child.set_clust_on(j);
      if(!child.fix_ok())
        child.set_clust_off(j);
    }

    // enforce Nmax limit
    //cout << "    enforce Nmax" << endl;
    while(child.get_Nclust_on() > Nmax) {
      j = mtrand.randInt(child.size() - 1);
      child.set_clust_off(j);
      if(!child.fix_ok())
        child.set_clust_on(j);
    }
    //cout << "finish GeneticAlgorithm:mate()" << endl;
    return child;
  };

  void set_population(BP_Vec< ECISet> &wannabeparents) {
    for(int i = 0; i < wannabeparents.size(); i++) {
      population.add(wannabeparents[i]);
    }
  };

  void set(int min, int max, int children, int mutations) {
    Nmin = min;
    Nmax = max;
    Nchildren = children;
    Nmutations = mutations;
  };

  void set_Nmin(int i) {
    Nmin = i;
  };
  void set_Nmax(int i) {
    Nmax = i;
  };
  void set_Nchildren(int i) {
    Nchildren = i;
  };
  void set_Nmutations(int i) {
    Nmutations = i;
  };

  BP_Vec<ECISet> get_population() {
    BP_Vec<ECISet> out_pop;

    for(int i = 0; i < population.size(); i++) {
      out_pop.add(population[i]);
    }
    return out_pop;
  };

  int get_Nmin() {
    return Nmin;
  };
  int get_Nmax() {
    return Nmax;
  };
  int get_Nchildren() {
    return Nchildren;
  };
  int get_Nmutations() {
    return Nmutations;
  };

  int get_best_index() {
    int best;
    double min = 1.0e20;
    for(int i = 0; i < population.size(); i++) {
      if(population[i].get_cv() < min) {
        best = i;
        min = population[i].get_cv();
      }
    }

    return best;
  };


  void run(MTRand &mtrand, Correlation &corr, EnergySet &nrg_set, bool print_steps) {
    run_ga(mtrand, corr, nrg_set, print_steps, 0, 0);
  };

  void run_dir(MTRand &mtrand, Correlation &corr, EnergySet &nrg_set, bool print_steps) {
    run_ga(mtrand, corr, nrg_set,  print_steps, 1, 0);
  };

  void run_dfs(long int Nstop, MTRand &mtrand, Correlation &corr, EnergySet &nrg_set, bool print_steps) {
    run_ga(mtrand, corr, nrg_set,  print_steps, 2, Nstop);
  };

private:
  void run_ga(MTRand &mtrand, Correlation &corr, EnergySet &nrg_set, bool print_steps, int min_option, long int Nstop) {
    //cout << "begin GeneticAlgorithm::run_ga() " << endl;
    int i;
    double fitness;
    bool singular;
    ECISet child;
    ECISet eci_best;
    BP_GVec_Member<ECISet> *Mom;
    BP_GVec_Member<ECISet> *Dad;
    BP_GVec_Member<ECISet> *worst_parent;
    BP_GVec_Member<ECISet> *new_addition;


    // add population to generator, make 'rate' = exp(1/cv) for weighting parent choice
    //cout << "add to generator" << endl;
    BP_RVG_tree<ECISet> generator;
    for(i = 0; i < population.size(); i++) {
      fit(corr, nrg_set, population[i], singular);
      if(eci_best.get_cv() == -1 || population[i].get_cv() < eci_best.get_cv())
        eci_best = population[i];
      //cout << " i: " << population[i].get_cv() << " bestsofar: " << eci_best.get_cv() << endl;
      //generator.add( population.member(i), exp(1.0/population[i].get_cv()));
      //generator.add( population.member(i), 1.0/population[i].get_cv());
      generator.add(population.member(i), 1.0);
    }

    // mate
    //cout << "run mate loop" << endl;
    fitness = max_cv(population, worst_parent);
    for(i = 0; i < Nchildren; i++) {
      //cout << "--pick parents" << endl;
      Mom = generator.pick(mtrand.randExc());
      do {
        Dad = generator.pick(mtrand.randExc());
      }
      while(Dad == Mom);

      //cout << "  mate parents" << endl;
      child = mate(Mom->get_val(), Dad->get_val(), mtrand);
      fit(corr, nrg_set, child, singular);

      //cout << "  check child fitness" << endl;




      //cout << "  child cv: " << child.get_cv() << endl;

      //cout << "  compare to pop" << endl;
      if(child.get_cv() != -1)
        if(child.get_cv() < fitness) {
          //new_addition = add_once( population, child);

          if(min_option == 0) {

          }
          else if(min_option == 1) {
            child = direct_min(Nmin, Nmax, child, corr, nrg_set, print_steps);

          }
          else if(min_option == 2) {
            child = dfs_min(Nstop, Nmin, Nmax, child, corr, nrg_set, print_steps);

          }

          if(eci_best.get_cv() == -1 || child.get_cv() < eci_best.get_cv()) {
            eci_best = child;
          }

          new_addition = add_once(population, child);

          if(new_addition != NULL) {


            //cout << "    add to pop" << endl;
            //cout << child.get_bit_string() << "    GA: " << i << "  cv: " << child.get_cv() << " rms: " << child.get_rms() << endl;
            population.remove(worst_parent);
            //generator.add(new_addition, exp( 1.0/child.get_cv()));
            //generator.add(new_addition, 1.0/child.get_cv() );
            generator.add(new_addition, 1.0);
            fitness = max_cv(population, worst_parent);
          }
        }


      cout << child.get_bit_string() << "  child: " << i << " Nclust: " << child.get_Nclust_on() << "  cv: " << child.get_cv() << " rms: " << child.get_rms() << " bestsofar: " << eci_best.get_cv() << " max_cv: " << fitness << endl;

    }

    //cout << "finish GeneticAlgorithm::run_ga() " << endl;

  };
};

void Correlation::write_covar(const ECISet &eci_in, const EnergySet &nrg_set, string filename) const {
  //cout << "begin fit()" << endl;

  int Nstruct = nrg_set.get_Nstruct_on();
  int Nclust = eci_in.get_Nclust_on();

  //cout << "fit Nstruct: " << Nstruct << "  Nclust: " << Nclust << endl;
  MatrixXd corr_matrix(Nstruct, Nclust);

  //cout << " set correlation matrix" << endl;
  set_correlation_matrix(corr_matrix, (*this), nrg_set, eci_in);

  //cout << " init svd" << endl;
  JacobiSVD<MatrixXd> svd(corr_matrix, ComputeThinU | ComputeThinV);

  BP_Write file(filename);
  file.newfile();
  file << "Singular values of correlation matrix used in fit:" << endl;
  for(int i = 0; i < Nclust; i++) {
    file << setiosflags(ios::right | ios::showpoint | ios::fixed) << setw(12) << setprecision(7) << svd.singularValues()[i] << ' ';
  }
  file << ' ' << endl;
  file << "Right singular vectors of correlation matrix used in fit:" << endl;
  for(int i = 0; i < Nclust; i++) {
    for(int j = 0; j < Nclust; j++) {
      file << setiosflags(ios::right | ios::showpoint | ios::fixed) << setw(12) << setprecision(7) << svd.matrixV()(i, j) << ' ';
    }
    file << ' ' << endl;
  }


}

void Correlation::set_cv_a(const MatrixXd &C, const EnergySet &nrg_set) {
  //cout << "begin set_cv_a()" << endl;
  // -- compute a_i = X_i*((X^T*X)^-1)*X_i^T ahead of time, once per corr

  MatrixXd XTXinv;
  unsigned long int i, ii;
  double a_i;
  cv_a.clear();
  cv_a.capacity(nrg_set.get_Nstruct_on());

  XTXinv = (C.transpose() * C).inverse();

  ii = 0;
  for(i = 0; i < nrg_set.size(); i++)
    if(nrg_set.get_weight(i) != 0) {
      a_i = C.row(ii) * XTXinv * C.row(ii).transpose();
      cv_a.add(a_i);
      ii++;
    }


  cv_a_ready = true;
  //cout << "finish set_cv_a()" << endl;

};

void EnergySet::calc_clex(const Correlation &corr, const ECISet &eci) {
  unsigned long int i, j;
  for(i = 0; i < size(); i++) {
    (*this)[i].Ef = 0;
    for(j = 0; j < corr[i].size(); j++) {
      (*this)[i].Ef += corr[i][j] * eci.get_value(j);
      (*this)[i].dist_from_hull = 0.0;
    }
  }
};

void EnergySet::calc_nrg_diff(BP_Vec<int> &index_list, const Correlation &corr, const ECISet &eci) {
  // subtract corr*eci from energies, where index_list gives the columns in corr that correspond to each eci

  if(index_list.size() != eci.size()) {
    cout << "Error.  In calc_nrg_diff(), index_list.size() (= " << index_list.size() << ") != eci.size() (= " << eci.size() << ")" << endl;
    exit(1);
  }

  unsigned long int i, j;
  for(i = 0; i < size(); i++) {
    //cout << "Ef: " << (*this)[i].Ef << endl;
    double tmp = 0;
    for(j = 0; j < index_list.size(); j++) {
      //cout << "j: " << j << "  corr: " << corr[i][index_list[j] ] << "  weight: " << eci.get_weight(j) << "  value: " << eci.get_value(j) << "  contribution: " << corr[i][index_list[j] ]*eci.get_weight(j)*eci.get_value(j) << endl;
      tmp += corr[i][index_list[j] ] * eci.get_weight(j) * eci.get_value(j);
      (*this)[i].Ef -= corr[i][index_list[j] ] * eci.get_weight(j) * eci.get_value(j);
      (*this)[i].dist_from_hull = 0.0;
    }

    //cout << endl << "   Ef(A): " << tmp << "   diff: " << (*this)[i].Ef << endl;
  }
};


////----------------------------
/// main functions
void calc_eci(string energy_filename, string eci_in_filename, string corr_in_filename, BP_Vec<ECISet> &population) {
  Correlation corr(corr_in_filename);
  EnergySet DFT_nrg(energy_filename);
  ECISet eci_in(eci_in_filename);


  bool singular;

  if(population.size() == 1)
    eci_in = population[0];

  fit(corr, DFT_nrg, eci_in, singular);

  if(!singular) {
    cout << eci_in.get_bit_string() << endl;
    cout << "cv: " << eci_in.get_cv() << " rms: " << eci_in.get_rms() << endl;
  }
  else {
    cout << "error, singular" << endl;
  }

  cout << endl << endl;

  // write corr.covar
  corr.write_covar(eci_in, DFT_nrg, "corr.covar");

  // write eci.in
  eci_in.write_ECIin("eci.in");
  cout << "Wrote 'eci.in'" << endl;

  // write eci.out
  eci_in.write_ECIout("eci.out", DFT_nrg);
  cout << "Wrote 'eci.out'" << endl;

  //cout << "Hull of: " << energy_filename << endl;
  //DFT_nrg.write_hull(cout);
  cout << endl << "Calculating hull of: " << energy_filename << endl;
  DFT_nrg.calc_hull(true);
  DFT_nrg.write_hull("hull");
  cout << "Wrote 'hull'" << endl;


  if(eci_in.get_cv() != 1e20) {

    // calc energy_clex
    EnergySet energy_clex = DFT_nrg;
    energy_clex.calc_clex(corr, eci_in);

    // write energy_clex
    cout << endl << "Calculating hull of: energy.clex " << endl;
    energy_clex.calc_hull(true);
    //energy_clex.write_hull(cout);
    energy_clex.write_hull("hull.clex");
    cout << "Wrote 'hull.clex'" << endl;
    energy_clex.write("energy.clex");
    cout << "Wrote 'energy.clex'" << endl;

    cout << endl;

    unsigned long int count;
    double rms;

    rms = energy_clex.calc_rms(DFT_nrg, 1, -1, count);
    cout << "weighted total rms: " <<  rms <<  "    #structures: " << count << endl;

    rms = energy_clex.calc_rms(DFT_nrg, 0, -1, count);
    cout << "non-weighted total rms: " <<  rms << "    #structures: " << count << endl;

    rms = energy_clex.calc_rms(DFT_nrg, 0, 0.2, count);
    cout << "non-weighted rms within 0.2 of hull: " << rms << "    #structures: " << count << endl;

    rms = energy_clex.calc_rms(DFT_nrg, 0, 0.1, count);
    cout << "non-weighted rms within 0.1 of hull: " << rms << "    #structures: " << count << endl;

    rms = energy_clex.calc_rms(DFT_nrg, 0, 0.05, count);
    cout << "non-weighted rms within 0.05 of hull: " << rms << "    #structures: " << count << endl;

    rms = energy_clex.calc_rms(DFT_nrg, 0, 0.01, count);
    cout << "non-weighted rms within 0.01 of hull: " << rms << "    #structures: " << count << endl;

    rms = energy_clex.calc_rms(DFT_nrg, 0, 0.005, count);
    cout << "non-weighted rms within 0.005 of hull: " << rms << "    #structures: " << count << endl;

    rms = energy_clex.calc_rms(DFT_nrg, 0, 0.001, count);
    cout << "non-weighted rms within 0.001 of hull: " << rms << "    #structures: " << count << endl;

    EnergySet energy_clex_of_DFT_hull;

    if(DFT_nrg.is_hull_found()) {
      energy_clex_of_DFT_hull.set(energy_clex.get(DFT_nrg.get_hull_indices()));
      cout << endl << "Calculating hull of: clex_of_DFT_hull" << endl;
      energy_clex_of_DFT_hull.calc_hull(true);
      energy_clex_of_DFT_hull.write_hull("hull.clex_of_DFT_hull");
      cout << "Wrote 'hull.clex_of_DFT_hull'" << endl;
      energy_clex_of_DFT_hull.write_below_hull("below.hull", energy_clex.get(DFT_nrg.get_non_hull_indices()));
      cout << "Wrote 'below.hull'" << endl;
      //energy_clex_of_DFT_hull.plot(plot, "DFT_hull_clex", "green", 1, 0.5, 0, 2, true, i, j);
    }

    cout << endl;
    for(int i = -1; i < (int) DFT_nrg.get_concentration_size(); i++) {
      for(int j = i + 1; j < (int) DFT_nrg.get_concentration_size(); j++) {
        // plot energy and energy.clex
        BP_Plot plot;
        //plot( BP_Plot &plot, string label, string s, int line_style, double line_width, int point_style, double point_size, bool point_face)
        DFT_nrg.plot(plot, "DFT", "red", 0, 1, 0, 5, false, i, j);
        energy_clex.plot(plot, "CLEX", "blue", 0, 0.5, 0, 3, true, i, j);

        // write below.hull
        if(DFT_nrg.is_hull_found()) {
          //EnergySet energy_clex_of_DFT_hull( energy_clex.get(DFT_nrg.get_hull_indices()));
          //cout << "Calculating hull of: clex_of_DFT_hull" << endl;
          //energy_clex_of_DFT_hull.calc_hull(true);
          //energy_clex_of_DFT_hull.write_hull("hull.clex_of_DFT_hull");
          //cout << "  wrote 'hull.clex_of_DFT_hull'" << endl;
          //energy_clex_of_DFT_hull.write_below_hull( "below.hull", energy_clex.get(DFT_nrg.get_non_hull_indices()));
          energy_clex_of_DFT_hull.plot(plot, "DFT_hull_clex", "green", 1, 0.5, 0, 2, true, i, j);
        }

        string si, sj;
        if(i == -1) si = "nrg";
        else si = "x" + itos(i);

        if(j == -1) sj = "nrg";
        else sj = "x" + itos(j);



        plot.write("clex_results_" + si + "_vs_" + sj);
        cout << "Wrote '" << "clex_results_" + si + "_vs_" + sj << "'" << endl;
      }
    }

  }

}

void calc_cs_eci(string energy_filename, string eci_in_filename, string corr_in_filename, const BP_Vec<double> &mu, int alg) {
  // solve E = Corr*ECI, for ECI using compressive sensing L1 norm minimization (plus L2 norm in practice)

  Correlation corr(corr_in_filename);
  EnergySet DFT_nrg(energy_filename);
  ECISet eci_in(eci_in_filename);
  ECISet eci_CS = eci_in;

  unsigned long int i;
  double prec_shrink = 1e-6;
  double prec_ECI = 1e-6;
  bool print_steps = true;

  for(i = 0; i < mu.size(); i++) {
    if(alg == 0) {
      // use Fixed-point continuation method, calc for single value of mu
      eci_CS = calc_FPC_eci(mu[i], prec_shrink, prec_ECI, eci_in, corr, DFT_nrg, print_steps);
    }
    else if(alg == 1) {
      // use Bregman iteration method, calc for single value of mu
      eci_CS = calc_BI_eci(mu[i], prec_shrink, prec_ECI, eci_in, corr, DFT_nrg, print_steps);
    }
    else {
      cout << "Error in calc_cs_eci().  alg == '" << alg << "' is not a valid option." << endl;
      cout << "  Options are:  0, Fixed-point continuation" << endl;
      cout << "                1, Bergman iteration" << endl;
      exit(1);
    }

    cout << endl << eci_CS.get_bit_string() << "    mu: " << mu[i] << "  Nclust: " << eci_CS.get_Nclust_on() << " rms: " << eci_CS.get_rms() << endl;
  }

  // write eci.in
  eci_CS.write_ECIin("eci.in");
  cout << "Wrote 'eci.in'" << endl;

  // write eci.out
  eci_CS.write_ECIout("eci.out", DFT_nrg);
  cout << "Wrote 'eci.out'" << endl;

  //cout << "Hull of: " << energy_filename << endl;
  //DFT_nrg.write_hull(cout);
  cout << endl << "Calculating hull of: " << energy_filename << endl;
  DFT_nrg.calc_hull(true);
  DFT_nrg.write_hull("hull");
  cout << "Wrote 'hull'" << endl;


  //if( eci_CS.get_cv() != 1e20)
  {

    // calc energy_clex
    EnergySet energy_clex = DFT_nrg;
    energy_clex.calc_clex(corr, eci_CS);

    // write energy_clex
    cout << endl << "Calculating hull of: energy.clex " << endl;
    energy_clex.calc_hull(true);
    //energy_clex.write_hull(cout);
    energy_clex.write_hull("hull.clex");
    cout << "Wrote 'hull.clex'" << endl;
    energy_clex.write("energy.clex");
    cout << "Wrote 'energy.clex'" << endl;


    cout << endl;

    unsigned long int count;
    double rms;

    rms = energy_clex.calc_rms(DFT_nrg, 1, -1, count);
    cout << "weighted total rms: " <<  rms <<  "    #structures: " << count << endl;

    rms = energy_clex.calc_rms(DFT_nrg, 0, -1, count);
    cout << "non-weighted total rms: " <<  rms << "    #structures: " << count << endl;

    rms = energy_clex.calc_rms(DFT_nrg, 0, 0.2, count);
    cout << "non-weighted rms within 0.2 of hull: " << rms << "    #structures: " << count << endl;

    rms = energy_clex.calc_rms(DFT_nrg, 0, 0.1, count);
    cout << "non-weighted rms within 0.1 of hull: " << rms << "    #structures: " << count << endl;

    rms = energy_clex.calc_rms(DFT_nrg, 0, 0.05, count);
    cout << "non-weighted rms within 0.05 of hull: " << rms << "    #structures: " << count << endl;

    rms = energy_clex.calc_rms(DFT_nrg, 0, 0.01, count);
    cout << "non-weighted rms within 0.01 of hull: " << rms << "    #structures: " << count << endl;

    rms = energy_clex.calc_rms(DFT_nrg, 0, 0.005, count);
    cout << "non-weighted rms within 0.005 of hull: " << rms << "    #structures: " << count << endl;

    rms = energy_clex.calc_rms(DFT_nrg, 0, 0.001, count);
    cout << "non-weighted rms within 0.001 of hull: " << rms << "    #structures: " << count << endl;


    EnergySet energy_clex_of_DFT_hull;

    if(DFT_nrg.is_hull_found()) {
      energy_clex_of_DFT_hull.set(energy_clex.get(DFT_nrg.get_hull_indices()));
      cout << endl << "Calculating hull of: clex_of_DFT_hull" << endl;
      energy_clex_of_DFT_hull.calc_hull(true);
      energy_clex_of_DFT_hull.write_hull("hull.clex_of_DFT_hull");
      cout << "Wrote 'hull.clex_of_DFT_hull'" << endl;
      energy_clex_of_DFT_hull.write_below_hull("below.hull", energy_clex.get(DFT_nrg.get_non_hull_indices()));
      cout << "Wrote 'below.hull'" << endl;
      //energy_clex_of_DFT_hull.plot(plot, "DFT_hull_clex", "green", 1, 0.5, 0, 2, true, i, j);
    }

    cout << endl;
    for(int i = -1; i < (int) DFT_nrg.get_concentration_size(); i++) {
      for(int j = i + 1; j < (int) DFT_nrg.get_concentration_size(); j++) {
        // plot energy and energy.clex
        BP_Plot plot;
        //plot( BP_Plot &plot, string label, string s, int line_style, double line_width, int point_style, double point_size, bool point_face)
        DFT_nrg.plot(plot, "DFT", "red", 0, 1, 0, 5, false, i, j);
        energy_clex.plot(plot, "CLEX", "blue", 0, 0.5, 0, 3, true, i, j);

        // write below.hull
        if(DFT_nrg.is_hull_found()) {
          //EnergySet energy_clex_of_DFT_hull( energy_clex.get(DFT_nrg.get_hull_indices()));
          //cout << "Calculating hull of: clex_of_DFT_hull" << endl;
          //energy_clex_of_DFT_hull.calc_hull(true);
          //energy_clex_of_DFT_hull.write_hull("hull.clex_of_DFT_hull");
          //cout << "  wrote 'hull.clex_of_DFT_hull'" << endl;
          //energy_clex_of_DFT_hull.write_below_hull( "below.hull", energy_clex.get(DFT_nrg.get_non_hull_indices()));
          energy_clex_of_DFT_hull.plot(plot, "DFT_hull_clex", "green", 1, 0.5, 0, 2, true, i, j);
        }

        string si, sj;
        if(i == -1) si = "nrg";
        else si = "x" + itos(i);

        if(j == -1) sj = "nrg";
        else sj = "x" + itos(j);



        plot.write("clex_results_" + si + "_vs_" + sj);
        cout << "Wrote '" << "clex_results_" + si + "_vs_" + sj << "'" << endl;
      }
    }

  }

}

void calc_ecistats(string energy_filename, string eci_in_filename, string corr_in_filename, BP_Vec<ECISet> &population) {
  Correlation corr(corr_in_filename);
  EnergySet DFT_nrg(energy_filename);
  //ECISet eci_in(eci_in_filename);
  BP_Vec<double> eci_value_list;
  BP_Vec<double> eci_nonzero_value_list;

  bool singular;

  for(int i = 0; i < population.size(); i++) {
    fit(corr, DFT_nrg, population[i], singular);
  }

  cout << "label" << "     " << "frac_on" << "     " << "mean(nonzero)" << "     "  << "min(nonzero)" << "     " << "max(nonzero)" << "     "  << "rms(nonzero)" << "     " << "mean" << "     "  << "min" << "     " << "max" << "     "  << "rms" << endl;
  for(int i = 0; i < population[0].size(); i++) {
    eci_value_list.clear();
    eci_nonzero_value_list.clear();
    for(int j = 0; j < population.size(); j++) {
      eci_value_list.add(population[j].get_value(i));
      if(population[j].get_weight(i)) {
        eci_nonzero_value_list.add(population[j].get_value(i));
      }
    }

    cout << i << " \t";
    cout << setprecision(6) << (1.0 * eci_nonzero_value_list.size()) / (1.0 * eci_value_list.size()) << " \t";
    if(eci_nonzero_value_list.size() != 0) {
      cout << mean(eci_nonzero_value_list) << " \t";
      cout << min(eci_nonzero_value_list) << " \t";
      cout << max(eci_nonzero_value_list) << " \t";
      cout << rms(eci_nonzero_value_list) << " \t";
      cout << mean(eci_value_list) << " \t";
      cout << min(eci_value_list) << " \t";
      cout << max(eci_value_list) << " \t";
      cout << rms(eci_value_list) << " \n";
    }
    else {
      cout << "-" << " \t";
      cout << "-" << " \t";
      cout << "-" << " \t";
      cout << "-" << " \t";
      cout << "-" << " \t";
      cout << "-" << " \t";
      cout << "-" << " \t";
      cout << "-" << " \n";
    }

  }

}

void calc_all_eci(int N, string energy_filename, string eci_in_filename, string corr_in_filename) {

  Correlation corr(corr_in_filename);
  EnergySet DFT_nrg(energy_filename);
  ECISet eci_in(eci_in_filename);
  ECISet eci_min_B = eci_in;
  bool singular;
  int bestsofar = 0;

  MTRand mtrand;

  BP_Comb combs(eci_in.size(), N);

  cout << "This calculates all ecisets with " << N << " eci." << endl;
  cout << "   For N = " << N << " and eci.in size = " << eci_in.size() << " that will be " << combs.total_combs() << " ecisets." << endl << endl;

  //for( int i=0; i<10; i++)
  //	combs.fix( i, 1);

  int count = 0;

  do {
    eci_in = combs;

    fit(corr, DFT_nrg, eci_in, singular);

    if(eci_min_B.get_cv() == -1 || eci_in.get_cv() < eci_min_B.get_cv()) {
      eci_min_B = eci_in;
      bestsofar = combs.get_count();
    }

    if(!singular) {
      cout << eci_in.get_bit_string() << "   i: " << combs.get_count() << "\t Nclust: " << eci_in.get_Nclust_on() << " cv: " << eci_in.get_cv() << " rms: " << eci_in.get_rms() << endl;
      //cout << "i: " << combs.get_count() << " Nclust: " << eci_in.get_Nclust_on() << " cv: " << cv_score << " rms: " << rms_score << " bit_string: " << combs.get_bit_string() << endl;
      //cout << "   eci: " ;
      //for( int i=0; i<eci_in.size(); i++)
      //	if( eci_in[i].weight != 0)
      //		cout << eci_in[i].value << " " ;
      //cout << endl;
    }
    else {
      cout << "error, singular" << endl;
    }

    combs.increment();

  }
  while(combs.complete() == false);

  cout << eci_min_B.get_bit_string() << "       BestSoFar: " << bestsofar << " Nclust: " << eci_min_B.get_Nclust_on() << " cv: " << eci_min_B.get_cv() << " rms: " << eci_min_B.get_rms() << endl;

}

void calc_directmin_eci(int Nrand, int Nmin, int Nmax, string energy_filename, string eci_in_filename, string corr_in_filename, BP_Vec<ECISet> &population) {
  Correlation corr(corr_in_filename);
  EnergySet DFT_nrg(energy_filename);
  ECISet eci_in(eci_in_filename);

  double cv_score;
  double rms_score;
  bool singular;

  ECISet eci_min_A = eci_in;
  ECISet eci_min_B = eci_in;

  MTRand mtrand;

  int count = 0;
  bool cont;
  int bestsofar = 1;

  // if no population input, create random population of size Nrand,
  //    if population is input, ignore Nrand
  if(population.size() == 0) {
    for(int i = 0; i < Nrand; i++) {
      population.add(eci_in);
      population[i].randomize(Nmin, mtrand);
    }
  }

  do {
    cout << endl << "---------------------------" << endl;
    // make random start with N clusters turned on
    count++;

    eci_min_A = population[count - 1];

    fit(corr, DFT_nrg, eci_min_A, singular);

    cout << eci_min_A.get_bit_string() << "            init: 0" << " Nclust: " << eci_min_A.get_Nclust_on() << " cv: " << eci_min_A.get_cv() << " rms: " << eci_min_A.get_rms() << endl;

    eci_min_A = direct_min(Nmin, Nmax, eci_min_A, corr, DFT_nrg, 1);

    cout << eci_min_A.get_bit_string() << "           final: " << count << " Nclust: " << eci_min_A.get_Nclust_on() << " cv: " << eci_min_A.get_cv() << " rms: " << eci_min_A.get_rms() << endl;

    population[count - 1] = eci_min_A;

    if(eci_min_B.get_cv() == -1 || eci_min_A.get_cv() < eci_min_B.get_cv()) {
      eci_min_B = eci_min_A;
      bestsofar = count;
    }

    cout << eci_min_B.get_bit_string() << "       bestsofar: " << bestsofar << " Nclust: " << eci_min_B.get_Nclust_on() << " cv: " << eci_min_B.get_cv() << " rms: " << eci_min_B.get_rms() << endl;

  }
  while(count < population.size());
  cout << endl << eci_min_B.get_bit_string() << "       BestSoFar: " << bestsofar << " Nclust: " << eci_min_B.get_Nclust_on() << " cv: " << eci_min_B.get_cv() << " rms: " << eci_min_B.get_rms() << endl;

}

void calc_dfsmin_eci(int Nrand, int Nstop, int Nmin, int Nmax, string energy_filename, string eci_in_filename, string corr_in_filename, BP_Vec<ECISet> &population) {
  Correlation corr(corr_in_filename);
  EnergySet DFT_nrg(energy_filename);
  ECISet eci_in(eci_in_filename);

  // set fix so only up to triplets are used
  /*
  if( population.size() == 0)
  {
  	for( int j=0; j<eci_in.size(); j++)
  	{
  		if( j == 2)
  			eci_in.fix(j,0);

  		if( eci_in.get_size(j) > 3)
  		{
  			cout << " fix: " << j << " to " << 0 << endl;
  			eci_in.fix(j,0);
  		}
  	}
  }
  else
  {
  	for( int k=0; k<population.size(); k++)
  		for( int j=0; j<population[k].size(); j++)
  		{
  			if( j == 2)
  				population[k].fix(j,0);

  			if( population[k].get_size(j) > 3)
  			{
  				if( k == 0) cout << " fix: " << j << " to " << 0 << endl;
  				population[k].fix(j,0);
  			}
  		}
  }
  */

  // set fix so only clusters less than some size are used
  //{
  //	eci_in.set_Rmax_fix(6.0);
  //
  //}


  double cv_score;
  double rms_score;
  bool singular;

  ECISet eci_min_A = eci_in;
  ECISet eci_min_B = eci_in;

  MTRand mtrand;

  if(TEST)
    mtrand.seed(1);

  int count = 0;
  bool cont;
  int bestsofar = 1;

  if(population.size() == 0) {
    for(int i = 0; i < Nrand; i++) {
      population.add(eci_in);
      population[i].randomize(Nmin, mtrand);
    }
  }

  do {
    cout << endl << "---------------------------" << endl;
    // make random start with N clusters turned on
    count++;

    eci_min_A = population[count - 1];

    fit(corr, DFT_nrg, eci_min_A, singular);

    cout << eci_min_A.get_bit_string() << "            init: 0" << " Nclust: " << eci_min_A.get_Nclust_on() << " cv: " << eci_min_A.get_cv() << " rms: " << eci_min_A.get_rms() << endl;

    eci_min_A = dfs_min(Nstop, Nmin, Nmax, eci_min_A, corr, DFT_nrg, 1);

    cout << eci_min_A.get_bit_string() << "           final: " << count << " Nclust: " << eci_min_A.get_Nclust_on() << " cv: " << eci_min_A.get_cv() << " rms: " << eci_min_A.get_rms() << endl;

    population[count - 1] = eci_min_A;

    if(eci_min_B.get_cv() == -1 || eci_min_A.get_cv() < eci_min_B.get_cv()) {
      eci_min_B = eci_min_A;
      bestsofar = count;
    }

    cout << eci_min_B.get_bit_string() << "       bestsofar: " << bestsofar << " Nclust: " << eci_min_B.get_Nclust_on() << " cv: " << eci_min_B.get_cv() << " rms: " << eci_min_B.get_rms() << endl;

  }
  while(count < population.size());

  cout << endl << eci_min_B.get_bit_string() << "       BestSoFar: " << bestsofar << " Nclust: " << eci_min_B.get_Nclust_on() << " cv: " << eci_min_B.get_cv() << " rms: " << eci_min_B.get_rms() << endl;

}

void calc_ga_eci(int Npopulation, int Nmin, int Nmax, int Nchildren, int Nmutations, string energy_filename, string eci_in_filename, string corr_in_filename, BP_Vec<ECISet> &population) {
  Correlation corr(corr_in_filename);
  EnergySet DFT_nrg(energy_filename);
  ECISet eci_in(eci_in_filename);
  MTRand mtrand;
  bool singular;
  int best;

  // set fix so that only up to triplets are used
  /*
  if( population.size() == 0)
  {
  	for( int j=0; j<eci_in.size(); j++)
  	{
  		if( j == 2)
  			eci_in.fix(j,0);

  		if( eci_in.get_size(j) > 3)
  		{
  			cout << " fix: " << j << " to " << 0 << endl;
  			eci_in.fix(j,0);
  		}
  	}
  }
  else
  {
  	for( int k=0; k<population.size(); k++)
  		for( int j=0; j<population[k].size(); j++)
  		{
  			if( j == 2)
  				population[k].fix(j,0);

  			if( population[k].get_size(j) > 3)
  			{
  				if( k == 0) cout << " fix: " << j << " to " << 0 << endl;
  				population[k].fix(j,0);
  			}
  		}
  }
  */

  // create initial random population
  //cout << "creating initial random population" << endl;
  //BP_Vec<ECISet> population;
  if(population.size() == 0) {
    for(int i = 0; i < Npopulation; i++) {
      population.add(eci_in);
      population[i].randomize(Nmin, mtrand);
    }
  }

  cout << endl << "Initial Population: " << endl;
  for(int i = 0; i < population.size(); i++) {
    fit(corr, DFT_nrg, population[i], singular);
    cout << population[i].get_bit_string() << "    i: " << i << " Nclust: " << population[i].get_Nclust_on() << "  cv: " << population[i].get_cv() << " rms: " << population[i].get_rms() << endl;
  }
  cout << endl;

  // create and run genetic algorithm
  //cout << "create gac" << endl;
  GeneticAlgorithm ga(population, Nmin, Nmax, Nchildren, Nmutations);
  //cout << "run ga" << endl;
  ga.run(mtrand, corr, DFT_nrg, 1);
  population = ga.get_population();

  // output results
  cout << endl << "Final Population: " << endl;
  for(int i = 0; i < ga.size(); i++) {
    cout << ga[i].get_bit_string() << "    i: " << i << " Nclust: " << ga[i].get_Nclust_on() << "  cv: " << ga[i].get_cv() << " rms: " << ga[i].get_rms() << endl;
  }

  best = ga.get_best_index();

  cout << endl << ga[best].get_bit_string() << "   BestSoFar: " << best << " Nclust: " << ga[best].get_Nclust_on() << "  cv: " << ga[best].get_cv() << " rms: " << ga[best].get_rms() << endl;

}

void calc_ga_dir_eci(int Npopulation, int Nmin, int Nmax, int Nchildren, int Nmutations, string energy_filename, string eci_in_filename, string corr_in_filename, BP_Vec<ECISet> &population) {
  Correlation corr(corr_in_filename);
  EnergySet DFT_nrg(energy_filename);
  ECISet eci_in(eci_in_filename);
  MTRand mtrand;
  bool singular;
  int best;
  bool add_ok;

  // create initial random population
  if(population.size() == 0) {
    cout << "Creating initial random population" << endl;
    for(int i = 0; i < Npopulation; i++) {
      do {
        eci_in.randomize(Nmin, mtrand);
        //eci_in = direct_min( Nmin, Nmax, eci_in, corr, DFT_nrg, 1);
        add_ok = add_once(population, eci_in);
      }
      while(!add_ok);
      //population.add(eci_in);
      cout << population[i].get_bit_string() << "    Added" << endl;
    }
  }

  cout << endl << "Initial Population: " << endl;
  for(int i = 0; i < population.size(); i++) {
    fit(corr, DFT_nrg, population[i], singular);
    cout << population[i].get_bit_string() << "    i: " << i << " Nclust: " << population[i].get_Nclust_on() << "  cv: " << population[i].get_cv() << " rms: " << population[i].get_rms() << endl;
  }
  cout << endl;

  // create and run genetic algorithm
  //cout << "create gac" << endl;
  GeneticAlgorithm ga(population, Nmin, Nmax, Nchildren, Nmutations);
  //cout << "run ga" << endl;
  ga.run_dir(mtrand, corr, DFT_nrg, 1);
  population = ga.get_population();

  // output results
  cout << endl << "Final Population: " << endl;
  for(int i = 0; i < ga.size(); i++) {
    cout << ga[i].get_bit_string() << "    i: " << i << " Nclust: " << ga[i].get_Nclust_on() << "  cv: " << ga[i].get_cv() << " rms: " << ga[i].get_rms() << endl;
  }

  best = ga.get_best_index();

  cout << endl << ga[best].get_bit_string() << "   BestSoFar: " << best << " Nclust: " << ga[best].get_Nclust_on() << "  cv: " << ga[best].get_cv() << " rms: " << ga[best].get_rms() << endl;

}

void calc_ga_dfs_eci(int Npopulation, int Nmin, int Nmax, int Nchildren, int Nmutations, int Nstop, string energy_filename, string eci_in_filename, string corr_in_filename, BP_Vec<ECISet> &population) {
  Correlation corr(corr_in_filename);
  EnergySet DFT_nrg(energy_filename);
  ECISet eci_in(eci_in_filename);
  MTRand mtrand;
  bool singular;
  int best;

  // create initial random population
  //cout << "creating initial random population" << endl;
  if(population.size() == 0) {
    for(int i = 0; i < Npopulation; i++) {
      population.add(eci_in);
      population[i].randomize(Nmin, mtrand);
    }
  }

  cout << endl << "Initial Population: " << endl;
  for(int i = 0; i < population.size(); i++) {
    fit(corr, DFT_nrg, population[i], singular);
    cout << population[i].get_bit_string() << "    i: " << i << " Nclust: " << population[i].get_Nclust_on() << "  cv: " << population[i].get_cv() << " rms: " << population[i].get_rms() << endl;
  }
  cout << endl;

  // create and run genetic algorithm
  //cout << "create gac" << endl;
  GeneticAlgorithm ga(population, Nmin, Nmax, Nchildren, Nmutations);
  //cout << "run ga" << endl;
  ga.run_dfs(Nstop, mtrand, corr, DFT_nrg, 1);
  population = ga.get_population();

  // output results
  cout << endl << "Final Population: " << endl;
  for(int i = 0; i < ga.size(); i++) {
    cout << ga[i].get_bit_string() << "    i: " << i << " Nclust: " << ga[i].get_Nclust_on() << "  cv: " << ga[i].get_cv() << " rms: " << ga[i].get_rms() << endl;
  }

  best = ga.get_best_index();

  cout << endl << ga[best].get_bit_string() << "   BestSoFar: " << best << " Nclust: " << ga[best].get_Nclust_on() << "  cv: " << ga[best].get_cv() << " rms: " << ga[best].get_rms() << endl;

}

void weight_nrg(double A, double B, double kT, string energy_filename, string out_filename) {
  EnergySet nrg_set(energy_filename);

  nrg_set.weight_nrg(A, B, kT);

  nrg_set.write(out_filename);
}

void weight_EMin(double A, double B, double kT, string energy_filename, string out_filename) {
  EnergySet nrg_set(energy_filename);
  double minEn = nrg_set.find_energy_min();
  nrg_set.weight_ERef(A, B, kT, minEn);
  nrg_set.write(out_filename);
}

void weight_ERef(double A, double B, double ERef, double kT, string energy_filename, string out_filename) {
  EnergySet nrg_set(energy_filename);
  nrg_set.weight_ERef(A, B, kT, ERef);
  nrg_set.write(out_filename);
}

////----------------------------
/// test functions

//eci_search -nrg_diff   tot_corr.index.A    AB_structs.energy   eci.in.tot   AB_structs.corr.in.tot   A_structs.energy  A_structs.eci.in  A_structs.corr.in
void calc_nrg_diff(string A_corr_index_filename, string AB_structs_energy_filename, string eci_in_tot_filename, string AB_structs_corr_tot_filename,
                   string A_structs_energy_filename, string A_structs_eci_in_filename, string A_structs_corr_in_filename) {
  unsigned long int i, j;
  Correlation AB_structs_corr_tot(AB_structs_corr_tot_filename);
  EnergySet AB_structs_energy(AB_structs_energy_filename);
  ECISet eci_in_tot(eci_in_tot_filename);

  bool singular;
  Correlation A_structs_corr_in(A_structs_corr_in_filename);
  EnergySet A_structs_energy(A_structs_energy_filename);
  ECISet A_structs_eci_in(A_structs_eci_in_filename);
  fit(A_structs_corr_in,  A_structs_energy, A_structs_eci_in, singular);

  A_structs_eci_in.write_ECIout("A_structs.eci.out", A_structs_energy);

  // read  tot_corr.index.A file, store in A_corr_index_list, and create AB_corr_index_list
  BP_Parse index_file(A_corr_index_filename);
  BP_Vec<int> A_corr_index_list;
  BP_Vec<int> AB_corr_index_list;
  BP_Vec<int> i_list;
  do {
    i_list = index_file.getline_int();
    if(i_list.size() == 1) {
      A_corr_index_list.add(i_list[0]);
    }
  }
  while(index_file.eof() == false);
  //cout << "A_corr_index_list: " << A_corr_index_list << endl;

  for(i = 0; i < eci_in_tot.size(); i++) {
    if(!A_corr_index_list.find_first(i, j)) {
      AB_corr_index_list.add(i);
    }
  }
  //cout << "AB_corr_index_list: " << AB_corr_index_list << endl;

  A_structs_eci_in.write_ECIout(A_corr_index_list, "A_structs.eci.out.part", A_structs_energy);

  // calculate and write the energy_diff: AB_DFT_nrg_diff = AB_DFT_nrg - CLEX(AB_strucs.A_corr, A_eci);
  EnergySet AB_structs_A_nrg_diff = AB_structs_energy;
  AB_structs_A_nrg_diff.calc_nrg_diff(A_corr_index_list, AB_structs_corr_tot, A_structs_eci_in);
  AB_structs_A_nrg_diff.write("AB_structs.energy.diff");

  // get subset of AB_structs.corr.in.tot corresponding to AB_strucs.corr.in.AB (AB clusters)
  Correlation AB_structs_corr_AB = AB_structs_corr_tot;
  AB_structs_corr_AB.cluster_subset(AB_corr_index_list);
  AB_structs_corr_AB.write("AB_structs.corr.in.AB");

  Correlation AB_structs_corr_A = AB_structs_corr_tot;
  AB_structs_corr_A.cluster_subset(A_corr_index_list);
  AB_structs_corr_A.write("AB_structs.corr.in.A");

  // write eci.in.AB
  eci_in_tot.write_ECIin(AB_corr_index_list, "eci.in.AB");
}

////----------------------------
/// sub functions
ECISet direct_min(int Nmin, int Nmax, const ECISet &eci_start, Correlation &corr, EnergySet &nrg_set, bool print_steps) {
  ECISet eci_min_A = eci_start;
  ECISet eci_in;
  bool cont;
  bool singular;
  int Nchoice;
  double last_cv;
  int last_flip;
  // make sure eci_min_A cv score is known
  fit(corr, nrg_set, eci_min_A, singular);

  // minimize
  int j_count = 1;
  do {
    // initialize for minimization step
    cont = 0;
    eci_in = eci_min_A;
    Nchoice = 0;
    last_cv = eci_min_A.get_cv();
    last_flip = -1;

    // try toggling each cluster on/off
    for(int i = 0; i < eci_in.size(); i++) {
      eci_in.toggle_clust(i);

      // check that the eci_set meets fix and Nmin/Nmax criteria
      if(eci_in.fix_ok())
        if(eci_in.get_Nclust_on() >= Nmin && eci_in.get_Nclust_on() <= Nmax) {
          // find fit/cv score
          fit(corr, nrg_set, eci_in, singular);

          if(eci_in.get_cv() < last_cv) {
            Nchoice++;
          }

          // keep the steepest descent option
          if(eci_in.get_cv() < eci_min_A.get_cv()) {

            cont = 1;
            last_flip = i;
            eci_min_A = eci_in;
          }
        }

      eci_in.toggle_clust(i);
    }


    if(print_steps) cout << eci_min_A.get_bit_string() << "      minimizing: " << j_count << " Nclust: " << eci_min_A.get_Nclust_on() << " cv: " << eci_min_A.get_cv() << " rms: " << eci_min_A.get_rms() << " Nchoice: " << Nchoice << " last_flip: " << last_flip << endl;

    j_count++;


  }
  while(cont);

  return eci_min_A;
}

class FitData {
public:

  ECISet eci_base;
  Correlation corr;
  EnergySet nrg_set;
  BP_Vec<string> *bit_string_list;
  int Nmin;
  int Nmax;
  int start;
  int stop;
  int mode;

  int Nchoice;
  int last_flip;
  bool cont;

  BP_Vec<string> new_bit_string;
  BP_Vec<ECISet> eci_results;

  FitData(const ECISet &_eci, const Correlation &_corr, const EnergySet &_nrg_set, BP_Vec<string> *_bit_string_list,
          int _Nstop, int _Nmin, int _Nmax, int _start, int _stop, int _mode)
    : corr(_corr), eci_base(_eci), nrg_set(_nrg_set), bit_string_list(_bit_string_list), Nmin(_Nmin), Nmax(_Nmax), start(_start), stop(_stop) {
    if(_mode == 0) {
      mode = 0;
    }
    else if(_mode == 1) {
      mode = 1;
    }
    else {
      cout << "Error in FitData::FitData. _mode = " << _mode << endl;
      cout << "  Options are: 0 (direct), 1 (dfs)" << endl;
      exit(1);
    }

  }


  void min() {
    new_bit_string.clear();
    eci_results.clear();
    Nchoice = 0;
    bool singular;
    double last_cv = eci_base.get_cv();

    if(mode == 0) {
      cont = false;
      eci_results.add(eci_base);
    }

    // find all changes that reduce the cv score, and add to RVG_tree with rate = (init_cv - curr_cv)
    for(int i = start; i < stop; i++) {
      eci_base.toggle_clust(i);

      if(eci_base.fix_ok()) {
        if(eci_base.get_Nclust_on() >= Nmin && eci_base.get_Nclust_on() <= Nmax) {
          if(!bit_string_list->contains(eci_base.get_bit_string())) {
            new_bit_string.add(eci_base.get_bit_string());
            fit(corr, nrg_set, eci_base, singular);

            if(mode == 0) {
              // direct
              if(eci_base.get_cv() < last_cv)
                Nchoice++;

              if(eci_base.get_cv() < eci_results[0].get_cv()) {
                cont = true;
                last_flip = i;
                eci_results[0] = eci_base;
              }
            }
            else {
              // dfs

              if(eci_base.get_cv() < last_cv) {
                Nchoice++;
                eci_results.add(eci_base);
              }
            }
          }
        }
      }

      eci_in.toggle_clust(i);

    }
  };

  static void *min_threaded(void *_arg) {
    ((FitData *) _arg)->min();
    return NULL;
  }


};


ECISet dfs_min(long int Nstop, int Nmin, int Nmax, const ECISet &eci_start, Correlation &corr, EnergySet &nrg_set, bool print_steps) {
  // search for the minimum by depth-first search
  // minimize
  ECISet eci_min_A = eci_start;
  ECISet eci_min_B;
  ECISet eci_in;
  //BP_RVG_tree<ECISet> generator;
  BP_GVec<ECISet> eci_list;
  BP_Vec<string> bit_string_list;
  int count_since_last_best = 0;
  int Nchoice;
  int bestsofar;
  bool singular;
  BP_GVec_Member<ECISet> *s_member;

  // For threading:
  BP_Vec<FitData> fit_data;

  if(THREADS > 1) {
    for(int i = 0; i < THREADS; i++)
      fit_data.add(FitData(eci_min_A, corr, nrg_set, &bit_string_list, Nmin, Nmax, i * (eci_min_A.size() / THREADS), i == THREADS - 1 ? eci_min_A.size() : eci_min_A.size() / THREADS, 1));

  }


  // make sure eci_min_A cv score is known
  fit(corr, nrg_set, eci_min_A, singular);
  eci_min_B = eci_min_A;

  // depth-first search
  int dfs_count = 0;
  do {
    //cout << endl << "  -+++-----------------------" << endl;


    // !!!Multi-thread this section!!!
    if(THREADS == 1) {
      // find all changes that reduce the cv score, and add to RVG_tree with rate = (init_cv - curr_cv)
      eci_in = eci_min_A;
      Nchoice = 0;
      for(int i = 0; i < eci_in.size(); i++) {
        eci_in.toggle_clust(i);

        if(eci_in.fix_ok())
          if(eci_in.get_Nclust_on() >= Nmin && eci_in.get_Nclust_on() <= Nmax)
            if(add_once(bit_string_list, eci_in.get_bit_string())) {
              fit(corr, nrg_set, eci_in, singular);

              if(eci_in.get_cv() < eci_min_A.get_cv()) {
                Nchoice++;
                //generator.add( eci_list.add(eci_in), exp(1.0/eci_in.get_cv()));
                eci_list.add(eci_in);
              }
            }

        eci_in.toggle_clust(i);

      }
    }
    else {
      for(int i = 0; i < THREADS; i++)
        fitpool.add_work(FitData::min_threaded, (void *) & (fit_data[i]));
      fitpool.finish();

      for(int i = 0; i < THREADS; i++) {
        // eci_list.add results
        // bit_string_list.add new_bit_string
        // Nchoice sum
      }
    }



    if(eci_list.size() > 0) {
      // pick a bitstring from the RVG_tree
      //s_member = generator.pick(mtrand.randExc());

      // pick the bitstring with lowest cv score
      s_member = get_best_member(eci_list);

      eci_min_A = eci_list[ s_member];
      eci_list.remove(s_member);

      if(eci_min_B.get_cv() == -1 || eci_min_A.get_cv() < eci_min_B.get_cv()) {
        eci_min_B = eci_min_A;
        bestsofar = dfs_count;
        count_since_last_best = 0;
      }

      if(print_steps) cout << eci_min_A.get_bit_string() << "      DFSchoice: " << dfs_count << " Nclust: " << eci_min_A.get_Nclust_on() << " cv: " << eci_min_A.get_cv() << " rms: " << eci_min_A.get_rms() << " Nchoice: " << Nchoice << " listsize: " << eci_list.size() << " Count: " << count_since_last_best << "/" << Nstop << "  bestofDFS: " << bestsofar << "  best_cv: " << eci_min_B.get_cv() << endl;

      dfs_count++;
    }

    count_since_last_best++;

  }
  while((eci_list.size() > 0) && (count_since_last_best < Nstop || Nstop == 0));

  return eci_min_B;
}

////----------------------------
/// eci fit
void set_correlation_matrix(QXu::Array &A, Correlation &corr, EnergySet &nrg_set, ECISet &eci_in) {
  int in_i, in_j;

  in_i = 0;
  for(int i = 0; i < nrg_set.size(); i++) {
    if(nrg_set.get_weight(i) != 0) {
      in_j = 0;
      for(int j = 0; j < corr[i].size(); j++) {
        if(eci_in.get_weight(j) != 0) {
          A.set_zelem(in_i, in_j, nrg_set.get_weight(i)*corr[i][j]);
          in_j++;
        }
      }
      in_i++;
    }
  }

  //cout << " in_i: " << in_i << " in_j: " << in_j << endl;
}

void set_energy_vector(QXu::Array &E_vec, EnergySet &nrg_set) {
  int in_i = 0;
  for(int i = 0; i < nrg_set.size(); i++)
    if(nrg_set.get_weight(i) != 0) {
      E_vec.set_zelem(in_i, 0, nrg_set.get_weight(i)*nrg_set.get_Ef(i));
      in_i++;
    }
};

bool check_if_singular(QXu::Array &W) {
  for(int i = 0; i < W.num_row(); i++)
    if(W.zelem(i, 0) < 1.0e-8)	// CONSTANT
      return true;
  return false;
}

void fit_Array(Correlation &corr, EnergySet &nrg_set, ECISet &eci_in, bool &singular) {
  //cout << "begin fit()" << endl;

  // Least Squares Fit of the energies given correlations & calculation of Leave One Out cv score
  //
  // Input:
  //   energies in 'nrg_set'
  //   correlations in 'corr'
  //   eci to fit in eci_in
  //
  // Output:
  //   eci value in eci_in
  //   cv_score
  //   rms_score
  //   singular if(check_singular), otherwise assumed not to be singular
  bool check_singular = 0;


  //cout << "begin fit()" << endl;
  double inf = 1.0e20;
  double cv_score;
  double rms_score;
  int Nstruct = nrg_set.get_Nstruct_on();
  eci_in.set_Nstruct(Nstruct);
  int Nclust = eci_in.get_Nclust_on();

  //cout << "fit Nstruct: " << Nstruct << "  Nclust: " << Nclust << endl;

  QXu::Array corr_matrix(Nstruct, Nclust);
  QXu::Array E_vec(Nstruct, 1);
  int Nstruct_in_cv;	//svdls the number of data used to calculate cv score

  //cout << " set correlation matrix" << endl;
  set_correlation_matrix(corr_matrix, corr, nrg_set, eci_in);

  //cout << " set energy vector" << endl;
  set_energy_vector(E_vec, nrg_set);


  if(check_singular) {
    QXu::Array U(Nstruct, Nclust);
    QXu::Array W(Nclust, 1);
    QXu::Array V(Nclust, Nclust);

    // do svd
    //cout << " do svd" << endl;
    corr_matrix.svd(U, W, V);

    //cout << " check if singular" << endl;
    singular = check_if_singular(W);

    if(singular) {
      cv_score = inf;
      rms_score = inf;
      return;
    }
  }
  else {
    singular = false;
  }

  QXu::Array covar(Nclust, Nclust);
  QXu::Array ECI(Nclust, 1);

  // Calculate rms fit and Leave One Out cv score
  //cout << " calc svdls" << endl;
  ECI = corr_matrix.svdls(E_vec, covar, rms_score, cv_score, Nstruct_in_cv);

  //for(int i=0;i<ECI.num_row(); i++)
  //	cout << "i: " << i << "  eci: " << ECI.zelem(i,0) << endl;

  //cout << " set eci values" << endl;
  eci_in.set_values(ECI);
  eci_in.set_cv(cv_score);
  eci_in.set_rms(rms_score);
  //cout << "finish fit()" << endl;

}

void set_correlation_matrix(MatrixXd &A, const Correlation &corr, const EnergySet &nrg_set, const ECISet &eci_in) {
  // set A to be the correlation matrix, including weights, and only the rows and columns being fit
  // assumes A is already the right size

  int i, j, in_i, in_j;

  in_i = 0;
  for(i = 0; i < nrg_set.size(); i++) {
    if(nrg_set.get_weight(i) != 0) {
      in_j = 0;
      for(j = 0; j < corr[i].size(); j++) {
        if(eci_in.get_weight(j) != 0) {
          A(in_i, in_j) = nrg_set.get_weight(i) * corr[i][j];
          in_j++;
        }
      }
      in_i++;
    }
  }

  //cout << " in_i: " << in_i << " in_j: " << in_j << endl;
}

bool check_if_singular(const VectorXd &S) {
  for(unsigned long int i = 0; i < S.size(); i++)
    if(fabs(S(i)) < 1.0e-8)	// CONSTANT
      return true;
  return false;
}

void fit(Correlation &corr, EnergySet &nrg_set, ECISet &eci_in, bool &singular) {
  // Least Squares Fit of the energies given correlations & calculation of Leave One Out cv score
  // Uses Eigen functions
  //
  // Input:
  //   energies in 'nrg_set'
  //   correlations in 'corr'
  //   eci to fit in eci_in
  //
  // Output:
  //   eci value in eci_in
  //   cv_score
  //   rms_score
  //   singular if(check_singular), otherwise assumed not to be singular
  bool check_singular = 0;

  //cout << "begin fit()" << endl;
  double inf = 1.0e20;
  int Nstruct = nrg_set.get_Nstruct_on();
  eci_in.set_Nstruct(Nstruct);
  int Nclust = eci_in.get_Nclust_on();

  //cout << "fit Nstruct: " << Nstruct << "  Nclust: " << Nclust << endl;
  MatrixXd corr_matrix(Nstruct, Nclust);

  //cout << " set correlation matrix" << endl;
  set_correlation_matrix(corr_matrix, corr, nrg_set, eci_in);

  //cout << " init svd" << endl;
  JacobiSVD<MatrixXd> svd(corr_matrix, ComputeThinU | ComputeThinV);

  if(check_singular) {
    VectorXd S;

    // do svd
    //cout << " do svd" << endl;
    S = svd.singularValues();

    //cout << " check if singular" << endl;
    singular = check_if_singular(S);

    if(singular) {
      eci_in.set_cv(inf);
      eci_in.set_rms(inf);
      return;
    }

  }
  else {
    singular = false;
  }

  VectorXd E_vec(Nstruct);
  VectorXd ECI(Nclust);
  VectorXd Err(Nstruct);
  double cv_score;
  double rms_score;
  unsigned long int i, ii;

  //cout << " set energy vector" << endl;
  if(!nrg_set.E_vec_is_ready()) {
    nrg_set.set_E_vec();
    corr.cv_a_reset();
  }

  // Calculate rms fit
  //cout << " rms fit" << endl;
  ECI = svd.solve(nrg_set.get_E_vec());

  //cout << "ECI: " << ECI << endl;

  // Calculate residuals
  //cout << " calc residuals" << endl;
  Err = corr_matrix * ECI - nrg_set.get_E_vec();

  // Calculate Leave One Out cv score
  // LOOCV = (1.0/Nnrg)*sum_i{ (e_i / 1 - X_i*((X^T*X)^-1)*X_i^T)^2 }
  // e_i is residual
  // X is corr_matrix
  // X_i is row i of X

  // -- compute cv_a_i = X_i*((X^T*X)^-1)*X_i^T ahead of time
  // this only changes if the structures you are fitting to change, or their weights, (?? or corr fitting to ?? Then need to recalc each time)
  //cout << " set cv_a" << endl;
  //if(!corr.cv_a_is_ready())
  {
    // trial ///////////////////////////////
    // include all corr ???

    //int Nstruct = nrg_set.get_Nstruct_on();
    //int Nclust_all = eci_in.size();

    //cout << "fit Nstruct: " << Nstruct << "  Nclust: " << Nclust << endl;
    //MatrixXd corr_matrix_all(Nstruct, Nclust_all);
    //ECISet all_eci = eci_in;
    //for( i=0; i<all_eci.size(); i++)
    //	all_eci.set_clust_on(i);

    //cout << " set correlation matrix" << endl;
    //set_correlation_matrix(corr_matrix_all, corr, nrg_set, all_eci);

    // this gives different results than using the corr_matrix with only corr being fit
    //corr.set_cv_a(corr_matrix_all, nrg_set);

    // trial ////////////////////////////////

    corr.set_cv_a(corr_matrix, nrg_set);
  }

  //cout << " calc rms and cv" << endl;
  rms_score = 0.0;
  cv_score = 0.0;
  ii = 0;
  for(i = 0; i < nrg_set.size(); i++)
    if(nrg_set.get_weight(i) != 0) {
      rms_score += Err(ii) * Err(ii);
      cv_score += sqr(Err(ii) / (1.0 - corr.get_cv_a(ii)));
      //cout << "cv_a " << ii << " :: " << corr.get_cv_a(ii) << endl;
      ii++;
    }

  rms_score = sqrt(rms_score / Nstruct);
  cv_score = sqrt(cv_score / Nstruct);

  //cout << " set eci values" << endl;
  eci_in.set_values(ECI);
  eci_in.set_cv(cv_score);
  eci_in.set_rms(rms_score);

  //cout << "finish fit_Eigen()" << endl;
  return;
}

////----------------------------
/// compressive sensing functions
ECISet calc_FPC_eci(double mu, double prec_shrink, double prec_ECI, const ECISet &eci_set, const Correlation &corr, const EnergySet &nrg_set, bool print_steps) {
  cout << "begin calc_FPC_eci()" << endl;

  // method:
  //		Start: ECI_0 = 0 vector
  //		Then:
  //			ECI_i+1 = shrink( ECI_i - tau*g_i, mu*tau)
  //		Where:
  //			g_i = Corr_transpose*(Corr*ECI_i - E) = M1*ECI_i - V1
  //			shrink(y,a) = sign(y)*max( fabs(y) - a, 0)
  //			tau = min( 1.999, -1.665*Corr.rows()/Corr.cols() + 2.665)
  //		Stop when:
  //			max(g)/mu - 1 < prec_shrink
  //		and
  //			2norm( ECI_i+1 - ECI_i)/2norm(ECI_i) < prec_ECI


  unsigned long int i, j, k, ii, jj, kk;
  ECISet eci_out = eci_set;

  unsigned long int Nnrg = nrg_set.get_Nstruct_on();
  unsigned long int Neci = eci_set.size();

  MatrixXd C(Nnrg, Neci);					// Correlation matrix
  MatrixXd Cn(Nnrg, Neci);				// normalized correlation matrix so that max eigenvalue of Cn.transpose*Cn <= 1
  MatrixXd M1(Neci, Neci);				// Cn.transpose()*Cn;
  VectorXd E(Nnrg);						// Enegry vector
  VectorXd En(Nnrg);						// normalized energy vector
  VectorXd V1(Neci);						// Cn.transpose()*En
  VectorXd ECI = VectorXd::Zero(Neci);	// current solution

  double tau = std::min(1.999, std::max(1.0, -1.665 * (1.0 * Nnrg) / (1.0 * Neci) + 2.665));
  double rms;

  //cout << "Nnrg: " << Nnrg << endl;
  //cout << "Neci: " << Neci << endl;
  //cout << "tau: " << tau << endl;

  // set Correlation matrix
  //cout << "set C" << endl;
  ii = 0;
  for(i = 0; i < nrg_set.size(); i++)
    if(nrg_set.get_weight(i) != 0) {
      jj = 0;
      for(j = 0; j < eci_set.size(); j++) {
        C(ii, jj) = nrg_set.get_weight(i) * corr[i][j];		// include weight!?
        jj++;
      }
      ii++;
    }

  // set Energy vector
  //cout << "set E" << endl;
  ii = 0;
  for(i = 0; i < nrg_set.size(); i++)
    if(nrg_set.get_weight(i) != 0) {
      E(ii) = nrg_set.get_weight(i) * nrg_set.get_Ef(i);		// include weight!?
      ii++;
    }

  // normalize C and E, so that largest eigenvalue of C.transpose*C is <= 1
  M1 = C.transpose() * C;
  SelfAdjointEigenSolver<MatrixXd> eigensolver(M1);
  if(eigensolver.info() != Success) {
    cout << "SelfAdjointEigenSolver failed!" << endl;
    exit(1);
  };
  //cout << "The eigenvalues of M1 are:\n" << eigensolver.eigenvalues() << endl;

  double max_eigenvalue = eigensolver.eigenvalues().maxCoeff();
  double a1 = sqrt(1.1 * max_eigenvalue);
  //cout << "The max eigenvalue of M1 is:\n" << max_eigenvalue << endl;

  // set normalized C & E
  Cn = C / a1;
  En = E / a1;

  //cout << "Cn: " << Cn << endl;
  //cout << "En: " << En << endl;


  // set M1 = corr.transpose() * corr
  //cout << "set M1" << endl;
  M1 = Cn.transpose() * Cn;

  // set V1 = corr.transpose() * nrg
  //cout << "set V1" << endl;
  V1 = Cn.transpose() * En;
  //cout << "V1: " << V1 << endl;

  // do Fixed-Point continuation algorithm to find ECI
  FPC(M1, V1, ECI, mu, tau, prec_shrink, prec_ECI, print_steps);

  //cout << "Final ECI:" << endl;
  //cout << ECI << endl;

  eci_out.set_values_and_weights(ECI);
  eci_out.set_cv(1e20);
  rms = (C * ECI - E).norm() / sqrt(1.0 * Nnrg);
  eci_out.set_rms(rms);

  cout << "finish calc_FPC_eci()" << endl;
  return eci_out;


};

ECISet calc_BI_eci(double mu, double prec_shrink, double prec_ECI, const ECISet &eci_set, const Correlation &corr, const EnergySet &nrg_set, bool print_steps) {
  cout << "begin calc_BI_eci()" << endl;

  // method:
  //		Start: ECI_0 = 0 vector, F_0 = 0 vector
  //		Then:
  //			F_i+1 = E + F_i - Corr*ECI_i
  //			ECI_i+1 = FPC( Corr, E, mu, tau)
  //		Where:
  //			tau = min( 1.999, -1.665*Corr.rows()/Corr.cols() + 2.665) // same as FPC
  //		Stop when:
  //			2norm( ECI_i+1 - ECI_i)/2norm(ECI_i) < prec_ECI


  unsigned long int i, j, k, ii, jj, kk;
  ECISet eci_out = eci_set;

  unsigned long int Nnrg = nrg_set.get_Nstruct_on();
  unsigned long int Neci = eci_set.size();

  MatrixXd C(Nnrg, Neci);					// Correlation matrix
  MatrixXd Cn(Nnrg, Neci);				// normalized correlation matrix so that max eigenvalue of Cn.transpose*Cn <= 1
  MatrixXd M1(Neci, Neci);				// Cn.transpose()*Cn;
  VectorXd E(Nnrg);						// Enegry vector
  VectorXd En(Nnrg);						// normalized energy vector
  VectorXd ECI = VectorXd::Zero(Neci);	// current solution

  double tau = std::min(1.999, std::max(1.0, -1.665 * (1.0 * Nnrg) / (1.0 * Neci) + 2.665));
  double rms;

  //cout << "Nnrg: " << Nnrg << endl;
  //cout << "Neci: " << Neci << endl;
  //cout << "tau: " << tau << endl;

  // set Correlation matrix
  //cout << "set C" << endl;
  ii = 0;
  for(i = 0; i < nrg_set.size(); i++)
    if(nrg_set.get_weight(i) != 0) {
      jj = 0;
      for(j = 0; j < eci_set.size(); j++) {
        C(ii, jj) = nrg_set.get_weight(i) * corr[i][j];		// include weight!?
        jj++;
      }
      ii++;
    }

  // set Energy vector
  //cout << "set E" << endl;
  ii = 0;
  for(i = 0; i < nrg_set.size(); i++)
    if(nrg_set.get_weight(i) != 0) {
      E(ii) = nrg_set.get_weight(i) * nrg_set.get_Ef(i);		// include weight!?
      ii++;
    }

  // normalize C and E, so that largest eigenvalue of C.transpose*C is <= 1
  M1 = C.transpose() * C;
  SelfAdjointEigenSolver<MatrixXd> eigensolver(M1);
  if(eigensolver.info() != Success) {
    cout << "SelfAdjointEigenSolver failed!" << endl;
    exit(1);
  };
  //cout << "The eigenvalues of M1 are:\n" << eigensolver.eigenvalues() << endl;

  double max_eigenvalue = eigensolver.eigenvalues().maxCoeff();
  double a1 = sqrt(1.1 * max_eigenvalue);
  //cout << "The max eigenvalue of M1 is:\n" << max_eigenvalue << endl;

  // set normalized C & E
  Cn = C / a1;
  En = E / a1;

  //cout << "Cn: " << Cn << endl;
  //cout << "En: " << En << endl;


  // set M1 = corr.transpose() * corr
  //cout << "set M1" << endl;
  M1 = Cn.transpose() * Cn;

  // set V1 = corr.transpose() * nrg
  //cout << "set V1" << endl;
  //V1 = Cn.transpose()*En;
  //cout << "V1: " << V1 << endl;

  // do Bergman Iteration algorithm to find ECI
  BI(Cn, En, M1, ECI, mu, tau, prec_shrink, prec_ECI, print_steps);

  //cout << "Final ECI:" << endl;
  //cout << ECI << endl;

  eci_out.set_values_and_weights(ECI);
  eci_out.set_cv(1e20);
  rms = (C * ECI - E).norm() / sqrt(1.0 * Nnrg);
  eci_out.set_rms(rms);

  cout << "finish calc_FPC_eci()" << endl;
  return eci_out;


};

void BI(const MatrixXd &Cn, const VectorXd &En, const MatrixXd &M1, VectorXd &ECI, double mu, double tau, double prec_shrink, double prec_ECI, bool print_steps) {
  // Bergman iteration for fitting ECI that minimize L1 norm

  VectorXd ECI_i;
  VectorXd F = VectorXd::Zero(En.size());
  VectorXd V1(ECI.size());						// Cn.transpose()*En

  unsigned long int step = 0;
  bool cont = true;
  double dECI, rms;

  do {
    step++;

    ECI_i = ECI;
    F = En + F - Cn * ECI;
    V1 = Cn.transpose() * F;
    FPC(M1, V1, ECI, mu, tau, prec_shrink, prec_ECI, false);

    dECI = (ECI - ECI_i).norm() / ECI_i.norm();

    if(print_steps) {
      if(step == 1) {
        cout << "    Step " << step << endl;
      }
      else
        cout << "    Step " << step << "  dECI: " << dECI << endl;
    }

    if(dECI < prec_ECI)
      cont = false;

  }
  while(cont);

}

void FPC(const MatrixXd &M1, const VectorXd &V1, VectorXd &ECI, double mu, double tau, double prec_shrink, double prec_ECI, bool print_steps) {
  //  Fixed-Point continuation algorithm for fitting ECI that minimize L1 norm
  //
  //

  //cout << "begin FPC()" << endl;

  unsigned long int step = 0;
  bool cont = true;
  double ECI_mag, delta_mag;

  VectorXd G = VectorXd::Zero(ECI.size());		// gradient of 2norm

  do {
    step++;

    //cout << " :1" << endl;
    ECI_mag = ECI.norm();
    //cout << " :2" << endl;
    G = M1 * ECI - V1;
    //cout << " :3" << endl;
    delta_mag = shrink(ECI, G, mu, tau);

    if(print_steps) {
      if(step % 1000 == 0) {
        //rms = (C*ECI - E).norm()/sqrt(1.0*Nnrg);
        //cout << ECI << endl;
        //cout << "  Step " << step << "  rms: " << rms << "  shrink: " << G.maxCoeff()/mu - 1.0 << "  dECI: " << (delta_mag/ECI_mag) << endl;
        cout << "  Step " << step << "  shrink: " << G.maxCoeff() / mu - 1.0 << "  dECI: " << (delta_mag / ECI_mag) << endl;

      }
    }

    //cout << " :4" << endl;
    if(G.maxCoeff() / mu - 1.0 < prec_shrink) {
      if((delta_mag / ECI_mag) < prec_ECI)
        cont = false;
    }

    //BP_pause();

  }
  while(cont);

};

double shrink(VectorXd &ECI, const VectorXd &G, double mu, double tau) {
  //cout << "begin shrink()" << endl;
  //	perform:
  //		ECI_i+1 = shrink( ECI_i - tau*g_i, mu*tau)
  //		shrink(y,a) = sign(y)*max( fabs(y) - a, 0)
  //	return:
  //		2norm(ECI_i+1 - ECI_i)

  double sqr_sum = 0.0;
  double y, a, init, delta;
  unsigned long int i;

  a = mu * tau;
  //cout << "a: " << a << endl;
  for(i = 0; i < ECI.size(); i++) {
    init = ECI(i);
    y = ECI(i) - tau * G(i);
    ECI(i) = sign(y) * std::max<double>(fabs(y) - a, 0.0);
    sqr_sum += sqr(ECI(i) - init);

    //cout << init << " " << ECI(i) << "  y: " << y <<  " delta: " << ECI(i) - init << endl;

  }

  //cout << "finish shrink()" << endl;
  //BP_pause();
  return sqrt(sqr_sum);

};



////----------------------------
/// misc.
bool is_bitstring(string s) {

  for(int i = 0; i < s.size(); i++) {
    if(s[i] != '0' && s[i] != '1') {
      return false;
    }
  }
  return true;
};

BP_Vec<ECISet> read_bit_strings_file(ECISet eci_in, string bit_strings_filename) {
  //cout << "begin read_bit_strings_file()" << endl;
  BP_Parse file(bit_strings_filename);

  ECISet eci_A = eci_in;
  BP_Vec<ECISet> eci_list;
  BP_Vec<string> entry;
  // data
  do {
    entry = file.getline_string();
    if(entry.size() != 0) {
      // check if the first entry is a bit_string
      if(is_bitstring(entry[0])) {
        eci_A.set_bit_string(entry[0]);
        eci_list.add(eci_A);
      }

    }
  }
  while(file.eof() == false);


  //cout << "Population from input file: " << endl;
  //for( int i=0; i<eci_list.size(); i++)
  //{
  //	cout << eci_list[i].get_bit_string() << "  i: " << i << endl;
  //}

  //cout << "finish read_bit_strings_file()" << endl;
  return eci_list;
};

double max_cv(BP_GVec<ECISet> &population, BP_GVec_Member<ECISet> *&worst_parent) {
  double max = -1;
  for(int i = 0; i < population.size(); i++) {
    if(population[i].get_cv() > max) {
      worst_parent = population.member(i);
      max = population[i].get_cv();
    }
  }

  return max;
};

double min_cv(BP_GVec<ECISet> &population, BP_GVec_Member<ECISet> *&worst_parent) {
  double max = -1;
  double min = 1.0e20;
  for(int i = 0; i < population.size(); i++) {
    if(population[i].get_cv() < min) {
      //worst_parent = population.member(i);
      min = population[i].get_cv();
    }

    if(population[i].get_cv() > max) {
      worst_parent = population.member(i);
      max = population[i].get_cv();
    }
  }

  return min;
};

int get_best_index(const BP_Vec<ECISet> &population) {
  int best;
  double min = 1.0e20;
  for(int i = 0; i < population.size(); i++) {
    if(population[i].get_cv() < min) {
      best = i;
      min = population[i].get_cv();
    }
  }

  return best;
};

int get_best_index(const BP_GVec<ECISet> &population) {
  int best;
  double min = 1.0e20;
  for(int i = 0; i < population.size(); i++) {
    if(population[i].get_cv() < min) {
      best = i;
      min = population[i].get_cv();
    }
  }

  return best;
};

BP_GVec_Member<ECISet> *get_best_member(BP_GVec<ECISet> &population) {
  int best;
  double min = 1.0e20;
  for(int i = 0; i < population.size(); i++) {
    if(population[i].get_cv() < min) {
      best = i;
      min = population[i].get_cv();
    }
  }

  return population.member(best);
};

#endif // eci_search_classes_HH
