/*
 *  EnergySet.cc
 */

#ifndef EnergySet_CC
#define EnergySet_CC

#include "EnergySet.hh"
#include "Functions.hh"

Energy::Energy(BP::BP_Vec<string> s_list) {
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

}

Energy::Energy(const CASM::jsonParser &json) {
  from_json(*this, json);
}



Energy::Energy(double i0, double i1, BP::BP_Vec<double> i2, double i3, string i4) {
  Ef = i0;
  weight = i1;
  concentration = i2;
  dist_from_hull = i3;
  dist_from_hull_found = true;
  name = i4;
}


CASM::jsonParser &to_json(const Energy &nrg, CASM::jsonParser &json) {
  json.put_obj();
  json["Ef"] = nrg.Ef;
  json["weight"] = nrg.weight;
  json["conc"] = nrg.concentration;
  if(nrg.dist_from_hull_found)
    json["dist_from_hull"] = nrg.dist_from_hull;
  json["name"] = nrg.name;
  return json;
}


void from_json(Energy &nrg, const CASM::jsonParser &json) {
  json["Ef"].get(nrg.Ef);
  json["weight"].get(nrg.weight);
  json["conc"].get(nrg.concentration);
  nrg.dist_from_hull_found = json.get_if(nrg.dist_from_hull, "dist_from_hull");
  json["name"].get(nrg.name);
}



EnergySet::EnergySet() {
  name = "-";
  E_vec_ready = false;
  hull_found = 0;
  Nstruct_set = 0;
  set_Nstruct_on();
  set_E_vec();
  m_format = "text";
}

EnergySet::EnergySet(BP::BP_Vec<Energy> energy_list) {
  name = "-";
  E_vec_ready = false;
  hull_found = 0;
  Nstruct_set = 0;

  BP::BP_Vec< Energy>::operator=(energy_list);

  set_Nstruct_on();
  set_E_vec();
  m_format = "text";
}

EnergySet::EnergySet(string energy_filename) {

  m_format = get_format_from_ext(energy_filename);

  if(m_format == "text") {
    name = energy_filename;
    E_vec_ready = false;
    hull_found = 0;
    Nstruct_set = 0;
    BP::BP_Parse file(energy_filename);

    // header line
    file.getline();

    BP::BP_Vec<string> entry;
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
    set_E_vec();
  }
  else if(m_format == "json") {
    CASM::jsonParser json(energy_filename);
    from_json(json);
  }
  else {
    std::cout << "Unexpected format option for EnergySet constructor" << std::endl;
    std::cout << "  Expected 'text' or 'json', but received: " << m_format << std::endl;
    exit(1);
  }
}

std::string EnergySet::format() const {
  return m_format;
}

double EnergySet::calc_clex(const Correlation &corr, const ECISet &eci, unsigned long int i) {
  unsigned long int j;
  double Ef = 0;
  for(j = 0; j < corr[i].size(); j++) {
    Ef += corr[i][j] * eci.get_value(j);
  }
  return Ef;
}

void EnergySet::calc_clex(const Correlation &corr, const ECISet &eci) {
  unsigned long int i, j;
  for(i = 0; i < size(); i++) {
    (*this)[i].Ef = 0;
    for(j = 0; j < corr[i].size(); j++) {
      (*this)[i].Ef += corr[i][j] * eci.get_value(j);
      (*this)[i].dist_from_hull = 0.0;
    }
  }
}

void EnergySet::calc_nrg_diff(BP::BP_Vec<int> &index_list, const Correlation &corr, const ECISet &eci) {
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
}

void EnergySet::set(BP::BP_Vec<Energy> energy_list) {
  E_vec_ready = false;
  hull_found = 0;
  Nstruct_set = 0;

  BP::BP_Vec< Energy>::operator=(energy_list);

  set_Nstruct_on();
}

bool EnergySet::E_vec_is_ready() const {
  return E_vec_ready;
}

void EnergySet::set_E_vec() {
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
}

const Eigen::VectorXd &EnergySet::get_E_vec() const {
  return E_vec;
}

BP::BP_Vec< Energy> EnergySet::get(BP::BP_Vec<int> indices) const {
  // return a BP::BP_Vec< Energy> for the input indices
  //    if any indices are bad, return an empty BP::BP_Vec

  BP::BP_Vec< Energy> energy_list;
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

int EnergySet::get_Nstruct_on() const {
  return Nstruct_on;
}

void EnergySet::set_Nstruct_on() {
  Nstruct_on = 0;
  for(int i = 0; i < size(); i++) {
    if((*this)[i].weight != 0)
      Nstruct_on++;
  }
  Nstruct_set = 1;
}

unsigned long int EnergySet::size() const {
  return BP::BP_Vec<Energy>::size();
}

double EnergySet::get_weight(unsigned long int i) const {
  return (*this)[i].weight;
}

void EnergySet::set_weight(unsigned long int i, double w) {
  (*this)[i].weight = w;
  Nstruct_set = 0;
  set_Nstruct_on();
}

double EnergySet::get_Ef(unsigned long int i) const {
  return (*this)[i].Ef;
}

double EnergySet::get_concentration(unsigned long int i, unsigned long int j) const {
  return (*this)[i].concentration[j];
}

unsigned long int EnergySet::get_concentration_size() const {
  return (*this)[0].concentration.size();
}

double EnergySet::get_dist_from_hull(unsigned long int i) const {
  return (*this)[i].dist_from_hull;
}

std::string EnergySet::get_name() const {
  return name;
}

double EnergySet::find_energy_min() {
  double minVal = get_Ef(0);
  for(int i = 0; i < BP::BP_Vec<Energy>::size(); i++) {
    if(get_Ef(i) <= minVal) {
      minVal = get_Ef(i);
    }
  }
  //        std:cout<<minVal<<"\n";
  return minVal;
}

bool EnergySet::concentration_is_less(unsigned long int i1, unsigned long int i2) const {

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

}

void EnergySet::calc_hull(bool calc_dist_from_hull, double hulltol) {
  // create matrix of points -------------------------
  unsigned long int i, j;
  unsigned long int Nrows = 1 + (*this)[0].concentration.size();	// dimensions are nrg + conc's
  unsigned long int Ncols = size();
  Eigen::MatrixXd m(Nrows, Ncols); //energy, composition1, composition2...

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
  hull.set_Geo_tol(hulltol);
  hull.set_CH_tol(hulltol);
  hull.reset_points(m, true);
  //hull.write_equivalent_points(std::cout);
  hull.set_verbosity(0);
  Eigen::VectorXd bottom_vector = Eigen::VectorXd::Zero(Nrows);
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
}

bool EnergySet::is_hull_found() const {
  return hull_found;
}

BP::BP_Vec<Energy> EnergySet::get_hull_points() {
  return get(get_hull_indices());
}

BP::BP_Vec<int> EnergySet::get_hull_indices() {
  if(!hull_found) {
    calc_hull(false);
  }

  return hull_indices;
}

BP::BP_Vec<Energy> EnergySet::get_non_hull_points() {
  return get(get_non_hull_indices());
}

BP::BP_Vec<int> EnergySet::get_non_hull_indices() {
  if(!hull_found) {
    calc_hull(false);
  }

  return non_hull_indices;
}

void EnergySet::write_hull(string s, std::string format) {

  if(format == "default")
    format = m_format;

  if(format == "text") {
    BP::BP_Write file(rm_json_ext(s));
    file.newfile();
    write_hull(file.get_ostream());
  }
  else if(format == "json") {
    if(!hull_found) {
      calc_hull(false);
      if(!hull_found)
        return;
    }

    CASM::jsonParser json;
    json.put_array();
    BP::BP_Vec< BP::BP_Vec< int> > equiv_points = hull.get_equivalent_points();
    for(int i = 0; i < hull_indices.size(); i++) {
      json.push_back((*this)[hull_indices[i]]);

      // add names of equivalent configurations
      unsigned long int k;
      json[json.size() - 1]["equivalents"].put_array();
      for(int j = 0; j < equiv_points.size(); j++) {
        if(equiv_points[j].find_first(hull_indices[i], k)) {
          for(int ii = 0; ii < equiv_points[j].size(); ii++)
            if(ii != k) {
              json[json.size() - 1]["equivalents"].push_back((*this)[ equiv_points[j][ii]].name);
            }
        }
      }

    }

    json.write(json_ext(s));

  }
  else {
    std::cout << "Unexpected format option for EnergySet::write" << std::endl;
    std::cout << "  Expected 'text' or 'json', but received: " << format << std::endl;
    exit(1);
  }

}

void EnergySet::write_hull(ostream &sout) {
  if(!hull_found) {
    calc_hull(false);
    if(!hull_found)
      return;
  }

  unsigned long int i, j, k, ii;
  BP::BP_Vec< BP::BP_Vec< int> > equiv_points = hull.get_equivalent_points();

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

}

void EnergySet::write(string s, std::string format) {

  if(format == "default")
    format = m_format;

  if(format == "text") {
    BP::BP_Write file(rm_json_ext(s));
    file.newfile();
    write(file.get_ostream());
  }
  else if(format == "json") {
    CASM::jsonParser json;
    to_json(json).write(json_ext(s));
  }
  else {
    std::cout << "Unexpected format option for EnergySet::write" << std::endl;
    std::cout << "  Expected 'text' or 'json', but received: " << format << std::endl;
    exit(1);
  }

}

void EnergySet::write(ostream &sout) {
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
}

void EnergySet::plot(BP::BP_Plot &plot, string label, string s, int line_style, double line_width, int point_style, double point_size, bool point_face, int component1, int component2) {
  BP::BP_Vec<double> x_list;
  BP::BP_Vec<double> y_list;
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

    BP::BP_Vec<int> nborlist;
    Eigen::VectorXd v, nbor_v;

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

}

void EnergySet::write_below_hull(string s, BP::BP_Vec<Energy> energy_list, std::string format) {

  if(!hull_found) {
    calc_hull(false);
    if(!hull_found)
      return;
  }

  Eigen::VectorXd v;
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

  if(format == "default")
    format = m_format;

  if(format == "text") {
    BP::BP_Write below_hull(rm_json_ext(s));
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
  }
  else if(format == "json") {
    CASM::jsonParser json;
    json.put_array();
    for(i = 0; i < energy_list.size(); i++)
      if(energy_list[i].dist_from_hull < 0.0)
        json.push_back(energy_list[i]);
    json.write(json_ext(s));
  }
  else {
    std::cout << "Unexpected format option for EnergySet::write_below_hull" << std::endl;
    std::cout << "  Expected 'text' or 'json', but received: " << format << std::endl;
    exit(1);
  }

}

void EnergySet::weight_nrg(double A, double B, double kT) {
  E_vec_ready = false;

  for(int i = 0; i < size(); i++)
    (*this)[i].weight = A * exp(-(*this)[i].dist_from_hull / kT) + B;

}

void EnergySet::weight_ERef(double A, double B, double kT, double ERef) {
  E_vec_ready = false;

  for(int i = 0; i < size(); i++)
    (*this)[i].weight = (((*this)[i].Ef <= ERef) ? 1.0 : (A * exp(-((*this)[i].Ef - ERef) / kT) + B));

}

double EnergySet::calc_rms(EnergySet &that, bool include_weights, double dist_from_hull_cutoff, unsigned long int &count) {
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
          rms += BP::sqr(get_Ef(i) * get_weight(i) - that.get_Ef(i) * that.get_weight(i));
          count++;
        }
      }
      else {
        if(dist_from_hull_cutoff <= 0 || that.get_dist_from_hull(i) < dist_from_hull_cutoff) {
          rms += BP::sqr(get_Ef(i) - that.get_Ef(i));
          count++;
        }
      }
    }

  //cout << "count: " << count << endl;

  rms /= count;
  return sqrt(rms);

}

CASM::jsonParser &EnergySet::to_json(CASM::jsonParser &json) const {
  json.put_array();
  for(int i = 0; i < size(); i++) {
    json.push_back((*this)[i]);
  }
  return json;
}

void EnergySet::from_json(const CASM::jsonParser &json) {
  name = "-";
  E_vec_ready = false;
  hull_found = 0;
  Nstruct_set = 0;

  clear();
  for(int i = 0; i < json.size(); i++) {
    add(json[i]);
  }
  set_Nstruct_on();
  set_E_vec();

}

CASM::jsonParser &to_json(const EnergySet &nrgset, CASM::jsonParser &json) {
  return nrgset.to_json(json);
}

void from_json(EnergySet &nrgset, const CASM::jsonParser &json) {
  nrgset.from_json(json);
}


#endif // EnergySet_CC
