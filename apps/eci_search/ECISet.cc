/*
 *  ECISet.CC
 */

#ifndef ECISet_CC
#define ECISet_CC

#include "ECISet.hh"
#include "Functions.hh"

ECI::ECI(const CASM::jsonParser &json) {
  from_json(*this, json);
}

CASM::jsonParser &to_json(const ECI &eci, CASM::jsonParser &json) {
  json.put_obj();
  json["label"] = eci.label;
  json["weight"] = eci.weight;
  json["mult"] = eci.mult;
  json["size"] = eci.size;  // 'branch'
  if(eci.orbit != -1)
    json["orbit"] = eci.orbit;
  if(eci.bfunc != -1)
    json["bfunc"] = eci.bfunc;
  json["length"] = eci.length;
  if(eci.weight == 1) {
    json["value"] = eci.value;
    json["value/mult"] = eci.value / eci.mult;
  }
  else {
    json["value"] = 0.0;
    json["value/mult"] = 0.0;
  }
  json["hierarchy"] = eci.hierarchy;
  return json;
}

void from_json(ECI &eci, const CASM::jsonParser &json) {
  json["label"].get(eci.label);
  if(json["weight"].is_string()) {
    if(json["weight"].get<std::string>() == "FixOn") {
      eci.weight = 1;
    }
    else if(json["weight"].get<std::string>() == "FixOff") {
      eci.weight = 0;
    }
  }
  else {
    json["weight"].get(eci.weight);
  }
  json["mult"].get(eci.mult);
  json["size"].get(eci.size);  // 'branch'
  json.get_else(eci.orbit, "orbit", -1);
  json.get_else(eci.bfunc, "bfunc", -1);
  json["length"].get(eci.length);
  json["value"].get(eci.value);
  json["hierarchy"].get(eci.hierarchy);
}



ECISet::ECISet(): cv(UK), rms(UK), Nmin(0), Nmax(-1) {}

ECISet::ECISet(const std::string &eci_in_filename) {

  m_format = get_format_from_ext(eci_in_filename);

  if(m_format == "text") {
    read(eci_in_filename);
  }
  else if(m_format == "json") {
    CASM::jsonParser json(eci_in_filename);
    from_json(json);
  }
  else {
    std::cout << "Unexpected format option for ECISet constructor" << std::endl;
    std::cout << "  Expected 'text' or 'json', but received: " << m_format << std::endl;
    exit(1);
  }
  set_Nmin(0);
  set_Nmax(size());
}

ECISet::ECISet(const std::string &eci_in_filename, int _Nmin, int _Nmax) {

  m_format = get_format_from_ext(eci_in_filename);

  if(m_format == "text") {
    read(eci_in_filename);
  }
  else if(m_format == "json") {
    CASM::jsonParser json(eci_in_filename);
    from_json(json);
  }
  else {
    std::cout << "Unexpected format option for ECISet constructor" << std::endl;
    std::cout << "  Expected 'text' or 'json', but received: " << m_format << std::endl;
    exit(1);
  }

  set_Nmin(_Nmin);
  set_Nmax(_Nmax);
}

std::string ECISet::format() const {
  return m_format;
}

void ECISet::read(const std::string &eci_in_filename) {
  clear();
  fix_pos_list.clear();
  fix_val_list.clear();
  fix_index.clear();

  Nclust_set = 0;
  cv = UK;
  rms = UK;

  BP::BP_Parse file(eci_in_filename);
  ECI eci;
  // header line
  file.getline();

  BP::BP_Vec<std::string> entry;
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

  while(fix_index.size() < size())
    fix_index.add(-1);

  get_Nclust_on();

  //std::cout << "eci_set.size(): " << eci_list.size() << std::endl;
  //for( int i=0; i<eci_list.size(); i++)
  //	std::cout << eci_list[i] << std::endl;
}

void ECISet::set_Nmin(int _Nmin) {
  Nmin = _Nmin;
}

int ECISet::get_Nmin() const {
  return Nmin;
}

void ECISet::set_Nmax(int _Nmax) {
  Nmax = _Nmax;
}

int ECISet::get_Nmax() const {
  return Nmax;
}

bool ECISet::toggle_allowed(int i) const {
  if((*this)[i].weight == 0) {
    if(fix_index[i] == 0) {
      return false;
    }

    if((get_Nclust_on() + 1 >= Nmin) && (get_Nclust_on() + 1 <= Nmax))
      return true;
    else
      return false;

  }
  else {
    if(fix_index[i] == 1) {
      return false;
    }

    if((get_Nclust_on() - 1 >= Nmin) && (get_Nclust_on() - 1 <= Nmax))
      return true;
    else
      return false;
  }
}

int ECISet::get_Nclust_on() const {
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
}

void ECISet::set_values(const Eigen::VectorXd &ECI) {
  int in_i = 0;
  for(int i = 0; i < size(); i++)
    if((*this)[i].weight != 0) {
      (*this)[i].value = ECI(in_i);
      in_i++;
    }
}

void ECISet::set_values_and_weights(const Eigen::MatrixXd &ECI) {
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
}

void ECISet::clear_weights() {
  for(int i = 0; i < size(); i++)
    (*this)[i].weight = 0;
  Nclust_on = 0;
  Nclust_set = 1;
}

void ECISet::set_clust_off(int i) {
  if((*this)[i].weight != 0) {
    Nclust_on--;
    (*this)[i].weight = 0;
  }
}

void ECISet::set_clust_on(int i) {
  if((*this)[i].weight == 0) {
    Nclust_on++;
    (*this)[i].weight = 1;
  }
}

void ECISet::toggle_clust(int i) {
  if((*this)[i].weight == 0) {
    Nclust_on++;
    (*this)[i].weight = 1;
  }
  else {
    Nclust_on--;
    (*this)[i].weight = 0;
  }
}

ECISet &ECISet::operator=(const BP::BP_Comb &comb) {
  if(comb.size() == size()) {
    for(int i = 0; i < comb.size(); i++)
      (*this)[i].weight = comb[i];
  }
  Nclust_set = 0;

  return (*this);
}

ECISetState ECISet::get_state() const {
  return ECISetState(get_bit_string(), get_Nclust_on(), get_cv(), get_rms());
};

void ECISet::set_state(const ECISetState &state) {
  set_bit_string(state.bit_string);
  cv = state.cv;
  rms = state.rms;
};


std::string ECISet::get_bit_string() const {
  std::stringstream ss;
  for(int i = 0; i < size(); i++)
    ss << (*this)[i].weight;
  return ss.str();
}

void ECISet::set_bit_string(const std::string &s) {
  if(s.size() != size()) {
    std::cout << "Error in ECISet::set_bit_string(). string length mismatch" << std::endl;
    std::cout << "  Attempting to set length " << s.size() << " bit_string: " << s << std::endl;
    std::cout << "  But this ECISet is size: " << size() << std::endl;
    exit(1);
  }

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

unsigned long int ECISet::size() const {
  return BP::BP_Vec<ECI>::size();
}

int ECISet::get_weight(unsigned long int i) const {
  return (*this)[i].weight;
}

int ECISet::get_size(unsigned long int i) const {
  return (*this)[i].size;
}

double ECISet::get_value(unsigned long int i) const {
  return (*this)[i].value;
}

double ECISet::get_cv() const {
  return cv;
}

double ECISet::get_rms() const {
  return rms;
}

int ECISet::get_Nstruct() const {
  return Nstruct;
}

bool ECISet::get_singular() const {
  return singular;
}

void ECISet::set_cv(double i1) {
  cv = i1;
}

void ECISet::set_rms(double i1) {
  rms = i1;
}

void ECISet::set_Nstruct(int _Nstruct) {
  Nstruct = _Nstruct;
}

void ECISet::set_data(const EnergySet &_nrg, const Correlation &_corr) {
  nrg = &_nrg;
  corr = &_corr;
};

bool ECISet::fix_ok() {
  for(int i = 0; i < fix_pos_list.size(); i++) {
    if((*this)[ fix_pos_list[i]].weight != fix_val_list[i])
      return false;
  }
  return true;
}

void ECISet::fix(int i1, bool i2) {
  fix_pos_list.add(i1);
  fix_val_list.add(i2);
  while(fix_index.size() <= i1)
    fix_index.add(-1);
  fix_index[i1] = (int) i2;
}

void ECISet::unfix(int i1) {
  while(fix_index.size() <= i1)
    fix_index.add(-1);
  fix_index[i1] = -1;
  for(int i = 0; i < fix_pos_list.size(); i++) {
    if(fix_pos_list[i] == i1) {
      fix_pos_list.remove(i);
      fix_val_list.remove(i);
    }
  }
}

void ECISet::unfix_all() {
  fix_pos_list.erase();
  fix_val_list.erase();
  fix_index = BP::BP_Vec<int>(size(), -1);
}

void ECISet::set_fix() {
  for(int i = 0; i < fix_pos_list.size(); i++) {
    (*this)[ fix_pos_list[i]].weight = fix_val_list[i];
  }
  Nclust_set = 0;
}

void ECISet::set_Rmax_fix(double Rmax) {
  for(int i = 0; i < size(); i++) {
    if((*this)[ i].length > Rmax)
      fix(i, 0);
  }
}

void ECISet::write_ECI(ostream &sout) {
  for(int i = 0; i < size(); i++) {
    if(get_weight(i) != 0) {
      sout << "i: " << i << " value: " << get_value(i) << std::endl;

    }

  }

}

void ECISet::write_ECIin(std::string filename, std::string format) {

  if(format == "default")
    format = m_format;

  if(format == "text") {
    BP::BP_Write file(rm_json_ext(filename));
    file.newfile();

    file << "label" << "     " << "weight" << "     " << "mult" << "     " << "size" << "     " << "length" << "     " << "heirarchy" << std::endl;
    for(int i = 0; i < size(); i++)
      file << (*this)[i] << std::endl;
  }
  else if(format == "json") {
    CASM::jsonParser json;
    to_json(json).write(json_ext(filename));
  }
  else {
    std::cout << "Unexpected format option for ECISET::write_ECIin" << std::endl;
    std::cout << "  Expected 'text' or 'json', but received: " << format << std::endl;
    exit(1);
  }
}

void ECISet::write_ECIin(BP::BP_Vec<int> &index_list, std::string filename, std::string format) {

  if(format == "default")
    format = m_format;

  if(format == "text") {
    BP::BP_Write file(rm_json_ext(filename));
    file.newfile();

    file << "label" << "     " << "weight" << "     " << "mult" << "     " << "size" << "     " << "length" << "     " << "heirarchy" << std::endl;
    for(int i = 0; i < index_list.size(); i++)
      file << (*this)[index_list[i] ] << std::endl;
  }
  else if(format == "json") {
    CASM::jsonParser json;
    to_json(json).write(json_ext(filename));
  }
  else {
    std::cout << "Unexpected format option for ECISET::write_ECIin" << std::endl;
    std::cout << "  Expected 'text' or 'json', but received: " << format << std::endl;
    exit(1);
  }
}

void ECISet::write_ECIout(std::string filename, EnergySet &nrg_set, std::string format) {

  if(format == "default")
    format = m_format;

  if(format == "text") {
    BP::BP_Write file(rm_json_ext(filename));
    file.newfile();

    file << "The number of clusters in fit is " << get_Nclust_on() << std::endl;
    file << "The number of structure in fit is " << nrg_set.get_Nstruct_on() << std::endl;
    file << "The rms is " << setprecision(9) << get_rms() << std::endl;
    file << "Number of structures in WCV calculation is Stru_wcv=" << nrg_set.get_Nstruct_on() << std::endl;
    file << "The WCV score is " << get_cv() << std::endl;
    file << "The rms based on the " << "N/A" << " structures having vasp E but not in the fit is: " << "N/A" << std::endl;
    file << setiosflags(ios::right) << setw(15) << "ECIs" << setiosflags(ios::right) << setw(15) << "ECI/mult" << setiosflags(ios::right) << setw(15) << "Cluster#" << std::endl;
    for(int i = 0; i < size(); i++) {
      if(get_weight(i) != 0) {
        file << setiosflags(ios::right) << setiosflags(ios::fixed) << setw(15) << setprecision(9) << (*this)[i].value ;
        file << setiosflags(ios::right) << setiosflags(ios::fixed) << setw(15) << setprecision(9) << (*this)[i].value / (*this)[i].mult ;
        file << setiosflags(ios::right) << setiosflags(ios::fixed) << setw(15) << setprecision(9) << (*this)[i].label << std::endl;
      }
    }
  }
  else if(format == "json") {
    CASM::jsonParser json;
    to_json(json).write(json_ext(filename));
  }
  else {
    std::cout << "Unexpected format option for ECISET::write_ECIin" << std::endl;
    std::cout << "  Expected 'text' or 'json', but received: " << format << std::endl;
    exit(1);
  }
}


void ECISet::write_ECIout(BP::BP_Vec<int> &index_list, std::string filename, EnergySet &nrg_set, std::string format) {

  if(format == "default")
    format = m_format;

  if(format == "text") {
    BP::BP_Write file(filename);
    file.newfile();

    file << "The number of clusters in fit is " << get_Nclust_on() << std::endl;
    file << "The number of structure in fit is " << nrg_set.get_Nstruct_on() << std::endl;
    file << "The rms is " << setprecision(9) << get_rms() << std::endl;
    file << "Number of structures in WCV calculation is Stru_wcv=" << nrg_set.get_Nstruct_on() << std::endl;
    file << "The WCV score is " << get_cv() << std::endl;
    file << "The rms based on the " << "N/A" << " structures having vasp E but not in the fit is: " << "N/A" << std::endl;
    file << setiosflags(ios::right) << setw(15) << "ECIs" << setiosflags(ios::right) << setw(15) << "ECI/mult" << setiosflags(ios::right) << setw(15) << "Cluster#" << std::endl;
    for(int i = 0; i < index_list.size(); i++) {
      if(get_weight(i) != 0) {
        file << setiosflags(ios::right) << setiosflags(ios::fixed) << setw(15) << setprecision(9) << (*this)[index_list[i]].value ;
        file << setiosflags(ios::right) << setiosflags(ios::fixed) << setw(15) << setprecision(9) << (*this)[index_list[i]].value / (*this)[index_list[i]].mult ;
        file << setiosflags(ios::right) << setiosflags(ios::fixed) << setw(15) << setprecision(9) << index_list[i] << std::endl;
      }
    }
  }
  else if(format == "json") {
    CASM::jsonParser json;
    json.put_obj();

    json["eci"].put_array();
    for(int i = 0; i < index_list.size(); i++) {
      json["eci"].push_back((*this)[index_list[i]]);
      if(fix_index[index_list[i]] == 0) {
        json["eci"][i]["weight"] = "FixOff";
      }
      else if(fix_index[index_list[i]] == 1) {
        json["eci"][i]["weight"] = "FixOn";
      }
    }
    json.write(json_ext(filename));
  }
  else {
    std::cout << "Unexpected format option for ECISET::write_ECIin" << std::endl;
    std::cout << "  Expected 'text' or 'json', but received: " << format << std::endl;
    exit(1);
  }

}


bool ECISet::operator==(const ECISet &E) const {
  if(size() != E.size())
    return false;

  for(int i = 0; i < size(); i++)
    if((*this)[i].weight != E.get_weight(i))
      return false;

  return true;
}

void ECISet::randomize(int Nclust, MTRand &mtrand) {
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

}

void ECISet::fit() {
  fit(*corr, *nrg, singular);
}

void ECISet::fit(const Correlation &_corr, const EnergySet &nrg_set, bool &_singular) {
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
  bool check_singular = 1;

  //std::cout << "begin fit()" << std::endl;
  double inf = 1.0e20;
  Nstruct = nrg_set.get_Nstruct_on();
  int Nclust = get_Nclust_on();

  //std::cout << "fit Nstruct: " << Nstruct << "  Nclust: " << Nclust << std::endl;
  Eigen::MatrixXd corr_matrix(Nstruct, Nclust);

  //std::cout << " set correlation matrix" << std::endl;
  set_correlation_matrix(corr_matrix, _corr, nrg_set);

  //std::cout << " init svd" << std::endl;
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(corr_matrix, Eigen::ComputeThinU | Eigen::ComputeThinV);

  if(check_singular) {
    Eigen::VectorXd S;

    // do svd
    //std::cout << " do svd" << std::endl;
    S = svd.singularValues();
    //std::cout << "singular values: " << S << "\n";

    //std::cout << " check if singular" << std::endl;
    _singular = check_if_singular(S);

    if(_singular) {
      set_cv(inf);
      set_rms(inf);
      return;
    }

  }
  else {
    _singular = false;
  }

  Eigen::VectorXd ECI(Nclust);
  Eigen::VectorXd Err(Nstruct);
  unsigned long int i, ii;

  //std::cout << " set energy vector" << std::endl;
  if(!nrg_set.E_vec_is_ready()) {
    std::cout << "Error in ECISet::fit.  nrg_set.E_vec is not ready." << std::endl;
    exit(1);
    //nrg_set.set_E_vec();
  }

  // Calculate rms fit
  //std::cout << " rms fit" << std::endl;
  ECI = svd.solve(nrg_set.get_E_vec());

  //std::cout << "ECI: " << ECI << std::endl;

  // Calculate residuals
  //std::cout << " calc residuals" << std::endl;
  Err = corr_matrix * ECI - nrg_set.get_E_vec();

  //std::cout << "n nrg Fit: " << nrg_set.get_E_vec().size() << std::endl;

  // Calculate Leave One Out cv score
  // LOOCV = (1.0/Nnrg)*sum_i{ (e_i / 1 - X_i*((X^T*X)^-1)*X_i^T)^2 }
  // e_i is residual
  // X is corr_matrix
  // X_i is row i of X

  // -- compute cv_a_i = X_i*((X^T*X)^-1)*X_i^T ahead of time
  BP::BP_Vec<double> cv_a = set_cv_a(corr_matrix, nrg_set);

  //std::cout << " calc rms and cv" << std::endl;
  rms = 0.0;
  cv = 0.0;
  ii = 0;
  for(i = 0; i < nrg_set.size(); i++)
    if(nrg_set.get_weight(i) != 0) {
      rms += Err(ii) * Err(ii);
      cv += BP::sqr(Err(ii) / (1.0 - cv_a[ii]));
      //std::cout << "cv_a " << ii << " :: " << cv_a[ii] << std::endl;
      ii++;
    }

  rms = sqrt(rms / Nstruct);
  cv = sqrt(cv / Nstruct);

  //std::cout << " set eci values" << std::endl;
  set_values(ECI);

  //std::cout << "finish fit_Eigen()" << std::endl;
  return;
}

void ECISet::check_cv(const Correlation &_corr, const EnergySet &_nrg_set, bool &_singular) {
  fit(_corr, _nrg_set, _singular);
  double fitted_cv = cv;
  std::cout << "fitted_cv: " << fitted_cv << "  fitted_rms: " << rms << std::endl;

  // for each data in nrg_set, leave it out, fit, calculate predicted value
  BP::BP_Vec<double> nrg_data;
  BP::BP_Vec<double> nrg_fit;
  double orig_weight;

  EnergySet nrg_set = _nrg_set;
  double sqr_sum = 0.0;
  double wsqr_sum = 0.0;

  for(int i = 0; i < nrg_set.size(); i++) {

    // remove
    orig_weight = nrg_set.get_weight(i);
    if(orig_weight != 0.0) {

      nrg_data.add(nrg_set.get_Ef(i));
      nrg_set.set_weight(i, 0.0);
      nrg_set.set_E_vec();

      // fit
      fit(_corr, nrg_set, _singular);

      // calc clex nrg
      nrg_fit.add(EnergySet::calc_clex(_corr, *this, i));

      // reset
      nrg_set.set_weight(i, orig_weight);

      std::cout << "i: " << i << " data: " << nrg_data.last() << "  fit: " << nrg_fit.last() << std::endl;
      sqr_sum += BP::sqr(nrg_data.last() - nrg_fit.last());
      wsqr_sum += BP::sqr(orig_weight * (nrg_data.last() - nrg_fit.last()));
    }
  }
  std::cout << "nFit: " << nrg_data.size() << std::endl;
  std::cout << "LOOCV: " << sqrt(sqr_sum / nrg_data.size()) << "  fitted_cv: " << fitted_cv << std::endl;
  std::cout << "wLOOCV: " << sqrt(wsqr_sum / nrg_data.size()) << "  fitted_cv: " << fitted_cv << std::endl;
}

void *ECISet::fit_threaded(void *arg) {
  ((ECISet *) arg)->fit();
  return NULL;
}

void ECISet::set_correlation_matrix(Eigen::MatrixXd &A, const Correlation &_corr, const EnergySet &nrg_set) const {
  // set A to be the correlation matrix, including weights, and only the rows and columns being fit
  // assumes A is already the right size

  int i, j, in_i, in_j;

  in_i = 0;
  for(i = 0; i < nrg_set.size(); i++) {
    if(nrg_set.get_weight(i) != 0) {
      in_j = 0;
      for(j = 0; j < _corr[i].size(); j++) {
        if(this->get_weight(j) != 0) {
          A(in_i, in_j) = nrg_set.get_weight(i) * _corr[i][j];
          in_j++;
        }
      }
      in_i++;
    }
  }

  //std::cout << " in_i: " << in_i << " in_j: " << in_j << std::endl;
}

// private:

bool ECISet::check_if_singular(const Eigen::VectorXd &S) const {
  for(unsigned long int i = 0; i < S.size(); i++)
    if(fabs(S(i)) < 1.0e-4)	// CONSTANT
      return true;
  return false;
}

BP::BP_Vec<double> ECISet::set_cv_a(const Eigen::MatrixXd &C, const EnergySet &nrg_set) const {
  //cout << "begin set_cv_a()" << endl;
  // -- compute a_i = X_i*((X^T*X)^-1)*X_i^T ahead of time, once per corr

  Eigen::MatrixXd XTXinv;
  unsigned long int i, ii;
  BP::BP_Vec<double> cv_a(0, nrg_set.get_Nstruct_on());

  XTXinv = (C.transpose() * C).inverse();

  ii = 0;
  for(i = 0; i < nrg_set.size(); i++)
    if(nrg_set.get_weight(i) != 0) {
      cv_a.add(C.row(ii) * XTXinv * C.row(ii).transpose());
      ii++;
    }
  return cv_a;
}


CASM::jsonParser &ECISet::to_json(CASM::jsonParser &json) const {

  json.put_obj();

  if(cv != UK) {
    json["Nstruct"] = Nstruct;
    json["Nclust"] = Nclust_on;
    json["wcv"] = cv;
    json["rms"] = rms;
  }


  json["eci"].put_array();
  for(int i = 0; i < size(); i++) {
    json["eci"].push_back((*this)[i]);
    if(fix_index[i] == 0) {
      json["eci"][i]["weight"] = "FixOff";
    }
    else if(fix_index[i] == 1) {
      json["eci"][i]["weight"] = "FixOn";
    }
  }
  return json;
}

void ECISet::from_json(const CASM::jsonParser &json) {
  clear();
  fix_pos_list.clear();
  fix_val_list.clear();
  fix_index.clear();

  Nstruct = 0;
  Nclust_set = 0;
  cv = UK;
  rms = UK;

  json.get_if(Nstruct, "Nstruct");
  json.get_if(cv, "wcv");
  json.get_if(rms, "rms");


  for(int i = 0; i < json["eci"].size(); i++) {
    if(json["eci"][i]["weight"].is_string()) {
      if(json["eci"][i]["weight"].get<std::string>() == "FixOn") {
        fix((*this).size(), 1);
      }
      else if(json["eci"][i]["weight"].get<std::string>() == "FixOff") {
        fix((*this).size(), 0);
      }
    }
    add(json["eci"][i]);

  }

  while(fix_index.size() < size())
    fix_index.add(-1);

  get_Nclust_on();

}

CASM::jsonParser &to_json(const ECISet &eciset, CASM::jsonParser &json) {
  return eciset.to_json(json);
}

void from_json(ECISet &eciset, const CASM::jsonParser &json) {
  eciset.from_json(json);
}

#endif // ECISet_CC
