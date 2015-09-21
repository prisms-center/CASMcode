/*
 *  Correlation.cc
 */

#ifndef Correlation_CC
#define Correlation_CC

#include <iostream>
#include "BP_Parse.hh"
#include "Correlation.hh"
#include "ECISet.hh"
#include "EnergySet.hh"
#include "Functions.hh"

// Construct from 'corr' file
Correlation::Correlation(std::string corr_in_filename) {

  m_format = get_format_from_ext(corr_in_filename);

  if(m_format == "text") {
    BP::BP_Parse file(corr_in_filename);

    int Nclust;
    int Ncon;
    std::string s1;
    BP::BP_Vec<double> list;

    s1 = file.getline();
    Nclust = BP::next_int(s1);

    s1 = file.getline();
    Ncon = BP::next_int(s1);

    //std::cout << "Reading " << corr_in_filename << std::endl;
    //std::cout << "  Nclust: " << Nclust << " Ncon: " << Ncon << std::endl;

    //"clusters" label
    s1 = file.getline();

    // data
    do {
      list = file.getline_double();
      if(list.size() != 0)
        add(list);
    }
    while(file.eof() == false);

    //std::cout << "val.size: " << val.size() << std::endl;
    if(size() != Ncon) {
      std::cout << "Error reading '" << corr_in_filename << "': stated #configurations == " << Ncon << ", but found #configurations == " << size() << std::endl;
      exit(1);
    }

    for(int i = 0; i < size(); i++) {
      //std::cout << "val[i].size(): " << val[i].size() << std::endl;
      if((*this)[i].size() != Nclust) {
        std::cout << "Error: the eci.in file stated #clusters == " << Nclust << ", but reading '" << corr_in_filename << "' found #clusters == " << (*this)[i].size() << " for configuration " << i << std::endl;
        exit(1);
      }
    }
  }
  else if(m_format == "json") {
    CASM::jsonParser json(corr_in_filename);
    from_json(*this, json);
  }
  else {
    std::cout << "Unexpected format option for Correlation constructor" << std::endl;
    std::cout << "  Expected 'text' or 'json', but received: " << m_format << std::endl;
    exit(1);
  }


}

std::string Correlation::format() const {
  return m_format;
}

// reduce this to only including the subset of clusters indicated by their indices in 'index_list'
void Correlation::cluster_subset(BP::BP_Vec<int> &index_list) {
  unsigned long int i, j;

  BP::BP_Vec< BP::BP_Vec< double> > tmp = (*this);

  for(i = 0; i < (*this).size(); i++)
    (*this)[i].clear();

  for(i = 0; i < tmp.size(); i++) {
    for(j = 0; j < index_list.size(); j++) {
      (*this)[i].add(tmp[i][ index_list[j]]);
    }
  }

}

// write a 'corr' file
void Correlation::write(const std::string &filename, std::string format) const {

  if(format == "default")
    format = m_format;

  if(format == "text") {
    unsigned long int i, j;
    BP::BP_Write file(rm_json_ext(filename));
    file.newfile();

    file << (*this)[0].size() << " # number of clusters" << std::endl;
    file << (*this).size() << " # number of configurations" << std::endl;
    file << "clusters" << std::endl;
    for(i = 0; i < size(); i++) {
      for(j = 0; j < (*this)[i].size(); j++) {
        file << "   " << (*this)[i][j] ;
      }
      file << "\n";
    }
  }
  else if(format == "json") {
    CASM::jsonParser json;
    to_json(*this, json).write(json_ext(filename));
  }
  else {
    std::cout << "Unexpected format option for Correlation constructor" << std::endl;
    std::cout << "  Expected 'text' or 'json', but received: " << format << std::endl;
    exit(1);
  }
}

void Correlation::write_covar(const ECISet &eci_in, const EnergySet &nrg_set, const std::string &filename, std::string format) const {
  //std::cout << "begin fit()" << std::endl;

  int Nstruct = nrg_set.get_Nstruct_on();
  int Nclust = eci_in.get_Nclust_on();

  //std::cout << "fit Nstruct: " << Nstruct << "  Nclust: " << Nclust << std::endl;
  Eigen::MatrixXd corr_matrix(Nstruct, Nclust);

  //std::cout << " set correlation matrix" << std::endl;
  eci_in.set_correlation_matrix(corr_matrix, (*this), nrg_set);

  //std::cout << " init svd" << std::endl;
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(corr_matrix, Eigen::ComputeThinU | Eigen::ComputeThinV);


  if(format == "default")
    format = m_format;

  if(format == "text") {

    BP::BP_Write file(rm_json_ext(filename));
    file.newfile();
    file << "Singular values of correlation matrix used in fit:" << std::endl;
    for(int i = 0; i < Nclust; i++) {
      file << setiosflags(ios::right | ios::showpoint | ios::fixed) << setw(12) << setprecision(7) << svd.singularValues()[i] << ' ';
    }
    file << ' ' << std::endl;
    file << "Right singular vectors of correlation matrix used in fit:" << std::endl;
    for(int i = 0; i < Nclust; i++) {
      for(int j = 0; j < Nclust; j++) {
        file << setiosflags(ios::right | ios::showpoint | ios::fixed) << setw(12) << setprecision(7) << svd.matrixV()(i, j) << ' ';
      }
      file << ' ' << std::endl;
    }
  }
  else if(format == "json") {
    CASM::jsonParser json;
    json.put_obj();
    json["singular_values"].put_array(Nclust);
    for(int i = 0; i < Nclust; i++) {
      json["singular_values"].push_back(svd.singularValues()[i]);
    }
    json["right_singular_vectors"].put_array(Nclust);
    for(int i = 0; i < Nclust; i++) {
      json["right_singular_vectors"].push_back(CASM::jsonParser::null());
      json["right_singular_vectors"][i].put_array(Nclust);
      for(int j = 0; j < Nclust; j++) {
        json["right_singular_vectors"][i][j] = svd.matrixV()(i, j);
      }
    }
    json.write(json_ext(filename));

  }
  else {
    std::cout << "Unexpected format option for Correlation constructor" << std::endl;
    std::cout << "  Expected 'text' or 'json', but received: " << format << std::endl;
    exit(1);
  }

}

CASM::jsonParser &to_json(const Correlation &corr, CASM::jsonParser &json) {
  return to_json((BP::BP_Vec< BP::BP_Vec< double> > &) corr, json);
}

void from_json(Correlation &corr, const CASM::jsonParser &json) {
  from_json((BP::BP_Vec< BP::BP_Vec< double> > &) corr, json);
}

#endif // Correlation_CC
