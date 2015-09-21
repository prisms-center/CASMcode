/*
 *  Correlation.hh
 */

#ifndef Correlation_HH
#define Correlation_HH

#include <string>
#include "jsonParser.hh"
#include "BP_Vec.hh"
#include "casm/external/Eigen/Dense"

class ECISet;
class EnergySet;

class Correlation : public BP::BP_Vec< BP::BP_Vec< double> > {

  /// detected input file format: "text" or "json"
  std::string m_format;

public:

  //Correlation() {}

  // Construct from 'corr' file
  Correlation(std::string corr_in_filename);

  std::string format() const;

  // reduce this to only including the subset of clusters indicated by their indices in 'index_list'
  void cluster_subset(BP::BP_Vec<int> &index_list);

  // write a 'corr' file
  void write(const std::string &filename, std::string format) const;

  // write covariance
  void write_covar(const ECISet &eci_in, const EnergySet &nrg_set, const std::string &fname, std::string format) const;

  // return a vector of 'a' values used for calculating LOOCV score
  BP::BP_Vec<double> cv_a(const Eigen::MatrixXd &C, const EnergySet &nrg_set) const;
};

CASM::jsonParser &to_json(const Correlation &corr, CASM::jsonParser &json);

void from_json(Correlation &corr, const CASM::jsonParser &json);

#endif // Correlation_HH
