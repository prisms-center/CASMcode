/*
 *  GeneticAlgorithm.hh
 */

#ifndef GeneticAlgorithm_HH
#define GeneticAlgorithm_HH

#include <iostream>
#include "BP_GVec.hh"
#include "BP_RVG.hh"

class ECISet;
class Correlation;
class EnergySet;
class MTRand;

class GeneticAlgorithm {
private:
  int Nmin;
  int Nmax;
  int Nchildren;
  //int Ngenerations;
  int Nmutations;
  BP::BP_GVec< ECISet> population;

public:

  GeneticAlgorithm();

  GeneticAlgorithm(BP::BP_Vec< ECISet> &wannabeparents, int min, int max, int children, int mutations);

  ECISet &operator[](unsigned long int i);
  unsigned long int size();



  ECISet mate(const ECISet &Mom, const ECISet &Dad, MTRand &mtrand);

  void set_population(BP::BP_Vec< ECISet> &wannabeparents);

  void set(int min, int max, int children, int mutations);

  void set_Nmin(int i);
  void set_Nmax(int i);
  void set_Nchildren(int i);
  void set_Nmutations(int i);

  BP::BP_Vec<ECISet> get_population() const;

  int get_Nmin() const;
  int get_Nmax() const;
  int get_Nchildren() const;
  int get_Nmutations() const;

  int get_best_index() const;


  void run(MTRand &mtrand, Correlation &corr, EnergySet &nrg_set, bool print_steps);

  void run_dir(MTRand &mtrand, Correlation &corr, EnergySet &nrg_set, bool print_steps);

  void run_dfs(long int Nstop, MTRand &mtrand, Correlation &corr, EnergySet &nrg_set, bool print_steps);

private:
  void run_ga(MTRand &mtrand, Correlation &corr, EnergySet &nrg_set, bool print_steps, int min_option, long int Nstop);
};


#endif // GeneticAlgorithm_HH
