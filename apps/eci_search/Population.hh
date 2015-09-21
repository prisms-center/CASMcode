/*
 *  Population.hh
 */

#ifndef Population_HH
#define Population_HH

#include <string>
#include <unistd.h>
#include "BP_GVec.hh"
#include "ECISet.hh"
#include "Functions.hh"
#include "Minimize.hh"


// This class contains a population of ECISets
//   It provides several methods to optimize that population:
//
//   direct: direct (steepest-descent) minimization
//   dfs: depth first search minimization
//   ga: genetic algorithm
//
class Population {
private:

  EnergySet nrg;
  ECISet eci_base;
  Correlation corr;

  MTRand mtrand;

  // the population of ECISets
  //   these should be added via the 'add' methods to ensure that ECISet.set_data is called
  BP::BP_Vec<ECISet> population;

  /// detected input bit_string file format: "text" or "json"
  std::string m_format;

public:

  Population(int Nmin,
             int Nmax,
             const std::string &nrgset_filename,
             const std::string &eciset_filename,
             const std::string &corr_filename):
    nrg(nrgset_filename), eci_base(eciset_filename, Nmin, Nmax), corr(corr_filename) {
    if(TEST)
      mtrand.seed(1);
    m_format = nrg.format();
  }

  std::string format() const;

  // Use the eci_base for the population
  void populate_with_base();

  // Fill in the population with randomly generated ECISets, with Nmin ECI
  void populate(int Npop);

  // Fill in the population with ECISets listed in a file
  void populate(const std::string &bit_string_filename);

  int size() const {
    return population.size();
  };

  // Add a random ECISet (with Nmin eci) to the population
  void add();

  // Add an ECISet with ECISetState to the population
  void add(const ECISetState &state);

  // Add an ECISet to the population with given bit_string
  void add(const std::string &bit_string);

  // Minimize each ECISet in the population using a direct (steepest-descent) minimization
  void direct(int pthreads = -1, int mthreads = -1);

  // Minimize each ECISet in the population using a depth-first search
  void dfs(int Nstop, int pthreads = -1, int mthreads = -1);

  // Genetic algorithm to optimize the population
  void ga(int Ngen, int Nmut, int pthreads = -1, int mode = 0, int Nstop = 1);

  // Genetic algorithm, with direct minimization of each child
  void ga_dir(int Ngen, int Nmut, int pthreads = -1);

  // Genetic algorithm, with dfs minimization of each child
  void ga_dfs(int Nstop, int Ngen, int Nmut, int pthreads = -1);

  // Calculate cv for all ECISet in population
  void calc(int pthreads = -1);

  // Calculate cv for all ECISets with N eci, based on eci_base
  void calc_all(int N, int pthreads);

  // Calculate & write hull, energy.clex, etc. files for population[index]
  void calc_details(int index, std::string out_format, double hulltol = 1.0e-14);

  // Print statistics on the value of the ECI across the population
  void calc_ecistats(int pthreads = -1);

  // Return a string containing the status of ECISets in population
  std::string population_status() const;

  // Write population bitstrings to file
  void write(std::string filename, std::string format) const;

private:

  // Add an ECISet to the population
  void add(const ECISet &eci);

  // Return a string containing the status of ECISets in population during a minimization
  std::string minimization_status(const BP::BP_Vec<Minimize> &minimization, BP::BP_Vec<bool> &completed) const;

  // Return a string containing the status of the gene pool
  std::string gene_pool_status(const BP::BP_Vec<ECISetState> &gene_pool) const;

};

#endif // Population_HH
