/*
 *  GeneticAlgorithm.cc
 */

#ifndef GeneticAlgorithm_CC
#define GeneticAlgorithm_CC


#include "ECISet.hh"
#include "GeneticAlgorithm.hh"
#include "casm/external/MersenneTwister/MersenneTwister.h"


GeneticAlgorithm::GeneticAlgorithm() {}

GeneticAlgorithm::GeneticAlgorithm(BP::BP_Vec< ECISet> &wannabeparents, int min, int max, int children, int mutations) {
  set_population(wannabeparents);
  Nmin = min;
  Nmax = max;
  Nchildren = children;
  Nmutations = mutations;
}

ECISet &GeneticAlgorithm::operator[](unsigned long int i) {
  return population[i];
}

unsigned long int GeneticAlgorithm::size() {
  return population.size();
}

ECISet GeneticAlgorithm::mate(const ECISet &Mom, const ECISet &Dad, MTRand &mtrand) {
  //std::cout << "begin GeneticAlgorithm:mate()" << std::endl;
  unsigned long int j;
  // mate
  //std::cout << "    mate" << std::endl;
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
  //std::cout << "    mutate" << std::endl;
  for(int i = 0; i < child.size(); i++) {
    if(mtrand.randExc() < (1.0 * Nmutations) / child.size()) {
      child.toggle_clust(i);
      if(!child.fix_ok())
        child.toggle_clust(i);
    }
  }

  // enforce Nmin limit
  //std::cout << "    enforce Nmin" << std::endl;
  while(child.get_Nclust_on() < Nmin) {
    j = mtrand.randInt(child.size() - 1);
    child.set_clust_on(j);
    if(!child.fix_ok())
      child.set_clust_off(j);
  }

  // enforce Nmax limit
  //std::cout << "    enforce Nmax" << std::endl;
  while(child.get_Nclust_on() > Nmax) {
    j = mtrand.randInt(child.size() - 1);
    child.set_clust_off(j);
    if(!child.fix_ok())
      child.set_clust_on(j);
  }
  //std::cout << "finish GeneticAlgorithm:mate()" << std::endl;
  return child;
}

void GeneticAlgorithm::set_population(BP::BP_Vec< ECISet> &wannabeparents) {
  for(int i = 0; i < wannabeparents.size(); i++) {
    population.add(wannabeparents[i]);
  }
}

void GeneticAlgorithm::set(int min, int max, int children, int mutations) {
  Nmin = min;
  Nmax = max;
  Nchildren = children;
  Nmutations = mutations;
}

void GeneticAlgorithm::set_Nmin(int i) {
  Nmin = i;
}

void GeneticAlgorithm::set_Nmax(int i) {
  Nmax = i;
}

void GeneticAlgorithm::set_Nchildren(int i) {
  Nchildren = i;
}

void GeneticAlgorithm::set_Nmutations(int i) {
  Nmutations = i;
}

BP::BP_Vec<ECISet> GeneticAlgorithm::get_population() const {
  BP::BP_Vec<ECISet> out_pop;

  for(int i = 0; i < population.size(); i++) {
    out_pop.add(population[i]);
  }
  return out_pop;
};

int GeneticAlgorithm::get_Nmin() const {
  return Nmin;
}

int GeneticAlgorithm::get_Nmax() const {
  return Nmax;
}

int GeneticAlgorithm::get_Nchildren() const {
  return Nchildren;
}

int GeneticAlgorithm::get_Nmutations() const {
  return Nmutations;
}

int GeneticAlgorithm::get_best_index() const {
  int best;
  double min = 1.0e20;
  for(int i = 0; i < population.size(); i++) {
    if(population[i].get_cv() < min) {
      best = i;
      min = population[i].get_cv();
    }
  }

  return best;
}

void GeneticAlgorithm::run(MTRand &mtrand, Correlation &corr, EnergySet &nrg_set, bool print_steps) {
  run_ga(mtrand, corr, nrg_set, print_steps, 0, 0);
}

void GeneticAlgorithm::run_dir(MTRand &mtrand, Correlation &corr, EnergySet &nrg_set, bool print_steps) {
  run_ga(mtrand, corr, nrg_set,  print_steps, 1, 0);
}

void GeneticAlgorithm::run_dfs(long int Nstop, MTRand &mtrand, Correlation &corr, EnergySet &nrg_set, bool print_steps) {
  run_ga(mtrand, corr, nrg_set,  print_steps, 2, Nstop);
}

void GeneticAlgorithm::run_ga(MTRand &mtrand, Correlation &corr, EnergySet &nrg_set, bool print_steps, int min_option, long int Nstop) {
  //std::cout << "begin GeneticAlgorithm::run_ga() " << std::endl;
  int i;
  double fitness;
  bool singular;
  ECISet child;
  ECISet eci_best;
  BP::BP_GVec_Member<ECISet> *Mom;
  BP::BP_GVec_Member<ECISet> *Dad;
  BP::BP_GVec_Member<ECISet> *worst_parent;
  BP::BP_GVec_Member<ECISet> *new_addition;


  // add population to generator, make 'rate' = exp(1/cv) for weighting parent choice
  //std::cout << "add to generator" << std::endl;
  BP::BP_RVG_tree<ECISet> generator;
  for(i = 0; i < population.size(); i++) {
    population[i].fit(corr, nrg_set, singular);
    if(eci_best.get_cv() == UK || population[i].get_cv() < eci_best.get_cv())
      eci_best = population[i];
    //std::cout << " i: " << population[i].get_cv() << " bestsofar: " << eci_best.get_cv() << std::endl;
    //generator.add( population.member(i), exp(1.0/population[i].get_cv()));
    //generator.add( population.member(i), 1.0/population[i].get_cv());
    generator.add(population.member(i), 1.0);
  }

  // mate
  //std::cout << "run mate loop" << std::endl;
  fitness = max_cv(population, worst_parent);
  for(i = 0; i < Nchildren; i++) {
    //std::cout << "--pick parents" << std::endl;
    Mom = generator.pick(mtrand.randExc());
    do {
      Dad = generator.pick(mtrand.randExc());
    }
    while(Dad == Mom);

    //std::cout << "  mate parents" << std::endl;
    child = mate(Mom->get_val(), Dad->get_val(), mtrand);
    child.fit(corr, nrg_set, singular);

    //std::cout << "  check child fitness" << std::endl;




    //std::cout << "  child cv: " << child.get_cv() << std::endl;

    //std::cout << "  compare to pop" << std::endl;
    if(child.get_cv() != UK)
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

        if(eci_best.get_cv() == UK || child.get_cv() < eci_best.get_cv()) {
          eci_best = child;
        }

        new_addition = add_once(population, child);

        if(new_addition != NULL) {


          //std::cout << "    add to pop" << std::endl;
          //std::cout << child.get_bit_string() << "    GA: " << i << "  cv: " << child.get_cv() << " rms: " << child.get_rms() << std::endl;
          population.remove(worst_parent);
          //generator.add(new_addition, exp( 1.0/child.get_cv()));
          //generator.add(new_addition, 1.0/child.get_cv() );
          generator.add(new_addition, 1.0);
          fitness = max_cv(population, worst_parent);
        }
        else {
          std::cout << "not new_addition" << std::endl;
        }
      }


    std::cout << child.get_bit_string() << "  child: " << i << " Nclust: " << child.get_Nclust_on() << "  cv: " << child.get_cv() << " rms: " << child.get_rms() << " bestsofar: " << eci_best.get_cv() << " max_cv: " << fitness << std::endl;

  }

  //std::cout << "finish GeneticAlgorithm::run_ga() " << std::endl;

}

#endif // GeneticAlgorithm_CC
