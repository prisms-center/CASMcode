/*
*  Population.cc
*/

#ifndef Population_CC
#define Population_CC

#include "jsonParser.hh"
#include "Population.hh"
#include "Functions.hh"
#include "BP_useful.hh"
#include <sstream>
#include <unistd.h>

std::string Population::format() const {
  return m_format;
}

void Population::populate_with_base() {
  population.clear();
  add(eci_base);
}

void Population::populate(int Npop) {

  // make population of Npop random ECISets

  population.clear();
  for(int i = 0; i < Npop; i++) {
    add();
  }
}

void Population::populate(const std::string &bit_string_filename) {

  m_format = get_format_from_ext(bit_string_filename);

  // make population from bit string file
  if(m_format == "text") {
    population.clear();
    BP::BP_Parse file(bit_string_filename);

    BP::BP_Vec<std::string> entry;

    do {
      entry = file.getline_string();
      if(entry.size() != 0) {
        // check if the first entry is a bit_string
        if(is_bitstring(entry[0])) {
          add(entry[0]);

        }

      }
    }
    while(file.eof() == false);
  }
  else if(m_format == "json") {

    population.clear();
    CASM::jsonParser json(bit_string_filename);

    for(int i = 0; i < json["population"].size(); i++) {
      add(json["population"][i]["eciset"].get<std::string>());
    }

  }
  else {
    std::cout << "Unexpected format option for bit strings population file reader" << std::endl;
    std::cout << "  Expected 'text' or 'json', but received: " << m_format << std::endl;
    exit(1);
  }
}

// Add a random ECISet to the population
void Population::add() {
  population.add(eci_base);
  population.last().set_data(nrg, corr);
  population.last().randomize(eci_base.get_Nmin(), mtrand);
};

// Add an ECISet with ECISetState to the population
void Population::add(const ECISetState &state) {
  population.add(eci_base);
  population.last().set_state(state);
  population.last().set_data(nrg, corr);

};

// Add an ECISet to the population with given bit_string
void Population::add(const std::string &bit_string) {
  population.add(eci_base);
  population.last().set_bit_string(bit_string);
  population.last().set_data(nrg, corr);

};

void Population::direct(int pthreads, int mthreads) {
  //std::cout << "begin Population::direct" << std::endl;

  // Direct minimization of the population
  std::cout << "\nBeginning direct minimization of " << population.size() << " samples." << std::endl << std::endl;

  if(population.size() == 0) {
    std::cout << "  Error. population size: " << population.size() << std::endl;
    exit(1);
  }

  // determine the number of threads
  int ncore = sysconf(_SC_NPROCESSORS_ONLN);

  if(pthreads == -1)
    pthreads = ncore < population.size() ? ncore : population.size() ;
  if(mthreads == -1)
    while(mthreads * pthreads < ncore)
      mthreads++;

  if(TEST) std::cout << "pthreads: " << pthreads << " mthreads: " << mthreads << std::endl << std::endl;

  // create a ThreadPool
  BP::ThreadPool pool(pthreads);

  // create minimization objects for each ECISet in the population & run the minimization
  BP::BP_Vec<Minimize> minimization;
  BP::BP_Vec<bool> completed;
  for(int i = 0; i < population.size(); i++) {
    completed.add(false);
    minimization.add(Minimize(nrg, population[i], corr, mthreads, false));
    pool.add_work(Minimize::direct_threaded, (void *) &minimization[i]);
  }

  // while the minimizations are ongoing, print results of completed fits
  do {
    sleep(1);

    std::cout << minimization_status(minimization, completed);
    std::cout << flush;
  }
  while(!pool.is_finished());

  std::cout << minimization_status(minimization, completed);
  std::cout << flush;

  // set the final, minimized population
  for(int i = 0; i < minimization.size(); i++) {
    population[i].set_state(minimization[i].get_eci().get_state());
  }

  //std::cout << "finish Population::direct" << std::endl;
}

void Population::dfs(int Nstop, int pthreads, int mthreads) {
  //std::cout << "begin Population::dfs" << std::endl;

  // Depth first search minimization of the population
  //   At each step, add to a queue all 1-eci perturbations of the ECISet which lower the cv score
  //   From the queue, choose the ECISet with the lowest cv score, and repeat
  //   Continue until a new overall best cv score has not been found for Nstop steps
  std::cout << "\nBeginning depth first search minimization of " << population.size() << " samples." << std::endl;
  std::cout << "  Nstop: " << Nstop << std::endl << std::endl;

  if(population.size() == 0) {
    std::cout << "  Error. population size: " << population.size() << std::endl;
    exit(1);
  }

  // determine the number of threads
  int ncore = sysconf(_SC_NPROCESSORS_ONLN);

  if(pthreads == -1)
    pthreads = ncore < population.size() ? ncore : population.size();
  if(mthreads == -1)
    while(mthreads * pthreads < ncore)
      mthreads++;

  if(TEST) std::cout << "pthreads: " << pthreads << " mthreads: " << mthreads << std::endl;

  // create a ThreadPool
  BP::ThreadPool pool(pthreads);

  // create minimization objects for each ECISet in the population & run the minimization
  BP::BP_Vec<Minimize> minimization;
  BP::BP_Vec<bool> completed;
  for(int i = 0; i < population.size(); i++) {
    completed.add(false);
    minimization.add(Minimize(nrg, population[i], corr, mthreads, false));
    minimization[i].set_Nstop(Nstop);
    pool.add_work(Minimize::dfs_threaded, (void *) &minimization[i]);
  }

  // while the minimizations are ongoing, print results of completed fits
  do {
    sleep(1);

    std::cout << minimization_status(minimization, completed);
    std::cout << flush;
  }
  while(!pool.is_finished());

  std::cout << minimization_status(minimization, completed);
  std::cout << flush;

  // set the final, minimized population
  for(int i = 0; i < minimization.size(); i++) {
    population[i].set_state(minimization[i].get_eci().get_state());
  }

  //std::cout << "finish Population::dfs" << std::endl;
}

void Population::ga(int Ngen, int Nmut, int pthreads, int mode, int Nstop) {
  //std::cout << "begin Population::ga" << std::endl;

  if(mode != 0 && mode != 1 && mode != 2) {
    std::cout << "\nError in Population::ga(). mode = " << mode << std::endl;
    std::cout << "Options are:" << std::endl;
    std::cout << "  0: normal" << std::endl;
    std::cout << "  1: direct minimization of each child" << std::endl;
    std::cout << "  2: dfs minimization of each child" << std::endl << std::endl;
  }
  else if(mode == 0) {
    std::cout << "\nBeginning genetic algorithm: " << std::endl;
    std::cout << "  population size: " << population.size() << std::endl;
    std::cout << "  number of generations: " << Ngen << std::endl;
    std::cout << "  avg number of mutations per child: " << Nmut << std::endl << std::endl;
  }
  else if(mode == 1) {
    std::cout << "\nBeginning genetic algorithm, with direct minimization: " << std::endl;
    std::cout << "  population size: " << population.size() << std::endl;
    std::cout << "  number of generations: " << Ngen << std::endl;
    std::cout << "  avg number of mutations per child: " << Nmut << std::endl << std::endl;
  }
  else if(mode == 2) {
    std::cout << "\nBeginning genetic algorithm, with depth first search minimization: " << std::endl;
    std::cout << "  population size: " << population.size() << std::endl;
    std::cout << "  number of generations: " << Ngen << std::endl;
    std::cout << "  avg number of mutations per child: " << Nmut << std::endl << std::endl;
    std::cout << "  Nstop: " << Nstop << std::endl << std::endl;
  }

  if(population.size() == 0) {
    std::cout << "  Error. population size: " << population.size() << std::endl;
    exit(1);
  }

  // determine the number of threads
  int ncore = sysconf(_SC_NPROCESSORS_ONLN);
  int ga_pthreads;

  if(pthreads <= 0)
    ga_pthreads = ncore;

  if(TEST) std::cout << "pthreads: " << pthreads << std::endl;

  // create ThreadPool
  BP::ThreadPool pool(ga_pthreads);

  // create gene pool
  BP::BP_Vec<ECISetState> gene_pool;
  for(int i = 0; i < population.size(); i++) {
    gene_pool.add(population[i].get_state());
  }
  std::cout << "---------------------------------------------" << std::endl;
  std::cout << "Initial gene pool:" << std::endl;
  std::cout << gene_pool_status(gene_pool);
  std::cout << "---------------------------------------------" << std::endl;

  for(int gen = 0; gen < Ngen; gen++) {

    // mate
    for(int i = 0; i < population.size(); i++) {
      int count = 0;
      bool found;
      do {
        mate(gene_pool, population[i], Nmut, mtrand);

        // only allow unique children
        found = false;
        for(int j = 0; j < i; j++) {
          if(population[j].get_bit_string() == population[i].get_bit_string()) {
            found = true;
            break;
          }
        }
        count++;
      }
      while(found);
    }

    // fit
    if(mode == 1) {
      std::cout << "\nDirect minimization of children, generation: " << gen << std::endl;
      direct(pthreads);
    }
    else if(mode == 2) {
      std::cout << "\nDepth first search minimization of children, generation: " << gen << std::endl;
      dfs(Nstop, pthreads);
    }
    else {
      for(int i = 0; i < population.size(); i++) {
        pool.add_work(ECISet::fit_threaded, (void *) &population[i]);
      }
      pool.finish();
    }

    // prune to best Npop unique ECISets of last 2 generations
    for(int i = 0; i < population.size(); i++) {
      BP::add_once(gene_pool, population[i].get_state());
    }
    while(gene_pool.size() > population.size()) {
      gene_pool.remove(BP::max_index(gene_pool));
    }

    // show current gene pool status periodically
    if(gen % 10 == 0) {
      std::cout << "\n---------------------------------------------" << std::endl;
      std::cout << "Current gene pool, generation: " << gen << std::endl;
      std::cout << gene_pool_status(gene_pool);
      std::cout << "---------------------------------------------" << std::endl;
    }
  }

  // set the final population
  for(int i = 0; i < population.size(); i++) {
    population[i].set_state(gene_pool[i]);
  }

  std::cout << "\nFinal Population:" << std::endl;
  std::cout << population_status();

  //std::cout << "finish Population::ga" << std::endl;

}

void Population::ga_dir(int Ngen, int Nmut, int pthreads) {
  ga(Ngen, Nmut, pthreads, 1);
}

void Population::ga_dfs(int Nstop, int Ngen, int Nmut, int pthreads) {
  ga(Ngen, Nmut, pthreads, 2, Nstop);
}

void Population::calc(int pthreads) {
  //std::cout << "begin Population::calc" << std::endl;

  // Calculate cv score for every ECISet in the population
  std::cout << "\nBeginning calculation of cv score for " << population.size() << " samples." << std::endl << std::endl;

  if(population.size() == 0) {
    std::cout << "  Error. population size: " << population.size() << std::endl;
    exit(1);
  }

  // determine the number of threads
  int ncore = sysconf(_SC_NPROCESSORS_ONLN);

  if(pthreads == -1)
    pthreads = ncore < population.size() ? ncore : population.size();

  if(TEST) std::cout << "pthreads: " << pthreads << std::endl;

  // create Threadpool
  BP::ThreadPool pool(pthreads);

  // fit each ECISet in the population
  for(int i = 0; i < population.size(); i++) {
    pool.add_work(ECISet::fit_threaded, &population[i]);
  }
  pool.finish();

  std::cout << population_status();

  //std::cout << "finish Population::calc" << std::endl;
}

void Population::calc_all(int N, int pthreads) {
  //std::cout << "begin Population::calc_all" << std::endl;

  // Find the optimal cv score considering all combinations of ECISets with N eci
  // There may be very many possible combinations, so don't set population.size() to number of combinations,
  //   instead fit Npop=(pthreads*10) combinations at a time. Find the best cv score.
  // Final population is size=1, with the best cv score.
  std::cout << "Beginning calculation of cv score for all combinations with " << N << " eci." << std::endl << std::endl;

  // determine the number of threads
  int ncore = sysconf(_SC_NPROCESSORS_ONLN);

  if(pthreads <= 0)
    pthreads = ncore;

  if(TEST) std::cout << "pthreads: " << pthreads << std::endl;

  // create the ThreadPool
  BP::ThreadPool pool(pthreads);

  bool singular;
  int bestsofar = 0;
  int zero = 0;
  ECISetState best_state;

  // create BP_Comb object to generate all combinations of ECISets with N eci
  BP::BP_Comb combs(eci_base.size(), N);

  std::cout << "This calculates all ecisets with " << N << " eci." << std::endl;
  std::cout << "   For N = " << N << " and ECISet size = " << eci_base.size() << " that will be " << combs.total_combs() << " ecisets." << std::endl << std::endl;



  // Create Npop ECISets
  int Npop = 10 * pthreads;
  int max = Npop;
  BP::BP_Vec<int> count;
  for(int i = 0; i < Npop; i++) {
    count.add();
    add(eci_base);
  }

  do {

    // Fit the current Npop ECISets
    for(int i = 0; i < Npop; i++) {
      population[i] = combs;
      count[i] = combs.get_count();
      pool.add_work(ECISet::fit_threaded, &population[i]);
      combs.increment();
      if(combs.complete()) {
        max = i;
        break;
      }
    }
    pool.finish();

    // Track the best cv score
    for(int i = 0; i < max; i++) {
      if(best_state.cv == UK || population[i].get_cv() < best_state.cv) {
        best_state = population[i].get_state();
        bestsofar = count[i];
      }

      if(!population[i].get_singular()) {
        std::cout << population[i].get_bit_string()
                  << "    i:" << std::setw(12) << count[i] << " "
                  << "   cv:" << std::setw(12) << population[i].get_cv() << " "
                  << "   rms:" << std::setw(12) << population[i].get_rms() << " " << std::endl;
      }
      else {
        std::cout << "error, singular" << std::endl;
      }
    }
  }
  while(combs.complete() == false);

  // Final population is just the 1 ECISet with the best cv score (no checking for equality...)
  population.clear();
  add(best_state);

  //std::cout << population[zero].get_bit_string() << "       BestSoFar: " << bestsofar << " Nclust: " << population[zero].get_Nclust_on() << " cv: " << population[zero].get_cv() << " rms: " << population[zero].get_rms() << std::endl;
  std::cout << population[zero].get_bit_string()
            << " Best:" << std::setw(12) << count[zero] << " "
            << "   cv:" << std::setw(12) << population[zero].get_cv() << " "
            << "   rms:" << std::setw(12) << population[zero].get_rms() << " " << std::endl;

  //std::cout << "finish Population::calc_all" << std::endl;
}

void Population::calc_details(int index, std::string out_format, double hulltol) {
  // Calculate and write files with details about ECISet at position 'index' in population
  //
  // Writes:
  //  corr.covar
  //  eci.in
  //  eci.out
  //  hull
  //  hull.clex
  //  energy.clex
  //  hull.clex_of_DFT_hull
  //  below.hull
  //  clex_results_X eps images

  if(out_format == "default")
    out_format = nrg.format();

  if(index >= population.size()) {
    std::cout << "Error in Population::calc_details(int index). index >= population.size()." << std::endl;
    std::cout << "  index: " << index << "  population.size(): " << population.size() << std::endl;
    exit(1);
  }

  ECISet &eci_in = population[index];
  EnergySet &DFT_nrg = nrg;
  bool singular;

  eci_in.fit(corr, DFT_nrg, singular);

  if(!singular) {
    std::cout << eci_in.get_bit_string() << std::endl;
    std::cout << "cv: " << eci_in.get_cv() << " rms: " << eci_in.get_rms() << std::endl;
  }
  else {
    std::cout << "error, singular" << std::endl;
  }

  std::cout << std::endl << std::endl;

  // write corr.covar
  corr.write_covar(eci_in, DFT_nrg, "corr.covar", out_format);

  // write eci.in
  eci_in.write_ECIin("eci.in", out_format);
  std::cout << "Wrote 'eci.in'" << std::endl;

  // write eci.out
  eci_in.write_ECIout("eci.out", DFT_nrg, out_format);
  std::cout << "Wrote 'eci.out'" << std::endl;

  //std::cout << "Hull of: " << energy_filename << std::endl;
  //DFT_nrg.write_hull(std::cout);
  std::cout << std::endl << "Calculating hull of: " << DFT_nrg.get_name() << std::endl;
  DFT_nrg.calc_hull(true, hulltol);
  DFT_nrg.write_hull("hull", out_format);
  std::cout << "Wrote 'hull'" << std::endl;


  if(eci_in.get_cv() != 1e20) {

    // calc energy_clex
    EnergySet energy_clex = DFT_nrg;
    energy_clex.calc_clex(corr, eci_in);

    // write energy_clex
    std::cout << std::endl << "Calculating hull of: energy.clex " << std::endl;
    energy_clex.calc_hull(true, hulltol);
    //energy_clex.write_hull(std::cout);
    energy_clex.write_hull("hull.clex", out_format);
    std::cout << "Wrote 'hull.clex'" << std::endl;
    energy_clex.write("energy.clex", out_format);
    std::cout << "Wrote 'energy.clex'" << std::endl;

    std::cout << std::endl;

    unsigned long int count;
    double rms;

    rms = energy_clex.calc_rms(DFT_nrg, 1, -1, count);
    std::cout << "weighted total rms: " <<  rms <<  "    #structures: " << count << std::endl;

    rms = energy_clex.calc_rms(DFT_nrg, 0, -1, count);
    std::cout << "non-weighted total rms: " <<  rms << "    #structures: " << count << std::endl;

    rms = energy_clex.calc_rms(DFT_nrg, 0, 0.2, count);
    std::cout << "non-weighted rms within 0.2 of hull: " << rms << "    #structures: " << count << std::endl;

    rms = energy_clex.calc_rms(DFT_nrg, 0, 0.1, count);
    std::cout << "non-weighted rms within 0.1 of hull: " << rms << "    #structures: " << count << std::endl;

    rms = energy_clex.calc_rms(DFT_nrg, 0, 0.05, count);
    std::cout << "non-weighted rms within 0.05 of hull: " << rms << "    #structures: " << count << std::endl;

    rms = energy_clex.calc_rms(DFT_nrg, 0, 0.01, count);
    std::cout << "non-weighted rms within 0.01 of hull: " << rms << "    #structures: " << count << std::endl;

    rms = energy_clex.calc_rms(DFT_nrg, 0, 0.005, count);
    std::cout << "non-weighted rms within 0.005 of hull: " << rms << "    #structures: " << count << std::endl;

    rms = energy_clex.calc_rms(DFT_nrg, 0, 0.001, count);
    std::cout << "non-weighted rms within 0.001 of hull: " << rms << "    #structures: " << count << std::endl;

    EnergySet energy_clex_of_DFT_hull;

    if(DFT_nrg.is_hull_found()) {
      energy_clex_of_DFT_hull.set(energy_clex.get(DFT_nrg.get_hull_indices()));
      std::cout << std::endl << "Calculating hull of: clex_of_DFT_hull" << std::endl;
      energy_clex_of_DFT_hull.calc_hull(true, hulltol);
      energy_clex_of_DFT_hull.write_hull("hull.clex_of_DFT_hull", out_format);
      std::cout << "Wrote 'hull.clex_of_DFT_hull'" << std::endl;
      energy_clex_of_DFT_hull.write_below_hull("below.hull", energy_clex.get(DFT_nrg.get_non_hull_indices()), out_format);
      std::cout << "Wrote 'below.hull'" << std::endl;
      //energy_clex_of_DFT_hull.plot(plot, "DFT_hull_clex", "green", 1, 0.5, 0, 2, true, i, j);
    }

    std::cout << std::endl;
    for(int i = -1; i < (int) DFT_nrg.get_concentration_size(); i++) {
      for(int j = i + 1; j < (int) DFT_nrg.get_concentration_size(); j++) {
        // plot energy and energy.clex
        BP::BP_Plot plot;
        //plot( BP::BP_Plot &plot, std::string label, std::string s, int line_style, double line_width, int point_style, double point_size, bool point_face)
        DFT_nrg.plot(plot, "DFT", "red", 0, 1, 0, 5, false, i, j);
        energy_clex.plot(plot, "CLEX", "blue", 0, 0.5, 0, 3, true, i, j);

        // write below.hull
        if(DFT_nrg.is_hull_found()) {
          //EnergySet energy_clex_of_DFT_hull( energy_clex.get(DFT_nrg.get_hull_indices()));
          //std::cout << "Calculating hull of: clex_of_DFT_hull" << std::endl;
          //energy_clex_of_DFT_hull.calc_hull(true);
          //energy_clex_of_DFT_hull.write_hull("hull.clex_of_DFT_hull", out_format);
          //std::cout << "  wrote 'hull.clex_of_DFT_hull'" << std::endl;
          //energy_clex_of_DFT_hull.write_below_hull( "below.hull", energy_clex.get(DFT_nrg.get_non_hull_indices()), out_format);
          energy_clex_of_DFT_hull.plot(plot, "DFT_hull_clex", "green", 1, 0.5, 0, 2, true, i, j);
        }

        std::string si, sj;
        if(i == -1) si = "nrg";
        else si = "x" + BP::itos(i);

        if(j == -1) sj = "nrg";
        else sj = "x" + BP::itos(j);



        plot.write("clex_results_" + si + "_vs_" + sj);
        std::cout << "Wrote '" << "clex_results_" + si + "_vs_" + sj << "'" << std::endl;
      }
    }

  }

  //eci_in.check_cv(corr, DFT_nrg, singular);
}

void Population::calc_ecistats(int pthreads) {

  // Print statistics for a population of ECI

  calc(pthreads);

  BP::BP_Vec<double> eci_value_list;
  BP::BP_Vec<double> eci_nonzero_value_list;
  int zero = 0;

  std::cout << std::endl;
  std::cout << std::setw(8) << "label" << " "
            << std::setw(12) << "frac_on" << " "
            << std::setw(14) << "mean(nonzero)" << " "
            << std::setw(14) << "min(nonzero)" << " "
            << std::setw(14) << "max(nonzero)" << " "
            << std::setw(14) << "rms(nonzero)" << " "
            << std::setw(12) << "mean" << " "
            << std::setw(12) << "min" << " "
            << std::setw(12) << "max" << " "
            << std::setw(12) << "rms" << " " << std::endl;
  for(int i = 0; i < population[zero].size(); i++) {
    eci_value_list.clear();
    eci_nonzero_value_list.clear();
    for(int j = 0; j < population.size(); j++) {
      eci_value_list.add(population[j].get_value(i));
      if(population[j].get_weight(i)) {
        eci_nonzero_value_list.add(population[j].get_value(i));
      }
    }

    std::cout << std::setw(8) << i << " ";
    std::cout << std::setw(12) << setprecision(6) << (1.0 * eci_nonzero_value_list.size()) / (1.0 * eci_value_list.size()) << " ";
    if(eci_nonzero_value_list.size() != 0) {
      std::cout << std::setw(14) << mean(eci_nonzero_value_list) << " ";
      std::cout << std::setw(14) << min(eci_nonzero_value_list) << " ";
      std::cout << std::setw(14) << max(eci_nonzero_value_list) << " ";
      std::cout << std::setw(14) << rms(eci_nonzero_value_list) << " ";
      std::cout << std::setw(12) << mean(eci_value_list) << " ";
      std::cout << std::setw(12) << min(eci_value_list) << " ";
      std::cout << std::setw(12) << max(eci_value_list) << " ";
      std::cout << std::setw(12) << rms(eci_value_list) << " " << "\n";
    }
    else {
      std::cout << std::setw(14) << "-" << " ";
      std::cout << std::setw(14) << "-" << " ";
      std::cout << std::setw(14) << "-" << " ";
      std::cout << std::setw(14) << "-" << " ";
      std::cout << std::setw(12) << "-" << " ";
      std::cout << std::setw(12) << "-" << " ";
      std::cout << std::setw(12) << "-" << " ";
      std::cout << std::setw(12) << "-" << " ""\n";
    }
  }
  std::cout << std::endl;
}
std::string Population::population_status() const {

  // return a string containing a summary of the population's Nclust, cv, and rms.
  //   highlights the best cv score with '***'

  double best_cv;
  int bestsofar = 0;
  best_cv = population[bestsofar].get_cv();

  stringstream ss;
  BP::BP_Vec<std::string> slist;

  for(int i = 0; i < population.size(); i++) {
    ss << population[i].get_bit_string() << "  set: " << std::setw(6) << std::left << i << "  Nclust: " << std::setw(6) << std::left << population[i].get_Nclust_on() << "  cv: " << std::setw(12) << population[i].get_cv() << "  rms: " << std::setw(12) << population[i].get_rms();
    slist.add(ss.str());
    if(population[i].get_cv() != UK) {
      if(population[i].get_cv() < best_cv) {
        best_cv = population[i].get_cv();
        bestsofar = i;
      }
    }
    ss.str("");
  }

  std::string result = "";
  for(int i = 0; i < population.size(); i++) {
    if(i == bestsofar) {
      slist[i] += " *** ";
    }

    result += slist[i] + "\n";
  }

  //ss << "\n";
  //ss << "BestSoFar: \n" << std::setw(6) << std::right << bestsofar << ": " << minimization[bestsofar].get_eci().get_bit_string() << "    step: " << std::setw(6) << std::left << minimization[bestsofar].get_step() << "  cv: " << std::setw(12) << minimization[bestsofar].get_eci().get_cv() << "  rms: " << std::setw(12) << minimization[bestsofar].get_eci().get_rms() << "\n";
  //result += ss.str();

  return result;
}

// Write population to file
void Population::write(std::string filename, std::string format) const {

  if(format == "default")
    format = m_format;

  if(format == "text") {
    BP::BP_Write file(rm_json_ext(filename));
    file.newfile();

    for(int i = 0; i < population.size(); i++) {
      file << population[i].get_bit_string();
      file << "  cv: " << population[i].get_cv();
      file << "  rms: " << population[i].get_rms() << "\n";
    }
    std::cout << "Wrote: " << rm_json_ext(filename) << std::endl;
  }
  else if(format == "json") {
    CASM::jsonParser json;
    json.put_obj();
    json["population"].put_array();
    for(int i = 0; i < population.size(); i++) {
      json["population"].push_back(CASM::jsonParser::object());
      json["population"][i]["eciset"] = population[i].get_bit_string();
      json["population"][i]["wcv"] = population[i].get_cv();
      json["population"][i]["rms"] = population[i].get_rms();
    }
    json.write(json_ext(filename));
    std::cout << "Wrote: " << json_ext(filename) << std::endl;
  }
  else {
    std::cout << "Unexpected format option for Population::write" << std::endl;
    std::cout << "  Expected 'text' or 'json', but received: " << format << std::endl;
    exit(1);
  }
}


// Add an ECISet to the population
void Population::add(const ECISet &eci) {
  population.add(eci);
  population.last().set_data(nrg, corr);

};

std::string Population::minimization_status(const BP::BP_Vec<Minimize> &minimization, BP::BP_Vec<bool> &completed) const {
  // during the course of a minimization, return a string containing newly completed fits
  //   highlights the best cv score with '*C*'

  double best_cv;
  int bestsofar;
  best_cv = minimization[0].get_eci().get_cv();
  bestsofar = 0;

  stringstream ss;
  BP::BP_Vec<std::string> slist;

  for(int i = 0; i < minimization.size(); i++) {
    //ss << minimization[i].get_eci().get_bit_string() << "  set: " << std::setw(6) << std::left << i << "  Nclust: " << std::setw(6) << std::left << minimization[i].get_eci().get_Nclust_on() << "  step: " << std::setw(6) << std::left << minimization[i].get_step() << "  cv: " << std::setw(12) << minimization[i].get_eci().get_cv() << "  rms: " << std::setw(12) << minimization[i].get_eci().get_rms();
    ss << minimization[i].get_eci().get_bit_string() << "  Nclust: " << std::setw(6) << std::left << minimization[i].get_eci().get_Nclust_on() << "  steps: " << std::setw(6) << std::left << minimization[i].get_step() << "  cv: " << std::setw(12) << minimization[i].get_eci().get_cv() << "  rms: " << std::setw(12) << minimization[i].get_eci().get_rms();
    slist.add(ss.str());
    if(minimization[i].get_eci().get_cv() != UK) {
      if(minimization[i].get_eci().get_cv() < best_cv) {
        best_cv = minimization[i].get_eci().get_cv();
        bestsofar = i;
      }
    }
    ss.str("");
  }

  std::string result = "";
  for(int i = 0; i < minimization.size(); i++) {
    if(!minimization[i].is_finished()) {
      /*(if(minimization[i].get_eci().get_cv() == UK) {
        slist[i] += "  Q  ";
        result += slist[i] + "\n";
      }
      else if(i != bestsofar) {
        slist[i] += "  R  ";
        result += slist[i] + "\n";
      }
      else {
        slist[i] += " *R* ";
        result += slist[i] + "\n";
      }*/
    }
    else if(!completed[i]) {
      if(i != bestsofar) {
        slist[i] += "  C  ";
        result += slist[i] + "\n";
      }
      else {
        slist[i] += " *C* ";
        result += slist[i] + "\n";
      }
      completed[i] = true;
    }
    //result += slist[i] + "\n";
  }

  //ss << "\n";
  //ss << "BestSoFar: \n" << std::setw(6) << std::right << bestsofar << ": " << minimization[bestsofar].get_eci().get_bit_string() << "    step: " << std::setw(6) << std::left << minimization[bestsofar].get_step() << "  cv: " << std::setw(12) << minimization[bestsofar].get_eci().get_cv() << "  rms: " << std::setw(12) << minimization[bestsofar].get_eci().get_rms() << "\n";
  //result += ss.str();

  return result;
}

std::string Population::gene_pool_status(const BP::BP_Vec<ECISetState> &gene_pool) const {

  // returns a string containing details about the curent gene_pool: Nclust, cv, rms
  //   highlights the best cv score with '***'

  double best_cv;
  int bestsofar = 0;
  best_cv = gene_pool[bestsofar].cv;

  stringstream ss;
  BP::BP_Vec<std::string> slist;

  for(int i = 0; i < gene_pool.size(); i++) {
    ss << gene_pool[i].bit_string << "  set: " << std::setw(6) << std::left << i << "  Nclust: " << std::setw(6) << std::left << gene_pool[i].Nclust << "  cv: " << std::setw(12) << gene_pool[i].cv << "  rms: " << std::setw(12) << gene_pool[i].rms;
    slist.add(ss.str());
    if(gene_pool[i].cv != UK) {
      if(gene_pool[i].cv < best_cv) {
        best_cv = gene_pool[i].cv;
        bestsofar = i;
      }
    }
    ss.str("");
  }

  std::string result = "";
  for(int i = 0; i < gene_pool.size(); i++) {
    if(i == bestsofar) {
      slist[i] += " *** ";
    }

    result += slist[i] + "\n";
  }

  //ss << "\n";
  //ss << "BestSoFar: \n" << std::setw(6) << std::right << bestsofar << ": " << minimization[bestsofar].get_eci().get_bit_string() << "    step: " << std::setw(6) << std::left << minimization[bestsofar].get_step() << "  cv: " << std::setw(12) << minimization[bestsofar].get_eci().get_cv() << "  rms: " << std::setw(12) << minimization[bestsofar].get_eci().get_rms() << "\n";
  //result += ss.str();

  return result;
}


#endif // Population_CC
