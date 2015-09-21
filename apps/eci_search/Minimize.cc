/*
 *  Minimize.cc
 */

#ifndef Minimize_CC
#define Minimize_CC

#include "Minimize.hh"
#include "BP_ThreadPool.hh"

void DirectMinStep::run() {
  Nchoice = 0;
  cont = false;
  bool singular;
  best_state = eci.get_state();
  double last_cv = best_state.cv;

  for(int i = 0; i < toggle.size(); i++) {
    eci.toggle_clust(toggle[i]);

    // find fit/cv score
    eci.fit(*corr, *nrg, singular);

    if(eci.get_cv() < last_cv) {
      Nchoice++;
    }

    // keep the steepest descent option
    if(eci.get_cv() < best_state.cv) {
      cont = true;
      best_state = eci.get_state();

    }

    eci.toggle_clust(toggle[i]);
  }
}

void *DirectMinStep::run_threaded(void *_this) {
  ((DirectMinStep *) _this)->run();
  return NULL;
}

void DFSMinStep::run() {
  bool singular;
  double last_cv = eci.get_cv();
  improved_state.clear();
  new_bit_string_list.clear();

  for(int i = 0; i < toggle.size(); i++) {
    eci.toggle_clust(toggle[i]);

    if(!bit_string_list->contains(eci.get_bit_string())) {
      new_bit_string_list.add(eci.get_bit_string());

      // find fit/cv score
      eci.fit(*corr, *nrg, singular);

      if(eci.get_cv() < last_cv) {
        improved_state.add(eci.get_state());
      }
    }

    eci.toggle_clust(toggle[i]);
  }
}

void *DFSMinStep::run_threaded(void *_this) {
  ((DFSMinStep *) _this)->run();
  return NULL;
}



void Minimize::direct() {
  bool cont, singular;
  int i, j, Nchoice;

  eci.fit(*corr, *nrg, singular);
  ECISetState best_state = eci.get_state();

  // create a threadpool
  BP::ThreadPool pool(Nthreads);

  // we'll be toggling each eci on/off
  //   so to parallelize, break eciset into portions
  BP::BP_Vec<DirectMinStep> portion;
  for(i = 0; i < Nthreads; i++) {
    portion.add(DirectMinStep(*nrg, eci, *corr));
  }

  std::stringstream ss;

  ss << eci.get_bit_string() << "            init:" << " Nclust: " << eci.get_Nclust_on() << " cv: " << eci.get_cv() << " rms: " << eci.get_rms() << std::endl;
  if(print_steps) std::cout << ss.str();
  sout += ss.str();
  ss.str("");

  // minimization loop
  step = 1;
  do {
    // reset the eci to flip
    for(i = 0; i < portion.size(); i++)
      portion[i].toggle.clear();

    j = 0;
    // portion out the eci to flip
    for(i = 0; i < portion[0].eci.size(); i++) {
      if(portion[0].eci.toggle_allowed(i)) {
        portion[j].toggle.add(i);
        j++;
        j = j % Nthreads;
      }
    }

    // run the step
    for(i = 0; i < portion.size(); i++)
      pool.add_work(DirectMinStep::run_threaded, (void *) &portion[i]);
    pool.finish();

    cont = false;
    Nchoice = 0;
    // find the best state
    for(i = 0; i < portion.size(); i++)
      if(portion[i].best_state.cv < best_state.cv) {
        cont = true;
        best_state = portion[i].best_state;
        Nchoice += portion[i].Nchoice;
      }

    // set the best state
    if(cont) {
      for(i = 0; i < portion.size(); i++)
        portion[i].eci.set_state(best_state);
      step++;
    }

    ss << portion[0].eci.get_bit_string() << "      minimizing: " << step << " Nclust: " << portion[0].eci.get_Nclust_on() << " cv: " << portion[0].eci.get_cv() << " rms: " << portion[0].eci.get_rms() << " Nchoice: " << Nchoice << "\n";
    if(print_steps) std::cout << ss.str();
    sout += ss.str();
    ss.str("");


    eci.set_state(best_state);
  }
  while(cont);

  ss << eci.get_bit_string() << "            init:" << " Nclust: " << eci.get_Nclust_on() << " cv: " << eci.get_cv() << " rms: " << eci.get_rms() << std::endl;
  if(print_steps) std::cout << ss.str();
  sout += ss.str();
  ss.str("");

  finished = true;

}


void Minimize::dfs() {
  bool cont, singular;
  int i, j, Nchoice;
  unsigned long int index = 0;
  BP::BP_Vec<std::string> bit_string_list;
  BP::BP_Vec<ECISetState> queue;
  int count_since_last_best = 0;
  int bestsofar = 0;
  ECISetState min_state;


  eci.fit(*corr, *nrg, singular);
  ECISetState best_state = eci.get_state();
  //queue.add(best_state);
  bit_string_list.add(eci.get_bit_string());

  // create a threadpool
  BP::ThreadPool pool(Nthreads);

  // we'll be toggling each eci on/off
  //   so to parallelize, break eciset into portions
  BP::BP_Vec<DFSMinStep> portion;
  for(i = 0; i < Nthreads; i++) {
    portion.add(DFSMinStep(*nrg, eci, *corr, bit_string_list));
  }

  std::stringstream ss;

  ss << eci.get_bit_string() << "            init:" << " Nclust: " << eci.get_Nclust_on() << " cv: " << eci.get_cv() << " rms: " << eci.get_rms() << std::endl;
  if(print_steps) std::cout << ss.str();
  sout += ss.str();
  ss.str("");

  // minimization loop
  step = 0;
  do {
    // reset the eci to flip
    for(i = 0; i < portion.size(); i++)
      portion[i].toggle.clear();

    j = 0;
    // portion out the eci to flip
    for(i = 0; i < portion[0].eci.size(); i++) {
      if(portion[0].eci.toggle_allowed(i)) {
        //ss << "toggle i: " << i << "\n";
        portion[j].toggle.add(i);
        j++;
        j = j % Nthreads;
      }
    }
    //sout += ss.str();
    //ss.str("");

    // run the step
    for(i = 0; i < portion.size(); i++)
      pool.add_work(DFSMinStep::run_threaded, (void *) &portion[i]);
    pool.finish();

    Nchoice = 0;
    // add the new_bit_string_list's to bit_string_list
    // add the improved_state's to the queue
    // count the number of possible choices
    for(i = 0; i < portion.size(); i++) {
      bit_string_list.append(portion[i].new_bit_string_list);
      queue.append(portion[i].improved_state);
      Nchoice += portion[i].improved_state.size();
    }

    if(queue.size() > 0) {
      // find the best state in the queue
      min_state = min(queue, index);
      queue.remove(index);

      if(min_state < best_state) {
        best_state = min_state;
        bestsofar = step;
        count_since_last_best = 0;
      }

      for(i = 0; i < portion.size(); i++)
        portion[i].eci.set_state(min_state);

      ss << portion[0].eci.get_bit_string() << "      DFSchoice: " << step << " Nclust: " << portion[0].eci.get_Nclust_on() << " cv: " << portion[0].eci.get_cv() << " rms: " << portion[0].eci.get_rms() << " Nchoice: " << Nchoice << " listsize: " << queue.size() << " Count: " << count_since_last_best << "/" << Nstop << "  bestofDFS: " << bestsofar << "  best_cv: " << best_state.cv << std::endl;
      if(print_steps) std::cout << ss.str();
      sout += ss.str();
      ss.str("");
    }

    step++;
    count_since_last_best++;
    eci.set_state(best_state);


  }
  while((queue.size() > 0) && (count_since_last_best < Nstop || Nstop == 0));

  ss << eci.get_bit_string() << "            final:" << " Nclust: " << eci.get_Nclust_on() << " cv: " << eci.get_cv() << " rms: " << eci.get_rms() << std::endl;
  if(print_steps) std::cout << ss.str();
  sout += ss.str();
  ss.str("");

  finished = true;

}

#endif // Minimize_CC
