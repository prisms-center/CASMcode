/*
 *  Minimize.hh
 */

#ifndef Minimize_HH
#define Minimize_HH

#include "EnergySet.hh"
#include "ECISet.hh"
#include "Correlation.hh"
#include "BP_Vec.hh"
#include <string>
#include <sstream>

class DirectMinStep {

public:

  const EnergySet *nrg;
  ECISet eci;
  const Correlation *corr;
  BP::BP_Vec<int> toggle;

  // for direct minimization
  int Nchoice;
  bool cont;  // continue?
  ECISetState best_state;

  DirectMinStep(const EnergySet &_nrg, const ECISet &_eci, const Correlation &_corr):
    nrg(&_nrg), eci(_eci), corr(&_corr) {
  }

  void run();
  static void *run_threaded(void *_this);
};


class DFSMinStep {

public:

  const EnergySet *nrg;
  ECISet eci;
  const Correlation *corr;
  const BP::BP_Vec<std::string> *bit_string_list;
  BP::BP_Vec<int> toggle;

  // for dfs minimization
  BP::BP_Vec<std::string> new_bit_string_list;
  BP::BP_Vec<ECISetState> improved_state;

  DFSMinStep(const EnergySet &_nrg,
             const ECISet &_eci,
             const Correlation &_corr,
             const BP::BP_Vec<std::string> &_bit_string_list):
    nrg(&_nrg), eci(_eci), corr(&_corr), bit_string_list(&_bit_string_list) {
  }

  void run();
  static void *run_threaded(void *_this);


};



class Minimize {

  const EnergySet *nrg;
  ECISet eci;
  const Correlation *corr;
  std::string sout;
  int Nthreads;

  int step;

  bool print_steps;
  int Nstop;

  int finished;

public:
  Minimize(const EnergySet &_nrg, const ECISet &_eci, const Correlation &_corr, int _Nthreads = 1, bool _print_steps = false):
    nrg(&_nrg), eci(_eci), corr(&_corr), Nthreads(_Nthreads), step(0), print_steps(_print_steps), finished(false) {
  }

  void direct();
  void dfs();

  void set_Nstop(int _Nstop) {
    Nstop = _Nstop;
  }

  static void *direct_threaded(void *_this) {
    ((Minimize *) _this)->direct();
    return NULL;
  }

  static void *dfs_threaded(void *_this) {
    ((Minimize *) _this)->dfs();
    return NULL;
  }

  const ECISet &get_eci() const {
    return eci;
  }

  const int &get_step() const {
    return step;
  }

  const std::string &get_sout() const {
    return sout;
  }

  bool is_finished() const {
    return finished;
  }

};

#endif // Minimize_HH
