/*
 *  eci_search.cpp
 *
 *
 *  Created by Brian Puchala on 4/17/2012.
 *	All rights reserved.
 *  eci optimization search code
 *
 */

#include "jsonParser.cc"
#include "BP_Vec.hh"

namespace CASM {
  template<typename T>
  CASM::jsonParser &to_json(const BP::BP_Vec<T> &value, CASM::jsonParser &json) {
    json.put_array();
    for(int i = 0; i < value.size(); i++)
      json.push_back(value[i]);
    return json;
  }

  /// This requires that 'T::T()' exists, if not, you must do this by hand
  template<typename T>
  void from_json(BP::BP_Vec<T> &value, const CASM::jsonParser &json) {
    try {
      value.capacity(json.size());
      for(int i = 0; i < json.size(); i++)
        value.add(json[i].get<T>());
    }
    catch(...) {
      /// re-throw exceptions
      throw;
    }
  }
}

#include "BP_basic.cc"
#include "BP_Dir.cc"
#include "BP_Vec.cc"
#include "BP_GVec.cc"
#include "BP_Plot.cc"
#include "BP_Geo.cc"
#include "BP_Parse.cc"
#include "BP_StopWatch.cc"
#include "BP_ThreadPool.cc"
#include "Correlation.cc"
#include "ECISet.cc"
#include "EnergySet.cc"
#include "GeneticAlgorithm.cc"
#include "Functions.cc"
#include "Minimize.cc"
#include "Population.cc"

void print_calc_man() {
  std::cout << "  eci_search -calc energy eci.in corr.in [bitstring | bitstring_file]" << std::endl;
}
void print_ecistats_man() {
  std::cout << "  eci_search -ecistats energy eci.in corr.in population_file" << std::endl;
}
void print_calc_all_man() {
  std::cout << "  eci_search -calc_all N energy eci.in corr.in" << std::endl;
}
void print_calc_directmin_man() {
  std::cout << "  eci_search -calc_directmin Nrand Nmin Nmax energy eci.in corr.in [population_file]" << std::endl;
}
void print_calc_dfsmin_man() {
  std::cout << "  eci_search -calc_dfsmin Nrand Nstop Nmin Nmax energy eci.in corr.in [population_file]" << std::endl;
}
void print_calc_ga_man() {
  std::cout << "  eci_search -calc_ga Npop Nmin Nmax Ngen Nmut energy eci.in corr.in [population_file]" << std::endl;
}
void print_calc_ga_dir_man() {
  std::cout << "  eci_search -calc_ga_dir Npop Nmin Nmax Ngen Nmut energy eci.in corr.in [population_file]" << std::endl;
}
void print_calc_ga_dfs_man() {
  std::cout << "  eci_search -calc_ga_dfs Npop Nmin Nmax Ngen Nmut Nstop energy eci.in corr.in [population_file]" << std::endl;
}
void print_weight_nrg_man() {
  std::cout << "  eci_search -weight_nrg A B kT energy weighted_energy" << std::endl;
}
void print_weight_emin_man() {
  std::cout << "  eci_search -weight_emin A B kT energy weighted_energy" << std::endl;
}
void print_weight_eref_man() {
  std::cout << "  eci_search -weight_eref A B E0 kT energy weighted_energy" << std::endl;
}
void print_convert_man() {
  std::cout << "  eci_search -convert-to-text energy.json eci.in.json corr.in.json" << std::endl;
  std::cout << "  eci_search -convert-to-json energy eci.in corr.in" << std::endl;
  std::cout << "  eci_search -convert-energy-to-text energy.json [...]" << std::endl;
  std::cout << "  eci_search -convert-energy-to-json energy [...]" << std::endl;
  std::cout << "  eci_search -convert-eci-to-text eci.in.json [...]" << std::endl;
  std::cout << "  eci_search -convert-eci-to-json eci.in [...]" << std::endl;
  std::cout << "  eci_search -convert-corr-to-text corr.in.json [...]" << std::endl;
  std::cout << "  eci_search -convert-corr-to-json corr.in [...]" << std::endl;
}
void print_calc_cs_fpc_man() {
  std::cout << "  eci_search -calc_cs_fpc energy eci.in corr.in mu" << std::endl;
}
void print_calc_cs_bi_man() {
  std::cout << "  eci_search -calc_cs_bi energy eci.in corr.in mu" << std::endl;
}

void print_calc_full_man() {
  std::cout << "  eci_search -calc energy eci.in corr.in [bitstring | bitstring_file]" << std::endl;
  std::cout << "      This calculates the eci, cv score, and rms score for a given eciset.   " << std::endl;
  std::cout << "      The eciset can be given as a bitstring ('11001000100...'), or as a     " << std::endl;
  std::cout << "      file which contains a bitstring. If neither of these is given, then    " << std::endl;
  std::cout << "      the eciset is read from eci.in.                                        " << std::endl;
  std::cout << "                                                                             " << std::endl;
  std::cout << "      This function then writes the files:  eci.in                           " << std::endl;
  std::cout << "                                            eci.out                          " << std::endl;
  std::cout << "                                            hull                             " << std::endl;
  std::cout << "                                            hull.clex                        " << std::endl;
  std::cout << "                                            below.hull                       " << std::endl;
  std::cout << "                                            energy.clex                      " << std::endl;
  std::cout << "                                            clex_results_i.eps               " << std::endl;
  std::cout << "                                                                             " << std::endl;
  std::cout << "      The file 'eci.in' is written with the weights for cluster basis functions   " << std::endl;
  std::cout << "      included in the eciset set to 1. The file 'eci.out' contains the calculated " << std::endl;
  std::cout << "      eci. The file 'hull' contains a list of the structures which make up the    " << std::endl;
  std::cout << "      convex hull of the structures included in the 'energy' file.  The file " << std::endl;
  std::cout << "      'hull.clex' contains a list of the structures which make up the convex " << std::endl;
  std::cout << "      hull of the cluster expanded energies, which are listed in the file    " << std::endl;
  std::cout << "      'energy.clex'. The file(s) 'clex_results_i_vs_j.eps' contains a 2d plot" << std::endl;
  std::cout << "       of 'energy' as red open circles, 'hull' as a red line, 'energy.clex'" << std::endl;
  std::cout << "      as blue filled circles, and 'hull.clex' as a blue line. Additionally a    " << std::endl;
  std::cout << "      green dashed line with green filled circles is drawn for the cluster   " << std::endl;
  std::cout << "      expanded energies of the structures which are listed in 'hull'. One    " << std::endl;
  std::cout << "      'clex_results_i_vs_j.eps' file is printed for 2d combination of energy " << std::endl;
  std::cout << "      and concentration, with 'i' corresponding to which energy/concentration" << std::endl;
  std::cout << "      is plotted along the y-axis, and 'j' corresponding to which            " << std::endl;
  std::cout << "      concentration is plotted along the x-axis. The file'below.hull' " << std::endl;
  std::cout << "      contains a list of the structures in 'energy.clex' which are below the " << std::endl;
  std::cout << "      green line. These correspond to the structures that the cluster        " << std::endl;
  std::cout << "      expansion wrongly predicts to be below the hull.                       " << std::endl;
  std::cout << std::endl << std::endl;
}
void print_ecistats_full_man() {
  std::cout << "  eci_search -ecistats energy eci.in corr.in population_file" << std::endl;
  std::cout << "      This calculates eci for each eciset in population_file.  Then for each eci it prints:" << std::endl;
  std::cout << "        frac_on: the fraction of ecisets in which this eci is non-zero" << std::endl;
  std::cout << "        mean, min, max, rms (nonzero): eci statistics for ecisets in which this eci is non-zero." << std::endl;
  std::cout << "        mean, min, max, rms : eci statistics including all ecisets" << std::endl;
  std::cout << std::endl << std::endl;

}
void print_calc_all_full_man() {
  std::cout << "  eci_search -calc_all N energy eci.in corr.in" << std::endl;
  std::cout << "      This calculates eci for all possible ecisets with N eci." << std::endl;
  std::cout << std::endl << std::endl;

}
void print_calc_directmin_full_man() {
  std::cout << "  eci_search -calc_directmin Nrand Nmin Nmax energy eci.in corr.in [population_file]" << std::endl;
  std::cout << "      This searches for an eciset which minimizes the cv score by doing a " << std::endl;
  std::cout << "      direct minimization (steepest descent) search. If a population_file " << std::endl;
  std::cout << "      is given, the search is done on each eciset in the population, and  " << std::endl;
  std::cout << "      Nrand is ignored. If no population_file is given, the search is done" << std::endl;
  std::cout << "      on Nrand ecisets in which Nmin eci are turned on randomly.          " << std::endl;
  std::cout << "                                                                          " << std::endl;
  std::cout << "      At each step in the search, one-by-one every eci is flipped, the cv " << std::endl;
  std::cout << "      score is calculated, and the eci is flipped back. The change which  " << std::endl;
  std::cout << "      results in the lower cv score is kept, and the process is repeated  " << std::endl;
  std::cout << "      until no eci flip lowers the cv score. During the search, only      " << std::endl;
  std::cout << "      flips which keep the total number of eci in the range               " << std::endl;
  std::cout << "                            Nmin <= Neci <= Nmax                          " << std::endl;
  std::cout << "      are allowed.                                                        " << std::endl;
  std::cout << std::endl << std::endl;
}
void print_calc_dfsmin_full_man() {
  std::cout << "  eci_search -calc_dfsmin Nrand Nstop Nmin Nmax energy eci.in corr.in [population_file]" << std::endl;
  std::cout << "      This searches for an eciset which minimizes the cv score by doing     " << std::endl;
  std::cout << "      a depth-first search. If a population_file is given, the search is    " << std::endl;
  std::cout << "      done on each eciset in the population, and Nrand is ignored. If no    " << std::endl;
  std::cout << "      population_file is given, the search is done on Nrand ecisets in      " << std::endl;
  std::cout << "      which Nmin eci are turned on randomly.                                " << std::endl;
  std::cout << "                                                                            " << std::endl;
  std::cout << "      At the beginning of the search a counter is set to zero. At each      " << std::endl;
  std::cout << "      step in the search, one-by-one every eci is flipped, the cv score is  " << std::endl;
  std::cout << "      calculated, and the eci is flipped back. Any eciset which results in  " << std::endl;
  std::cout << "      a lower cv score is stored in a queue. The eciset in the queue with   " << std::endl;
  std::cout << "      the lowest cv score is chosen to begin the next step. If this eciset  " << std::endl;
  std::cout << "      results in a new best cv score, the counter is reset to zero. If this " << std::endl;
  std::cout << "      eciset does not have a new best cv score the counter is incremented.  " << std::endl;
  std::cout << "      When the counter reaches Nstop the process stops. During the search,  " << std::endl;
  std::cout << "      only flips which keep the total number of eci in the range            " << std::endl;
  std::cout << "                            Nmin <= Neci <= Nmax                            " << std::endl;
  std::cout << "      are allowed.                                                          " << std::endl;
  std::cout << std::endl << std::endl;


}
void print_calc_ga_full_man() {
  std::cout << "  eci_search -calc_ga Npop Nmin Nmax Ngen Nmut energy eci.in corr.in [population_file]" << std::endl;
  std::cout << "      This searches for an eciset which minimizes the cv score by using    " << std::endl;
  std::cout << "      a genetic algorithm. If a population_file is given, that population  " << std::endl;
  std::cout << "      is used. If no population_file is given, the starting population     " << std::endl;
  std::cout << "      consists of Npop ecisets in which Nmin eci are turned on randomly.   " << std::endl;
  std::cout << "                                                                           " << std::endl;
  std::cout << "      At each step in the search, a new generation of Npop children are    " << std::endl;
  std::cout << "      created. Two 'parent' ecisets are randomly chosen from the           " << std::endl;
  std::cout << "      population and are mated to create each 'child' eciset, and the      " << std::endl;
  std::cout << "      cv score is calculated. The Npop children or parents with the lowest " << std::endl;
  std::cout << "      cv scores become the parents for the next generation. The process    " << std::endl;
  std::cout << "      process continues for Ngen generations.                              " << std::endl;
  std::cout << "                                                                           " << std::endl;
  std::cout << "      Mating is done by comparing each eci of the two parents. Any eci that" << std::endl;
  std::cout << "      is turned on in both or turned off in both is kept in that same state" << std::endl;
  std::cout << "      in the child. Any eci which is on for one parent and off for the     " << std::endl;
  std::cout << "      other is randomly chosen in the child. Then mutations are performed  " << std::endl;
  std::cout << "      by flipping each eci in the child with probability equal to          " << std::endl;
  std::cout << "      (Nmut/eciset.size()). After mutations, the limits on the total number" << std::endl;
  std::cout << "      of eci on are enforced by randomly choosing eci to flip on until     " << std::endl;
  std::cout << "      Neci >= Nmin, and then randomly choosing eci to flip off until       " << std::endl;
  std::cout << "      Neci <= Nmax.                                                        " << std::endl;
  std::cout << "                                                                           " << std::endl;
  std::cout << "                                                                           " << std::endl;
  std::cout << "      *** -old depracted behavior, input Nchild instead of Ngen ***        " << std::endl;
  std::cout << "      At each step in the search, two 'parent' ecisets are randomly chosen " << std::endl;
  std::cout << "      from the population and are mated to create a 'child' eciset, and the" << std::endl;
  std::cout << "      cv score is calculated. If the child has a lower cv score than the   " << std::endl;
  std::cout << "      eciset in the population with the highest cv score, the child is     " << std::endl;
  std::cout << "      added to the population and the worst eciset is discarded. The       " << std::endl;
  std::cout << "      process continues for Nchild children.                               " << std::endl;
  std::cout << "                                                                           " << std::endl;
  std::cout << std::endl << std::endl;
}
void print_calc_ga_dir_full_man() {
  std::cout << "  eci_search -calc_ga_dir Npop Nmin Nmax Nchild Nmut energy eci.in corr.in [population_file]" << std::endl;
  std::cout << "      This searches for an eciset which minimizes the cv score by using     " << std::endl;
  std::cout << "      a genetic algorithm combined with a direct minimization search. If    " << std::endl;
  std::cout << "      a population_file is given, that population is used. If no            " << std::endl;
  std::cout << "      population_file is given, the starting population consists of Npop    " << std::endl;
  std::cout << "      ecisets in which Nmin eci are turned on randomly.                     " << std::endl;
  std::cout << "                                                                            " << std::endl;
  std::cout << "      At each step in the search, a new generation of Npop children are     " << std::endl;
  std::cout << "      created. Two 'parent' ecisets are randomly chosen from the            " << std::endl;
  std::cout << "      population and are mated to create each 'child' eciset. Each child is " << std::endl;
  std::cout << "      minimized using a direct minimization. The Npop children or parents   " << std::endl;
  std::cout << "      with the lowest cv scores become the parents for the next generation. " << std::endl;
  std::cout << "      The process process continues for Ngen generations.                   " << std::endl;
  std::cout << "                                                                            " << std::endl;
  std::cout << "      *** -old depracted behavior, input Nchild instead of Ngen ***         " << std::endl;
  std::cout << "      At each step in the search, two 'parent' ecisets are randomly chosen  " << std::endl;
  std::cout << "      from the population and are mated to create an initial 'child' eciset," << std::endl;
  std::cout << "      a direct minimization search is performed to optimize the child, and  " << std::endl;
  std::cout << "      the cv score is calculated. If the optimized child has a lower cv     " << std::endl;
  std::cout << "      score than the eciset in the population with the highest cv score, the" << std::endl;
  std::cout << "      optimized child is  added to the population and the worst eciset is   " << std::endl;
  std::cout << "      discarded. The process continues for Nchild children.                 " << std::endl;
  std::cout << "                                                                            " << std::endl;
  std::cout << "      Mating and the direct minimization search are done as in -calc_ga and " << std::endl;
  std::cout << "      -calc_directmin.                                                      " << std::endl;
  std::cout << std::endl << std::endl;

}
void print_calc_ga_dfs_full_man() {
  std::cout << "  eci_search -calc_ga_dfs Npop Nmin Nmax Nchild Nmut Nstop energy eci.in corr.in [population_file]" << std::endl;
  std::cout << "      This searches for an eciset which minimizes the cv score by using     " << std::endl;
  std::cout << "      a genetic algorithm combined with a depth-first search. If a          " << std::endl;
  std::cout << "      population_file is given, that population is used. If no              " << std::endl;
  std::cout << "      population_file is given, the starting population consists of Npop    " << std::endl;
  std::cout << "      ecisets in which Nmin eci are turned on randomly.                     " << std::endl;
  std::cout << "                                                                            " << std::endl;
  std::cout << "      At each step in the search, a new generation of Npop children are     " << std::endl;
  std::cout << "      created. Two 'parent' ecisets are randomly chosen from the            " << std::endl;
  std::cout << "      population and are mated to create each 'child' eciset. Each child is " << std::endl;
  std::cout << "      minimized using a depth first search. The Npop children or parents    " << std::endl;
  std::cout << "      with the lowest cv scores become the parents for the next generation. " << std::endl;
  std::cout << "      The process process continues for Ngen generations.                   " << std::endl;
  std::cout << "                                                                            " << std::endl;
  std::cout << "      *** -old depracted behavior, input Nchild instead of Ngen ***         " << std::endl;
  std::cout << "      At each step in the search, two 'parent' ecisets are randomly chosen  " << std::endl;
  std::cout << "      from the population and are mated to create an initial 'child' eciset," << std::endl;
  std::cout << "      a depth-first search is performed to optimize the child, and the cv   " << std::endl;
  std::cout << "      score is calculated. If the optimized child has a lower cv score than " << std::endl;
  std::cout << "      the eciset in the population with the highest cv score, the optimized " << std::endl;
  std::cout << "      child is added to the population and the worst eciset is discarded.   " << std::endl;
  std::cout << "      The process continues for Nchild children.                            " << std::endl;
  std::cout << "                                                                            " << std::endl;
  std::cout << "      Mating and the depth-first search are done as in -calc_ga and         " << std::endl;
  std::cout << "      -calc_dfsmin.                                                         " << std::endl;
  std::cout << std::endl << std::endl;
}
void print_weight_nrg_full_man() {
  std::cout << "  eci_search -weight_nrg A B kT energy weighted_energy" << std::endl;
  std::cout << "      This reads in the file 'energy', and writes the file 'weighted_energy'" << std::endl;
  std::cout << "      with the weights set to: weight = A*exp(-dist_from_hull/kT) + B." << std::endl;
  std::cout << std::endl << std::endl;
}
void print_weight_emin_full_man() {
  std::cout << "  eci_search -weight_emin A B kT energy weighted_energy" << std::endl;
  std::cout << "      This reads in the file 'energy', and writes the file 'weighted_energy'" << std::endl;
  std::cout << "      with the weights set to: weight = A*exp(-dist_from_minEnergy/kT) + B." << std::endl;
  std::cout << std::endl << std::endl;
}
void print_weight_eref_full_man() {
  std::cout << "  eci_search -weight_eref A B E0 kT energy weighted_energy" << std::endl
            << "      This reads in the file 'energy', and writes the file 'weighted_energy'" << std::endl
            << "      with the weights set to: weight = A*exp[-(E-E0)/kT] + B if E >  E0 and" << std::endl
            << "                               weight = 1.0                   if E <= E0" << std::endl
            << std::endl << std::endl;
}
void print_convert_full_man() {
  std::cout << "  eci_search -convert-to-text energy.json eci.in.json corr.in.json" << std::endl;
  std::cout << "  eci_search -convert-to-json energy eci.in corr.in" << std::endl;
  std::cout << "  eci_search -convert-energy-to-text energy.json [...]" << std::endl;
  std::cout << "  eci_search -convert-energy-to-json energy [...]" << std::endl;
  std::cout << "  eci_search -convert-eci-to-text eci.in.json [...]" << std::endl;
  std::cout << "  eci_search -convert-eci-to-json eci.in [...]" << std::endl;
  std::cout << "  eci_search -convert-corr-to-text corr.in.json [...]" << std::endl;
  std::cout << "  eci_search -convert-corr-to-json corr.in [...]" << std::endl;
  std::cout << "      Convert 'energy', 'eci.in', and 'corr.in' files to/from json          " << std::endl;

}
void print_calc_cs_fpc_full_man() {
  std::cout << "  eci_search -calc_cs_fpc energy eci.in corr.in mu" << std::endl;
  std::cout << "      This searches for an eciset using compressive sensing, with the       " << std::endl;
  std::cout << "      fixed-point continuation algorithm. Decreasing mu results in using    " << std::endl;
  std::cout << "      more ECI. See Nelson, Hart, Zhou, and Ozolins, PRB, 87, 035125 (2013)." << std::endl;
  std::cout << std::endl << std::endl;
}
void print_calc_cs_bi_full_man() {
  std::cout << "  eci_search -calc_cs_fpc energy eci.in corr.in mu" << std::endl;
  std::cout << "      This searches for an eciset using compressive sensing, with the       " << std::endl;
  std::cout << "      Bregman iteration algorithm. Decreasing mu results in using more ECI. " << std::endl;
  std::cout << "      See Nelson, Hart, Zhou, and Ozolins, PRB, 87, 035125 (2013).          " << std::endl;
  std::cout << std::endl << std::endl;
}

void print_eci_search_quick_man() {
  std::cout << "*** eci_search quick manual ***" << std::endl;
  std::cout << "  eci_search -help" << std::endl;
  print_calc_man();
  print_ecistats_man();
  print_calc_all_man();
  print_calc_directmin_man();
  print_calc_dfsmin_man();
  print_calc_ga_man();
  print_calc_ga_dir_man();
  print_calc_ga_dfs_man();
  print_weight_nrg_man();
  print_weight_emin_man();
  print_weight_eref_man();
  print_convert_man();

  std::cout << std::endl << "  * under development *" << std::endl;
  print_calc_cs_fpc_man();
  print_calc_cs_bi_man();

  std::cout << "\n\n  Note: Use 'FixOn' or 'FixOff' for the 'weight' in the 'eci.in' file to set particular eci on/off manually." << std::endl << std::endl;

  std::cout << "  Note: Use '-pthreads X' to specify number of threads used for parallel processing ECISets in a population. Default is number of cores." << std::endl << std::endl;

  std::cout << "  Note: Use '-mthreads X' to specify number of threads to use for each minimization step. Default is generally 1, but more if population size is less than pthreads." << std::endl << std::endl;

  std::cout << "  Note: Use '-tol X' to set hull finding tolerances. Default is 1.0e-14." << std::endl << std::endl;

  std::cout << "  Note: Use '-old' to use deprecated serial functions." << std::endl << std::endl;

  std::cout << "  Note: Use '-maxstep X' to specify maximum number of FPC iterations during calc_cs_fpc and calc_cs_bi" << std::endl << std::endl;


};

void print_eci_search_man() {

  print_eci_search_quick_man();

  std::cout << std::endl << std::endl << "*** eci_search full manual ***" << std::endl << std::endl;

  std::cout << "  Note: Use 'FixOn' or 'FixOff' for the 'weight' in the 'eci.in' file to set particular eci on/off manually." << std::endl << std::endl;

  std::cout << "  Note: Use '-pthreads X' to specify number of threads used for parallel processing ECISets in a population. Default is number of cores." << std::endl << std::endl;

  std::cout << "  Note: Use '-mthreads X' to specify number of threads to use for each minimization step. Default is generally 1, but more if population size is less than pthreads." << std::endl << std::endl;

  std::cout << "  Note: Use '-tol X' to set hull finding tolerances. Default is 1.0e-14." << std::endl << std::endl;

  std::cout << "  Note: Use '-old' to use deprecated serial functions." << std::endl << std::endl;

  std::cout << "  Note: Use '-maxstep X' to specify maximum number of FPC iterations during calc_cs_fpc and calc_cs_bi" << std::endl << std::endl;

  print_calc_full_man();
  print_ecistats_full_man();
  print_calc_all_full_man();
  print_calc_directmin_full_man();
  print_calc_dfsmin_full_man();
  print_calc_ga_full_man();
  print_calc_ga_dir_full_man();
  print_calc_ga_dfs_full_man();
  print_weight_nrg_full_man();
  print_weight_emin_full_man();
  print_weight_eref_full_man();
  print_convert_full_man();

  std::cout << "  * under development *" << std::endl;
  print_calc_cs_fpc_full_man();
  std::cout << "  * under development *" << std::endl;
  print_calc_cs_bi_full_man();

};

int main(int argc, char *argv[]) {
  //std::cout << "BEGIN eci_search" << std::endl;


  //BP::BP_pause();
  int argc_adjustment = 0;
  unsigned long int max_step;

  NEW = true;

  BP::BP_Vec<string> args;
  for(int i = 0; i < argc; i++) {
    std::cout << argv[i] << " ";
    if(std::string(argv[i]) == "-test") {
      TEST = true;
      argc_adjustment++;
    }
    else if(std::string(argv[i]) == "-mthreads") {
      // threading for minimization step (toggle parallelization)
      i++;
      std::cout << argv[i] << " ";
      MTHREADS = BP::stoi(string(argv[i]));
      argc_adjustment += 2;

    }
    else if(std::string(argv[i]) == "-pthreads") {
      // threading for population (parellelize minimization of each ECISet in population)
      i++;
      std::cout << argv[i] << " ";
      PTHREADS = BP::stoi(string(argv[i]));
      argc_adjustment += 2;
    }
    else if(std::string(argv[i]) == "-tol") {
      // threading for population (parellelize minimization of each ECISet in population)
      i++;
      std::cout << argv[i] << " ";
      HULLTOL = BP::stod(string(argv[i]));
      argc_adjustment += 2;
    }
    else if(std::string(argv[i]) == "-old") {
      // use original functions
      NEW = false;
      argc_adjustment += 1;
    }
    else if(std::string(argv[i]) == "-maxstep") {
      // use original functions
      i++;
      std::cout << argv[i] << " ";
      max_step = static_cast<unsigned long int>(BP::stod(string(argv[i])));
      argc_adjustment += 2;
    }
    else {
      args.add(string(argv[i]));
    }
  }
  std::cout << std::endl;
  argc -= argc_adjustment;

  //std::cout << "argc: " << argc << std::endl;
  //std::cout << "args: " << args << std::endl;

  BP::BP_StopWatch timer;
  timer.set_start();

  if(TEST) {
    BP::BP_pause();
    std::cout << "WARNING!!! THIS IS A TEST!!! THE RANDOM NUMBERS MAY NOT BE RANDOM!!!" << std::endl;
  }



  try {
    if(argc == 1) {
      //std::cout << "print_eci_search_man()" << std::endl;
      print_eci_search_quick_man();
      return 1;
    }
    else if(argc >= 2) {
      if(args[1] == "-calc") {
        if(!NEW) {
          //eci_search -calc energy eci.in corr.in [bitstring | bitstring_file]
          // writes eci.in with bitstring
          // writes eci.out
          if(argc == 5) {
            BP::BP_Vec<ECISet> population;

            calc_eci(args[2], args[3], args[4], population, HULLTOL);
          }
          else if(argc == 6) {
            BP::BP_Vec<ECISet> population;
            if(is_bitstring(args[5])) {
              ECISet eci_in(args[3]);
              eci_in.set_bit_string(args[5]);
              population.add(eci_in);
            }
            else {
              population = read_bit_strings_file(ECISet(args[3]), args[5]);
              if(population.size() != 1) {
                std::cout << "error, the file " << args[5] << " has multiple bit strings." << std::endl;
                exit(1);
              }
            }

            calc_eci(args[2], args[3], args[4], population);
          }
          else {
            //std::cout << "print_eci_search_man()" << std::endl;
            print_calc_man();
            return 1;
          }
        }
        else {
          if(argc == 5 || argc == 6) {
            Population population(0, 0, args[2], args[3], args[4]);

            if(argc == 5) {
              population.populate_with_base();
            }
            else {
              if(is_bitstring(args[5]))
                population.add(args[5]);
              else
                population.populate(args[5]);
            }

            if(population.size() == 1)
              population.calc_details(0, "default", HULLTOL);
            else
              population.calc(PTHREADS);
          }
          else {
            //std::cout << "print_eci_search_man()" << std::endl;
            print_calc_man();

            return 1;
          }
        }

      }
      else if(args[1] == "-convert-to-text") {
        //eci_search -convert-to-text energy.json eci.in.json corr.in.json
        if(argc == 5) {

          EnergySet nrgset(args[2]);
          nrgset.write(args[2], "text");

          ECISet eciset(args[3]);
          eciset.write_ECIin(args[3], "text");

          Correlation corr(args[4]);
          corr.write(args[4], "text");

        }
        else {
          print_convert_man();
        }

      }
      else if(args[1] == "-convert-to-json") {
        //eci_search -convert-to-json energy eci.in corr.in
        if(argc == 5) {

          EnergySet nrgset(args[2]);
          nrgset.write(args[2], "json");

          ECISet eciset(args[3]);
          eciset.write_ECIin(args[3], "json");

          Correlation corr(args[4]);
          corr.write(args[4], "json");

        }
        else {
          print_convert_man();
        }

      }
      else if(args[1] == "-convert-energy-to-text") {
        //eci_search -convert-energy-to-text energy [...]
        if(argc > 2) {

          for(int i = 2; i < args.size(); i++) {
            EnergySet nrgset(args[i]);
            nrgset.write(args[i], "text");
          }
        }
        else {
          print_convert_man();
        }

      }
      else if(args[1] == "-convert-energy-to-json") {
        //eci_search -convert-energy-to-json energy [...]
        if(argc > 2) {

          for(int i = 2; i < args.size(); i++) {
            EnergySet nrgset(args[i]);
            nrgset.write(args[i], "json");
          }
        }
        else {
          print_convert_man();
        }

      }
      else if(args[1] == "-convert-eci-to-text") {
        //eci_search -convert-eci-to-text energy [...]
        if(argc > 2) {

          for(int i = 2; i < args.size(); i++) {
            ECISet eciset(args[i]);
            eciset.write_ECIin(args[i], "text");
          }
        }
        else {
          print_convert_man();
        }

      }
      else if(args[1] == "-convert-eci-to-json") {
        //eci_search -convert-eci-to-json energy [...]
        if(argc > 2) {

          for(int i = 2; i < args.size(); i++) {
            ECISet eciset(args[i]);
            eciset.write_ECIin(args[i], "json");
          }
        }
        else {
          print_convert_man();
        }

      }
      else if(args[1] == "-convert-corr-to-text") {
        //eci_search -convert-corr-to-text energy [...]
        if(argc > 2) {

          for(int i = 2; i < args.size(); i++) {
            Correlation corr(args[i]);
            corr.write(args[i], "text");
          }
        }
        else {
          print_convert_man();
        }

      }
      else if(args[1] == "-convert-corr-to-json") {
        //eci_search -convert-corr-to-json energy [...]
        if(argc > 2) {

          for(int i = 2; i < args.size(); i++) {
            Correlation corr(args[i]);
            corr.write(args[i], "json");
          }
        }
        else {
          print_convert_man();
        }

      }
      else if(args[1] == "-calc_cs_fpc") {
        //eci_search -calc_cs_fpc energy eci.in corr.in mu
        if(argc == 6) {
          BP::BP_Vec<double> mu;
          mu.add(BP::stod(args[5]));

          calc_cs_eci(args[2], args[3], args[4], mu, 0, HULLTOL, max_step);
        }
        else {
          //std::cout << "print_eci_search_man()" << std::endl;
          //print_calc_man();
        }

      }
      else if(args[1] == "-calc_cs_bi") {
        //eci_search -calc_cs_bi energy eci.in corr.in mu
        if(argc == 6) {
          BP::BP_Vec<double> mu;
          mu.add(BP::stod(args[5]));

          calc_cs_eci(args[2], args[3], args[4], mu, 1, HULLTOL, max_step);
        }
        else {
          //std::cout << "print_eci_search_man()" << std::endl;
          //print_calc_man();
        }

      }
      else if(args[1] == "-ecistats") {
        //eci_search -ecistats energy eci.in corr.in population_file

        if(!NEW) {
          // writes eci.in with bitstring
          // writes eci.out
          if(argc == 6) {
            BP::BP_Vec<ECISet> population = read_bit_strings_file(ECISet(args[3]), args[5]);

            calc_ecistats(args[2], args[3], args[4], population);
          }
          else {
            //std::cout << "print_ecistats_man()" << std::endl;
            print_ecistats_man();

            return 1;
          }
        }
        else {
          if(argc == 6) {
            Population population(0, 0, args[2], args[3], args[4]);
            population.populate(args[5]);
            population.calc_ecistats(PTHREADS);
          }
          else {
            //std::cout << "print_ecistats_man()" << std::endl;
            print_ecistats_man();

            return 1;
          }
        }

      }
      else if(args[1] == "-calc_all") {
        if(!NEW) {
          if(argc == 6) {
            //eci_search -calc_all N energy eci.in corr.in
            calc_all_eci(BP::stoi(args[2]), args[3], args[4], args[5]);
          }
          else {
            //std::cout << "print_eci_search_man()" << std::endl;
            print_calc_all_man();

            return 1;
          }
        }
        else {
          if(argc == 6) {
            //eci_search -calc_all N energy eci.in corr.in
            Population population(0, 0, args[3], args[4], args[5]);
            population.calc_all(BP::stoi(args[2]), PTHREADS);
          }
          else {
            //std::cout << "print_eci_search_man()" << std::endl;
            print_calc_all_man();

            return 1;
          }
        }
      }
      else if(args[1] == "-calc_directmin") {

        //eci_search -calc_directmin Npop Nmin Nmax energy eci.in corr.in [population_file]
        if(!NEW) {
          if(argc == 8) {
            int Nrand = BP::stoi(args[2]);
            int Nmin = BP::stoi(args[3]);
            int Nmax = BP::stoi(args[4]);
            BP::BP_Vec<ECISet> population;

            if(Nmin > 0 && Nmax > Nmin)
              calc_directmin_eci(Nrand, Nmin, Nmax, args[5], args[6], args[7], population);
            else {
              //std::cout << "print_eci_search_man()" << std::endl;
              print_calc_directmin_man();

              return 1;
            }
          }
          else if(argc == 9) {
            int Nrand = BP::stoi(args[2]);
            int Nmin = BP::stoi(args[3]);
            int Nmax = BP::stoi(args[4]);
            BP::BP_Vec<ECISet> population = read_bit_strings_file(ECISet(args[6]), args[8]);

            if(Nmin > 0 && Nmax > Nmin)
              calc_directmin_eci(Nrand, Nmin, Nmax, args[5], args[6], args[7], population);
            else {
              //std::cout << "print_eci_search_man()" << std::endl;
              print_calc_directmin_man();

              return 1;
            }
          }
          else {
            //std::cout << "print_eci_search_man()" << std::endl;
            print_calc_directmin_man();

            return 1;
          }
        }
        else {
          if(argc == 8 || argc == 9) {
            int Npop = BP::stoi(args[2]);
            int Nmin = BP::stoi(args[3]);
            int Nmax = BP::stoi(args[4]);

            Population population(Nmin, Nmax, args[5], args[6], args[7]);

            if(argc == 8)
              population.populate(Npop);
            else
              population.populate(args[8]);

            population.direct(PTHREADS, MTHREADS);

            population.write(unique("population"), "default");
          }
          else {
            //std::cout << "print_eci_search_man()" << std::endl;
            print_calc_directmin_man();

            return 1;
          }
        }

      }
      else if(args[1] == "-calc_dfsmin") {
        if(!NEW) {
          if(argc == 9) {
            //eci_search -calc_dfsmin Npop Nstop Nmin Nmax energy eci.in corr.in [population_file]
            int Nrand = BP::stoi(args[2]);
            int Nstop = BP::stoi(args[3]);
            int Nmin = BP::stoi(args[4]);
            int Nmax = BP::stoi(args[5]);
            BP::BP_Vec<ECISet> population;

            if(Nmin > 0 && Nmax > Nmin)
              calc_dfsmin_eci(Nrand, Nstop, Nmin, Nmax, args[6], args[7], args[8], population);
            else {
              //std::cout << "print_eci_search_man()" << std::endl;
              print_calc_dfsmin_man();

              return 1;
            }
          }
          else if(argc == 10) {
            //eci_search -calc_dfsmin Npop Nstop Nmin Nmax energy eci.in corr.in [population_file]
            int Nrand = BP::stoi(args[2]);
            int Nstop = BP::stoi(args[3]);
            int Nmin = BP::stoi(args[4]);
            int Nmax = BP::stoi(args[5]);
            BP::BP_Vec<ECISet> population = read_bit_strings_file(ECISet(args[7]), args[9]);

            if(Nmin > 0 && Nmax > Nmin)
              calc_dfsmin_eci(Nrand, Nstop, Nmin, Nmax, args[6], args[7], args[8], population);
            else {
              //std::cout << "print_eci_search_man()" << std::endl;
              print_calc_dfsmin_man();

              return 1;
            }
          }
          else {
            //std::cout << "print_eci_search_man()" << std::endl;
            print_calc_dfsmin_man();

            return 1;
          }
        }
        else {

          if(argc == 9 || argc == 10) {
            //eci_search -calc_dfsmin Npop Nstop Nmin Nmax energy eci.in corr.in [population_file]

            int Npop = BP::stoi(args[2]);
            int Nstop = BP::stoi(args[3]);
            int Nmin = BP::stoi(args[4]);
            int Nmax = BP::stoi(args[5]);

            Population population(Nmin, Nmax, args[6], args[7], args[8]);

            if(argc == 9)
              population.populate(Npop);
            else
              population.populate(args[9]);

            population.dfs(Nstop, PTHREADS, MTHREADS);

            population.write(unique("population"), "default");

          }
          else {
            print_calc_dfsmin_man();

            return 1;
          }
        }

      }
      else if(args[1] == "-calc_ga") {
        //eci_search -calc_ga Npop Nmin Nmax Nchild Nmut energy eci.in corr.in [population_file]
        if(!NEW) {

          if(argc == 10) {
            int Npopulation = BP::stoi(args[2]);
            int Nmin = BP::stoi(args[3]);
            int Nmax = BP::stoi(args[4]);
            int Nchildren = BP::stoi(args[5]);
            int Nmutations = BP::stoi(args[6]);

            BP::BP_Vec<ECISet> population;

            if(Nmin > 0 && Nmax > Nmin)
              calc_ga_eci(Npopulation, Nmin, Nmax, Nchildren, Nmutations, args[7], args[8], args[9], population);
            else {
              //std::cout << "print_eci_search_man()" << std::endl;
              print_calc_ga_man();

              return 1;
            }
          }
          else if(argc == 11) {
            int Npopulation = BP::stoi(args[2]);
            int Nmin = BP::stoi(args[3]);
            int Nmax = BP::stoi(args[4]);
            int Nchildren = BP::stoi(args[5]);
            int Nmutations = BP::stoi(args[6]);
            BP::BP_Vec<ECISet> population = read_bit_strings_file(ECISet(args[8]), args[10]);

            if(Nmin > 0 && Nmax > Nmin)
              calc_ga_eci(Npopulation, Nmin, Nmax, Nchildren, Nmutations, args[7], args[8], args[9], population);
            else {
              //std::cout << "print_eci_search_man()" << std::endl;
              print_calc_ga_man();

              return 1;
            }
          }
          else {
            //std::cout << "print_eci_search_man()" << std::endl;
            print_calc_ga_man();

            return 1;
          }
        }
        else {
          //eci_search -calc_ga Npop Nmin Nmax Ngen Nmut energy eci.in corr.in [population_file]
          if(argc == 10 || argc == 11) {
            int Npopulation = BP::stoi(args[2]);
            int Nmin = BP::stoi(args[3]);
            int Nmax = BP::stoi(args[4]);
            int Ngenerations = BP::stoi(args[5]);
            int Nmutations = BP::stoi(args[6]);
            Population population(Nmin, Nmax, args[7], args[8], args[9]);

            if(argc == 10)
              population.populate(Npopulation);
            else
              population.populate(args[10]);

            population.ga(Ngenerations, Nmutations, PTHREADS);

            population.write(unique("population"), "default");
          }
          else {
            //std::cout << "print_eci_search_man()" << std::endl;
            print_calc_ga_man();

            return 1;
          }
        }
      }
      else if(args[1] == "-calc_ga_dir") {
        //eci_search -calc_ga_dir Npop Nmin Nmax Nchild Nmut energy eci.in corr.in [population_file]
        if(!NEW) {
          if(argc == 10) {
            int Npopulation = BP::stoi(args[2]);
            int Nmin = BP::stoi(args[3]);
            int Nmax = BP::stoi(args[4]);
            int Nchildren = BP::stoi(args[5]);
            int Nmutations = BP::stoi(args[6]);

            BP::BP_Vec<ECISet> population;

            if(Nmin > 0 && Nmax > Nmin)
              calc_ga_dir_eci(Npopulation, Nmin, Nmax, Nchildren, Nmutations, args[7], args[8], args[9], population);
            else {
              //std::cout << "print_eci_search_man()" << std::endl;
              print_calc_ga_dir_man();

              return 1;
            }
          }
          else if(argc == 11) {
            int Npopulation = BP::stoi(args[2]);
            int Nmin = BP::stoi(args[3]);
            int Nmax = BP::stoi(args[4]);
            int Nchildren = BP::stoi(args[5]);
            int Nmutations = BP::stoi(args[6]);
            BP::BP_Vec<ECISet> population = read_bit_strings_file(ECISet(args[8]), args[10]);

            if(Nmin > 0 && Nmax > Nmin)
              calc_ga_dir_eci(Npopulation, Nmin, Nmax, Nchildren, Nmutations, args[7], args[8], args[9], population);
            else {
              //std::cout << "print_eci_search_man()" << std::endl;
              print_calc_ga_dir_man();

              return 1;
            }
          }
          else {
            //std::cout << "print_eci_search_man()" << std::endl;
            print_calc_ga_dir_man();

            return 1;
          }
        }
        else {
          if(argc == 10 || argc == 11) {
            int Npopulation = BP::stoi(args[2]);
            int Nmin = BP::stoi(args[3]);
            int Nmax = BP::stoi(args[4]);
            int Ngenerations = BP::stoi(args[5]);
            int Nmutations = BP::stoi(args[6]);
            Population population(Nmin, Nmax, args[7], args[8], args[9]);

            if(argc == 10)
              population.populate(Npopulation);
            else
              population.populate(args[10]);

            population.ga_dir(Ngenerations, Nmutations, PTHREADS);

            population.write(unique("population"), "default");
          }
          else {
            //std::cout << "print_eci_search_man()" << std::endl;
            print_calc_ga_dir_man();

            return 1;
          }
        }
      }
      else if(args[1] == "-calc_ga_dfs") {
        //eci_search -calc_ga_dfs Npop Nmin Nmax Nchild Nmut Nstop energy eci.in corr.in [population_file]
        if(!NEW) {
          if(argc == 11) {
            int Npopulation = BP::stoi(args[2]);
            int Nmin = BP::stoi(args[3]);
            int Nmax = BP::stoi(args[4]);
            int Nchildren = BP::stoi(args[5]);
            int Nmutations = BP::stoi(args[6]);
            int Nstop = BP::stoi(args[7]);

            BP::BP_Vec<ECISet> population;

            if(Nmin > 0 && Nmax > Nmin)
              calc_ga_dfs_eci(Npopulation, Nmin, Nmax, Nchildren, Nmutations, Nstop, args[8], args[9], args[10], population);
            else {
              //std::cout << "print_eci_search_man()" << std::endl;
              print_calc_ga_dfs_man();

              return 1;
            }
          }
          else if(argc == 12) {
            int Npopulation = BP::stoi(args[2]);
            int Nmin = BP::stoi(args[3]);
            int Nmax = BP::stoi(args[4]);
            int Nchildren = BP::stoi(args[5]);
            int Nmutations = BP::stoi(args[6]);
            int Nstop = BP::stoi(args[7]);

            BP::BP_Vec<ECISet> population = read_bit_strings_file(ECISet(args[9]), args[11]);

            if(Nmin > 0 && Nmax > Nmin)
              calc_ga_dfs_eci(Npopulation, Nmin, Nmax, Nchildren, Nmutations, Nstop, args[8], args[9], args[10], population);
            else {
              //std::cout << "print_eci_search_man()" << std::endl;
              print_calc_ga_dfs_man();

              return 1;
            }
          }
          else {
            //std::cout << "print_eci_search_man()" << std::endl;
            print_calc_ga_dfs_man();

            return 1;
          }
        }
        else {
          if(argc == 11 || argc == 12) {
            int Npopulation = BP::stoi(args[2]);
            int Nmin = BP::stoi(args[3]);
            int Nmax = BP::stoi(args[4]);
            int Ngenerations = BP::stoi(args[5]);
            int Nmutations = BP::stoi(args[6]);
            int Nstop = BP::stoi(args[7]);
            Population population(Nmin, Nmax, args[8], args[9], args[10]);

            if(argc == 11)
              population.populate(Npopulation);
            else
              population.populate(args[11]);

            population.ga_dfs(Nstop, Ngenerations, Nmutations, PTHREADS);

            population.write(unique("population"), "default");
          }
          else {
            //std::cout << "print_eci_search_man()" << std::endl;
            print_calc_ga_dfs_man();

            return 1;
          }
        }
      }
      else if(args[1] == "-weight_nrg") {
        //eci_search -weight_nrg A B kT energy weighted_energy
        if(argc == 7) {
          weight_nrg(BP::stod(args[2]), BP::stod(args[3]), BP::stod(args[4]), args[5], args[6]);

        }
        else {
          //std::cout << "print_eci_search_man()" << std::endl;
          print_weight_nrg_man();

          return 1;

        }
      }
      else if(args[1] == "-nrg_diff") {
        //eci_search -nrg_diff tot_corr.index.A AB_structs.energy eci.in.tot AB_structs.corr.in.tot eciset.A
        //eci_search -nrg_diff   tot_corr.index.A    AB_structs.energy   eci.in.tot   AB_structs.corr.in.tot   A_structs.energy  A_structs.eci.in  A_structs.corr.in
        if(argc == 9) {
          calc_nrg_diff(args[2], args[3], args[4], args[5], args[6],  args[7],  args[8]);

        }
        else {
          //std::cout << "print_eci_search_man()" << std::endl;
          std::cout << "eci_search -nrg_diff    tot_corr.index.A    AB_structs.energy   eci.in.tot   AB_structs.corr.in.tot   A_structs.energy  A_structs.eci.in  A_structs.corr.in" << std::endl;
        }
      }
      else if(args[1] == "-weight_emin") {
        //eci_search -weight_emin A B kT energy weighted_energy
        if(argc == 7) {
          weight_EMin(BP::stod(args[2]), BP::stod(args[3]), BP::stod(args[4]), args[5], args[6]);
        }
        else {
          //std::cout << "print_eci_search_man()" << std::endl;
          print_weight_emin_man();

          return 1;
        }
      }
      else if(args[1] == "-weight_eref") {
        //eci_search -weight_eref A B E0 kT energy weighted_energy
        if(argc == 8) {
          weight_ERef(BP::stod(args[2]), BP::stod(args[3]), BP::stod(args[4]), BP::stod(args[5]), args[6], args[7]);
        }
        else {
          //std::cout << "print_eci_search_man()" << std::endl;
          print_weight_eref_man();

          return 1;
        }
      }
      else {
        //std::cout << "print_eci_search_man()" << std::endl;
        print_eci_search_man();

        return 1;
      }
    }
    else {
      //std::cout << "print_eci_search_man()" << std::endl;
      print_eci_search_man();

      return 1;
    }
  }
  catch(std::exception &e) {
    std::cerr << "Error in eci_search: " << e.what() << std::endl;
    return 1;
  }

  std::cout << "TIME: " << timer.total_time_s() << std::endl;

  if(TEST) {
    std::cout << "WARNING!!! THIS IS A TEST!!! THE RANDOM NUMBERS MAY NOT BE RANDOM!!!" << std::endl;
  }

  return 0;
}
