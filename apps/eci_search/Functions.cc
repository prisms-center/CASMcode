/*
 *  Functions.cc
 */

#ifndef Functions_CC
#define Functions_CC

#include <cstdlib>
#include "jsonParser.hh"
#include "BP_Dir.hh"
#include "Functions.hh"
#include "Correlation.hh"
#include "ECISet.hh"
#include "EnergySet.hh"
#include "GeneticAlgorithm.hh"

////----------------------------
/// main functions
void calc_eci(std::string energy_filename, std::string eci_in_filename, std::string corr_in_filename, BP::BP_Vec<ECISet> &population, double hulltol) {
  Correlation corr(corr_in_filename);
  EnergySet DFT_nrg(energy_filename);
  ECISet eci_in(eci_in_filename);


  bool singular;

  if(population.size() == 1)
    eci_in = population[0];

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
  corr.write_covar(eci_in, DFT_nrg, "corr.covar", "default");

  // write eci.in
  eci_in.write_ECIin("eci.in", "default");
  std::cout << "Wrote 'eci.in'" << std::endl;

  // write eci.out
  eci_in.write_ECIout("eci.out", DFT_nrg, "default");
  std::cout << "Wrote 'eci.out'" << std::endl;

  //std::cout << "Hull of: " << energy_filename << std::endl;
  //DFT_nrg.write_hull(std::cout);
  std::cout << std::endl << "Calculating hull of: " << energy_filename << std::endl;
  DFT_nrg.calc_hull(true, hulltol);
  DFT_nrg.write_hull("hull", "default");
  std::cout << "Wrote 'hull'" << std::endl;


  if(eci_in.get_cv() != 1e20) {

    // calc energy_clex
    EnergySet energy_clex = DFT_nrg;
    energy_clex.calc_clex(corr, eci_in);

    // write energy_clex
    std::cout << std::endl << "Calculating hull of: energy.clex " << std::endl;
    energy_clex.calc_hull(true, hulltol);
    //energy_clex.write_hull(std::cout);
    energy_clex.write_hull("hull.clex", "default");
    std::cout << "Wrote 'hull.clex'" << std::endl;
    energy_clex.write("energy.clex", "default");
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
      energy_clex_of_DFT_hull.write_hull("hull.clex_of_DFT_hull", "default");
      std::cout << "Wrote 'hull.clex_of_DFT_hull'" << std::endl;
      energy_clex_of_DFT_hull.write_below_hull("below.hull", energy_clex.get(DFT_nrg.get_non_hull_indices()), "default");
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
          //energy_clex_of_DFT_hull.calc_hull(true, hulltol);
          //energy_clex_of_DFT_hull.write_hull("hull.clex_of_DFT_hull");
          //std::cout << "  wrote 'hull.clex_of_DFT_hull'" << std::endl;
          //energy_clex_of_DFT_hull.write_below_hull( "below.hull", energy_clex.get(DFT_nrg.get_non_hull_indices()));
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

}

void calc_cs_eci(std::string energy_filename, std::string eci_in_filename, std::string corr_in_filename, const BP::BP_Vec<double> &mu, int alg, double hulltol, unsigned long int max_step) {
  // solve E = Corr*ECI, for ECI using compressive sensing L1 norm minimization (plus L2 norm in practice)

  Correlation corr(corr_in_filename);
  EnergySet DFT_nrg(energy_filename);
  ECISet eci_in(eci_in_filename);
  ECISet eci_CS = eci_in;

  unsigned long int i;
  double prec_shrink = 1e-6;
  double prec_ECI = 1e-6;
  bool print_steps = true;

  for(i = 0; i < mu.size(); i++) {
    if(alg == 0) {
      // use Fixed-point continuation method, calc for single value of mu
      eci_CS = calc_FPC_eci(mu[i], prec_shrink, prec_ECI, eci_in, corr, DFT_nrg, print_steps, max_step);
    }
    else if(alg == 1) {
      // use Bregman iteration method, calc for single value of mu
      eci_CS = calc_BI_eci(mu[i], prec_shrink, prec_ECI, eci_in, corr, DFT_nrg, print_steps, max_step);
    }
    else {
      std::cout << "Error in calc_cs_eci().  alg == '" << alg << "' is not a valid option." << std::endl;
      std::cout << "  Options are:  0, Fixed-point continuation" << std::endl;
      std::cout << "                1, Bergman iteration" << std::endl;
      exit(1);
    }

    std::cout << std::endl << eci_CS.get_bit_string() << "    mu: " << mu[i] << "  Nclust: " << eci_CS.get_Nclust_on() << " rms: " << eci_CS.get_rms() << std::endl;
  }

  // write eci.in
  eci_CS.write_ECIin("eci.in", "default");
  std::cout << "Wrote 'eci.in'" << std::endl;

  // write eci.out
  eci_CS.write_ECIout("eci.out", DFT_nrg, "default");
  std::cout << "Wrote 'eci.out'" << std::endl;

  //std::cout << "Hull of: " << energy_filename << std::endl;
  //DFT_nrg.write_hull(std::cout);
  std::cout << std::endl << "Calculating hull of: " << energy_filename << std::endl;
  DFT_nrg.calc_hull(true, hulltol);
  DFT_nrg.write_hull("hull", "default");
  std::cout << "Wrote 'hull'" << std::endl;


  //if( eci_CS.get_cv() != 1e20)
  {

    // calc energy_clex
    EnergySet energy_clex = DFT_nrg;
    energy_clex.calc_clex(corr, eci_CS);

    // write energy_clex
    std::cout << std::endl << "Calculating hull of: energy.clex " << std::endl;
    energy_clex.calc_hull(true, hulltol);
    //energy_clex.write_hull(std::cout);
    energy_clex.write_hull("hull.clex", "default");
    std::cout << "Wrote 'hull.clex'" << std::endl;
    energy_clex.write("energy.clex", "default");
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
      energy_clex_of_DFT_hull.write_hull("hull.clex_of_DFT_hull", "default");
      std::cout << "Wrote 'hull.clex_of_DFT_hull'" << std::endl;
      energy_clex_of_DFT_hull.write_below_hull("below.hull", energy_clex.get(DFT_nrg.get_non_hull_indices()), "default");
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
          //energy_clex_of_DFT_hull.calc_hull(true, hulltol);
          //energy_clex_of_DFT_hull.write_hull("hull.clex_of_DFT_hull");
          //std::cout << "  wrote 'hull.clex_of_DFT_hull'" << std::endl;
          //energy_clex_of_DFT_hull.write_below_hull( "below.hull", energy_clex.get(DFT_nrg.get_non_hull_indices()));
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

}

void calc_ecistats(std::string energy_filename, std::string eci_in_filename, std::string corr_in_filename, BP::BP_Vec<ECISet> &population) {
  Correlation corr(corr_in_filename);
  EnergySet DFT_nrg(energy_filename);
  //ECISet eci_in(eci_in_filename);
  BP::BP_Vec<double> eci_value_list;
  BP::BP_Vec<double> eci_nonzero_value_list;

  bool singular;

  for(int i = 0; i < population.size(); i++) {
    population[i].fit(corr, DFT_nrg, singular);
  }

  std::cout << "label" << "     " << "frac_on" << "     " << "mean(nonzero)" << "     "  << "min(nonzero)" << "     " << "max(nonzero)" << "     "  << "rms(nonzero)" << "     " << "mean" << "     "  << "min" << "     " << "max" << "     "  << "rms" << std::endl;
  for(int i = 0; i < population[0].size(); i++) {
    eci_value_list.clear();
    eci_nonzero_value_list.clear();
    for(int j = 0; j < population.size(); j++) {
      eci_value_list.add(population[j].get_value(i));
      if(population[j].get_weight(i)) {
        eci_nonzero_value_list.add(population[j].get_value(i));
      }
    }

    std::cout << i << " \t";
    std::cout << setprecision(6) << (1.0 * eci_nonzero_value_list.size()) / (1.0 * eci_value_list.size()) << " \t";
    if(eci_nonzero_value_list.size() != 0) {
      std::cout << mean(eci_nonzero_value_list) << " \t";
      std::cout << min(eci_nonzero_value_list) << " \t";
      std::cout << max(eci_nonzero_value_list) << " \t";
      std::cout << rms(eci_nonzero_value_list) << " \t";
      std::cout << mean(eci_value_list) << " \t";
      std::cout << min(eci_value_list) << " \t";
      std::cout << max(eci_value_list) << " \t";
      std::cout << rms(eci_value_list) << " \n";
    }
    else {
      std::cout << "-" << " \t";
      std::cout << "-" << " \t";
      std::cout << "-" << " \t";
      std::cout << "-" << " \t";
      std::cout << "-" << " \t";
      std::cout << "-" << " \t";
      std::cout << "-" << " \t";
      std::cout << "-" << " \n";
    }

  }

}

void calc_all_eci(int N, std::string energy_filename, std::string eci_in_filename, std::string corr_in_filename) {

  Correlation corr(corr_in_filename);
  EnergySet DFT_nrg(energy_filename);
  ECISet eci_in(eci_in_filename);
  ECISet eci_min_B = eci_in;
  bool singular;
  int bestsofar = 0;

  MTRand mtrand;

  BP::BP_Comb combs(eci_in.size(), N);

  std::cout << "This calculates all ecisets with " << N << " eci." << std::endl;
  std::cout << "   For N = " << N << " and eci.in size = " << eci_in.size() << " that will be " << combs.total_combs() << " ecisets." << std::endl << std::endl;

  //for( int i=0; i<10; i++)
  //	combs.fix( i, 1);

  int count = 0;

  do {
    eci_in = combs;

    eci_in.fit(corr, DFT_nrg, singular);

    if(eci_min_B.get_cv() == UK || eci_in.get_cv() < eci_min_B.get_cv()) {
      eci_min_B = eci_in;
      bestsofar = combs.get_count();
    }

    if(!singular) {
      std::cout << eci_in.get_bit_string() << "   i: " << combs.get_count() << "\t Nclust: " << eci_in.get_Nclust_on() << " cv: " << eci_in.get_cv() << " rms: " << eci_in.get_rms() << std::endl;
      //std::cout << "i: " << combs.get_count() << " Nclust: " << eci_in.get_Nclust_on() << " cv: " << cv_score << " rms: " << rms_score << " bit_string: " << combs.get_bit_string() << std::endl;
      //std::cout << "   eci: " ;
      //for( int i=0; i<eci_in.size(); i++)
      //	if( eci_in[i].weight != 0)
      //		std::cout << eci_in[i].value << " " ;
      //std::cout << std::endl;
    }
    else {
      std::cout << "error, singular" << std::endl;
    }

    combs.increment();

  }
  while(combs.complete() == false);

  std::cout << eci_min_B.get_bit_string() << "       BestSoFar: " << bestsofar << " Nclust: " << eci_min_B.get_Nclust_on() << " cv: " << eci_min_B.get_cv() << " rms: " << eci_min_B.get_rms() << std::endl;

}

void calc_directmin_eci(int Nrand, int Nmin, int Nmax, std::string energy_filename, std::string eci_in_filename, std::string corr_in_filename, BP::BP_Vec<ECISet> &population) {
  Correlation corr(corr_in_filename);
  EnergySet DFT_nrg(energy_filename);
  ECISet eci_in(eci_in_filename);

  double cv_score;
  double rms_score;
  bool singular;

  ECISet eci_min_A = eci_in;
  ECISet eci_min_B = eci_in;

  MTRand mtrand;

  if(TEST)
    mtrand.seed(1);

  int count = 0;
  bool cont;
  int bestsofar = 1;

  // if no population input, create random population of size Nrand,
  //    if population is input, ignore Nrand
  if(population.size() == 0) {
    for(int i = 0; i < Nrand; i++) {
      population.add(eci_in);
      population[i].randomize(Nmin, mtrand);
    }
  }

  do {
    std::cout << std::endl << "---------------------------" << std::endl;
    // make random start with N clusters turned on
    count++;

    eci_min_A = population[count - 1];

    eci_min_A.fit(corr, DFT_nrg, singular);

    std::cout << eci_min_A.get_bit_string() << "            init: 0" << " Nclust: " << eci_min_A.get_Nclust_on() << " cv: " << eci_min_A.get_cv() << " rms: " << eci_min_A.get_rms() << std::endl;

    eci_min_A = direct_min(Nmin, Nmax, eci_min_A, corr, DFT_nrg, 1);

    std::cout << eci_min_A.get_bit_string() << "           final: " << count << " Nclust: " << eci_min_A.get_Nclust_on() << " cv: " << eci_min_A.get_cv() << " rms: " << eci_min_A.get_rms() << std::endl;

    population[count - 1] = eci_min_A;

    if(eci_min_B.get_cv() == UK || eci_min_A.get_cv() < eci_min_B.get_cv()) {
      eci_min_B = eci_min_A;
      bestsofar = count;
    }

    std::cout << eci_min_B.get_bit_string() << "       bestsofar: " << bestsofar << " Nclust: " << eci_min_B.get_Nclust_on() << " cv: " << eci_min_B.get_cv() << " rms: " << eci_min_B.get_rms() << std::endl;

  }
  while(count < population.size());
  std::cout << std::endl << eci_min_B.get_bit_string() << "       BestSoFar: " << bestsofar << " Nclust: " << eci_min_B.get_Nclust_on() << " cv: " << eci_min_B.get_cv() << " rms: " << eci_min_B.get_rms() << std::endl;

}

void calc_dfsmin_eci(int Nrand, int Nstop, int Nmin, int Nmax, std::string energy_filename, std::string eci_in_filename, std::string corr_in_filename, BP::BP_Vec<ECISet> &population) {
  Correlation corr(corr_in_filename);
  EnergySet DFT_nrg(energy_filename);
  ECISet eci_in(eci_in_filename);

  // set fix so only up to triplets are used
  /*
  if( population.size() == 0)
  {
  	for( int j=0; j<eci_in.size(); j++)
  	{
  		if( j == 2)
  			eci_in.fix(j,0);

  		if( eci_in.get_size(j) > 3)
  		{
  			std::cout << " fix: " << j << " to " << 0 << std::endl;
  			eci_in.fix(j,0);
  		}
  	}
  }
  else
  {
  	for( int k=0; k<population.size(); k++)
  		for( int j=0; j<population[k].size(); j++)
  		{
  			if( j == 2)
  				population[k].fix(j,0);

  			if( population[k].get_size(j) > 3)
  			{
  				if( k == 0) std::cout << " fix: " << j << " to " << 0 << std::endl;
  				population[k].fix(j,0);
  			}
  		}
  }
  */

  // set fix so only clusters less than some size are used
  //{
  //	eci_in.set_Rmax_fix(6.0);
  //
  //}


  double cv_score;
  double rms_score;
  bool singular;

  ECISet eci_min_A = eci_in;
  ECISet eci_min_B = eci_in;

  MTRand mtrand;

  if(TEST)
    mtrand.seed(1);

  int count = 0;
  bool cont;
  int bestsofar = 1;

  if(population.size() == 0) {
    for(int i = 0; i < Nrand; i++) {
      population.add(eci_in);
      population[i].randomize(Nmin, mtrand);
    }
  }

  std::cout << "Nstop: " << Nstop << std::endl;
  do {
    std::cout << std::endl << "---------------------------" << std::endl;
    // make random start with N clusters turned on
    count++;

    eci_min_A = population[count - 1];

    eci_min_A.fit(corr, DFT_nrg, singular);

    std::cout << eci_min_A.get_bit_string() << "            init: 0" << " Nclust: " << eci_min_A.get_Nclust_on() << " cv: " << eci_min_A.get_cv() << " rms: " << eci_min_A.get_rms() << std::endl;

    eci_min_A = dfs_min(Nstop, Nmin, Nmax, eci_min_A, corr, DFT_nrg, 1);

    std::cout << eci_min_A.get_bit_string() << "           final: " << count << " Nclust: " << eci_min_A.get_Nclust_on() << " cv: " << eci_min_A.get_cv() << " rms: " << eci_min_A.get_rms() << std::endl;

    population[count - 1] = eci_min_A;

    if(eci_min_B.get_cv() == UK || eci_min_A.get_cv() < eci_min_B.get_cv()) {
      eci_min_B = eci_min_A;
      bestsofar = count;
    }

    std::cout << eci_min_B.get_bit_string() << "       bestsofar: " << bestsofar << " Nclust: " << eci_min_B.get_Nclust_on() << " cv: " << eci_min_B.get_cv() << " rms: " << eci_min_B.get_rms() << std::endl;

  }
  while(count < population.size());

  std::cout << std::endl << eci_min_B.get_bit_string() << "       BestSoFar: " << bestsofar << " Nclust: " << eci_min_B.get_Nclust_on() << " cv: " << eci_min_B.get_cv() << " rms: " << eci_min_B.get_rms() << std::endl;

}

void calc_ga_eci(int Npopulation, int Nmin, int Nmax, int Nchildren, int Nmutations, std::string energy_filename, std::string eci_in_filename, std::string corr_in_filename, BP::BP_Vec<ECISet> &population) {
  Correlation corr(corr_in_filename);
  EnergySet DFT_nrg(energy_filename);
  ECISet eci_in(eci_in_filename);
  MTRand mtrand;
  bool singular;
  int best;

  // set fix so that only up to triplets are used
  /*
  if( population.size() == 0)
  {
  	for( int j=0; j<eci_in.size(); j++)
  	{
  		if( j == 2)
  			eci_in.fix(j,0);

  		if( eci_in.get_size(j) > 3)
  		{
  			std::cout << " fix: " << j << " to " << 0 << std::endl;
  			eci_in.fix(j,0);
  		}
  	}
  }
  else
  {
  	for( int k=0; k<population.size(); k++)
  		for( int j=0; j<population[k].size(); j++)
  		{
  			if( j == 2)
  				population[k].fix(j,0);

  			if( population[k].get_size(j) > 3)
  			{
  				if( k == 0) std::cout << " fix: " << j << " to " << 0 << std::endl;
  				population[k].fix(j,0);
  			}
  		}
  }
  */

  // create initial random population
  //std::cout << "creating initial random population" << std::endl;
  //BP::BP_Vec<ECISet> population;
  if(population.size() == 0) {
    for(int i = 0; i < Npopulation; i++) {
      population.add(eci_in);
      population[i].randomize(Nmin, mtrand);
    }
  }

  std::cout << std::endl << "Initial Population: " << std::endl;
  for(int i = 0; i < population.size(); i++) {
    population[i].fit(corr, DFT_nrg, singular);
    std::cout << population[i].get_bit_string() << "    i: " << i << " Nclust: " << population[i].get_Nclust_on() << "  cv: " << population[i].get_cv() << " rms: " << population[i].get_rms() << std::endl;
  }
  std::cout << std::endl;

  // create and run genetic algorithm
  //std::cout << "create gac" << std::endl;
  GeneticAlgorithm ga(population, Nmin, Nmax, Nchildren, Nmutations);
  //std::cout << "run ga" << std::endl;
  ga.run(mtrand, corr, DFT_nrg, 1);
  population = ga.get_population();

  // output results
  std::cout << std::endl << "Final Population: " << std::endl;
  for(int i = 0; i < ga.size(); i++) {
    std::cout << ga[i].get_bit_string() << "    i: " << i << " Nclust: " << ga[i].get_Nclust_on() << "  cv: " << ga[i].get_cv() << " rms: " << ga[i].get_rms() << std::endl;
  }

  best = ga.get_best_index();

  std::cout << std::endl << ga[best].get_bit_string() << "   BestSoFar: " << best << " Nclust: " << ga[best].get_Nclust_on() << "  cv: " << ga[best].get_cv() << " rms: " << ga[best].get_rms() << std::endl;

}

void calc_ga_dir_eci(int Npopulation, int Nmin, int Nmax, int Nchildren, int Nmutations, std::string energy_filename, std::string eci_in_filename, std::string corr_in_filename, BP::BP_Vec<ECISet> &population) {
  Correlation corr(corr_in_filename);
  EnergySet DFT_nrg(energy_filename);
  ECISet eci_in(eci_in_filename);
  MTRand mtrand;
  bool singular;
  int best;
  bool add_ok;

  // create initial random population
  if(population.size() == 0) {
    std::cout << "Creating initial random population" << std::endl;
    for(int i = 0; i < Npopulation; i++) {
      do {
        eci_in.randomize(Nmin, mtrand);
        //eci_in = direct_min( Nmin, Nmax, eci_in, corr, DFT_nrg, 1);
        add_ok = add_once(population, eci_in);
      }
      while(!add_ok);
      //population.add(eci_in);
      std::cout << population[i].get_bit_string() << "    Added" << std::endl;
    }
  }

  std::cout << std::endl << "Initial Population: " << std::endl;
  for(int i = 0; i < population.size(); i++) {
    population[i].fit(corr, DFT_nrg, singular);
    std::cout << population[i].get_bit_string() << "    i: " << i << " Nclust: " << population[i].get_Nclust_on() << "  cv: " << population[i].get_cv() << " rms: " << population[i].get_rms() << std::endl;
  }
  std::cout << std::endl;

  // create and run genetic algorithm
  //std::cout << "create gac" << std::endl;
  GeneticAlgorithm ga(population, Nmin, Nmax, Nchildren, Nmutations);
  //std::cout << "run ga" << std::endl;
  ga.run_dir(mtrand, corr, DFT_nrg, 1);
  population = ga.get_population();

  // output results
  std::cout << std::endl << "Final Population: " << std::endl;
  for(int i = 0; i < ga.size(); i++) {
    std::cout << ga[i].get_bit_string() << "    i: " << i << " Nclust: " << ga[i].get_Nclust_on() << "  cv: " << ga[i].get_cv() << " rms: " << ga[i].get_rms() << std::endl;
  }

  best = ga.get_best_index();

  std::cout << std::endl << ga[best].get_bit_string() << "   BestSoFar: " << best << " Nclust: " << ga[best].get_Nclust_on() << "  cv: " << ga[best].get_cv() << " rms: " << ga[best].get_rms() << std::endl;

}

void calc_ga_dfs_eci(int Npopulation, int Nmin, int Nmax, int Nchildren, int Nmutations, int Nstop, std::string energy_filename, std::string eci_in_filename, std::string corr_in_filename, BP::BP_Vec<ECISet> &population) {
  Correlation corr(corr_in_filename);
  EnergySet DFT_nrg(energy_filename);
  ECISet eci_in(eci_in_filename);
  MTRand mtrand;
  bool singular;
  int best;

  // create initial random population
  //std::cout << "creating initial random population" << std::endl;
  if(population.size() == 0) {
    for(int i = 0; i < Npopulation; i++) {
      population.add(eci_in);
      population[i].randomize(Nmin, mtrand);
    }
  }

  std::cout << std::endl << "Initial Population: " << std::endl;
  for(int i = 0; i < population.size(); i++) {
    population[i].fit(corr, DFT_nrg, singular);
    std::cout << population[i].get_bit_string() << "    i: " << i << " Nclust: " << population[i].get_Nclust_on() << "  cv: " << population[i].get_cv() << " rms: " << population[i].get_rms() << std::endl;
  }
  std::cout << std::endl;

  // create and run genetic algorithm
  //std::cout << "create gac" << std::endl;
  GeneticAlgorithm ga(population, Nmin, Nmax, Nchildren, Nmutations);
  //std::cout << "run ga" << std::endl;
  ga.run_dfs(Nstop, mtrand, corr, DFT_nrg, 1);
  population = ga.get_population();

  // output results
  std::cout << std::endl << "Final Population: " << std::endl;
  for(int i = 0; i < ga.size(); i++) {
    std::cout << ga[i].get_bit_string() << "    i: " << i << " Nclust: " << ga[i].get_Nclust_on() << "  cv: " << ga[i].get_cv() << " rms: " << ga[i].get_rms() << std::endl;
  }

  best = ga.get_best_index();

  std::cout << std::endl << ga[best].get_bit_string() << "   BestSoFar: " << best << " Nclust: " << ga[best].get_Nclust_on() << "  cv: " << ga[best].get_cv() << " rms: " << ga[best].get_rms() << std::endl;

}

void weight_nrg(double A, double B, double kT, std::string energy_filename, std::string out_filename) {
  EnergySet nrg_set(energy_filename);

  nrg_set.weight_nrg(A, B, kT);

  nrg_set.write(out_filename, "default");
}

void weight_EMin(double A, double B, double kT, std::string energy_filename, std::string out_filename) {
  EnergySet nrg_set(energy_filename);
  double minEn = nrg_set.find_energy_min();
  nrg_set.weight_ERef(A, B, kT, minEn);
  nrg_set.write(out_filename, "default");
}

void weight_ERef(double A, double B, double ERef, double kT, std::string energy_filename, std::string out_filename) {
  EnergySet nrg_set(energy_filename);
  nrg_set.weight_ERef(A, B, kT, ERef);
  nrg_set.write(out_filename, "default");
}

////----------------------------
/// test functions

//eci_search -nrg_diff   tot_corr.index.A    AB_structs.energy   eci.in.tot   AB_structs.corr.in.tot   A_structs.energy  A_structs.eci.in  A_structs.corr.in
void calc_nrg_diff(std::string A_corr_index_filename, std::string AB_structs_energy_filename, std::string eci_in_tot_filename, std::string AB_structs_corr_tot_filename,
                   std::string A_structs_energy_filename, std::string A_structs_eci_in_filename, std::string A_structs_corr_in_filename) {
  unsigned long int i, j;
  Correlation AB_structs_corr_tot(AB_structs_corr_tot_filename);
  EnergySet AB_structs_energy(AB_structs_energy_filename);
  ECISet eci_in_tot(eci_in_tot_filename);

  bool singular;
  Correlation A_structs_corr_in(A_structs_corr_in_filename);
  EnergySet A_structs_energy(A_structs_energy_filename);
  ECISet A_structs_eci_in(A_structs_eci_in_filename);
  A_structs_eci_in.fit(A_structs_corr_in,  A_structs_energy, singular);

  A_structs_eci_in.write_ECIout("A_structs.eci.out", A_structs_energy, "default");

  // read  tot_corr.index.A file, store in A_corr_index_list, and create AB_corr_index_list
  BP::BP_Parse index_file(A_corr_index_filename);
  BP::BP_Vec<int> A_corr_index_list;
  BP::BP_Vec<int> AB_corr_index_list;
  BP::BP_Vec<int> i_list;
  do {
    i_list = index_file.getline_int();
    if(i_list.size() == 1) {
      A_corr_index_list.add(i_list[0]);
    }
  }
  while(index_file.eof() == false);
  //std::cout << "A_corr_index_list: " << A_corr_index_list << std::endl;

  for(i = 0; i < eci_in_tot.size(); i++) {
    if(!A_corr_index_list.find_first(i, j)) {
      AB_corr_index_list.add(i);
    }
  }
  //std::cout << "AB_corr_index_list: " << AB_corr_index_list << std::endl;

  A_structs_eci_in.write_ECIout(A_corr_index_list, "A_structs.eci.out.part", A_structs_energy, "default");

  // calculate and write the energy_diff: AB_DFT_nrg_diff = AB_DFT_nrg - CLEX(AB_strucs.A_corr, A_eci);
  EnergySet AB_structs_A_nrg_diff = AB_structs_energy;
  AB_structs_A_nrg_diff.calc_nrg_diff(A_corr_index_list, AB_structs_corr_tot, A_structs_eci_in);
  AB_structs_A_nrg_diff.write("AB_structs.energy.diff", "default");

  // get subset of AB_structs.corr.in.tot corresponding to AB_strucs.corr.in.AB (AB clusters)
  Correlation AB_structs_corr_AB = AB_structs_corr_tot;
  AB_structs_corr_AB.cluster_subset(AB_corr_index_list);
  AB_structs_corr_AB.write("AB_structs.corr.in.AB", "default");

  Correlation AB_structs_corr_A = AB_structs_corr_tot;
  AB_structs_corr_A.cluster_subset(A_corr_index_list);
  AB_structs_corr_A.write("AB_structs.corr.in.A", "default");

  // write eci.in.AB
  eci_in_tot.write_ECIin(AB_corr_index_list, "eci.in.AB", "default");
}

////----------------------------
/// sub functions
ECISet direct_min(int Nmin, int Nmax, const ECISet &eci_start, Correlation &corr, EnergySet &nrg_set, bool print_steps) {
  ECISet eci_min_A = eci_start;
  ECISet eci_in;
  bool cont;
  bool singular;
  int Nchoice;
  double last_cv;
  int last_flip;
  // make sure eci_min_A cv score is known
  eci_min_A.fit(corr, nrg_set, singular);

  // minimize
  int j_count = 1;
  do {
    // initialize for minimization step
    cont = 0;
    eci_in = eci_min_A;
    Nchoice = 0;
    last_cv = eci_min_A.get_cv();
    last_flip = -1;

    // try toggling each cluster on/off
    for(int i = 0; i < eci_in.size(); i++) {
      eci_in.toggle_clust(i);

      // check that the eci_set meets fix and Nmin/Nmax criteria
      if(eci_in.fix_ok())
        if(eci_in.get_Nclust_on() >= Nmin && eci_in.get_Nclust_on() <= Nmax) {
          // find fit/cv score
          eci_in.fit(corr, nrg_set, singular);

          if(eci_in.get_cv() < last_cv) {
            Nchoice++;
          }

          // keep the steepest descent option
          if(eci_in.get_cv() < eci_min_A.get_cv()) {

            cont = 1;
            last_flip = i;
            eci_min_A = eci_in;
          }
        }

      eci_in.toggle_clust(i);
    }


    if(print_steps) std::cout << eci_min_A.get_bit_string() << "      minimizing: " << j_count << " Nclust: " << eci_min_A.get_Nclust_on() << " cv: " << eci_min_A.get_cv() << " rms: " << eci_min_A.get_rms() << " Nchoice: " << Nchoice << " last_flip: " << last_flip << std::endl;

    j_count++;


  }
  while(cont);

  return eci_min_A;
}

ECISet dfs_min(long int Nstop, int Nmin, int Nmax, const ECISet &eci_start, Correlation &corr, EnergySet &nrg_set, bool print_steps) {
  // search for the minimum by depth-first search
  // minimize
  ECISet eci_min_A = eci_start;
  ECISet eci_min_B;
  ECISet eci_in;
  //BP::BP_RVG_tree<ECISet> generator;
  BP::BP_GVec<ECISet> eci_list;
  BP::BP_Vec<std::string> bit_string_list;
  int count_since_last_best = 0;
  int Nchoice;
  int bestsofar;
  bool singular;
  BP::BP_GVec_Member<ECISet> *s_member;

  // make sure eci_min_A cv score is known
  eci_min_A.fit(corr, nrg_set, singular);
  eci_min_B = eci_min_A;

  // depth-first search
  int dfs_count = 0;
  do {
    //std::cout << std::endl << "  -+++-----------------------" << std::endl;


    // find all changes that reduce the cv score, and add to RVG_tree with rate = (init_cv - curr_cv)
    eci_in = eci_min_A;
    Nchoice = 0;
    for(int i = 0; i < eci_in.size(); i++) {
      eci_in.toggle_clust(i);

      if(eci_in.fix_ok())
        if(eci_in.get_Nclust_on() >= Nmin && eci_in.get_Nclust_on() <= Nmax) {
          //std::cout << "toggle i: " << i << "\n";
          if(add_once(bit_string_list, eci_in.get_bit_string())) {
            eci_in.fit(corr, nrg_set, singular);

            if(eci_in.get_cv() < eci_min_A.get_cv()) {
              Nchoice++;
              //generator.add( eci_list.add(eci_in), exp(1.0/eci_in.get_cv()));
              eci_list.add(eci_in);
            }
          }
        }

      eci_in.toggle_clust(i);

    }




    if(eci_list.size() > 0) {
      // pick a bitstring from the RVG_tree
      //s_member = generator.pick(mtrand.randExc());

      // pick the bitstring with lowest cv score
      s_member = get_best_member(eci_list);

      eci_min_A = eci_list[ s_member];
      eci_list.remove(s_member);

      if(eci_min_B.get_cv() == UK || eci_min_A.get_cv() < eci_min_B.get_cv()) {
        eci_min_B = eci_min_A;
        bestsofar = dfs_count;
        count_since_last_best = 0;
      }

      if(print_steps) std::cout << eci_min_A.get_bit_string() << "      DFSchoice: " << dfs_count << " Nclust: " << eci_min_A.get_Nclust_on() << " cv: " << eci_min_A.get_cv() << " rms: " << eci_min_A.get_rms() << " Nchoice: " << Nchoice << " listsize: " << eci_list.size() << " Count: " << count_since_last_best << "/" << Nstop << "  bestofDFS: " << bestsofar << "  best_cv: " << eci_min_B.get_cv() << std::endl;

      dfs_count++;
    }

    count_since_last_best++;

  }
  while((eci_list.size() > 0) && (count_since_last_best < Nstop || Nstop == 0));

  return eci_min_B;
}


////----------------------------
/// compressive sensing functions
ECISet calc_FPC_eci(double mu, double prec_shrink, double prec_ECI, const ECISet &eci_set, const Correlation &corr, const EnergySet &nrg_set, bool print_steps, unsigned long int max_step) {
  std::cout << "begin calc_FPC_eci()" << std::endl;

  // method:
  //		Start: ECI_0 = 0 vector
  //		Then:
  //			ECI_i+1 = shrink( ECI_i - tau*g_i, mu*tau)
  //		Where:
  //			g_i = Corr_transpose*(Corr*ECI_i - E) = M1*ECI_i - V1
  //			shrink(y,a) = sign(y)*max( fabs(y) - a, 0)
  //			tau = min( 1.999, -1.665*Corr.rows()/Corr.cols() + 2.665)
  //		Stop when:
  //			max(g)/mu - 1 < prec_shrink
  //		and
  //			2norm( ECI_i+1 - ECI_i)/2norm(ECI_i) < prec_ECI


  unsigned long int i, j, k, ii, jj, kk;
  ECISet eci_out = eci_set;

  unsigned long int Nnrg = nrg_set.get_Nstruct_on();
  unsigned long int Neci = eci_set.size();

  Eigen::MatrixXd C(Nnrg, Neci);					// Correlation matrix
  Eigen::MatrixXd Cn(Nnrg, Neci);				// normalized correlation matrix so that max eigenvalue of Cn.transpose*Cn <= 1
  Eigen::MatrixXd M1(Neci, Neci);				// Cn.transpose()*Cn;
  Eigen::VectorXd E(Nnrg);						// Enegry vector
  Eigen::VectorXd En(Nnrg);						// normalized energy vector
  Eigen::VectorXd V1(Neci);						// Cn.transpose()*En
  Eigen::VectorXd ECI = Eigen::VectorXd::Zero(Neci);	// current solution

  double tau = std::min(1.999, std::max(1.0, -1.665 * (1.0 * Nnrg) / (1.0 * Neci) + 2.665));
  double rms;

  //std::cout << "Nnrg: " << Nnrg << std::endl;
  //std::cout << "Neci: " << Neci << std::endl;
  //std::cout << "tau: " << tau << std::endl;

  // set Correlation matrix
  //std::cout << "set C" << std::endl;
  ii = 0;
  for(i = 0; i < nrg_set.size(); i++)
    if(nrg_set.get_weight(i) != 0) {
      jj = 0;
      for(j = 0; j < eci_set.size(); j++) {
        C(ii, jj) = nrg_set.get_weight(i) * corr[i][j];		// include weight!?
        jj++;
      }
      ii++;
    }

  // set Energy vector
  //std::cout << "set E" << std::endl;
  ii = 0;
  for(i = 0; i < nrg_set.size(); i++)
    if(nrg_set.get_weight(i) != 0) {
      E(ii) = nrg_set.get_weight(i) * nrg_set.get_Ef(i);		// include weight!?
      ii++;
    }

  // normalize C and E, so that largest eigenvalue of C.transpose*C is <= 1
  M1 = C.transpose() * C;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(M1);
  if(eigensolver.info() != Eigen::Success) {
    std::cout << "SelfAdjointEigenSolver failed!" << std::endl;
    exit(1);
  };
  //std::cout << "The eigenvalues of M1 are:\n" << eigensolver.eigenvalues() << std::endl;

  double max_eigenvalue = eigensolver.eigenvalues().maxCoeff();
  double a1 = sqrt(1.1 * max_eigenvalue);
  //std::cout << "The max eigenvalue of M1 is:\n" << max_eigenvalue << std::endl;

  // set normalized C & E
  Cn = C / a1;
  En = E / a1;

  //std::cout << "Cn: " << Cn << std::endl;
  //std::cout << "En: " << En << std::endl;


  // set M1 = corr.transpose() * corr
  //std::cout << "set M1" << std::endl;
  M1 = Cn.transpose() * Cn;

  // set V1 = corr.transpose() * nrg
  //std::cout << "set V1" << std::endl;
  V1 = Cn.transpose() * En;
  //std::cout << "V1: " << V1 << std::endl;

  // do Fixed-Point continuation algorithm to find ECI
  FPC(M1, V1, ECI, mu, tau, prec_shrink, prec_ECI, print_steps, max_step);

  //std::cout << "Final ECI:" << std::endl;
  //std::cout << ECI << std::endl;

  eci_out.set_values_and_weights(ECI);
  eci_out.set_cv(1e20);
  rms = (C * ECI - E).norm() / sqrt(1.0 * Nnrg);
  eci_out.set_rms(rms);

  std::cout << "finish calc_FPC_eci()" << std::endl;
  return eci_out;


};

ECISet calc_BI_eci(double mu, double prec_shrink, double prec_ECI, const ECISet &eci_set, const Correlation &corr, const EnergySet &nrg_set, bool print_steps, unsigned long int max_step) {
  std::cout << "begin calc_BI_eci()" << std::endl;

  // method:
  //		Start: ECI_0 = 0 vector, F_0 = 0 vector
  //		Then:
  //			F_i+1 = E + F_i - Corr*ECI_i
  //			ECI_i+1 = FPC( Corr, E, mu, tau)
  //		Where:
  //			tau = min( 1.999, -1.665*Corr.rows()/Corr.cols() + 2.665) // same as FPC
  //		Stop when:
  //			2norm( ECI_i+1 - ECI_i)/2norm(ECI_i) < prec_ECI


  unsigned long int i, j, k, ii, jj, kk;
  ECISet eci_out = eci_set;

  unsigned long int Nnrg = nrg_set.get_Nstruct_on();
  unsigned long int Neci = eci_set.size();

  Eigen::MatrixXd C(Nnrg, Neci);					// Correlation matrix
  Eigen::MatrixXd Cn(Nnrg, Neci);				// normalized correlation matrix so that max eigenvalue of Cn.transpose*Cn <= 1
  Eigen::MatrixXd M1(Neci, Neci);				// Cn.transpose()*Cn;
  Eigen::VectorXd E(Nnrg);						// Enegry vector
  Eigen::VectorXd En(Nnrg);						// normalized energy vector
  Eigen::VectorXd ECI = Eigen::VectorXd::Zero(Neci);	// current solution

  double tau = std::min(1.999, std::max(1.0, -1.665 * (1.0 * Nnrg) / (1.0 * Neci) + 2.665));
  double rms;

  //std::cout << "Nnrg: " << Nnrg << std::endl;
  //std::cout << "Neci: " << Neci << std::endl;
  //std::cout << "tau: " << tau << std::endl;

  // set Correlation matrix
  //std::cout << "set C" << std::endl;
  ii = 0;
  for(i = 0; i < nrg_set.size(); i++)
    if(nrg_set.get_weight(i) != 0) {
      jj = 0;
      for(j = 0; j < eci_set.size(); j++) {
        C(ii, jj) = nrg_set.get_weight(i) * corr[i][j];		// include weight!?
        jj++;
      }
      ii++;
    }

  // set Energy vector
  //std::cout << "set E" << std::endl;
  ii = 0;
  for(i = 0; i < nrg_set.size(); i++)
    if(nrg_set.get_weight(i) != 0) {
      E(ii) = nrg_set.get_weight(i) * nrg_set.get_Ef(i);		// include weight!?
      ii++;
    }

  // normalize C and E, so that largest eigenvalue of C.transpose*C is <= 1
  M1 = C.transpose() * C;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(M1);
  if(eigensolver.info() != Eigen::Success) {
    std::cout << "SelfAdjointEigenSolver failed!" << std::endl;
    exit(1);
  };
  //std::cout << "The eigenvalues of M1 are:\n" << eigensolver.eigenvalues() << std::endl;

  double max_eigenvalue = eigensolver.eigenvalues().maxCoeff();
  double a1 = sqrt(1.1 * max_eigenvalue);
  //std::cout << "The max eigenvalue of M1 is:\n" << max_eigenvalue << std::endl;

  // set normalized C & E
  Cn = C / a1;
  En = E / a1;

  //std::cout << "Cn: " << Cn << std::endl;
  //std::cout << "En: " << En << std::endl;


  // set M1 = corr.transpose() * corr
  //std::cout << "set M1" << std::endl;
  M1 = Cn.transpose() * Cn;

  // set V1 = corr.transpose() * nrg
  //std::cout << "set V1" << std::endl;
  //V1 = Cn.transpose()*En;
  //std::cout << "V1: " << V1 << std::endl;

  // do Bergman Iteration algorithm to find ECI
  BI(Cn, En, M1, ECI, mu, tau, prec_shrink, prec_ECI, print_steps, max_step);

  //std::cout << "Final ECI:" << std::endl;
  //std::cout << ECI << std::endl;

  eci_out.set_values_and_weights(ECI);
  eci_out.set_cv(1e20);
  rms = (C * ECI - E).norm() / sqrt(1.0 * Nnrg);
  eci_out.set_rms(rms);

  std::cout << "finish calc_FPC_eci()" << std::endl;
  return eci_out;


};

void BI(const Eigen::MatrixXd &Cn, const Eigen::VectorXd &En, const Eigen::MatrixXd &M1, Eigen::VectorXd &ECI, double mu, double tau, double prec_shrink, double prec_ECI, bool print_steps, unsigned long int max_step) {
  // Bergman iteration for fitting ECI that minimize L1 norm

  Eigen::VectorXd ECI_i;
  Eigen::VectorXd F = Eigen::VectorXd::Zero(En.size());
  Eigen::VectorXd V1(ECI.size());						// Cn.transpose()*En

  unsigned long int step = 0;
  bool cont = true;
  double dECI, rms;

  do {
    step++;

    ECI_i = ECI;
    F = En + F - Cn * ECI;
    V1 = Cn.transpose() * F;
    FPC(M1, V1, ECI, mu, tau, prec_shrink, prec_ECI, false, max_step);

    dECI = (ECI - ECI_i).norm() / ECI_i.norm();

    if(print_steps) {

      if(dECI != dECI) {
        throw std::runtime_error("Error in BI iteration: dECI = nan\nTry a smaller mu value.");
      }

      if(step == 1) {
        std::cout << "    Step " << step << std::endl;
      }
      else
        std::cout << "    Step " << step << "  dECI: " << dECI << std::endl;
    }

    if(dECI < prec_ECI)
      cont = false;

  }
  while(cont);

}

void FPC(const Eigen::MatrixXd &M1, const Eigen::VectorXd &V1, Eigen::VectorXd &ECI, double mu, double tau, double prec_shrink, double prec_ECI, bool print_steps, unsigned long int max_step) {
  //  Fixed-Point continuation algorithm for fitting ECI that minimize L1 norm
  //
  //

  //std::cout << "begin FPC()" << std::endl;

  unsigned long int step = 0;
  bool cont = true;
  double ECI_mag, delta_mag;

  Eigen::VectorXd G = Eigen::VectorXd::Zero(ECI.size());		// gradient of 2norm

  do {
    step++;

    //std::cout << " :1" << std::endl;
    ECI_mag = ECI.norm();
    //std::cout << " :2" << std::endl;
    G = M1 * ECI - V1;
    //std::cout << " :3" << std::endl;
    delta_mag = shrink(ECI, G, mu, tau);

    if(print_steps) {
      if(step % 1000 == 0) {

        double _shrink = G.maxCoeff() / mu - 1.0;
        double _dECI = (delta_mag / ECI_mag);

        if(_dECI != _dECI) {
          throw std::runtime_error("Error in FPC iteration: dECI = nan\nTry a smaller mu value.");
        }
        if(std::isinf(_shrink)) {
          throw std::runtime_error("Error in FPC iteration: shrink = inf\nTry a larger mu value.");
        }

        //rms = (C*ECI - E).norm()/sqrt(1.0*Nnrg);
        //std::cout << ECI << std::endl;
        //std::cout << "  Step " << step << "  rms: " << rms << "  shrink: " << G.maxCoeff()/mu - 1.0 << "  dECI: " << (delta_mag/ECI_mag) << std::endl;
        std::cout << "  Step " << step << "  shrink: " << _shrink << "  dECI: " << _dECI << std::endl;

      }
    }

    //std::cout << " :4" << std::endl;
    if(G.maxCoeff() / mu - 1.0 < prec_shrink) {
      if((delta_mag / ECI_mag) < prec_ECI)
        cont = false;
    }

    if(step > max_step) {
      std::cout << "STOPPING: Reached FPC iteration max step: " << max_step << std::endl;
      cont = false;
    }

    //BP::BP_pause();

  }
  while(cont);

};

double shrink(Eigen::VectorXd &ECI, const Eigen::VectorXd &G, double mu, double tau) {
  //std::cout << "begin shrink()" << std::endl;
  //	perform:
  //		ECI_i+1 = shrink( ECI_i - tau*g_i, mu*tau)
  //		shrink(y,a) = sign(y)*max( fabs(y) - a, 0)
  //	return:
  //		2norm(ECI_i+1 - ECI_i)

  double sqr_sum = 0.0;
  double y, a, init, delta;
  unsigned long int i;

  a = mu * tau;
  //std::cout << "a: " << a << std::endl;
  for(i = 0; i < ECI.size(); i++) {
    init = ECI(i);
    y = ECI(i) - tau * G(i);
    ECI(i) = BP::sign(y) * std::max<double>(fabs(y) - a, 0.0);
    sqr_sum += BP::sqr(ECI(i) - init);

    //std::cout << init << " " << ECI(i) << "  y: " << y <<  " delta: " << ECI(i) - init << std::endl;

  }

  //std::cout << "finish shrink()" << std::endl;
  //BP::BP_pause();
  return sqrt(sqr_sum);

};



////----------------------------
/// misc.
bool is_bitstring(std::string s) {

  for(int i = 0; i < s.size(); i++) {
    if(s[i] != '0' && s[i] != '1') {
      return false;
    }
  }
  return true;
};

BP::BP_Vec<ECISet> read_bit_strings_file(ECISet eci_in, string bit_strings_filename) {

  std::string format = get_format_from_ext(bit_strings_filename);

  if(format == "text") {
    //cout << "begin read_bit_strings_file()" << endl;
    BP::BP_Parse file(bit_strings_filename);

    ECISet eci_A = eci_in;
    BP::BP_Vec<ECISet> eci_list;
    BP::BP_Vec<string> entry;
    // data
    do {
      entry = file.getline_string();
      if(entry.size() != 0) {
        // check if the first entry is a bit_string
        if(is_bitstring(entry[0])) {
          eci_A.set_bit_string(entry[0]);
          eci_list.add(eci_A);
        }

      }
    }
    while(file.eof() == false);


    //cout << "Population from input file: " << endl;
    //for( int i=0; i<eci_list.size(); i++)
    //{
    //	cout << eci_list[i].get_bit_string() << "  i: " << i << endl;
    //}

    //cout << "finish read_bit_strings_file()" << endl;
    return eci_list;
  }
  else if(format == "json") {
    CASM::jsonParser json(bit_strings_filename);

    ECISet eci_A = eci_in;
    BP::BP_Vec<ECISet> eci_list;
    for(int i = 0; i < json["population"].size(); i++) {
      eci_A.set_bit_string(json["population"][i]["eciset"].get<std::string>());
      eci_list.add(eci_A);
    }

    //cout << "Population from input file: " << endl;
    //for( int i=0; i<eci_list.size(); i++)
    //{
    //	cout << eci_list[i].get_bit_string() << "  i: " << i << endl;
    //}

    //cout << "finish read_bit_strings_file()" << endl;
    return eci_list;
  }
  else {
    std::cout << "Unexpected format option for read_bit_strings_file" << std::endl;
    std::cout << "  Expected 'text' or 'json', but received: " << format << std::endl;
    exit(1);
  }
};

double max_cv(BP::BP_GVec<ECISet> &population, BP::BP_GVec_Member<ECISet> *&worst_parent) {
  double max = -1;
  for(int i = 0; i < population.size(); i++) {
    if(population[i].get_cv() > max) {
      worst_parent = population.member(i);
      max = population[i].get_cv();
    }
  }

  return max;
};

double min_cv(BP::BP_GVec<ECISet> &population, BP::BP_GVec_Member<ECISet> *&worst_parent) {
  double max = -1;
  double min = 1.0e20;
  for(int i = 0; i < population.size(); i++) {
    if(population[i].get_cv() < min) {
      //worst_parent = population.member(i);
      min = population[i].get_cv();
    }

    if(population[i].get_cv() > max) {
      worst_parent = population.member(i);
      max = population[i].get_cv();
    }
  }

  return min;
};

int get_best_index(const BP::BP_Vec<ECISet> &population) {
  int best;
  double min = 1.0e20;
  for(int i = 0; i < population.size(); i++) {
    if(population[i].get_cv() < min) {
      best = i;
      min = population[i].get_cv();
    }
  }

  return best;
};

int get_best_index(const BP::BP_GVec<ECISet> &population) {
  int best;
  double min = 1.0e20;
  for(int i = 0; i < population.size(); i++) {
    if(population[i].get_cv() < min) {
      best = i;
      min = population[i].get_cv();
    }
  }

  return best;
};

BP::BP_GVec_Member<ECISet> *get_best_member(BP::BP_GVec<ECISet> &population) {
  int best;
  double min = 1.0e20;
  for(int i = 0; i < population.size(); i++) {
    if(population[i].get_cv() < min) {
      best = i;
      min = population[i].get_cv();
    }
  }

  return population.member(best);
};

// randomly pick two parents from gene_pool
//   mate them and set the bit_string in child
void mate(const BP::BP_Vec<ECISetState> &gene_pool, ECISet &child, int Nmut, MTRand &mtrand) {
  int Mom, Dad, j;

  if(gene_pool.size() <= 1) {
    std::cout << "Error in mate(): gene_pool.size() <= 1" << std::endl;
    exit(1);
  }

  Mom = mtrand.randInt(gene_pool.size() - 1);
  do {
    Dad = mtrand.randInt(gene_pool.size() - 1);
  }
  while(Dad == Mom);

  child.set_bit_string(gene_pool[Mom].bit_string);

  for(int i = 0; i < child.size(); i++) {

    if(gene_pool[Mom].bit_string[i] != gene_pool[Dad].bit_string[i]) {
      if(mtrand.randExc() < 0.5)
        child.set_clust_off(i);
      else
        child.set_clust_on(i);
    }
  }

  // mutate (random with probability 1/Nmutations)
  for(int i = 0; i < child.size(); i++) {
    if(mtrand.randExc() < (1.0 * Nmut) / child.size()) {
      child.toggle_clust(i);
      if(!child.fix_ok())
        child.toggle_clust(i);
    }
  }

  // enforce Nmin limit
  while(child.get_Nclust_on() < child.get_Nmin()) {
    j = mtrand.randInt(child.size() - 1);
    child.set_clust_on(j);
    if(!child.fix_ok())
      child.set_clust_off(j);
  }

  // enforce Nmax limit
  while(child.get_Nclust_on() > child.get_Nmax()) {
    j = mtrand.randInt(child.size() - 1);
    child.set_clust_off(j);
    if(!child.fix_ok())
      child.set_clust_on(j);
  }

}

std::string unique(std::string filename) {
  BP::BP_Dir dir;
  if(!dir.is_file(filename) && !dir.is_file(json_ext(filename))) {
    return filename;
  }
  int i = 0;
  while(dir.is_file(filename + "." + BP::itos(i)) || dir.is_file(json_ext(filename + "." + BP::itos(i)))) {
    i++;
  }
  return filename + "." + BP::itos(i);
}

std::string json_ext(std::string filename) {
  if(filename.size() <= 5 || filename.compare(filename.size() - 5, 5, ".json") != 0) {
    return filename + ".json";
  }
  return filename;
}

std::string rm_json_ext(std::string filename) {
  if(filename.size() > 5 && filename.compare(filename.size() - 5, 5, ".json") == 0) {
    return filename.substr(0, filename.size() - 5);
  }
  return filename;
}

std::string get_format_from_ext(std::string filename) {
  if(rm_json_ext(filename) == filename)
    return "text";
  return "json";
}


#endif // Functions_CC
