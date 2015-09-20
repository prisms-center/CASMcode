#ifndef BP_useful_HH
#define BP_useful_HH

#include <iostream>
#include <math.h>
#include "casm/BP_C++/BP_Vec.hh"
#include "casm/BP_C++/BP_basic.hh"
#include "casm/BP_C++/BP_Parse.hh"
#include "casm/BP_C++/BP_StopWatch.hh"
#include "casm/external/MersenneTwister/MersenneTwister.h"

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
///  Useful functions.  These use BP:: functions, but no BP:: class depends on these

namespace BP {

  //int min(int a, int b);
  //unsigned long int min(unsigned long int a, unsigned long int b);
  //int max(int a, int b);
  //unsigned long int max(unsigned long int a, unsigned long int b);

  template<class T> int sign(T val) {
    return (T(0) < val) - (val < T(0));
  }

  ///		return greatest common factor, using Euclid alg
  unsigned long int gcf(BP_Vec<unsigned long int> i_list);

  ///		return greatest common factor, using Euclid alg
  long int gcf(BP_Vec<long int> i_list);

  ///		return greatest common factor, using Euclid alg
  int gcf(BP_Vec<int> i_list);

  ///		return least common multiple, using lcm(a, b, c) = lcm( lcm(a,b),c) & gcf(a,b)*lcm(a,b) = a*b
  unsigned long int lcm(BP_Vec<unsigned long int> i_list);

  ///		return all factors of a number
  BP_Vec<unsigned long int> factors(unsigned long int n);

  ///		return prime factorization of a number
  BP_Vec<unsigned long int> prime_factorization(unsigned long int n);

  ///		return true if prime, false if not
  bool is_prime(unsigned long int n);

  ///		seed a MTRand random number generator from a file of seeds
  void read_rand_seed(int seed_index, std::string seed_list, MTRand &mtrand);

  ///		write to a stream a list of random numbers which can be used to seed a MTRand object
  void write_rand_seed(int N, std::ostream &sout);

  ///		write the current state of the generator to a file
  void write_rand_state(std::string rand_filename, MTRand &mtrand);

  ///		read the current state of the generator from a file
  bool read_rand_state(std::string rand_filename, MTRand &mtrand);

}

#endif // BP_useful_HH
