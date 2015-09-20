#ifndef BP_useful_CC
#define BP_useful_CC

#include "casm/BP_C++/BP_useful.hh"

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
///  Useful functions

namespace BP {

  //int min(int a, int b) {
  //  if(a < b) return a;
  //  else return b;
  //};
  //unsigned long int min(unsigned long int a, unsigned long int b) {
  //  if(a < b) return a;
  //  else return b;
  //};
  //int max(int a, int b) {
  //  if(a > b) return a;
  //  else return b;
  //};
  //unsigned long int max(unsigned long int a, unsigned long int b) {
  //  if(a > b) return a;
  //  else return b;
  //};

  ///		return greatest common factor, using Euclid alg
  unsigned long int gcf(BP_Vec<unsigned long int> i_list) {

    while(i_list.size() > 1) {
      // find gcf of i_list[0] & i_list[1]
      if(i_list[0] == 0) {
        i_list.remove(0);
      }
      else {
        while(i_list[1] != 0) {
          if(i_list[0] > i_list[1])
            i_list[0] = i_list[0] - i_list[1];
          else
            i_list[1] = i_list[1] - i_list[0];

        }
        i_list.remove(1);
      }

    }

    return i_list[0];
  };

  ///		return greatest common factor, using Euclid alg
  long int gcf(BP_Vec<long int> i_list) {
    // make all positive
    i_list = abs(i_list);

    while(i_list.size() > 1) {
      // find gcf of i_list[0] & i_list[1]
      if(i_list[0] == 0) {
        i_list.remove(0);
      }
      else {
        while(i_list[1] != 0) {
          if(i_list[0] > i_list[1])
            i_list[0] = i_list[0] - i_list[1];
          else
            i_list[1] = i_list[1] - i_list[0];

        }
        i_list.remove(1);
      }

    }

    return i_list[0];
  };

  ///		return greatest common factor, using Euclid alg
  int gcf(BP_Vec<int> i_list) {
    // make all positive
    i_list = abs(i_list);

    while(i_list.size() > 1) {
      // find gcf of i_list[0] & i_list[1]
      if(i_list[0] == 0) {
        i_list.remove(0);
      }
      else {
        while(i_list[1] != 0) {
          if(i_list[0] > i_list[1])
            i_list[0] = i_list[0] - i_list[1];
          else
            i_list[1] = i_list[1] - i_list[0];

        }
        i_list.remove(1);
      }

    }

    return i_list[0];
  };

  ///		return least common multiple, using lcm(a, b, c) = lcm( lcm(a,b),c) & gcf(a,b)*lcm(a,b) = a*b
  unsigned long int lcm(BP_Vec<unsigned long int> i_list) {
    if(i_list.size() == 0)
      return 0;
    else if(i_list.size() == 1)
      return i_list[0];
    else if(i_list.size() == 2)
      return (i_list[0] * i_list[1]) / gcf(i_list);
    else if(i_list.size() > 2) {
      BP_Vec<unsigned long int> new_list;
      new_list.add(lcm(i_list.range(0, 1)));
      new_list.append(i_list.range(2, i_list.size() - 1));
      return lcm(new_list);
    }
    return -1;
  }

  ///		return all factors of a number
  BP_Vec<unsigned long int> factors(unsigned long int n) {
    BP_Vec< unsigned long int> _factors;
    _factors.add(1);
    _factors.add(n);
    unsigned long int max = n;
    unsigned long int q;
    for(unsigned long int i = 2; i < max; i++) {
      q = n / i;
      if((i * q) == n) {
        _factors.add(i);
        if(q != i) _factors.add(q);
        max = q;
      }
    }

    _factors.sort();
    return _factors;
  }

  ///		return prime factorization of a number
  BP_Vec<unsigned long int> prime_factorization(unsigned long int n) {
    BP_Vec< unsigned long int> prime_factors;

    if(n == 0) return prime_factors;
    else if(n == 1) return prime_factors;

    unsigned long int q;
    bool cont = true;

    for(unsigned long int i = 2; i <= n; i++) {
      do {
        q = n / i;

        if((i * q) == n) {
          prime_factors.add(i);
          n = q;
          cont = true;
        }
        else {
          cont = false;
        }
      }
      while(cont);
    }

    return prime_factors;
  }

  ///		return true if prime, false if not
  bool is_prime(unsigned long int n) {
    if(n == 0)
      return false;
    else if(n == 1)
      return false;
    else if(prime_factorization(n).size() == 1)
      return true;
    else
      return false;
  }

  ///		seed a MTRand random number generator from a file of seeds
  void read_rand_seed(int seed_index, std::string seed_list, MTRand &mtrand) {
    //std::cout << "begin read_rand_seed()" << std::endl;


    if(seed_list == "") {
      std::cout << "Error, 'seed_index = " << seed_index << "', but no 'seed_list' given." << std::endl;
      exit(1);
    }

    if(seed_index < 0) {
      std::cout << "Error. 'seed_index = '" << seed_index << "' < 0." << std::endl;
      exit(1);
    }

    BP_Parse file(seed_list);
    MTRand::uint32 seed_vals[ MTRand::N ];


    unsigned long int n1 = 0, n2 = 0, number;
    while(n1 < seed_index && file.get_istream()) {
      file >> number;
      n1++;
    }

    while(file >> number && n2 < MTRand::N) {
      seed_vals[n2] = number;
      //std::cout << "n2: " << n2 << "  val: " << number << std::endl;
      n2++;
    }

    if(n2 < MTRand::N) {
      std::cout << "Error. 'seed_index = '" <<  seed_index << "' is too large." << std::endl;
      std::cout << "  In 'seed_list = " <<  seed_list << "' the maximum possible seed_index is " << n1 + n2 - MTRand::N << std::endl;
      exit(1);
    }

    mtrand.seed(seed_vals);


    //std::cout << "finish read rand seed" << std::endl;
    //BP_pause();
  }

  ///		write to a stream a list of random numbers which can be used to seed a MTRand object
  void write_rand_seed(int N, std::ostream &sout) {
    BP_StopWatch timer;
    MTRand mtrand((unsigned long int) floor(timer.gettime()));
    for(int i = 0; i < N + MTRand::N; i++)
      sout << mtrand.randInt() << " ";
    sout << std::endl;
  }

  ///		write the current state of the generator to a file
  void write_rand_state(std::string rand_filename, MTRand &mtrand) {
    BP_Write file(rand_filename);
    file.newfile();
    file << mtrand;
  };

  ///		write the current state of the generator to a file
  bool read_rand_state(std::string rand_filename, MTRand &mtrand) {
    BP_Parse file;
    if(file.try_open(rand_filename)) {
      file.get_istream() >> mtrand;
      return true;
    }
    return false;
  };

}

#endif // BP_useful_CC
