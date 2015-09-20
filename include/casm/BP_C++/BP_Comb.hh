#ifndef BP_Comb_HH
#define BP_Comb_HH


#include "casm/BP_C++/BP_Vec.hh"
#include <cstring>
#include <sstream>
#include <iostream>

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
///  Coordinate classes

namespace BP {

  class BP_Comb {
  private:
    BP_Vec<bool> val_list;
    BP_Vec<unsigned int> pos_list;
    BP_Vec<unsigned int> fix_pos_list;
    BP_Vec<bool> fix_val_list;
    int curr_index;
    bool is_complete;
    int N;
    int K;
    unsigned long int count;

  public:
    BP_Comb(int n, int k) {
      reset(n, k);
    };

    void reset(int n, int k) {
      N = n;
      K = k;
      val_list = BP_Vec<bool>(N, 0);

      fix_pos_list.erase();
      fix_val_list.erase();

      pos_list.clear();
      for(int i = 0; i < K; i++) {
        pos_list.add(i);
        val_list[i] = 1;
      }

      curr_index = K - 1;

      is_complete = false;
      count = 0;
    };

    double total_combs() {
      double tot = 1;
      unsigned long int i;
      unsigned long int j = K;
      for(i = N; i > N - K; i--) {
        tot *= i;
      }

      for(i = K; i > 0; i--) {
        tot /= i;
      }

      return tot;

    };

    bool get(int i) const {
      return val_list[i];
    };

    bool operator[](int i) const {
      return val_list[i];
    };

    std::string get_bit_string() const {
      std::stringstream ss;
      for(int i = 0; i < size(); i++)
        ss << val_list[i];
      return ss.str();
    }

    BP_Vec<bool> get_all() const {
      return val_list;
    };

    void increment() {
      //cout << "begin increment()" << endl;

      bool ok = 1;
      do {
        //cout << " here b1" << endl;

        int curr_index = K - 1;
        if(pos_list[curr_index] == N - 1) {
          curr_index--;
          if(curr_index == -1) {
            is_complete = 1;
            //cout << "finish increment() 1a" << endl;
            return;
          }

          bool cont = 1;

          //cout << " here b2" << endl;
          //cout << " curr_index: " << curr_index << endl;

          do {
            if(pos_list[curr_index] < pos_list[curr_index + 1] - 1) {
              cont = 0;
            }
            else {
              curr_index--;
              if(curr_index == -1) {
                is_complete = 1;
                //cout << "finish increment() 1" << endl;
                return;
              }
            }

          }
          while(cont == 1);

          //cout << " curr_index: " << curr_index << endl;
          for(int i = curr_index; i < pos_list.size(); i++) {
            //cout << "i: " << i << endl;
            val_list[ pos_list[ i]] = 0;
          }
          //cout << " here b3" << endl;
          int new_start = pos_list[curr_index] + 1;
          //cout << " new_start: " << new_start << endl;
          for(int i = 0; i < pos_list.size() - curr_index; i++) {
            //cout << " here b3 i: "<< i << endl;
            pos_list[ curr_index + i] = new_start + i;
            //cout << " here b3 pos_list: "<< pos_list[ curr_index+ i] << endl;
            val_list[ pos_list[ curr_index + i]] = 1;
          }
          //cout << " here b4" << endl;

          //count++;
          curr_index = K - 1;
          //cout << " here b5" << endl;

        }
        else {
          //cout << " here a" << endl;
          val_list[ pos_list[ curr_index]] = 0;
          pos_list[curr_index]++;
          val_list[ pos_list[ curr_index]] = 1;
          //count++;
          //cout << " here a2" << endl;

        }

        //cout << " here c" << endl;

        ok = 1;
        for(int i = 0; i < fix_pos_list.size(); i++) {
          if(val_list[ fix_pos_list[i]] != fix_val_list[i]) {
            ok = 0;
            continue;
          }
        }

        //cout << " here d" << endl;

      }
      while(ok == 0);
      count++;

      //cout << "finish increment() 2" << endl;

    };

    int get_N() const {
      return N;
    };
    int get_K() const {
      return K;
    };
    int size() const {
      return N;
    };
    unsigned long int get_count() const {
      return count;
    };

    bool complete() const {
      return is_complete;
    };

    void fix(int i1, bool i2) {
      fix_pos_list.add(i1);
      fix_val_list.add(i2);
    };

    void unfix(int i1) {
      for(int i = 0; i < fix_pos_list.size(); i++) {
        if(fix_pos_list[i] == i1) {
          fix_pos_list.remove(i);
          fix_val_list.remove(i);
        }
      }
    };

    void unfix_all() {
      fix_pos_list.erase();
      fix_val_list.erase();
    };

    friend std::ostream &operator<<(std::ostream &outstream, const BP_Comb &c) {
      for(int i = 0; i < c.size(); i++)
        outstream << c[i];
      return outstream;
    };
  };


}

#endif // BP_Comb_HH
