#ifndef BP_RVG_HH
#define BP_RVG_HH

#include "casm/BP_C++/BP_GVec.hh"
#include "casm/BP_C++/BP_basic.hh"

namespace BP {

  /// \ingroup BP BP_RVG
  template< class T> class BP_RVG_base : public BP_Group<T> {
  public:
    virtual BP_GVec_Member<T> *add(BP_Gen_GVec_Member *, double) = 0;
    virtual BP_GVec_Member<T> *add_with_update(BP_Gen_GVec_Member *, double) = 0;
    virtual void remove(unsigned long int) = 0;
    virtual void remove(BP_Gen_GVec_Member *) = 0;
    virtual BP_GVec_Member<T> *pick(double i1) = 0;
    virtual double total() = 0 ;
    virtual double min_rate()const = 0;
    virtual double max_rate()const = 0;
    virtual BP_GVec_Member<T> *min_rate_member() = 0;
    virtual BP_GVec_Member<T> *max_rate_member() = 0;
    virtual void refresh_total() = 0;
    virtual unsigned long int size() const = 0;
    virtual void capacity(unsigned long int) = 0;
    virtual void print() = 0;
    virtual void clear() = 0;
    virtual void crazy_clear() = 0;
    virtual void free() = 0;
    virtual double get_rate(unsigned long int) = 0;

  public:
    BP_RVG_base() {
    };

    virtual ~BP_RVG_base() {
    };

  };


  /// \ingroup BP BP_RVG
  template <class T> class BP_RVG_linear : public BP_RVG_base<T> {
  private:

    BP_Vec<double>	rates;
    double			rates_total;
    bool			rates_total_ready;

  public:

    BP_RVG_linear();

    void capacity(unsigned long int i1);
    BP_GVec_Member<T> *add(BP_Gen_GVec_Member *i1, double i2);
    BP_GVec_Member<T> *add_with_update(BP_Gen_GVec_Member *i1, double i2);
    void remove(unsigned long int i1);
    void remove(BP_Gen_GVec_Member *i1);
    void clear();
    void crazy_clear();
    void free();
    BP_GVec_Member<T> *pick(double i1);   // pass a random number on range [0,1)  (excluding 1), for example: mtrand.randExc()
    double total();
    double curr_total() const;
    double min_rate() const;
    double max_rate() const;
    BP_GVec_Member<T> *min_rate_member();
    BP_GVec_Member<T> *max_rate_member();
    void refresh_total();
    unsigned long int size() const;
    void print();
    double get_rate(unsigned long int i1);

  };

  /// \ingroup BP BP_RVG
  template <class T> class BP_RVG_tree : public BP_RVG_base<T> {

  private:

    BP_Vec< BP_Vec<double> >	rate_sums;	// rate_sums[0][j] contains rate for event j, rate; capacity is always power of 2, and extra spots are set to 0.0
    // rate_sums[i>0][j] contains the sum: rate_sums[i-1][j*2] + rate_sums[i-1][j*2+1]

    bool rate_sums_ready;

    void increase_max_mvs();
    void add_new_rate(unsigned long int j, double r);
    void clear_to_zero();


  public:

    BP_RVG_tree();

    double get_rate(unsigned long int i1);
    void print();
    void capacity(unsigned long int i1);
    BP_GVec_Member<T> *add(BP_Gen_GVec_Member *i1, double i2);
    BP_GVec_Member<T> *add_with_update(BP_Gen_GVec_Member *i1, double i2);
    void remove(unsigned long int i1);
    void remove(BP_Gen_GVec_Member *i1);
    void clear();
    void crazy_clear();
    void free();
    BP_GVec_Member<T> *pick(double i1);
    double total();
    double curr_total() const;
    double min_rate() const;
    double max_rate() const;
    BP_GVec_Member<T> *min_rate_member();
    BP_GVec_Member<T> *max_rate_member();
    void refresh_total();
    unsigned long int size() const;

  };

  // to do: composition rejection algorithm
  /*template <class T> class BP_RVG_cr : public BP_RVG_base<T>
  {

  	public:
  	BP_RVG_cr()
  	{

  	};
  };
  */

  // /////////////////////////////////////
  // /////////////////////////////////////
  // Functions for (linear):
  //     template <class T> class BP_RVG_linear : public BP_RVG_base<T>

  //public:
  template<class T> BP_RVG_linear<T>::BP_RVG_linear() {
    rates_total = 0;
    rates_total_ready = true;
  };

  template<class T> void BP_RVG_linear<T>::capacity(unsigned long int i1) {
    BP_Gen_Group::capacity(i1);
    rates.capacity(i1);
  };

  template<class T> BP_GVec_Member<T> *BP_RVG_linear<T>::add(BP_Gen_GVec_Member *i1, double i2) {
    //std::cout << "begin linear::add()" << std::endl;
    rates_total += i2;
    rates.add(i2);

    return BP_Group<T>::add(i1);
  };

  template<class T> BP_GVec_Member<T> *BP_RVG_linear<T>::add_with_update(BP_Gen_GVec_Member *i1, double i2) {
    //std::cout << "begin linear::add()" << std::endl;
    rates_total += i2;
    rates.add(i2);

    return BP_Group<T>::add(i1);
  };

  template<class T> void BP_RVG_linear<T>::remove(unsigned long int i1) {
    //std::cout << "BP_RVG_linear remove" << std::endl;
    //if( i1 >= size())
    //{
    //	std::cout << "error: remove: " << i1 << " size: " << size() << std::endl;
    //	exit(1);
    //}

    //rates_total -= rates[i1];
    rates_total_ready = false;
    rates.remove(i1);



    BP_Gen_Group::remove(i1);

  };

  template<class T> double BP_RVG_linear<T>::curr_total() const {
    return sum(rates);
  };


  template<class T> void BP_RVG_linear<T>::remove(BP_Gen_GVec_Member *i1) {
    remove(BP_Gen_Group::get_index(i1));
  };

  template<class T> void BP_RVG_linear<T>::clear() {

    rates.clear();
    rates_total = 0;
    rates_total_ready = true;

    while(size() > 0)
      BP_Gen_Group::remove((unsigned long int) 0);

  };

  template<class T> void BP_RVG_linear<T>::crazy_clear() {
    rates.clear();
    rates_total = 0;
    rates_total_ready = true;

    BP_Gen_Group::crazy_clear();
  };

  template<class T> void BP_RVG_linear<T>::free() {
    rates.free();
    BP_Gen_Group::free();
  };

  // pass a random number on range [0,1)  (excluding 1)
  //   for example from mtrand.randExc();
  template<class T> BP_GVec_Member<T> *BP_RVG_linear<T>::pick(double i1) {
    if(!rates_total_ready)
      refresh_total();

    double rand = rates_total * i1;
    double sum = rates[0];
    unsigned long int index = 0;
    while(rand >= sum) {
      index++;
      sum += rates[index];
    }

    return BP_Group<T>::member(index);

  };

  template<class T> double BP_RVG_linear<T>::total() {
    if(!rates_total_ready) refresh_total();
    return rates_total;
  };

  template<class T> double BP_RVG_linear<T>::min_rate() const {
    return min(rates);
  };

  template<class T> double BP_RVG_linear<T>::max_rate() const {
    return max(rates);
  };

  template<class T> BP_GVec_Member<T> *BP_RVG_linear<T>::min_rate_member() {
    unsigned long int index;
    min(rates, index);
    return BP_Group<T>::member(index);
  };

  template<class T> BP_GVec_Member<T> *BP_RVG_linear<T>::max_rate_member() {
    unsigned long int index;
    max(rates, index);
    return BP_Group<T>::member(index);
  };


  template<class T> void BP_RVG_linear<T>::refresh_total() {
    rates_total = sum(rates);
    rates_total_ready = true;
  };

  template<class T> unsigned long int BP_RVG_linear<T>::size() const {
    return BP_Gen_Group::size();
  };

  template<class T> void BP_RVG_linear<T>::print() {
    std::cout << "Linear: " ;
    for(unsigned long int i = 0; i < size(); i++) {
      std::cout << rates[i] << " " ;
    }
    std::cout << std::endl;
  };

  template<class T> double BP_RVG_linear<T>::get_rate(unsigned long int i1) {
    return rates[i1];
  };


  // /////////////////////////////////////
  // /////////////////////////////////////
  // Functions for (tree):
  //     template <class T> class BP_RVG_tree : public BP_RVG_base<T>

  //private:
  template<class T> void  BP_RVG_tree<T>::increase_max_mvs() {
    // always powers of two
    //std::cout << "begin increase_max_mvs(): " << size() << std::endl;

    unsigned long int i, j;
    for(i = 0; i < rate_sums.size(); i++) {
      if(i == 0) {
        BP_Gen_Group::capacity(2 * rate_sums[i].size_capacity());
      }

      rate_sums[i].capacity(2 * rate_sums[i].size_capacity());
      while(rate_sums[i].size() < rate_sums[i].size_capacity())
        rate_sums[i].add(0.0);
    }
    rate_sums.add();
    rate_sums[ rate_sums.size() - 1].capacity(1);
    rate_sums[ rate_sums.size() - 1].add(rate_sums[ rate_sums.size() - 2][0]);



    //std::cout << "finish increase_max_mvs(): " << size() << std::endl;


    //print();



  };

  template<class T> void  BP_RVG_tree<T>::add_new_rate(unsigned long int j, double r) {
    //std::cout << "begin add_new_rate()" << std::endl;

    unsigned long int i = 0;
    rate_sums[i][j] = r;

    i++;
    j /= 2;


    while(i < rate_sums.size()) {
      //std::cout << "  i: " << i << " j: " << j << std::endl;
      rate_sums[i][j] = rate_sums[i - 1][2 * j] + rate_sums[i - 1][2 * j + 1];
      i++;
      j /= 2;
    }

    //std::cout << "finish add_new_rate()" << std::endl;

  };

  //public:
  template<class T> BP_RVG_tree<T>::BP_RVG_tree() {
    //std::cout << "begin BP_RVG_tree()" << std::endl;

    rate_sums.add();

    rate_sums[0].capacity(1);
    rate_sums[0].add(0.0);

    rate_sums_ready = true;

    //std::cout << "finish BP_RVG_tree()" << std::endl;

  };

  template<class T> double  BP_RVG_tree<T>::get_rate(unsigned long int i1) {
    return rate_sums[0][i1];
  };

  template<class T> void  BP_RVG_tree<T>::print() {
    std::cout << std::endl;
    std::cout << "  rate_sums.size(): " << rate_sums.size() << std::endl;
    for(int i = 0; i < rate_sums.size(); i++) {
      std::cout << "rate_sums[" << i << "]: " << rate_sums[i] << std::endl;

    }
    std::cout << std::endl;
    std::cout << std::endl;
  };



  template<class T> void BP_RVG_tree<T>::capacity(unsigned long int i1) {
    if(i1 > BP_Gen_Group::size()) {
      //BP_Gen_Group::capacity(i1);
      while(BP_Gen_Group::size_capacity() < i1)
        increase_max_mvs();
    }
    else if(i1 < BP_Gen_Group::size()) {
      unsigned long int j = 1;
      long int p = 0;
      while(j < BP_Gen_Group::size()) {
        j *= 2;
        p++;
      }

      BP_Gen_Group::capacity(ulint_pow(2, p));

      for(unsigned long int i = 0; i < rate_sums.size(); i++) {
        if(p >= 0) {
          rate_sums[i].capacity(ulint_pow(2, p));
          p--;
        }
        else {
          rate_sums[i].clear();
          rate_sums[i].free();
          rate_sums.remove(i);
        }
      }

      rate_sums.free();
    }

    return;
  };

  template<class T> BP_GVec_Member<T> *BP_RVG_tree<T>::add(BP_Gen_GVec_Member *i1, double i2) {
    //std::cout << "begin tree::add()" << std::endl;
    unsigned long int j = BP_Gen_Group::size();
    if(j == rate_sums[0].size_capacity()) {
      increase_max_mvs();
    }

    rate_sums[0][j] = i2;
    rate_sums_ready = false;

    //add_new_rate(j,i2);

    //std::cout << "finish tree::add()" << std::endl;

    return BP_Group<T>::add(i1);

  };

  template<class T> BP_GVec_Member<T> *BP_RVG_tree<T>::add_with_update(BP_Gen_GVec_Member *i1, double i2) {
    //std::cout << "begin tree::add_with_update()" << std::endl;

    unsigned long int j = BP_Gen_Group::size();
    if(j == rate_sums[0].size_capacity()) {
      increase_max_mvs();
    }

    add_new_rate(j, i2);

    //std::cout << "finish tree::add_with_update()" << std::endl;

    return BP_Group<T>::add(i1);

  };

  template<class T> void BP_RVG_tree<T>::remove(unsigned long int i1) {
    //std::cout << "begin tree::remove()" << std::endl;
    add_new_rate(i1, rate_sums[0][BP_Gen_Group::size() - 1]);
    add_new_rate(BP_Gen_Group::size() - 1, 0.0);
    BP_Gen_Group::remove(i1);
    //std::cout << "finish tree::remove()" << std::endl;
    return;
  };

  template<class T> void BP_RVG_tree<T>::remove(BP_Gen_GVec_Member *i1) {
    remove(BP_Gen_Group::get_index(i1));

    return;
  };

  template<class T> void BP_RVG_tree<T>::clear_to_zero() {
    unsigned long int i, j;
    unsigned long int Nj = BP_Gen_Group::size() - 1;

    /*for( i=0; i<rate_sums.size(); i++)
    {
    	for( j=0; j<Nj; j++)
    	{
    		rate_sums[i][j] = 0.0;

    	}
    	Nj /= 2;
    }*/
    //std::cout << "size: " <<  BP_Gen_Group::size() << std::endl;
    //std::cout << "rate_sums.size: " <<  rate_sums.size() << std::endl;
    //std::cout << "Nj: " << Nj << std::endl;

    for(i = 0; i < rate_sums.size(); i++) {
      for(j = Nj; j != -1; j--) {
        //std::cout << "i: " << i << " j: " << j << std::endl;
        rate_sums[i][j] = 0.0;

      }
      if(Nj != -1) Nj /= 2;
      //std::cout << "Nj: " << Nj << std::endl;

    }


  };

  template<class T> void BP_RVG_tree<T>::clear() {
    //while(size() > 0)
    //	remove((unsigned long int) 0);
    //
    //return;

    //std::cout << "begin tree::clear()" << std::endl;

    clear_to_zero();

    //std::cout << "rate_sums.size(): " << rate_sums.size() << std::endl;
    //std::cout << "size: " << size() << std::endl;
    //for( i=0; i<rate_sums.size(); i++)
    //{
    //	for( j=0; j<rate_sums[i].size(); j++)
    //	{
    //		if( rate_sums[i][j] != 0.0)
    //		{
    //			std::cout << "i: " << i << " j: " << j << std::endl;
    //		}
    //
    //	}
    //}


    //for( i=Nj-1; i!=-1; i--)
    //{
    //add_new_rate( i, 0.0);
    //BP_Gen_Group::remove(i);

    //std::cout << "i: " << i << std::endl;
    //remove(i);

    //add_new_rate( i-1, 0.0);
    //}

    while(size() > 0)
      BP_Gen_Group::remove((unsigned long int) 0);

    rate_sums_ready = true;

    //std::cout << "finish tree::clear() size: " << size() << std::endl;

  };

  template<class T> void BP_RVG_tree<T>::crazy_clear() {
    //std::cout << "begin tree::crazy_clear()" << std::endl;
    clear_to_zero();
    //std::cout << "  tree::crazy_clear() here 1" << std::endl;
    BP_Gen_Group::crazy_clear();
    rate_sums_ready = true;
    //std::cout << "finish tree::crazy_clear()" << std::endl;

  };

  template<class T> void BP_RVG_tree<T>::free() {
    capacity(BP_Gen_Group::size());
  };

  template<class T> BP_GVec_Member<T> *BP_RVG_tree<T>::pick(double i1) {
    //std::cout << "begin tree::pick()" << std::endl;
    if(!rate_sums_ready) {
      //std::cout << "  refresh" << std::endl;
      refresh_total();
    }

    unsigned long int i, j;

    i = rate_sums.size() - 1;
    j = 0;

    //std::cout << " init i: " << i << " j: " << j << std::endl;
    //std::cout << " pick size: " << size() << std::endl;
    //std::cout << " pick rate_sums.size: " << rate_sums.size() << std::endl;
    //std::cout << " pick total: " << total() << std::endl;
    //std::cout << " pick rate_sums[i][j]: " << rate_sums[i][j] << std::endl;
    //std::cout << " pick rand: " << i1 << std::endl;


    double rand = rate_sums[i][j] * i1;
    //std::cout << " pick rand*rate_sum: " << rand << std::endl;

    if(i > 0) i--;

    //std::cout << " start i: " << i << " j: " << j << std::endl;
    while(i > 0) {
      if(rand >= rate_sums[i][j]) {
        rand -= rate_sums[i][j];
        j = (j + 1) * 2;
      }
      else {
        j *= 2;
      }
      i--;

      //std::cout << "i: " << i << " j: " << j << std::endl;
    }

    //std::cout << " final i: " << i << " j: " << j << std::endl;

    if(rand >= rate_sums[0][j]) {
      //std::cout << " pick r1: " << j+1 << std::endl;
      return BP_Group<T>::member(j + 1);
    }
    else {
      //std::cout << " pick r2: " << j << std::endl;
      return BP_Group<T>::member(j);
    }
  };

  template<class T> double BP_RVG_tree<T>::total() {
    if(!rate_sums_ready) refresh_total();
    return rate_sums[ rate_sums.size() - 1][0];
  };

  template<class T> double BP_RVG_tree<T>::curr_total() const {
    return sum(rate_sums[0]);
  };

  template<class T> double BP_RVG_tree<T>::min_rate() const {
    return min(rate_sums[0]);
  };

  template<class T> double BP_RVG_tree<T>::max_rate() const {
    return max(rate_sums[0]);
  };

  template<class T> BP_GVec_Member<T> *BP_RVG_tree<T>::min_rate_member() {
    unsigned long int index;
    min(rate_sums[0], index);
    return BP_Group<T>::member(index);
  };

  template<class T> BP_GVec_Member<T> *BP_RVG_tree<T>::max_rate_member() {
    unsigned long int index;
    max(rate_sums[0], index);
    return BP_Group<T>::member(index);
  };

  template<class T> void BP_RVG_tree<T>::refresh_total() {
    //std::cout << "begin tree::refresh_total()" << std::endl;
    if(BP_Gen_Group::size() == 0) return;

    unsigned long int i, j;
    unsigned long int Nj = (BP_Gen_Group::size() - 1) / 2;

    //std::cout << "BP_Gen_Group::size(): " << BP_Gen_Group::size() << std::endl;
    //std::cout << "Nj: " << Nj << std::endl;

    for(i = 1; i < rate_sums.size(); i++) {
      //std::cout << "i: " << i << " rate_sums[i].size(): " << rate_sums[i].size() << std::endl;
      for(j = 0; j <= Nj; j++) {
        //std::cout << "i: " << i << " j: " << j << std::endl;
        rate_sums[i][j] = rate_sums[i - 1][2 * j] + rate_sums[i - 1][2 * j + 1];

      }
      Nj /= 2;
      //std::cout << "Nj: " << Nj << std::endl;

    }


    rate_sums_ready = true;
    //std::cout << "finish tree::refresh_total()" << std::endl;

  };

  template<class T> unsigned long int BP_RVG_tree<T>::size() const {
    return BP_Gen_Group::size();
  };


  // /////////////////////////////////////
  // /////////////////////////////////////
  // to do: Functions for (cr = composition rejection):
  //     composition rejection algorithm
  //     template <class T> class BP_RVG_cr : public BP_RVG_base<T>


}

#endif // BP_RVG_HH

