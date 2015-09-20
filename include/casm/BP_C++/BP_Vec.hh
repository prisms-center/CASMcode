#ifndef BP_Vec_HH
#define BP_Vec_HH

#include <iostream>
#include <cmath>
#include <cstdlib>
#include "casm/BP_C++/BP_basic.hh"

namespace BP {


  /*
  // A template class for a Vector-like container:
  //    0 1 2 ...
  //   |---------------N--------------------------------N_alloc----------------------N_max|
  //   |  active       |  allocated but inactive           |         unallocated          |
  //
  //
  //	 T must have copy constructor T::T( T t){};  (so that "= new T(t) works)
  //	 T must be assignable: T=T;
  //
  //	 Contains an array of pointers to T objects, so all operations actually done moving pointers around instead of objects
  //				- objects do not move
  //				- pointers to objects always remain valid, until the object is 'freed'
  //
  //		access obect "BP_Vec<T> V" with V[i]
  //
  //		V1 = V2, assigns V2 to V1		// V1 now points to copies of the V2 objects, not the V2 objects themselves
  //
  //		'void V.add(U u)' adds an object
  //			if allocated space is available
  //				if ( U == T) uses 'V[N] = u'				(using overloaded void V.add(T t))
  //				else  uses 'V[N] = T(u)'
  //			else
  //				'V[N] = new T(u)'
  //
  //
  //		'void V.add_in_place(U u, int i)' adds an object 'u' at position 'i', maintaining order
  //
  //		'void V.remove( int i)' removes object at [i] and puts object at [N-1] at [i], and reduces N = N-1
  //
  //		'void V.ordered_remove( int i)' removes object at [i], shifts all objects at [>i] down 1 place, and reduces N = N-1
  //
  //		'void V.swap( int i, int j)' swaps objects at [i] and [j]
  //
  //		'int V.size()' returns N;
  //
  //		'void V.clear()' sets N = 0
  //
  //		'void V.capacity( int i)' sets N_max = i (and N and N_alloc if necessary)
  //
  //
  //	Memory is only deallocated by ~V() or V.free() ( or "=" or V.capacity() if it becomes smaller)
  //
  //		'void V.free()' calls 'delete V[i]' for all N <= i < N_alloc
  //
  //
  //	!! use functions V.put_in() and V.take_out() with caution: improper use could cause memory leakage or double deletion
  //
  //			'void V.put_in(T* t)' makes V[N] = t, and increments N = N+1... !!! it does not make a copy of *T !!!
  //
  //			'void V.put_in_place(T* t, int i)' puts t in position i, maintaining order !!! it does not make a copy of *T !!!
  //
  //			'T* V.take_out(int i)' returns V[i], sets V[i] = V[N-1], sets V[N-1] = V[N_alloc-1],
  //				and decrements N = N-1 and N_alloc = N_alloc-1 !!! it does not make a copy of *T !!!
  //
  //			'T* V.ordered_take_out(int i)' returns V[i], shifts all objects at [>i] down 1 space,
  //				and decrements N = N-1 and N_alloc = N_alloc-1 !!! it does not make a copy of *T !!!
  */

  /// \ingroup BP BP_Vec
  template<class T> class BP_Vec {
  protected:

    unsigned long int N_max;			//lenth of T **val
    unsigned long int N_alloc;		//max i, for which val[i] = new T has been called
    unsigned long int N;				//current max i, (val[i], for N <= i < N_alloc, point to existing objects, but are not currently in use)
    int cap_incr;				//amount that capacity is increased if N == N_alloc == N_max when add() or put_in() is called
    T **val;

  public:
    BP_Vec &operator=(const BP_Vec &);
    inline T &operator[](unsigned long int i1) {
      return *val[i1];
    };
    inline const T &operator[](unsigned long int i1) const {
      return *val[i1];
    };
    inline T &last() {
      return *val[ N - 1];
    };
    inline const T &last() const {
      return *val[ N - 1];
    };
    void capacity(unsigned long int);
    void min_capacity() {
      capacity(N);
    };
    void cap_increment(unsigned long int);
    T *add();
    T *add(const T &);
    template <class U> T *add(const U &);
    template <class U> T *add_in_place(const U &, unsigned long int);
    void remove(unsigned long int);
    void ordered_remove(unsigned long int);
    void swap(unsigned long int, unsigned long int);
    void ordered_swap(unsigned long int i1, unsigned long int i2);
    void clear();
    void free();
    void erase();
    inline unsigned long int size() const {
      return N;
    };
    inline unsigned long int size_alloc() const {
      return N_alloc;
    };
    inline unsigned long int size_capacity() const {
      return N_max;
    };
    void put_in(T *);
    void put_in_place(T *, unsigned long int);
    T *take_out(unsigned long int);
    T *ordered_take_out(unsigned long int);
    BP_Vec range(unsigned long int, unsigned long int);
    void append(const BP_Vec &);
    //unsigned long int find_first( const T&);
    bool contains(const T &i1) const;
    bool find_first(const T &i1, unsigned long int &index) const;
    BP_Vec<unsigned long int> find_all(const T &) const;

    void sort();
    template<class U> void sort(BP_Vec<U> &values);
    void flip();

    template<class U> BP_Vec &operator+=(const U &);
    template<class U> BP_Vec &operator-=(const U &);
    template<class U> BP_Vec &operator*=(const U &);
    template<class U> BP_Vec &operator/=(const U &);

    template<class U> BP_Vec &operator+=(const BP_Vec<U> &);
    template<class U> BP_Vec &operator-=(const BP_Vec<U> &);
    template<class U> BP_Vec &operator*=(const BP_Vec<U> &);
    template<class U> BP_Vec &operator/=(const BP_Vec<U> &);

    template<class U> friend std::ostream &operator<<(std::ostream &outstream, const BP_Vec<U> &vec);
    void write_multiline(std::ostream &sout);
    void write_inline(std::ostream &sout);

    void init();
    BP_Vec();
    BP_Vec(const BP_Vec &);
    BP_Vec(unsigned long int i1, const T &i2);
    ~BP_Vec();

  private:
    void qsort();
    template<class U> void qsort(BP_Vec<U> &values);
    void isort();
    template<class U> void isort(BP_Vec<U> &values);
    void hsort();

    void heap(unsigned long int N);
    void sift_down(unsigned long int start, unsigned long int end);

  };



  //  template<typename T>
  //  jsonParser &to_json(const BP_Vec<T> &value, jsonParser &json) {
  //    json.put_array();
  //    for(int i = 0; i < value.size(); i++)
  //      json.push_back(value[i]);
  //    return json;
  //  }
  //
  //  template<typename T>
  //  void from_json(BP_Vec<T> &value, const jsonParser &json) {
  //    try {
  //      value.clear();
  //      T jvalue;
  //      for(int i = 0; i < json.size(); i++)
  //        from_json(jvalue, json[i]);
  //        value.add(jvalue);
  //    }
  //    catch(...) {
  //      /// re-throw exceptions
  //      throw;
  //    }
  //  }

  /*
  template<class T> unsigned long int BP_Vec<T>::find_first( const T &i1)
  {
  	for( unsigned long int i=0; i<size(); i++)
  		if( *val[i] == i1)
  			return i;
  	return -1;
  };
  */

  template<class T> bool BP_Vec<T>::contains(const T &i1) const {
    for(unsigned long int i = 0; i < size(); i++)
      if(*val[i] == i1) {
        return true;
      }
    return false;
  };

  template<class T> bool BP_Vec<T>::find_first(const T &i1, unsigned long int &index) const {
    for(unsigned long int i = 0; i < size(); i++)
      if(*val[i] == i1) {
        index = i;
        return true;
      }
    return false;
  };

  template<class T> BP_Vec<unsigned long int> BP_Vec<T>::find_all(const T &i1) const {
    BP_Vec<unsigned long int> list;
    for(unsigned long int i = 0; i < size(); i++)
      if(*val[i] == i1)
        list.add(i);
    return list;
  };


  template<class U>  std::ostream &operator<<(std::ostream &sout, const BP_Vec<U> &vec) {
    for(unsigned long int i = 0; i < vec.size(); i++)
      sout << "'" << vec[i] << "' ";
    return sout;
  };

  template<class T> void BP_Vec<T>::write_multiline(std::ostream &sout) {
    for(unsigned long int i = 0; i < size(); i++)
      sout << i << ": " << *val[i] << std::endl;
  };

  template<class T> void BP_Vec<T>::write_inline(std::ostream &sout) {
    for(unsigned long int i = 0; i < size(); i++)
      sout << *val[i] << " ";
  };


  // exact assignment, including the allocated but unused (N <= i < N_alloc)
  template<class T> BP_Vec<T> &BP_Vec<T>::operator=(const BP_Vec &i1) {
    unsigned long int i;
    cap_incr = i1.cap_incr;
    capacity(i1.N_max);

    if(i1.N_alloc > N_alloc) {
      // Objects exist, just assign from i1.val
      for(i = 0; i < N_alloc; i++)
        *(val[i]) = *(i1.val[i]);

      // Objects must be created from copy of i1.val
      for(i = N_alloc; i < i1.N_alloc; i++) {
        //cout << "new 1" << endl;
        val[i] = new T(*(i1.val[i]));
        //*(val[i]) = *(i1.val[i]);
      }

      N = i1.N;
      N_alloc = i1.N_alloc;
    }
    else {
      // Objects exist, just assign from i1.val
      for(i = 0; i < i1.N_alloc; i++)
        *(val[i]) = *(i1.val[i]);

      for(i = i1.N_alloc; i < N_alloc; i++)
        delete val[i];

      N = i1.N;
      N_alloc = i1.N_alloc;

    }

    return *this;
  };

  // add/subtract/multiply/divide a scalar into all elements
  template<class T> template<class U>  BP_Vec<T>  &BP_Vec<T>::operator+=(const U &B) {
    for(unsigned long int i = 0; i < size(); i++)
      (*this)[i] += B;

    return *this;
  };
  template<class T> template<class U>  BP_Vec<T>  &BP_Vec<T>::operator-=(const U &B) {
    for(unsigned long int i = 0; i < size(); i++)
      (*this)[i] -= B;

    return *this;
  };

  template<class T> template<class U>  BP_Vec<T>  &BP_Vec<T>::operator*=(const U &B) {
    for(unsigned long int i = 0; i < size(); i++)
      (*this)[i] *= B;

    return *this;
  };
  template<class T> template<class U>  BP_Vec<T>  &BP_Vec<T>::operator/=(const U &B) {
    for(unsigned long int i = 0; i < size(); i++)
      (*this)[i] /= B;

    return *this;
  };

  // elementwise vector-vector operations
  template<class T> template<class U>  BP_Vec<T>  &BP_Vec<T>::operator+=(const BP_Vec<U> &B) {
    for(unsigned long int i = 0; i < size(); i++)
      (*this)[i] += B[i];

    return *this;
  };
  template<class T> template<class U>  BP_Vec<T>  &BP_Vec<T>::operator-=(const BP_Vec<U> &B) {
    for(unsigned long int i = 0; i < size(); i++)
      (*this)[i] -= B[i];

    return *this;
  };

  template<class T> template<class U>  BP_Vec<T>  &BP_Vec<T>::operator*=(const BP_Vec<U> &B) {
    for(unsigned long int i = 0; i < size(); i++)
      (*this)[i] *= B[i];

    return *this;
  };
  template<class T> template<class U>  BP_Vec<T>  &BP_Vec<T>::operator/=(const BP_Vec<U> &B) {
    for(unsigned long int i = 0; i < size(); i++)
      (*this)[i] /= B[i];

    return *this;
  };


  // Change N_max to i1
  //  uses one "new T*[i1]", no "new T"
  //  calls "delete val[i]" for i>=i1
  //  val[i], where i<i1 remain unchanged
  template<class T> void BP_Vec<T>::capacity(unsigned long int i1) {
    if(i1 == N_max) return;

    unsigned long int i;
    T **tmp;
    //cout << "new 2: " << i1 << endl;
    tmp = new T*[i1];

    if(i1 > N_alloc) {
      for(i = 0; i < N_alloc; i++)
        tmp[i] = val[i];
    }
    else {
      for(i = 0; i < i1; i++)
        tmp[i] = val[i];

      for(i = i1; i < N_alloc; i++)
        delete val[i];

      if(i1 < N) N = i1;
      N_alloc = i1;
    }

    if(N_max != 0) delete [] val;
    val = tmp;
    N_max = i1;

  };

  template<class T> void BP_Vec<T>::cap_increment(unsigned long int i1) {
    cap_incr = i1;
    if(cap_incr == 0) cap_incr = 1;
  };

  // Put object pointed to into BP_Vec without making a copy
  // could be dangerous... results in two pointers to same object, and now only this one can be used to delete
  template<class T> void BP_Vec<T>::put_in(T *i1) {
    if(N == N_max) {
      capacity(N_max + cap_incr);
    }

    if(N == N_alloc) {
      val[N] = i1;
      N++;
      N_alloc++;
    }
    else {
      if(N_alloc == N_max) {
        delete val[N];
        val[N] = i1;
        N++;
      }
      else {
        val[N_alloc] = i1;
        swap(N, N_alloc);
        N_alloc++;
        N++;
      }
    }

    return;
  };

  // Put object pointed to into BP_Vec at position i2 without making a copy
  // could be dangerous... results in two pointers to same object, and now only this one can be used to delete
  template<class T> void BP_Vec<T>::put_in_place(T *i1, unsigned long int i2) {
    put_in(i1);

    ordered_swap(N - 1, i2);

    return;
  };


  // Return pointer to object at val[i1], eliminate pointer
  // could be dangerous... removes pointer without freeing memory, which must now be done elsewhere
  template<class T> T *BP_Vec<T>::take_out(unsigned long int i1) {
    T *t = val[i1];
    val[i1] = val[N - 1];
    val[N - 1] = val[N_alloc - 1];
    N--;
    N_alloc--;

    return t;
  };

  // Return pointer to object at val[i1], eliminate pointer, shift all down one space
  // could be dangerous... removes pointer without freeing memory, which must now be done elsewhere
  template<class T> T *BP_Vec<T>::ordered_take_out(unsigned long int i1) {
    T *t = val[i1];
    for(unsigned long int i = i1; i < N - 1; i++)
      val[i] = val[i + 1];
    val[N - 1] = val[N_alloc - 1];
    N--;
    N_alloc--;

    return t;
  };

  // add an object (as is if memory already allocated, or initialized with T() if memory needs to be allocated) and return a pointer to it
  template<class T> T *BP_Vec<T>::add() {
    if(N == N_max) {
      capacity(N_max + cap_incr);
    }

    if(N == N_alloc) {
      //cout << "new 3" << endl;
      val[N] = new T;
      N++;
      N_alloc++;
    }
    else {
      N++;
    }

    return val[N - 1];
  };

  // add to BP_Vec (like push_back())
  template<class T> T *BP_Vec<T>::add(const T &i1) {
    if(N == N_max) {
      capacity(N_max + cap_incr);
    }

    if(N == N_alloc) {
      //cout << "new 4" << endl;

      val[N] = new T(i1);
      N++;
      N_alloc++;
    }
    else {
      *(val[N]) = i1;
      N++;
    }

    return val[N - 1];
  };

  // add to BP_Vec (like push_back())
  template <class T> template <class U> T *BP_Vec<T>::add(const U &i1) {
    if(N == N_max) {
      capacity(N_max + cap_incr);
    }

    if(N == N_alloc) {
      val[N] = new T(i1);
      N++;
      N_alloc++;
    }
    else {
      *(val[N]) = T(i1);
      N++;
    }

    return val[N - 1];
  };


  // add to BP_Vec in position i2
  template <class T> template <class U> T *BP_Vec<T>::add_in_place(const U &i1, unsigned long int i2) {
    add(i1);

    ordered_swap(N - 1, i2);

    return val[i2];
  };

  template<class T> void BP_Vec<T>::ordered_swap(unsigned long int i1, unsigned long int i2) {
    if(i1 < i2) {
      T *tmp = val[i1];
      for(unsigned long int i = i1; i < i2; i++)
        val[i] = val[i + 1];
      val[i2] = tmp;
    }
    else {
      T *tmp = val[i1];
      for(unsigned long int i = i1; i > i2; i--)
        val[i] = val[i - 1];
      val[i2] = tmp;
    }

  };

  // Removes object val[i1] from list, moves val[N-1] to val[i1].  Does not actually free memory.
  template<class T> void BP_Vec<T>::remove(unsigned long int i1) {
    swap(i1, N - 1);
    N--;

  };

  // Removes object val[i1] from list, moves val[i1+1] to val[N-1] down one index.  Does not actually free memory.
  template<class T> void BP_Vec<T>::ordered_remove(unsigned long int i1) {
    T *tmp = val[i1];
    for(unsigned long int i = i1; i < N - 1; i++)
      val[i] = val[i + 1];
    val[N - 1] = tmp;
    N--;

  };

  // Deletes objects stored in val[i] i>=N.  Frees unused memory.
  template<class T> void BP_Vec<T>::free() {
    for(unsigned long int i = N; i < N_alloc; i++) {
      delete val[i];

    }

    N_alloc = N;

  };

  // Swap val[a] and val[b]
  template<class T> void BP_Vec<T>::swap(unsigned long int a, unsigned long int b) {
    T *tmp = val[a];
    val[a] = val[b];
    val[b] = tmp;
  };

  // Set N=0. Does not actually clear memory.
  template<class T> void BP_Vec<T>::clear() {
    N = 0;
  };

  template<class T> void BP_Vec<T>::erase() {
    clear();
    free();

  };


  template<class T> BP_Vec<T> BP_Vec<T>::range(unsigned long int i1, unsigned long int i2) {
    BP_Vec<T> subvec;
    subvec.capacity(i2 - i1 + 1);
    for(unsigned long int i = i1; i <= i2; i++)
      subvec.add(*val[i]);
    return subvec;
  };

  template<class T> void BP_Vec<T>::append(const BP_Vec<T> &v2) {
    if(size() + v2.size() < N_max)
      capacity(size() + v2.size());

    for(unsigned long int i = 0; i < v2.size(); i++)
      add(v2[i]);
  };

  template<class T> void BP_Vec<T>::flip() {
    // flip from '1, 2, 3' to '3, 2, 1'

    unsigned long int i, j;
    i = 0;
    j = size() - 1;
    while(i < j) {
      swap(i, j);
      i++;
      j--;
    }
  };

  template<class T> void BP_Vec<T>::sort() {
    if(size() < 50)
      isort();
    else
      qsort();
  }

  template<class T> template<class U> void BP_Vec<T>::sort(BP_Vec<U> &values) {
    if(size() < 50)
      isort(values);
    else
      qsort(values);
  }

  template<class T> void BP_Vec<T>::qsort() {
    /// quicksort sorting algorithm (in ascending order)
    ///  - assumes that T::operator<(T) exists


    unsigned long int left, right, middle, pivot, curr, i;
    BP_Vec<unsigned long int> queue_left, queue_right;

    // estimate queue size
    unsigned long int cap = (unsigned long int) 3 * (log(size()) / log(2));
    //cout << "cap: " << cap << endl;
    queue_left.capacity(cap);
    queue_right.capacity(cap);
    queue_left.cap_increment(cap);
    queue_right.cap_increment(cap);

    queue_left.add(0);
    queue_right.add(size() - 1);

    unsigned long int max_q = 0;


    do {
      left = queue_left[0];
      right = queue_right[0];
      queue_left.remove(0);
      queue_right.remove(0);

      if(left < right) {

        // avoid using (right+left) since this could be larger than max allowed value
        middle = left + (right - left) / 2;

        // choose median of (left, middle, right) for pivot
        if(*val[left] < *val[middle]) {
          if(*val[middle] < *val[right]) pivot = middle;
          else {
            if(*val[left] < *val[right]) pivot = right;
            else pivot = left;
          }
        }
        else {
          if(*val[right] < *val[middle]) pivot = middle;
          else {
            if(*val[right] < *val[left]) pivot = right;
            else pivot = left;
          }
        }


        // move pivot to end
        swap(pivot, right);

        // in the current region (left:right-1), seperate the values less than the pivot value (which is now at 'right')
        curr = left;
        for(i = left; i < right; i++) {
          if(*val[i] < *val[right]) {
            swap(curr, i);
            curr++;
          }
        }

        // move the pivot value to the 'curr' position
        swap(curr, right);

        // now everything in 'left:curr-1' is < *val[curr], and everything in 'curr+1:right' is >= *val[curr]
        //   so add those two regions to the queue

        if(curr != 0) {
          if(left != curr - 1) {
            //cout << "  add region: " << left << ":" << curr-1 << endl;
            queue_left.add(left);
            queue_right.add(curr - 1);
          }
        }

        if(curr + 1 != right) {
          //cout << "  add region: " << curr+1 << ":" << right << endl;
          queue_left.add(curr + 1);
          queue_right.add(right);
        }
        // repeat until there is nothing left
      }

      if(queue_left.size() > max_q)
        max_q = queue_left.size();


    }
    while(queue_left.size() > 0);

    //cout << "max_q: " << max_q << endl;
    //cout << "max_q/cap: " << (1.0*max_q)/(1.0*cap) << endl;
  };

  template<class T> template<class U> void BP_Vec<T>::qsort(BP_Vec<U> &values) {
    /// quicksort sorting algorithm (in ascending order)
    ///   - sort *this by the values in argument 'values'
    ///   - 'values' also becomes sorted
    ///   - assumes that U::operator<(U) exists

    if(size() != values.size()) {
      //cout << "Error in qsort(values):  *this.size() != values.size()" << endl;
      return;
    }

    //cout << "values.size(): " << values.size() << endl;
    //cout << "qsort: " << values << endl;

    unsigned long int left, right, middle, pivot, curr, i;
    BP_Vec<unsigned long int> queue_left, queue_right;

    // estimate queue size
    unsigned long int cap = (unsigned long int) ceil(3 * (log(size()) / log(2)));
    //cout << "cap: " << cap << endl;
    queue_left.capacity(cap);
    queue_right.capacity(cap);
    queue_left.cap_increment(cap);
    queue_right.cap_increment(cap);

    queue_left.add(0);
    queue_right.add(size() - 1);

    //unsigned long int count = 0;


    do {
      //cout << "  queue_left.size(): " << queue_left.size() << endl;
      left = queue_left[0];
      right = queue_right[0];
      queue_left.remove(0);
      queue_right.remove(0);

      //cout << "count: " << count << "  left: " << left << "  right: " << right << endl;
      //count++;

      if(left < right) {
        // avoid using (right+left) since this could be larger than max allowed value
        middle = left + (right - left) / 2;

        //cout << "  left: " << left << "  middle: " << middle << "  right: " << right << endl;
        //cout << "          vals: " << values[left] << "  " << values[middle] << "  " << values[right] << endl;

        // choose median of (left, middle, right) for pivot
        if(values[left] < values[middle]) {
          if(values[middle] < values[right]) pivot = middle;
          else {
            if(values[left] < values[right]) pivot = right;
            else pivot = left;
          }
        }
        else {
          if(values[right] < values[middle]) pivot = middle;
          else {
            if(values[right] < values[left]) pivot = right;
            else pivot = left;
          }
        }

        //cout << "  pivot: " << pivot << "  val: " << values[pivot] << endl;

        // move pivot to end
        swap(pivot, right);
        values.swap(pivot, right);

        // in the current region (left:right-1), seperate the values less than the pivot value (which is now at 'right')
        curr = left;
        for(i = left; i < right; i++) {
          //cout << "  i: " << i << "  curr: " << curr << endl;

          if(values[i] < values[right]) {
            swap(curr, i);
            values.swap(curr, i);
            curr++;
            //cout << "    swap" << endl;
          }

        }

        //cout << "  move pivot" << endl;
        // move the pivot value to the 'curr' position
        swap(curr, right);
        values.swap(curr, right);

        // now everything in 'left:curr-1' is < values[curr], and everything in 'curr+1:right' is >= values[curr]
        //   so add those two regions to the queue

        if(curr != 0) {
          if(left != curr - 1) {
            //cout << "  add region: " << left << ":" << curr-1 << endl;
            queue_left.add(left);
            queue_right.add(curr - 1);
          }
        }

        if(curr + 1 != right) {
          //cout << "  add region: " << curr+1 << ":" << right << endl;
          queue_left.add(curr + 1);
          queue_right.add(right);

        }
        // repeat until there is nothing left

      }


    }
    while(queue_left.size() > 0);

  };

  template<class T> void BP_Vec<T>::isort() {
    /// insertion sorting algorithm (in ascending order)
    ///  - assumes that T::operator<(T) exists


    unsigned long int i, j;

    for(i = 1; i < size(); i++) {
      j = i;
      while(*val[i] < *val[j - 1]) {
        j--;
        if(j == 0) break;
      }
      ordered_swap(i, j);
    }
  };

  template<class T> template<class U> void BP_Vec<T>::isort(BP_Vec<U> &values) {
    /// insertion sorting algorithm (in ascending order)
    ///   - sort *this by the values in argument 'values'
    ///   - 'values' also becomes sorted
    ///   - assumes that U::operator<(U) exists

    if(size() != values.size()) {
      //cout << "Error in isort(values):  *this.size() != values.size()" << endl;
      return;
    }


    unsigned long int i, j;

    for(i = 1; i < size(); i++) {
      j = i;
      while(values[i] < values[j - 1]) {
        j--;
        if(j == 0) break;
      }
      ordered_swap(i, j);
      values.ordered_swap(i, j);
    }
  };

  template<class T> void BP_Vec<T>::hsort() {
    /// heapsort sorting algorithm (in ascending order)
    ///  - assumes that T::operator<(T) exists

    heap(size());

    unsigned long int end = size() - 1;

    while(end > 0) {
      //cout << "hsort: " << end << endl;
      swap(0, end);
      end--;
      sift_down(0, end);
    }


  };

  template<class T> void BP_Vec<T>::heap(unsigned long int N) {
    /// makes heap for hsort()

    unsigned long int start = (N - 2) / 2;
    bool stop = false;

    while(!stop) {
      //cout << "heap start: " << start << endl;
      sift_down(start, N - 1);
      if(start == 0) stop = true;
      start--;
    }
  };

  template<class T> void BP_Vec<T>::sift_down(unsigned long int start, unsigned long int end) {
    /// gets the heap in order for heap()
    unsigned long int root, child, i;

    root = start;

    while(root * 2 + 1 <= end) {
      //cout << "root: " << root << endl;
      child = root * 2 + 1;
      i = root;

      if(*val[i] < *val[child])
        i = child;

      if((child + 1 <= end) && (*val[i] < *val[child + 1]))
        i = child + 1;

      if(i != root) {
        swap(root, i);
        root = i;
      }
      else
        return;
    }

  };


  ///////////////////////////////////
  /// Constructors & Destructor

  template<class T> void BP_Vec<T>::init() {
    N_max = 0;
    N_alloc = 0;
    N = 0;
    cap_incr = 10;
  };

  template<class T> BP_Vec<T>::BP_Vec() {
    init();

  };

  template<class T> BP_Vec<T>::BP_Vec(const BP_Vec &i1) {
    init();
    *this = i1;

  };

  template<class T> BP_Vec<T>::BP_Vec(unsigned long int i1, const T &i2) {
    init();
    capacity(i1);
    for(unsigned long int i = 0; i < i1; i++)
      add(i2);
  };

  template<class T> BP_Vec<T>::~BP_Vec() {
    for(unsigned long int i = 0; i < N_alloc; i++) {
      delete val[i];
    }

    if(N_max != 0) delete [] val;

  };


  //////////////////////////////////////
  /// Bonus functions

  bool is_integer(const BP_Vec< double > &m, double tol);

  bool is_integer(const BP_Vec< BP_Vec<double> > &m, double tol);

  double mean(const BP_Vec<double> &i_list);

  double rms(const BP_Vec<double> &i_list);


  //////////////////////////////////////
  /// Bonus template functions

  template<class T> T sum(const BP_Vec<T> &i_list) {
    T s = 0;
    for(unsigned long int i = 0; i < i_list.size(); i++)
      s += i_list[i];
    return s;
  }

  template<class T> T min(const BP_Vec<T> &i_list) {
    unsigned long int min_index = 0;
    for(unsigned long int i = 1; i < i_list.size(); i++)
      if(i_list[i] < i_list[min_index])
        min_index = i;

    return i_list[min_index];
  };

  template<class T> T min(const BP_Vec<T> &i_list, unsigned long int &min_index) {
    min_index = 0;
    for(unsigned long int i = 1; i < i_list.size(); i++)
      if(i_list[i] < i_list[min_index]) {
        min_index = i;
      }

    return i_list[min_index];
  };

  template<class T> T min(const BP_Vec<T> &i_list, BP_Vec<unsigned long int> &min_index) {
    min_index.clear();
    min_index.add(0);
    for(unsigned long int i = 1; i < i_list.size(); i++)
      if(i_list[i] < i_list[min_index[0]]) {
        min_index.clear();
        min_index.add(i);
      }
      else if(i_list[i] == min) {
        min_index.add(i);
      }

    return i_list[min_index[0]];
  };

  template<class T> unsigned long int min_index(const BP_Vec<T> &i_list) {
    unsigned long int min_index = 0;
    for(unsigned long int i = 1; i < i_list.size(); i++)
      if(i_list[i] < i_list[min_index]) {
        min_index = i;
      }

    return min_index;
  };

  template<class T> BP_Vec<unsigned long int> min_indices(const BP_Vec<T> &i_list) {
    BP_Vec<unsigned long int> min_index;
    min_index.add(0);
    for(unsigned long int i = 1; i < i_list.size(); i++)
      if(i_list[i] < i_list[min_index[0]]) {
        min_index.clear();
        min_index.add(i);
      }
      else if(i_list[i] == min) {
        min_index.add(i);
      }

    return min_index;
  };

  template<class T> T max(const BP_Vec<T> &i_list) {
    unsigned long int max_index = 0;
    for(int i = 1; i < i_list.size(); i++)
      if(i_list[i] > i_list[max_index])
        max_index = i;

    return i_list[max_index];
  };

  template<class T> T max(const BP_Vec<T> &i_list, unsigned long int &max_index) {
    max_index = 0;
    for(unsigned long int i = 1; i < i_list.size(); i++)
      if(i_list[i] > i_list[max_index]) {
        max_index = i;
      }

    return i_list[max_index];
  };

  template<class T> T max(const BP_Vec<T> &i_list, BP_Vec<unsigned long int> &max_index) {
    max_index.clear();
    max_index.add(0);
    for(unsigned long int i = 1; i < i_list.size(); i++)
      if(i_list[i] > i_list[max_index[0]]) {
        max_index.clear();
        max_index.add(i);
      }
      else if(i_list[i] == max) {
        max_index.add(i);
      }

    return i_list[max_index[0]];
  };

  template<class T> unsigned long int max_index(const BP_Vec<T> &i_list) {
    unsigned long int max_index = 0;
    for(unsigned long int i = 1; i < i_list.size(); i++)
      if(i_list[i] > i_list[max_index]) {
        max_index = i;
      }

    return max_index;
  };

  template<class T> BP_Vec<unsigned long int> max_indices(const BP_Vec<T> &i_list) {
    BP_Vec<unsigned long int> max_index;
    max_index.add(0);
    for(unsigned long int i = 1; i < i_list.size(); i++)
      if(i_list[i] > i_list[max_index[0]]) {
        max_index.clear();
        max_index.add(i);
      }
      else if(i_list[i] == max) {
        max_index.add(i);
      }

    return max_index;
  };



  template<class T> double mag(const BP_Vec<T> &vec) {
    double mag = 0.0;
    for(unsigned long int i = 0; i < vec.size(); i++)
      mag += vec[i] * vec[i];

    return sqrt(mag);
  }

  template<class T> bool add_once(BP_Vec<T> &vec, const T &i1) {

    for(unsigned long int i = 0; i < vec.size(); i++) {
      if(vec[i] == i1) {
        return false;
      }
    }
    vec.add(i1);
    return true;
  }

  template<class T> bool add_once(BP_Vec<T> &vec, const T &i1, double eps) {

    for(unsigned long int i = 0; i < vec.size(); i++) {
      if(std::fabs(vec[i] - i1) < eps) {
        return false;
      }
    }
    vec.add(i1);
    return true;
  }

  template<class T> void add_once(BP_Vec<T> &vec, const BP_Vec<T> &vec2) {
    for(unsigned long int i = 0; i < vec2.size(); i++)
      add_once(vec, vec2[i]);
  }

  template<class T> unsigned long int mem_size(BP_Vec< T> &vec) {
    unsigned long int mem_sum = 0;

    // add size of val lists
    mem_sum += vec.size_capacity() * sizeof(T *);

    // add size of data
    mem_sum += sizeof(BP_Vec<T>);
    mem_sum += vec.size_alloc() * sizeof(T);

    return mem_sum;
  }

  template<class T> unsigned long int mem_size(BP_Vec< BP_Vec< T> > &vec) {
    unsigned long int i;
    unsigned long int mem_sum = 0;

    // add size of val lists
    mem_sum += vec.size_capacity() * sizeof(BP_Vec<T> *);
    for(i = 0; i < vec.size(); i++)
      mem_sum += vec[i].size_capacity() * sizeof(T *);

    // add size of data
    mem_sum += sizeof(BP_Vec< BP_Vec<T> >);
    mem_sum += vec.size_alloc() * sizeof(BP_Vec<T>);
    for(i = 0; i < vec.size(); i++)
      mem_sum += vec[i].size_alloc() * sizeof(T);

    return mem_sum;
  }

  template<class T> unsigned long int mem_size(BP_Vec< BP_Vec< BP_Vec< T> > > &vec) {
    unsigned long int i, j;
    unsigned long int mem_sum = 0;

    // add size of val lists
    mem_sum += vec.size_capacity() * sizeof(BP_Vec< BP_Vec<T> > *);
    for(i = 0; i < vec.size(); i++) {
      mem_sum += vec[i].size_capacity() * sizeof(BP_Vec<T> *);
      for(j = 0; j < vec[i].size(); j++)
        mem_sum += vec[i][j].size_capacity() * sizeof(T *);
    }

    // add size of data
    mem_sum += sizeof(BP_Vec< BP_Vec< BP_Vec<T> > >);
    mem_sum += vec.size_alloc() * sizeof(BP_Vec< BP_Vec<T> >);
    for(i = 0; i < vec.size(); i++) {
      mem_sum += vec[i].size_alloc() * sizeof(BP_Vec<T>);
      for(j = 0; j < vec[i].size(); j++)
        mem_sum += vec[i][j].size_alloc() * sizeof(T);

    }

    return mem_sum;
  }

  template<class T> void set(BP_Vec<T> &vec, const T &val) {
    for(unsigned long int i = 0; i < vec.size(); i++)
      vec[i] = val;
  }

  template<class T> void set(BP_Vec< BP_Vec<T> > &vec, const T &val) {
    for(unsigned long int i = 0; i < vec.size(); i++)
      set(vec[i], val);
  }

  template<class T> void set(BP_Vec< BP_Vec< BP_Vec<T> > > &vec, const T &val) {
    for(unsigned long int i = 0; i < vec.size(); i++)
      set(vec[i], val);
  }

  template<class T> void set(BP_Vec< BP_Vec< BP_Vec< BP_Vec<T> > > > &vec, const T &val) {
    for(unsigned long int i = 0; i < vec.size(); i++)
      set(vec[i], val);
  }

  template<class T> void set(BP_Vec< BP_Vec< BP_Vec< BP_Vec< BP_Vec<T> > > > > &vec, const T &val) {
    for(unsigned long int i = 0; i < vec.size(); i++)
      set(vec[i], val);
  }

  // matrix multiplication (not elementwise)
  template<class T> const BP_Vec< BP_Vec<T> > operator*(const BP_Vec< BP_Vec<T> > &mA, const BP_Vec< BP_Vec<T> > &mB) {
    unsigned long int r, c, n, i, j, k;
    r = mA.size();
    c = mB[0].size();
    n = mA[0].size();

    BP_Vec< BP_Vec<T> > m = BP_Vec< BP_Vec<T> >(r, BP_Vec<T>(c, 0.0));

    for(i = 0; i < r; i++)
      for(j = 0; j < c; j++)
        for(k = 0; k < n; k++) {
          m[i][j] += mA[i][k] * mB[k][j];
        }

    return m;
  }

  // matrix*vector
  template<class T> const BP_Vec<T> operator*(const BP_Vec< BP_Vec<T> > &mA, const BP_Vec<T> &vB) {
    unsigned long int r, c, n, i, k;
    r = mA.size();
    c = 1;
    n = mA[0].size();

    BP_Vec<T> v = BP_Vec<T>(r, 0.0);

    for(i = 0; i < r; i++)
      for(k = 0; k < n; k++)
        v[i] += mA[i][k] * vB[k];

    return v;
  }

  // vector*scalar
  template<class T> const BP_Vec<T> operator*(const BP_Vec<T> &vA, const T &B) {
    BP_Vec<T> v = vA;
    v *= B;
    return v;
  }

  // matrix*scalar
  template<class T> const BP_Vec<T> operator*(const BP_Vec< BP_Vec<T> > &mA, const T &B) {
    BP_Vec< BP_Vec<T> >  m = mA;
    m *= B;
    return m;
  }

  // vector+scalar
  template<class T> const BP_Vec<T> operator+(const BP_Vec<T> &vA, const T &B) {
    BP_Vec<T> v = vA;
    v += B;
    return v;
  }

  // vector-scalar
  template<class T> const BP_Vec<T> operator-(const BP_Vec<T> &vA, const T &B) {
    BP_Vec<T> v = vA;
    v -= B;
    return v;
  }

  // vector/scalar
  template<class T> const BP_Vec<T> operator/(const BP_Vec<T> &vA, const T &B) {
    BP_Vec<T> v = vA;
    v /= B;
    return v;
  }

  template<class T> const BP_Vec<T> abs(const BP_Vec<T> &vA) {
    BP_Vec<T> v = vA;
    for(unsigned long int i = 0; i < vA.size(); i++)
      v[i] = std::abs(vA[i]);
    return v;
  }

  template<class T> const BP_Vec<T> fabs(const BP_Vec<T> &vA) {
    BP_Vec<T> v = vA;
    for(unsigned long int i = 0; i < vA.size(); i++)
      v[i] = std::fabs(vA[i]);
    return v;
  }

  template<class T> const BP_Vec<T> get_unique(const BP_Vec<T> &v) {
    BP_Vec<T> unique_values;
    unsigned long int index;

    for(unsigned long int i = 0; i < v.size(); i++) {
      if(!unique_values.find_first(v[i], index)) {
        unique_values.add(v[i]);
      }
    }

    return unique_values;
  }

  template<class T> const BP_Vec<T> get_unique(const BP_Vec<T> &v, BP_Vec<unsigned long int> &count) {
    BP_Vec<T> unique_values;
    count.erase();
    unsigned long int index;

    for(unsigned long int i = 0; i < v.size(); i++) {
      if(!unique_values.find_first(v[i], index)) {
        unique_values.add(v[i]);
        count.add(1);
      }
      else {
        count[index]++;
      }
    }

    return unique_values;
  }

}

#endif // BP_Vec_HH

