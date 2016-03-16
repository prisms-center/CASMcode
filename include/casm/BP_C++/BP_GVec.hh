#ifndef BP_GVec_HH
#define BP_GVec_HH

#include <limits>
#include <cstdlib>
#include <iostream>
#include "casm/BP_C++/BP_Vec.hh"

// vocab usage
//
// remove:			remove data from active part of list, does not free memory
// clear:			remove all data from active parts, does not free memory
// reset:			reset all data to initial settings, does not free memory
// erase:			(Member): remove GVec_Member from all lists it is a part of, does not free memory
//					(BP_Vec or BP_GVec): clear() & free() the list
// erase_if_last:	(Member): remove GVec_Member from Group or Graph; if after removal it is no longer in any
//								other Group or Graph, remove it from it's GVec home
// free:			frees unused memory, but does not change capacity
// capacity:		changes capacity, freeing memory if necessary

/// \ingroup BP
namespace BP {

  //GVec
  class BP_Gen_GVec;
  template<class T> class BP_GVec; //: public BP_Gen_GVec

  //GVec_Member
  class BP_Group_Data_class;
  class BP_Graph_Data_class;
  class BP_Gen_GVec_Member;
  template<class T> class BP_GVec_Member; //: public BP_Gen_GVec_Member

  //Group
  class BP_Gen_Group;
  template<class T> class BP_Group; //: public BP_Gen_Group

  //Graph
  class BP_Edge_Data_class;
  class BP_Gen_Vertex;
  class BP_Gen_Edge;
  class BP_Gen_Graph;
  template<class T, class U> class BP_Graph; //: public BP_Gen_Graph


  // //////////////////////
  //  GVec

  /// \ingroup BP BP_GVec
  class BP_Gen_GVec {
  protected:

    unsigned long int N_max;			//lenth of T **val
    unsigned long int N_alloc;		//max i, for which val[i] = new T has been called
    unsigned long int N;				//current max i, (val[i], for N <= i < N_alloc, point to existing objects, but are not currently in use)
    unsigned long int cap_incr;				//amount that capacity is increased if N == N_alloc == N_max when add() or put_in() is called

    BP_Gen_GVec_Member **val;
  public:
    unsigned long int get_index(BP_Gen_GVec_Member *);
    void remove(BP_Gen_GVec_Member *);
    void ordered_remove(BP_Gen_GVec_Member *);
    void remove(unsigned long int);
    void ordered_remove(unsigned long int);
    void swap(unsigned long int, unsigned long int);
    void ordered_swap(unsigned long int i1, unsigned long int i2);
    void clear();
    void crazy_clear();				// potentially unsafe clear, doesn't remove members from groups and graphs,
    //	so you must call crazy_clear() for each group and graph separately
    //  but if you do this right it might be faster
    void free();
    void erase();
    unsigned long int  size() const {
      return N;
    };
    void cap_increment(unsigned long int);


    void init();
    BP_Gen_GVec();
    ~BP_Gen_GVec();

  };

  /// \ingroup BP BP_GVec
  template<class T> class BP_GVec : public BP_Gen_GVec {
  public:


    //BP_GVec<T>& operator=(const BP_GVec<T>&);
    BP_GVec_Member<T> *member(unsigned long int);
    BP_GVec_Member<T> *member(BP_Gen_Vertex *);
    BP_GVec_Member<T> *member(BP_Gen_Edge *);
    BP_GVec_Member<T> *last_member();
    inline T &operator[](unsigned long int);
    inline const T &operator[](unsigned long int) const;
    inline T &operator[](BP_Gen_GVec_Member *);
    inline const T &operator[](BP_Gen_GVec_Member *) const;
    inline T &operator[](BP_Gen_Vertex *);
    inline const T &operator[](BP_Gen_Vertex *) const;
    inline T &operator[](BP_Gen_Edge *);
    inline const T &operator[](BP_Gen_Edge *) const;
    BP_GVec_Member<T> *add();
    BP_GVec_Member<T> *add(const T &);
    template <class U> BP_GVec_Member<T> *add(const U &);
    template <class U> BP_GVec_Member<T> *add_in_place(const U &, unsigned long int);
    inline T &last() {
      return (*this)[N - 1];
    };
    inline const T &last() const {
      return (*this)[N - 1];
    };
    BP_GVec_Member<T> *put_in(BP_GVec_Member<T> *);
    BP_GVec_Member<T> *put_in_place(BP_GVec_Member<T> *, unsigned long int);
    BP_GVec_Member<T> *take_out(BP_Gen_GVec_Member *);
    BP_GVec_Member<T> *take_out(unsigned long int i1);
    BP_GVec_Member<T> *ordered_take_out(BP_Gen_GVec_Member *);
    BP_GVec_Member<T> *ordered_take_out(unsigned long int);

    void capacity(unsigned long int);



  };

  ////////////////////////
  /// GVec_Member

  /// \ingroup BP BP_GVec
  class BP_Group_Data_class {
  public:

    BP_Gen_Group *address;
    unsigned long int index;	//index into group_members list
  };

  /// \ingroup BP BP_GVec
  class BP_Graph_Data_class {
  public:

    BP_Gen_Graph *address;
    int type;	//0: vert; 1: edge
    unsigned long int index;	//index into Vertex_list or Edge_list
  };

  /// \ingroup BP BP_GVec
  class BP_Gen_GVec_Member {
    friend class BP_Gen_GVec;
    template<class> friend class BP_GVec;
    friend class BP_Gen_Group;
    template<class> friend class BP_Group;
    friend class BP_Gen_Graph;
    template<class, class> friend class BP_Graph;

  private:

    BP_Gen_GVec *BP_GVec_home;
    unsigned long int BP_GVec_index;
    BP_Vec< BP_Group_Data_class> groups;
    BP_Vec< BP_Graph_Data_class> graphs;

  public:
    BP_Gen_GVec *GVec_home() const {
      return BP_GVec_home;
    };
    unsigned long int GVec_index() const {
      return BP_GVec_index;
    };
    void set_GVec_home(BP_Gen_GVec *i1) {
      BP_GVec_home = i1;
    };
    void set_GVec_index(unsigned long int i1) {
      BP_GVec_index = i1;
    };

    unsigned long int group_size() const {
      return groups.size();
    };
    unsigned long int graph_size() const {
      return graphs.size();
    };
    BP_Gen_Group *group(unsigned long int i1) const {
      return groups[i1].address;
    };
    BP_Gen_Graph *graph(unsigned long int i1) const {
      return graphs[i1].address;
    };

    void erase();
    void erase_if_last(BP_Gen_Group *i1);
    void erase_if_last(BP_Gen_Graph *i1);
    void crazy_clear();
    BP_Gen_GVec_Member() {
      BP_GVec_home = 0;
    };
    BP_Gen_GVec_Member(BP_Gen_GVec *i1, unsigned long int i2) {
      BP_GVec_home = i1;
      BP_GVec_index = i2;
    };
    virtual ~BP_Gen_GVec_Member();

    void set_group_index(BP_Gen_Group *addy, unsigned long int indy) {
      if(groups[0].address == addy) {
        //std::cout << "addy 0" << std::endl;
        groups[0].index = indy;
        return;
      }
      else if(groups[1].address == addy) {
        //std::cout << "addy 1" << std::endl;
        groups[1].index = indy;
        return;
      }
      else if(groups[2].address == addy) {
        //std::cout << "addy 2" << std::endl;
        groups[2].index = indy;
        return;
      }
      else {

        //std::cout << "addy >2" << std::endl;
        for(unsigned long int i = 3; i < groups.size(); i++) {
          if(groups[i].address == addy) {
            groups[i].index = indy;
            return;
          }
        }

      }
    };

    bool is_in_group(BP_Gen_Group *addy) {
      for(unsigned long int i = 0; i < groups.size(); i++)
        if(groups[i].address == addy)
          return true;

      return false;
    }

    bool is_in_graph(BP_Gen_Graph *addy) {
      for(unsigned long int i = 0; i < graphs.size(); i++)
        if(graphs[i].address == addy)
          return true;

      return false;
    }

  };

  /// \ingroup BP BP_GVec
  template<class T> class BP_GVec_Member : public BP_Gen_GVec_Member {
    friend class BP_Gen_GVec;
    template<class> friend class BP_GVec;

  private:

    bool alloc;
    T *val;

  public:

    BP_GVec_Member<T>() {
      alloc = 0;
    };
    BP_GVec_Member<T>(BP_Gen_GVec *i1, unsigned long int i2) : BP_Gen_GVec_Member(i1, i2) {
      val = new T();
      alloc = 1;
    };
    BP_GVec_Member<T>(BP_Gen_GVec *i1, unsigned long int i2, const T &i3) : BP_Gen_GVec_Member(i1, i2) {
      val = new T(i3);
      alloc = 1;
    };
    template <class U> BP_GVec_Member<T>(BP_Gen_GVec *i1, unsigned long int i2, const U &i3) : BP_Gen_GVec_Member(i1, i2) {
      val = new T(i3);
      alloc = 1;
    };

    BP_GVec_Member<T>(BP_Gen_GVec *i1, unsigned long int i2, T *i3) : BP_Gen_GVec_Member(i1, i2) {
      val = i3;
      alloc = 1;
    };

    virtual ~BP_GVec_Member<T>() {
      //std::cout << "-- Deconstrucor: ~BP_GVec_Member<T>" << std::endl;
      if(alloc) delete val;
    };

    T &get_val() {
      return *val;
    };
    T *getp_val() {
      return val;
    };

  };



  ////////////////////////
  /// Group

  /// \ingroup BP BP_GVec
  class BP_Gen_Group {

  public:
    std::string Name;

  protected:
    BP_Vec<BP_Gen_GVec_Member *> group_members;
    BP_Vec<unsigned long int> group_member_index;


  public:
    virtual inline unsigned long int size() const {
      return group_members.size();
    };
    virtual inline unsigned long int size_alloc() const {
      return group_members.size_alloc();
    };
    virtual inline unsigned long int size_capacity() const {
      return group_members.size_capacity();
    };
    virtual void capacity(unsigned long int);
    virtual void cap_increment(unsigned long int);
    virtual void set_Name(std::string i1) {
      Name = i1;
    };
    virtual std::string get_Name() {
      return Name;
    };
    virtual void erase(unsigned long int i1) {
      group_members[i1]->erase();
    };
    virtual void erase_if_last(unsigned long int i1) {
      group_members[i1]->erase_if_last(this);
    };
    virtual void remove(unsigned long int);
    virtual void remove(BP_Gen_GVec_Member *);
    virtual void remove(const BP_Vec<  BP_Gen_GVec_Member *> &);
    virtual void clear();
    virtual void crazy_clear();
    virtual void free();
    virtual void write(std::ostream &sout);
    virtual unsigned long int get_index(BP_Gen_GVec_Member *);
    virtual bool contains(BP_Gen_GVec_Member *i1);


    BP_Gen_Group();
    ~BP_Gen_Group();

    // testing
    unsigned long int remove_part_one(BP_Gen_GVec_Member *i1);
    void remove_part_two(unsigned long int vi);

  };

  //return *( ((BP_GVec_Member<T>*) val[i1])->val);
  /// \ingroup BP BP_GVec
  template <class T> class BP_Group : public BP_Gen_Group {
  public:
    virtual BP_GVec_Member<T> *add(BP_Gen_GVec_Member *);
    virtual BP_GVec_Member<T> *member(unsigned long int);
    virtual BP_Vec< BP_GVec_Member<T>* > all_members();
    virtual inline T *getp(unsigned long int i1) {
      return ((BP_GVec_Member<T> *) group_members[i1])->getp_val();
    };
    virtual inline T &operator[](unsigned long int i1) {
      return ((BP_GVec_Member<T> *) group_members[i1])->get_val();
    };
    virtual inline const T &operator[](unsigned long int i1)const {
      return ((BP_GVec_Member<T> *) group_members[i1])->get_val();
    };

  };

  ////////////////////////
  /// Graph

  /// \ingroup BP BP_GVec
  class BP_Edge_Data_class {
  public:
    BP_Gen_Edge *address;
    unsigned long int vert_index;	//1: BP_Gen_Edge->a points to this Vertex, 0: BP_Gen_Edge->b points to this Vertex
  };

  /// \ingroup BP BP_GVec
  class BP_Gen_Vertex {
  public:
    unsigned long int vertex_list_index;
    BP_Gen_GVec_Member *val;
    BP_Vec< BP_Edge_Data_class > edges;
  };

  /// \ingroup BP BP_GVec
  class BP_Gen_Edge {
  public:

    unsigned long int edge_list_index;	//index in BP_Gen_Graph.Edge_list

    bool dir;	//1: a->b, 0:a-b (undirected)		('a' being vert_list[0], 'b' being vert_list[1])

    BP_Gen_Vertex *vert_list[2];
    unsigned long int edge_index[2];	//BP_Gen_Vertex vert_list[i]->edges[edge_index[i]] points to this BP_Gen_Edge
    BP_Gen_GVec_Member *val;

  };

  /// \ingroup BP BP_GVec
  class BP_Gen_Graph {
  public:

    std::string Name;

  protected:
    BP_Vec< BP_Gen_Vertex > Vertex_list;
    BP_Vec< BP_Gen_Edge > Edge_list;

    //BP_Vec< int > Vertex_list_index;
    //BP_Vec< int > Edge_list_index;

  public:

    inline unsigned long int vert_size() const {
      return Vertex_list.size();
    };		//number of vertices in graph
    inline unsigned long int edge_size() const {
      return Edge_list.size();
    };			//number of edges in graph
    unsigned long int num_incident_edges(unsigned long int i1);					//number of edges connected to vertex i1
    unsigned long int num_incident_edges(unsigned long int i1, int i2);			//number of edges connected to vertex i1 ( -1:out edges, 0:bidir, 1:in)
    unsigned long int num_incident_edges(BP_Gen_GVec_Member *i1);					//number of edges connected to vertex i1
    unsigned long int num_incident_edges(BP_Gen_GVec_Member *i1, int i2);			//number of edges connected to vertex i1 ( -1:out edges, 0:bidir, 1:in)
    unsigned long int num_incident_edges(BP_Gen_Vertex *i1);					//number of edges connected to vertex i1
    unsigned long int num_incident_edges(BP_Gen_Vertex *i1, int i2);			//number of edges connected to vertex i1 ( -1:out edges, 0:bidir, 1:in)


    /// these virtual functions could be specialized in derived classes

    virtual BP_Gen_Vertex *add_vertex(BP_Gen_GVec_Member *);

    virtual void remove(BP_Gen_GVec_Member *);

    virtual void erase_vert(unsigned long int i1) {
      Vertex_list[i1].val->erase();
    };
    virtual void erase_edge(unsigned long int i1) {
      Edge_list[i1].val->erase();
    };
    virtual void erase(BP_Gen_GVec_Member *i1) {
      i1->erase();
    };
    virtual void erase(BP_Gen_Vertex *i1) {
      i1->val->erase();
    };
    virtual void erase(BP_Gen_Edge *i1) {
      i1->val->erase();
    };

    virtual void erase_if_last_vert(unsigned long int i1) {
      Vertex_list[i1].val->erase_if_last(this);
    };
    virtual void erase_if_last_edge(unsigned long int i1) {
      Edge_list[i1].val->erase_if_last(this);
    };
    virtual void erase_if_last(BP_Gen_GVec_Member *i1) {
      i1->erase_if_last(this);
    };
    virtual void erase_if_last(BP_Gen_Vertex *i1) {
      i1->val->erase_if_last(this);
    };
    virtual void erase_if_last(BP_Gen_Edge *i1) {
      i1->val->erase_if_last(this);
    };


    virtual void remove_vertex(unsigned long int);
    virtual void remove_vertex(BP_Gen_GVec_Member *);
    virtual void remove_vertex(BP_Gen_Vertex *);

    virtual void disconnect_vertex(unsigned long int);
    virtual void disconnect_vertex(BP_Gen_GVec_Member *);
    virtual void disconnect_vertex(BP_Gen_Vertex *);

    virtual void add_halfedge(unsigned long int, BP_Gen_GVec_Member *, int);
    virtual void add_halfedge(BP_Gen_GVec_Member *, BP_Gen_GVec_Member *, int);
    virtual void add_halfedge(BP_Gen_Vertex *, BP_Gen_GVec_Member *, int);

    virtual void connect_vertices(unsigned long int, unsigned long int, BP_Gen_GVec_Member *, int);
    virtual void connect_vertices(BP_Gen_GVec_Member *, BP_Gen_GVec_Member *, BP_Gen_GVec_Member *, int);
    virtual void connect_vertices(BP_Gen_Vertex *, BP_Gen_Vertex *, BP_Gen_GVec_Member *, int);

    virtual void disconnect_vertices(unsigned long int, unsigned long int);
    virtual void disconnect_vertices(BP_Gen_GVec_Member *, BP_Gen_GVec_Member *);
    virtual void disconnect_vertices(BP_Gen_Vertex *, BP_Gen_Vertex *);

    virtual void remove_edge(unsigned long int);
    virtual void remove_edge(BP_Gen_GVec_Member *);
    virtual void remove_edge(BP_Gen_Edge *);

    virtual void clear_edges();
    virtual void clear();
    virtual void crazy_clear();
    virtual void free();
    virtual void capacity(unsigned long int);
    virtual void cap_increment(unsigned long int);

    /// end comment


    unsigned long int vert_to_index(BP_Gen_Vertex *i1);
    unsigned long int edge_to_index(BP_Gen_Edge *i1);
    BP_Gen_Vertex *index_to_vert(unsigned long int i1);
    BP_Gen_Edge *index_to_edge(unsigned long int i1);
    unsigned long int member_vert_to_index(BP_Gen_GVec_Member *i1);
    unsigned long int member_edge_to_index(BP_Gen_GVec_Member *i1);
    BP_Gen_Vertex *member_vert_to_vert(BP_Gen_GVec_Member *i1);
    BP_Gen_Edge *member_edge_to_edge(BP_Gen_GVec_Member *i1);

    unsigned long int get_neighbor(unsigned long int i1, unsigned long int i2);					//vertex i2 connected to vertex i1 considering all edges
    unsigned long int get_neighbor(unsigned long int i1, unsigned long int i2, int i3);			//vertices connected to vertex i1 considering ( -1:out edges, 0:bidir, 1:in)
    BP_Gen_Vertex *get_neighbor(BP_Gen_Vertex *i1, unsigned long int i2);					//vertex i2 connected to vertex i1 considering all edges
    BP_Gen_Vertex *get_neighbor(BP_Gen_Vertex *i1, unsigned long int i2, int i3);			//vertices connected to vertex i1 considering ( -1:out edges, 0:bidir, 1:in)

    unsigned long int get_incident_edge(unsigned long int i1, unsigned long int i2);				//edge i2 connected to vertex i1 considring all edges
    unsigned long int get_incident_edge(unsigned long int i1, unsigned long int i2, int i3);		//edge i2 connected to vertex i1 considering (i3=-1:out edges, 0:bidir, 1:in)
    BP_Gen_Edge *get_incident_edge(BP_Gen_Vertex *i1, unsigned long int i2);				//edge i2 connected to vertex i1 considring all edges
    BP_Gen_Edge *get_incident_edge(BP_Gen_Vertex *i1, unsigned long int i2, int i3);		//edge i2 connected to vertex i1 considering (i3=-1:out edges, 0:bidir, 1:in)

    bool are_1NN(BP_Gen_Vertex *i1, BP_Gen_Vertex *i2);
    bool are_1NN(unsigned long int i1, unsigned long int i2);
    bool are_1NN(BP_Gen_GVec_Member *i1, BP_Gen_GVec_Member *i2);

    // give an edge, return startpoint
    unsigned long int edge_startvert(unsigned long int i1);
    BP_Gen_Vertex *edge_startvert(BP_Gen_Edge *i1);

    // give an edge, return endpoint
    unsigned long int edge_endvert(unsigned long int i1);
    BP_Gen_Vertex *edge_endvert(BP_Gen_Edge *i1);

    bool contains(BP_Gen_GVec_Member *i1);

    virtual void write(std::ostream &sout);

    BP_Gen_Graph();
    ~BP_Gen_Graph();

  protected:
    unsigned long int get_index(BP_Gen_GVec_Member *);	//  get the index in Vertex_list or Edge_list for this member

  };

  //        Vertex, Edge
  /// \ingroup BP BP_GVec
  template<class T, class U> class BP_Graph: public BP_Gen_Graph {
  public:
    BP_GVec_Member<T> *index_to_member_vert(unsigned long int i1);
    BP_GVec_Member<U> *index_to_member_edge(unsigned long int i1);
    BP_GVec_Member<T> *vert_to_member_vert(BP_Gen_Vertex *i1);
    BP_GVec_Member<U> *edge_to_member_edge(BP_Gen_Edge *i1);


    // give an edge, return startpoint
    BP_GVec_Member<T> *edge_startvert(BP_Gen_GVec_Member *i1);

    // give an edge, return endpoint
    BP_GVec_Member<T> *edge_endvert(BP_Gen_GVec_Member *i1);


    inline T &vert_val(unsigned long int i1) {
      return ((BP_GVec_Member<T> *) Vertex_list[i1].val)->get_val();
    };
    inline U &edge_val(unsigned long int i1) {
      return ((BP_GVec_Member<U> *) Edge_list[i1].val)->get_val();
    };
    inline T &vert_val(BP_Gen_Vertex *i1) {
      return ((BP_GVec_Member<T> *) i1->val)->get_val();
    };
    inline U &edge_val(BP_Gen_Edge *i1) {
      return ((BP_GVec_Member<U> *) i1->val)->get_val();
    };
    inline T &vert_val(BP_Gen_GVec_Member *i1) {
      return ((BP_GVec_Member<T> *) i1)->get_val();
    };
    inline U &edge_val(BP_Gen_GVec_Member *i1) {
      return ((BP_GVec_Member<U> *) i1)->get_val();
    };

    inline T *getp_vert_val(unsigned long int i1) {
      return ((BP_GVec_Member<T> *) Vertex_list[i1].val)->val;
    };
    inline U *getp_edge_val(unsigned long int i1) {
      return ((BP_GVec_Member<U> *) Edge_list[i1].val)->val;
    };
    inline T *getp_vert_val(BP_Gen_Vertex *i1) {
      return ((BP_GVec_Member<T> *) i1->val)->val;
    };
    inline U *getp_edge_val(BP_Gen_Edge *i1) {
      return ((BP_GVec_Member<U> *) i1->val)->val;
    };
    inline T *getp_vert_val(BP_Gen_GVec_Member *i1) {
      return ((BP_GVec_Member<T> *) i1)->val;
    };
    inline U *getp_edge_val(BP_Gen_GVec_Member *i1) {
      return ((BP_GVec_Member<U> *) i1)->val;
    };

    using BP_Gen_Graph::get_neighbor;
    BP_GVec_Member<T> *get_neighbor(BP_Gen_GVec_Member *i1, unsigned long int i2);					//vertex i2 connected to vertex i1 considering all edges
    BP_GVec_Member<T> *get_neighbor(BP_Gen_GVec_Member *i1, unsigned long int i2, int i3);			//vertices connected to vertex i1 considering ( -1:out edges, 0:bidir, 1:in)

    using BP_Gen_Graph::get_incident_edge;
    BP_GVec_Member<U> *get_incident_edge(BP_Gen_GVec_Member *i1, unsigned long int i2);				//edge i2 connected to vertex i1 considring all edges
    BP_GVec_Member<U> *get_incident_edge(BP_Gen_GVec_Member *i1, unsigned long int i2, int i3);		//edge i2 connected to vertex i1 considering (i3=-1:out edges, 0:bidir, 1:in)
    //T* operator[](const BP_Vert_Iter &i1){return ((T*) Vertex_list[i1.val]);};
    //U* operator[](const BP_Edge_Iter &i1){return ((U*) Edge_list[i1.val]);};

  };



  ////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////
  //// Functions /////////////////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////////////////////////////
  /////// BP_GVec<T> ////////////////////////////////////////////////////////////////////////////

  template<class T> BP_GVec_Member<T> *BP_GVec<T>::member(unsigned long int i1) {
    //if( i1 >= N) std::cout << "BP_GVec index >= N" << std::endl;
    //if( i1 >= N_alloc) std::cout << "BP_GVec index >= N_alloc" << std::endl;
    return (BP_GVec_Member<T> *) val[i1];
    //return *(val[i1]->val);
  };

  template<class T> BP_GVec_Member<T> *BP_GVec<T>::member(BP_Gen_Vertex *i1) {
    BP_GVec_Member<T> *gm = (BP_GVec_Member<T> *) i1->val;
    if(gm->BP_GVec_home != this) {
      std::cout << " home:  " << gm->BP_GVec_home << " != this: " << this << ", return 0" << std::endl;
      return 0;
    }
    else return gm;
  };

  template<class T> BP_GVec_Member<T> *BP_GVec<T>::member(BP_Gen_Edge *i1) {
    BP_GVec_Member<T> *gm = (BP_GVec_Member<T> *) i1->val;
    if(gm->BP_GVec_home != this) return 0;
    else return gm;
  };

  template<class T> BP_GVec_Member<T> *BP_GVec<T>::last_member() {
    //if( i1 >= N) std::cout << "BP_GVec index >= N" << std::endl;
    //if( i1 >= N_alloc) std::cout << "BP_GVec index >= N_alloc" << std::endl;
    return (BP_GVec_Member<T> *) val[N - 1];
    //return *(val[i1]->val);
  };

  template<class T> inline T &BP_GVec<T>::operator[](unsigned long int i1) {
    //if( i1 >= N) std::cout << "BP_GVec index >= N" << std::endl;
    //if( i1 >= N_alloc) std::cout << "BP_GVec index >= N_alloc" << std::endl;
    return ((BP_GVec_Member<T> *) val[i1])->get_val();
    //return *(val[i1]->val);
  };

  template<class T> inline const T &BP_GVec<T>::operator[](unsigned long int i1) const {
    //if( i1 >= N) std::cout << "BP_GVec index >= N" << std::endl;
    //if( i1 >= N_alloc) std::cout << "BP_GVec index >= N_alloc" << std::endl;
    return *(((BP_GVec_Member<T> *) val[i1])->val);
    //return *(val[i1]->val);
  };

  template<class T> inline T &BP_GVec<T>::operator[](BP_Gen_GVec_Member *i1) {
    //if( i1->BP_GVec_index >= N) std::cout << "BP_GVec index >= N" << std::endl;
    //if( i1->BP_GVec_index >= N_alloc) std::cout << "BP_GVec index >= N_alloc" << std::endl;
    return *(((BP_GVec_Member<T> *) val[i1->BP_GVec_index])->val);
    //return *(val[i1]->val);
  };

  template<class T> inline const T &BP_GVec<T>::operator[](BP_Gen_GVec_Member *i1) const {
    //if( i1->BP_GVec_index >= N) std::cout << "BP_GVec index >= N" << std::endl;
    //if( i1->BP_GVec_index >= N_alloc) std::cout << "BP_GVec index >= N_alloc" << std::endl;
    return *(((BP_GVec_Member<T> *) val[i1->BP_GVec_index])->val);
    //return *(val[i1]->val);
  };

  template<class T> inline T &BP_GVec<T>::operator[](BP_Gen_Vertex *i1) {
    BP_GVec_Member<T> *gm = (BP_GVec_Member<T> *) i1->val;
    //if( gm->BP_GVec_index >= N) std::cout << "BP_GVec index: " << gm->BP_GVec_index << " >= N: " << N << std::endl;
    //if( gm->BP_GVec_index >= N_alloc) std::cout << "BP_GVec index: " << gm->BP_GVec_index << " >= N_alloc: " << N_alloc << std::endl;
    return *(gm->val);
    //return *(val[i1]->val);
  };

  template<class T> inline const T &BP_GVec<T>::operator[](BP_Gen_Vertex *i1) const {
    BP_GVec_Member<T> *gm = (BP_GVec_Member<T> *) i1->val;
    //if( gm->BP_GVec_index >= N) std::cout << "BP_GVec index >= N" << std::endl;
    //if( gm->BP_GVec_index >= N_alloc) std::cout << "BP_GVec index >= N_alloc" << std::endl;
    return *(gm->val);
    //return *(val[i1]->val);
  };

  template<class T> inline T &BP_GVec<T>::operator[](BP_Gen_Edge *i1) {
    BP_GVec_Member<T> *gm = (BP_GVec_Member<T> *) i1->val;
    //if( gm->BP_GVec_index >= N) std::cout << "BP_GVec index >= N" << std::endl;
    //if( gm->BP_GVec_index >= N_alloc) std::cout << "BP_GVec index >= N_alloc" << std::endl;
    return *(gm->val);
    //return *(val[i1]->val);
  };

  template<class T> inline const T &BP_GVec<T>::operator[](BP_Gen_Edge *i1) const {
    BP_GVec_Member<T> *gm = (BP_GVec_Member<T> *) i1->val;
    //if( gm->BP_GVec_index >= N) std::cout << "BP_GVec index >= N" << std::endl;
    //if( gm->BP_GVec_index >= N_alloc) std::cout << "BP_GVec index >= N_alloc" << std::endl;
    return *(gm->val);
    //return *(val[i1]->val);
  };

  // Change N_max to i1
  //  uses one "new T*[i1]", no "new T"
  //  calls "delete val[i]" for i>=i1
  //  val[i], where i<i1 remain unchanged
  template<class T> void BP_GVec<T>::capacity(unsigned long int i1) {
    //std::cout << "begin capacity(" << i1 << ")" << "  N_max: " << N_max << std::endl;
    if(i1 == N_max) return;
    //std::cout << "  increase" << std::endl;

    unsigned long int i;
    BP_Gen_GVec_Member **tmp;
    tmp = new BP_Gen_GVec_Member*[i1];


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

    //std::cout << "finish capacity(" << i1 << ")" << std::endl;


  };


  // Put object pointed to into BP_GVec without making a copy
  // could be dangerous... results in two pointers to same object, and now only this one can be used to delete
  template<class T> BP_GVec_Member<T> *BP_GVec<T>::put_in(BP_GVec_Member<T> *i1) {
    if(N == N_max) {
      capacity(N_max + cap_incr);
    }

    i1->BP_GVec_home = this;
    i1->BP_GVec_index = N;

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

    return (BP_GVec_Member<T> *) val[N - 1];
  };

  // Put object pointed to into BP_GVec at position i2 without making a copy
  // could be dangerous... results in two pointers to same object, and now only this one can be used to delete
  template<class T> BP_GVec_Member<T> *BP_GVec<T>::put_in_place(BP_GVec_Member<T> *i1, unsigned long int i2) {
    put_in(i1);

    ordered_swap(N - 1, i2);

    return (BP_GVec_Member<T> *) val[i2];
  };

  // Return pointer to object at val[i1], eliminate pointer
  // could be dangerous... removes pointer without freeing memory, which must now be done elsewhere
  template<class T> BP_GVec_Member<T> *BP_GVec<T>::take_out(BP_Gen_GVec_Member *i1) {
    if(i1->GVec_home() == this) return take_out(i1->GVec_index());
    else return 0;
  };

  // Return pointer to object at val[i1], eliminate pointer
  // could be dangerous... removes pointer without freeing memory, which must now be done elsewhere
  template<class T> BP_GVec_Member<T> *BP_GVec<T>::take_out(unsigned long int i1) {
    //std::cout << "begin take_out(" << i1 << ")" << std::endl;
    BP_GVec_Member<T> *tM = (BP_GVec_Member<T> *) val[i1];
    tM->BP_GVec_home = 0;
    tM->BP_GVec_index = -1;

    swap(i1, N - 1);
    swap(N - 1, N_alloc - 1);
    val[N_alloc - 1] = nullptr;

    N--;
    N_alloc--;


    return tM;
  };


  // Return pointer to object at val[i1], eliminate pointer
  // could be dangerous... removes pointer without freeing memory, which must now be done elsewhere
  template<class T> BP_GVec_Member<T> *BP_GVec<T>::ordered_take_out(BP_Gen_GVec_Member *i1) {
    if(i1->GVec_home() == this) return ordered_take_out(i1->GVec_index());
    else return 0;
  };

  // Return pointer to object at val[i1], eliminate pointer, shift all down one space
  // could be dangerous... removes pointer without freeing memory, which must now be done elsewhere
  template<class T> BP_GVec_Member<T> *BP_GVec<T>::ordered_take_out(unsigned long int i1) {
    ordered_swap(i1, N - 1);
    return take_out(N - 1);
  };

  // add an object (as is if memory already allocated, or initialized with T() if memory needs to be allocated) and return a pointer to it
  template<class T> BP_GVec_Member<T> *BP_GVec<T>::add() {
    if(N == N_max) {
      capacity(N_max + cap_incr);
    }

    if(N == N_alloc) {
      val[N] = new BP_GVec_Member<T>(this, N);

      N++;
      N_alloc++;
    }
    else {
      N++;
    }

    return (BP_GVec_Member<T> *) val[N - 1];
  };

  // add to BP_GVec (like push_back())
  template<class T> BP_GVec_Member<T> *BP_GVec<T>::add(const T &i1) {
    //std::cout << "begin BP_GVec<T>::add()" << std::endl;

    if(N == N_max) {
      capacity(N_max + cap_incr);
    }

    if(N == N_alloc) {
      val[N] = new BP_GVec_Member<T>(this, N, i1);

      N++;
      N_alloc++;
    }
    else {
      ((BP_GVec_Member<T> *) val[N])->get_val() = i1;
      val[N]->set_GVec_home(this);
      val[N]->set_GVec_index(N);

      N++;
    }

    //std::cout << "finish BP_GVec<T>::add()" << std::endl;

    return (BP_GVec_Member<T> *) val[N - 1];
  };

  // add to BP_GVec (like push_back())
  template <class T> template <class U> BP_GVec_Member<T> *BP_GVec<T>::add(const U &i1) {
    if(N == N_max) {
      capacity(N_max + cap_incr);
    }

    if(N == N_alloc) {
      val[N] = new BP_GVec_Member<T>(this, N, i1);

      N++;
      N_alloc++;
    }
    else {
      *(val[N]) = BP_GVec_Member<T>(this, N, i1);

      N++;
    }

    return (BP_GVec_Member<T> *) val[N - 1];
  };

  // add to BP_GVec in position i2
  template <class T> template <class U> BP_GVec_Member<T> *BP_GVec<T>::add_in_place(const U &i1, unsigned long int i2) {
    add(i1);

    ordered_swap(N - 1, i2);

    return (BP_GVec_Member<T> *) val[i2];
  };

  ////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////
  /// BP_Group<T> member functions

  template<class T> BP_GVec_Member<T> *BP_Group<T>::member(unsigned long int i1) {
    return (BP_GVec_Member<T> *) group_members[i1];
  };

  template<class T> BP_Vec< BP_GVec_Member<T>* > BP_Group<T>::all_members() {
    //return (BP_Vec< BP_GVec_Member<T>* >) group_members;

    BP_Vec< BP_GVec_Member<T>* > group_members_list;
    for(unsigned long int i = 0; i < size(); i++)
      group_members_list.add(member(i));

    return group_members_list;

  };

  template<class T> BP_GVec_Member<T> *BP_Group<T>::add(BP_Gen_GVec_Member *i1) {
    for(unsigned long int i = 0; i < i1->group_size(); i++)
      if(i1->group(i) == this) return 0;

    BP_Group_Data_class *gd = i1->groups.add();

    gd->address = this;
    gd->index = group_members.size();

    group_members.add(i1);

    return (BP_GVec_Member<T> *) i1;
  };

  /////////////////////////////////////////////////////////////////////////
  /// BP_Graph<T,U> member functions

  template<class T, class U> BP_GVec_Member<T> *BP_Graph<T, U>::get_neighbor(BP_Gen_GVec_Member *i1, unsigned long int i2) {
    //return get_neighbor( get_index(i1), i2);
    //return (BP_GVec_Member<T>*) Vertex_list[get_neigbhor( get_index(i1), i2)].val;
    return (BP_GVec_Member<T> *) get_neighbor(&Vertex_list[get_index(i1)], i2)->val;

  };

  template<class T, class U> BP_GVec_Member<T> *BP_Graph<T, U>::get_neighbor(BP_Gen_GVec_Member *i1, unsigned long int i2, int i3) {
    //return get_neighbor( get_index(i1), i2, i3);
    return (BP_GVec_Member<T> *) get_neighbor(&Vertex_list[get_index(i1)], i2, i3)->val;
  };

  template<class T, class U> BP_GVec_Member<U> *BP_Graph<T, U>::get_incident_edge(BP_Gen_GVec_Member *i1, unsigned long int i2) {
    return (BP_GVec_Member<U> *) get_incident_edge(&Vertex_list[get_index(i1)], i2)->val;
  };

  template<class T, class U> BP_GVec_Member<U> *BP_Graph<T, U>::get_incident_edge(BP_Gen_GVec_Member *i1, unsigned long int i2, int i3) {
    return (BP_GVec_Member<U> *) get_incident_edge(&Vertex_list[get_index(i1)], i2, i3)->val;
  };

  template<class T, class U> BP_GVec_Member<T> *BP_Graph<T, U>::index_to_member_vert(unsigned long int i1) {
    return (BP_GVec_Member<T> *) Vertex_list[i1].val;
  };

  template<class T, class U> BP_GVec_Member<U> *BP_Graph<T, U>::index_to_member_edge(unsigned long int i1) {
    return (BP_GVec_Member<U> *) Edge_list[i1].val;
  };

  template<class T, class U> BP_GVec_Member<T> *BP_Graph<T, U>::vert_to_member_vert(BP_Gen_Vertex *i1) {
    return (BP_GVec_Member<T> *) i1->val;
  };

  template<class T, class U> BP_GVec_Member<U> *BP_Graph<T, U>::edge_to_member_edge(BP_Gen_Edge *i1) {
    return (BP_GVec_Member<U> *) i1->val;
  };

  template<class T, class U> BP_GVec_Member<T> *BP_Graph<T, U>::edge_startvert(BP_Gen_GVec_Member *i1) {
    // give an edge, return startpoint

    //std::cout << "begin edge_startvert( BP_Gen_GVec_Member*)" << std::endl;
    //std::cout << "  get_index: " << get_index(i1) << std::endl;
    //std::cout << "  member_edge_to_edge(i1): " << member_edge_to_edge(i1) << std::endl;
    //std::cout << "  BP_Gen_Graph::edge_startvert(): " << BP_Gen_Graph::edge_startvert( member_edge_to_edge(i1) ) << std::endl;
    //std::cout << "  vert_to_member_vert(): " << vert_to_member_vert( BP_Gen_Graph::edge_startvert( member_edge_to_edge(i1) ) ) << std::endl;

    return vert_to_member_vert(BP_Gen_Graph::edge_startvert(member_edge_to_edge(i1)));
  };

  template<class T, class U> BP_GVec_Member<T> *BP_Graph<T, U>::edge_endvert(BP_Gen_GVec_Member *i1) {
    // give an edge, return endpoint
    return vert_to_member_vert(BP_Gen_Graph::edge_endvert(member_edge_to_edge(i1)));
  };


  //////////////////////////////////////
  /// Bonus functions
  template<class T> BP_GVec_Member<T> *add_once(BP_GVec<T> &gvec, const T &i1) {
    for(unsigned long int i = 0; i < gvec.size(); i++) {
      if(gvec[i] == i1) {
        return nullptr;
      }
    }
    return gvec.add(i1);
  }


}

#endif // BP_GVec_HH

