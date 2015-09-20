#ifndef BP_GVec_CC
#define BP_GVec_CC

#include "casm/BP_C++/BP_GVec.hh"

namespace BP {


  ////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////
  //// Functions /////////////////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////////////////////////////
  //// BP_Gen_GVec_Member /////////////////////////////////////////////////////////////////////////////

  void BP_Gen_GVec_Member::erase() {
    //std::cout << "begin BP_Gen_GVec_Member::erase" << std::endl;
    while(groups.size() > 0)
      groups[0].address->remove(this);
    while(graphs.size() > 0)
      graphs[0].address->remove(this);
    if(BP_GVec_home != NULL) BP_GVec_home->remove(BP_GVec_index);
  };

  void BP_Gen_GVec_Member::erase_if_last(BP_Gen_Graph *i1) {
    i1->remove(this);

    if(groups.size() == 0 && graphs.size() == 0 && BP_GVec_home != NULL)
      BP_GVec_home->remove(BP_GVec_index);
  };

  void BP_Gen_GVec_Member::erase_if_last(BP_Gen_Group *i1) {
    i1->remove(this);

    if(groups.size() == 0 && graphs.size() == 0 && BP_GVec_home != NULL)
      BP_GVec_home->remove(BP_GVec_index);
  };

  void BP_Gen_GVec_Member::crazy_clear() {
    //std::cout << "begin BP_Gen_GVec_Member::crazy_clear()" << std::endl;
    groups.clear();
    graphs.clear();

  };

  BP_Gen_GVec_Member::~BP_Gen_GVec_Member() {
    //std::cout << "**Deconstructor: ~BP_Gen_GVec_Member" << std::endl;
    while(groups.size() > 0)
      groups[0].address->remove(this);
    while(graphs.size() > 0)
      graphs[0].address->remove(this);

  };

  ////////////////////////////////////////////////////////////////////////////////////////////
  //// BP_Gen_GVec /////////////////////////////////////////////////////////////////////////////

  //  get the index in Vertex_list or Edge_list for i1
  unsigned long int BP_Gen_GVec::get_index(BP_Gen_GVec_Member *i1) {
    if(i1->GVec_home() == this) return i1->GVec_index();
    else return std::numeric_limits<unsigned long int>::max();

  };

  void BP_Gen_GVec::ordered_swap(unsigned long int i1, unsigned long int i2) {
    if(i1 < i2) {
      BP_Gen_GVec_Member *tmp = val[i1];
      for(unsigned long int i = i1; i < i2; i++) {
        val[i] = val[i + 1];
        val[i]->set_GVec_index(i);
      }
      val[i2] = tmp;
      val[i2]->set_GVec_index(i2);
    }
    else {
      BP_Gen_GVec_Member *tmp = val[i1];
      for(unsigned long int i = i1; i > i2; i--) {
        val[i] = val[i - 1];
        val[i]->set_GVec_index(i);
      }
      val[i2] = tmp;
      val[i2]->set_GVec_index(i2);
    }

  };

  void BP_Gen_GVec::remove(BP_Gen_GVec_Member *i1) {
    if(i1->GVec_home() == this) remove(i1->GVec_index());
  };

  void BP_Gen_GVec::ordered_remove(BP_Gen_GVec_Member *i1) {
    if(i1->GVec_home() == this) ordered_remove(i1->GVec_index());
  };

  // Removes object val[i1] from list, moves val[N-1] to val[i1].  Does not actually free memory.
  void BP_Gen_GVec::remove(unsigned long int i1) {
    //std::cout << "begin BP_Gen_GVec::remove()" << std::endl;
    BP_Gen_GVec_Member *gm = val[i1];
    while(gm->groups.size() > 0)
      gm->groups[0].address->remove(gm);
    while(gm->graphs.size() > 0)
      gm->graphs[0].address->remove(gm);

    //val[i1]->BP_GVec_home = 0;
    //val[i1]->BP_GVec_index = -1;
    swap(i1, N - 1);
    N--;

  };

  // Removes object val[i1] from list, moves val[i1+1] to val[N-1] down one index.  Does not actually free memory.
  void BP_Gen_GVec::ordered_remove(unsigned long int i1) {
    BP_Gen_GVec_Member *tmp = val[i1];
    while(tmp->groups.size() > 0)
      tmp->groups[0].address->remove(tmp);
    while(tmp->graphs.size() > 0)
      tmp->graphs[0].address->remove(tmp);

    //tmp->BP_GVec_home = 0;
    //tmp->BP_GVec_index = -1;
    for(unsigned long int i = i1; i < N - 1; i++) {
      val[i] = val[i + 1];
      val[i]->BP_GVec_index = i;
    }
    val[N - 1] = tmp;
    N--;

  };

  // Deletes objects stored in val[i] i>=N.  Frees unused memory.
  void BP_Gen_GVec::free() {
    for(unsigned long int i = N; i < N_alloc; i++)
      delete val[i];
    N_alloc = N;

  };

  // Swap val[a] and val[b]
  void BP_Gen_GVec::swap(unsigned long int a, unsigned long int b) {
    BP_Gen_GVec_Member *tmp = val[a];

    val[a] = val[b];
    val[a]->BP_GVec_index = a;

    val[b] = tmp;
    val[b]->BP_GVec_index = b;

  };

  // Set N=0. Does not actually delete memory.
  void BP_Gen_GVec::clear() {
    while(N > 0)
      remove((unsigned long int) 0);
  };

  void BP_Gen_GVec::crazy_clear() {
    for(unsigned long int i = 0; i < N; i++)
      val[i]->crazy_clear();

    N = 0;
  };

  void BP_Gen_GVec::erase() {
    clear();
    free();

  };

  void BP_Gen_GVec::cap_increment(unsigned long int i1) {
    if(i1 > 0)
      cap_incr = i1;
  };

  ///////////////////////////////////
  /// Constructors & Destructor

  void BP_Gen_GVec::init() {
    N_max = 0;
    N_alloc = 0;
    N = 0;
    cap_incr = 10;
  };

  BP_Gen_GVec::BP_Gen_GVec() {
    init();

  };

  BP_Gen_GVec::~BP_Gen_GVec() {
    //std::cout << "deconstructor: ~BP_Gen_GVec" << std::endl;
    //std::cout << "  N_alloc: " << N_alloc << std::endl;
    for(unsigned long int i = 0; i < N_alloc; i++) {
      delete val[i];
    }


    if(N_max != 0) delete [] val;

  };

  ////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////
  /// BP_Gen_Group functions

  BP_Gen_Group::BP_Gen_Group() {

  };

  BP_Gen_Group::~BP_Gen_Group() {
    //std::cout << "begin ~BP_Gen_Group()" << std::endl;
    while(group_members.size() > 0) {
      BP_Gen_Group::remove(group_members[0]);
    };
    //std::cout << "finish ~BP_Gen_Group()" << std::endl;

  };

  //  get the index in Vertex_list or Edge_list for i1
  unsigned long int BP_Gen_Group::get_index(BP_Gen_GVec_Member *i1) {
    for(unsigned long int i = 0; i < i1->groups.size(); i++)
      if(i1->groups[i].address == this) {
        return i1->groups[i].index;
      }

    return std::numeric_limits<unsigned long int>::max();

  };

  bool BP_Gen_Group::contains(BP_Gen_GVec_Member *i1) {
    for(unsigned long int i = 0; i < i1->groups.size(); i++)
      if(i1->groups[i].address == this) {
        return true;
      }

    return false;
  };

  void BP_Gen_Group::remove(unsigned long int i1) {
    //std::cout << "begin BP_Gen_Group::remove(int)" << std::endl;
    BP_Gen_Group::remove(group_members[i1]);
  };



  unsigned long int BP_Gen_Group::remove_part_one(BP_Gen_GVec_Member *i1) {
    unsigned long int i, vi;

    //remove this group from BP_Gen_GVec_Member
    for(i = 0; i < i1->groups.size(); i++) {
      if(i1->groups[i].address == this) {
        vi = i1->groups[i].index;
        i1->groups.remove(i);
        return vi;

      }
    }


    return -1;
  };

  void BP_Gen_Group::remove_part_two(unsigned long int vi) {
    //update vertex_list_index of vertex now at vi
    if(vi < group_members.size()) {


      //std::cout << "group_members[vi]->groups.size(): " << group_members[vi]->groups.size() << std::endl;
      for(unsigned long int i = 0; i < group_members[vi]->groups.size(); i++)
        if(group_members[vi]->groups[i].address == this) {
          group_members[vi]->groups[i].index = vi;
          return;
          //break;
        }

    }
  };

  void BP_Gen_Group::remove(BP_Gen_GVec_Member *i1) {
    //std::cout << "begin BP_Gen_Group::remove(member)" << std::endl;
    unsigned long int i, vi;

    //vi = remove_part_one(i1);

    //remove this group from BP_Gen_GVec_Member
    for(i = 0; i < i1->groups.size(); i++) {
      if(i1->groups[i].address == this) {
        vi = i1->groups[i].index;
        i1->groups.remove(i);
        break;
      }
    }


    //remove this member from group_members
    group_members.remove(vi);

    if(vi < group_members.size()) group_members[vi]->set_group_index(this, vi);


    //remove_part_two( vi);
    /*
    //update vertex_list_index of vertex now at vi
    if( vi < group_members.size())
    {
    	for( i=0;i<group_members[vi]->groups.size();i++)
    		if( group_members[vi]->groups[i].address == this)
    		{
    			group_members[vi]->groups[i].index = vi;
    			break;
    		}
    }
    */

  };

  void BP_Gen_Group::remove(const BP_Vec<BP_Gen_GVec_Member *> &i1) {
    for(unsigned long int i = 0; i < i1.size(); i++) {
      BP_Gen_Group::remove(i1[i]);
    }

  };

  void BP_Gen_Group::clear() {
    while(group_members.size() > 0) {
      remove(group_members[group_members.size() - 1]);
    }
  };

  void BP_Gen_Group::crazy_clear() {
    group_members.clear();
    //group_member_index.clear();
  };

  void BP_Gen_Group::free() {
    group_members.free();
    //group_member_index.free();
  };

  void BP_Gen_Group::capacity(unsigned long int i1) {
    group_members.capacity(i1);
    //group_member_index.capacity(i1);
  };

  void BP_Gen_Group::cap_increment(unsigned long int i1) {
    group_members.cap_increment(i1);
    //group_member_index.cap_increment(i1);
  };

  void BP_Gen_Group::write(std::ostream &sout) {
    sout << Name << " contains: " << std::endl;
    for(unsigned long int i = 0; i < group_members.size(); i++)
      sout << group_members[i] << " ";
    sout << std::endl;
  }


  ////////////////////////////////////////////
  /// BP_Gen_Graph member functions

  BP_Gen_Graph::BP_Gen_Graph() {

  };

  BP_Gen_Graph::~BP_Gen_Graph() {
    //std::cout << "begin ~BP_Gen_Graph()" << std::endl;
    while(Edge_list.size() > 0)
      remove_edge((unsigned long int) 0);
    while(Vertex_list.size() > 0)
      remove_vertex((unsigned long int) 0);
    //std::cout << "finish ~BP_Gen_Graph()" << std::endl;

  };

  //vertex i2 connected to vertex i1 considering all edges
  unsigned long int BP_Gen_Graph::get_neighbor(unsigned long int i1, unsigned long int i2) {
    //BP_Gen_Vertex *v = &Vertex_list[i1];
    //BP_Gen_Edge *e = v->edges[i2].address;

    //if( e->vert_list[0] == v) return ((BP_GVec_Member<T>*) e->vert_list[1]->val)->val;
    //else return ((BP_GVec_Member<T>*) e->vert_list[0]->val)->val;
    //if( e->vert_list[0] == v) return e->vert_list[1]->vertex_list_index;
    //else return e->vert_list[0]->vertex_list_index;

    return get_neighbor(&Vertex_list[i1], i2)->vertex_list_index;

  };

  //vertex i2 connected to vertex i1 considering (i3= -1:out edges, 0:bidir, 1:in)
  unsigned long int BP_Gen_Graph::get_neighbor(unsigned long int i1, unsigned long int i2, int i3) {
    return get_neighbor(&Vertex_list[i1], i2, i3)->vertex_list_index;
  };

  //vertex i2 connected to vertex i1 considering all edges
  BP_Gen_Vertex  *BP_Gen_Graph::get_neighbor(BP_Gen_Vertex *i1, unsigned long int i2) {
    BP_Gen_Edge *e = i1->edges[i2].address;

    if(e->vert_list[0] == i1) return e->vert_list[1];
    else return e->vert_list[0];

  };

  //vertex i2 connected to vertex i1 considering particular types of edges (i3= -1:out edges, 0:bidir, 1:in)
  BP_Gen_Vertex *BP_Gen_Graph::get_neighbor(BP_Gen_Vertex *i1, unsigned long int i2, int i3) {
    //return get_neighbor( i1->vertex_list_index, i2, i3);

    //BP_Gen_Vertex *v = &Vertex_list[i1];
    BP_Gen_Edge *e;
    unsigned long int i;
    long int sum = -1;
    if(i3 == 0) {
      for(i = 0; i < i1->edges.size(); i++) {
        e = i1->edges[i].address;
        if(e->dir == 0) {
          sum++;
          if(sum == i2) {
            //if( e->vert_list[0] == v) return ((BP_GVec_Member<T>*) e->vert_list[1]->val)->val;
            //else return ((BP_GVec_Member<T>*) e->vert_list[0]->val)->val;
            if(e->vert_list[0] == i1) return e->vert_list[1];
            else return e->vert_list[0];

          }
        }
      }

    }
    else if(i3 == -1) {
      for(i = 0; i < i1->edges.size(); i++) {
        e = i1->edges[i].address;
        if(e->dir == 1 & e->vert_list[0] == i1) {
          sum++;
          if(sum == i2) {
            //if( e->vert_list[0] == v) return ((BP_GVec_Member<T>*) e->vert_list[1]->val)->val;
            //else return ((BP_GVec_Member<T>*) e->vert_list[0]->val)->val;
            if(e->vert_list[0] == i1) return e->vert_list[1];
            else return e->vert_list[0];

          }
        }
      }

    }
    else if(i3 == 1) {
      for(i = 0; i < i1->edges.size(); i++) {
        e = i1->edges[i].address;
        if(e->dir == 1 & e->vert_list[1] == i1) {
          sum++;
          if(sum == i2) {
            //if( e->vert_list[0] == v) return ((BP_GVec_Member<T>*) e->vert_list[1]->val)->val;
            //else return ((BP_GVec_Member<T>*) e->vert_list[0]->val)->val;
            if(e->vert_list[0] == i1) return e->vert_list[1];
            else return e->vert_list[0];

          }
        }
      }

    }
    else
      return 0;

    return 0;
  };

  //return edge i2 connected to vertex i1 considring all edges
  unsigned long int BP_Gen_Graph::get_incident_edge(unsigned long int i1, unsigned long int i2) {
    //return ((BP_GVec_Member<U>*) Vertex_list[i1].edges[i2].address->val)->val;
    return Vertex_list[i1].edges[i2].address->edge_list_index;
    //return get_incident_edge( &Vertex_list[i1], i2)->edge_list_index;
  };

  //return edge i2 connected to vertex i1 considering (i3=-1:out edges, 0:bidir, 1:in)
  unsigned long int BP_Gen_Graph::get_incident_edge(unsigned long int i1, unsigned long int i2, int i3) {
    return get_incident_edge(&Vertex_list[i1], i2, i3)->edge_list_index;

  };

  //return edge i2 connected to vertex i1 considring all edges
  BP_Gen_Edge *BP_Gen_Graph::get_incident_edge(BP_Gen_Vertex *i1, unsigned long int i2) {
    //return get_incident_edge( i1->vertex_list_index, i2);
    return i1->edges[i2].address;
  };

  //return edge i2 connected to vertex i1 considring all edges
  BP_Gen_Edge *BP_Gen_Graph::get_incident_edge(BP_Gen_Vertex *i1, unsigned long int i2, int i3) {
    //return get_incident_edge( i1->vertex_list_index, i2, i3);
    //BP_Gen_Vertex *v = &Vertex_list[i1];
    BP_Gen_Edge *e;
    unsigned long int i;
    long int sum = -1;
    if(i3 == 0) {
      for(i = 0; i < i1->edges.size(); i++) {
        e = i1->edges[i].address;
        if(e->dir == 0) {
          sum++;
          if(sum == i2) {
            //return ((BP_GVec_Member<U>*) e->val)->val;
            return e;
          }
        }
      }

    }
    else if(i3 == -1) {
      for(i = 0; i < i1->edges.size(); i++) {
        e = i1->edges[i].address;
        if(e->dir == 1 & e->vert_list[0] == i1) {
          sum++;
          if(sum == i2) {
            //return ((BP_GVec_Member<U>*) e->val)->val;
            return e;
          }
        }
      }

    }
    else if(i3 == 1) {
      for(i = 0; i < i1->edges.size(); i++) {
        e = i1->edges[i].address;
        if(e->dir == 1 & e->vert_list[1] == i1) {
          sum++;
          if(sum == i2) {
            //return ((BP_GVec_Member<U>*) e->val)->val;
            return e;
          }
        }
      }

    }
    else
      return 0;

    return 0;
  };

  BP_Gen_Vertex *BP_Gen_Graph::index_to_vert(unsigned long int i1) {
    return &Vertex_list[i1];
  };

  BP_Gen_Edge *BP_Gen_Graph::index_to_edge(unsigned long int i1) {
    return &Edge_list[i1];
  };

  unsigned long int BP_Gen_Graph::vert_to_index(BP_Gen_Vertex *i1) {
    return i1->vertex_list_index;
  };

  unsigned long int BP_Gen_Graph::edge_to_index(BP_Gen_Edge *i1) {
    return i1->edge_list_index;
  };

  unsigned long int BP_Gen_Graph::member_vert_to_index(BP_Gen_GVec_Member *i1) {
    return get_index(i1);
  };

  unsigned long int BP_Gen_Graph::member_edge_to_index(BP_Gen_GVec_Member *i1) {
    return get_index(i1);
  };

  BP_Gen_Vertex *BP_Gen_Graph::member_vert_to_vert(BP_Gen_GVec_Member *i1) {
    return &Vertex_list[get_index(i1)];
  };

  BP_Gen_Edge *BP_Gen_Graph::member_edge_to_edge(BP_Gen_GVec_Member *i1) {
    return &Edge_list[get_index(i1)];
  };

  //  get the index in Vertex_list or Edge_list for i1
  unsigned long int BP_Gen_Graph::get_index(BP_Gen_GVec_Member *i1) {
    for(unsigned long int i = 0; i < i1->graphs.size(); i++)
      if(i1->graphs[i].address == this) {
        return i1->graphs[i].index;
      }

    return std::numeric_limits<unsigned long int>::max();

  };



  unsigned long int BP_Gen_Graph::num_incident_edges(BP_Gen_GVec_Member *i1) {
    return num_incident_edges(get_index(i1));
  };

  unsigned long int BP_Gen_Graph::num_incident_edges(BP_Gen_GVec_Member *i1, int i2) {
    return num_incident_edges(get_index(i1), i2);
  };

  //return number of edges connected to vertex i1 (any edge)
  unsigned long int BP_Gen_Graph::num_incident_edges(unsigned long int i1) {
    return Vertex_list[i1].edges.size();

  };

  //return number of edges connected to vertex i1 considering only(i2= -1:out edges, 0:bidir, 1:in)
  unsigned long int BP_Gen_Graph::num_incident_edges(unsigned long int i1, int i2) {

    BP_Gen_Vertex *v = &Vertex_list[i1];

    unsigned long int i, sum = 0;

    if(i2 == 0) {
      for(i = 0; i < v->edges.size(); i++) {
        if(v->edges[i].address->dir == 0)
          sum++;
      }
      return sum;
    }
    else if(i2 == -1) {
      BP_Gen_Edge *e;
      for(i = 0; i < v->edges.size(); i++) {
        e = v->edges[i].address;
        if(e->dir == 1 & e->vert_list[0] == v)
          sum++;
      }
      return sum;
    }
    else if(i2 == 1) {
      BP_Gen_Edge *e;
      for(i = 0; i < v->edges.size(); i++) {
        e = v->edges[i].address;
        if(e->dir == 1 & e->vert_list[1] == v)
          sum++;
      }
      return sum;
    }
    //else

    return -1;
  };

  //return number of edges connected to vertex i1 (any edge)
  unsigned long int BP_Gen_Graph::num_incident_edges(BP_Gen_Vertex *i1) {
    return i1->edges.size();

  };

  unsigned long int BP_Gen_Graph::num_incident_edges(BP_Gen_Vertex *i1, int i2) {
    return num_incident_edges(i1->vertex_list_index, i2);
  };



  BP_Gen_Vertex *BP_Gen_Graph::add_vertex(BP_Gen_GVec_Member *i1) {
    //if( get_graph_index(i1) != -1) return;
    //check if already in graph
    for(unsigned long int i = 0; i < i1->graphs.size(); i++)
      if(i1->graphs[i].address == this) return 0;


    //Vertex_list_index.add(i1->graphs.size());
    //i1->graph_index.add(Vertex_list.size());
    //Vertex_list.add(i1);
    //i1->graphs.add(this);

    BP_Gen_Vertex *v = Vertex_list.add();
    v->val = i1;
    v->vertex_list_index = Vertex_list.size() - 1;

    BP_Graph_Data_class *gd = i1->graphs.add();
    gd->address = this;
    gd->type = 0;
    gd->index = v->vertex_list_index;

    return v;

  };



  void BP_Gen_Graph::remove(BP_Gen_GVec_Member *i1) {
    for(unsigned long int i = 0; i < i1->graphs.size(); i++)
      if(i1->graphs[i].address == this) {
        if(i1->graphs[i].type == 0)
          remove_vertex(i1->graphs[i].index);
        else
          remove_edge(i1->graphs[i].index);
      }

  };

  void BP_Gen_Graph::remove_vertex(BP_Gen_GVec_Member *i1) {
    remove_vertex(get_index(i1));
  };

  void BP_Gen_Graph::remove_vertex(unsigned long int i1) {
    remove_vertex(&Vertex_list[i1]);
  };

  void BP_Gen_Graph::remove_vertex(BP_Gen_Vertex *i1) {
    unsigned long int i, vi;
    //remove this graph from BP_Gen_GVec_Member
    BP_Gen_GVec_Member *gm = i1->val;
    for(i = 0; i < gm->graphs.size(); i++) {
      if(gm->graphs[i].address == this) {
        gm->graphs.remove(i);
        break;
      }
    }

    //remove edges connecting to this vertex
    while(i1->edges.size() > 0)
      remove_edge(i1->edges[0].address);

    //remove this vertex from Vertex_list
    vi = i1->vertex_list_index;
    Vertex_list.remove(vi);

    //update vertex_list_index of vertex now at vi
    if(vi < Vertex_list.size()) {
      Vertex_list[vi].vertex_list_index = vi;
      gm = Vertex_list[vi].val;
      for(i = 0; i < gm->graphs.size(); i++) {
        if(gm->graphs[i].address == this) {
          gm->graphs[i].index = vi;
          break;
        }
      }
    }
  };



  void BP_Gen_Graph::add_halfedge(BP_Gen_GVec_Member *i1,  BP_Gen_GVec_Member *i3, int i4) {
    add_halfedge(get_index(i1), i3, i4);
  };

  //									Vertex1,  Edge                  , dir (0:undirected, 1: i1->i2)
  void BP_Gen_Graph::add_halfedge(unsigned long int i1,  BP_Gen_GVec_Member *i3, int i4) {
    add_halfedge(&Vertex_list[i1], i3, i4);
  };

  //									Vertex1		,		Edge        , dir (0:undirected, 1: i1->i2)
  void BP_Gen_Graph::add_halfedge(BP_Gen_Vertex *i1, BP_Gen_GVec_Member *i3, int i4) {
    if((i4 != 0) && (i4 != 1)) {
      std::cout << "connect_vertices direction error: " << i4 << std::endl;
      abort();
    }

    //check if edge already in graph
    for(unsigned long int i = 0; i < i3->graphs.size(); i++)
      if(i3->graphs[i].address == this) return;

    // add edge into the Edge_list, then add it
    BP_Gen_Edge *e = Edge_list.add();
    e->edge_list_index = Edge_list.size() - 1;
    e->val = i3;
    e->dir = i4;
    e->vert_list[0] = i1;
    e->vert_list[1] = NULL;
    e->edge_index[0] = i1->edges.size();
    e->edge_index[1] = std::numeric_limits<unsigned long int>::max();

    //add info on this edge to vertex i1
    BP_Edge_Data_class *ed = i1->edges.add();
    ed->address = e;
    ed->vert_index = 0;

    // add this graph to i3
    BP_Graph_Data_class *gd = i3->graphs.add();
    gd->address = this;
    gd->type = 1;
    gd->index = e->edge_list_index;


  };


  void BP_Gen_Graph::connect_vertices(BP_Gen_GVec_Member *i1,  BP_Gen_GVec_Member *i2,    BP_Gen_GVec_Member *i3, int i4) {
    connect_vertices(get_index(i1), get_index(i2), i3, i4);
  };

  //									Vertex1, Vertex2,  Edge                  , dir (0:undirected, 1: i1->i2)
  void BP_Gen_Graph::connect_vertices(unsigned long int i1,  unsigned long int i2,    BP_Gen_GVec_Member *i3, int i4) {
    connect_vertices(&Vertex_list[i1], &Vertex_list[i2], i3, i4);
  };

  //									Vertex1		,		Vertex2		,   Edge        , dir (0:undirected, 1: i1->i2)
  void BP_Gen_Graph::connect_vertices(BP_Gen_Vertex *i1, BP_Gen_Vertex *i2, BP_Gen_GVec_Member *i3, int i4) {
    if((i4 != 0) && (i4 != 1)) {
      std::cout << "connect_vertices direction error: " << i4 << std::endl;
      abort();
    }

    //check if edge already in graph
    for(unsigned long int i = 0; i < i3->graphs.size(); i++)
      if(i3->graphs[i].address == this) return;

    // add edge into the Edge_list, then add it
    BP_Gen_Edge *e = Edge_list.add();
    e->edge_list_index = Edge_list.size() - 1;
    e->val = i3;
    e->dir = i4;
    e->vert_list[0] = i1;
    e->vert_list[1] = i2;
    e->edge_index[0] = i1->edges.size();
    e->edge_index[1] = i2->edges.size();

    //add info on this edge to vertices i1 and i2
    BP_Edge_Data_class *ed = i1->edges.add();
    ed->address = e;
    ed->vert_index = 0;

    ed = i2->edges.add();
    ed->address = e;
    ed->vert_index = 1;

    // add this graph to i3
    BP_Graph_Data_class *gd = i3->graphs.add();
    gd->address = this;
    gd->type = 1;
    gd->index = e->edge_list_index;


  };



  void BP_Gen_Graph::disconnect_vertex(BP_Gen_GVec_Member *i1) {
    disconnect_vertex(get_index(i1));
  };

  void BP_Gen_Graph::disconnect_vertex(unsigned long int i1) {
    disconnect_vertex(&Vertex_list[i1]);
  };

  void BP_Gen_Graph::disconnect_vertex(BP_Gen_Vertex *i1) {
    while(i1->edges.size() > 0)
      remove_edge(i1->edges[0].address);

  };



  void BP_Gen_Graph::disconnect_vertices(BP_Gen_GVec_Member *i1, BP_Gen_GVec_Member *i2) {
    disconnect_vertices(get_index(i1), get_index(i2));
  };

  void BP_Gen_Graph::disconnect_vertices(unsigned long int i1, unsigned long int i2) {
    disconnect_vertices(&Vertex_list[i1], &Vertex_list[i2]);
  };

  void BP_Gen_Graph::disconnect_vertices(BP_Gen_Vertex *i1, BP_Gen_Vertex *i2) {
    // First collect all the edges going between i1 and i2, store in v_edge
    BP_Vec< BP_Gen_Edge *> v_edge;
    BP_Gen_Edge *e;
    for(unsigned long int i = 0; i < i1->edges.size(); i++) {
      e = i1->edges[i].address;
      if(((e->vert_list[0] == i1) & (e->vert_list[1] == i2)) || ((e->vert_list[0] == i2) & (e->vert_list[1] == i1))) {
        v_edge.add(e);
      }
    }

    // Then remove all those edges
    for(unsigned long int i = 0; i < v_edge.size(); i++)
      remove_edge(v_edge[i]);
  };



  void BP_Gen_Graph::remove_edge(BP_Gen_GVec_Member *i1) {
    remove_edge(get_index(i1));
  };

  void BP_Gen_Graph::remove_edge(unsigned long int i1) {
    remove_edge(&Edge_list[i1]);
  };

  void BP_Gen_Graph::remove_edge(BP_Gen_Edge *i1) {
    unsigned long int i, ei, vi;
    //remove this graph from BP_Gen_GVec_Member
    BP_Gen_GVec_Member *gm = i1->val;
    for(i = 0; i < gm->graphs.size(); i++) {
      if(gm->graphs[i].address == this) {
        gm->graphs.remove(i);
        break;
      }
    }

    //remove this edge from it's vertices' "edges" lists
    for(i = 0; i < 2; i++) {
      BP_Gen_Vertex *v = i1->vert_list[i];
      if(v != NULL) {	// allow for halfedges
        vi = i1->edge_index[i];

        v->edges.remove(vi);

        //update edge_index for edge now at vi
        if(vi < v->edges.size()) {
          v->edges[vi].address->edge_index[ v->edges[vi].vert_index] = vi;
        }
      }
    }
    //remove this edge from Edge_list
    ei = i1->edge_list_index;
    Edge_list.remove(ei);

    //update edge_list_index of edge now at ei
    if(ei < Edge_list.size()) {
      Edge_list[ei].edge_list_index = ei;
      gm = Edge_list[ei].val;
      for(i = 0; i < gm->graphs.size(); i++) {
        if(gm->graphs[i].address == this) {
          gm->graphs[i].index = ei;
          break;
        }
      }
    }


  };



  bool BP_Gen_Graph::are_1NN(unsigned long int i1, unsigned long int i2) {
    return are_1NN(&Vertex_list[i1], &Vertex_list[i2]);
  };

  bool BP_Gen_Graph::are_1NN(BP_Gen_GVec_Member *i1, BP_Gen_GVec_Member *i2) {
    return are_1NN(get_index(i1), get_index(i2));
  };

  bool BP_Gen_Graph::are_1NN(BP_Gen_Vertex *i1, BP_Gen_Vertex *i2) {
    unsigned long int i;
    for(i = 0; i < num_incident_edges(i1); i++) {
      if(i2 == get_neighbor(i1, i)) {
        return 1;
      }
    }
    return 0;
  };



  unsigned long int BP_Gen_Graph::edge_startvert(unsigned long int i1) {
    // give an edge, return startpoint
    return vert_to_index(Edge_list[i1].vert_list[0]);
  };

  BP_Gen_Vertex *BP_Gen_Graph::edge_startvert(BP_Gen_Edge *i1) {
    // give an edge, return startpoint
    return i1->vert_list[0];
  };



  unsigned long int BP_Gen_Graph::edge_endvert(unsigned long int i1) {
    // give an edge, return endpoint
    return vert_to_index(Edge_list[i1].vert_list[1]);
  };

  BP_Gen_Vertex *BP_Gen_Graph::edge_endvert(BP_Gen_Edge *i1) {
    // give an edge, return endpoint
    return i1->vert_list[1];
  };



  void BP_Gen_Graph::clear_edges() {
    while(Edge_list.size() > 0) {
      remove_edge(&Edge_list[0]);
    }
  };

  void BP_Gen_Graph::clear() {
    clear_edges();
    while(Vertex_list.size() > 0) {
      remove_vertex(&Vertex_list[0]);
    }
  };

  void BP_Gen_Graph::crazy_clear() {
    Vertex_list.clear();
    Edge_list.clear();
  };

  void BP_Gen_Graph::free() {
    Vertex_list.free();
    Edge_list.free();
  };

  void BP_Gen_Graph::capacity(unsigned long int i1) {
    Vertex_list.capacity(i1);
    Edge_list.capacity(i1);
  };

  void BP_Gen_Graph::cap_increment(unsigned long int i1) {
    Vertex_list.cap_increment(i1);
    Edge_list.cap_increment(i1);
  };

  void BP_Gen_Graph::write(std::ostream &sout) {
    sout << "BP_Graph: " << Name << " contains: " << std::endl;
    sout << "  Vertex: Neighbor_list" << std::endl;
    for(unsigned long int i = 0; i < vert_size(); i++) {
      sout << "    " << i << ": ";
      for(unsigned long int j = 0; j < num_incident_edges(i); j++)
        sout << get_neighbor(i, j) << " ";
      sout << std::endl;
    }

  }

  bool BP_Gen_Graph::contains(BP_Gen_GVec_Member *i1) {
    for(unsigned long int i = 0; i < i1->graphs.size(); i++)
      if(i1->graphs[i].address == this) {
        return true;
      }

    return false;
  };

}

#endif // BP_GVec_CC

